#include "methane_fit.hpp"
#include "methane_fit_internal.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_regression {
namespace {

constexpr double gas_constant = 8.31446261815324;
using internal::doubles;
using internal::parse_payload;
using internal::parse_row;
using internal::Payload;
using internal::positive_finite;
using internal::Row;
using internal::parameter_count;
using internal::residual_count;
using internal::residuals_per_row;
using internal::row_count;
using internal::variable_count;

struct Phase final {
    double volume;
    double pressure;
    double chemical_potential;
    double stability_slope;
    std::array<double, 25> hessian;
    std::string fingerprint;
};

struct RowDiagnostic final {
    Row source;
    Phase liquid;
    Phase vapor;
    std::array<double, residuals_per_row> raw;
    std::array<double, residuals_per_row> scaled;
};

struct Evaluation final {
    std::array<double, residual_count> residuals{};
    std::array<double, residual_count * variable_count> jacobian{};
    std::array<RowDiagnostic, row_count> diagnostics;
    std::string fingerprint;
};

struct MatrixDiagnostics final {
    std::vector<double> singular_values;
    int rank;
    double condition_number;
};

struct SolveOutcome final {
    ceres::Solver::Summary summary;
    std::array<double, variable_count> variables;
    Evaluation evaluation;
    MatrixDiagnostics full_jacobian;
    MatrixDiagnostics parameter_jacobian;
    bool complete_columns;
};

struct ReportingOutcome final {
    Row row;
    double predicted_pressure;
    double predicted_liquid_density;
    Phase liquid;
    Phase vapor;
    std::array<double, 3> raw_residuals;
    std::string termination;
    bool usable;
};

Phase evaluate_phase(
    const epcsaft_native_sdk_v1& table,
    double temperature,
    double amount,
    double volume,
    const std::array<double, 3>& parameters
) {
    epcsaft_parameterized_phase_block_result_v1 result{};
    result.struct_size = sizeof(result);
    const int status = table.evaluate_pure_phase_parameters(
        table.model_context,
        temperature,
        amount,
        volume,
        parameters[0],
        parameters[1],
        parameters[2],
        &result
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
        throw std::runtime_error(std::string("provider evaluation failed: ") + result.error);
    }
    Phase phase{};
    phase.volume = volume;
    phase.pressure = result.pressure_pa;
    phase.chemical_potential = result.chemical_potential_over_rt;
    std::copy(std::begin(result.hessian), std::end(result.hessian), phase.hessian.begin());
    phase.stability_slope = gas_constant * temperature * result.hessian[6] * volume * volume / amount;
    phase.fingerprint.assign(
        result.parameter_fingerprint,
        strnlen(result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE)
    );
    if (!std::isfinite(phase.pressure) || !std::isfinite(phase.chemical_potential)
        || !positive_finite(phase.stability_slope)
        || !std::all_of(phase.hessian.begin(), phase.hessian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::runtime_error("provider phase result is nonfinite or mechanically unstable");
    }
    return phase;
}

Evaluation evaluate_problem(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const std::array<double, variable_count>& variables
) {
    Evaluation evaluation{};
    evaluation.jacobian.fill(std::numeric_limits<double>::quiet_NaN());
    std::array<double, parameter_count> parameters{};
    for (std::size_t index = 0; index < parameter_count; ++index) {
        parameters[index] = payload.start[index] + payload.parameter_scale[index] * variables[index];
        if (!std::isfinite(parameters[index]) || parameters[index] < payload.lower[index]
            || parameters[index] > payload.upper[index]) {
            throw std::runtime_error("transformed parameter is outside its physical bounds");
        }
    }
    for (std::size_t row_index = 0; row_index < row_count; ++row_index) {
        const Row& row = payload.rows[row_index];
        const double liquid_start = payload.molar_mass / row.liquid_density;
        const double vapor_start = gas_constant * row.temperature / row.pressure;
        const double liquid_volume = liquid_start * std::exp(variables[3 + 2 * row_index]);
        const double vapor_volume = vapor_start * std::exp(variables[4 + 2 * row_index]);
        if (!std::isfinite(liquid_volume) || !std::isfinite(vapor_volume)
            || liquid_volume < payload.liquid_volume_bounds[0]
            || liquid_volume > payload.liquid_volume_bounds[1]
            || vapor_volume < payload.vapor_volume_bounds[0]
            || vapor_volume > payload.vapor_volume_bounds[1]
            || liquid_volume >= vapor_volume
            || (vapor_volume - liquid_volume) / vapor_volume <= payload.topology_separation) {
            throw std::runtime_error("phase volume bounds, ordering, or topology separation failed");
        }
        Phase liquid = evaluate_phase(table, row.temperature, payload.amount, liquid_volume, parameters);
        Phase vapor = evaluate_phase(table, row.temperature, payload.amount, vapor_volume, parameters);
        if (liquid.fingerprint.empty() || liquid.fingerprint != vapor.fingerprint
            || (!evaluation.fingerprint.empty() && evaluation.fingerprint != liquid.fingerprint)) {
            throw std::runtime_error("provider source fingerprint changed within the compiled problem");
        }
        evaluation.fingerprint = liquid.fingerprint;
        const double density = payload.molar_mass / liquid_volume;
        const std::array<double, 4> raw{
            liquid.pressure - row.pressure,
            vapor.pressure - row.pressure,
            liquid.chemical_potential - vapor.chemical_potential,
            density - row.liquid_density,
        };
        const std::array<double, 4> scales{row.pressure, row.pressure, 1.0, row.liquid_density};
        std::array<double, 4> scaled{};
        for (std::size_t residual = 0; residual < residuals_per_row; ++residual) {
            scaled[residual] = std::sqrt(payload.weights[residual]) * raw[residual] / scales[residual];
            evaluation.residuals[row_index * residuals_per_row + residual] = scaled[residual];
        }

        auto set_jacobian = [&](std::size_t local_row, std::size_t column, double value) {
            evaluation.jacobian[
                (row_index * residuals_per_row + local_row) * variable_count + column
            ] = value;
        };
        for (std::size_t local_row = 0; local_row < residuals_per_row; ++local_row) {
            for (std::size_t column = 0; column < variable_count; ++column) {
                set_jacobian(local_row, column, 0.0);
            }
        }
        const double pressure_factor = std::sqrt(payload.weights[0]) / row.pressure;
        const double vapor_pressure_factor = std::sqrt(payload.weights[1]) / row.pressure;
        const double mu_factor = std::sqrt(payload.weights[2]);
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            const std::size_t coordinate = 2 + parameter;
            const double liquid_dp = -gas_constant * row.temperature
                * liquid.hessian[5 + coordinate];
            const double vapor_dp = -gas_constant * row.temperature
                * vapor.hessian[5 + coordinate];
            const double liquid_dmu = liquid.hessian[coordinate];
            const double vapor_dmu = vapor.hessian[coordinate];
            set_jacobian(0, parameter, pressure_factor * liquid_dp * payload.parameter_scale[parameter]);
            set_jacobian(1, parameter, vapor_pressure_factor * vapor_dp * payload.parameter_scale[parameter]);
            set_jacobian(
                2,
                parameter,
                mu_factor * (liquid_dmu - vapor_dmu) * payload.parameter_scale[parameter]
            );
        }
        const std::size_t liquid_column = 3 + 2 * row_index;
        const std::size_t vapor_column = liquid_column + 1;
        const double liquid_dpdv = -gas_constant * row.temperature * liquid.hessian[6];
        const double vapor_dpdv = -gas_constant * row.temperature * vapor.hessian[6];
        set_jacobian(0, liquid_column, pressure_factor * liquid_dpdv * liquid_volume);
        set_jacobian(1, vapor_column, vapor_pressure_factor * vapor_dpdv * vapor_volume);
        set_jacobian(2, liquid_column, mu_factor * liquid.hessian[1] * liquid_volume);
        set_jacobian(2, vapor_column, -mu_factor * vapor.hessian[1] * vapor_volume);
        set_jacobian(
            3,
            liquid_column,
            std::sqrt(payload.weights[3]) * (-density) / row.liquid_density
        );
        evaluation.diagnostics[row_index] = RowDiagnostic{
            row, std::move(liquid), std::move(vapor), raw, scaled
        };
    }
    if (!std::all_of(evaluation.residuals.begin(), evaluation.residuals.end(), [](double value) {
            return std::isfinite(value);
        }) || !std::all_of(evaluation.jacobian.begin(), evaluation.jacobian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::runtime_error("assembled residual or Jacobian is nonfinite");
    }
    if (evaluation.fingerprint != payload.identity[36]) {
        throw std::runtime_error("provider source fingerprint differs from the compiled problem");
    }
    return evaluation;
}

class TrainingCost final : public ceres::CostFunction {
public:
    TrainingCost(const epcsaft_native_sdk_v1& table, const Payload& payload)
        : table_(table), payload_(payload) {
        set_num_residuals(static_cast<int>(residual_count));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(parameter_count));
        for (std::size_t index = 0; index < row_count; ++index) {
            mutable_parameter_block_sizes()->push_back(2);
        }
    }

    bool Evaluate(double const* const* blocks, double* residuals, double** jacobians) const override {
        try {
            std::array<double, variable_count> variables{};
            std::copy(blocks[0], blocks[0] + parameter_count, variables.begin());
            for (std::size_t row = 0; row < row_count; ++row) {
                variables[3 + 2 * row] = blocks[row + 1][0];
                variables[4 + 2 * row] = blocks[row + 1][1];
            }
            const Evaluation evaluation = evaluate_problem(table_, payload_, variables);
            std::copy(evaluation.residuals.begin(), evaluation.residuals.end(), residuals);
            if (jacobians != nullptr) {
                if (jacobians[0] != nullptr) {
                    for (std::size_t row = 0; row < residual_count; ++row) {
                        for (std::size_t column = 0; column < parameter_count; ++column) {
                            jacobians[0][row * parameter_count + column] =
                                evaluation.jacobian[row * variable_count + column];
                        }
                    }
                }
                for (std::size_t block = 0; block < row_count; ++block) {
                    if (jacobians[block + 1] == nullptr) continue;
                    for (std::size_t row = 0; row < residual_count; ++row) {
                        jacobians[block + 1][row * 2] =
                            evaluation.jacobian[row * variable_count + 3 + 2 * block];
                        jacobians[block + 1][row * 2 + 1] =
                            evaluation.jacobian[row * variable_count + 4 + 2 * block];
                    }
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    const epcsaft_native_sdk_v1& table_;
    const Payload& payload_;
};

std::string termination_name(ceres::TerminationType termination) {
    switch (termination) {
        case ceres::CONVERGENCE: return "CONVERGENCE";
        case ceres::NO_CONVERGENCE: return "NO_CONVERGENCE";
        case ceres::FAILURE: return "FAILURE";
        case ceres::USER_SUCCESS: return "USER_SUCCESS";
        case ceres::USER_FAILURE: return "USER_FAILURE";
    }
    return "UNKNOWN";
}

MatrixDiagnostics matrix_diagnostics(const Eigen::MatrixXd& matrix) {
    const Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto singular = decomposition.singularValues();
    MatrixDiagnostics diagnostics{};
    diagnostics.singular_values.reserve(static_cast<std::size_t>(singular.size()));
    for (Eigen::Index index = 0; index < singular.size(); ++index) {
        diagnostics.singular_values.push_back(singular(index));
    }
    const double maximum = singular.size() == 0 ? 0.0 : singular(0);
    const double threshold = maximum * static_cast<double>(std::max(matrix.rows(), matrix.cols()))
        * std::numeric_limits<double>::epsilon() * 100.0;
    diagnostics.rank = 0;
    double minimum_accepted = std::numeric_limits<double>::infinity();
    for (double value : diagnostics.singular_values) {
        if (value > threshold) {
            ++diagnostics.rank;
            minimum_accepted = std::min(minimum_accepted, value);
        }
    }
    diagnostics.condition_number = diagnostics.rank == 0
        ? std::numeric_limits<double>::infinity()
        : maximum / minimum_accepted;
    return diagnostics;
}

SolveOutcome solve_training(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    double liquid_start_offset,
    double vapor_start_offset
) {
    SolveOutcome outcome{};
    outcome.variables.fill(0.0);
    std::array<double, 3> parameters{};
    std::array<std::array<double, 2>, row_count> volumes{};
    for (std::size_t row = 0; row < row_count; ++row) {
        volumes[row] = {liquid_start_offset, vapor_start_offset};
    }
    ceres::Problem problem;
    auto* cost = new TrainingCost(table, payload);
    std::vector<double*> blocks{parameters.data()};
    for (auto& volume : volumes) blocks.push_back(volume.data());
    problem.AddResidualBlock(cost, nullptr, blocks);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        problem.SetParameterLowerBound(
            parameters.data(), static_cast<int>(parameter),
            (payload.lower[parameter] - payload.start[parameter]) / payload.parameter_scale[parameter]
        );
        problem.SetParameterUpperBound(
            parameters.data(), static_cast<int>(parameter),
            (payload.upper[parameter] - payload.start[parameter]) / payload.parameter_scale[parameter]
        );
    }
    for (std::size_t row = 0; row < row_count; ++row) {
        const double liquid_reference = payload.molar_mass / payload.rows[row].liquid_density;
        const double vapor_reference = gas_constant * payload.rows[row].temperature / payload.rows[row].pressure;
        problem.SetParameterLowerBound(
            volumes[row].data(), 0, std::log(payload.liquid_volume_bounds[0] / liquid_reference)
        );
        problem.SetParameterUpperBound(
            volumes[row].data(), 0, std::log(payload.liquid_volume_bounds[1] / liquid_reference)
        );
        problem.SetParameterLowerBound(
            volumes[row].data(), 1, std::log(payload.vapor_volume_bounds[0] / vapor_reference)
        );
        problem.SetParameterUpperBound(
            volumes[row].data(), 1, std::log(payload.vapor_volume_bounds[1] / vapor_reference)
        );
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = payload.num_threads;
    ceres::Solve(options, &problem, &outcome.summary);
    std::copy(parameters.begin(), parameters.end(), outcome.variables.begin());
    for (std::size_t row = 0; row < row_count; ++row) {
        outcome.variables[3 + 2 * row] = volumes[row][0];
        outcome.variables[4 + 2 * row] = volumes[row][1];
    }
    outcome.evaluation = evaluate_problem(table, payload, outcome.variables);
    Eigen::MatrixXd full(static_cast<Eigen::Index>(residual_count), static_cast<Eigen::Index>(variable_count));
    Eigen::MatrixXd parameter_matrix(
        static_cast<Eigen::Index>(residual_count), static_cast<Eigen::Index>(parameter_count)
    );
    outcome.complete_columns = true;
    for (std::size_t row = 0; row < residual_count; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            const double value = outcome.evaluation.jacobian[row * variable_count + column];
            full(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = value;
            if (column < parameter_count) {
                parameter_matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = value;
            }
        }
    }
    outcome.complete_columns = std::all_of(
        outcome.evaluation.jacobian.begin(),
        outcome.evaluation.jacobian.end(),
        [](double value) { return std::isfinite(value); }
    );
    outcome.full_jacobian = matrix_diagnostics(full);
    outcome.parameter_jacobian = matrix_diagnostics(parameter_matrix);
    return outcome;
}

class ReportingCost final : public ceres::SizedCostFunction<3, 3> {
public:
    ReportingCost(
        const epcsaft_native_sdk_v1& table,
        const Payload& payload,
        Row row,
        std::array<double, 3> parameters
    ) : table_(table), payload_(payload), row_(row), parameters_(parameters) {}

    bool Evaluate(double const* const* blocks, double* residuals, double** jacobians) const override {
        try {
            const double liquid_reference = payload_.molar_mass / row_.liquid_density;
            const double vapor_reference = gas_constant * row_.temperature / row_.pressure;
            const double liquid_volume = liquid_reference * std::exp(blocks[0][0]);
            const double vapor_volume = vapor_reference * std::exp(blocks[0][1]);
            const double pressure = row_.pressure * std::exp(blocks[0][2]);
            if (liquid_volume < payload_.liquid_volume_bounds[0]
                || liquid_volume > payload_.liquid_volume_bounds[1]
                || vapor_volume < payload_.vapor_volume_bounds[0]
                || vapor_volume > payload_.vapor_volume_bounds[1]
                || liquid_volume >= vapor_volume
                || (vapor_volume - liquid_volume) / vapor_volume <= payload_.topology_separation) {
                return false;
            }
            const Phase liquid = evaluate_phase(
                table_, row_.temperature, payload_.amount, liquid_volume, parameters_
            );
            const Phase vapor = evaluate_phase(
                table_, row_.temperature, payload_.amount, vapor_volume, parameters_
            );
            residuals[0] = (liquid.pressure - pressure) / row_.pressure;
            residuals[1] = (vapor.pressure - pressure) / row_.pressure;
            residuals[2] = liquid.chemical_potential - vapor.chemical_potential;
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                jacobians[0][0] = -gas_constant * row_.temperature * liquid.hessian[6]
                    * liquid_volume / row_.pressure;
                jacobians[0][1] = 0.0;
                jacobians[0][2] = -pressure / row_.pressure;
                jacobians[0][3] = 0.0;
                jacobians[0][4] = -gas_constant * row_.temperature * vapor.hessian[6]
                    * vapor_volume / row_.pressure;
                jacobians[0][5] = -pressure / row_.pressure;
                jacobians[0][6] = liquid.hessian[1] * liquid_volume;
                jacobians[0][7] = -vapor.hessian[1] * vapor_volume;
                jacobians[0][8] = 0.0;
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    const epcsaft_native_sdk_v1& table_;
    const Payload& payload_;
    Row row_;
    std::array<double, 3> parameters_;
};

ReportingOutcome solve_reporting(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const Row& row,
    const std::array<double, 3>& parameters
) {
    std::array<double, 3> variables{};
    ceres::Problem problem;
    problem.AddResidualBlock(new ReportingCost(table, payload, row, parameters), nullptr, variables.data());
    const double liquid_reference = payload.molar_mass / row.liquid_density;
    const double vapor_reference = gas_constant * row.temperature / row.pressure;
    problem.SetParameterLowerBound(variables.data(), 0, std::log(payload.liquid_volume_bounds[0] / liquid_reference));
    problem.SetParameterUpperBound(variables.data(), 0, std::log(payload.liquid_volume_bounds[1] / liquid_reference));
    problem.SetParameterLowerBound(variables.data(), 1, std::log(payload.vapor_volume_bounds[0] / vapor_reference));
    problem.SetParameterUpperBound(variables.data(), 1, std::log(payload.vapor_volume_bounds[1] / vapor_reference));
    problem.SetParameterLowerBound(
        variables.data(), 2, std::log(payload.reporting_pressure_bounds[0] / row.pressure)
    );
    problem.SetParameterUpperBound(
        variables.data(), 2, std::log(payload.reporting_pressure_bounds[1] / row.pressure)
    );
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = payload.num_threads;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    const double liquid_volume = liquid_reference * std::exp(variables[0]);
    const double vapor_volume = vapor_reference * std::exp(variables[1]);
    const double pressure = row.pressure * std::exp(variables[2]);
    Phase liquid = evaluate_phase(table, row.temperature, payload.amount, liquid_volume, parameters);
    Phase vapor = evaluate_phase(table, row.temperature, payload.amount, vapor_volume, parameters);
    return ReportingOutcome{
        row,
        pressure,
        payload.molar_mass / liquid_volume,
        liquid,
        vapor,
        {
            liquid.pressure - pressure,
            vapor.pressure - pressure,
            liquid.chemical_potential - vapor.chemical_potential,
        },
        termination_name(summary.termination_type),
        summary.IsSolutionUsable(),
    };
}

PyObject* tuple_from_values(const double* values, std::size_t size) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(size));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < size; ++index) {
        PyObject* value = PyFloat_FromDouble(values[index]);
        if (value == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* diagnostics_to_python(const std::array<RowDiagnostic, row_count>& diagnostics) {
    PyObject* rows = PyTuple_New(static_cast<Py_ssize_t>(row_count));
    if (rows == nullptr) return nullptr;
    for (std::size_t index = 0; index < row_count; ++index) {
        const RowDiagnostic& row = diagnostics[index];
        PyObject* raw = tuple_from_values(row.raw.data(), row.raw.size());
        PyObject* scaled = tuple_from_values(row.scaled.data(), row.scaled.size());
        if (raw == nullptr || scaled == nullptr) {
            Py_XDECREF(raw);
            Py_XDECREF(scaled);
            Py_DECREF(rows);
            return nullptr;
        }
        PyObject* item = Py_BuildValue(
            "(ssdddddddddNN)",
            row.source.row_id.c_str(),
            row.source.source_id.c_str(),
            row.source.temperature,
            row.liquid.volume,
            row.vapor.volume,
            row.liquid.pressure,
            row.vapor.pressure,
            row.liquid.chemical_potential,
            row.vapor.chemical_potential,
            row.liquid.stability_slope,
            row.vapor.stability_slope,
            raw,
            scaled
        );
        if (item == nullptr) {
            Py_DECREF(rows);
            return nullptr;
        }
        PyTuple_SET_ITEM(rows, static_cast<Py_ssize_t>(index), item);
    }
    return rows;
}

PyObject* vector_to_tuple(const std::vector<double>& values) {
    return tuple_from_values(values.data(), values.size());
}

PyObject* strings_to_tuple(const std::vector<std::string>& values) {
    PyObject* result = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (result == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyUnicode_FromStringAndSize(
            values[index].data(), static_cast<Py_ssize_t>(values[index].size())
        );
        if (value == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), value);
    }
    return result;
}

PyObject* reporting_to_python(const std::vector<ReportingOutcome>& outcomes) {
    PyObject* rows = PyTuple_New(static_cast<Py_ssize_t>(outcomes.size()));
    if (rows == nullptr) return nullptr;
    for (std::size_t index = 0; index < outcomes.size(); ++index) {
        const ReportingOutcome& outcome = outcomes[index];
        PyObject* raw = tuple_from_values(outcome.raw_residuals.data(), outcome.raw_residuals.size());
        if (raw == nullptr) {
            Py_DECREF(rows);
            return nullptr;
        }
        PyObject* item = Py_BuildValue(
            "(ssdddddddddNsO)",
            outcome.row.row_id.c_str(),
            outcome.row.source_id.c_str(),
            outcome.row.temperature,
            outcome.row.pressure,
            outcome.row.liquid_density,
            outcome.predicted_pressure,
            outcome.predicted_liquid_density,
            outcome.liquid.volume,
            outcome.vapor.volume,
            outcome.liquid.stability_slope,
            outcome.vapor.stability_slope,
            raw,
            outcome.termination.c_str(),
            outcome.usable ? Py_True : Py_False
        );
        if (item == nullptr) {
            Py_DECREF(rows);
            return nullptr;
        }
        PyTuple_SET_ITEM(rows, static_cast<Py_ssize_t>(index), item);
    }
    return rows;
}

}  // namespace

PyObject* evaluate_python(PyObject* capsule, PyObject* payload_object, PyObject* variables_object) {
    const epcsaft_native_sdk_v1* table = checked_provider_table(capsule);
    if (table == nullptr) return nullptr;
    try {
        const Payload payload = parse_payload(payload_object);
        const std::vector<double> parsed_variables = doubles(
            variables_object, variable_count, "transformed variables"
        );
        std::array<double, variable_count> variables{};
        std::copy(parsed_variables.begin(), parsed_variables.end(), variables.begin());
        const Evaluation evaluation = evaluate_problem(*table, payload, variables);
        PyObject* residuals = tuple_from_values(evaluation.residuals.data(), residual_count);
        PyObject* jacobian = tuple_from_values(
            evaluation.jacobian.data(), residual_count * variable_count
        );
        PyObject* diagnostics = diagnostics_to_python(evaluation.diagnostics);
        PyObject* fingerprint = PyUnicode_FromStringAndSize(
            evaluation.fingerprint.data(), static_cast<Py_ssize_t>(evaluation.fingerprint.size())
        );
        if (residuals == nullptr || jacobian == nullptr || diagnostics == nullptr
            || fingerprint == nullptr) {
            Py_XDECREF(residuals);
            Py_XDECREF(jacobian);
            Py_XDECREF(diagnostics);
            Py_XDECREF(fingerprint);
            return nullptr;
        }
        PyObject* result = PyTuple_New(4);
        if (result == nullptr) {
            Py_DECREF(residuals);
            Py_DECREF(jacobian);
            Py_DECREF(diagnostics);
            Py_DECREF(fingerprint);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, 0, residuals);
        PyTuple_SET_ITEM(result, 1, jacobian);
        PyTuple_SET_ITEM(result, 2, diagnostics);
        PyTuple_SET_ITEM(result, 3, fingerprint);
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_ValueError, error.what());
        return nullptr;
    }
}

PyObject* solve_python(PyObject* capsule, PyObject* payload_object, PyObject* reporting_rows_object) {
    const epcsaft_native_sdk_v1* table = checked_provider_table(capsule);
    if (table == nullptr) return nullptr;
    try {
        const Payload payload = parse_payload(payload_object);
        const SolveOutcome primary = solve_training(*table, payload, 0.0, 0.0);
        const SolveOutcome confirmation = solve_training(
            *table,
            payload,
            std::log(payload.confirmation_liquid_start_multiplier),
            std::log(payload.confirmation_vapor_start_multiplier)
        );
        double parameter_delta = 0.0;
        for (std::size_t index = 0; index < parameter_count; ++index) {
            parameter_delta = std::max(
                parameter_delta, std::abs(primary.variables[index] - confirmation.variables[index])
            );
        }
        const double cost_delta = std::abs(primary.summary.final_cost - confirmation.summary.final_cost)
            / std::max(1.0, std::abs(primary.summary.final_cost));
        std::array<double, 3> final_parameters{};
        for (std::size_t index = 0; index < parameter_count; ++index) {
            final_parameters[index] = payload.start[index]
                + payload.parameter_scale[index] * primary.variables[index];
        }
        PyObject* reporting_sequence = PySequence_Fast(
            reporting_rows_object, "reporting rows must be a sequence"
        );
        if (reporting_sequence == nullptr || PySequence_Fast_GET_SIZE(reporting_sequence) != 9) {
            Py_XDECREF(reporting_sequence);
            throw std::invalid_argument("reporting rows must contain the ordered 100-180 K table");
        }
        std::vector<ReportingOutcome> reporting;
        reporting.reserve(9);
        for (Py_ssize_t index = 0; index < 9; ++index) {
            const Row row = parse_row(
                PySequence_Fast_GET_ITEM(reporting_sequence, index),
                static_cast<std::size_t>(index)
            );
            reporting.push_back(solve_reporting(
                *table, payload, row, final_parameters
            ));
        }
        Py_DECREF(reporting_sequence);
        PyObject* variables = tuple_from_values(primary.variables.data(), primary.variables.size());
        PyObject* residuals = tuple_from_values(primary.evaluation.residuals.data(), residual_count);
        PyObject* jacobian = tuple_from_values(
            primary.evaluation.jacobian.data(), residual_count * variable_count
        );
        PyObject* diagnostics = diagnostics_to_python(primary.evaluation.diagnostics);
        PyObject* full_singular = vector_to_tuple(primary.full_jacobian.singular_values);
        PyObject* parameter_singular = vector_to_tuple(primary.parameter_jacobian.singular_values);
        PyObject* reporting_python = reporting_to_python(reporting);
        PyObject* compiled_identity = strings_to_tuple(payload.identity);
        if (variables == nullptr || residuals == nullptr || jacobian == nullptr
            || diagnostics == nullptr || full_singular == nullptr
            || parameter_singular == nullptr || reporting_python == nullptr
            || compiled_identity == nullptr) {
            Py_XDECREF(variables);
            Py_XDECREF(residuals);
            Py_XDECREF(jacobian);
            Py_XDECREF(diagnostics);
            Py_XDECREF(full_singular);
            Py_XDECREF(parameter_singular);
            Py_XDECREF(reporting_python);
            Py_XDECREF(compiled_identity);
            return nullptr;
        }
        PyObject* result = PyTuple_New(24);
        if (result == nullptr) {
            Py_DECREF(variables);
            Py_DECREF(residuals);
            Py_DECREF(jacobian);
            Py_DECREF(diagnostics);
            Py_DECREF(full_singular);
            Py_DECREF(parameter_singular);
            Py_DECREF(reporting_python);
            Py_DECREF(compiled_identity);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, 0, PyUnicode_FromString(termination_name(primary.summary.termination_type).c_str()));
        PyTuple_SET_ITEM(result, 1, Py_NewRef(primary.summary.IsSolutionUsable() ? Py_True : Py_False));
        PyTuple_SET_ITEM(result, 2, PyFloat_FromDouble(primary.summary.initial_cost));
        PyTuple_SET_ITEM(result, 3, PyFloat_FromDouble(primary.summary.final_cost));
        PyTuple_SET_ITEM(result, 4, PyLong_FromSize_t(primary.summary.iterations.size()));
        PyTuple_SET_ITEM(result, 5, variables);
        PyTuple_SET_ITEM(result, 6, residuals);
        PyTuple_SET_ITEM(result, 7, jacobian);
        PyTuple_SET_ITEM(result, 8, diagnostics);
        PyTuple_SET_ITEM(result, 9, full_singular);
        PyTuple_SET_ITEM(result, 10, PyLong_FromLong(primary.full_jacobian.rank));
        PyTuple_SET_ITEM(result, 11, PyFloat_FromDouble(primary.full_jacobian.condition_number));
        PyTuple_SET_ITEM(result, 12, parameter_singular);
        PyTuple_SET_ITEM(result, 13, PyLong_FromLong(primary.parameter_jacobian.rank));
        PyTuple_SET_ITEM(result, 14, PyFloat_FromDouble(primary.parameter_jacobian.condition_number));
        PyTuple_SET_ITEM(result, 15, Py_NewRef(primary.complete_columns ? Py_True : Py_False));
        PyTuple_SET_ITEM(result, 16, PyFloat_FromDouble(parameter_delta));
        PyTuple_SET_ITEM(result, 17, PyFloat_FromDouble(cost_delta));
        PyTuple_SET_ITEM(result, 18, PyUnicode_FromString(termination_name(confirmation.summary.termination_type).c_str()));
        PyTuple_SET_ITEM(result, 19, Py_NewRef(confirmation.summary.IsSolutionUsable() ? Py_True : Py_False));
        PyTuple_SET_ITEM(result, 20, reporting_python);
        PyTuple_SET_ITEM(result, 21, PyUnicode_FromString(primary.evaluation.fingerprint.c_str()));
        PyTuple_SET_ITEM(result, 22, Py_True);
        Py_INCREF(Py_True);
        PyTuple_SET_ITEM(result, 23, compiled_identity);
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_regression
