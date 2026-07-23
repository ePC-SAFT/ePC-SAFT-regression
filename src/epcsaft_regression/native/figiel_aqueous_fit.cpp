#include "figiel_aqueous_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_regression {
namespace {

constexpr std::size_t salt_count = 6;

struct Observation final {
    std::string row_id;
    std::string salt;
    std::size_t salt_index;
    double molality;
    double observed;
    std::vector<std::size_t> columns;
};

struct Start final {
    std::string name;
    std::vector<double> values;
};

struct Payload final {
    std::vector<std::string> identity;
    std::string stage;
    std::vector<Observation> observations;
    std::array<std::string, salt_count> expected_fingerprints;
    std::size_t parameter_count;
    std::array<double, 2> bounds;
    std::vector<Start> starts;
    double temperature;
    double pressure;
    int max_iterations;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    double rank_multiplier;
};

struct Row final {
    std::string row_id;
    std::string salt;
    double molality;
    double observed;
    double log_modeled;
    double modeled;
    double residual;
    std::vector<double> local_derivative;
    double reference_molality;
    double reference_convergence;
    double derivative_convergence;
    std::string fingerprint;
};

struct Evaluation final {
    std::vector<double> residuals;
    std::vector<double> jacobian;
    std::vector<Row> rows;
};

struct SolveOutcome final {
    std::string name;
    ceres::Solver::Summary summary;
    std::vector<double> parameters;
    Evaluation evaluation;
    std::vector<double> singular_values;
    double rank_threshold{0.0};
    int rank{0};
    double condition_number{std::numeric_limits<double>::infinity()};
    std::vector<double> least_sensitive_direction;
    bool complete_columns{false};
    std::string failure_reason;
};

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    Py_ssize_t size = 0;
    const char* value = PyUnicode_AsUTF8AndSize(object, &size);
    if (value == nullptr) throw std::invalid_argument(std::string(label) + " must be UTF-8");
    return std::string(value, static_cast<std::size_t>(size));
}

double number(PyObject* object, const char* label) {
    const double value = PyFloat_AsDouble(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(label) + " must be numeric");
    }
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(label) + " must be finite");
    }
    return value;
}

std::size_t index_value(PyObject* object, const char* label) {
    const unsigned long long value = PyLong_AsUnsignedLongLong(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(label) + " must be a nonnegative integer");
    }
    return static_cast<std::size_t>(value);
}

std::vector<double> doubles(PyObject* object, std::size_t expected, const char* label) {
    PyObject* sequence = PySequence_Fast(object, label);
    if (sequence == nullptr) throw std::invalid_argument(label);
    if (PySequence_Fast_GET_SIZE(sequence) != static_cast<Py_ssize_t>(expected)) {
        Py_DECREF(sequence);
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::vector<double> result(expected);
    for (std::size_t index = 0; index < expected; ++index) {
        result[index] = number(
            PySequence_Fast_GET_ITEM(sequence, static_cast<Py_ssize_t>(index)), label
        );
    }
    Py_DECREF(sequence);
    return result;
}

Payload parse_payload(PyObject* object) {
    PyObject* sequence = PySequence_Fast(object, "aqueous payload must be a sequence");
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 14) {
        Py_XDECREF(sequence);
        throw std::invalid_argument("aqueous payload has the wrong length");
    }
    Payload payload{};
    PyObject* identity = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 0), "aqueous identity must be a sequence"
    );
    if (identity == nullptr || PySequence_Fast_GET_SIZE(identity) == 0) {
        Py_XDECREF(identity);
        Py_DECREF(sequence);
        throw std::invalid_argument("aqueous identity must be nonempty");
    }
    for (Py_ssize_t index = 0; index < PySequence_Fast_GET_SIZE(identity); ++index) {
        payload.identity.push_back(text(PySequence_Fast_GET_ITEM(identity, index), "identity"));
    }
    Py_DECREF(identity);
    payload.stage = text(PySequence_Fast_GET_ITEM(sequence, 1), "stage");
    if (payload.stage != "solvation_factor" && payload.stage != "aqueous_kij") {
        Py_DECREF(sequence);
        throw std::invalid_argument("aqueous stage is unsupported");
    }
    PyObject* observations = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 2), "observations must be a sequence"
    );
    if (observations == nullptr) {
        Py_DECREF(sequence);
        throw std::invalid_argument("observations must be a sequence");
    }
    const Py_ssize_t row_count = PySequence_Fast_GET_SIZE(observations);
    const Py_ssize_t expected_rows = payload.stage == "solvation_factor" ? 21 : 164;
    if (row_count != expected_rows) {
        Py_DECREF(observations);
        Py_DECREF(sequence);
        throw std::invalid_argument("aqueous stage has the wrong row count");
    }
    payload.observations.reserve(static_cast<std::size_t>(row_count));
    for (Py_ssize_t row_index = 0; row_index < row_count; ++row_index) {
        PyObject* row = PySequence_Fast(
            PySequence_Fast_GET_ITEM(observations, row_index),
            "aqueous observation must be a sequence"
        );
        if (row == nullptr || PySequence_Fast_GET_SIZE(row) != 6) {
            Py_XDECREF(row);
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument("aqueous observation has the wrong length");
        }
        PyObject* columns = PySequence_Fast(
            PySequence_Fast_GET_ITEM(row, 5), "aqueous columns must be a sequence"
        );
        const Py_ssize_t expected_columns = payload.stage == "solvation_factor" ? 1 : 3;
        if (columns == nullptr || PySequence_Fast_GET_SIZE(columns) != expected_columns) {
            Py_XDECREF(columns);
            Py_DECREF(row);
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument("aqueous observation column map is invalid");
        }
        Observation observed{
            text(PySequence_Fast_GET_ITEM(row, 0), "row id"),
            text(PySequence_Fast_GET_ITEM(row, 1), "salt"),
            index_value(PySequence_Fast_GET_ITEM(row, 2), "salt index"),
            number(PySequence_Fast_GET_ITEM(row, 3), "molality"),
            number(PySequence_Fast_GET_ITEM(row, 4), "observed MIAC"),
            {},
        };
        for (Py_ssize_t column = 0; column < expected_columns; ++column) {
            observed.columns.push_back(index_value(
                PySequence_Fast_GET_ITEM(columns, column), "column index"
            ));
        }
        Py_DECREF(columns);
        Py_DECREF(row);
        if (observed.salt_index >= salt_count || observed.observed <= 0.0
            || observed.molality <= 0.0) {
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument("aqueous observation is outside the frozen domain");
        }
        payload.observations.push_back(std::move(observed));
    }
    Py_DECREF(observations);
    PyObject* fingerprints = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 3), "fingerprints must be a sequence"
    );
    if (fingerprints == nullptr || PySequence_Fast_GET_SIZE(fingerprints) != 6) {
        Py_XDECREF(fingerprints);
        Py_DECREF(sequence);
        throw std::invalid_argument("fingerprints must contain six entries");
    }
    for (std::size_t index = 0; index < salt_count; ++index) {
        payload.expected_fingerprints[index] = text(
            PySequence_Fast_GET_ITEM(fingerprints, static_cast<Py_ssize_t>(index)),
            "fingerprint"
        );
    }
    Py_DECREF(fingerprints);
    payload.parameter_count = index_value(
        PySequence_Fast_GET_ITEM(sequence, 4), "parameter count"
    );
    const std::size_t expected_parameters = payload.stage == "solvation_factor" ? 1 : 11;
    if (payload.parameter_count != expected_parameters) {
        Py_DECREF(sequence);
        throw std::invalid_argument("aqueous stage has the wrong parameter count");
    }
    const auto bounds = doubles(PySequence_Fast_GET_ITEM(sequence, 5), 2, "bounds");
    payload.bounds = {bounds[0], bounds[1]};
    PyObject* starts = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 6), "starts must be a sequence"
    );
    if (starts == nullptr || PySequence_Fast_GET_SIZE(starts) == 0) {
        Py_XDECREF(starts);
        Py_DECREF(sequence);
        throw std::invalid_argument("starts must be nonempty");
    }
    for (Py_ssize_t index = 0; index < PySequence_Fast_GET_SIZE(starts); ++index) {
        PyObject* item = PySequence_Fast(
            PySequence_Fast_GET_ITEM(starts, index), "start must be a pair"
        );
        if (item == nullptr || PySequence_Fast_GET_SIZE(item) != 2) {
            Py_XDECREF(item);
            Py_DECREF(starts);
            Py_DECREF(sequence);
            throw std::invalid_argument("start must be a name/value pair");
        }
        payload.starts.push_back(Start{
            text(PySequence_Fast_GET_ITEM(item, 0), "start name"),
            doubles(PySequence_Fast_GET_ITEM(item, 1), payload.parameter_count, "start values"),
        });
        Py_DECREF(item);
    }
    Py_DECREF(starts);
    payload.temperature = number(PySequence_Fast_GET_ITEM(sequence, 7), "temperature");
    payload.pressure = number(PySequence_Fast_GET_ITEM(sequence, 8), "pressure");
    payload.max_iterations = static_cast<int>(index_value(
        PySequence_Fast_GET_ITEM(sequence, 9), "maximum iterations"
    ));
    payload.function_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 10), "function tolerance"
    );
    payload.gradient_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 11), "gradient tolerance"
    );
    payload.parameter_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 12), "parameter tolerance"
    );
    payload.rank_multiplier = number(PySequence_Fast_GET_ITEM(sequence, 13), "rank multiplier");
    Py_DECREF(sequence);
    if (payload.temperature != 298.15 || payload.pressure != 100000.0
        || payload.bounds != (payload.stage == "solvation_factor"
            ? std::array<double, 2>{1.0, 2.0}
            : std::array<double, 2>{-1.0, 1.0})
        || payload.max_iterations != 500 || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10
        || payload.parameter_tolerance != 1.0e-10 || payload.rank_multiplier != 100.0) {
        throw std::invalid_argument("aqueous numerical contract does not match the frozen design");
    }
    return payload;
}

const epcsaft_native_sdk_v1* checked_table(
    PyObject* capsule, const Payload& payload, std::size_t salt_index
) {
    if (!PyCapsule_CheckExact(capsule)) {
        throw std::invalid_argument("provider transport must be an exact CPython capsule");
    }
    void* pointer = PyCapsule_GetPointer(capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME);
    if (pointer == nullptr) throw std::runtime_error("provider capability unavailable");
    struct Prefix final { std::uint32_t abi_version; std::size_t table_size; } prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    const std::size_t minimum_size = payload.stage == "solvation_factor"
        ? offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_solvation_factor)
            + sizeof(epcsaft_evaluate_aqueous_miac_solvation_factor_v1)
        : offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_kij)
            + sizeof(epcsaft_evaluate_aqueous_miac_kij_v1);
    if (prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || prefix.table_size < minimum_size) {
        throw std::runtime_error("provider capability unavailable");
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    constexpr std::array<const char*, salt_count> cations{
        "lithium-cation", "sodium-cation", "potassium-cation",
        "lithium-cation", "sodium-cation", "potassium-cation",
    };
    constexpr std::array<const char*, salt_count> anions{
        "chloride-anion", "chloride-anion", "chloride-anion",
        "bromide-anion", "bromide-anion", "bromide-anion",
    };
    if (table->component_count != 3 || table->component_ids == nullptr
        || table->component_charges == nullptr || table->component_ids[0] == nullptr
        || table->component_ids[1] == nullptr || table->component_ids[2] == nullptr
        || std::string(table->component_ids[0]) != "water"
        || std::string(table->component_ids[1]) != cations[salt_index]
        || std::string(table->component_ids[2]) != anions[salt_index]
        || table->component_charges[0] != 0 || table->component_charges[1] != 1
        || table->component_charges[2] != -1) {
        throw std::invalid_argument("provider component order does not match the aqueous salt");
    }
    if (payload.stage == "solvation_factor") {
        if (table->aqueous_miac_solvation_factor_result_size
                != sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1)
            || table->evaluate_aqueous_miac_solvation_factor == nullptr) {
            throw std::runtime_error("provider capability unavailable");
        }
    } else if (table->aqueous_miac_kij_result_size
            != sizeof(epcsaft_aqueous_miac_kij_result_v1)
        || table->evaluate_aqueous_miac_kij == nullptr) {
        throw std::runtime_error("provider capability unavailable");
    }
    return table;
}

std::array<const epcsaft_native_sdk_v1*, salt_count> parse_tables(
    PyObject* capsules, const Payload& payload
) {
    PyObject* sequence = PySequence_Fast(capsules, "provider capsules must be a sequence");
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 6) {
        Py_XDECREF(sequence);
        throw std::invalid_argument("provider capsules must contain six salt models");
    }
    std::array<const epcsaft_native_sdk_v1*, salt_count> result{};
    for (std::size_t index = 0; index < salt_count; ++index) {
        result[index] = checked_table(
            PySequence_Fast_GET_ITEM(sequence, static_cast<Py_ssize_t>(index)),
            payload,
            index
        );
    }
    Py_DECREF(sequence);
    return result;
}

Evaluation evaluate(
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
    const Payload& payload,
    const std::vector<double>& parameters
) {
    Evaluation evaluation{
        std::vector<double>(payload.observations.size()),
        std::vector<double>(payload.observations.size() * payload.parameter_count),
        {},
    };
    evaluation.rows.reserve(payload.observations.size());
    for (std::size_t row_index = 0; row_index < payload.observations.size(); ++row_index) {
        const Observation& observed = payload.observations[row_index];
        const auto* table = tables[observed.salt_index];
        double log_modeled = 0.0;
        double reference_molality = 0.0;
        double reference_convergence = 0.0;
        double derivative_convergence = 0.0;
        std::string fingerprint;
        std::vector<double> local_derivative(observed.columns.size());
        if (payload.stage == "solvation_factor") {
            epcsaft_aqueous_miac_solvation_factor_result_v1 result{};
            result.struct_size = sizeof(result);
            const int status = table->evaluate_aqueous_miac_solvation_factor(
                table->model_context,
                payload.expected_fingerprints[observed.salt_index].c_str(),
                payload.temperature,
                payload.pressure,
                observed.molality,
                parameters[0],
                &result
            );
            if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
                throw std::runtime_error(
                    std::string("provider solvent-factor evaluation failed: ")
                    + result.error
                );
            }
            log_modeled = result.log_mean_ionic_activity_coefficient_molality;
            local_derivative[0] = result.derivative;
            reference_molality = result.reference_molality_mol_per_kg;
            reference_convergence = result.reference_convergence_error;
            derivative_convergence = result.reference_derivative_convergence_error;
            fingerprint.assign(
                result.parameter_fingerprint,
                strnlen(result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE)
            );
        } else {
            std::array<double, 3> local_parameters{};
            for (std::size_t local = 0; local < 3; ++local) {
                local_parameters[local] = parameters[observed.columns[local]];
            }
            epcsaft_aqueous_miac_kij_result_v1 result{};
            result.struct_size = sizeof(result);
            const int status = table->evaluate_aqueous_miac_kij(
                table->model_context,
                payload.expected_fingerprints[observed.salt_index].c_str(),
                payload.temperature,
                payload.pressure,
                observed.molality,
                local_parameters.data(),
                local_parameters.size(),
                &result
            );
            if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
                throw std::runtime_error(
                    std::string("provider aqueous-kij evaluation failed: ")
                    + result.error
                );
            }
            log_modeled = result.log_mean_ionic_activity_coefficient_molality;
            std::copy(result.derivative, result.derivative + 3, local_derivative.begin());
            reference_molality = result.reference_molality_mol_per_kg;
            reference_convergence = result.reference_convergence_error;
            derivative_convergence = result.reference_derivative_convergence_error;
            fingerprint.assign(
                result.parameter_fingerprint,
                strnlen(result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE)
            );
        }
        if (fingerprint != payload.expected_fingerprints[observed.salt_index]
            || !std::isfinite(log_modeled) || !std::isfinite(reference_molality)
            || !std::isfinite(reference_convergence)
            || !std::isfinite(derivative_convergence)) {
            throw std::runtime_error("provider aqueous diagnostics failed");
        }
        const double modeled = std::exp(log_modeled);
        const double ratio = modeled / observed.observed;
        const double residual = 1.0 - ratio;
        evaluation.residuals[row_index] = residual;
        for (std::size_t local = 0; local < observed.columns.size(); ++local) {
            const std::size_t column = observed.columns[local];
            if (column >= payload.parameter_count) {
                throw std::invalid_argument(
                    "aqueous observation column is outside the parameter block"
                );
            }
            evaluation.jacobian[row_index * payload.parameter_count + column]
                = -ratio * local_derivative[local];
        }
        evaluation.rows.push_back(Row{
            observed.row_id,
            observed.salt,
            observed.molality,
            observed.observed,
            log_modeled,
            modeled,
            residual,
            std::move(local_derivative),
            reference_molality,
            reference_convergence,
            derivative_convergence,
            fingerprint,
        });
    }
    return evaluation;
}

class AqueousCost final : public ceres::CostFunction {
public:
    AqueousCost(
        const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
        const Payload& payload
    ) : tables_(tables), payload_(payload) {
        set_num_residuals(static_cast<int>(payload.observations.size()));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(payload.parameter_count));
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            const std::vector<double> parameters(
                values[0], values[0] + payload_.parameter_count
            );
            if (parameters != cached_parameters_) {
                cached_evaluation_ = evaluate(tables_, payload_, parameters);
                cached_parameters_ = parameters;
            }
            const Evaluation& result = *cached_evaluation_;
            std::copy(result.residuals.begin(), result.residuals.end(), residuals);
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(result.jacobian.begin(), result.jacobian.end(), jacobians[0]);
            }
            failure_reason_.clear();
            return true;
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
    }

    const std::string& failure_reason() const noexcept { return failure_reason_; }

private:
    std::array<const epcsaft_native_sdk_v1*, salt_count> tables_;
    const Payload& payload_;
    mutable std::string failure_reason_;
    mutable std::vector<double> cached_parameters_;
    mutable std::optional<Evaluation> cached_evaluation_;
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

SolveOutcome solve_one(
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
    const Payload& payload,
    const Start& start
) {
    SolveOutcome outcome{};
    outcome.name = start.name;
    outcome.parameters = start.values;
    AqueousCost cost(tables, payload);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, outcome.parameters.data());
    for (std::size_t index = 0; index < payload.parameter_count; ++index) {
        problem.SetParameterLowerBound(outcome.parameters.data(), index, payload.bounds[0]);
        problem.SetParameterUpperBound(outcome.parameters.data(), index, payload.bounds[1]);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    ceres::Solve(options, &problem, &outcome.summary);
    if (!cost.failure_reason().empty()) outcome.failure_reason = cost.failure_reason();
    outcome.evaluation = evaluate(tables, payload, outcome.parameters);
    outcome.complete_columns = std::all_of(
        outcome.evaluation.jacobian.begin(), outcome.evaluation.jacobian.end(),
        [](double value) { return std::isfinite(value); }
    );
    for (std::size_t column = 0; column < payload.parameter_count; ++column) {
        bool nonzero = false;
        for (std::size_t row = 0; row < payload.observations.size(); ++row) {
            nonzero = nonzero
                || outcome.evaluation.jacobian[
                    row * payload.parameter_count + column
                ] != 0.0;
        }
        outcome.complete_columns = outcome.complete_columns && nonzero;
    }
    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    > row_major(
        outcome.evaluation.jacobian.data(),
        static_cast<Eigen::Index>(payload.observations.size()),
        static_cast<Eigen::Index>(payload.parameter_count)
    );
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        Eigen::MatrixXd(row_major), Eigen::ComputeFullV
    );
    outcome.singular_values.resize(payload.parameter_count);
    for (std::size_t index = 0; index < payload.parameter_count; ++index) {
        outcome.singular_values[index] = svd.singularValues()[static_cast<Eigen::Index>(index)];
    }
    outcome.rank_threshold = outcome.singular_values.front()
        * static_cast<double>(std::max(payload.observations.size(), payload.parameter_count))
        * std::numeric_limits<double>::epsilon() * payload.rank_multiplier;
    outcome.rank = static_cast<int>(std::count_if(
        outcome.singular_values.begin(), outcome.singular_values.end(),
        [&](double value) { return value > outcome.rank_threshold; }
    ));
    if (outcome.singular_values.back() > 0.0) {
        outcome.condition_number = outcome.singular_values.front()
            / outcome.singular_values.back();
    }
    outcome.least_sensitive_direction.resize(payload.parameter_count);
    const auto least = svd.matrixV().col(static_cast<Eigen::Index>(payload.parameter_count - 1));
    for (std::size_t index = 0; index < payload.parameter_count; ++index) {
        outcome.least_sensitive_direction[index] = least[static_cast<Eigen::Index>(index)];
    }
    return outcome;
}

PyObject* doubles_to_tuple(const std::vector<double>& values) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyFloat_FromDouble(values[index]);
        if (value == nullptr) { Py_DECREF(tuple); return nullptr; }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* strings_to_tuple(const std::vector<std::string>& values) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyUnicode_FromString(values[index].c_str());
        if (value == nullptr) { Py_DECREF(tuple); return nullptr; }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* rows_to_tuple(const std::vector<Row>& rows) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(rows.size()));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const Row& row = rows[index];
        PyObject* derivative = doubles_to_tuple(row.local_derivative);
        PyObject* item = PyTuple_New(12);
        if (derivative == nullptr || item == nullptr) {
            Py_XDECREF(derivative); Py_XDECREF(item); Py_DECREF(tuple); return nullptr;
        }
        PyTuple_SET_ITEM(item, 0, PyUnicode_FromString(row.row_id.c_str()));
        PyTuple_SET_ITEM(item, 1, PyUnicode_FromString(row.salt.c_str()));
        PyTuple_SET_ITEM(item, 2, PyFloat_FromDouble(row.molality));
        PyTuple_SET_ITEM(item, 3, PyFloat_FromDouble(row.observed));
        PyTuple_SET_ITEM(item, 4, PyFloat_FromDouble(row.log_modeled));
        PyTuple_SET_ITEM(item, 5, PyFloat_FromDouble(row.modeled));
        PyTuple_SET_ITEM(item, 6, PyFloat_FromDouble(row.residual));
        PyTuple_SET_ITEM(item, 7, derivative);
        PyTuple_SET_ITEM(item, 8, PyFloat_FromDouble(row.reference_molality));
        PyTuple_SET_ITEM(item, 9, PyFloat_FromDouble(row.reference_convergence));
        PyTuple_SET_ITEM(item, 10, PyFloat_FromDouble(row.derivative_convergence));
        PyTuple_SET_ITEM(item, 11, PyUnicode_FromString(row.fingerprint.c_str()));
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), item);
    }
    return tuple;
}

PyObject* evaluation_to_tuple(const Evaluation& evaluation, const Payload& payload) {
    PyObject* residuals = doubles_to_tuple(evaluation.residuals);
    PyObject* jacobian = doubles_to_tuple(evaluation.jacobian);
    PyObject* rows = rows_to_tuple(evaluation.rows);
    PyObject* identity = strings_to_tuple(payload.identity);
    if (residuals == nullptr || jacobian == nullptr || rows == nullptr || identity == nullptr) {
        Py_XDECREF(residuals); Py_XDECREF(jacobian); Py_XDECREF(rows); Py_XDECREF(identity);
        return nullptr;
    }
    PyObject* result = PyTuple_New(4);
    if (result == nullptr) {
        Py_DECREF(residuals); Py_DECREF(jacobian); Py_DECREF(rows); Py_DECREF(identity);
        return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, residuals);
    PyTuple_SET_ITEM(result, 1, jacobian);
    PyTuple_SET_ITEM(result, 2, rows);
    PyTuple_SET_ITEM(result, 3, identity);
    return result;
}

PyObject* solution_to_tuple(const SolveOutcome& outcome) {
    PyObject* parameters = doubles_to_tuple(outcome.parameters);
    PyObject* residuals = doubles_to_tuple(outcome.evaluation.residuals);
    PyObject* jacobian = doubles_to_tuple(outcome.evaluation.jacobian);
    PyObject* rows = rows_to_tuple(outcome.evaluation.rows);
    PyObject* singular = doubles_to_tuple(outcome.singular_values);
    PyObject* least = doubles_to_tuple(outcome.least_sensitive_direction);
    if (parameters == nullptr || residuals == nullptr || jacobian == nullptr
        || rows == nullptr || singular == nullptr || least == nullptr) {
        Py_XDECREF(parameters); Py_XDECREF(residuals); Py_XDECREF(jacobian);
        Py_XDECREF(rows); Py_XDECREF(singular); Py_XDECREF(least); return nullptr;
    }
    PyObject* result = PyTuple_New(17);
    if (result == nullptr) {
        Py_DECREF(parameters); Py_DECREF(residuals); Py_DECREF(jacobian);
        Py_DECREF(rows); Py_DECREF(singular); Py_DECREF(least); return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, PyUnicode_FromString(outcome.name.c_str()));
    PyTuple_SET_ITEM(
        result,
        1,
        PyUnicode_FromString(
            termination_name(outcome.summary.termination_type).c_str()
        )
    );
    PyTuple_SET_ITEM(result, 2, Py_NewRef(outcome.summary.IsSolutionUsable() ? Py_True : Py_False));
    PyTuple_SET_ITEM(result, 3, PyFloat_FromDouble(outcome.summary.initial_cost));
    PyTuple_SET_ITEM(result, 4, PyFloat_FromDouble(outcome.summary.final_cost));
    PyTuple_SET_ITEM(result, 5, PyLong_FromSize_t(outcome.summary.iterations.size()));
    PyTuple_SET_ITEM(result, 6, parameters);
    PyTuple_SET_ITEM(result, 7, residuals);
    PyTuple_SET_ITEM(result, 8, jacobian);
    PyTuple_SET_ITEM(result, 9, rows);
    PyTuple_SET_ITEM(result, 10, singular);
    PyTuple_SET_ITEM(result, 11, PyFloat_FromDouble(outcome.rank_threshold));
    PyTuple_SET_ITEM(result, 12, PyLong_FromLong(outcome.rank));
    PyTuple_SET_ITEM(result, 13, PyFloat_FromDouble(outcome.condition_number));
    PyTuple_SET_ITEM(result, 14, least);
    PyTuple_SET_ITEM(result, 15, Py_NewRef(outcome.complete_columns ? Py_True : Py_False));
    PyTuple_SET_ITEM(result, 16, PyUnicode_FromString(outcome.failure_reason.c_str()));
    return result;
}

}  // namespace

PyObject* evaluate_figiel_aqueous_python(
    PyObject* capsules, PyObject* payload_object, PyObject* parameters_object
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules, payload);
        const auto parameters = doubles(
            parameters_object, payload.parameter_count, "aqueous parameters"
        );
        return evaluation_to_tuple(evaluate(tables, payload, parameters), payload);
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* solve_figiel_aqueous_python(PyObject* capsules, PyObject* payload_object) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules, payload);
        PyObject* solutions = PyTuple_New(static_cast<Py_ssize_t>(payload.starts.size()));
        if (solutions == nullptr) return nullptr;
        for (std::size_t index = 0; index < payload.starts.size(); ++index) {
            PyObject* item = solution_to_tuple(
                solve_one(tables, payload, payload.starts[index])
            );
            if (item == nullptr) { Py_DECREF(solutions); return nullptr; }
            PyTuple_SET_ITEM(solutions, static_cast<Py_ssize_t>(index), item);
        }
        PyObject* identity = strings_to_tuple(payload.identity);
        if (identity == nullptr) { Py_DECREF(solutions); return nullptr; }
        PyObject* result = PyTuple_New(2);
        if (result == nullptr) { Py_DECREF(solutions); Py_DECREF(identity); return nullptr; }
        PyTuple_SET_ITEM(result, 0, solutions);
        PyTuple_SET_ITEM(result, 1, identity);
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_regression
