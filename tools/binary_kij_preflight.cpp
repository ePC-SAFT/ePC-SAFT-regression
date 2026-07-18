#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <ceres/ceres.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <epcsaft/native_sdk_v1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t row_count = 17;
constexpr std::size_t residual_count = 68;
constexpr std::size_t variable_count = 35;
constexpr std::size_t coordinate_count = 4;
constexpr double gas_constant = 8.31446261815324;
constexpr double weight = 0.25;
constexpr double kij_scale = 0.01;
constexpr double liquid_reference = 6.5e-5;
constexpr std::array<double, 2> liquid_bounds{2.0e-5, 1.0e-4};
constexpr std::array<double, 2> vapor_bounds{1.0e-4, 1.0e-2};

struct PyDeleter final { void operator()(PyObject* value) const noexcept { Py_XDECREF(value); } };
using OwnedPy = std::unique_ptr<PyObject, PyDeleter>;

struct Row final {
    std::string id;
    double temperature{};
    double pressure{};
    double x{};
    double y{};
};

struct Phase final {
    double phi{};
    std::array<double, 4> gradient{};
    std::array<double, 16> hessian{};
    double pressure{};
    std::string fingerprint;
};

struct RowResult final {
    Row row;
    double liquid_volume{};
    double vapor_volume{};
    Phase liquid;
    Phase vapor;
    std::array<double, 4> raw{};
    double relative_separation{};
};

struct Evaluation final {
    std::array<double, residual_count> residuals{};
    std::array<double, residual_count * variable_count> jacobian{};
    std::array<RowResult, row_count> rows{};
    std::string fingerprint;
};

struct MatrixResult final {
    std::vector<double> singular;
    int rank{};
    double condition{std::numeric_limits<double>::infinity()};
    double threshold{};
};

struct SolveResult final {
    ceres::Solver::Summary summary;
    std::array<double, variable_count> variables{};
    Evaluation evaluation;
    MatrixResult full;
    double projected_singular{};
    int projected_rank{};
    bool complete{};
    double max_pressure{};
    double min_stability{std::numeric_limits<double>::infinity()};
    double min_separation{std::numeric_limits<double>::infinity()};
    std::size_t worst_row{};
    double max_mu{};
    std::string failure;
};

double finite_double(PyObject* value, const char* label) {
    const double result = PyFloat_AsDouble(value);
    if (PyErr_Occurred() != nullptr || !std::isfinite(result)) {
        throw std::invalid_argument(std::string(label) + " must be finite");
    }
    return result;
}

std::array<Row, row_count> parse_rows(PyObject* object) {
    OwnedPy sequence{PySequence_Fast(object, "rows must be a sequence")};
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence.get()) != 17) {
        throw std::invalid_argument("preflight requires all 17 source rows");
    }
    std::array<Row, row_count> result{};
    for (std::size_t index = 0; index < row_count; ++index) {
        OwnedPy item{PySequence_Fast(
            PySequence_Fast_GET_ITEM(sequence.get(), static_cast<Py_ssize_t>(index)),
            "row must be a sequence"
        )};
        if (item == nullptr || PySequence_Fast_GET_SIZE(item.get()) != 5) {
            throw std::invalid_argument("row must contain id, T, P, x, y");
        }
        PyObject** values = PySequence_Fast_ITEMS(item.get());
        const char* id = PyUnicode_AsUTF8(values[0]);
        if (id == nullptr) throw std::invalid_argument("row id must be text");
        result[index] = {
            id,
            finite_double(values[1], "temperature"),
            finite_double(values[2], "pressure"),
            finite_double(values[3], "x"),
            finite_double(values[4], "y"),
        };
        const std::string expected = "may2015-ch4-c2h6-"
            + std::string(index < 9 ? "00" : "0") + std::to_string(index + 1);
        if (result[index].id != expected || result[index].temperature <= 0.0
            || result[index].pressure <= 0.0 || result[index].x <= 0.0
            || result[index].x >= 1.0 || result[index].y <= 0.0
            || result[index].y >= 1.0) {
            throw std::invalid_argument("source row identity or domain changed");
        }
    }
    return result;
}

const epcsaft_native_sdk_v1& table_from_capsule(PyObject* capsule) {
    if (!PyCapsule_CheckExact(capsule)) throw std::invalid_argument("provider must be a capsule");
    void* pointer = PyCapsule_GetPointer(capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME);
    if (pointer == nullptr) throw std::runtime_error("provider capsule pointer is unavailable");
    struct Prefix final { std::uint32_t abi; std::size_t size; } prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    constexpr std::size_t required = offsetof(epcsaft_native_sdk_v1, evaluate_mixture_phase_kij)
        + sizeof(epcsaft_evaluate_mixture_phase_kij_v1);
    if (prefix.abi != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION || prefix.size < required) {
        throw std::runtime_error("provider lacks the active-kij ABI tail");
    }
    const auto& table = *static_cast<const epcsaft_native_sdk_v1*>(pointer);
    if (table.component_count != 2
        || table.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)
        || table.model_context == nullptr || table.evaluate_mixture_phase_kij == nullptr) {
        throw std::runtime_error("provider active-kij contract is unavailable");
    }
    return table;
}

Phase phase(
    const epcsaft_native_sdk_v1& table,
    double temperature,
    std::array<double, 2> amounts,
    double volume,
    double kij
) {
    Phase result{};
    epcsaft_mixture_phase_block_result_v1 native{};
    native.struct_size = sizeof(native);
    native.coordinate_count = 4;
    native.gradient_capacity = 4;
    native.hessian_capacity = 16;
    native.gradient = result.gradient.data();
    native.hessian = result.hessian.data();
    const int status = table.evaluate_mixture_phase_kij(
        table.model_context, temperature, amounts.data(), 2, volume, kij, &native
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || native.status != status) {
        throw std::runtime_error(std::string("provider evaluation failed: ") + native.error);
    }
    result.phi = native.helmholtz_over_rt_reference_amount;
    result.pressure = native.pressure_pa;
    result.fingerprint.assign(
        native.parameter_fingerprint,
        strnlen(native.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE)
    );
    const bool finite = std::isfinite(result.phi) && std::isfinite(result.pressure)
        && std::all_of(result.gradient.begin(), result.gradient.end(), [](double value) {
            return std::isfinite(value);
        })
        && std::all_of(result.hessian.begin(), result.hessian.end(), [](double value) {
            return std::isfinite(value);
        });
    if (!finite || result.fingerprint.empty()) throw std::runtime_error("nonfinite provider result");
    return result;
}

double h(const Phase& value, std::size_t row, std::size_t column) {
    return value.hessian[row * coordinate_count + column];
}

Evaluation evaluate(
    const epcsaft_native_sdk_v1& table,
    const std::array<Row, row_count>& rows,
    const std::array<double, variable_count>& variables
) {
    Evaluation result{};
    const double kij = kij_scale * variables[0];
    for (std::size_t index = 0; index < row_count; ++index) {
        const Row& source = rows[index];
        const double liquid_volume = liquid_reference * std::exp(variables[1 + 2 * index]);
        const double vapor_reference = gas_constant * source.temperature / source.pressure;
        const double vapor_volume = vapor_reference * std::exp(variables[2 + 2 * index]);
        if (liquid_volume < liquid_bounds[0] || liquid_volume > liquid_bounds[1]
            || vapor_volume < vapor_bounds[0] || vapor_volume > vapor_bounds[1]
            || liquid_volume >= vapor_volume) {
            throw std::runtime_error("lifted volume bounds or ordering failed");
        }
        const Phase liquid = phase(
            table, source.temperature, {source.x, 1.0 - source.x}, liquid_volume, kij
        );
        const Phase vapor = phase(
            table, source.temperature, {source.y, 1.0 - source.y}, vapor_volume, kij
        );
        if (liquid.fingerprint != vapor.fingerprint
            || (!result.fingerprint.empty() && result.fingerprint != liquid.fingerprint)) {
            throw std::runtime_error("provider fingerprint changed within the problem");
        }
        result.fingerprint = liquid.fingerprint;
        RowResult& row = result.rows[index];
        row.row = source;
        row.liquid_volume = liquid_volume;
        row.vapor_volume = vapor_volume;
        row.liquid = liquid;
        row.vapor = vapor;
        row.relative_separation = (vapor_volume - liquid_volume) / vapor_volume;
        row.raw = {
            liquid.pressure - source.pressure,
            vapor.pressure - source.pressure,
            liquid.gradient[0] - vapor.gradient[0],
            liquid.gradient[1] - vapor.gradient[1],
        };
        const std::size_t residual = 4 * index;
        result.residuals[residual] = weight * row.raw[0] / source.pressure;
        result.residuals[residual + 1] = weight * row.raw[1] / source.pressure;
        result.residuals[residual + 2] = weight * row.raw[2];
        result.residuals[residual + 3] = weight * row.raw[3];
        const std::size_t liquid_column = 1 + 2 * index;
        const std::size_t vapor_column = liquid_column + 1;
        const auto set = [&result](std::size_t r, std::size_t c, double value) {
            result.jacobian[r * variable_count + c] = value;
        };
        set(residual, 0, weight / source.pressure
            * (-gas_constant * source.temperature * h(liquid, 2, 3)) * kij_scale);
        set(residual, liquid_column, weight / source.pressure
            * (-gas_constant * source.temperature * h(liquid, 2, 2)) * liquid_volume);
        set(residual + 1, 0, weight / source.pressure
            * (-gas_constant * source.temperature * h(vapor, 2, 3)) * kij_scale);
        set(residual + 1, vapor_column, weight / source.pressure
            * (-gas_constant * source.temperature * h(vapor, 2, 2)) * vapor_volume);
        for (std::size_t component = 0; component < 2; ++component) {
            set(residual + 2 + component, 0, weight
                * (h(liquid, component, 3) - h(vapor, component, 3)) * kij_scale);
            set(residual + 2 + component, liquid_column,
                weight * h(liquid, component, 2) * liquid_volume);
            set(residual + 2 + component, vapor_column,
                -weight * h(vapor, component, 2) * vapor_volume);
        }
    }
    return result;
}

class Cost final : public ceres::CostFunction {
public:
    Cost(const epcsaft_native_sdk_v1& table, const std::array<Row, row_count>& rows, std::string* error)
        : table_(table), rows_(rows), error_(error) {
        set_num_residuals(68);
        mutable_parameter_block_sizes()->push_back(35);
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        try {
            std::array<double, variable_count> variables{};
            std::copy_n(parameters[0], variable_count, variables.begin());
            const Evaluation result = evaluate(table_, rows_, variables);
            std::copy(result.residuals.begin(), result.residuals.end(), residuals);
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(result.jacobian.begin(), result.jacobian.end(), jacobians[0]);
            }
            return true;
        } catch (const std::exception& exception) {
            *error_ = exception.what();
            return false;
        }
    }

private:
    const epcsaft_native_sdk_v1& table_;
    const std::array<Row, row_count>& rows_;
    std::string* error_;
};

std::string termination(ceres::TerminationType value) {
    switch (value) {
        case ceres::CONVERGENCE: return "CONVERGENCE";
        case ceres::NO_CONVERGENCE: return "NO_CONVERGENCE";
        case ceres::FAILURE: return "FAILURE";
        case ceres::USER_SUCCESS: return "USER_SUCCESS";
        case ceres::USER_FAILURE: return "USER_FAILURE";
    }
    return "UNKNOWN";
}

MatrixResult matrix_result(const Eigen::MatrixXd& matrix) {
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    MatrixResult result{};
    for (Eigen::Index index = 0; index < svd.singularValues().size(); ++index) {
        result.singular.push_back(svd.singularValues()(index));
    }
    const double maximum = result.singular.empty() ? 0.0 : result.singular.front();
    result.threshold = maximum * static_cast<double>(std::max(matrix.rows(), matrix.cols()))
        * std::numeric_limits<double>::epsilon();
    double minimum = std::numeric_limits<double>::infinity();
    for (double value : result.singular) {
        if (value > result.threshold) {
            ++result.rank;
            minimum = std::min(minimum, value);
        }
    }
    if (result.rank > 0) result.condition = maximum / minimum;
    return result;
}

SolveResult solve(
    const epcsaft_native_sdk_v1& table,
    const std::array<Row, row_count>& rows,
    double kij_start,
    double liquid_multiplier,
    double vapor_multiplier
) {
    SolveResult result{};
    result.variables[0] = kij_start / kij_scale;
    for (std::size_t index = 0; index < row_count; ++index) {
        result.variables[1 + 2 * index] = std::log(liquid_multiplier);
        result.variables[2 + 2 * index] = std::log(vapor_multiplier);
    }
    ceres::Problem problem;
    problem.AddResidualBlock(new Cost(table, rows, &result.failure), nullptr, result.variables.data());
    problem.SetParameterLowerBound(result.variables.data(), 0, -15.0);
    problem.SetParameterUpperBound(result.variables.data(), 0, 10.0);
    for (std::size_t index = 0; index < row_count; ++index) {
        const double vapor_reference = gas_constant * rows[index].temperature / rows[index].pressure;
        problem.SetParameterLowerBound(result.variables.data(), static_cast<int>(1 + 2 * index),
            std::log(liquid_bounds[0] / liquid_reference));
        problem.SetParameterUpperBound(result.variables.data(), static_cast<int>(1 + 2 * index),
            std::log(liquid_bounds[1] / liquid_reference));
        problem.SetParameterLowerBound(result.variables.data(), static_cast<int>(2 + 2 * index),
            std::log(vapor_bounds[0] / vapor_reference));
        problem.SetParameterUpperBound(result.variables.data(), static_cast<int>(2 + 2 * index),
            std::log(vapor_bounds[1] / vapor_reference));
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 500;
    options.function_tolerance = 1.0e-10;
    options.gradient_tolerance = 1.0e-10;
    options.parameter_tolerance = 1.0e-10;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    ceres::Solve(options, &problem, &result.summary);
    result.evaluation = evaluate(table, rows, result.variables);
    result.complete = std::all_of(
        result.evaluation.jacobian.begin(), result.evaluation.jacobian.end(),
        [](double value) { return std::isfinite(value); }
    );
    Eigen::MatrixXd full(68, 35);
    for (std::size_t row = 0; row < residual_count; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            full(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                result.evaluation.jacobian[row * variable_count + column];
        }
    }
    result.full = matrix_result(full);
    const Eigen::MatrixXd volumes = full.rightCols(34);
    const Eigen::JacobiSVD<Eigen::MatrixXd> volume_svd(
        volumes, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    const double volume_threshold = volume_svd.singularValues()(0) * 68.0
        * std::numeric_limits<double>::epsilon();
    Eigen::Index volume_rank = 0;
    while (volume_rank < volume_svd.singularValues().size()
        && volume_svd.singularValues()(volume_rank) > volume_threshold) {
        ++volume_rank;
    }
    Eigen::VectorXd projected = full.col(0);
    if (volume_rank > 0) {
        const Eigen::MatrixXd basis = volume_svd.matrixU().leftCols(volume_rank);
        projected -= basis * (basis.transpose() * full.col(0));
    }
    result.projected_singular = projected.norm();
    result.projected_rank = result.projected_singular > result.full.threshold ? 1 : 0;
    for (std::size_t index = 0; index < row_count; ++index) {
        const RowResult& row = result.evaluation.rows[index];
        const double row_pressure = std::max(
            std::abs(row.raw[0] / row.row.pressure), std::abs(row.raw[1] / row.row.pressure)
        );
        if (row_pressure > result.max_pressure) {
            result.max_pressure = row_pressure;
            result.worst_row = index;
        }
        result.max_mu = std::max({result.max_mu, std::abs(row.raw[2]), std::abs(row.raw[3])});
        result.min_stability = std::min({result.min_stability, h(row.liquid, 2, 2), h(row.vapor, 2, 2)});
        result.min_separation = std::min(result.min_separation, row.relative_separation);
    }
    if (result.summary.termination_type == ceres::CONVERGENCE && result.summary.IsSolutionUsable()) {
        result.failure.clear();
    } else if (result.failure.empty()) {
        result.failure = "Ceres did not return a usable converged solution";
    }
    return result;
}

double scaled_error(double observed, double expected, double atol, double rtol) {
    return std::abs(observed - expected) / (atol + rtol * std::abs(expected));
}

struct CheckResult final {
    double value{};
    double gradient{};
    double symmetry{};
    double pressure{};
    double residual_abs{};
    double residual_scaled{};
    Phase representative_liquid;
    Phase representative_vapor;
    double liquid_volume{};
    double vapor_volume{};
};

void phase_check(
    const epcsaft_native_sdk_v1& table,
    double temperature,
    std::array<double, 2> amounts,
    double volume,
    CheckResult& result
) {
    constexpr double step = 1.0e-5;
    const std::array<double, 4> base{amounts[0], amounts[1], volume, 0.0};
    const std::array<double, 4> direction{0.07, -0.05, 0.12 * volume, 0.015};
    const auto shifted = [&base, &direction](double sign) {
        std::array<double, 4> point{};
        for (std::size_t index = 0; index < 4; ++index) {
            point[index] = base[index] + sign * step * direction[index];
        }
        return point;
    };
    const auto plus_point = shifted(1.0);
    const auto minus_point = shifted(-1.0);
    const Phase center = phase(table, temperature, amounts, volume, 0.0);
    const Phase plus = phase(
        table, temperature, {plus_point[0], plus_point[1]}, plus_point[2], plus_point[3]
    );
    const Phase minus = phase(
        table, temperature, {minus_point[0], minus_point[1]}, minus_point[2], minus_point[3]
    );
    double value_exact = 0.0;
    for (std::size_t index = 0; index < 4; ++index) value_exact += center.gradient[index] * direction[index];
    const double value_fd = (plus.phi - minus.phi) / (2.0 * step);
    result.value = std::max(result.value, scaled_error(value_fd, value_exact, 5.0e-7, 5.0e-5));
    for (std::size_t row = 0; row < 4; ++row) {
        double exact = 0.0;
        for (std::size_t column = 0; column < 4; ++column) exact += h(center, row, column) * direction[column];
        const double finite_difference = (plus.gradient[row] - minus.gradient[row]) / (2.0 * step);
        result.gradient = std::max(
            result.gradient, scaled_error(finite_difference, exact, 5.0e-7, 5.0e-5)
        );
        for (std::size_t column = 0; column < row; ++column) {
            result.symmetry = std::max(
                result.symmetry,
                scaled_error(h(center, row, column), h(center, column, row), 1.0e-12, 5.0e-13)
            );
        }
    }
    result.pressure = std::max(
        result.pressure,
        scaled_error(center.pressure, -gas_constant * temperature * center.gradient[2], 1.0e-6, 5.0e-13)
    );
}

CheckResult checks(
    const epcsaft_native_sdk_v1& table,
    const std::array<Row, row_count>& rows
) {
    CheckResult result{};
    const Row& source = rows[0];
    result.liquid_volume = liquid_reference;
    result.vapor_volume = gas_constant * source.temperature / source.pressure;
    result.representative_liquid = phase(
        table, source.temperature, {source.x, 1.0 - source.x}, result.liquid_volume, 0.0
    );
    result.representative_vapor = phase(
        table, source.temperature, {source.y, 1.0 - source.y}, result.vapor_volume, 0.0
    );
    phase_check(table, source.temperature, {source.x, 1.0 - source.x}, result.liquid_volume, result);
    phase_check(table, source.temperature, {source.y, 1.0 - source.y}, result.vapor_volume, result);

    std::array<double, variable_count> variables{};
    std::array<double, variable_count> direction{};
    direction[0] = 0.2;
    for (std::size_t index = 0; index < row_count; ++index) {
        direction[1 + 2 * index] = 0.012 - 0.0007 * static_cast<double>(index);
        direction[2 + 2 * index] = -0.009 + 0.0005 * static_cast<double>(index);
    }
    constexpr double step = 1.0e-6;
    std::array<double, variable_count> plus_variables{};
    std::array<double, variable_count> minus_variables{};
    for (std::size_t index = 0; index < variable_count; ++index) {
        plus_variables[index] = variables[index] + step * direction[index];
        minus_variables[index] = variables[index] - step * direction[index];
    }
    const Evaluation center = evaluate(table, rows, variables);
    const Evaluation plus = evaluate(table, rows, plus_variables);
    const Evaluation minus = evaluate(table, rows, minus_variables);
    for (std::size_t row = 0; row < residual_count; ++row) {
        double exact = 0.0;
        for (std::size_t column = 0; column < variable_count; ++column) {
            exact += center.jacobian[row * variable_count + column] * direction[column];
        }
        const double finite_difference = (plus.residuals[row] - minus.residuals[row]) / (2.0 * step);
        result.residual_abs = std::max(result.residual_abs, std::abs(finite_difference - exact));
        result.residual_scaled = std::max(
            result.residual_scaled, scaled_error(finite_difference, exact, 2.0e-7, 2.0e-6)
        );
    }
    return result;
}

PyObject* doubles(const double* values, std::size_t size) {
    PyObject* result = PyTuple_New(static_cast<Py_ssize_t>(size));
    if (result == nullptr) return nullptr;
    for (std::size_t index = 0; index < size; ++index) {
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), PyFloat_FromDouble(values[index]));
    }
    return result;
}

PyObject* phase_tuple(const Phase& value, double volume, double n1) {
    PyObject* gradient = doubles(value.gradient.data(), value.gradient.size());
    if (gradient == nullptr) return nullptr;
    return Py_BuildValue(
        "(ddddNs)", volume, n1, 1.0 - n1, value.pressure, gradient, value.fingerprint.c_str()
    );
}

PyObject* solve_tuple(const SolveResult& value) {
    const RowResult& worst = value.evaluation.rows[value.worst_row];
    PyObject* singular = doubles(value.full.singular.data(), value.full.singular.size());
    PyObject* worst_row = Py_BuildValue(
        "(sddddd)",
        worst.row.id.c_str(),
        worst.raw[0] / worst.row.pressure,
        worst.raw[1] / worst.row.pressure,
        worst.raw[2],
        worst.raw[3],
        value.max_mu
    );
    if (singular == nullptr || worst_row == nullptr) {
        Py_XDECREF(singular);
        Py_XDECREF(worst_row);
        return nullptr;
    }
    return Py_BuildValue(
        "(sOddKdNidddidddONs)",
        termination(value.summary.termination_type).c_str(),
        value.summary.IsSolutionUsable() ? Py_True : Py_False,
        value.summary.initial_cost,
        value.summary.final_cost,
        value.summary.iterations.size(),
        kij_scale * value.variables[0],
        singular,
        value.full.rank,
        value.full.condition,
        value.full.threshold,
        value.projected_singular,
        value.projected_rank,
        value.max_pressure,
        value.min_stability,
        value.min_separation,
        value.complete ? Py_True : Py_False,
        worst_row,
        value.failure.c_str()
    );
}

PyObject* py_run(PyObject*, PyObject* arguments) {
    PyObject* capsule = nullptr;
    PyObject* rows_object = nullptr;
    if (!PyArg_ParseTuple(arguments, "OO:run", &capsule, &rows_object)) return nullptr;
    try {
        const auto rows = parse_rows(rows_object);
        const auto& table = table_from_capsule(capsule);
        const CheckResult check = checks(table, rows);
        const std::array<SolveResult, 3> solves{
            solve(table, rows, 0.0, 1.0, 1.0),
            solve(table, rows, -0.05, 1.01, 0.98),
            solve(table, rows, 0.05, 1.01, 0.98),
        };
        PyObject* callback = Py_BuildValue(
            "(dddd)", check.value, check.gradient, check.symmetry, check.pressure
        );
        PyObject* residual = Py_BuildValue(
            "(ddiii)", check.residual_abs, check.residual_scaled, 68, 35, 2380
        );
        PyObject* liquid = phase_tuple(check.representative_liquid, check.liquid_volume, rows[0].x);
        PyObject* vapor = phase_tuple(check.representative_vapor, check.vapor_volume, rows[0].y);
        PyObject* solve_values = PyTuple_New(3);
        if (callback == nullptr || residual == nullptr || liquid == nullptr || vapor == nullptr
            || solve_values == nullptr) {
            Py_XDECREF(callback);
            Py_XDECREF(residual);
            Py_XDECREF(liquid);
            Py_XDECREF(vapor);
            Py_XDECREF(solve_values);
            return nullptr;
        }
        for (std::size_t index = 0; index < solves.size(); ++index) {
            PyObject* item = solve_tuple(solves[index]);
            if (item == nullptr) {
                Py_DECREF(callback);
                Py_DECREF(residual);
                Py_DECREF(liquid);
                Py_DECREF(vapor);
                Py_DECREF(solve_values);
                return nullptr;
            }
            PyTuple_SET_ITEM(solve_values, static_cast<Py_ssize_t>(index), item);
        }
        PyObject* result = PyTuple_New(5);
        PyTuple_SET_ITEM(result, 0, callback);
        PyTuple_SET_ITEM(result, 1, residual);
        PyTuple_SET_ITEM(result, 2, liquid);
        PyTuple_SET_ITEM(result, 3, vapor);
        PyTuple_SET_ITEM(result, 4, solve_values);
        return result;
    } catch (const std::exception& exception) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, exception.what());
        return nullptr;
    }
}

PyMethodDef methods[] = {
    {"run", py_run, METH_VARARGS, "Run the frozen installed-artifact binary-kij preflight."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_binary_kij_preflight",
    "Non-installed binary-kij preflight evidence target.",
    -1,
    methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

}  // namespace

PyMODINIT_FUNC PyInit__binary_kij_preflight() { return PyModule_Create(&module); }
