#include "general_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_regression {
namespace {

constexpr double gas_constant = 8.31446261815324;

struct Row final {
    std::string row_id;
    std::string partition;
    double temperature;
    double pressure;
    double liquid_first;
    double vapor_first;
    double pressure_scale;
    std::array<double, 2> chemical_potential_scales;
    double liquid_volume_origin;
    double liquid_volume_start;
    std::array<double, 2> liquid_volume_bounds;
    double vapor_volume_origin;
    double vapor_volume_start;
    std::array<double, 2> vapor_volume_bounds;
};

struct Payload final {
    std::string capability_id;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::array<std::string, 2> component_ids;
    double parameter_origin;
    double parameter_scale;
    double parameter_lower_bound;
    double parameter_upper_bound;
    std::vector<double> starts;
    double maximum_condition_number;
    int maximum_iterations;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    double confirmation_parameter_delta;
    double confirmation_cost_delta;
    std::vector<Row> training_rows;
    std::vector<Row> reporting_rows;
};

struct Phase final {
    double pressure;
    std::array<double, 4> gradient;
    std::array<double, 16> hessian;
};

struct Evaluation final {
    std::vector<double> residuals;
    std::vector<double> jacobian;
};

struct MatrixDiagnostics final {
    std::vector<double> singular_values;
    int rank{0};
    double condition_number{std::numeric_limits<double>::infinity()};
};

struct SolveOutcome final {
    ceres::Solver::Summary summary;
    std::vector<double> variables;
    Evaluation evaluation;
    MatrixDiagnostics full_jacobian;
    MatrixDiagnostics projected_parameter_jacobian;
    std::string failure_reason;
};

struct RowOutcome final {
    Row row;
    double liquid_volume;
    double vapor_volume;
    std::array<double, 4> residuals;
    bool usable;
    std::string failure_reason;
};

class OwnedPyObject final {
public:
    explicit OwnedPyObject(PyObject* object) : object_(object) {}
    ~OwnedPyObject() { Py_XDECREF(object_); }
    OwnedPyObject(const OwnedPyObject&) = delete;
    OwnedPyObject& operator=(const OwnedPyObject&) = delete;
    PyObject* get() const noexcept { return object_; }
    explicit operator bool() const noexcept { return object_ != nullptr; }

private:
    PyObject* object_;
};

double number(PyObject* object, const char* name) {
    const double value = PyFloat_AsDouble(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(name) + " must be numeric");
    }
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
    return value;
}

std::string text(PyObject* object, const char* name) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(name) + " must be text");
    }
    const char* value = PyUnicode_AsUTF8(object);
    if (value == nullptr || value[0] == '\0') {
        PyErr_Clear();
        throw std::invalid_argument(std::string(name) + " must be nonempty");
    }
    return value;
}

std::vector<double> doubles(PyObject* object, const char* name) {
    OwnedPyObject sequence{PySequence_Fast(object, name)};
    if (!sequence) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(name) + " must be a sequence");
    }
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence.get());
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(count));
    for (Py_ssize_t index = 0; index < count; ++index) {
        values.push_back(number(
            PySequence_Fast_GET_ITEM(sequence.get(), index), name
        ));
    }
    return values;
}

Row parse_row(PyObject* object) {
    OwnedPyObject sequence{
        PySequence_Fast(object, "observation payload must be a sequence")
    };
    if (!sequence || PySequence_Fast_GET_SIZE(sequence.get()) != 17) {
        PyErr_Clear();
        throw std::invalid_argument(
            "observation payload must contain exactly 17 fields"
        );
    }
    auto item = [&](Py_ssize_t index) {
        return PySequence_Fast_GET_ITEM(sequence.get(), index);
    };
    Row row{};
    row.row_id = text(item(0), "row id");
    row.partition = text(item(1), "partition");
    row.temperature = number(item(2), "temperature");
    row.pressure = number(item(3), "pressure");
    row.liquid_first = number(item(4), "liquid composition");
    row.vapor_first = number(item(5), "vapor composition");
    row.pressure_scale = number(item(6), "pressure scale");
    row.chemical_potential_scales = {
        number(item(7), "first chemical-potential scale"),
        number(item(8), "second chemical-potential scale"),
    };
    row.liquid_volume_origin = number(item(9), "liquid volume origin");
    row.liquid_volume_start = number(item(10), "liquid volume start");
    row.liquid_volume_bounds = {
        number(item(11), "liquid volume lower bound"),
        number(item(12), "liquid volume upper bound"),
    };
    row.vapor_volume_origin = number(item(13), "vapor volume origin");
    row.vapor_volume_start = number(item(14), "vapor volume start");
    row.vapor_volume_bounds = {
        number(item(15), "vapor volume lower bound"),
        number(item(16), "vapor volume upper bound"),
    };
    return row;
}

std::vector<Row> parse_rows(PyObject* object, const char* name) {
    OwnedPyObject sequence{PySequence_Fast(object, name)};
    if (!sequence) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(name) + " must be a sequence");
    }
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence.get());
    std::vector<Row> rows;
    rows.reserve(static_cast<std::size_t>(count));
    for (Py_ssize_t index = 0; index < count; ++index) {
        rows.push_back(parse_row(
            PySequence_Fast_GET_ITEM(sequence.get(), index)
        ));
    }
    return rows;
}

Payload parse_payload(PyObject* object) {
    OwnedPyObject sequence{
        PySequence_Fast(object, "general regression payload must be a sequence")
    };
    if (!sequence || PySequence_Fast_GET_SIZE(sequence.get()) != 18) {
        PyErr_Clear();
        throw std::invalid_argument(
            "general regression payload must contain exactly 18 fields"
        );
    }
    auto item = [&](Py_ssize_t index) {
        return PySequence_Fast_GET_ITEM(sequence.get(), index);
    };
    OwnedPyObject components{
        PySequence_Fast(item(3), "component ids must be a sequence")
    };
    if (!components || PySequence_Fast_GET_SIZE(components.get()) != 2) {
        PyErr_Clear();
        throw std::invalid_argument("component ids must contain exactly two values");
    }
    Payload payload{};
    payload.capability_id = text(item(0), "capability id");
    payload.parameter_fingerprint = text(item(1), "parameter fingerprint");
    payload.topology_fingerprint = text(item(2), "topology fingerprint");
    payload.component_ids = {
        text(PySequence_Fast_GET_ITEM(components.get(), 0), "first component id"),
        text(PySequence_Fast_GET_ITEM(components.get(), 1), "second component id"),
    };
    payload.parameter_origin = number(item(4), "parameter origin");
    payload.parameter_scale = number(item(5), "parameter scale");
    payload.parameter_lower_bound = number(item(6), "parameter lower bound");
    payload.parameter_upper_bound = number(item(7), "parameter upper bound");
    payload.starts = doubles(item(8), "parameter starts");
    payload.maximum_condition_number = number(
        item(9), "maximum condition number"
    );
    const long maximum_iterations = PyLong_AsLong(item(10));
    if (PyErr_Occurred() != nullptr || maximum_iterations <= 0
        || maximum_iterations > std::numeric_limits<int>::max()) {
        PyErr_Clear();
        throw std::invalid_argument("maximum iterations must be a positive integer");
    }
    payload.maximum_iterations = static_cast<int>(maximum_iterations);
    payload.function_tolerance = number(item(11), "function tolerance");
    payload.gradient_tolerance = number(item(12), "gradient tolerance");
    payload.parameter_tolerance = number(item(13), "parameter tolerance");
    payload.confirmation_parameter_delta = number(
        item(14), "confirmation parameter delta"
    );
    payload.confirmation_cost_delta = number(
        item(15), "confirmation cost delta"
    );
    payload.training_rows = parse_rows(item(16), "training rows");
    payload.reporting_rows = parse_rows(item(17), "reporting rows");
    if (payload.training_rows.empty()) {
        throw std::invalid_argument("at least one training row is required");
    }
    return payload;
}

const epcsaft_native_sdk_v1* capability_table(PyObject* capsule) {
    if (!PyCapsule_CheckExact(capsule)) {
        throw std::invalid_argument("provider transport must be an exact CPython capsule");
    }
    void* pointer = PyCapsule_GetPointer(
        capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME
    );
    if (pointer == nullptr) {
        throw std::invalid_argument("provider transport capsule is invalid");
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    constexpr std::size_t required_size =
        offsetof(epcsaft_native_sdk_v1, capabilities)
        + sizeof(table->capabilities);
    if (table->abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || table->table_size < required_size) {
        throw std::runtime_error(
            "provider native SDK lacks the capability descriptor tail"
        );
    }
    if (table->capability_count > 0 && table->capabilities == nullptr) {
        throw std::runtime_error("provider capability descriptor pointer is null");
    }
    return table;
}

bool bounded_field_equal(
    const std::string& expected,
    const char field[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE]
) {
    const std::size_t length = strnlen(
        field, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    return length == expected.size()
        && std::memcmp(field, expected.data(), length) == 0;
}

const char* capability_id(std::uint32_t value);

void validate_descriptor(
    const epcsaft_native_capability_descriptor_v1& descriptor
);

const epcsaft_native_capability_descriptor_v1& checked_descriptor(
    const epcsaft_native_sdk_v1& table, const Payload& payload
) {
    if (table.capability_count != 1 || table.capabilities == nullptr) {
        throw std::runtime_error(
            "provider does not advertise exactly one supported capability"
        );
    }
    const auto& descriptor = table.capabilities[0];
    validate_descriptor(descriptor);
    if (payload.capability_id != capability_id(descriptor.capability)
        || !bounded_field_equal(
            payload.parameter_fingerprint, descriptor.parameter_fingerprint
        )
        || !bounded_field_equal(
            payload.topology_fingerprint, descriptor.topology_fingerprint
        )
        || descriptor.component_count != 2
        || descriptor.component_ids == nullptr
        || payload.component_ids[0] != descriptor.component_ids[0]
        || payload.component_ids[1] != descriptor.component_ids[1]
        || table.evaluate_mixture_phase_kij == nullptr
        || table.mixture_result_size
            != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw std::runtime_error(
            "regression problem does not match the installed Provider capability"
        );
    }
    return descriptor;
}

const char* capability_id(std::uint32_t value) {
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1) {
        return "neutral_binary_phase_kij_v1";
    }
    throw std::runtime_error("provider advertised an unknown capability");
}

const char* parameter_family(std::uint32_t value) {
    if (value
        == EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1) {
        return "k_ij";
    }
    throw std::runtime_error("provider advertised an unknown parameter family");
}

const char* coordinate_kind(std::uint32_t value) {
    switch (value) {
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1:
            return "amount";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1:
            return "volume";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1:
            return "k_ij";
        default:
            throw std::runtime_error(
                "provider advertised an unknown capability coordinate"
            );
    }
}

void validate_descriptor(
    const epcsaft_native_capability_descriptor_v1& descriptor
) {
    if (descriptor.struct_size
            != sizeof(epcsaft_native_capability_descriptor_v1)
        || descriptor.schema_version
            != EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1
        || descriptor.capability
            != EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1
        || descriptor.parameter_family
            != EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1
        || descriptor.parameter_identity
            != EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1
        || descriptor.observation_contract
            != EPCSAFT_NATIVE_OBSERVATION_FIXED_COMPOSITION_HELMHOLTZ_PHASE_V1
        || descriptor.model_domain
            != EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_BINARY_V1
        || descriptor.tensor_layout
            != EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1
        || descriptor.derivative_order != 2
        || descriptor.maturity
            != EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1
        || descriptor.authority_effect
            != EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1
        || descriptor.unsupported_status
            != EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1
        || descriptor.domain_status != EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1
        || descriptor.state_coordinate_count != 3
        || descriptor.active_parameter_count != 1
        || descriptor.coordinate_count != 4
        || descriptor.component_count != 2
        || descriptor.coordinates == nullptr
        || descriptor.component_ids == nullptr
        || descriptor.component_ids[0] == nullptr
        || descriptor.component_ids[1] == nullptr
        || !std::isfinite(descriptor.temperature_min_k)
        || !std::isfinite(descriptor.temperature_max_k)
        || descriptor.temperature_min_k >= descriptor.temperature_max_k
        || descriptor.helmholtz_basis_id
            != std::string(EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID)) {
        throw std::runtime_error(
            "provider capability descriptor does not match the supported v1 contract"
        );
    }
    const std::array<std::uint32_t, 4> kinds = {
        EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
        EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
        EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1,
        EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1,
    };
    const std::array<int, 4> components = {0, 1, -1, -1};
    const std::array<int, 4> pair_a = {-1, -1, -1, 0};
    const std::array<int, 4> pair_b = {-1, -1, -1, 1};
    const std::array<const char*, 4> units = {
        "mol", "mol", "m3", "dimensionless"
    };
    for (std::size_t index = 0; index < 4; ++index) {
        const auto& coordinate = descriptor.coordinates[index];
        if (coordinate.struct_size
                != sizeof(epcsaft_native_capability_coordinate_v1)
            || coordinate.kind != kinds[index]
            || coordinate.component_index != components[index]
            || coordinate.pair_component_index_a != pair_a[index]
            || coordinate.pair_component_index_b != pair_b[index]
            || coordinate.unit == nullptr
            || coordinate.unit != std::string(units[index])) {
            throw std::runtime_error(
                "provider capability coordinate order does not match the supported v1 contract"
            );
        }
    }
    const auto valid_fingerprint = [](const char* value) {
        return value != nullptr
            && strnlen(value, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE) == 71
            && std::strncmp(value, "sha256:", 7) == 0;
    };
    if (!valid_fingerprint(descriptor.parameter_fingerprint)
        || !valid_fingerprint(descriptor.topology_fingerprint)) {
        throw std::runtime_error(
            "provider capability fingerprints are incomplete"
        );
    }
}

PyObject* string_tuple(
    const char* const* values, std::size_t count
) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(count));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < count; ++index) {
        if (values[index] == nullptr) {
            Py_DECREF(tuple);
            throw std::runtime_error("provider descriptor contains a null string");
        }
        PyObject* item = PyUnicode_FromString(values[index]);
        if (item == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), item);
    }
    return tuple;
}

PyObject* descriptor_to_python(
    const epcsaft_native_capability_descriptor_v1& descriptor
) {
    validate_descriptor(descriptor);

    PyObject* components = string_tuple(
        descriptor.component_ids, descriptor.component_count
    );
    PyObject* kinds = PyTuple_New(
        static_cast<Py_ssize_t>(descriptor.coordinate_count)
    );
    PyObject* units = PyTuple_New(
        static_cast<Py_ssize_t>(descriptor.coordinate_count)
    );
    if (components == nullptr || kinds == nullptr || units == nullptr) {
        Py_XDECREF(components);
        Py_XDECREF(kinds);
        Py_XDECREF(units);
        return nullptr;
    }
    for (std::size_t index = 0; index < descriptor.coordinate_count; ++index) {
        const auto& coordinate = descriptor.coordinates[index];
        if (coordinate.struct_size
                != sizeof(epcsaft_native_capability_coordinate_v1)
            || coordinate.unit == nullptr) {
            Py_DECREF(components);
            Py_DECREF(kinds);
            Py_DECREF(units);
            throw std::runtime_error(
                "provider capability coordinate metadata is incomplete"
            );
        }
        PyObject* kind = PyUnicode_FromString(coordinate_kind(coordinate.kind));
        PyObject* unit = PyUnicode_FromString(coordinate.unit);
        if (kind == nullptr || unit == nullptr) {
            Py_XDECREF(kind);
            Py_XDECREF(unit);
            Py_DECREF(components);
            Py_DECREF(kinds);
            Py_DECREF(units);
            return nullptr;
        }
        PyTuple_SET_ITEM(kinds, static_cast<Py_ssize_t>(index), kind);
        PyTuple_SET_ITEM(units, static_cast<Py_ssize_t>(index), unit);
    }

    const std::size_t parameter_length = strnlen(
        descriptor.parameter_fingerprint,
        EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    const std::size_t topology_length = strnlen(
        descriptor.topology_fingerprint,
        EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    const std::size_t basis_length = strnlen(
        descriptor.helmholtz_basis_id,
        EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    PyObject* result = Py_BuildValue(
        "(ssNNNs#s#issddssssnns#ss)",
        capability_id(descriptor.capability),
        parameter_family(descriptor.parameter_family),
        components,
        kinds,
        units,
        descriptor.parameter_fingerprint,
        static_cast<Py_ssize_t>(parameter_length),
        descriptor.topology_fingerprint,
        static_cast<Py_ssize_t>(topology_length),
        static_cast<int>(descriptor.derivative_order),
        "DERIVATIVE_READY",
        "NONE",
        descriptor.temperature_min_k,
        descriptor.temperature_max_k,
        "unordered_component_pair",
        "fixed_composition_helmholtz_phase",
        "neutral_nonassociating_binary",
        "row_major",
        static_cast<Py_ssize_t>(descriptor.state_coordinate_count),
        static_cast<Py_ssize_t>(descriptor.active_parameter_count),
        descriptor.helmholtz_basis_id,
        static_cast<Py_ssize_t>(basis_length),
        "UNSUPPORTED_MODEL",
        "DOMAIN_ERROR"
    );
    return result;
}

Phase evaluate_phase(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const Row& row,
    double first_fraction,
    double volume,
    double parameter
) {
    std::array<double, 2> amounts = {
        first_fraction, 1.0 - first_fraction
    };
    Phase phase{};
    epcsaft_mixture_phase_block_result_v1 result{};
    result.struct_size = sizeof(result);
    result.coordinate_count = phase.gradient.size();
    result.gradient_capacity = phase.gradient.size();
    result.hessian_capacity = phase.hessian.size();
    result.gradient = phase.gradient.data();
    result.hessian = phase.hessian.data();
    const int status = table.evaluate_mixture_phase_kij(
        table.model_context,
        row.temperature,
        amounts.data(),
        amounts.size(),
        volume,
        parameter,
        &result
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
        throw std::runtime_error(
            std::string("Provider phase evaluation failed: ") + result.error
        );
    }
    if (!bounded_field_equal(
            payload.parameter_fingerprint, result.parameter_fingerprint
        )) {
        throw std::runtime_error(
            "Provider evaluation parameter fingerprint changed"
        );
    }
    phase.pressure = result.pressure_pa;
    return phase;
}

Evaluation evaluate_problem(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const std::vector<double>& variables
) {
    const std::size_t row_count = payload.training_rows.size();
    const std::size_t variable_count = 1 + 2 * row_count;
    if (variables.size() != variable_count) {
        throw std::invalid_argument(
            "solver variables do not match the training-row dimension"
        );
    }
    Evaluation evaluation{};
    evaluation.residuals.resize(4 * row_count);
    evaluation.jacobian.assign(4 * row_count * variable_count, 0.0);
    const double parameter =
        payload.parameter_origin + payload.parameter_scale * variables[0];
    for (std::size_t row_index = 0; row_index < row_count; ++row_index) {
        const Row& row = payload.training_rows[row_index];
        const std::size_t liquid_column = 1 + 2 * row_index;
        const std::size_t vapor_column = liquid_column + 1;
        const double liquid_volume =
            row.liquid_volume_origin * std::exp(variables[liquid_column]);
        const double vapor_volume =
            row.vapor_volume_origin * std::exp(variables[vapor_column]);
        const Phase liquid = evaluate_phase(
            table, payload, row, row.liquid_first, liquid_volume, parameter
        );
        const Phase vapor = evaluate_phase(
            table, payload, row, row.vapor_first, vapor_volume, parameter
        );
        const std::size_t residual_offset = 4 * row_index;
        evaluation.residuals[residual_offset] =
            (liquid.pressure - row.pressure) / row.pressure_scale;
        evaluation.residuals[residual_offset + 1] =
            (vapor.pressure - row.pressure) / row.pressure_scale;
        evaluation.residuals[residual_offset + 2] =
            (liquid.gradient[0] - vapor.gradient[0])
            / row.chemical_potential_scales[0];
        evaluation.residuals[residual_offset + 3] =
            (liquid.gradient[1] - vapor.gradient[1])
            / row.chemical_potential_scales[1];

        auto jacobian = [&](std::size_t residual, std::size_t column)
            -> double& {
            return evaluation.jacobian[
                residual * variable_count + column
            ];
        };
        jacobian(residual_offset, 0) =
            -gas_constant * row.temperature * liquid.hessian[2 * 4 + 3]
            * payload.parameter_scale / row.pressure_scale;
        jacobian(residual_offset, liquid_column) =
            -gas_constant * row.temperature * liquid.hessian[2 * 4 + 2]
            * liquid_volume / row.pressure_scale;
        jacobian(residual_offset + 1, 0) =
            -gas_constant * row.temperature * vapor.hessian[2 * 4 + 3]
            * payload.parameter_scale / row.pressure_scale;
        jacobian(residual_offset + 1, vapor_column) =
            -gas_constant * row.temperature * vapor.hessian[2 * 4 + 2]
            * vapor_volume / row.pressure_scale;
        for (std::size_t component = 0; component < 2; ++component) {
            const std::size_t residual = residual_offset + 2 + component;
            const double scale = row.chemical_potential_scales[component];
            jacobian(residual, 0) =
                (liquid.hessian[component * 4 + 3]
                 - vapor.hessian[component * 4 + 3])
                * payload.parameter_scale / scale;
            jacobian(residual, liquid_column) =
                liquid.hessian[component * 4 + 2] * liquid_volume / scale;
            jacobian(residual, vapor_column) =
                -vapor.hessian[component * 4 + 2] * vapor_volume / scale;
        }
    }
    return evaluation;
}

class GeneralKijCost final : public ceres::CostFunction {
public:
    GeneralKijCost(
        const epcsaft_native_sdk_v1* table, const Payload& payload
    ) : table_(table), payload_(payload) {
        set_num_residuals(static_cast<int>(4 * payload.training_rows.size()));
        mutable_parameter_block_sizes()->push_back(
            static_cast<int>(1 + 2 * payload.training_rows.size())
        );
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            const std::size_t count = static_cast<std::size_t>(
                parameter_block_sizes()[0]
            );
            std::vector<double> variables(values[0], values[0] + count);
            const Evaluation evaluation = evaluate_problem(
                *table_, payload_, variables
            );
            std::copy(
                evaluation.residuals.begin(),
                evaluation.residuals.end(),
                residuals
            );
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(
                    evaluation.jacobian.begin(),
                    evaluation.jacobian.end(),
                    jacobians[0]
                );
            }
            failure_reason_.clear();
            return true;
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
    }

    const std::string& failure_reason() const noexcept {
        return failure_reason_;
    }

private:
    const epcsaft_native_sdk_v1* table_;
    const Payload& payload_;
    mutable std::string failure_reason_;
};

MatrixDiagnostics matrix_diagnostics(const Eigen::MatrixXd& matrix) {
    MatrixDiagnostics diagnostics{};
    const Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(
        matrix, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    const Eigen::VectorXd singular = decomposition.singularValues();
    diagnostics.singular_values.assign(
        singular.data(), singular.data() + singular.size()
    );
    if (singular.size() == 0 || !std::isfinite(singular[0])) {
        return diagnostics;
    }
    const double threshold =
        100.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(matrix.rows(), matrix.cols()))
        * singular[0];
    for (Eigen::Index index = 0; index < singular.size(); ++index) {
        if (singular[index] > threshold) {
            ++diagnostics.rank;
        }
    }
    if (diagnostics.rank > 0) {
        diagnostics.condition_number =
            singular[0] / singular[diagnostics.rank - 1];
    }
    return diagnostics;
}

void diagnose_jacobian(SolveOutcome& outcome) {
    const Eigen::Index residual_count = static_cast<Eigen::Index>(
        outcome.evaluation.residuals.size()
    );
    const Eigen::Index variable_count = static_cast<Eigen::Index>(
        outcome.variables.size()
    );
    const Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    > full(
        outcome.evaluation.jacobian.data(), residual_count, variable_count
    );
    outcome.full_jacobian = matrix_diagnostics(full);
    const Eigen::MatrixXd nuisance = full.rightCols(variable_count - 1);
    const Eigen::JacobiSVD<Eigen::MatrixXd> nuisance_svd(
        nuisance, Eigen::ComputeThinU
    );
    const Eigen::VectorXd nuisance_singular = nuisance_svd.singularValues();
    Eigen::Index nuisance_rank = 0;
    if (nuisance_singular.size() > 0 && std::isfinite(nuisance_singular[0])) {
        const double threshold =
            100.0 * std::numeric_limits<double>::epsilon()
            * static_cast<double>(
                std::max(nuisance.rows(), nuisance.cols())
            )
            * nuisance_singular[0];
        while (nuisance_rank < nuisance_singular.size()
               && nuisance_singular[nuisance_rank] > threshold) {
            ++nuisance_rank;
        }
    }
    Eigen::MatrixXd projected = full.leftCols(1);
    if (nuisance_rank > 0) {
        const Eigen::MatrixXd basis =
            nuisance_svd.matrixU().leftCols(nuisance_rank);
        projected -= basis * (basis.transpose() * projected);
    }
    outcome.projected_parameter_jacobian = matrix_diagnostics(projected);
}

ceres::Solver::Options solver_options(const Payload& payload) {
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.maximum_iterations;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    return options;
}

SolveOutcome solve_training(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    double physical_start
) {
    SolveOutcome outcome{};
    const std::size_t variable_count = 1 + 2 * payload.training_rows.size();
    outcome.variables.resize(variable_count);
    outcome.variables[0] =
        (physical_start - payload.parameter_origin) / payload.parameter_scale;
    for (std::size_t index = 0; index < payload.training_rows.size(); ++index) {
        const Row& row = payload.training_rows[index];
        outcome.variables[1 + 2 * index] =
            std::log(row.liquid_volume_start / row.liquid_volume_origin);
        outcome.variables[2 + 2 * index] =
            std::log(row.vapor_volume_start / row.vapor_volume_origin);
    }
    GeneralKijCost cost(&table, payload);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, outcome.variables.data());
    const double parameter_lower =
        (payload.parameter_lower_bound - payload.parameter_origin)
        / payload.parameter_scale;
    const double parameter_upper =
        (payload.parameter_upper_bound - payload.parameter_origin)
        / payload.parameter_scale;
    problem.SetParameterLowerBound(
        outcome.variables.data(), 0, std::min(parameter_lower, parameter_upper)
    );
    problem.SetParameterUpperBound(
        outcome.variables.data(), 0, std::max(parameter_lower, parameter_upper)
    );
    for (std::size_t index = 0; index < payload.training_rows.size(); ++index) {
        const Row& row = payload.training_rows[index];
        problem.SetParameterLowerBound(
            outcome.variables.data(),
            static_cast<int>(1 + 2 * index),
            std::log(
                row.liquid_volume_bounds[0] / row.liquid_volume_origin
            )
        );
        problem.SetParameterUpperBound(
            outcome.variables.data(),
            static_cast<int>(1 + 2 * index),
            std::log(
                row.liquid_volume_bounds[1] / row.liquid_volume_origin
            )
        );
        problem.SetParameterLowerBound(
            outcome.variables.data(),
            static_cast<int>(2 + 2 * index),
            std::log(
                row.vapor_volume_bounds[0] / row.vapor_volume_origin
            )
        );
        problem.SetParameterUpperBound(
            outcome.variables.data(),
            static_cast<int>(2 + 2 * index),
            std::log(
                row.vapor_volume_bounds[1] / row.vapor_volume_origin
            )
        );
    }
    ceres::Solve(
        solver_options(payload), &problem, &outcome.summary
    );
    outcome.failure_reason = cost.failure_reason();
    try {
        outcome.evaluation = evaluate_problem(table, payload, outcome.variables);
        diagnose_jacobian(outcome);
    } catch (const std::exception& error) {
        if (outcome.failure_reason.empty()) {
            outcome.failure_reason = error.what();
        }
    }
    return outcome;
}

class ReportingCost final : public ceres::CostFunction {
public:
    ReportingCost(
        const epcsaft_native_sdk_v1* table,
        const Payload& payload,
        double parameter_solver_value
    ) : table_(table), payload_(payload), parameter_(parameter_solver_value) {
        set_num_residuals(4);
        mutable_parameter_block_sizes()->push_back(2);
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            const std::vector<double> variables = {
                parameter_, values[0][0], values[0][1]
            };
            const Evaluation evaluation = evaluate_problem(
                *table_, payload_, variables
            );
            std::copy(
                evaluation.residuals.begin(),
                evaluation.residuals.end(),
                residuals
            );
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                for (std::size_t row = 0; row < 4; ++row) {
                    jacobians[0][2 * row] =
                        evaluation.jacobian[3 * row + 1];
                    jacobians[0][2 * row + 1] =
                        evaluation.jacobian[3 * row + 2];
                }
            }
            failure_reason_.clear();
            return true;
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
    }

    const std::string& failure_reason() const noexcept {
        return failure_reason_;
    }

private:
    const epcsaft_native_sdk_v1* table_;
    const Payload& payload_;
    double parameter_;
    mutable std::string failure_reason_;
};

RowOutcome solve_reporting(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const Row& row,
    double parameter_solver_value
) {
    Payload row_payload = payload;
    row_payload.training_rows = {row};
    row_payload.reporting_rows.clear();
    std::array<double, 2> variables = {
        std::log(row.liquid_volume_start / row.liquid_volume_origin),
        std::log(row.vapor_volume_start / row.vapor_volume_origin),
    };
    ReportingCost cost(&table, row_payload, parameter_solver_value);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, variables.data());
    problem.SetParameterLowerBound(
        variables.data(), 0,
        std::log(row.liquid_volume_bounds[0] / row.liquid_volume_origin)
    );
    problem.SetParameterUpperBound(
        variables.data(), 0,
        std::log(row.liquid_volume_bounds[1] / row.liquid_volume_origin)
    );
    problem.SetParameterLowerBound(
        variables.data(), 1,
        std::log(row.vapor_volume_bounds[0] / row.vapor_volume_origin)
    );
    problem.SetParameterUpperBound(
        variables.data(), 1,
        std::log(row.vapor_volume_bounds[1] / row.vapor_volume_origin)
    );
    ceres::Solver::Summary summary;
    ceres::Solve(solver_options(payload), &problem, &summary);
    RowOutcome outcome{
        row,
        row.liquid_volume_origin * std::exp(variables[0]),
        row.vapor_volume_origin * std::exp(variables[1]),
        {},
        summary.IsSolutionUsable(),
        cost.failure_reason(),
    };
    try {
        const Evaluation evaluation = evaluate_problem(
            table,
            row_payload,
            {parameter_solver_value, variables[0], variables[1]}
        );
        std::copy(
            evaluation.residuals.begin(),
            evaluation.residuals.end(),
            outcome.residuals.begin()
        );
    } catch (const std::exception& error) {
        outcome.usable = false;
        outcome.failure_reason = error.what();
    }
    return outcome;
}

std::string termination_name(ceres::TerminationType termination) {
    switch (termination) {
        case ceres::CONVERGENCE:
            return "CONVERGENCE";
        case ceres::NO_CONVERGENCE:
            return "NO_CONVERGENCE";
        case ceres::FAILURE:
            return "FAILURE";
        case ceres::USER_SUCCESS:
            return "USER_SUCCESS";
        case ceres::USER_FAILURE:
            return "USER_FAILURE";
    }
    return "UNKNOWN";
}

PyObject* doubles_to_tuple(const std::vector<double>& values) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyFloat_FromDouble(values[index]);
        if (value == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* row_to_python(
    const Row& row,
    double liquid_volume,
    double vapor_volume,
    const std::vector<double>& residuals,
    bool usable,
    const std::string& failure_reason
) {
    PyObject* residual_tuple = doubles_to_tuple(residuals);
    if (residual_tuple == nullptr) return nullptr;
    return Py_BuildValue(
        "(ssddNOs)",
        row.row_id.c_str(),
        row.partition.c_str(),
        liquid_volume,
        vapor_volume,
        residual_tuple,
        usable ? Py_True : Py_False,
        failure_reason.c_str()
    );
}

bool complete_evaluation(
    const SolveOutcome& outcome, const Payload& payload
) {
    const std::size_t residual_count = 4 * payload.training_rows.size();
    return outcome.evaluation.residuals.size() == residual_count
        && outcome.evaluation.jacobian.size()
            == residual_count * outcome.variables.size();
}

PyObject* rows_to_python(
    const Payload& payload,
    const SolveOutcome& primary,
    const epcsaft_native_sdk_v1& table
) {
    const std::size_t count =
        payload.training_rows.size() + payload.reporting_rows.size();
    PyObject* result = PyTuple_New(static_cast<Py_ssize_t>(count));
    if (result == nullptr) return nullptr;
    const bool primary_evaluation_available =
        complete_evaluation(primary, payload);
    const std::string unavailable_reason = primary.failure_reason.empty()
        ? "primary exact evaluation is unavailable"
        : primary.failure_reason;
    for (std::size_t index = 0; index < payload.training_rows.size(); ++index) {
        const Row& row = payload.training_rows[index];
        const double liquid_volume = primary_evaluation_available
            ? row.liquid_volume_origin
                * std::exp(primary.variables[1 + 2 * index])
            : row.liquid_volume_start;
        const double vapor_volume = primary_evaluation_available
            ? row.vapor_volume_origin
                * std::exp(primary.variables[2 + 2 * index])
            : row.vapor_volume_start;
        const std::vector<double> residuals = primary_evaluation_available
            ? std::vector<double>(
                primary.evaluation.residuals.begin() + 4 * index,
                primary.evaluation.residuals.begin() + 4 * index + 4
            )
            : std::vector<double>(
                4, std::numeric_limits<double>::quiet_NaN()
            );
        PyObject* item = row_to_python(
            row,
            liquid_volume,
            vapor_volume,
            residuals,
            primary_evaluation_available,
            primary_evaluation_available ? "" : unavailable_reason
        );
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), item);
    }
    for (std::size_t index = 0; index < payload.reporting_rows.size(); ++index) {
        const Row& row = payload.reporting_rows[index];
        const RowOutcome outcome = primary_evaluation_available
            ? solve_reporting(table, payload, row, primary.variables[0])
            : RowOutcome{
                row,
                row.liquid_volume_start,
                row.vapor_volume_start,
                {
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                },
                false,
                unavailable_reason,
            };
        PyObject* item = row_to_python(
            outcome.row,
            outcome.liquid_volume,
            outcome.vapor_volume,
            std::vector<double>(
                outcome.residuals.begin(), outcome.residuals.end()
            ),
            outcome.usable,
            outcome.failure_reason
        );
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(
            result,
            static_cast<Py_ssize_t>(payload.training_rows.size() + index),
            item
        );
    }
    return result;
}

}  // namespace

PyObject* parameter_capabilities_python(PyObject* capsule) {
    try {
        const auto* table = capability_table(capsule);
        PyObject* result = PyTuple_New(
            static_cast<Py_ssize_t>(table->capability_count)
        );
        if (result == nullptr) return nullptr;
        for (std::size_t index = 0; index < table->capability_count; ++index) {
            PyObject* item = descriptor_to_python(table->capabilities[index]);
            if (item == nullptr) {
                Py_DECREF(result);
                return nullptr;
            }
            PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), item);
        }
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* evaluate_general_kij_python(
    PyObject* capsule, PyObject* payload_object, PyObject* variables_object
) {
    try {
        const auto* table = capability_table(capsule);
        const Payload payload = parse_payload(payload_object);
        checked_descriptor(*table, payload);
        const std::vector<double> variables = doubles(
            variables_object, "solver variables"
        );
        const Evaluation evaluation = evaluate_problem(
            *table, payload, variables
        );
        PyObject* residuals = doubles_to_tuple(evaluation.residuals);
        PyObject* jacobian = doubles_to_tuple(evaluation.jacobian);
        if (residuals == nullptr || jacobian == nullptr) {
            Py_XDECREF(residuals);
            Py_XDECREF(jacobian);
            return nullptr;
        }
        PyObject* result = PyTuple_New(2);
        if (result == nullptr) {
            Py_DECREF(residuals);
            Py_DECREF(jacobian);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, 0, residuals);
        PyTuple_SET_ITEM(result, 1, jacobian);
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* solve_general_kij_python(
    PyObject* capsule, PyObject* payload_object
) {
    try {
        const auto* table = capability_table(capsule);
        const Payload payload = parse_payload(payload_object);
        checked_descriptor(*table, payload);
        std::vector<SolveOutcome> outcomes;
        outcomes.reserve(payload.starts.size());
        for (const double start : payload.starts) {
            outcomes.push_back(solve_training(*table, payload, start));
        }
        if (outcomes.empty()) {
            throw std::invalid_argument("at least one parameter start is required");
        }
        const SolveOutcome& primary = outcomes.front();
        double maximum_parameter_delta = 0.0;
        double maximum_cost_delta = 0.0;
        bool confirmations_usable = true;
        for (std::size_t index = 1; index < outcomes.size(); ++index) {
            const SolveOutcome& confirmation = outcomes[index];
            confirmations_usable =
                confirmations_usable
                && confirmation.summary.IsSolutionUsable()
                && confirmation.failure_reason.empty()
                && complete_evaluation(confirmation, payload);
            maximum_parameter_delta = std::max(
                maximum_parameter_delta,
                std::abs(primary.variables[0] - confirmation.variables[0])
            );
            maximum_cost_delta = std::max(
                maximum_cost_delta,
                std::abs(
                    primary.summary.final_cost - confirmation.summary.final_cost
                ) / std::max({
                    std::abs(primary.summary.final_cost),
                    std::abs(confirmation.summary.final_cost),
                    std::numeric_limits<double>::min(),
                })
            );
        }
        PyObject* residuals = doubles_to_tuple(primary.evaluation.residuals);
        PyObject* jacobian = doubles_to_tuple(primary.evaluation.jacobian);
        PyObject* full_singular = doubles_to_tuple(
            primary.full_jacobian.singular_values
        );
        PyObject* projected_singular = doubles_to_tuple(
            primary.projected_parameter_jacobian.singular_values
        );
        PyObject* rows = rows_to_python(payload, primary, *table);
        if (residuals == nullptr || jacobian == nullptr
            || full_singular == nullptr || projected_singular == nullptr
            || rows == nullptr) {
            Py_XDECREF(residuals);
            Py_XDECREF(jacobian);
            Py_XDECREF(full_singular);
            Py_XDECREF(projected_singular);
            Py_XDECREF(rows);
            return nullptr;
        }
        const double physical_parameter =
            payload.parameter_origin
            + payload.parameter_scale * primary.variables[0];
        const double bound_distance = std::min(
            physical_parameter - payload.parameter_lower_bound,
            payload.parameter_upper_bound - physical_parameter
        );
        const double active_tolerance =
            std::sqrt(std::numeric_limits<double>::epsilon())
            * std::max(
                1.0,
                payload.parameter_upper_bound - payload.parameter_lower_bound
            );
        const char* active_bound = "";
        if (std::abs(
                physical_parameter - payload.parameter_lower_bound
            ) <= active_tolerance) {
            active_bound = "lower";
        } else if (std::abs(
                       physical_parameter - payload.parameter_upper_bound
                   ) <= active_tolerance) {
            active_bound = "upper";
        }
        PyObject* result = PyTuple_New(24);
        if (result == nullptr) {
            Py_DECREF(residuals);
            Py_DECREF(jacobian);
            Py_DECREF(full_singular);
            Py_DECREF(projected_singular);
            Py_DECREF(rows);
            return nullptr;
        }
        PyTuple_SET_ITEM(
            result, 0,
            PyUnicode_FromString(
                termination_name(primary.summary.termination_type).c_str()
            )
        );
        PyTuple_SET_ITEM(
            result, 1,
            Py_NewRef(primary.summary.IsSolutionUsable() ? Py_True : Py_False)
        );
        PyTuple_SET_ITEM(
            result, 2, PyFloat_FromDouble(primary.summary.initial_cost)
        );
        PyTuple_SET_ITEM(
            result, 3, PyFloat_FromDouble(primary.summary.final_cost)
        );
        PyTuple_SET_ITEM(
            result, 4,
            PyLong_FromSize_t(primary.summary.iterations.size())
        );
        PyTuple_SET_ITEM(result, 5, PyFloat_FromDouble(physical_parameter));
        PyTuple_SET_ITEM(result, 6, PyFloat_FromDouble(bound_distance));
        PyTuple_SET_ITEM(result, 7, PyUnicode_FromString(active_bound));
        PyTuple_SET_ITEM(result, 8, residuals);
        PyTuple_SET_ITEM(result, 9, jacobian);
        PyTuple_SET_ITEM(result, 10, full_singular);
        PyTuple_SET_ITEM(
            result, 11, PyLong_FromLong(primary.full_jacobian.rank)
        );
        PyTuple_SET_ITEM(
            result, 12,
            PyFloat_FromDouble(primary.full_jacobian.condition_number)
        );
        PyTuple_SET_ITEM(result, 13, projected_singular);
        PyTuple_SET_ITEM(
            result, 14,
            PyLong_FromLong(primary.projected_parameter_jacobian.rank)
        );
        PyTuple_SET_ITEM(
            result, 15,
            PyFloat_FromDouble(
                primary.projected_parameter_jacobian.condition_number
            )
        );
        PyTuple_SET_ITEM(
            result, 16, PyLong_FromSize_t(outcomes.size() - 1)
        );
        PyTuple_SET_ITEM(
            result, 17, PyFloat_FromDouble(maximum_parameter_delta)
        );
        PyTuple_SET_ITEM(
            result, 18, PyFloat_FromDouble(maximum_cost_delta)
        );
        PyTuple_SET_ITEM(
            result, 19,
            Py_NewRef(confirmations_usable ? Py_True : Py_False)
        );
        PyTuple_SET_ITEM(result, 20, rows);
        PyTuple_SET_ITEM(
            result, 21,
            PyUnicode_FromString(primary.failure_reason.c_str())
        );
        PyTuple_SET_ITEM(
            result, 22, PyLong_FromSize_t(primary.variables.size())
        );
        PyTuple_SET_ITEM(
            result, 23,
            PyLong_FromSize_t(primary.evaluation.residuals.size())
        );
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_regression
