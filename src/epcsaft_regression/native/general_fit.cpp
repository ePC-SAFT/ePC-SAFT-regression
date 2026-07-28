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

enum class ObservationKind {
    pair_phase,
    pure_phase,
    mean_ionic_activity,
    solvation_gibbs,
};

struct Row final {
    ObservationKind kind{ObservationKind::pair_phase};
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
    double molar_mass{0.0};
    double liquid_density{0.0};
    double liquid_density_scale{0.0};
    double direct_state{0.0};
    double observed_value{0.0};
    double direct_scale{0.0};
};

struct Payload final {
    std::string capability_id;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::vector<std::string> component_ids;
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
    std::size_t coordinate_count;
};

struct Evaluation final {
    std::vector<double> residuals;
    std::vector<double> jacobian;
    std::vector<double> modeled_values;
    std::vector<double> provider_derivatives;
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

bool direct_observation(const Payload& payload) {
    return payload.capability_id == "ion_solvation_born_v1"
        || payload.capability_id == "aqueous_solvation_factor_miac_v1";
}

std::size_t residual_count(const Payload& payload) {
    return (direct_observation(payload) ? 1u : 4u)
        * payload.training_rows.size();
}

std::size_t variable_count(const Payload& payload) {
    return direct_observation(payload)
        ? 1u
        : 1u + 2u * payload.training_rows.size();
}

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

Row parse_row(PyObject* object, ObservationKind kind) {
    OwnedPyObject sequence{
        PySequence_Fast(object, "observation payload must be a sequence")
    };
    const bool direct = kind == ObservationKind::mean_ionic_activity
        || kind == ObservationKind::solvation_gibbs;
    const Py_ssize_t expected_size = direct ? 7 : 17;
    if (!sequence
        || PySequence_Fast_GET_SIZE(sequence.get()) != expected_size) {
        PyErr_Clear();
        throw std::invalid_argument(
            "observation payload has the wrong field count"
        );
    }
    auto item = [&](Py_ssize_t index) {
        return PySequence_Fast_GET_ITEM(sequence.get(), index);
    };
    Row row{};
    row.kind = kind;
    row.row_id = text(item(0), "row id");
    row.partition = text(item(1), "partition");
    row.temperature = number(item(2), "temperature");
    row.pressure = number(item(3), "pressure");
    if (direct) {
        row.direct_state = number(item(4), "direct observation state");
        row.observed_value = number(item(5), "observed value");
        row.direct_scale = number(item(6), "direct residual scale");
        const bool valid = row.temperature > 0.0 && row.pressure > 0.0
            && row.direct_scale > 0.0
            && (kind != ObservationKind::mean_ionic_activity
                || (row.direct_state > 0.0 && row.observed_value > 0.0));
        if (!valid) {
            throw std::invalid_argument(
                "direct observation state, value, and scale are invalid"
            );
        }
        return row;
    }
    const bool pure = kind == ObservationKind::pure_phase;
    if (pure) {
        row.pressure_scale = number(item(4), "pressure scale");
        row.chemical_potential_scales = {
            number(item(5), "chemical-potential scale"),
            number(item(8), "liquid-density scale"),
        };
        row.molar_mass = number(item(6), "molar mass");
        row.liquid_density = number(item(7), "liquid density");
        row.liquid_density_scale = row.chemical_potential_scales[1];
        row.liquid_first = 1.0;
        row.vapor_first = 1.0;
    } else {
        row.liquid_first = number(item(4), "liquid composition");
        row.vapor_first = number(item(5), "vapor composition");
        row.pressure_scale = number(item(6), "pressure scale");
        row.chemical_potential_scales = {
            number(item(7), "first chemical-potential scale"),
            number(item(8), "second chemical-potential scale"),
        };
    }
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
    const auto valid_volume = [](double origin, double start, const auto& bounds) {
        return origin > 0.0 && start > 0.0 && bounds[0] > 0.0
            && bounds[0] < bounds[1] && start >= bounds[0]
            && start <= bounds[1];
    };
    const bool common_valid = row.temperature > 0.0 && row.pressure > 0.0
        && row.pressure_scale > 0.0
        && valid_volume(
            row.liquid_volume_origin,
            row.liquid_volume_start,
            row.liquid_volume_bounds
        )
        && valid_volume(
            row.vapor_volume_origin,
            row.vapor_volume_start,
            row.vapor_volume_bounds
        );
    const bool observation_valid = pure
        ? row.chemical_potential_scales[0] > 0.0
            && row.molar_mass > 0.0 && row.liquid_density > 0.0
            && row.liquid_density_scale > 0.0
        : row.liquid_first > 0.0 && row.liquid_first < 1.0
            && row.vapor_first > 0.0 && row.vapor_first < 1.0
            && row.chemical_potential_scales[0] > 0.0
            && row.chemical_potential_scales[1] > 0.0;
    if (!common_valid || !observation_valid) {
        throw std::invalid_argument(
            "observation scales, state values, and volume contracts must be "
            "positive and ordered"
        );
    }
    return row;
}

std::vector<Row> parse_rows(
    PyObject* object, const char* name, ObservationKind kind
) {
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
            PySequence_Fast_GET_ITEM(sequence.get(), index), kind
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
    if (!components
        || (PySequence_Fast_GET_SIZE(components.get()) < 1
            || PySequence_Fast_GET_SIZE(components.get()) > 3)) {
        PyErr_Clear();
        throw std::invalid_argument(
            "component ids must contain one to three values"
        );
    }
    Payload payload{};
    payload.capability_id = text(item(0), "capability id");
    payload.parameter_fingerprint = text(item(1), "parameter fingerprint");
    payload.topology_fingerprint = text(item(2), "topology fingerprint");
    const Py_ssize_t component_count = PySequence_Fast_GET_SIZE(components.get());
    payload.component_ids.reserve(static_cast<std::size_t>(component_count));
    for (Py_ssize_t index = 0; index < component_count; ++index) {
        payload.component_ids.push_back(text(
            PySequence_Fast_GET_ITEM(components.get(), index), "component id"
        ));
    }
    const ObservationKind observation_kind =
        payload.capability_id == "ion_solvation_born_v1"
        ? ObservationKind::solvation_gibbs
        : payload.capability_id == "aqueous_solvation_factor_miac_v1"
            ? ObservationKind::mean_ionic_activity
            : component_count == 1
                ? ObservationKind::pure_phase
                : ObservationKind::pair_phase;
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
    payload.training_rows = parse_rows(
        item(16), "training rows", observation_kind
    );
    payload.reporting_rows = parse_rows(
        item(17), "reporting rows", observation_kind
    );
    if (payload.training_rows.empty()) {
        throw std::invalid_argument("at least one training row is required");
    }
    const bool starts_valid = payload.starts.size() >= 2
        && std::all_of(
            payload.starts.begin(),
            payload.starts.end(),
            [&](double start) {
                return start >= payload.parameter_lower_bound
                    && start <= payload.parameter_upper_bound;
            }
        );
    if (payload.parameter_scale == 0.0
        || payload.parameter_lower_bound >= payload.parameter_upper_bound
        || !starts_valid || payload.maximum_condition_number <= 0.0
        || payload.function_tolerance <= 0.0
        || payload.gradient_tolerance <= 0.0
        || payload.parameter_tolerance <= 0.0
        || payload.confirmation_parameter_delta <= 0.0
        || payload.confirmation_cost_delta <= 0.0) {
        throw std::invalid_argument(
            "parameter, solver, and confirmation contracts are invalid"
        );
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
    const epcsaft_native_capability_descriptor_v1* selected = nullptr;
    for (std::size_t index = 0; index < table.capability_count; ++index) {
        const auto& candidate = table.capabilities[index];
        const bool supported =
            candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1
            || candidate.capability
                == EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1;
        if (!supported) {
            continue;
        }
        validate_descriptor(candidate);
        if (payload.capability_id != capability_id(candidate.capability)) {
            continue;
        }
        if (selected != nullptr) {
            throw std::runtime_error(
                "provider advertises a duplicate requested capability"
            );
        }
        selected = &candidate;
    }
    if (selected == nullptr) {
        throw std::runtime_error(
            "provider does not advertise the requested capability"
        );
    }
    const auto& descriptor = *selected;
    const bool kij = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1;
    const bool lij = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1;
    const bool born = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1;
    const bool solvation_factor = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1;
    const bool binary = kij || lij;
    const bool direct = born || solvation_factor;
    const bool callback_available = kij
        ? table.evaluate_mixture_phase_kij != nullptr
        : lij
            ? (table.table_size
                    >= offsetof(
                        epcsaft_native_sdk_v1,
                        evaluate_mixture_phase_lij
                    ) + sizeof(table.evaluate_mixture_phase_lij)
                && table.evaluate_mixture_phase_lij != nullptr)
            : born
                ? table.evaluate_ion_solvation_born != nullptr
                    && table.ion_solvation_born_result_size
                        == sizeof(epcsaft_ion_solvation_born_result_v1)
                : solvation_factor
                    ? table.evaluate_aqueous_miac_solvation_factor_batch
                            != nullptr
                        && table.aqueous_miac_solvation_factor_result_size
                            == sizeof(
                                epcsaft_aqueous_miac_solvation_factor_result_v1
                            )
                    : (table.table_size
                    >= offsetof(
                        epcsaft_native_sdk_v1,
                        evaluate_pure_phase_parameter
                    ) + sizeof(table.evaluate_pure_phase_parameter)
                && table.evaluate_pure_phase_parameter != nullptr);
    const std::size_t expected_component_count =
        direct ? 3u : binary ? 2u : 1u;
    const bool components_match =
        payload.component_ids.size() == expected_component_count
        && descriptor.component_count == expected_component_count
        && descriptor.component_ids != nullptr
        && std::equal(
            payload.component_ids.begin(),
            payload.component_ids.end(),
            descriptor.component_ids,
            [](const std::string& expected, const char* observed) {
                return observed != nullptr && expected == observed;
            }
        );
    if (!bounded_field_equal(
            payload.parameter_fingerprint, descriptor.parameter_fingerprint
        )
        || !bounded_field_equal(
            payload.topology_fingerprint, descriptor.topology_fingerprint
        )
        || !components_match
        || !callback_available
        || (!direct
            && table.mixture_result_size
                != sizeof(epcsaft_mixture_phase_block_result_v1))) {
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
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1) {
        return "neutral_binary_phase_lij_v1";
    }
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1) {
        return "neutral_pure_segment_count_v1";
    }
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1) {
        return "neutral_pure_segment_diameter_v1";
    }
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1) {
        return "neutral_pure_dispersion_energy_over_k_v1";
    }
    if (value == EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1) {
        return "ion_solvation_born_v1";
    }
    if (value
        == EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1) {
        return "aqueous_solvation_factor_miac_v1";
    }
    throw std::runtime_error("provider advertised an unknown capability");
}

const char* parameter_family(std::uint32_t value) {
    if (value
        == EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1) {
        return "k_ij";
    }
    if (value
        == EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_LIJ_V1) {
        return "l_ij";
    }
    if (value == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1) {
        return "segment_count";
    }
    if (value == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1) {
        return "segment_diameter";
    }
    if (value
        == EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1) {
        return "dispersion_energy_over_k";
    }
    if (value == EPCSAFT_NATIVE_PARAMETER_FAMILY_BORN_DIAMETER_V1) {
        return "born_diameter";
    }
    if (value == EPCSAFT_NATIVE_PARAMETER_FAMILY_SOLVATION_FACTOR_V1) {
        return "solvation_factor";
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
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_LIJ_V1:
            return "l_ij";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_COUNT_V1:
            return "segment_count";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1:
            return "segment_diameter";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1:
            return "dispersion_energy_over_k";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BORN_DIAMETER_V1:
            return "born_diameter";
        case EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SOLVATION_FACTOR_V1:
            return "solvation_factor";
        default:
            throw std::runtime_error(
                "provider advertised an unknown capability coordinate"
            );
    }
}

void validate_descriptor(
    const epcsaft_native_capability_descriptor_v1& descriptor
) {
    const bool kij = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1;
    const bool lij = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1;
    const bool segment_count = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1;
    const bool segment_diameter = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1;
    const bool dispersion_energy = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1;
    const bool born = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1;
    const bool solvation_factor = descriptor.capability
        == EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1;
    const bool binary = kij || lij;
    const bool pure = segment_count || segment_diameter || dispersion_energy;
    const bool direct = born || solvation_factor;
    const bool matching_family =
        (kij
         && descriptor.parameter_family
             == EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1)
        || (lij
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_LIJ_V1)
        || (segment_count
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1)
        || (segment_diameter
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1)
        || (dispersion_energy
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1)
        || (born
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_BORN_DIAMETER_V1)
        || (solvation_factor
            && descriptor.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_SOLVATION_FACTOR_V1);
    const std::uint32_t expected_observation = born
        ? EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1
        : solvation_factor
            ? EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1
            : EPCSAFT_NATIVE_OBSERVATION_FIXED_COMPOSITION_HELMHOLTZ_PHASE_V1;
    const std::uint32_t expected_domain = binary
        ? EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_BINARY_V1
        : pure
            ? EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_PURE_V1
            : born
                ? EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_WATER_SINGLE_ION_V1
                : EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_AQUEOUS_NABR_V1;
    const std::size_t expected_component_count =
        binary ? 2u : pure ? 1u : 3u;
    const std::size_t expected_state_count =
        binary ? 3u : pure ? 2u : 0u;
    const std::size_t expected_coordinate_count =
        binary ? 4u : pure ? 3u : 1u;
    if (descriptor.struct_size
            != sizeof(epcsaft_native_capability_descriptor_v1)
        || descriptor.schema_version
            != EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1
        || (!binary && !pure && !direct)
        || !matching_family
        || descriptor.parameter_identity
            != (binary
                    ? EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1
                    : EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1)
        || descriptor.observation_contract != expected_observation
        || descriptor.model_domain != expected_domain
        || descriptor.tensor_layout
            != EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1
        || descriptor.derivative_order != (direct ? 1u : 2u)
        || descriptor.maturity
            != EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1
        || descriptor.authority_effect
            != EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1
        || descriptor.unsupported_status
            != EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1
        || descriptor.domain_status != EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1
        || descriptor.state_coordinate_count != expected_state_count
        || descriptor.active_parameter_count != 1
        || descriptor.coordinate_count != expected_coordinate_count
        || descriptor.component_count != expected_component_count
        || descriptor.coordinates == nullptr
        || descriptor.component_ids == nullptr
        || descriptor.component_ids[0] == nullptr
        || (expected_component_count > 1
            && descriptor.component_ids[1] == nullptr)
        || (expected_component_count > 2
            && descriptor.component_ids[2] == nullptr)
        || !std::isfinite(descriptor.temperature_min_k)
        || !std::isfinite(descriptor.temperature_max_k)
        || descriptor.temperature_min_k > descriptor.temperature_max_k
        || descriptor.helmholtz_basis_id
            != std::string(EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID)) {
        throw std::runtime_error(
            "provider capability descriptor does not match the supported v1 contract"
        );
    }
    std::vector<std::uint32_t> kinds;
    std::vector<int> components;
    std::vector<int> pair_a;
    std::vector<int> pair_b;
    std::vector<const char*> units;
    if (direct) {
        kinds.push_back(
            born
                ? EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BORN_DIAMETER_V1
                : EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SOLVATION_FACTOR_V1
        );
        components.push_back(born ? 1 : 0);
        pair_a.push_back(-1);
        pair_b.push_back(-1);
        units.push_back(born ? "angstrom" : "dimensionless");
    } else {
        kinds.push_back(EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1);
        components.push_back(0);
        pair_a.push_back(-1);
        pair_b.push_back(-1);
        units.push_back("mol");
    }
    if (binary) {
        kinds.push_back(EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1);
        components.push_back(1);
        pair_a.push_back(-1);
        pair_b.push_back(-1);
        units.push_back("mol");
    }
    if (!direct) {
        kinds.push_back(EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1);
        components.push_back(-1);
        pair_a.push_back(-1);
        pair_b.push_back(-1);
        units.push_back("m3");
    }
    if (binary) {
        kinds.push_back(kij
            ? EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1
            : EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_LIJ_V1);
        components.push_back(-1);
        pair_a.push_back(0);
        pair_b.push_back(1);
        units.push_back("dimensionless");
    } else if (pure) {
        kinds.push_back(
            segment_count
                ? EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_COUNT_V1
                : segment_diameter
                    ? EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1
                    : EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1
        );
        components.push_back(0);
        pair_a.push_back(-1);
        pair_b.push_back(-1);
        units.push_back(
            segment_count ? "dimensionless"
                          : segment_diameter ? "angstrom" : "kelvin"
        );
    }
    for (std::size_t index = 0; index < kinds.size(); ++index) {
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
    const bool binary = descriptor.capability
            == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1
        || descriptor.capability
            == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1;
    const char* observation_contract =
        descriptor.observation_contract
                == EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1
            ? "ion_solvation_gibbs"
            : descriptor.observation_contract
                    == EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1
                ? "aqueous_mean_ionic_activity"
                : "fixed_composition_helmholtz_phase";
    const char* model_domain =
        descriptor.model_domain
                == EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_WATER_SINGLE_ION_V1
            ? "figiel_water_single_ion"
            : descriptor.model_domain
                    == EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_AQUEOUS_NABR_V1
                ? "figiel_aqueous_nabr"
                : binary
                    ? "neutral_nonassociating_binary"
                    : "neutral_nonassociating_pure";
    const auto& parameter_coordinate =
        descriptor.coordinates[descriptor.coordinate_count - 1];
    const std::size_t active_component_count = binary ? 2u : 1u;
    PyObject* active_components = PyTuple_New(
        static_cast<Py_ssize_t>(active_component_count)
    );
    if (active_components == nullptr) {
        Py_DECREF(components);
        Py_DECREF(kinds);
        Py_DECREF(units);
        return nullptr;
    }
    for (std::size_t index = 0; index < active_component_count; ++index) {
        const int component_index = binary
            ? (index == 0
                   ? parameter_coordinate.pair_component_index_a
                   : parameter_coordinate.pair_component_index_b)
            : parameter_coordinate.component_index;
        PyObject* component = PyUnicode_FromString(
            descriptor.component_ids[component_index]
        );
        if (component == nullptr) {
            Py_DECREF(components);
            Py_DECREF(kinds);
            Py_DECREF(units);
            Py_DECREF(active_components);
            return nullptr;
        }
        PyTuple_SET_ITEM(
            active_components, static_cast<Py_ssize_t>(index), component
        );
    }
    PyObject* result = Py_BuildValue(
        "(ssNNNs#s#issddssssnns#ssN)",
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
        binary ? "unordered_component_pair" : "component",
        observation_contract,
        model_domain,
        "row_major",
        static_cast<Py_ssize_t>(descriptor.state_coordinate_count),
        static_cast<Py_ssize_t>(descriptor.active_parameter_count),
        descriptor.helmholtz_basis_id,
        static_cast<Py_ssize_t>(basis_length),
        "UNSUPPORTED_MODEL",
        "DOMAIN_ERROR",
        active_components
    );
    return result;
}

PyObject* unsupported_descriptor_to_python(
    const epcsaft_native_capability_descriptor_v1& descriptor
) {
    if (descriptor.struct_size
        != sizeof(epcsaft_native_capability_descriptor_v1)) {
        throw std::runtime_error(
            "provider capability descriptor size is unsupported"
        );
    }
    return Py_BuildValue(
        "(III)",
        descriptor.capability,
        descriptor.schema_version,
        descriptor.parameter_family
    );
}

Phase evaluate_phase(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const Row& row,
    double first_fraction,
    double volume,
    double parameter
) {
    Phase phase{};
    const bool pure = payload.component_ids.size() == 1;
    phase.coordinate_count = pure ? 3 : 4;
    phase.gradient.fill(std::numeric_limits<double>::quiet_NaN());
    phase.hessian.fill(std::numeric_limits<double>::quiet_NaN());
    epcsaft_mixture_phase_block_result_v1 result{};
    result.struct_size = sizeof(result);
    result.coordinate_count = phase.coordinate_count;
    result.gradient_capacity = phase.gradient.size();
    result.hessian_capacity = phase.hessian.size();
    result.gradient = phase.gradient.data();
    result.hessian = phase.hessian.data();
    int status = EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
    if (pure) {
        std::uint32_t family = 0;
        if (payload.capability_id == "neutral_pure_segment_count_v1") {
            family = EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1;
        } else if (
            payload.capability_id == "neutral_pure_segment_diameter_v1"
        ) {
            family = EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1;
        } else if (
            payload.capability_id
            == "neutral_pure_dispersion_energy_over_k_v1"
        ) {
            family =
                EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1;
        }
        status = table.evaluate_pure_phase_parameter(
            table.model_context,
            row.temperature,
            1.0,
            volume,
            family,
            parameter,
            &result
        );
    } else {
        const std::array<double, 2> amounts = {
            first_fraction, 1.0 - first_fraction
        };
        status = payload.capability_id == "neutral_binary_phase_kij_v1"
            ? table.evaluate_mixture_phase_kij(
                  table.model_context,
                  row.temperature,
                  amounts.data(),
                  amounts.size(),
                  volume,
                  parameter,
                  &result
              )
            : table.evaluate_mixture_phase_lij(
                  table.model_context,
                  row.temperature,
                  amounts.data(),
                  amounts.size(),
                  volume,
                  parameter,
                  &result
              );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
        const std::size_t error_length = strnlen(
            result.error, EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE
        );
        throw std::runtime_error(
            std::string("Provider phase evaluation failed: ")
            + std::string(result.error, error_length)
        );
    }
    const bool finite_gradient = std::all_of(
        phase.gradient.begin(),
        phase.gradient.begin()
            + static_cast<std::ptrdiff_t>(phase.coordinate_count),
        [](double value) { return std::isfinite(value); }
    );
    const bool finite_hessian = std::all_of(
        phase.hessian.begin(),
        phase.hessian.begin()
            + static_cast<std::ptrdiff_t>(
                phase.coordinate_count * phase.coordinate_count
            ),
        [](double value) { return std::isfinite(value); }
    );
    if (result.coordinate_count != phase.coordinate_count
        || result.gradient_capacity < phase.coordinate_count
        || result.hessian_capacity
            < phase.coordinate_count * phase.coordinate_count
        || result.gradient != phase.gradient.data()
        || result.hessian != phase.hessian.data()
        || !std::isfinite(result.pressure_pa)
        || !finite_gradient || !finite_hessian) {
        throw std::runtime_error(
            "Provider phase derivative result is incomplete or nonfinite"
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

void evaluate_direct_problem(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    double parameter,
    Evaluation& evaluation
) {
    const std::size_t row_count = payload.training_rows.size();
    if (payload.capability_id == "ion_solvation_born_v1") {
        for (std::size_t index = 0; index < row_count; ++index) {
            const Row& row = payload.training_rows[index];
            epcsaft_ion_solvation_born_result_v1 result{};
            result.struct_size = sizeof(result);
            const int status = table.evaluate_ion_solvation_born(
                table.model_context,
                row.temperature,
                row.pressure,
                parameter,
                &result
            );
            if (status != EPCSAFT_NATIVE_STATUS_OK_V1
                || result.status != status
                || !bounded_field_equal(
                    payload.parameter_fingerprint,
                    result.parameter_fingerprint
                )
                || !std::isfinite(result.solvation_gibbs_j_per_mol)
                || !std::isfinite(
                    result.derivative_j_per_mol_per_angstrom
                )
                || !std::isfinite(result.reference_molality_mol_per_kg)
                || !std::isfinite(result.reference_convergence_error)) {
                throw std::runtime_error(
                    std::string("Provider Born evaluation failed: ")
                    + result.error
                );
            }
            evaluation.modeled_values[index] =
                result.solvation_gibbs_j_per_mol;
            evaluation.provider_derivatives[index] =
                result.derivative_j_per_mol_per_angstrom;
            evaluation.residuals[index] =
                (result.solvation_gibbs_j_per_mol - row.observed_value)
                / row.direct_scale;
            evaluation.jacobian[index] =
                result.derivative_j_per_mol_per_angstrom
                * payload.parameter_scale / row.direct_scale;
        }
        return;
    }

    const Row& first = payload.training_rows.front();
    std::vector<double> molalities;
    molalities.reserve(row_count);
    for (const Row& row : payload.training_rows) {
        if (row.temperature != first.temperature
            || row.pressure != first.pressure) {
            throw std::invalid_argument(
                "batched direct observations must share temperature and pressure"
            );
        }
        molalities.push_back(row.direct_state);
    }
    std::vector<epcsaft_aqueous_miac_solvation_factor_result_v1> results(
        row_count
    );
    for (auto& result : results) {
        result.struct_size = sizeof(result);
    }
    const int status = table.evaluate_aqueous_miac_solvation_factor_batch(
        table.model_context,
        payload.parameter_fingerprint.c_str(),
        first.temperature,
        first.pressure,
        molalities.data(),
        molalities.size(),
        parameter,
        results.data(),
        results.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::runtime_error(
            std::string("Provider solvation-factor batch failed: ")
            + results.front().error
        );
    }
    for (std::size_t index = 0; index < row_count; ++index) {
        const Row& row = payload.training_rows[index];
        const auto& result = results[index];
        if (result.status != status
            || !bounded_field_equal(
                payload.parameter_fingerprint,
                result.parameter_fingerprint
            )
            || !std::isfinite(
                result.log_mean_ionic_activity_coefficient_molality
            )
            || !std::isfinite(result.derivative)
            || !std::isfinite(result.reference_molality_mol_per_kg)
            || !std::isfinite(result.reference_convergence_error)
            || !std::isfinite(
                result.reference_derivative_convergence_error
            )) {
            throw std::runtime_error(
                std::string("Provider solvation-factor row failed: ")
                + result.error
            );
        }
        const double modeled = std::exp(
            result.log_mean_ionic_activity_coefficient_molality
        );
        const double ratio = modeled / row.observed_value;
        evaluation.modeled_values[index] = modeled;
        evaluation.provider_derivatives[index] = result.derivative;
        evaluation.residuals[index] =
            (1.0 - ratio) / row.direct_scale;
        evaluation.jacobian[index] =
            -ratio * result.derivative * payload.parameter_scale
            / row.direct_scale;
    }
}

Evaluation make_evaluation(const Payload& payload) {
    const std::size_t row_count = payload.training_rows.size();
    const std::size_t variables = variable_count(payload);
    const std::size_t residuals = residual_count(payload);
    return Evaluation{
        std::vector<double>(residuals),
        std::vector<double>(residuals * variables),
        std::vector<double>(
            row_count, std::numeric_limits<double>::quiet_NaN()
        ),
        std::vector<double>(
            row_count, std::numeric_limits<double>::quiet_NaN()
        ),
    };
}

void evaluate_problem(
    const epcsaft_native_sdk_v1& table,
    const Payload& payload,
    const double* variables,
    std::size_t variable_size,
    Evaluation& evaluation
) {
    const std::size_t row_count = payload.training_rows.size();
    const std::size_t variable_total = variable_count(payload);
    const std::size_t residuals = residual_count(payload);
    if (variable_size != variable_total) {
        throw std::invalid_argument(
            "solver variables do not match the training-row dimension"
        );
    }
    if (evaluation.residuals.size() != residuals
        || evaluation.jacobian.size() != residuals * variable_total) {
        throw std::logic_error("evaluation scratch dimensions are invalid");
    }
    std::fill(evaluation.jacobian.begin(), evaluation.jacobian.end(), 0.0);
    const double parameter =
        payload.parameter_origin + payload.parameter_scale * variables[0];
    if (!std::isfinite(parameter)
        || parameter < payload.parameter_lower_bound
        || parameter > payload.parameter_upper_bound) {
        throw std::invalid_argument(
            "active parameter is outside its declared bounds"
        );
    }
    if (direct_observation(payload)) {
        evaluate_direct_problem(table, payload, parameter, evaluation);
        if (!std::all_of(
                evaluation.residuals.cbegin(),
                evaluation.residuals.cend(),
                [](double value) { return std::isfinite(value); }
            )
            || !std::all_of(
                evaluation.jacobian.cbegin(),
                evaluation.jacobian.cend(),
                [](double value) { return std::isfinite(value); }
            )) {
            throw std::runtime_error(
                "assembled residual or Jacobian is nonfinite"
            );
        }
        return;
    }
    for (std::size_t row_index = 0; row_index < row_count; ++row_index) {
        const Row& row = payload.training_rows[row_index];
        const std::size_t liquid_column = 1 + 2 * row_index;
        const std::size_t vapor_column = liquid_column + 1;
        const double liquid_volume =
            row.liquid_volume_origin * std::exp(variables[liquid_column]);
        const double vapor_volume =
            row.vapor_volume_origin * std::exp(variables[vapor_column]);
        if (!std::isfinite(liquid_volume) || !std::isfinite(vapor_volume)
            || liquid_volume < row.liquid_volume_bounds[0]
            || liquid_volume > row.liquid_volume_bounds[1]
            || vapor_volume < row.vapor_volume_bounds[0]
            || vapor_volume > row.vapor_volume_bounds[1]
            || liquid_volume >= vapor_volume) {
            throw std::invalid_argument(
                "phase volumes violate their declared bounds or topology"
            );
        }
        const Phase liquid = evaluate_phase(
            table, payload, row, row.liquid_first, liquid_volume, parameter
        );
        const Phase vapor = evaluate_phase(
            table, payload, row, row.vapor_first, vapor_volume, parameter
        );
        const bool pure = row.kind == ObservationKind::pure_phase;
        const std::size_t phase_volume_coordinate = pure ? 1 : 2;
        if (!(liquid.hessian[
                  phase_volume_coordinate * liquid.coordinate_count
                  + phase_volume_coordinate
              ] > 0.0)
            || !(vapor.hessian[
                  phase_volume_coordinate * vapor.coordinate_count
                  + phase_volume_coordinate
              ] > 0.0)) {
            throw std::runtime_error(
                "phase state is not mechanically stable"
            );
        }
        const std::size_t residual_offset = 4 * row_index;
        evaluation.residuals[residual_offset] =
            (liquid.pressure - row.pressure) / row.pressure_scale;
        evaluation.residuals[residual_offset + 1] =
            (vapor.pressure - row.pressure) / row.pressure_scale;
        evaluation.residuals[residual_offset + 2] =
            (liquid.gradient[0] - vapor.gradient[0])
            / row.chemical_potential_scales[0];
        evaluation.residuals[residual_offset + 3] = pure
            ? (row.molar_mass / liquid_volume - row.liquid_density)
                / row.liquid_density_scale
            : (liquid.gradient[1] - vapor.gradient[1])
                / row.chemical_potential_scales[1];

        auto jacobian = [&](std::size_t residual, std::size_t column)
            -> double& {
            return evaluation.jacobian[
                residual * variable_total + column
            ];
        };
        const std::size_t coordinate_count = liquid.coordinate_count;
        const std::size_t volume_coordinate = pure ? 1 : 2;
        const std::size_t parameter_coordinate = pure ? 2 : 3;
        jacobian(residual_offset, 0) =
            -gas_constant * row.temperature
            * liquid.hessian[
                volume_coordinate * coordinate_count + parameter_coordinate
            ]
            * payload.parameter_scale / row.pressure_scale;
        jacobian(residual_offset, liquid_column) =
            -gas_constant * row.temperature
            * liquid.hessian[
                volume_coordinate * coordinate_count + volume_coordinate
            ]
            * liquid_volume / row.pressure_scale;
        jacobian(residual_offset + 1, 0) =
            -gas_constant * row.temperature
            * vapor.hessian[
                volume_coordinate * coordinate_count + parameter_coordinate
            ]
            * payload.parameter_scale / row.pressure_scale;
        jacobian(residual_offset + 1, vapor_column) =
            -gas_constant * row.temperature
            * vapor.hessian[
                volume_coordinate * coordinate_count + volume_coordinate
            ]
            * vapor_volume / row.pressure_scale;
        const std::size_t chemical_residual_count = pure ? 1 : 2;
        for (
            std::size_t component = 0;
            component < chemical_residual_count;
            ++component
        ) {
            const std::size_t residual = residual_offset + 2 + component;
            const double scale = row.chemical_potential_scales[component];
            jacobian(residual, 0) =
                (liquid.hessian[
                    component * coordinate_count + parameter_coordinate
                 ] - vapor.hessian[
                    component * coordinate_count + parameter_coordinate
                 ])
                * payload.parameter_scale / scale;
            jacobian(residual, liquid_column) =
                liquid.hessian[
                    component * coordinate_count + volume_coordinate
                ] * liquid_volume / scale;
            jacobian(residual, vapor_column) =
                -vapor.hessian[
                    component * coordinate_count + volume_coordinate
                ] * vapor_volume / scale;
        }
        if (pure) {
            jacobian(residual_offset + 3, liquid_column) =
                -(row.molar_mass / liquid_volume)
                / row.liquid_density_scale;
        }
    }
    if (!std::all_of(
            evaluation.residuals.cbegin(),
            evaluation.residuals.cend(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            evaluation.jacobian.cbegin(),
            evaluation.jacobian.cend(),
            [](double value) { return std::isfinite(value); }
        )) {
        throw std::runtime_error(
            "assembled residual or Jacobian is nonfinite"
        );
    }
}

class GeneralCost final : public ceres::CostFunction {
public:
    GeneralCost(
        const epcsaft_native_sdk_v1* table, const Payload& payload
    ) : table_(table), payload_(payload), scratch_(make_evaluation(payload)) {
        set_num_residuals(static_cast<int>(residual_count(payload)));
        mutable_parameter_block_sizes()->push_back(
            static_cast<int>(variable_count(payload))
        );
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            evaluate_problem(
                *table_,
                payload_,
                values[0],
                static_cast<std::size_t>(parameter_block_sizes()[0]),
                scratch_
            );
            std::copy(
                scratch_.residuals.begin(),
                scratch_.residuals.end(),
                residuals
            );
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(
                    scratch_.jacobian.begin(),
                    scratch_.jacobian.end(),
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
    mutable Evaluation scratch_;
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
    if (variable_count == 1) {
        outcome.projected_parameter_jacobian = outcome.full_jacobian;
        return;
    }
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
    outcome.variables.resize(variable_count(payload));
    outcome.variables[0] =
        (physical_start - payload.parameter_origin) / payload.parameter_scale;
    if (!direct_observation(payload)) {
        for (
            std::size_t index = 0;
            index < payload.training_rows.size();
            ++index
        ) {
            const Row& row = payload.training_rows[index];
            outcome.variables[1 + 2 * index] =
                std::log(
                    row.liquid_volume_start / row.liquid_volume_origin
                );
            outcome.variables[2 + 2 * index] =
                std::log(
                    row.vapor_volume_start / row.vapor_volume_origin
                );
        }
    }
    GeneralCost cost(&table, payload);
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
    if (!direct_observation(payload)) {
        for (
            std::size_t index = 0;
            index < payload.training_rows.size();
            ++index
        ) {
            const Row& row = payload.training_rows[index];
            problem.SetParameterLowerBound(
                outcome.variables.data(),
                static_cast<int>(1 + 2 * index),
                std::log(
                    row.liquid_volume_bounds[0]
                    / row.liquid_volume_origin
                )
            );
            problem.SetParameterUpperBound(
                outcome.variables.data(),
                static_cast<int>(1 + 2 * index),
                std::log(
                    row.liquid_volume_bounds[1]
                    / row.liquid_volume_origin
                )
            );
            problem.SetParameterLowerBound(
                outcome.variables.data(),
                static_cast<int>(2 + 2 * index),
                std::log(
                    row.vapor_volume_bounds[0]
                    / row.vapor_volume_origin
                )
            );
            problem.SetParameterUpperBound(
                outcome.variables.data(),
                static_cast<int>(2 + 2 * index),
                std::log(
                    row.vapor_volume_bounds[1]
                    / row.vapor_volume_origin
                )
            );
        }
    }
    ceres::Solve(
        solver_options(payload), &problem, &outcome.summary
    );
    outcome.failure_reason = cost.failure_reason();
    try {
        outcome.evaluation = make_evaluation(payload);
        evaluate_problem(
            table,
            payload,
            outcome.variables.data(),
            outcome.variables.size(),
            outcome.evaluation
        );
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
    ) : table_(table),
        payload_(payload),
        parameter_(parameter_solver_value),
        scratch_(make_evaluation(payload)) {
        set_num_residuals(4);
        mutable_parameter_block_sizes()->push_back(2);
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            const std::array<double, 3> variables = {
                parameter_, values[0][0], values[0][1]
            };
            evaluate_problem(
                *table_,
                payload_,
                variables.data(),
                variables.size(),
                scratch_
            );
            std::copy(
                scratch_.residuals.begin(),
                scratch_.residuals.end(),
                residuals
            );
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                for (std::size_t row = 0; row < 4; ++row) {
                    jacobians[0][2 * row] =
                        scratch_.jacobian[3 * row + 1];
                    jacobians[0][2 * row + 1] =
                        scratch_.jacobian[3 * row + 2];
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
    mutable Evaluation scratch_;
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
        Evaluation evaluation = make_evaluation(row_payload);
        const std::array<double, 3> final_variables = {
            parameter_solver_value, variables[0], variables[1]
        };
        evaluate_problem(
            table,
            row_payload,
            final_variables.data(),
            final_variables.size(),
            evaluation
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
    const std::string& failure_reason,
    double modeled_value,
    double provider_derivative
) {
    PyObject* residual_tuple = doubles_to_tuple(residuals);
    if (residual_tuple == nullptr) return nullptr;
    return Py_BuildValue(
        "(ssddNOsdd)",
        row.row_id.c_str(),
        row.partition.c_str(),
        liquid_volume,
        vapor_volume,
        residual_tuple,
        usable ? Py_True : Py_False,
        failure_reason.c_str(),
        modeled_value,
        provider_derivative
    );
}

bool complete_evaluation(
    const SolveOutcome& outcome, const Payload& payload
) {
    const std::size_t residuals = residual_count(payload);
    return outcome.failure_reason.empty()
        && outcome.evaluation.residuals.size() == residuals
        && outcome.evaluation.jacobian.size()
            == residuals * outcome.variables.size();
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
    if (direct_observation(payload)) {
        for (
            std::size_t index = 0;
            index < payload.training_rows.size();
            ++index
        ) {
            const Row& row = payload.training_rows[index];
            const std::vector<double> residuals = {
                primary_evaluation_available
                    ? primary.evaluation.residuals[index]
                    : std::numeric_limits<double>::quiet_NaN()
            };
            PyObject* item = row_to_python(
                row,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                residuals,
                primary_evaluation_available,
                primary_evaluation_available ? "" : unavailable_reason,
                primary_evaluation_available
                    ? primary.evaluation.modeled_values[index]
                    : std::numeric_limits<double>::quiet_NaN(),
                primary_evaluation_available
                    ? primary.evaluation.provider_derivatives[index]
                    : std::numeric_limits<double>::quiet_NaN()
            );
            if (item == nullptr) {
                Py_DECREF(result);
                return nullptr;
            }
            PyTuple_SET_ITEM(
                result, static_cast<Py_ssize_t>(index), item
            );
        }
        for (
            std::size_t index = 0;
            index < payload.reporting_rows.size();
            ++index
        ) {
            const Row& row = payload.reporting_rows[index];
            Payload row_payload = payload;
            row_payload.training_rows = {row};
            row_payload.reporting_rows.clear();
            Evaluation evaluation = make_evaluation(row_payload);
            bool usable = primary_evaluation_available;
            std::string failure_reason =
                usable ? "" : unavailable_reason;
            if (usable) {
                try {
                    evaluate_problem(
                        table,
                        row_payload,
                        primary.variables.data(),
                        primary.variables.size(),
                        evaluation
                    );
                } catch (const std::exception& error) {
                    usable = false;
                    failure_reason = error.what();
                }
            }
            const std::vector<double> residuals = {
                usable
                    ? evaluation.residuals[0]
                    : std::numeric_limits<double>::quiet_NaN()
            };
            PyObject* item = row_to_python(
                row,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                residuals,
                usable,
                failure_reason,
                usable
                    ? evaluation.modeled_values[0]
                    : std::numeric_limits<double>::quiet_NaN(),
                usable
                    ? evaluation.provider_derivatives[0]
                    : std::numeric_limits<double>::quiet_NaN()
            );
            if (item == nullptr) {
                Py_DECREF(result);
                return nullptr;
            }
            PyTuple_SET_ITEM(
                result,
                static_cast<Py_ssize_t>(
                    payload.training_rows.size() + index
                ),
                item
            );
        }
        return result;
    }
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
            primary_evaluation_available ? "" : unavailable_reason,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
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
            outcome.failure_reason,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
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
            const auto& descriptor = table->capabilities[index];
            PyObject* item = (
                descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1
                || descriptor.capability
                    == EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1
            )
                ? descriptor_to_python(descriptor)
                : unsupported_descriptor_to_python(descriptor);
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

PyObject* evaluate_general_python(
    PyObject* capsule, PyObject* payload_object, PyObject* variables_object
) {
    try {
        const auto* table = capability_table(capsule);
        const Payload payload = parse_payload(payload_object);
        checked_descriptor(*table, payload);
        const std::vector<double> variables = doubles(
            variables_object, "solver variables"
        );
        Evaluation evaluation = make_evaluation(payload);
        evaluate_problem(
            *table,
            payload,
            variables.data(),
            variables.size(),
            evaluation
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

PyObject* solve_general_python(
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
                    payload.function_tolerance,
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
