#include "evaluator_fit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {

struct AnalyticContext final {
    std::size_t value_calls{0};
    std::size_t jacobian_calls{0};
};

template <std::size_t Capacity>
void copy_text(char (&target)[Capacity], const char* value) {
    const std::size_t length = std::strlen(value);
    if (length >= Capacity) {
        throw std::runtime_error("analytic evaluator text is too long");
    }
    std::fill(target, target + Capacity, '\0');
    std::copy(value, value + length, target);
}

int evaluate_analytic(
    void* opaque,
    const double* parameters,
    std::size_t parameter_count,
    std::int32_t request_mode,
    epcsaft_regression_evaluator_result_v1* result
) {
    if (opaque == nullptr || parameters == nullptr || parameter_count != 2
        || result == nullptr || result->struct_size != sizeof(*result)
        || result->row_count != 2 || result->parameter_count != 2
        || result->value_capacity < 2 || result->values == nullptr
        || result->row_result_capacity < 2 || result->row_results == nullptr) {
        return EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
    }
    const bool jacobian =
        request_mode
        == EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1;
    if ((!jacobian
         && request_mode
             != EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_ONLY_V1)
        || (jacobian
            && (result->jacobian == nullptr
                || result->jacobian_capacity < 4))) {
        return EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
    }
    auto& context = *static_cast<AnalyticContext*>(opaque);
    jacobian ? ++context.jacobian_calls : ++context.value_calls;
    const double p0 = parameters[0];
    const double p1 = parameters[1];
    result->values[0] = std::exp(p0 + 0.5 * p1);
    result->values[1] = std::exp(0.5 * p0 - p1);
    if (jacobian) {
        result->jacobian[0] = result->values[0];
        result->jacobian[1] = 0.5 * result->values[0];
        result->jacobian[2] = 0.5 * result->values[1];
        result->jacobian[3] = -result->values[1];
    }
    for (std::size_t index = 0; index < 2; ++index) {
        auto& row = result->row_results[index];
        if (row.struct_size != sizeof(row)) {
            return EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
        }
        row.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_OK_V1;
        row.reason[0] = '\0';
        copy_text(row.solver_status, "CONVERGED");
        copy_text(row.numerical_status, "CONDITIONED");
        copy_text(row.physical_status, "ANALYTIC_FIXTURE");
        copy_text(row.derivative_status, "EXACT_ANALYTIC");
        copy_text(row.chart_topology, "fixed");
        copy_text(
            row.provider_topology_fingerprint,
            "sha256:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
        );
        row.kkt_dimension = 2;
        row.kkt_rank = 2;
        row.kkt_condition_number_inf = 10.0;
    }
    copy_text(
        result->provider_parameter_fingerprint,
        "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    );
    copy_text(
        result->artifact_identity,
        "sha256:cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
    );
    result->error[0] = '\0';
    result->request_mode = request_mode;
    result->status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_OK_V1;
    return result->status;
}

}  // namespace

void test_evaluator_fit() {
    constexpr const char* rows[] = {"log-row", "identity-row"};
    constexpr const char* states[] = {"state", "state"};
    constexpr const char* state_schemas[] = {"fixed-state-v1", "fixed-state-v1"};
    constexpr const char* sources[] = {"source", "source"};
    constexpr const char* primitives[] = {"primitive-a", "primitive-b"};
    constexpr const char* primitive_units[] = {"1", "1"};
    constexpr const char* transforms[] = {"natural_log", "identity"};
    constexpr const char* references[] = {"reference", "reference"};
    constexpr const char* reference_fingerprints[] = {
        "sha256:1111111111111111111111111111111111111111111111111111111111111111",
        "sha256:2222222222222222222222222222222222222222222222222222222222222222",
    };
    constexpr const char* parameter_ids[] = {"parameter-a", "parameter-b"};
    constexpr const char* parameter_units[] = {"1", "1"};
    AnalyticContext context;
    epcsaft_regression_evaluator_sdk_v1 sdk{};
    sdk.abi_version = EPCSAFT_REGRESSION_EVALUATOR_V1_ABI_VERSION;
    sdk.table_size = sizeof(sdk);
    sdk.model_context = &context;
    sdk.row_count = 2;
    sdk.parameter_count = 2;
    sdk.row_ids = rows;
    sdk.state_ids = states;
    sdk.state_schema_ids = state_schemas;
    sdk.observation_source_ids = sources;
    sdk.primitive_ids = primitives;
    sdk.primitive_units = primitive_units;
    sdk.transform_ids = transforms;
    sdk.reference_ids = references;
    sdk.reference_fingerprints = reference_fingerprints;
    sdk.parameter_ids = parameter_ids;
    sdk.parameter_units = parameter_units;
    sdk.evaluator_identity = "analytic.evaluator.v1";
    sdk.capability_id = "positive-scalars-v1";
    sdk.capability_fingerprint =
        "sha256:3333333333333333333333333333333333333333333333333333333333333333";
    sdk.provider_artifact_identity = "provider-artifact";
    sdk.owner_artifact_identity = "owner-artifact";
    sdk.contract_fingerprint =
        "sha256:8888888888888888888888888888888888888888888888888888888888888888";
    sdk.model_fingerprint =
        "sha256:9999999999999999999999999999999999999999999999999999999999999999";
    sdk.provider_parameter_fingerprint =
        "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    sdk.expected_provider_topology_fingerprint =
        "sha256:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
    sdk.provider_sdk_capsule_name = "epcsaft.native_sdk.v1";
    sdk.provider_sdk_abi_version = 1;
    sdk.provider_sdk_table_size = 1;
    sdk.provider_sdk_result_size = 1;
    sdk.provider_sdk_mixture_result_size = 1;
    sdk.provider_sdk_neutral_reference_result_size = 1;
    sdk.provider_sdk_neutral_reference_derivative_result_size = 1;
    sdk.provider_sdk_reacting_phase_parameter_result_size = 1;
    sdk.artifact_identity =
        "sha256:cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc";
    sdk.single_thread_non_reentrant = 1;
    sdk.value_only_avoids_derivative_work = 1;
    sdk.result_size = sizeof(epcsaft_regression_evaluator_result_v1);
    sdk.row_result_size =
        sizeof(epcsaft_regression_evaluator_row_result_v1);
    sdk.evaluate = evaluate_analytic;

    epcsaft_regression::evaluator::Problem problem{
        {
            sdk.evaluator_identity,
            sdk.capability_id,
            sdk.capability_fingerprint,
            sdk.provider_artifact_identity,
            sdk.owner_artifact_identity,
            sdk.contract_fingerprint,
            sdk.model_fingerprint,
            sdk.provider_parameter_fingerprint,
            sdk.expected_provider_topology_fingerprint,
            sdk.artifact_identity,
        },
        {
            {"parameter-a", "1", 0.0, 1.0, -4.0, 4.0},
            {"parameter-b", "1", 0.0, 1.0, -4.0, 4.0},
        },
        {{0.0, 0.0}, {2.0, -1.0}},
        {
            {
                "log-row",
                "training",
                "state",
                "fixed-state-v1",
                "source",
                "primitive-a",
                "1",
                "natural_log",
                "reference",
                reference_fingerprints[0],
                std::exp(2.0),
                1.0,
            },
            {
                "identity-row",
                "training",
                "state",
                "fixed-state-v1",
                "source",
                "primitive-b",
                "1",
                "identity",
                "reference",
                reference_fingerprints[1],
                std::exp(-1.5),
                1.0,
            },
        },
        {100, 10.0, 1.0e-14, 1.0e-14, 1.0e-14},
        1.0e6,
        1.0e-8,
        1.0e-8,
    };
    const auto fit = epcsaft_regression::evaluator::solve(sdk, problem);
    if (fit.solves.size() != 2
        || !fit.solves.front().summary.IsSolutionUsable()
        || fit.solves.front().full_jacobian.rank != 2
        || fit.solves.front().projected_parameter_jacobian.rank != 2
        || std::abs(fit.solves.front().variables[0] - 1.0) > 1.0e-8
        || std::abs(fit.solves.front().variables[1] - 2.0) > 1.0e-8
        || !fit.confirmations_usable
        || fit.confirmation_parameter_delta > 1.0e-8
        || context.value_calls == 0 || context.jacobian_calls == 0) {
        throw std::runtime_error(
            "exact analytic evaluator did not satisfy the fit contract"
        );
    }
    const char* reversed_ids[] = {"parameter-b", "parameter-a"};
    auto invalid = sdk;
    invalid.parameter_ids = reversed_ids;
    bool rejected = false;
    try {
        epcsaft_regression::evaluator::validate_contract(invalid, problem);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    if (!rejected) {
        throw std::runtime_error(
            "reordered evaluator parameter columns were not rejected"
        );
    }
}
