from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path
import sys

from run_candidate import (
    _canonical_json_sha256,
    _canonical_receipt_bytes,
    _require_hash,
    _require_import_origin,
    _require_installed_distribution_matches_wheel,
    _sha256,
)


PROVIDER_COMMIT = "45d5764f61729d387100348a8ff91792f6e0a395"
PROVIDER_TREE = "271d4848faf73afd4cf0683efe5c855053df7d01"
PROVIDER_WHEEL_SHA256 = "95d2292b052ab74657931f2dec97c3ea4160d9b17812956515c7195a853e6c5b"
PROVIDER_HEADER_SHA256 = "6f3a186bf5359f32449a31c544e8b8525be6c594804151d3e74a74e411ded8f4"
SOURCE_COMMIT = "73a37f5935e919a34d1e4fa3af285951d6fac8e7"
SOURCE_SHA256 = "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f"
SOURCE_PACKET_SHA256 = "d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e"
PAIR_RECORD_SHA256 = "747e8281c7a1e4240ee4badbc0bedd047521fb303726699b38c99fccf7f74c2a"
MODEL_FINGERPRINT = "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
HEADER_MEMBER = "epcsaft/include/epcsaft/native_sdk_v1.h"
CORE_MEMBER = "epcsaft/_core.cpython-313-x86_64-linux-gnu.so"


def _load_rows(path: Path) -> list[tuple[str, float, float, float, float]]:
    _require_hash(path, SOURCE_SHA256, "May 2015 source CSV")
    with path.open(newline="", encoding="utf-8") as stream:
        source_rows = list(csv.DictReader(stream))
    rows = [
        (
            row["row_id"],
            float(row["T_K"]),
            float(row["P_Pa"]),
            float(row["x_methane"]),
            float(row["y_methane"]),
        )
        for row in source_rows
    ]
    if len(rows) != 17:
        raise SystemExit("the frozen preflight requires all 17 May 2015 rows")
    return rows


def _load_native(path: Path) -> object:
    if not path.is_file():
        raise SystemExit(f"preflight module not found: {path}")
    spec = importlib.util.spec_from_file_location("_binary_kij_preflight", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load preflight module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _phase_semantics(
    model: object,
    row: tuple[str, float, float, float, float],
    phase: tuple[object, ...],
) -> tuple[object, tuple[float, ...]]:
    from epcsaft import unit_registry

    volume, n_methane, n_ethane, callback_pressure, gradient, fingerprint = phase
    evaluated = model.evaluate(
        temperature=row[1] * unit_registry.kelvin,
        molar_density=(1.0 / volume) * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(n_methane, n_ethane),
    )
    pressure_error = (
        abs(float(evaluated.pressure.magnitude) - callback_pressure)
        / callback_pressure
    )
    return evaluated, (pressure_error, *gradient, volume, n_methane, n_ethane, fingerprint)


def _solve_record(start: float, values: tuple[object, ...]) -> dict[str, object]:
    (
        termination,
        usable,
        initial_cost,
        final_cost,
        iterations,
        fitted_kij,
        singular_values,
        full_rank,
        condition,
        rank_threshold,
        projected_singular,
        projected_rank,
        max_pressure_closure,
        min_stability,
        min_phase_separation,
        complete_jacobian,
        worst,
        failure,
    ) = values
    return {
        "start_kij": start,
        "termination": termination,
        "solution_usable": usable,
        "initial_cost": initial_cost,
        "final_cost": final_cost,
        "iterations": iterations,
        "fitted_kij": fitted_kij,
        "largest_scaled_jacobian_singular_value": singular_values[0],
        "smallest_scaled_jacobian_singular_value": singular_values[-1],
        "full_rank": full_rank,
        "condition_number": condition,
        "rank_threshold": rank_threshold,
        "projected_parameter_singular_value": projected_singular,
        "projected_parameter_rank": projected_rank,
        "max_pressure_relative_closure": max_pressure_closure,
        "minimum_volume_stability_slope": min_stability,
        "minimum_relative_phase_separation": min_phase_separation,
        "complete_finite_jacobian": complete_jacobian,
        "worst_pressure_row": {
            "row_id": worst[0],
            "liquid_pressure_relative_closure": worst[1],
            "vapor_pressure_relative_closure": worst[2],
            "methane_mu_over_rt_equality_residual": worst[3],
            "ethane_mu_over_rt_equality_residual": worst[4],
            "max_abs_mu_over_rt_equality_residual_all_rows": worst[5],
        },
        "failure": failure or None,
    }


def _symmetric_relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), sys.float_info.min)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Replay the frozen binary-kij numerical preflight."
    )
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--native-module", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    repository = Path(__file__).resolve().parents[1]
    source_path = repository / "evidence" / "may-2015-methane-ethane-vle.csv"
    rows = _load_rows(source_path)
    _require_hash(arguments.provider_wheel, PROVIDER_WHEEL_SHA256, "provider wheel")
    binding = _require_installed_distribution_matches_wheel(arguments.provider_wheel, "epcsaft")
    verified_paths = binding["verified_paths"]
    if not isinstance(verified_paths, dict):
        raise SystemExit("invalid installed-provider binding")
    header_path = Path(verified_paths[HEADER_MEMBER])
    _require_hash(header_path, PROVIDER_HEADER_SHA256, "installed provider public header")

    import epcsaft
    from epcsaft import EPCSAFT, ParameterBundle, native_sdk

    _require_import_origin(epcsaft.__file__, binding, "epcsaft/__init__.py")
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    if (
        parameters.component_ids != ("methane", "ethane")
        or parameters.fingerprint != MODEL_FINGERPRINT
    ):
        raise SystemExit("provider component order or parameter fingerprint changed")
    pair_records = [record for record in parameters.records if record.family == "k_ij"]
    if len(pair_records) != 1 or float(pair_records[0].value) != 0.0:
        raise SystemExit("the source kij=0 start is unavailable")
    model = EPCSAFT(parameters)
    if (
        model.component_ids != ("methane", "ethane")
        or model.parameter_fingerprint != MODEL_FINGERPRINT
    ):
        raise SystemExit("model-bound provider identity changed")

    native = _load_native(arguments.native_module.resolve())
    callback_raw, jacobian_raw, liquid_raw, vapor_raw, solve_raw = native.run(
        native_sdk(model), rows
    )
    callback = {
        "value_directional_scaled_error": callback_raw[0],
        "gradient_hessian_directional_scaled_error": callback_raw[1],
        "hessian_symmetry_scaled_error": callback_raw[2],
        "pressure_identity_scaled_error": callback_raw[3],
        "central_difference_step": 1.0e-5,
        "value_and_gradient_scale": "abs(error)/(5e-7 + 5e-5*abs(expected))",
        "symmetry_scale": "abs(error)/(1e-12 + 5e-13*abs(expected))",
        "pressure_identity_scale_pa": "abs(error)/(1e-6 + 5e-13*abs(expected))",
        "acceptance": "each scaled error <= 1",
    }
    jacobian = {
        "max_abs_directional_error": jacobian_raw[0],
        "max_scaled_directional_error": jacobian_raw[1],
        "residual_count": jacobian_raw[2],
        "variable_count": jacobian_raw[3],
        "jacobian_entry_count": jacobian_raw[4],
        "central_difference_step": 1.0e-6,
        "scale": "abs(error)/(2e-7 + 2e-6*abs(expected))",
        "acceptance": "max scaled directional error <= 1 and exact 68 x 35 shape",
    }
    liquid_public, liquid = _phase_semantics(model, rows[0], liquid_raw)
    vapor_public, vapor = _phase_semantics(model, rows[0], vapor_raw)
    chemical_semantic_errors = []
    for component in range(2):
        callback_difference = liquid[1 + component] - vapor[1 + component]
        public_difference = (
            liquid_public.residual_chemical_potential_over_rt[component]
            - vapor_public.residual_chemical_potential_over_rt[component]
            + math.log(
                (liquid[6 + component] / liquid[5])
                / (vapor[6 + component] / vapor[5])
            )
        )
        chemical_semantic_errors.append(abs(callback_difference - public_difference))
    semantics = {
        "component_order": ["methane", "ethane"],
        "representative_row_id": rows[0][0],
        "max_pressure_relative_difference_from_public_value_api": max(liquid[0], vapor[0]),
        "max_mu_over_rt_phase_difference_from_public_value_api": max(chemical_semantic_errors),
        "acceptance": "both maxima <= 1e-12",
    }

    starts = (0.0, -0.05, 0.05)
    solves = [_solve_record(start, values) for start, values in zip(starts, solve_raw, strict=True)]
    primary = solves[0]
    confirmation = []
    for solve in solves[1:]:
        confirmation.append(
            {
                "start_kij": solve["start_kij"],
                "scaled_kij_delta_from_primary": abs(
                    solve["fitted_kij"] - primary["fitted_kij"]
                ) / 0.01,
                "symmetric_relative_cost_delta_from_primary": _symmetric_relative_difference(
                    solve["final_cost"], primary["final_cost"]
                ),
            }
        )

    convergence = all(
        solve["termination"] == "CONVERGENCE"
        and solve["solution_usable"]
        and solve["complete_finite_jacobian"]
        and math.isfinite(solve["final_cost"])
        and solve["final_cost"] < solve["initial_cost"]
        for solve in solves
    )
    fitted_kij = primary["fitted_kij"]
    gates = {
        "installed_artifact_identity": True,
        "provider_value_gradient_hessian_contract": all(
            math.isfinite(value) and 0.0 <= value <= 1.0 for value in callback_raw
        ),
        "methane_ethane_semantic_coordinate_order": max(
            liquid[0], vapor[0], *chemical_semantic_errors
        ) <= 1.0e-12,
        "exact_scaled_residual_jacobian": math.isfinite(jacobian_raw[1])
        and 0.0 <= jacobian_raw[1] <= 1.0
        and jacobian_raw[2:] == (68, 35, 2380),
        "ceres_convergence_and_cost_reduction_all_starts": convergence,
        "full_scaled_jacobian_rank_35": primary["full_rank"] == 35,
        "projected_parameter_rank_1": primary["projected_parameter_rank"] == 1,
        "finite_moved_non_bound_kij": math.isfinite(fitted_kij)
        and abs(fitted_kij) > 1.0e-8
        and fitted_kij > -0.15 + 1.0e-8
        and fitted_kij < 0.10 - 1.0e-8,
        "positive_stability_and_separated_phases": primary[
            "minimum_volume_stability_slope"
        ]
        > 0.0
        and primary["minimum_relative_phase_separation"] > 1.0e-3,
        "perturbed_start_confirmation": max(
            item["scaled_kij_delta_from_primary"] for item in confirmation
        ) <= 1.0e-5
        and max(item["symmetric_relative_cost_delta_from_primary"] for item in confirmation)
        <= 1.0e-8,
        "scaled_pressure_closure_le_1e-8": primary["max_pressure_relative_closure"]
        <= 1.0e-8,
    }
    unexpected_failures = [
        name
        for name, passed in gates.items()
        if not passed and name != "scaled_pressure_closure_le_1e-8"
    ]
    if unexpected_failures:
        raise SystemExit(
            f"preflight failed outside the frozen pressure gate: {unexpected_failures}"
        )
    if gates["scaled_pressure_closure_le_1e-8"]:
        raise SystemExit("the expected frozen pressure-closure blocker was not reproduced")

    payload: dict[str, object] = {
        "receipt_id": "binary-kij-all-may2015-installed-artifact-preflight-v1",
        "status": "BLOCKED_FROZEN_PRESSURE_CLOSURE",
        "authority_effect": "none",
        "provider": {
            "commit": PROVIDER_COMMIT,
            "tree": PROVIDER_TREE,
            "wheel_sha256": PROVIDER_WHEEL_SHA256,
            "installed_public_header_sha256": PROVIDER_HEADER_SHA256,
            "installed_wheel_binding": {
                "distribution": binding["distribution"],
                "version": binding["version"],
                "verified_member_count": binding["verified_member_count"],
                "imported_core_sha256": _sha256(Path(verified_paths[CORE_MEMBER])),
            },
            "model_fingerprint": MODEL_FINGERPRINT,
            "component_order": ["methane", "ethane"],
            "callback_coordinates": ["n_methane_mol", "n_ethane_mol", "volume_m3", "kij"],
            "callback_derivative_order": 2,
            "callback_derivative_backend": "provider CppAD value, gradient, and Hessian",
        },
        "source": {
            "validation_commit": SOURCE_COMMIT,
            "csv_sha256": SOURCE_SHA256,
            "packet_sha256": SOURCE_PACKET_SHA256,
            "pair_record_sha256": PAIR_RECORD_SHA256,
            "row_count": 17,
            "row_ids": [row[0] for row in rows],
            "use": "all rows are training data",
        },
        "problem": {
            "variables": "[z_kij, u_L_1, u_V_1, ..., u_L_17, u_V_17]",
            "variable_count": 35,
            "kij_transform": "kij = 0 + 0.01*z_kij",
            "kij_start": 0.0,
            "kij_bounds": [-0.15, 0.10],
            "liquid_volume": "V_L = 6.5e-5 m3*exp(u_L)",
            "vapor_volume": "V_V = (R*T/P) m3*exp(u_V)",
            "liquid_volume_bounds_m3": [2.0e-5, 1.0e-4],
            "vapor_volume_bounds_m3": [1.0e-4, 1.0e-2],
            "residuals_per_row": [
                "0.25*(P_L-P_observed)/P_observed",
                "0.25*(P_V-P_observed)/P_observed",
                "0.25*((mu_methane/RT)_L-(mu_methane/RT)_V)",
                "0.25*((mu_ethane/RT)_L-(mu_ethane/RT)_V)",
            ],
            "residual_count": 68,
            "jacobian": "exact provider Hessian chain rule; no third derivatives or density roots",
            "ceres": {
                "linear_solver": "DENSE_QR",
                "threads": 1,
                "max_iterations": 500,
                "function_gradient_parameter_tolerances": [1.0e-10, 1.0e-10, 1.0e-10],
            },
            "starts": [
                {"kij": 0.0, "liquid_volume_multiplier": 1.0, "vapor_volume_multiplier": 1.0},
                {"kij": -0.05, "liquid_volume_multiplier": 1.01, "vapor_volume_multiplier": 0.98},
                {"kij": 0.05, "liquid_volume_multiplier": 1.01, "vapor_volume_multiplier": 0.98},
            ],
        },
        "implementation": {
            "installed": False,
            "public_api": False,
            "native_source_sha256": _sha256(repository / "tools" / "binary_kij_preflight.cpp"),
            "runner_sha256": _sha256(Path(__file__).resolve()),
            "cmake_sha256": _sha256(repository / "CMakeLists.txt"),
            "native_module_sha256": _sha256(arguments.native_module),
        },
        "checks": {
            "callback_derivatives": callback,
            "semantic_coordinate_order": semantics,
            "residual_jacobian": jacobian,
            "gates": gates,
        },
        "solves": solves,
        "confirmation": confirmation,
        "status_split": {
            "solver_convergence": "PASS",
            "numerical_derivatives_ranks_and_confirmation": "PASS",
            "physical_pressure_closure": "BLOCKED",
        },
        "blocking_gate": {
            "name": "scaled_pressure_closure_le_1e-8",
            "threshold": 1.0e-8,
            "observed": primary["max_pressure_relative_closure"],
            "worst_row": primary["worst_pressure_row"],
            "interpretation": (
                "the frozen fixed-composition lifted-volume fit cannot satisfy "
                "the approved physical closure gate"
            ),
        },
        "claim": "installed-artifact derivative and in-sample numerical preflight feasibility only",
        "predictive_status": "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
        "excluded": [
            "runtime or package API",
            "changed residual weights, bounds, or tolerances",
            "predicted compositions or held-out evidence",
            "uncertainty, covariance, or global identifiability",
            "density roots, equilibrium dependencies, or copied provider equations",
            "provider catalog persistence, promotion, publication, or authority change",
        ],
    }
    receipt = _canonical_receipt_bytes(payload)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_bytes(receipt)
    print(f"wrote {arguments.output}")
    print(f"receipt subject SHA-256: {_canonical_json_sha256(payload)}")
    print(
        "BLOCKED_FROZEN_PRESSURE_CLOSURE: "
        f"{primary['max_pressure_relative_closure']:.16g} > 1e-8 at "
        f"{primary['worst_pressure_row']['row_id']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
