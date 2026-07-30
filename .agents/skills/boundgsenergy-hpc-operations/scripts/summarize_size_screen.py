#!/usr/bin/env python3

import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


REQUIRED_GATES = (
    "all_three_optimal",
    "raw_monotonicity_passed",
    "fusion_dominates_rg_only",
    "all_raw_objectives_below_finite_reference",
    "local_rdm_exactness_passed",
    "physical_relaxation_compatibility_proven",
)


def format_number(value: Any) -> str:
    return f"{value:.12g}" if isinstance(value, (int, float)) else ""


def replay_correction(result: Dict[str, Any], route: str) -> str:
    replay = result.get(route)
    if not isinstance(replay, dict):
        return ""
    correction = replay.get("total_replay_correction")
    if correction is None:
        raw = replay.get("raw_objective")
        strict = replay.get("strict_lower_endpoint")
        if isinstance(raw, (int, float)) and isinstance(strict, (int, float)):
            correction = strict - raw
    return format_number(correction)


def replay_endpoint(result: Dict[str, Any], route: str) -> str:
    replay = result.get(route)
    return (
        format_number(replay.get("strict_lower_endpoint"))
        if isinstance(replay, dict)
        else ""
    )


def exact_replay_endpoint(result: Dict[str, Any], route: str) -> Optional[Fraction]:
    replay = result.get(route)
    if not isinstance(replay, dict):
        return None
    rational = replay.get("strict_lower_endpoint_rational")
    if isinstance(rational, str):
        numerator, denominator = rational.split("//", maxsplit=1)
        return Fraction(int(numerator), int(denominator))
    endpoint = replay.get("strict_lower_endpoint")
    return Fraction.from_float(endpoint) if isinstance(endpoint, (int, float)) else None


def shifted_blocks(result: Dict[str, Any], route: str) -> str:
    replay = result.get(route)
    if not isinstance(replay, dict):
        return ""
    repairs = replay.get("psd_repair_diagnostics")
    return str(repairs.get("shifted_block_count", "")) if isinstance(repairs, dict) else ""


def load_summary(
    path: Path,
) -> Tuple[List[str], bool, Optional[Fraction], Optional[Fraction]]:
    with path.open(encoding="utf-8") as handle:
        result = json.load(handle)

    meta = result["meta"]
    rg_only = result.get("rg_only") or {}
    baseline = result.get("baseline") or {}
    hybrid = result.get("hybrid") or {}
    comparison = result.get("comparison") or {}
    gates = result.get("gates") or {}
    required_gates = REQUIRED_GATES
    if "omega_normalization" in meta:
        required_gates += ("omega_normalization_exact",)
    structural_ok = all(gates.get(gate) is True for gate in required_gates)
    row = [
        str(path),
        str(meta["system_size"]),
        str(meta["n_super"]),
        str(meta.get("formal_certificate", False)).lower(),
        format_number(rg_only.get("energy")),
        format_number(baseline.get("energy")),
        format_number(hybrid.get("energy")),
        format_number(comparison.get("hybrid_over_rg_only")),
        format_number(comparison.get("raw_improvement")),
        format_number(comparison.get("strict_improvement")),
        str(gates.get("strict_improvement_passed", False)).lower(),
        str(structural_ok).lower(),
        str(meta.get("coefficient_scale", "")),
        str(meta.get("omega_normalization", "")),
        str(meta.get("solve_form", "")),
        str((result.get("failure") or {}).get("phase", "")),
        str(rg_only.get("status", "")),
        str(baseline.get("status", "")),
        str(hybrid.get("status", "")),
        replay_endpoint(result, "baseline_replay"),
        replay_endpoint(result, "hybrid_replay"),
        replay_correction(result, "baseline_replay"),
        replay_correction(result, "hybrid_replay"),
        shifted_blocks(result, "baseline_replay"),
        shifted_blocks(result, "hybrid_replay"),
    ]
    return (
        row,
        structural_ok,
        exact_replay_endpoint(result, "baseline_replay"),
        exact_replay_endpoint(result, "hybrid_replay"),
    )


def main(paths: List[str]) -> int:
    if not paths:
        print("usage: summarize_size_screen.py <size_screen.json> ...", file=sys.stderr)
        return 2

    header = (
        "path",
        "N",
        "n_super",
        "formal_certificate",
        "rg_only",
        "structured_npa",
        "fusion",
        "fusion_over_rg",
        "fusion_over_structured",
        "strict_fusion_over_structured",
        "strict_gain",
        "structural_gates",
        "coefficient_scale",
        "omega_normalization",
        "solve_form",
        "failure_phase",
        "rg_status",
        "baseline_status",
        "hybrid_status",
        "baseline_strict_endpoint",
        "hybrid_strict_endpoint",
        "baseline_replay_correction",
        "hybrid_replay_correction",
        "baseline_psd_shifts",
        "hybrid_psd_shifts",
        "robust_strict_gain",
        "robust_strict_gain_passed",
    )
    summaries = []
    all_structural_ok = True
    for value in paths:
        try:
            summary = load_summary(Path(value))
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            print(f"{value}: {error}", file=sys.stderr)
            return 2
        summaries.append(summary)
        all_structural_ok &= summary[1]

    best_baselines: Dict[Tuple[str, str], Fraction] = {}
    for row, _, baseline_endpoint, _ in summaries:
        if baseline_endpoint is not None:
            key = row[1], row[2]
            best_baselines[key] = max(
                best_baselines.get(key, baseline_endpoint),
                baseline_endpoint,
            )
    rows = []
    for row, _, _, hybrid_endpoint in summaries:
        key = row[1], row[2]
        robust_gain = (
            hybrid_endpoint - best_baselines[key]
            if hybrid_endpoint is not None and key in best_baselines
            else None
        )
        row.append(format_number(float(robust_gain)) if robust_gain is not None else "")
        row.append(str(robust_gain is not None and robust_gain > 0).lower())
        rows.append(row)

    print("\t".join(header))
    for row in rows:
        print("\t".join(row))
    return 0 if all_structural_ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
