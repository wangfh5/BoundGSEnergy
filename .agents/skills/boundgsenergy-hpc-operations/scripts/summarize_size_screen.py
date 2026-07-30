#!/usr/bin/env python3

import json
import sys
from pathlib import Path
from typing import List, Tuple


REQUIRED_GATES = (
    "all_three_optimal",
    "raw_monotonicity_passed",
    "fusion_dominates_rg_only",
    "all_raw_objectives_below_finite_reference",
    "local_rdm_exactness_passed",
    "physical_relaxation_compatibility_proven",
)


def load_summary(path: Path) -> Tuple[List[str], bool]:
    with path.open(encoding="utf-8") as handle:
        result = json.load(handle)

    meta = result["meta"]
    rg_only = result["rg_only"]
    baseline = result["baseline"]
    hybrid = result["hybrid"]
    comparison = result["comparison"]
    gates = result["gates"]
    structural_ok = all(gates.get(gate) is True for gate in REQUIRED_GATES)
    row = [
        str(path),
        str(meta["system_size"]),
        str(meta["n_super"]),
        str(meta.get("formal_certificate", False)).lower(),
        f"{rg_only['energy']:.12g}",
        f"{baseline['energy']:.12g}",
        f"{hybrid['energy']:.12g}",
        f"{comparison['hybrid_over_rg_only']:.12g}",
        f"{comparison['raw_improvement']:.12g}",
        f"{comparison['strict_improvement']:.12g}",
        str(gates.get("strict_improvement_passed", False)).lower(),
        str(structural_ok).lower(),
    ]
    return row, structural_ok


def main(paths: List[str]) -> int:
    if not paths:
        print("usage: summarize_size_screen.py <size_screen.json> ...", file=sys.stderr)
        return 2

    print(
        "\t".join(
            (
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
            )
        )
    )
    all_structural_ok = True
    for value in paths:
        try:
            row, structural_ok = load_summary(Path(value))
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            print(f"{value}: {error}", file=sys.stderr)
            return 2
        print("\t".join(row))
        all_structural_ok &= structural_ok

    return 0 if all_structural_ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
