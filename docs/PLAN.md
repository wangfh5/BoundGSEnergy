# PLAN — BoundGSEnergy / challenge #49

Ratified 2026-07-27. These four answers are the skeleton of the final report.

The staged execution plan and current decision gates are maintained in [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

## 1. What we compute

Implement the NPA + coarse-graining combination: block neighbouring spins
into super-sites following PRX 14, 021008 (Sec. III-D-2), feed the blocked
model into the structured-SDP certification pipeline already working from the
Track-2 reproduction (QMBCertify.jl + Mosek). Start at ladder rung 1:
certified lower bounds for the 1D Heisenberg chain.

## 2. Sizes

Ladder: N = 20 (align with known values) → 50 → 100 → 200 (the #49 target).
Each rung must validate (precision and cost) before the next is attempted.

The existing compressed-LTI MVE bounds the energy density of an infinite translation-invariant chain. Its `2n` physical-site support is not a finite `N = 2n` system, so it does not by itself complete this ladder. The detailed research plan treats a finite periodic chain as the default interpretation of the challenge target unless the challenge owner confirms a hierarchy-support interpretation.

## 3. Success criteria (tiered)

- **Floor (teaching success):** coarse-grained SDP reproduces known bounds at
  N = 20–50 AND is visibly cheaper (SDP size / wall time) than the
  non-coarse-grained version at the same accuracy — proves coarse-graining
  actually buys something.
- **Target:** N = 200 lower bound meeting the #49 criterion of 10⁻⁵ accuracy.
- **Stretch:** ladder rung 2 (1D J₁–J₂, 100 spins @ 10⁻³) if time allows.

## 4. Kill criteria — stop and report honestly if either fires

- Coarse-graining shows no size/cost advantage at N = 50 (implementation bug
  or no real benefit — both are reportable findings).
- A single SDP point exceeds local memory (54 GB) — would need cluster /
  redesign, outside this time box.
