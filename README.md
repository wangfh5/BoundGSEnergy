# BoundGSEnergy

Certified ground-state energy bounds for quantum spin systems via structured NPA and renormalization-group coarse-graining.

This repository is team **BoundGSEnergy**'s entry for [QuantumBFS/quantum.harness#49](https://github.com/QuantumBFS/quantum.harness/issues/49), registered as [PR #219](https://github.com/QuantumBFS/quantum.harness/pull/219).

## Current result

The challenge target is a certified energy-density gap of at most `1e-5` for the periodic spin-1/2 Heisenberg chain at `N = 200`. The target system size has not yet been reached.

At `N = 20`, the unified structured-NPA + exact-compatible RG certificate gives

```text
certified lower bound  = -0.44522454084723984
certified upper bound  = -0.44521932649373910
certified gap          =  5.214353500759827e-6
```

This passes the `1e-5` accuracy gate and is strictly tighter than both matched parent relaxations. The matched structured-NPA endpoint is `-0.4452272566443788`, while the independent exact-compatible RG/RDM8 route gives `-0.44528862186638035`. The unified model improves the matched structured-NPA certificate by `2.715797138952425e-6`.

The accepted model uses order 4, `extra = 1`, RDM9, PSD state optimality level 3, linear state optimality, an audited eight-spin RDM link, and a 13-bit dyadic `D = 4`, `n_super = 5` RG hierarchy. Multiplying only the `omega4` link equations by the exact power of two `16` preserves the feasible set while reducing the strict replay loss enough to expose the RG tightening.

The concise [Challenge Report](docs/challenge-report/report.html) summarizes the current stage. The complete scientific status, retained evidence, compute scaling, and next decision are in [docs/PROGRESS.md](docs/PROGRESS.md).

## Repository layout

- `scripts/heisenberg1d/` contains the retained 1D Heisenberg implementations, entrypoints, and solver-free regression tests.
- `results/` contains generated run artifacts and is ignored by Git. New runs must use a fresh dated folder.
- `artifacts/coarse_grainers/` contains generated coarse-graining caches and is ignored by Git.
- `docs/` contains the Challenge Report, ratified plan, progress report, and mathematical contracts behind the certified claims.
- `notes/` is a pedagogical guide to the method and repository.
- `knowledge/` is the survey library; `knowledge/ref.bib` is the bibliography source of truth.

## Reproduce

Install the Juliaup `1.11.5` channel, then instantiate the pinned environment with that explicit runtime:

```bash
juliaup add 1.11.5
julia +1.11.5 --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run all retained solver-free checks:

```bash
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_metadata.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h4.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_rg_rdm8.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_structured_rg_hybrid.jl
```

Recompute the accepted full structured-NPA certificate:

```bash
run_stamp=$(date +%Y%m%d-%H%M%S)
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/h11_full_sohs_certificate.jl --run-name="${run_stamp}-full-sohs-n20" --mosek-threads=16
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl --result="results/${run_stamp}-full-sohs-n20/full_sohs_certificate.json"
```

The command builds its exact finite-energy upper certificate, solves the structured-NPA SDP, writes the SOHS transcript, and immediately replays both endpoints. Generated results are not part of the commit.

Recompute the accepted unified certificate and independently rebuild and replay it:

```bash
run_stamp=$(date +%Y%m%d-%H%M%S)
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/structured_rg_local.jl --run-name="${run_stamp}-structured-rg-n20" --link-mode=rdm8 --n-super=5 --mosek-threads=16 --feasibility-tolerance=1e-7
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_structured_rg_hybrid.jl --result="results/${run_stamp}-structured-rg-n20/structured_rg.json"
```

The longer RG/RDM8 reproduction chain, including every generated predecessor, is in [scripts/heisenberg1d/README.md](scripts/heisenberg1d/README.md). The SDP runs require a Mosek license. A clean `OPTIMAL` status is necessary but not sufficient for a certified claim: the validators also replay the exact arithmetic, verify the model or RG map provenance, and check the rigorous sandwich against an exact-rational upper bound.

See [scripts/heisenberg1d/README.md](scripts/heisenberg1d/README.md) for the full script map and finite-size ladder.

Regenerate the standalone Challenge Report after editing its source:

```bash
python3 skills/report/render_report.py docs/challenge-report
```
