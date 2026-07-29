# BoundGSEnergy

Certified ground-state energy bounds for quantum spin systems via structured NPA and renormalization-group coarse-graining.

This repository is team **BoundGSEnergy**'s entry for [QuantumBFS/quantum.harness#49](https://github.com/QuantumBFS/quantum.harness/issues/49), registered as [PR #219](https://github.com/QuantumBFS/quantum.harness/pull/219).

## Current result

The challenge target is a certified energy-density gap of at most `1e-5` for the periodic spin-1/2 Heisenberg chain at `N = 200`. The target system size has not yet been reached.

At `N = 20`, the retained full structured-NPA certificate gives

```text
certified lower bound  = -0.44522841380294176
certified upper bound  = -0.44521932649373910
certified gap          =  9.087309202669595e-6
```

This passes the `1e-5` accuracy gate. The independent exact-compatible RG/RDM8 route gives a certified `N = 20` gap of `6.967350626092331e-5` and demonstrates that the RDM8 lift improves the compressed certificate by `3.772029497473017e-5`.

The complete scientific status, retained evidence, discarded directions, compute scaling, and next decision are in [docs/PROGRESS.md](docs/PROGRESS.md).

## Repository layout

- `scripts/heisenberg1d/` contains the retained 1D Heisenberg implementations, entrypoints, and solver-free regression tests.
- `results/` contains generated run artifacts and is ignored by Git. New runs must use a fresh dated folder.
- `artifacts/coarse_grainers/` contains generated coarse-graining caches and is ignored by Git.
- `docs/` contains the ratified plan and the mathematical contracts behind the finite-size, hybrid-RDM, dual-replay, and full-SOHS claims.
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
```

Recompute the accepted full structured-NPA certificate:

```bash
run_stamp=$(date +%Y%m%d-%H%M%S)
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/h11_full_sohs_certificate.jl --run-name="${run_stamp}-full-sohs-n20" --mosek-threads=16
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl --result="results/${run_stamp}-full-sohs-n20/full_sohs_certificate.json"
```

The command builds its exact finite-energy upper certificate, solves the structured-NPA SDP, writes the SOHS transcript, and immediately replays both endpoints. Generated results are not part of the commit.

The longer RG/RDM8 reproduction chain, including every generated predecessor, is in [scripts/heisenberg1d/README.md](scripts/heisenberg1d/README.md). The SDP runs require a Mosek license. A clean `OPTIMAL` status is necessary but not sufficient for a certified claim: the validators also replay the exact arithmetic, verify the model or RG map provenance, and check the rigorous sandwich against an exact-rational upper bound.

See [scripts/heisenberg1d/README.md](scripts/heisenberg1d/README.md) for the full script map and finite-size ladder.
