# Progress report — challenge #49

Updated 2026-07-29. This is the authoritative status report for [QuantumBFS/quantum.harness#49](https://github.com/QuantumBFS/quantum.harness/issues/49) and PR #219.

## Executive status

The challenge asks for certified energy-density lower bounds for increasingly large spin-1/2 Heisenberg systems using NPA constraints with RG coarse-graining. The first target is a periodic chain with `N = 200` and certified gap at most `1e-5`.

The project has established two independently replayable certificate routes at `N = 20`:

1. An exact-compatible RG/RDM8 relaxation uses one 13-bit dyadic `D = 4` tensor to assemble every recursive map and SDP coefficient. It certifies a gap of `6.967350626092331e-5`.
2. A full structured-NPA/SOHS relaxation without RG certifies a gap of `9.087309202669595e-6`, passing the requested accuracy threshold at the validation size.

The issue is not complete. `N = 50` and `N = 100` currently have numerical screens, not accepted lower certificates, and no `N = 200` solve is claimed.

## Claim discipline

This report uses “certified lower bound” only when all of the following hold:

- the solver termination status is `OPTIMAL`;
- the stored dual or SOHS transcript is replayed independently in exact arithmetic;
- every PSD correction has a rigorous eigenvalue enclosure and trace or operator-norm bound;
- the transcript is bound to the physical model or to an exactly compatible RG construction;
- the finite-energy upper endpoint is an exact-rational Rayleigh quotient;
- the lower endpoint does not exceed the upper endpoint.

Residual-adjusted and downward-rounded values from ordinary floating-point solves are numerical diagnostics. They are useful for selecting methods and estimating scale, but they are not described as certified bounds.

## Retained methods

### Exact-compatible RG/RDM8

The compressed-LTI hierarchy uses a two-spin super-site and a VUMPS-derived coarse-grainer. Exact U(1) charge blocks reduce the SDP, while an eight-physical-site RDM lift adds a structured-NPA-style local positivity constraint linked to the compressed variables.

The physical compatibility issue is resolved by quantizing the selected `D = 4` tensor to 13 dyadic bits and deriving `W2`, `MR`, `ML`, the Kronecker maps, the local Hamiltonian, and every quadratic SDP coefficient from that same tensor. Exact-rational replay then certifies the SDP that is proven to contain every physical recursive witness.

The retained implementation is split by responsibility:

- `rdm8_common.jl` contains shared constants and plotting only.
- `rg_rdm8_local.jl` certifies the local `D = 4`, `n_super = 6` comparison.
- `rg_rdm8_finite_n20.jl` certifies the finite periodic `N = 20`, `n_super = 10` comparison and validates the complete H0/H1/H3 predecessor chain.
- `finite_upper_certificate.jl` supplies the exact-rational Rayleigh upper endpoint.

### Full structured NPA and SOHS replay

The selected QMBCertify configuration is order 4, `extra = 1`, RDM9, PSD state optimality level 3, and linear state optimality. Its transcript covers the base Gram blocks, RDM terms, PSD state-optimality terms, linear state-optimality equalities, and Pauli-word reduction.

Validation deterministically rebuilds the source-pinned model without solving and compares the objective expression and `MAX_SENSE`, every equality coefficient, every PSD cone and variable mapping, variable indexing, and canonical word support. Only then does it replay the exact dyadic SOHS certificate and rigorous residual correction.

The retained implementation is `h10_qmbcertify_n20.jl` for configuration construction and `h11_full_sohs_certificate.jl` for the accepted solve and replay.

## Certified results

All numbers in this table are energy per physical site.

| Route | System | Certified lower bound | Certified upper bound | Rigorous gap | Reference run |
| --- | ---: | ---: | ---: | ---: | --- |
| Exact-compatible RG, compressed baseline | periodic `N = 20` | `-0.4453263421613551` | `-0.4452193264937391` | `1.0767350626092332e-4` | `20260729-h7-finite-n20-dyadic-k13` |
| Exact-compatible RG + RDM8 | periodic `N = 20` | `-0.44528862186638035` | `-0.4452193264937391` | `6.967350626092331e-5` | `20260729-h7-finite-n20-dyadic-k13` |
| Full structured NPA + SOHS | periodic `N = 20` | `-0.44522841380294176` | `-0.4452193264937391` | `9.087309202669595e-6` | `20260729-repro-julia115-full-sohs-n20` |

The RDM8 lift raises the exact-compatible RG certificate by `3.772029497473017e-5`, reducing its finite gap by about 35%. The full structured-NPA route independently passes the `1e-5` validation gate, but it does not yet include RG coarse-graining.

The H0 baseline also retains a rationally post-certified order-2 structured-NPA lower bound of `-0.44993915712553295` at `N = 20`. It is a correctness reference rather than a competitive endpoint.

## Numerical ladder and scaling evidence

| System and method | Endpoint | Status | Interpretation |
| --- | ---: | --- | --- |
| `N = 20`, `D = 4`, `n_super = 10` compressed | raw objective `-0.44532643488209284` | `OPTIMAL`, numerical | Establishes the predecessor used by the exact-compatible RG certificate. |
| `N = 50`, `D = 4`, `n_super = 25` compressed | raw objective `-0.443902769427835` | `OPTIMAL`, numerical | Improves the matched `D = 3` objective by `1.1728327994980914e-3`; diagnostic finite gap `4.258676536669781e-4`. |
| `N = 100`, `D = 4`, `n_super = 50` compressed | raw objective `-0.4438211577378427` | `OPTIMAL`, numerical | Improves the matched `D = 3` objective by `1.3336436432175303e-3`; diagnostic finite gap `5.924292490814831e-4`. |
| `N = 100`, `D = 5`, `n_super = 40` compressed | raw objective `-0.4436360816706201` | `OPTIMAL`, numerical | Improves the deeper `D = 4`, `n_super = 50` objective by `1.850760672226226e-4`; total solve time was 911.5 s. |
| `N = 50`, full structured-NPA configuration | raw objective `-0.4435171724808517` | `SLOW_PROGRESS`, rejected | 133,108 scalar variables, 7,902 equations, 109 PSD blocks, largest block 126, 214.4 s; numerical reference gap `4.0040134519e-5`. |

The `N = 100`, `D = 5` run recorded a live RSS snapshot of `43,763,876 KiB`; this is not a scheduler-measured maximum. This route therefore needs a high-memory system before an `N = 200` numerical ladder is sensible. In contrast, the structured-NPA `N = 50` screen remains locally tractable: its immediate blockers are accuracy and `SLOW_PROGRESS`, not compute capacity.

## Symmetry and constraint findings

- Exact U(1) blocks reduce each `D = 4` compressed `omega` PSD matrix from dimension 256 to a largest block of 80. At `n_super = 4`, matched dense and blocked objectives differ by `4.06e-12`, while solver wall time falls from 53.45 s to 1.52 s. At `n_super = 6`, the difference is `1.67e-9` and wall time falls from 175.03 s to 6.59 s; both differences are interpreted at the recorded residual scale.
- The shallow RDM8 probe at `n_super = 4` changes raw objectives by less than `4e-11` and is inactive there. At the exact-compatible active depths it becomes useful, producing certified improvements of `2.127205949921928e-5` locally and `3.772029497473017e-5` at finite `N = 20`.
- Six linear stationarity equations change the tested `D = 3` objective by less than `1e-13`.
- The `15 × 15` PSD state-optimality matrix changes the tested `D = 3` objective by about `1e-12`, below its residual scale.
- Direct SU(2) invariance changes the tested `D = 3` objective by about `1.46e-10`, also below its residual scale.
- The `D = 5` SU(2)-symmetric map passes dense/U(1) equivalence at `n_super = 4` with raw objective difference `4.39e-13`. At `n_super = 6`, its raw improvement over `D = 4` is `2.1534726996674802e-5`.
- Active-depth stationarity and PSD-optimality screens at exact-dyadic `N = 20` do not produce an accepted numerical candidate, so no strict certificate was attempted.

## Directions not retained as certificate routes

- Independently rounded floating-point RG maps can replay a dual certificate for the rounded SDP while failing the exact recursive witness identities. Those artifacts are not physical energy certificates and are not retained.
- A local coarse-grainer perturbation improved the compressed objective but did not resolve against the RDM8 residual scale. It is not retained as an accepted optimizer.
- An RDM10 lift increased memory without a residual-resolved gain. Its implementation and large transcript are not retained.
- Increasing the structured-NPA range to `extra = 2` returned `SLOW_PROGRESS` at `N = 20` for only about `8.7e-7` raw improvement.
- Increasing PSD state optimality to level 4 returned `SLOW_PROGRESS` at `N = 20`, increased the scalar-variable count by about 59%, and improved the raw objective by only about `2.28e-6`.

The final structured-NPA screens used run names matching `20260729-h12-*`. They are explicitly rejected numerical experiments; their conclusions are recorded here rather than retaining generated outputs.

## Reproduction and provenance

Generated results and coarse-grainer caches are intentionally excluded from version control. The commit contains the Julia `1.11.5` pin, package manifest, source code, tests, mathematical contracts, reference numbers, and complete commands needed to rebuild the two accepted `N = 20` certificate chains.

- `h11_full_sohs_certificate.jl` independently constructs its exact finite-energy upper endpoint, solves the selected structured-NPA model, writes the full transcript, and validates the generated result.
- The exact-compatible RG/RDM8 route is a generated H0 → H1 → H3 → H7 chain. Every stage records its explicit input paths and validates their metadata and numerical gates.
- A fresh Julia `1.11.5` H11 run reproduced the reported numerical objective, certified lower endpoint, exact upper endpoint, and `9.087309202669595e-6` gap; its full replay and tamper suite passed 32 tests.
- A fresh fixed-seed D=3 VUMPS run and H0 ladder also reproduced all four OPTIMAL compressed points and the rationally post-certified structured-NPA baseline.
- The historical run names in the tables identify the runs from which the quoted numbers were taken; they are report labels, not committed directories.
- [scripts/heisenberg1d/README.md](../scripts/heisenberg1d/README.md) contains the clean-clone reproduction commands and the optional post-generation replay commands.

The `N = 50` and `N = 100` entries are historical numerical screening evidence, not accepted certificate routes. Their scripts remain available with explicit predecessor arguments, but this commit does not claim a documented clean-clone reproduction chain for those exploratory ladders.

## Challenge status and next decision

| Rung | Current state | Gate |
| --- | --- | --- |
| `N = 20` | certified by two routes | Full structured NPA passes `1e-5`; RG/RDM8 is certified at `6.97e-5`. |
| `N = 50` | numerical screens only | Need `OPTIMAL` termination, a stricter model with predicted gap near `1e-5`, and exact post-certification. |
| `N = 100` | high-memory RG numerical screen only | Need a strict certificate before its gap can be claimed. |
| `N = 200` | not run | Launch only after a smaller rung demonstrates both accuracy and a credible resource projection. |

The next work should remain local: improve conditioning and constraint selection for the `N = 50` structured-NPA model, or optimize an exactly representable RG map using the cheap local and `N = 20` certificate gates. High-performance computing becomes justified for the RG ladder only after one of those method gates predicts that `N = 200` can meet the accuracy target.

## Validation commands

The result-independent suite is:

```bash
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_metadata.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h4.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_rg_rdm8.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl
```

After generating a certificate, pass its result path to the corresponding test entrypoint for full transcript replay and tamper checks.

The mathematical contracts are [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md), [HYBRID_RDM_CONTRACT.md](HYBRID_RDM_CONTRACT.md), [DUAL_POSTCERTIFICATION_CONTRACT.md](DUAL_POSTCERTIFICATION_CONTRACT.md), [HIGHER_D_CONTRACT.md](HIGHER_D_CONTRACT.md), and [FULL_SOHS_CERTIFICATE.md](FULL_SOHS_CERTIFICATE.md).
