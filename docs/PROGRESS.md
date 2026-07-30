# Progress report — challenge #49

Updated 2026-07-30. This is the authoritative status report for [QuantumBFS/quantum.harness#49](https://github.com/QuantumBFS/quantum.harness/issues/49) and PR #219.

For a concise four-section summary, see the standalone [Challenge Report](challenge-report/report.html).

## Executive status

The challenge asks for certified energy-density lower bounds for increasingly large spin-1/2 Heisenberg systems using NPA constraints with RG coarse-graining. Its four official targets are:

| Official target | Accepted state |
| --- | --- |
| 1D Heisenberg, `N = 200`, gap at most `1e-5` | Not complete. The accepted fused certificate is at `N = 20`; larger screens have not passed the strict certificate gate. |
| 1D J1-J2 Heisenberg, `N = 100`, gap at most `1e-3` | Not attempted as a fused certificate target. |
| 2D Heisenberg, `16 × 16`, gap at most `1e-3` | Not attempted as a fused certificate target. |
| 2D J1-J2 Heisenberg, `10 × 10`, gap at most `1e-2` | Not attempted as a fused certificate target. |

Formal target completion is therefore `0 / 4`. The `N = 20` result is a completed validation-scale proof of concept for the fused algorithm, not completion of any full-size challenge target. The achieved system size is `20 / 200` of the first target, but this ratio must not be interpreted as a compute- or research-effort percentage because the SDP cost grows superlinearly.

The project has established three independently replayable certificate routes at `N = 20`:

1. An exact-compatible RG/RDM8 relaxation uses one 13-bit dyadic `D = 4` tensor to assemble every recursive map and SDP coefficient. Its exact certificate endpoints give a gap of `6.929537264127523e-5`; rounding the lower bound downward to six decimal places gives the more conservative published gap `6.967350626092331e-5`.
2. A full structured-NPA/SOHS relaxation without RG certifies a gap of `9.040106778343054e-6`, passing the requested accuracy threshold at the validation size.
3. A unified structured-NPA + RG relaxation reconstructs an eight-spin RDM from the full moment basis, audits the moment-to-RDM coefficients exactly, and attaches the exact dyadic RG hierarchy in the same SDP. At `D = 4`, `n_super = 5`, it certifies a gap of `5.214353500759827e-6`, passing the accuracy gate and strictly improving both parent relaxations.

The issue is not complete. At `N = 30`, matched three-route screens solve cleanly and show a positive raw fusion gain, but exact replay gives a negative strict gain. In the first large-size HPC campaign, the `N = 50` fused solves stalled, the `N = 100` structured-NPA parent stalled before a fused result was available, and the queued `N = 200` job was canceled without consuming compute. These are diagnostic failures, not accepted lower certificates.

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

### Unified structured NPA + RG

The unified model dualizes the accepted SOHS formulation into physical moments and reconstructs an eight-spin RDM for four two-spin super-sites. An independent U(1)-blocked PSD variable is linked to that RDM by coefficients whose complete dyadic accumulation is audited below the Float64 exact-integer limit. The already audited 13-bit dyadic `W2`, `MR`, and `ML` maps then attach the recursive RG marginal constraints before the model is dualized back to a stable SOHS-like maximization. This is a genuine intersection of the two relaxations, not a comparison of separately optimized bounds.

The transcript binds the local-RDM coefficient audit, exact coarse-grainer, U(1) charges, complete structured-NPA stationarity system, all SOHS PSD blocks, every RG affine PSD slack, and the exact equality-family scales `(1, 16, 1)` to a deterministically rebuilt model. Exact replay combines Pauli-word residual correction with analytic trace corrections for `rho4` and every `omega` group. The retained implementation is `structured_rg_hybrid.jl`, `structured_rg_certificate.jl`, and `structured_rg_local.jl`; the mathematical contract is [STRUCTURED_RG_CONTRACT.md](STRUCTURED_RG_CONTRACT.md).

## Certified results

All numbers in this table are energy per physical site.

| Route | System | Certified lower bound | Certified upper bound | Rigorous gap | Reference run |
| --- | ---: | ---: | ---: | ---: | --- |
| Exact-compatible RG, compressed baseline | periodic `N = 20` | `-0.4453263421613551` | `-0.4452193264937391` | `1.070156676160054e-4` | `20260729-h7-finite-n20-dyadic-k13` |
| Exact-compatible RG + RDM8 | periodic `N = 20` | `-0.44528862186638035` | `-0.4452193264937391` | `6.929537264127523e-5` | `20260729-h7-finite-n20-dyadic-k13` |
| Full structured NPA + SOHS | periodic `N = 20` | `-0.44522836660051746` | `-0.4452193264937391` | `9.040106778343054e-6` | `20260730-stage-report-full-sohs-n20` |
| Matched full structured NPA parent | periodic `N = 20` | `-0.4452272566443788` | `-0.4452193264937391` | `7.930150639712252e-6` | `20260730-structured-rg-rdm8-scaled-n20-d4-n5` |
| Full structured NPA + exact-compatible RG/RDM8 | periodic `N = 20` | `-0.44522454084723984` | `-0.4452193264937391` | `5.214353500759827e-6` | `20260730-structured-rg-rdm8-scaled-n20-d4-n5` |

The RDM8 lift raises the exact-compatible RG certificate by `3.772029497473017e-5`, reducing its finite gap by about 35%. The full structured-NPA route independently passes the `1e-5` validation gate, but it does not yet include RG coarse-graining.

The unified model proves that all retained structured constraints and RG constraints can be certified together through an exact physical witness. Its raw objective improves by `3.9178868482814266e-7`; exact mixed replay improves the matched structured-NPA endpoint by `2.715797138952425e-6`. Multiplying only the `omega4` link equalities by the exact factor `16` leaves the feasible set unchanged and exposes the strict improvement by improving the numerical conditioning. Both solves are `OPTIMAL / FEASIBLE_POINT`, and the independent validator rebuilds the complete model and rejects altered coefficients or equality scales after recomputed hashes.

The H0 baseline also retains a rationally post-certified order-2 structured-NPA lower bound of `-0.44993915712553295` at `N = 20`. It is a correctness reference rather than a competitive endpoint.

## Numerical ladder and scaling evidence

| System and method | Endpoint | Status | Interpretation |
| --- | ---: | --- | --- |
| `N = 20`, `D = 4`, `n_super = 10` compressed | raw objective `-0.44532643488209284` | `OPTIMAL`, numerical | Establishes the predecessor used by the exact-compatible RG certificate. |
| `N = 20`, unified structured NPA + RG/RDM8, `n_super = 5` | certified lower bound `-0.44522454084723984` | `OPTIMAL`, certified | 76,501 scalar variables; rigorous gap `5.214353500759827e-6`; strictly improves both parent routes. |
| `N = 30`, unified structured NPA + RG/RDM8, `n_super = 5` | raw gain `2.1437405645086116e-7` | `OPTIMAL`, exploratory | Fixed-depth raw tightening persists, but same-process strict replay loses `1.3332406440951713e-5` because moment residuals grow; not a formal finite-size certificate. |
| `N = 30`, normalized unified model, `n_super = 8`, coefficient scale `512` | raw gain `1.6160551729593742e-6` | all three routes `OPTIMAL`, exploratory | The deeper RG hierarchy solves, but exact replay gives strict gain `-3.2971641651619926e-5`; the strict-improvement gate fails. |
| `N = 30`, normalized unified model, `n_super = 8`, coefficient scale `16384` | raw gain `1.5211146759397387e-6` | all three routes `OPTIMAL`, exploratory | The equivalent scaling control also fails strict replay, with strict gain `-3.580298304150444e-5`. |
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

Generated results and coarse-grainer caches are intentionally excluded from version control. The commit contains the Julia `1.11.5` pin, package manifest, source code, tests, mathematical contracts, reference numbers, and complete commands needed to rebuild the three accepted `N = 20` certificate chains.

- `h11_full_sohs_certificate.jl` independently constructs its exact finite-energy upper endpoint, solves the selected structured-NPA model, writes the full transcript, and validates the generated result.
- `structured_rg_local.jl --link-mode=rdm8` solves both the full structured-NPA baseline and the fused `n_super = 5` model, writes their exact transcripts, and validates the unified result against a rebuilt model.
- The exact-compatible RG/RDM8 route is a generated H0 → H1 → H3 → H7 chain. Every stage records its explicit input paths and validates their metadata and numerical gates.
- A fresh Julia `1.11.5` H11 run reproduced the numerical objective, current energy-aware certified lower endpoint, exact upper endpoint, and `9.040106778343054e-6` gap. The generator immediately rebuilt and validated the complete model; the optional replay command exercises the additional tamper suite.
- A fresh fixed-seed D=3 VUMPS run and H0 ladder also reproduced all four OPTIMAL compressed points and the rationally post-certified structured-NPA baseline.
- The historical run names in the tables identify the runs from which the quoted numbers were taken; they are report labels, not committed directories.
- [scripts/heisenberg1d/README.md](../scripts/heisenberg1d/README.md) contains the clean-clone reproduction commands and the optional post-generation replay commands.

The `N = 50` and `N = 100` entries are historical numerical screening evidence, not accepted certificate routes. Their scripts remain available with explicit predecessor arguments, but this commit does not claim a documented clean-clone reproduction chain for those exploratory ladders.

## Challenge status and next decision

| Rung | Current state | Gate |
| --- | --- | --- |
| `N = 20` | certified by three routes | The unified route is the tightest, passes `1e-5` at `5.21e-6`, and strictly improves both parent methods. |
| `N = 30` | three-route `n_super = 8` screens solve, strict gate fails | Raw fusion gains remain positive, but ordinary moment residuals make the strict gains negative at both tested equivalent coefficient scales. |
| `N = 50` | compressed numerical screens exist; matched fused HPC jobs failed | The structured parent completed in the first campaign, but the fused SDP returned `MSK_RES_TRM_STALL`; no strict fused certificate exists. |
| `N = 100` | high-memory compressed screens exist; matched structured baseline failed | The first matched HPC job stalled in the structured-NPA parent, before an accepted fused result. |
| `N = 200` | not run | The queued job was canceled after the smaller-size failures showed that the same formulation was not ready. |

The local fusion gate is complete, but the size trend is not. Exact dyadic normalization and RG equality rescaling remove the exponential `omega_k` magnitude growth and allow the normalized `N = 30`, `n_super = 8` formulation to solve. They do not yet control the accumulated ordinary moment residuals well enough for strict replay. The next technical blocker is therefore the numerical conditioning of the structured-NPA coefficient-matching equations, not additional CPU parallelism. Further `N = 50`, `N = 100`, or `N = 200` runs are scientifically justified only after an `N = 30` formulation passes the same strict-improvement gate without post-processing away the residual loss.

## Validation commands

The result-independent suite is:

```bash
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_metadata.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h4.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_rg_rdm8.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_structured_rg_hybrid.jl
```

After generating a certificate, pass its result path to the corresponding test entrypoint for full transcript replay and tamper checks.

The mathematical contracts are [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md), [HYBRID_RDM_CONTRACT.md](HYBRID_RDM_CONTRACT.md), [DUAL_POSTCERTIFICATION_CONTRACT.md](DUAL_POSTCERTIFICATION_CONTRACT.md), [HIGHER_D_CONTRACT.md](HIGHER_D_CONTRACT.md), [FULL_SOHS_CERTIFICATE.md](FULL_SOHS_CERTIFICATE.md), and [STRUCTURED_RG_CONTRACT.md](STRUCTURED_RG_CONTRACT.md).
