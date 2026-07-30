# H2 hybrid RDM contract

This note defines the first structured-NPA strengthening of the compressed-LTI relaxation. It is an infinite-chain feasibility probe; it becomes evidence for a finite challenge size only after the same construction is applied under the finite-size contract in [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md).

## Variables and map

Each super-site contains two physical spins and has Hilbert-space dimension four. Let `rho3` be the existing three-super-site density matrix, let `W2: C^16 -> C^(D^2)` be the existing two-super-site coarse-graining map, and introduce a four-super-site density matrix

```math
\rho_4 \succeq 0,\qquad \operatorname{Tr}(\rho_4)=1.
```

Translation compatibility requires

```math
\operatorname{Tr}_{1}(\rho_4)=\rho_3,\qquad \operatorname{Tr}_{4}(\rho_4)=\rho_3.
```

In the established most-significant-factor-first index order, define

```math
C=I_4\otimes W_2\otimes I_4.
```

The existing compressed boundary variable is linked exactly by

```math
\omega_4=C\rho_4 C^\mathsf{T}.
```

No Hamiltonian coefficient or objective normalization changes.

## Why every physical state remains feasible

For any translation-invariant infinite-chain state, take its reduced density matrices on three and four consecutive super-sites. They are positive semidefinite, normalized, and satisfy the two partial-trace identities above. Applying the same real coarse-graining map to the middle two super-sites gives `omega4 = C rho4 C'`.

The original `omega4` marginal constraints follow directly:

```math
\operatorname{Tr}_{4}(\omega_4)=(I_4\otimes W_2)\rho_3(I_4\otimes W_2)^\mathsf{T},
```

```math
\operatorname{Tr}_{1}(\omega_4)=(W_2\otimes I_4)\rho_3(W_2\otimes I_4)^\mathsf{T}.
```

The strengthened feasible set therefore still contains the image of every physical translation-invariant state. The exact SDP optimum remains a lower bound; a floating-point solver value becomes a strict numerical certificate only after the dual post-step in [DUAL_POSTCERTIFICATION_CONTRACT.md](DUAL_POSTCERTIFICATION_CONTRACT.md).

The implementation checks both identities above without a solver on a deterministic cyclically translation-invariant four-super-site test state. The recorded maximum residual is `3.47e-18`.

## Relation to structured NPA

`rho4` is a positive semidefinite reduced-density-matrix block on eight contiguous physical spins. It is the direct coarse-grained analogue of the `rdm=8` strengthening in QMBCertify, with an additional linear link to the compressed hierarchy.

The existing relaxation constrains the two marginals of `omega4` but does not require a common positive semidefinite four-super-site lift. This lift is potentially strict whenever `W2` loses information. At the recorded coarse-grainers, `rank(W2)=9<16` for the ordinary `D=3` tensor and `rank(W2)=11<16` for the symmetry-preserving `D=4` tensor.

## U(1) blocks and acceptance gate

In integer `2Sz` units, the four-super-site charge multiplicities are `[1, 8, 28, 56, 70, 56, 28, 8, 1]`, so the 256-dimensional RDM becomes nine PSD blocks with largest dimension 70.

The H2 feasibility gate is:

- baseline and hybrid solves terminate `OPTIMAL`;
- `E_compressed <= E_hybrid <= 1/4 - ln(2)` within the configured numerical tolerance;
- `E_hybrid - E_compressed > 2e-7`;
- the raw objective separation exceeds the sum of the two independent solver-residual scales;
- both solves pass strict dual post-certification;
- solver residuals, wall time, block inventory, system scope, support, and coarse-grainer fingerprint are recorded.

## Depth-dependent results

The lift is inactive at the shallow `n = 4` probe: the D = 3 and D = 4 raw objectives change by less than `4e-11`. This does not establish redundancy after the recursive compressed constraints become active.

The retained RG/RDM8 implementation reconstructs the D = 4 model from one exactly compatible 13-bit dyadic tensor. At n = 6, its physical compressed and RDM8 certificates differ by `2.127e-5`. At finite N = 20 and n = 10, they differ by `3.772e-5`; the RDM8 endpoint `-0.44528862186638035` and exact dyadic Rayleigh upper bound `-0.4452193264937391` give the rigorous finite gap `6.9295372641e-5`. Rounding the certified lower endpoint downward to six decimal places gives the separately recorded conservative published gap `6.9673506261e-5`.

The runnable entrypoints are `scripts/heisenberg1d/rg_rdm8_local.jl` and `scripts/heisenberg1d/rg_rdm8_finite_n20.jl`. Independently rounded-map experiments are excluded from the retained certificate chain.
