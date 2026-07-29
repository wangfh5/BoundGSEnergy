# H2 positive-semidefinite ground-state optimality contract

This note defines the third structured-NPA strengthening tested against the compressed-LTI relaxation. It is an infinite-chain method probe and does not become evidence for a finite challenge size unless the same constraint is applied under [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md).

## Matrix ground-state condition

For local operators `a1, ..., am`, every infinite-system ground state satisfies

```math
C_{ij}=\omega(a_i^\dagger[H,a_j]),\qquad C\succeq0.
```

This is the zero-temperature limit of the matrix ground-state condition in Fawzi, Fawzi, and Scalet, *Certified algorithms for equilibrium states of local quantum Hamiltonians*. It strengthens the linear equations in [STATE_OPTIMALITY_CONTRACT.md](STATE_OPTIMALITY_CONTRACT.md).

## Minimal local basis

Each `ai` acts on the middle of three consecutive four-dimensional super-sites. As before, only

```math
H_{\mathrm{touch}}=H_{B,12}+h_2+H_{B,23}
```

contributes to the commutator, so every entry of `C` is a linear function of the existing 64-dimensional `rho3`.

The identity direction is removed because its commutator vanishes. The implemented real basis contains 15 traceless matrices: three diagonal differences, six symmetric off-diagonal generators, and six antisymmetric off-diagonal generators. Complex linear combinations span the complete traceless operator algebra on the central super-site. Six linear stationarity equations are imposed alongside the 15 by 15 positive-semidefinite matrix so that adding an identity component does not change the condition.

## Physical-state regression

On the exact periodic six-spin ground state, the constructed matrix is symmetric to numerical precision and positive semidefinite. This simultaneously checks the dagger, commutator sign, local Hamiltonian support, and most-significant-factor-first index order before any SDP is solved.

## Acceptance gate

- Baseline and PSD-optimality solves terminate `OPTIMAL`.
- The PSD-optimality result stays below the Bethe energy within the configured numerical tolerance.
- The raw objective improvement exceeds `2e-7`.
- Solver residuals, wall time, operator basis, exact-state eigenvalue check, support, and coarse-grainer fingerprint are recorded.

The first point is dense D = 3, n = 4. A symmetry-aware D = 4 point and finite-size port are attempted only if this gate passes.

## Probe result

**Rejected as a tightening on 2026-07-28.** The reference run used the name `20260728-h2-psd-optimality-d3`; generated outputs are not versioned.

The exact periodic N = 6 regression produces a symmetric optimality matrix with minimum eigenvalue `0.0562`, confirming the implemented dagger, sign, support, and index order. The baseline and PSD-optimality formulations both terminate `OPTIMAL` and publish the same six-decimal bound `-0.456387`. Their raw solver objectives change by `-1.02e-12`, below the hybrid primal feasibility residual `1.37e-8` and far below the `2e-7` gate. Mosek time changes from 2.05 s to 2.02 s. The D = 3 gate fails, so no D = 4 or finite-size run is made.
