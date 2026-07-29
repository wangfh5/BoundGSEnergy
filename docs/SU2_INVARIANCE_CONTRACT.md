# H2 physical SU(2)-invariance contract

This note defines the final structured-NPA strengthening tested against the compressed-LTI relaxation. It is an infinite-chain method probe and does not become evidence for a finite challenge size unless the same constraint is applied under [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md).

## Why an invariant ground state may be selected

The spin-1/2 Heisenberg Hamiltonian is invariant under global SU(2) rotations. Averaging any ground state over SU(2) produces an invariant mixed ground state with the same energy. The relaxation may therefore be restricted to local marginals of invariant states without excluding the ground energy.

This argument does not require the fixed D = 3 coarse-grainer to be covariant. The symmetry-averaged physical ground state is itself a valid physical state, and applying any fixed coarse-graining maps to its marginals produces feasible compressed variables.

## Implemented equations

On the six physical spins represented by `rho3`, define the total-spin generators

```math
S^x_{\mathrm{tot}}=\sum_{i=1}^{6}S_i^x,\qquad S^z_{\mathrm{tot}}=\sum_{i=1}^{6}S_i^z.
```

The real symmetric marginal is constrained by

```math
[\rho_3,S^x_{\mathrm{tot}}]=0,\qquad [\rho_3,S^z_{\mathrm{tot}}]=0.
```

Commutation with these two noncommuting generators implies full SU(2) invariance. Each commutator is antisymmetric, so only its strict upper triangle is imposed. The ambient real symmetric matrix space has dimension 2080; for six spin-1/2 sites the SU(2)-invariant real symmetric subspace has dimension 76, giving 2004 independent equations despite 4032 generated upper-triangle equations.

## Physical-state regression and gate

The exact periodic N = 6 ground-state density matrix commutes with both generators to numerical precision. The D = 3, n = 4 baseline and symmetry-constrained solves must terminate `OPTIMAL`, stay below the Bethe value, and improve the raw objective by more than `2e-7`.

If this final candidate fails, H2 closes as a documented negative feasibility stage and no tested strengthening is ported to the finite-size solver.

## Probe result

**Rejected as a tightening on 2026-07-28.** The reference run used the name `20260728-h2-su2-invariance-d3`; generated outputs are not versioned.

The exact periodic N = 6 density matrix commutes with total `Sx` and `Sz` with maximum residual `3.05e-16`. The baseline and SU(2)-invariant formulations both terminate `OPTIMAL` and publish the same six-decimal bound `-0.456387`. Their raw solver objectives change by `-1.46e-10`, below the hybrid primal feasibility residual `1.12e-8` and far below the `2e-7` gate. Mosek time rises from 1.95 s to 5.16 s. The D = 3 gate fails and H2 closes without a strengthening to port.
