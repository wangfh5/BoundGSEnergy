# H2 linear ground-state stationarity contract

This note defines the second structured-NPA strengthening tested against the compressed-LTI relaxation. It is an infinite-chain method probe and does not become evidence for a finite challenge size unless the same constraint family is applied under [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md).

## Ground-state condition

For an infinite-system ground state `omega` and every local operator `a`, linear stationarity requires

```math
\omega([H,a])=0.
```

This is the zero-temperature stationarity constraint in Eq. 45d of Fawzi, Fawzi, and Scalet, *Certified algorithms for equilibrium states of local quantum Hamiltonians*. It is weaker than their full positive-semidefinite ground-state condition `[\omega(a_i^\dagger[H,a_j])]_{ij} >= 0`, but it is linear and introduces no new PSD variable.

## Local support and real formulation

Let `a` act on the middle of three consecutive super-sites. Only the two inter-super-site bonds and the middle intra-super-site bond fail to commute with it, so

```math
H_{\mathrm{touch}}=H_{B,12}+h_2+H_{B,23}
```

is supported on the existing 64-dimensional `rho3` variable. The implemented constraint is

```math
\operatorname{Tr}\!\left(\rho_3[H_{\mathrm{touch}},a]\right)=0.
```

The SDP is real and `rho3` is symmetric. A real symmetric generator produces a real antisymmetric commutator whose trace against `rho3` vanishes identically. The nontrivial minimal family therefore uses the six real antisymmetric generators `Aij = |i><j| - |j><i|` of one four-dimensional super-site. Each `[H_touch,Aij]` is real symmetric, and the six resulting matrices have linear rank six.

## Physical-state and symmetry checks

Terms outside `H_touch` commute with the local generator, so every infinite-chain ground state satisfies the constraint. As an index and sign regression, the six matrices are evaluated on the exact periodic six-spin ground state; the largest residual is `3.5e-16`.

With super-site charges `[-2,0,0,2]`, only `A23` acts within a fixed U(1) charge sector. Projection of the six matrices onto the existing `rho3` charge blocks therefore leaves one nonzero independent constraint. The first gate is tested only in the dense D = 3, n = 4 formulation; the U(1)-blocked D = 4 point is run only if D = 3 passes.

## Acceptance gate

- Baseline and stationarity solves terminate `OPTIMAL`.
- The stationarity result stays below the Bethe energy within the configured numerical tolerance.
- The raw objective improvement exceeds `2e-7`.
- Solver residuals, wall time, support, operator-family rank, exact-state residuals, and coarse-grainer fingerprint are recorded.

If the D = 3 gate fails, the candidate is retained as a negative result and is not expanded to D = 4 or ported to the finite-size solver.

## Probe result

**Rejected as a tightening on 2026-07-28.** The reference run used the name `20260728-h2-linear-stationarity-d3`; generated outputs are not versioned.

The baseline and six-equation stationarity formulations both terminated `OPTIMAL`, stayed below the Bethe value, and publish the same six-decimal bound `-0.456387`. Their raw solver objectives change by `-9.48e-14`, below the largest recorded primal feasibility residual `5.28e-11` and far below the `2e-7` acceptance threshold. Mosek time changes from 1.95 s to 1.84 s. The D = 3 gate fails, so the single surviving U(1) equation is not tested at D = 4 and this candidate is not ported to the finite-size solver.
