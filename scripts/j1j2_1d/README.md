# 1D J1-J2 Heisenberg chain

This folder is the staged workspace for the frustrated spin-1/2 chain

```math
H = J_1 \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1} + J_2 \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+2}.
```

The challenge target is a certified energy-density lower bound for 100 spins with absolute accuracy `1e-3`. Work starts only after the 1D Heisenberg pipeline passes its finite-size `N = 50` gate.

The Majumdar-Ghosh point `J2/J1 = 0.5` is the first exact end-to-end anchor. Other coupling points must be fixed before compute and recorded with the Hamiltonian normalization and boundary conditions.

This directory intentionally contains no runnable scripts yet. The staged work and acceptance gates are defined in [`docs/RESEARCH_PLAN.md`](../../docs/RESEARCH_PLAN.md).
