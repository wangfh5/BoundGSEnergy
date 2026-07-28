# 2D J1-J2 Heisenberg model

This folder is the staged workspace for the frustrated spin-1/2 square-lattice model

```math
H = J_1 \sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j + J_2 \sum_{\langle\langle i,j \rangle\rangle} \mathbf{S}_i \cdot \mathbf{S}_j.
```

The challenge target is a certified energy-density lower bound for a `10 × 10` lattice with absolute accuracy `1e-2`. Work starts only after the 2D Heisenberg coarse-graining gate passes.

An energy bound alone cannot resolve the debated phase structure. Any phase claim requires a separate observable-certification target and validation protocol.

This directory intentionally contains no runnable scripts yet. The staged work and acceptance gates are defined in [`docs/RESEARCH_PLAN.md`](../../docs/RESEARCH_PLAN.md).
