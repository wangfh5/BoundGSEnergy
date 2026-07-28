# Ideas Session — 2026-07-27 19:20

## Context (carried in from challenge registration)

- Challenge #49 (BoundGSEnergy, PR #219): combine NPA SDP hierarchy with RG
  coarse-graining (PRX 14, 021008 Sec. III-D-2) for certified ground-state
  energy lower bounds, spin-1/2 Heisenberg systems.
- PLAN.md ratified: ladder N=20→50→100→200; success = coarse-grained SDP
  visibly cheaper at matched accuracy; kill criteria written.
- Assets: working QMBCertify+Mosek pipeline (Track 2), 34-ref survey library
  incl. rendered Kull-Schuch-Dive-Navascués (2212.03014).
- User: Fo-Hong Wang, QMC researcher, learning mindset, skipped Track 4.

## Phase 2, exchange 1 — the key conceptual finding

Mentor read the local rendered 2212.03014 (Kull et al., PRX 14, 021008).
Key finding: their coarse-graining acts on the *relaxation constraints*, NOT
on the Hamiltonian. Architecture:
- LTI problem: hierarchy of compatibility constraints between reduced states
  ρ^(2) ← ρ^(3) ← ... ← ρ^(n).
- Coarse-graining maps W_m (from MPS/TTN, e.g. variationally optimized ones)
  compress each ρ^(m) into a small variable ω^(m); compatibility of the maps
  (extra maps L_m, R_m, Eq. 9) lets the chain of constraints close on the
  compressed variables.
- Result: rigorous lower bounds for infinite 1D chains, 1–2 orders of
  magnitude better than un-compressed truncations; same MPS that give
  variational upper bounds work well as coarse-grainers for the lower bound.

Consequence surfaced to user: naive "block 2 spins → effective Heisenberg
model → QMBCertify" would be *variational* (Wilson-style truncation of H
keeps no certificate); the lower-bound property only survives if
coarse-graining happens at the constraint/relaxation level à la Kull-Schuch.

Fork presented:
1. (Recommended) MVE = reimplement compressed-LTI (Kull-Schuch) at small
   scale for the 1D Heisenberg chain with an MPS coarse-grainer (ITensors
   available), reproduce a published bound, THEN add NPA-style
   moment/symmetry constraints and measure the tightening.
2. Keep QMBCertify moment-SDP architecture, swap in a coarse-grained
   monomial basis — cheaper coding but conceptually murky (the W-map
   compatibility issue is exactly what Kull-Schuch solves).

Open question for user: check whether Kull et al. published code before
writing our own.

## Phase 2, exchange 2 — MVE ratified, route 1 chosen

User picked route 1: reimplement Kull-Schuch compressed LTI before combining
with NPA. Mentor read Section II of 2212.03014 in depth and extracted the
exact relaxation (Eq. 14):

- Variables: ρ^(3) (8x8 state) + compressed ω^(4..n), each d²D² × d²D²
  (D = MPS bond dim; D=2 → 16×16).
- Constraints: ρ^(3) ⪰ 0, normalized; LTI TrL = TrR; compressed chain via
  L_A/R_A maps built from the MPS tensor A; ω ⪰ 0. NOTE: ω are NOT
  trace-1 (W is not an isometry) — normalization only on ρ^(3).
- Objective: min Tr(h ρ^(2)), nearest-neighbour bond energy density.
- Coarse-grainer: VUMPS uniform MPS of the SAME Hamiltonian (paper's
  heuristic; works well per their Fig. 2).
- Order of bounds (self-check): E_relax(n,D) ≤ E_LTI(n) ≤ E_TI = e₀.

Verification is free: Bethe-ansatz exact e₀ = 1/4 − ln 2 ≈ −0.443147
(analytic anchor, no ED needed) + monotonicity + E_LTI(n) computable
uncompressed for n ≤ 8.

Only hard part flagged honestly: Eqs. 13/14 are images in the paper; the
precise index order of the L_A/R_A contractions must be derived from the
text. Mitigation: the e₀ ceiling exposes a wrong map direction immediately.

MVE: reproduce Fig. 2b for the Heisenberg chain (D = 2, 3; n ladder to
~150). Confirmed by user ("可以").

## Implementation task split (test-first)

- T1: VUMPS uniform MPS for 1D Heisenberg, D=2,3 (MPSKit).
  Test: variational energy must approach e₀ from ABOVE, tighter with D.
- T2: Build W₂ / L_A / R_A maps + compressed SDP (JuMP + Mosek).
  Test: E_relax(n,D) ≤ e₀ always (violation ⇒ index bug).
- T3: n ladder 4→150, D=2,3; ΔE vs n curve.
  Test: E_relax(n,D) ≤ E_LTI(n) with E_LTI solved uncompressed for n ≤ 8.
- T4: Results + figure to tracks/polyopt/results/<run>/, commit to PR branch.

## Phase 2, exchange 3 — MVE EXECUTED AND PASSED (2026-07-27 late session)

Route 1 implemented in-session (subagent delegation hit an API quota 403, so
the mentor coded it directly — noted to user).

- T1: VUMPS uniform MPS, D=2 (−0.42791) and D=3 (−0.43578), both variational
  (above e₀), ‖B‖ ~ 1e-13. Had to switch to a 2-site unit cell (1-site
  oscillated) and to Float64 eltype (ComplexF64 gauge blocks the real SDP).
- T2/T3: compressed-LTI SDP (super-spin d̂=4 formulation) with JuMP+Mosek.
  Debugging caught by test-first asserts: (i) kron/MSF vs reshape/LSF index
  mismatch in ptrace helpers and W/M maps; (ii) Mosek SLOW_PROGRESS at
  stalled residuals 2e-8 — fixed with honest 1e-7 tolerances.
- Results (all OPTIMAL, all below e₀ = 1/4−ln2, monotone in n):
  D=2 saturates at gap 4.69e-3 from super-n=15 (30 physical) onward, up to
  super-120 (240 physical) in 24 s. D=3 saturates at gap 2.01e-3 by super-50
  (100 physical), 306 s. Exact LTI computed for n=4, 6 (n=8 infeasible
  uncompressed — itself the demonstration of the compression advantage).
  Compressed D=3 at 30 physical spins: gap 2.1e-3 vs LTI(6) gap 2.4e-2 —
  >10× accuracy at comparable cost. Fig. 2b shape reproduced (n⁻²
  continuation + D-dependent saturation).
- PLAN floor criterion MET: coarse-grained SDP is visibly cheaper at matched
  accuracy. Scripts committed to PR branch (f647e72); figure + JSON in
  tracks/polyopt/results/20260727-kullschuch-mve/.

Next candidates (open): (a) the actual #49 combination — add NPA
moment/symmetry constraints to the compressed relaxation or use NCTSSoS with
coarse-grained variables; (b) climb D (5, 7) toward the paper's n_eff≈60;
(c) ladder rung 2 (J1-J2); (d) wrap session.
