# NOTES — SDP certification of quantum spin-system ground states

Built by /survey on 2026-07-27 (7 parallel strategies; themes A/C/E/F kept).
Anchor reproduction: Wang et al. 2026 Table 8 (beginner-training Track 2).

## Field landscape

**Certification main line (NPA → many-body).** The field converged from the
quantum-information NPA hierarchy [@navascues_2008_convergent;
@0903.4368 in this folder] to many-body ground states through the PRX 2024
reformulation [@wang_2024_certifying in this folder] and the 2026 structured
NPA paper reaching 16×16 lattices [@wang_2026_scalable in this folder].
Active groups: ICFO Barcelona (Acín, Renou, Frérot, Wang, Jansen) +
LAAS-CNRS Toulouse (Magron, Legat). Direct descendants: certified phase
diagrams via energy-window observable bounds [@jansen_2025_mapping],
steady-state certification of open systems [@mortimer_2024_certifying],
certification from partial experimental data [@mortimer_2026_bounding],
finite-level convergence-rate theory [@klep_2026_quantitative], systematic
moment optimization [@flora_2026_moment], and the certified-algorithm line
for equilibrium states [@fawzi_2023_certified; @fawzi_2023_entropy].

**Gap & observable certification (adjacent cluster).** Smaller and newer:
spectral-gap certificates for frustration-free systems [@rai_2024_hierarchy]
and gap bootstrapping of spin chains [@nancarrow_2022_bootstrapping]. Beyond
these two, gap bounds are dominated by analytic finite-size criteria, not SDP.

**Cross-method benchmark infrastructure (upper-bound sources).** Every SDP
lower bound is validated against these: Sandvik's QMC tables for the 2D
Heisenberg model [@sandvik_1997_finitesize; @sandvik_2026_highprecision],
DMRG at maximal frustration [@gong_2014_plaquette; @wang_2018_critical],
VMC/NQS/PEPS variational energies [@morita_2015_quantum;
@choo_2019_twodimensional; @liu_2022_gapless], the community variational
compendium [@wu_2024_variational], and the Hubbard-model benchmark suites
[@leblanc_2015_solutions; @qin_2016_benchmark; @sawaya_2024_hamlib].

**Known failure modes (negative results).** The bosonic non-convergence
paradox [@navascues_2013_paradox]; 2-RDM SDP walls on the Hubbard model
[@rubin_2013_v2dm]; convergence is fast only for the targeted observable
[@berenstein_2025_goldilocks]; level-dependent ceilings for Quantum Max Cut
[@takahashi_2023_su2; @watts_2024_relaxations]; product-state approximation
limits [@brandao_2013_productstate]; floating-point solver output is not a
rigorous certificate without rational post-processing
[@naceur_2025_certified].

## Key open problems

- **Tightening 2D bounds at scale.** The 2026 structured NPA reaches 16×16
  but gaps vs QMC grow with L (0.00% → 0.69%); better monomial bases for
  non-local observables are explicitly open.
- **Long-range observables.** Hamiltonian-local bases leave C(L/2, L/2)-type
  bounds uselessly wide at large L; basis tailoring per observable.
- **Long-range order at the thermodynamic limit.** LRO certified only for
  L ≤ 8; finite-size certified bounds do not extrapolate.
- **Frustrated regimes.** Certification where no exact method exists is at
  few-percent looseness; solver stability needs care (PSD state-optimality
  destabilizes).
- **Rigorous end-to-end certificates.** Rational post-processing exists for
  1D chains; 2D has no packaged exact certification.
- **The 2026 "new column" reproduction gap.** Public QMBCertify reproduces
  the PRX-2024-era (old) column exactly; the 2026 extended reductions are
  not yet in the shipped code (found during beginner-training Track 2).

## Key bottlenecks

- **Memory, not time.** Mosek factorizations; 128 GB carried N=100 (1D) and
  L=10 (2D); 16×16 took 1 TB.
- **Solver trust.** Only clean OPTIMAL statuses certify; trustworthy digits
  are capped by residuals; rational rounding needed for frontier claims.
- **Basis design is manual.** Tightness per cost lives in the monomial
  basis; automate/optimize it (cf. moment-optimization framework
  [@flora_2026_moment]).
- **Benchmark asymmetry.** No leaderboard for *certified* bounds; validation
  is ad-hoc comparison against QMC/VMC tables.
