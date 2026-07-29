# H4 higher-dimensional SU(2) coarse-grainer contract

## Representation

The smallest odd-dimensional SU(2)-symmetric outer virtual space compatible with the two-site Heisenberg unit cell is `0⊕0⊕1`, with dense dimension five. The alternating inner virtual space is `1/2⊕1/2`, with dense dimension four. The two physical spin-1/2 tensors therefore have dense shapes `5×2×4` and `4×2×5`; their contraction is a real `5×4×5` super-site tensor.

The outer virtual U(1) charges in TensorKit dense order are `[0, 0, -2, 0, 2]`, and the two-spin super-site charges are `[-2, 0, 0, 2]`. The exactified tensor has zero forbidden-charge entries. `W2`, `MR`, `ML`, both Kronecker link maps, and `HLOC` preserve their derived charges exactly in the stored run.

The finite-size lower-bound argument is unchanged from [FINITE_SIZE_CONTRACT.md](FINITE_SIZE_CONTRACT.md): global U(1) twirling preserves the Heisenberg objective, and a charge-covariant coarse-grainer maps every twirled physical witness into the declared blocks. Changing the virtual representation changes the relaxation but not the certificate direction.

## Acceptance gates

1. VUMPS must converge with Galerkin norm at most `1e-12`, and its variational energy must remain above the Bethe energy.
2. At `n_super = 4`, dense and U(1)-blocked SDPs must both terminate `OPTIMAL`; their raw objectives must agree within `2e-7`, their residual-adjusted bounds must remain below the Bethe ceiling, and the largest PSD block must shrink by at least a factor of two.
3. Before a finite-size run, the D = 5 raw objective must improve over the validated D = 4 objective by more than both `2e-7` and the sum of their achieved solver-residual scales.
4. Before an N = 200 attempt, a D = 5 point at N = 100 must improve over the validated D = 4 N = 100 endpoint by more than both thresholds.

Every run consumes explicit predecessor paths, verifies stored objectives and residual diagnostics, and rejects stale snapshots. Residual-adjusted fields are not strict certificates; a published rigorous lower bound additionally requires [DUAL_POSTCERTIFICATION_CONTRACT.md](DUAL_POSTCERTIFICATION_CONTRACT.md).

## Validated results

The D = 5 VUMPS coarse-grainer has energy density `-0.440303008286` and Galerkin norm `7.32e-13`. Relative to the validated D = 4 map, the variational energy improves by `7.76e-6`.

At `n_super = 4`, dense and U(1) objectives differ by `4.39e-13`; both are `OPTIMAL`, with numerical downward diagnostic `-0.456387`. U(1) blocking reduces the 400-dimensional dense `omega` matrix to a largest block of 116, a `3.45×` reduction, and reduces measured total wall time from 95.0 s to 3.64 s in the matched run.

At `n_super = 6`, the D = 5 raw objective is `-0.448950917510`, compared with `-0.448972452237` for D = 4. The raw improvement is `2.153e-5`; this is a numerical map-selection result, not a strict bound comparison.

At finite periodic N = 100 and `n_super = 40`, the D = 5 solve is `OPTIMAL` with raw objective `-0.443636081671`. Its residual-adjusted value `-0.443636287511` and downward-rounded `-0.443637` are numerical diagnostics. The raw improvement over the D = 4 N = 100, `n_super = 50` endpoint is `1.851e-4`, so the numerical N = 200 accuracy screen passes, but a strict N = 100 finite gap is not claimed.

This is not the challenge target. The N = 100 point contains 622,464 scalar variables, takes 911.5 s, and live process monitoring observes an RSS snapshot of 43,763,876 KiB during the solve. Reusing its raw objective as an N = 200 numerical diagnostic would leave a gap of about `4.69e-4`, so a substantially deeper cluster run or a certificate-preserving coarse-grainer optimization is still required to reach `1e-5`.

The prepared `h4_finite_n200_d5.jl` entrypoint admits only `n_super = 50, 75, 100`. The n = 50 point consumes the numerical N = 100 gate; n = 75 and 100 additionally require the preceding N = 200 result. Every point uses a fresh directory, but the `1e-5` target can pass only after strict dual post-certification is added to that entrypoint.
