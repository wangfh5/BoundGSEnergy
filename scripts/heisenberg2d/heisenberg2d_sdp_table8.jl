# Reproduction of arXiv:2604.01555 Table 8 / Fig. 7 (beginner tier)
# Certified SDP lower bounds on E0/N, square-lattice AFM Heisenberg (PBC), L = 4, 6, 8.
# Settings follow the paper's official 2D example: d = 4, rdm/pso/lso = 0, extra = 0.
# Usage: julia --project=. scripts/heisenberg2d/heisenberg2d_sdp_table8.jl <run-dir>

using QMBCertify
using JSON
using Plots
using Printf

# WORKAROUND for an upstream bug (QMBCertify master, commits 70c9a4f/4d82be4b):
# `resort` is called in src/basic_function.jl (eigen_circmat, lines 316/341) but
# its definition was lost when src/sdp.jl was split into the current file layout.
# Body ported verbatim from the pre-split sdp.jl; bfind adapted to its current
# 2-argument form (same binary-search semantics). @eval defines it inside the
# module (plain `function QMBCertify.resort` is rejected for new names).
@eval QMBCertify function resort(Locb, coef)
    nLocb = copy(Locb)
    sort!(nLocb)
    unique!(nLocb)
    ncoef = zeros(typeof(coef[1]), length(nLocb))
    for i = 1:length(Locb)
        ind = bfind(nLocb, Locb[i])
        ncoef[ind] += coef[i]
    end
    ind = ncoef .!= 0
    return nLocb[ind], ncoef[ind]
end

rundir = length(ARGS) >= 1 ? ARGS[1] : "results/20260727-183425-wang2026-table8-heisenberg2d"
mkpath(joinpath(rundir, "figs"))

# 2D square Heisenberg: nearest-neighbour bond, per-spin normalization (3/2 = 2 bonds/site x 3 components x 1/4)
supp = [[1, 4]]
coe  = [3/2]

# Paper Table 8 reference values (ED for L=4,6; QMC for L=8)
refs = Dict(4 => -0.701780, 6 => -0.678872, 8 => -0.673487)
paper_sdp = Dict(4 => -0.701783, 6 => -0.680886, 8 => -0.676370)

results = Dict{String, Any}()

for L in [4, 6, 8]
    println("=" ^ 60); flush(stdout)
    println("L = $L  ($(L*L) spins) — assembling structured SDP (d=4)..."); flush(stdout)
    t0 = time()
    opt, data = GSB(supp, coe, L, 4; lattice="square", rdm=0, pso=0, lso=0,
                    extra=0, QUIET=false)
    wall = time() - t0
    ref = refs[L]
    gap = (ref - opt) / abs(ref) * 100
    println(@sprintf("L = %d done: E_SDP/N = %.6f  (paper: %.6f, ref: %.6f, gap vs ref: %.2f%%, wall: %.1f s)",
                     L, opt, paper_sdp[L], ref, gap, wall)); flush(stdout)
    @assert opt <= ref + 1e-6 "sandwich violated at L=$L: bound $opt above reference $ref"
    results[string(L)] = Dict("E_SDP" => opt, "paper_SDP" => paper_sdp[L],
                              "E_ref" => ref, "gap_pct" => gap, "wall_s" => wall)
    # incremental write after each size — a kill loses at most one point
    open(joinpath(rundir, "sdp_bounds.json"), "w") do io
        JSON.print(io, results, 2)
    end
    println("checkpoint written for L = $L"); flush(stdout)
end

Ls = sort!(parse.(Int, collect(keys(results))))
bounds = [results[string(l)]["E_SDP"] for l in Ls]
refvals = [results[string(l)]["E_ref"] for l in Ls]
papervals = [results[string(l)]["paper_SDP"] for l in Ls]

p = plot(title="Certified lower bounds on E0/N — square-lattice Heisenberg (PBC)",
         xlabel="L (LxL lattice)", ylabel="E / N", legend=:bottomleft,
         size=(700, 450))
plot!(p, Ls, refvals, marker=:circle, label="reference (ED L<=6, QMC L=8)")
plot!(p, Ls, papervals, marker=:square, ls=:dash, label="paper SDP (Table 8)")
plot!(p, Ls, bounds, marker=:diamond, label="this run SDP (d=4)")
savefig(p, joinpath(rundir, "figs", "table8_bounds.png"))
println("figure saved; all sizes complete."); flush(stdout)
