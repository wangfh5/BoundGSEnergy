# L=4 setting sweep for the square-lattice Heisenberg model: which knob lifts the SDP bound from the old column
# (-0.703051) toward the paper's new Table 8 value (-0.701783)?
# Baseline = official 2D example (d=4, rdm=0, pso=0, lso=0, extra=0, SU2 off).

using QMBCertify
using Printf

@eval QMBCertify function resort(Locb, coef)
    nLocb = copy(Locb); sort!(nLocb); unique!(nLocb)
    ncoef = zeros(typeof(coef[1]), length(nLocb))
    for i = 1:length(Locb)
        ind = bfind(nLocb, Locb[i])
        ncoef[ind] += coef[i]
    end
    ind = ncoef .!= 0
    return nLocb[ind], ncoef[ind]
end

supp = [[1, 4]]
coe  = [3/2]
L, d = 4, 4
TARGET = -0.701783
OLD    = -0.703051

configs = [
    ("extra=2",                     Dict(:rdm=>0, :pso=>0, :lso=>0, :extra=>2, :SU2=>false)),
    ("rdm=8",                       Dict(:rdm=>8, :pso=>0, :lso=>0, :extra=>0, :SU2=>false)),
    ("lso=1",                       Dict(:rdm=>0, :pso=>0, :lso=>1, :extra=>0, :SU2=>false)),
    ("pso=3",                       Dict(:rdm=>0, :pso=>3, :lso=>0, :extra=>0, :SU2=>false)),
    ("lso=1 + pso=3 (GSB defaults)",Dict(:rdm=>0, :pso=>3, :lso=>1, :extra=>0, :SU2=>false)),
    ("rdm=8 + lso=1 + pso=3",       Dict(:rdm=>8, :pso=>3, :lso=>1, :extra=>0, :SU2=>false)),
]

println(@sprintf("%-28s %12s %10s %8s", "config", "E_SDP/N", "vs target", "wall(s)")); flush(stdout)
for (name, c) in configs
    t0 = time()
    try
        opt, _ = GSB(supp, coe, L, d; lattice="square", rdm=c[:rdm], pso=c[:pso],
                     lso=c[:lso], lol=L, extra=c[:extra], SU2_symmetry=c[:SU2],
                     three_type=[1,1], QUIET=true)
        wall = time() - t0
        println(@sprintf("%-28s %12.6f %+10.6f %8.1f", name, opt, opt - TARGET, wall)); flush(stdout)
    catch err
        println(@sprintf("%-28s   CRASH: %s", name, sprint(showerror, err) |> x -> first(split(x, "\n")))); flush(stdout)
    end
end
println("reference: old column = $OLD, target (new column) = $TARGET"); flush(stdout)
