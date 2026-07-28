# T1 — VUMPS uniform MPS for the 1D Heisenberg AFM (the coarse-grainer for
# the Kull-Schuch compressed-LTI relaxation, challenge #49).
# Test-first: the variational energy must approach e0 = 1/4 - ln2 from ABOVE
# and tighten with D. Saves the MPS tensors for T2 via Serialization.

using MPSKit, MPSKitModels, TensorKit, LinearAlgebra, Serialization, Printf, JSON

const E0_EXACT = 1/4 - log(2)   # Bethe ansatz, H = sum S_i . S_{i+1}, J = 1
const DEFAULT_DS = [2, 3, 4, 5, 6, 7]
const DEFAULT_TOL = 1e-12
const DEFAULT_MAXITER = 400

parse_int_list(spec::AbstractString) = [parse(Int, strip(x)) for x in split(spec, ",") if !isempty(strip(x))]

function parse_args(args)
    Ds = copy(DEFAULT_DS)
    tol = DEFAULT_TOL
    maxiter = DEFAULT_MAXITER
    for arg in args
        if startswith(arg, "--Ds=")
            Ds = parse_int_list(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--tol=")
            tol = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--maxiter=")
            maxiter = parse(Int, split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    sort!(unique!(Ds))
    Ds, tol, maxiter
end

function write_outputs(psis, summary)
    serialize(joinpath(@__DIR__, "vumps_mps.jls"), psis)
    open(joinpath(@__DIR__, "vumps_summary.json"), "w") do io
        JSON.print(io, summary, 2)
    end
end

function main(args)
    Ds, tol, maxiter = parse_args(args)
    L = 2   # 2-site unit cell: AFM Néel period; VUMPS is far better conditioned than 1-site
    H = heisenberg_XXX(Float64, Trivial, InfiniteChain(L); spin=1//2)
    psis_path = joinpath(@__DIR__, "vumps_mps.jls")
    summary_path = joinpath(@__DIR__, "vumps_summary.json")
    psis = isfile(psis_path) ? deserialize(psis_path) : Dict{Int, Any}()
    summary = isfile(summary_path) ? JSON.parsefile(summary_path) : Dict(
        "exact_energy" => E0_EXACT,
        "unit_cell" => L,
        "tol" => tol,
        "maxiter" => maxiter,
        "Ds_requested" => Ds,
        "runs" => Dict{String, Any}(),
    )
    summary["exact_energy"] = E0_EXACT
    summary["unit_cell"] = L
    summary["tol"] = tol
    summary["maxiter"] = maxiter
    summary["Ds_requested"] = Ds
    summary["runs"] = get(summary, "runs", Dict{String, Any}())

    for D in Ds
        println("Running VUMPS for D=$D"); flush(stdout)
        psi0 = InfiniteMPS(Float64, fill(physicalspace(H, 1), L), fill(ComplexSpace(D), L))
        ok = false
        record = Dict{String, Any}("D" => D)
        delete!(psis, D)
        try
            psi, envs, delta = find_groundstate(psi0, H, VUMPS(; tol=tol, maxiter=maxiter))
            E = real(expectation_value(psi, H, envs) / L)   # energy density per site
            gap = E - E0_EXACT
            converged = delta <= tol
            record["energy"] = E
            record["gap"] = gap
            record["bnorm"] = delta
            record["converged"] = converged
            @printf("D = %d:  E/N = %.10f   (exact %.10f, variational gap %.2e, ||B|| = %.1e, converged=%s)\n",
                    D, E, E0_EXACT, gap, delta, converged)
            flush(stdout)
            @assert E > E0_EXACT - 1e-12 "variational principle violated at D=$D"
            if converged
                psis[D] = psi
                ok = true
            else
                println("   WARNING: VUMPS did not meet tol at D=$D; skipping this coarse-grainer")
                flush(stdout)
            end
        catch err
            record["error"] = sprint(showerror, err)
            println("   WARNING: VUMPS failed at D=$D: $(record["error"])")
            flush(stdout)
        end
        record["saved"] = ok
        summary["runs"]["D$D"] = record
        write_outputs(psis, summary)
    end

    println("Saved converged MPS tensors to $psis_path"); flush(stdout)
    println("Saved VUMPS summary to $summary_path"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
