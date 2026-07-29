# H4 — smallest D=5 SU(2)-symmetric coarse-grainer and dense/U(1) correctness gate.

include(joinpath(@__DIR__, "symmetry_probe.jl"))

const H4_D = 5
const H4_N_SUPER = 4
const H4_MOSEK_THREADS = 16
const H4_SEED = 20260729
const H4_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const H4_COARSE_GRAINER_PATH = joinpath(H4_REPO_ROOT, "artifacts", "coarse_grainers", "heisenberg1d", "symmetry_probe_d5.jls")

function h4_input_path(spec)
    normpath(isabspath(spec) ? spec : joinpath(H4_REPO_ROOT, spec))
end

function h4_recorded_input(result, key)
    inputs = get(get(result, "meta", Dict{String, Any}()), "inputs", Dict{String, Any}())
    haskey(inputs, key) || error("H4 predecessor does not record input $key")
    path = h4_input_path(inputs[key])
    isfile(path) || error("recorded H4 input does not exist: $path")
    path
end

function parse_h4_args(args)
    run_name = nothing
    mosek_threads = H4_MOSEK_THREADS
    d4_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d4-result=")
            d4_result = h4_input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h4-su2-d5-n4")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(d4_result) && error("--d4-result is required")
    isfile(d4_result) || error("D=4 predecessor does not exist: $d4_result")
    run_name, mosek_threads, d4_result
end

function h4_result_paths(run_name)
    result_dir = joinpath(H4_REPO_ROOT, "results", run_name)
    result_dir, joinpath(result_dir, "symmetry_probe_d5.json"), joinpath(result_dir, "figs", "symmetry_probe_d5.png")
end

function build_d5_coarse_grainer()
    Random.seed!(H4_SEED)
    H = heisenberg_XXX(Float64, SU2Irrep, InfiniteChain(2); spin=1 // 2)
    physical = Rep[SU₂](1 // 2 => 1)
    half_integer_bond = Rep[SU₂](1 // 2 => 2)
    integer_bond = Rep[SU₂](0 => 2, 1 => 1)
    psi0 = InfiniteMPS(Float64, [physical, physical], [half_integer_bond, integer_bond])
    t0 = time_ns()
    psi, envs, delta = find_groundstate(psi0, H, VUMPS(; tol=H1_TOL, maxiter=H1_MAXITER, verbosity=1))
    wall = elapsed_seconds(t0)
    bnorm = norm(delta)
    energy = real(sum(expectation_value(psi, H, envs)) / 2)
    @assert bnorm <= H1_TOL "D=5 VUMPS did not converge: ||B||=$bnorm"
    @assert energy >= E0 - H1_TOL "D=5 variational energy fell below the Bethe value"
    B = super_tensor(psi)
    @assert size(B) == (H4_D, 4, H4_D)
    B, Dict(
        "algorithm" => "VUMPS",
        "tensor_symmetry" => "SU(2)",
        "block_symmetry" => "U(1) subgroup",
        "outer_bond_representation" => "0=>2,1=>1",
        "inner_bond_representation" => "1/2=>2",
        "outer_bond_dimension" => H4_D,
        "inner_bond_dimension" => 4,
        "energy_density" => energy,
        "bethe_gap" => energy - E0,
        "galerkin_norm" => bnorm,
        "tolerance" => H1_TOL,
        "maxiter" => H1_MAXITER,
        "seed" => H4_SEED,
        "wall_seconds" => wall,
        "converged" => true,
        "julia_version" => string(VERSION),
    )
end

function load_or_build_d5_coarse_grainer()
    if isfile(H4_COARSE_GRAINER_PATH)
        B, metadata = deserialize(H4_COARSE_GRAINER_PATH)
        size(B) == (H4_D, 4, H4_D) || error("cached D=5 coarse-grainer has the wrong shape")
        metadata["outer_bond_dimension"] == H4_D || error("cached D=5 coarse-grainer metadata is inconsistent")
        return B, metadata, true
    end
    B, metadata = build_d5_coarse_grainer()
    mkpath(dirname(H4_COARSE_GRAINER_PATH))
    mktemp(dirname(H4_COARSE_GRAINER_PATH)) do tmp, io
        serialize(io, (B, metadata))
        close(io)
        mv(tmp, H4_COARSE_GRAINER_PATH; force=true)
    end
    B, metadata, false
end

function validate_d4_probe(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H1" || error("D=4 predecessor is not an H1 result")
    meta["bond_dimension"] == H1_D || error("D=4 predecessor has the wrong bond dimension")
    meta["n_super"] == H4_N_SUPER || error("D=4 predecessor has the wrong hierarchy depth")
    result["vumps"]["converged"] || error("D=4 predecessor VUMPS did not converge")
    result["vumps"]["galerkin_norm"] <= result["vumps"]["tolerance"] || error("D=4 predecessor VUMPS residual exceeds tolerance")
    validated_solver_bounds(result["solves"]["dense"], "D=4 dense predecessor")
    validated_solver_bounds(result["solves"]["u1_blocks"], "D=4 U(1) predecessor")
    result
end

function plot_h4_probe(results, figure_path)
    omega_inventory = results["block_inventory"]["omega_u1"]
    charges = sort(parse.(Int, collect(keys(omega_inventory))))
    dimensions = [omega_inventory[string(charge)] for charge in charges]
    dense_dimension = only(results["block_inventory"]["omega_dense"])
    size_plot = bar(string.(charges), dimensions; xlabel="total charge 2Sz", ylabel="PSD dimension", title="D=5 U(1) omega blocks (dense: $dense_dimension)", legend=false, color=:steelblue, left_margin=5mm, bottom_margin=4mm)
    times = [results["solves"]["dense"]["total_wall_seconds"], results["solves"]["u1_blocks"]["total_wall_seconds"]]
    time_plot = bar(["dense", "U(1) blocks"], times; ylabel="total wall time (s)", title="D=5, n=4 after warm-up", legend=false, color=[:gray50, :darkorange], ylims=(0, 1.08 * maximum(times)), bottom_margin=4mm)
    mkpath(dirname(figure_path))
    savefig(plot(size_plot, time_plot; layout=(1, 2), size=(1050, 460)), figure_path)
end

function h4_main(args)
    run_name, mosek_threads, d4_result_path = parse_h4_args(args)
    result_dir, result_path, figure_path = h4_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    d4 = validate_d4_probe(d4_result_path)
    virtual_charges = [0, 0, -2, 0, 2]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_d5_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))

    rho_inventory = block_inventory(rho3_charges(super_charges))
    omega_inventory = block_inventory(omega_charges(super_charges, virtual_charges))
    dense_omega_dimension = 16 * H4_D^2
    largest_omega = maximum(values(omega_inventory))
    block_reduction = dense_omega_dimension / largest_omega
    @printf("D=5 rho3 blocks: %s\n", rho_inventory)
    @printf("D=5 omega blocks: %s (dense %d, largest %d)\n", omega_inventory, dense_omega_dimension, largest_omega)
    flush(stdout)

    warmup = warm_up(mosek_threads)
    GC.gc()
    dense = solve_formulation(B, H4_N_SUPER, virtual_charges, super_charges; blocked=false, mosek_threads)
    GC.gc()
    blocked = solve_formulation(B, H4_N_SUPER, virtual_charges, super_charges; blocked=true, mosek_threads)

    both_optimal = dense["status"] == "OPTIMAL" && blocked["status"] == "OPTIMAL"
    both_optimal || error("D=5 dense and U(1) solves must both terminate OPTIMAL")
    validated_solver_bounds(dense, "D=5 dense solve")
    validated_solver_bounds(blocked, "D=5 U(1) solve")
    energy_difference = abs(dense["energy"] - blocked["energy"])
    correctness = energy_difference <= 2e-7 && dense["residual_adjusted_lower_bound"] < E0 && blocked["residual_adjusted_lower_bound"] < E0
    correctness || error("D=5 dense/U(1) energy comparison failed")
    block_reduction >= 2 || error("D=5 largest PSD block did not shrink by a factor of two")
    wall_nonregression = blocked["total_wall_seconds"] <= dense["total_wall_seconds"]
    variational_improvement = d4["vumps"]["energy_density"] - vumps["energy_density"]

    results = Dict(
        "meta" => Dict(
            "stage" => "H4",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 2 * H4_N_SUPER,
            "n_super" => H4_N_SUPER,
            "bond_dimension" => H4_D,
            "energy_precision_digits" => 6,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "charge_units" => "2Sz",
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
            "inputs" => Dict("d4_result" => relpath(d4_result_path, H4_REPO_ROOT)),
        ),
        "vumps" => vumps,
        "D4_predecessor" => Dict(
            "source" => relpath(d4_result_path, H4_REPO_ROOT),
            "coarse_grainer_fingerprint" => d4["meta"]["coarse_grainer_fingerprint"],
            "energy_density" => d4["vumps"]["energy_density"],
            "galerkin_norm" => d4["vumps"]["galerkin_norm"],
        ),
        "charge_check" => Dict("largest_forbidden_entry" => largest_forbidden, "exactified_to_zero" => true, "map_violations" => map_violations),
        "block_inventory" => Dict(
            "rho3_dense" => [64],
            "rho3_u1" => rho_inventory,
            "omega_dense" => [dense_omega_dimension],
            "omega_u1" => omega_inventory,
            "omega_variable_count" => H4_N_SUPER - 3,
            "largest_rho_reduction" => 64 / maximum(values(rho_inventory)),
            "largest_omega_reduction" => block_reduction,
        ),
        "warmup" => warmup,
        "solves" => Dict("dense" => dense, "u1_blocks" => blocked),
        "gates" => Dict(
            "both_optimal" => both_optimal,
            "energy_difference" => energy_difference,
            "energy_agreement_atol" => 2e-7,
            "correctness_passed" => correctness,
            "largest_block_reduction" => block_reduction,
            "block_reduction_passed" => block_reduction >= 2,
            "wall_nonregression_passed" => wall_nonregression,
            "variational_improvement_over_D4" => variational_improvement,
            "variational_improvement_passed" => variational_improvement > 2e-7,
            "usefulness_passed" => block_reduction >= 2 && wall_nonregression && variational_improvement > 2e-7,
        ),
    )
    write_h1_json(result_path, results)
    plot_h4_probe(results, figure_path)
    @printf("D=5 VUMPS improvement over D=4: %.3e\n", variational_improvement)
    @printf("dense E %.9f, U(1) E %.9f, |Delta E| %.3e\n", dense["energy"], blocked["energy"], energy_difference)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h4_main(ARGS)
end
