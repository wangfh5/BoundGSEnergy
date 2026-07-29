# H4 — residual-resolved D=5 versus D=4 certificate comparison at n=6.

include(joinpath(@__DIR__, "h4_symmetry_probe_d5.jl"))

function parse_h4_n6_args(args)
    run_name = nothing
    mosek_threads = H4_MOSEK_THREADS
    d5_n4_result = nothing
    d4_n6_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d5-n4-result=")
            d5_n4_result = h4_input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d4-n6-result=")
            d4_n6_result = h4_input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h4-su2-d5-n6")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(d5_n4_result) && error("--d5-n4-result is required")
    isnothing(d4_n6_result) && error("--d4-n6-result is required")
    isfile(d5_n4_result) || error("D=5 n=4 result does not exist: $d5_n4_result")
    isfile(d4_n6_result) || error("D=4 n=6 result does not exist: $d4_n6_result")
    run_name, mosek_threads, d5_n4_result, d4_n6_result
end

function validate_d5_n4_probe(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H4" || error("D=5 n=4 predecessor is not an H4 result")
    meta["bond_dimension"] == H4_D || error("D=5 n=4 predecessor has the wrong bond dimension")
    meta["n_super"] == H4_N_SUPER || error("D=5 n=4 predecessor has the wrong hierarchy depth")
    result["gates"]["correctness_passed"] || error("D=5 dense/U(1) correctness gate did not pass")
    validated_solver_bounds(result["solves"]["dense"], "D=5 n=4 dense predecessor")
    validated_solver_bounds(result["solves"]["u1_blocks"], "D=5 n=4 U(1) predecessor")
    d4_path = h4_recorded_input(result, "d4_result")
    d4 = validate_d4_probe(d4_path)
    snapshot = result["D4_predecessor"]
    snapshot["source"] == result["meta"]["inputs"]["d4_result"] || error("D=5 n=4 D=4 source snapshot is stale")
    snapshot["coarse_grainer_fingerprint"] == d4["meta"]["coarse_grainer_fingerprint"] || error("D=5 n=4 D=4 fingerprint snapshot is stale")
    snapshot["energy_density"] == d4["vumps"]["energy_density"] || error("D=5 n=4 D=4 energy snapshot is stale")
    result["vumps"]["converged"] && result["vumps"]["galerkin_norm"] <= result["vumps"]["tolerance"] || error("D=5 n=4 VUMPS gate did not pass")
    maximum(values(result["charge_check"]["map_violations"])) <= H1_TOL || error("D=5 n=4 charge checks did not pass")
    result
end

function validate_d4_n6_probe(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H1" || error("D=4 n=6 predecessor is not an H1 result")
    meta["bond_dimension"] == H1_D || error("D=4 n=6 predecessor has the wrong bond dimension")
    meta["n_super"] == 6 || error("D=4 predecessor has the wrong hierarchy depth")
    validated_solver_bounds(result["solves"]["u1_blocks"], "D=4 n=6 U(1) predecessor")
    result
end

function plot_h4_n6(results, figure_path)
    gates = results["gates"]
    values = [gates["raw_improvement_over_D4"], gates["largest_solver_residual_scale"], gates["minimum_improvement"]]
    labels = ["raw D=5 improvement", "largest solver residual", "acceptance threshold"]
    scale_plot = scatter(1:3, values; yscale=:log10, xticks=(1:3, labels), ylabel="energy-density scale", title="D=5 certificate signal at n=6", label=false, markersize=8, markercolor=[:steelblue, :darkorange, :black], ylims=(minimum(values) / 2, maximum(values) * 2), left_margin=6mm, bottom_margin=12mm)
    bounds = [results["D4_predecessor"]["published_lower_bound"], results["D5"]["published_lower_bound"]]
    bound_plot = scatter(1:2, bounds; xticks=(1:2, ["D=4", "D=5"]), ylabel="published lower bound", title="Conservative n=6 bounds", label=false, markersize=9, markercolor=[:gray55, :steelblue], ylims=(minimum(bounds) - 5e-6, maximum(bounds) + 5e-6), left_margin=7mm, bottom_margin=5mm)
    mkpath(dirname(figure_path))
    savefig(plot(scale_plot, bound_plot; layout=(1, 2), size=(1100, 470)), figure_path)
end

function h4_n6_main(args)
    run_name, mosek_threads, d5_n4_path, d4_n6_path = parse_h4_n6_args(args)
    result_dir = joinpath(H4_REPO_ROOT, "results", run_name)
    result_path = joinpath(result_dir, "d5_n6.json")
    figure_path = joinpath(result_dir, "figs", "d5_n6.png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    d5_n4 = validate_d5_n4_probe(d5_n4_path)
    d4_n6 = validate_d4_n6_probe(d4_n6_path)
    virtual_charges = [0, 0, -2, 0, 2]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_d5_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
    fingerprint == d5_n4["meta"]["coarse_grainer_fingerprint"] || error("cached D=5 coarse-grainer does not match the n=4 gate")

    d5 = solve_formulation(B, 6, virtual_charges, super_charges; blocked=true, mosek_threads)
    validated_solver_bounds(d5, "D=5 n=6 U(1) solve")
    d5["residual_adjusted_lower_bound"] < E0 || error("D=5 n=6 residual-adjusted bound violates the Bethe ceiling")

    d4 = d4_n6["solves"]["u1_blocks"]
    raw_improvement = d5["energy"] - d4["energy"]
    conservative_improvement = d5["residual_adjusted_lower_bound"] - d4["residual_adjusted_lower_bound"]
    largest_residual = max(solver_residual_scale(d4), solver_residual_scale(d5))
    minimum_improvement = 2e-7
    residual_resolved = raw_improvement > max(minimum_improvement, largest_residual)
    residual_resolved || error("D=5 n=6 improvement is not residual-resolved")

    results = Dict(
        "meta" => Dict(
            "stage" => "H4",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 12,
            "n_super" => 6,
            "bond_dimension" => H4_D,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "charge_units" => "2Sz",
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
            "inputs" => Dict(
                "d5_n4_result" => relpath(d5_n4_path, H4_REPO_ROOT),
                "d4_n6_result" => relpath(d4_n6_path, H4_REPO_ROOT),
            ),
        ),
        "vumps" => vumps,
        "charge_check" => Dict("largest_forbidden_entry" => largest_forbidden, "map_violations" => map_violations),
        "D4_predecessor" => Dict(
            "source" => relpath(d4_n6_path, H4_REPO_ROOT),
            "coarse_grainer_fingerprint" => d4_n6["meta"]["coarse_grainer_fingerprint"],
            "energy" => d4["energy"],
            "residual_adjusted_lower_bound" => d4["residual_adjusted_lower_bound"],
            "published_lower_bound" => d4["published_lower_bound"],
            "residual_scale" => solver_residual_scale(d4),
        ),
        "D5" => d5,
        "gates" => Dict(
            "minimum_improvement" => minimum_improvement,
            "raw_improvement_over_D4" => raw_improvement,
            "conservative_improvement_over_D4" => conservative_improvement,
            "largest_solver_residual_scale" => largest_residual,
            "improvement_to_residual_ratio" => raw_improvement / largest_residual,
            "residual_resolved_improvement_passed" => residual_resolved,
            "timing_comparison_valid" => false,
        ),
    )
    write_h1_json(result_path, results)
    plot_h4_n6(results, figure_path)
    @printf("D=5 n=6 raw improvement %.3e, conservative improvement %.3e, residual ratio %.1fx\n", raw_improvement, conservative_improvement, raw_improvement / largest_residual)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h4_n6_main(ARGS)
end
