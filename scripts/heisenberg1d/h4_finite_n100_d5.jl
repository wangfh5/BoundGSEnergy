# H4 — finite N=100 D=5 probe against the validated D=4 endpoint.

include(joinpath(@__DIR__, "h4_d5_n6.jl"))

const H4_N100 = 100

function parse_h4_n100_args(args)
    run_name = nothing
    n_super = 40
    mosek_threads = H4_MOSEK_THREADS
    d5_n6_result = nothing
    d4_n100_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--n-super=")
            n_super = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d5-n6-result=")
            d5_n6_result = h4_input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d4-n100-result=")
            d4_n100_result = h4_input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h4-finite-n100-d5-n$(n_super)")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    n_super in (30, 40, 50) || error("n_super must match a validated D=4 depth: 30, 40, or 50")
    2 * n_super <= H4_N100 || error("constraint support exceeds N=100")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(d5_n6_result) && error("--d5-n6-result is required")
    isnothing(d4_n100_result) && error("--d4-n100-result is required")
    isfile(d5_n6_result) || error("D=5 n=6 result does not exist: $d5_n6_result")
    isfile(d4_n100_result) || error("D=4 N=100 result does not exist: $d4_n100_result")
    run_name, n_super, mosek_threads, d5_n6_result, d4_n100_result
end

function validate_d5_n6_result(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H4" || error("D=5 n=6 predecessor is not an H4 result")
    meta["bond_dimension"] == H4_D || error("D=5 n=6 predecessor has the wrong bond dimension")
    meta["n_super"] == 6 || error("D=5 predecessor has the wrong hierarchy depth")
    result["gates"]["residual_resolved_improvement_passed"] || error("D=5 n=6 accuracy gate did not pass")
    validated_solver_bounds(result["D5"], "D=5 n=6 predecessor")
    d5_n4 = validate_d5_n4_probe(h4_recorded_input(result, "d5_n4_result"))
    d4_n6 = validate_d4_n6_probe(h4_recorded_input(result, "d4_n6_result"))
    meta["coarse_grainer_fingerprint"] == d5_n4["meta"]["coarse_grainer_fingerprint"] || error("D=5 n=6 fingerprint does not match the n=4 gate")
    d4 = d4_n6["solves"]["u1_blocks"]
    snapshot = result["D4_predecessor"]
    snapshot["source"] == result["meta"]["inputs"]["d4_n6_result"] || error("D=5 n=6 D=4 source snapshot is stale")
    snapshot["energy"] == d4["energy"] || error("D=5 n=6 D=4 energy snapshot is stale")
    snapshot["residual_adjusted_lower_bound"] == d4["residual_adjusted_lower_bound"] || error("D=5 n=6 D=4 bound snapshot is stale")
    raw_improvement = result["D5"]["energy"] - d4["energy"]
    conservative_improvement = result["D5"]["residual_adjusted_lower_bound"] - d4["residual_adjusted_lower_bound"]
    largest_residual = max(solver_residual_scale(result["D5"]), solver_residual_scale(d4))
    gates = result["gates"]
    gates["raw_improvement_over_D4"] == raw_improvement || error("D=5 n=6 raw-improvement gate is stale")
    gates["conservative_improvement_over_D4"] == conservative_improvement || error("D=5 n=6 conservative-improvement gate is stale")
    gates["largest_solver_residual_scale"] == largest_residual || error("D=5 n=6 residual gate is stale")
    raw_improvement > max(gates["minimum_improvement"], largest_residual) || error("D=5 n=6 improvement is not residual-resolved")
    result
end

function validate_d4_n100_result(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H3" || error("D=4 N=100 predecessor is not an H3 result")
    meta["system_scope"] == "finite_periodic_chain" || error("D=4 predecessor has the wrong system scope")
    meta["system_size"] == H4_N100 || error("D=4 predecessor has the wrong system size")
    meta["boundary"] == "periodic" || error("D=4 predecessor has the wrong boundary")
    meta["bond_dimension"] == H1_D || error("D=4 predecessor has the wrong bond dimension")
    [point["n_super"] for point in result["points"]] == [30, 40, 50] || error("D=4 predecessor has unexpected depths")
    expected_reference = bethe_ground_reference(H4_N100)
    isapprox(result["exact_reference"]["energy_per_site"], expected_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("D=4 predecessor has the wrong exact reference")
    exact_energy = result["exact_reference"]["energy_per_site"]
    for point in result["points"]
        adjusted, published = validated_solver_bounds(point["solve"], "D=4 N=100 point")
        point["energy_per_site"] == point["solve"]["energy"] || error("D=4 point energy snapshot is stale")
        point["residual_adjusted_lower_bound"] == adjusted || error("D=4 point residual-adjusted snapshot is stale")
        point["published_lower_bound"] == published || error("D=4 point published snapshot is stale")
        point["system_size"] == H4_N100 && point["boundary"] == "periodic" && point["bond_dimension"] == H1_D || error("D=4 point metadata is inconsistent")
        point["support_sites"] == 2 * point["n_super"] <= H4_N100 || error("D=4 point support is inconsistent")
        point["raw_finite_gap"] == exact_energy - point["energy_per_site"] || error("D=4 raw finite gap is stale")
        point["finite_gap"] == exact_energy - published || error("D=4 published finite gap is stale")
        published <= adjusted <= point["energy_per_site"] < exact_energy || error("D=4 point violates the finite ceiling")
    end
    result["gates"]["all_optimal"] && result["gates"]["all_below_exact"] || error("D=4 predecessor gates did not pass")
    result
end

function validate_h4_n100_result(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H4" || error("D=5 N=100 result is not an H4 result")
    meta["system_scope"] == "finite_periodic_chain" && meta["system_size"] == H4_N100 && meta["boundary"] == "periodic" || error("D=5 N=100 result metadata is inconsistent")
    meta["bond_dimension"] == H4_D || error("D=5 N=100 result has the wrong bond dimension")

    d5_n6 = validate_d5_n6_result(h4_recorded_input(result, "d5_n6_result"))
    d4_n100 = validate_d4_n100_result(h4_recorded_input(result, "d4_n100_result"))
    meta["coarse_grainer_fingerprint"] == d5_n6["meta"]["coarse_grainer_fingerprint"] || error("D=5 N=100 fingerprint does not match its predecessor")
    result["D4_predecessor"]["points"] == d4_n100["points"] || error("D=5 N=100 D=4 ladder snapshot is stale")

    exact_energy = result["exact_reference"]["energy_per_site"]
    expected_reference = bethe_ground_reference(H4_N100)
    isapprox(exact_energy, expected_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("D=5 N=100 exact reference is inconsistent")
    point = result["D5_point"]
    adjusted, published = validated_solver_bounds(point["solve"], "D=5 N=100 stored solve")
    point["energy_per_site"] == point["solve"]["energy"] || error("D=5 N=100 energy snapshot is stale")
    point["residual_adjusted_lower_bound"] == adjusted || error("D=5 N=100 residual-adjusted snapshot is stale")
    point["published_lower_bound"] == published || error("D=5 N=100 published snapshot is stale")
    point["n_super"] == meta["n_super"] && point["support_sites"] == meta["support_sites"] == 2 * meta["n_super"] || error("D=5 N=100 support metadata is inconsistent")
    point["raw_finite_gap"] == exact_energy - point["energy_per_site"] || error("D=5 N=100 raw finite gap is stale")
    point["finite_gap"] == exact_energy - published || error("D=5 N=100 published finite gap is stale")
    published <= adjusted <= point["energy_per_site"] < exact_energy || error("D=5 N=100 point violates the finite ceiling")

    endpoint = last(d4_n100["points"])
    matched = only(candidate for candidate in d4_n100["points"] if candidate["n_super"] == point["n_super"])
    raw_improvement = point["energy_per_site"] - endpoint["energy_per_site"]
    conservative_improvement = adjusted - endpoint["residual_adjusted_lower_bound"]
    largest_residual = max(solver_residual_scale(point["solve"]), solver_residual_scale(endpoint["solve"]))
    gates = result["gates"]
    gates["raw_improvement_over_D4_endpoint"] == raw_improvement || error("D=5 N=100 raw endpoint-improvement gate is stale")
    gates["conservative_improvement_over_D4_endpoint"] == conservative_improvement || error("D=5 N=100 conservative endpoint-improvement gate is stale")
    gates["raw_improvement_over_matched_D4_depth"] == point["energy_per_site"] - matched["energy_per_site"] || error("D=5 N=100 matched-depth gate is stale")
    gates["largest_solver_residual_scale"] == largest_residual || error("D=5 N=100 residual gate is stale")
    gates["N200_entry_passed"] == (raw_improvement > max(gates["minimum_improvement"], largest_residual)) || error("D=5 N=100 entry gate is stale")
    result
end

function plot_h4_n100(results, figure_path)
    exact_energy = results["exact_reference"]["energy_per_site"]
    d4_points = results["D4_predecessor"]["points"]
    d5_point = results["D5_point"]
    d4_gaps = [exact_energy - point["published_lower_bound"] for point in d4_points]
    gap_plot = plot([point["support_sites"] for point in d4_points], d4_gaps; marker=:circle, yscale=:log10, xlabel="constraint support (physical sites)", ylabel="residual-adjusted diagnostic gap", title="N=100 D=5 entry probe", label="D=4", legend=:topright, left_margin=10mm, bottom_margin=7mm)
    scatter!(gap_plot, [d5_point["support_sites"]], [d5_point["finite_gap"]]; marker=:diamond, markersize=8, label="D=5")

    gates = results["gates"]
    values = [abs(gates["raw_improvement_over_D4_endpoint"]), gates["largest_solver_residual_scale"], gates["minimum_improvement"]]
    labels = ["raw improvement", "solver residual", "2e-7 threshold"]
    signal_plot = scatter(1:3, values; yscale=:log10, xticks=(1:3, labels), ylabel="energy-density scale", title="N=200 entry signal", label=false, markersize=8, markercolor=[:steelblue, :darkorange, :black], ylims=(minimum(values) / 2, maximum(values) * 2), left_margin=8mm, right_margin=5mm, bottom_margin=10mm)

    mkpath(dirname(figure_path))
    savefig(plot(gap_plot, signal_plot; layout=(1, 2), size=(1250, 500)), figure_path)
end

function h4_n100_main(args)
    run_name, n_super, mosek_threads, d5_n6_path, d4_n100_path = parse_h4_n100_args(args)
    result_dir = joinpath(H4_REPO_ROOT, "results", run_name)
    result_path = joinpath(result_dir, "finite_n100_d5_n$(n_super).json")
    figure_path = joinpath(result_dir, "figs", "finite_n100_d5_n$(n_super).png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    d5_n6 = validate_d5_n6_result(d5_n6_path)
    d4_n100 = validate_d4_n100_result(d4_n100_path)
    exact_reference = bethe_ground_reference(H4_N100)
    exact_energy = exact_reference["energy_per_site"]
    virtual_charges = [0, 0, -2, 0, 2]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_d5_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
    fingerprint == d5_n6["meta"]["coarse_grainer_fingerprint"] || error("cached D=5 coarse-grainer does not match the n=6 gate")

    solve = solve_formulation(B, n_super, virtual_charges, super_charges; blocked=true, mosek_threads)
    adjusted, published = validated_solver_bounds(solve, "D=5 finite N=100 solve")
    energy = solve["energy"]
    energy < exact_energy || error("D=5 raw objective violates the finite N=100 ceiling")
    adjusted < exact_energy || error("D=5 residual-adjusted bound violates the finite N=100 ceiling")
    energy >= d5_n6["D5"]["energy"] - 2e-7 || error("D=5 hierarchy is not monotone from n=6")

    d4_endpoint = last(d4_n100["points"])
    d4_matched = only(point for point in d4_n100["points"] if point["n_super"] == n_super)
    raw_improvement = energy - d4_endpoint["energy_per_site"]
    conservative_improvement = adjusted - d4_endpoint["residual_adjusted_lower_bound"]
    matched_improvement = energy - d4_matched["energy_per_site"]
    largest_residual = max(solver_residual_scale(solve), solver_residual_scale(d4_endpoint["solve"]))
    minimum_improvement = 2e-7
    entry_passed = raw_improvement > max(minimum_improvement, largest_residual)

    d5_point = Dict(
        "system_scope" => "finite_periodic_chain",
        "system_size" => H4_N100,
        "boundary" => "periodic",
        "support_sites" => 2 * n_super,
        "n_super" => n_super,
        "bond_dimension" => H4_D,
        "energy_per_site" => energy,
        "residual_adjusted_lower_bound" => adjusted,
        "published_lower_bound" => published,
        "raw_finite_gap" => exact_energy - energy,
        "finite_gap" => exact_energy - published,
        "solve" => solve,
    )
    results = Dict(
        "meta" => Dict(
            "stage" => "H4",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H4_N100,
            "boundary" => "periodic",
            "contract" => "docs/FINITE_SIZE_CONTRACT.md",
            "support_sites" => 2 * n_super,
            "n_super" => n_super,
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
                "d5_n6_result" => relpath(d5_n6_path, H4_REPO_ROOT),
                "d4_n100_result" => relpath(d4_n100_path, H4_REPO_ROOT),
            ),
        ),
        "exact_reference" => exact_reference,
        "vumps" => vumps,
        "charge_check" => Dict("largest_forbidden_entry" => largest_forbidden, "map_violations" => map_violations),
        "D4_predecessor" => Dict(
            "source" => relpath(d4_n100_path, H4_REPO_ROOT),
            "coarse_grainer_fingerprint" => d4_n100["meta"]["coarse_grainer_fingerprint"],
            "points" => d4_n100["points"],
            "endpoint" => d4_endpoint,
        ),
        "D5_point" => d5_point,
        "gates" => Dict(
            "minimum_improvement" => minimum_improvement,
            "raw_improvement_over_D4_endpoint" => raw_improvement,
            "conservative_improvement_over_D4_endpoint" => conservative_improvement,
            "raw_improvement_over_matched_D4_depth" => matched_improvement,
            "largest_solver_residual_scale" => largest_residual,
            "improvement_to_residual_ratio" => raw_improvement / largest_residual,
            "N200_entry_passed" => entry_passed,
            "timing_comparison_valid" => false,
        ),
    )
    write_h1_json(result_path, results)
    validate_h4_n100_result(result_path)
    plot_h4_n100(results, figure_path)
    @printf("N=100 D=5 n=%d raw E %.9f, published %.6f, endpoint improvement %.3e, entry gate %s\n", n_super, energy, published, raw_improvement, entry_passed)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h4_n100_main(ARGS)
end
