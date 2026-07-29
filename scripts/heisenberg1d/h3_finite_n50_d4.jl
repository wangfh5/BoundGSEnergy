# H3 — finite periodic N=50 U(1)-blocked D=4 ladder with a Bethe reference.

include(joinpath(@__DIR__, "h3_finite_n20_d4.jl"))

const H3_N50 = 50
const H3_N50_DEPTHS = [15, 20, 25]

function parse_h3_n50_args(args)
    run_name = nothing
    h0_result = nothing
    n20_result = nothing
    legacy_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--h0-result=")
            h0_result = input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--n20-result=")
            n20_result = input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d3-result=")
            legacy_result = input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h3-finite-n50-d4")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(h0_result) && error("--h0-result is required")
    isnothing(n20_result) && error("--n20-result is required")
    isnothing(legacy_result) && error("--d3-result is required")
    isfile(h0_result) || error("H0 result does not exist: $h0_result")
    isfile(n20_result) || error("N=20 result does not exist: $n20_result")
    isfile(legacy_result) || error("D=3 result does not exist: $legacy_result")
    run_name, h0_result, n20_result, legacy_result
end

function h3_n50_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "finite_n50_d4.json"), joinpath(result_dir, "figs", "finite_n50_d4.png")
end

function validate_n50_chain(n50, h0, n20, legacy)
    expected_reference = bethe_ground_reference(H3_N50)
    isapprox(n50["exact_reference"]["energy_per_site"], expected_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("N=50 exact reference does not match a fresh Bethe calculation")
    reference_validation = n50["reference_validation"]
    expected_n20_reference = bethe_ground_reference(H3_N)
    isapprox(reference_validation["N20_bethe"]["energy_per_site"], expected_n20_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("N=50 embedded N=20 Bethe reference does not match a fresh calculation")
    reference_validation["N20_xdiag_energy_per_site"] == h0["exact_reference"]["energy_per_site"] || error("N=50 exact-reference snapshot does not match H0")
    expected_difference = abs(reference_validation["N20_bethe"]["energy_per_site"] - h0["exact_reference"]["energy_per_site"])
    reference_validation["absolute_difference"] == expected_difference || error("N=50 reference cross-check difference is inconsistent")
    reference_validation["passed"] && expected_difference <= 2e-12 || error("N=50 reference-validation gate did not pass")

    comparison = n50["comparison"]
    previous_point = last(n20["points"])
    prior_snapshot = comparison["D4_n10_N20_provenance"]
    for key in ("support_sites", "n_super", "energy_per_site", "residual_adjusted_lower_bound")
        prior_snapshot[key] == previous_point[key] || error("N=50 prior-point snapshot field $key does not match N=20")
    end
    prior_snapshot["source"] == n50["meta"]["inputs"]["n20_result"] || error("N=50 prior-point source does not match its recorded input")
    n50["meta"]["coarse_grainer_fingerprint"] == n20["meta"]["coarse_grainer_fingerprint"] || error("N=50 coarse-grainer fingerprint does not match N=20")

    haskey(legacy, "D3-n20") || error("the D=3 comparison does not contain D3-n20")
    legacy_point = legacy["D3-n20"]
    legacy_point["status"] == "OPTIMAL" || error("the D=3 comparison is not OPTIMAL")
    isfinite(legacy_point["E"]) && isfinite(legacy_point["wall"]) || error("the D=3 comparison contains non-finite values")
    d3_snapshot = comparison["D3_n20"]
    d3_snapshot["source"] == n50["meta"]["inputs"]["d3_result"] || error("N=50 D=3 source does not match its recorded input")
    d3_snapshot["status"] == legacy_point["status"] || error("N=50 D=3 status snapshot does not match its source")
    d3_snapshot["energy_per_site"] == legacy_point["E"] || error("N=50 D=3 energy snapshot does not match its source")
    d3_snapshot["wall_seconds_reported"] == legacy_point["wall"] || error("N=50 D=3 timing snapshot does not match its source")
    n50
end

function plot_h3_n50(results, figure_path)
    exact_energy = results["exact_reference"]["energy_per_site"]
    d4 = vcat([results["comparison"]["D4_n10_N20_provenance"]], results["points"])
    d3 = results["comparison"]["D3_n20"]
    d4_bounds = [get(point, "residual_adjusted_lower_bound", point["energy_per_site"]) for point in d4]
    plot_result = plot(
        [point["support_sites"] for point in d4],
        exact_energy .- d4_bounds;
        marker=:diamond,
        yscale=:log10,
        xlabel="constraint support (physical sites)",
        ylabel="finite-N gap: exact energy - lower bound",
        title="N=50 periodic Heisenberg chain",
        label="D=4 U(1)",
        legend=:topright,
        size=(800, 500),
    )
    scatter!(
        plot_result,
        [d3["support_sites"]],
        [exact_energy - d3["energy_per_site"]];
        marker=:circle,
        markersize=7,
        label="D=3 dense",
    )
    mkpath(dirname(figure_path))
    savefig(plot_result, figure_path)
end

function h3_n50_main(args)
    run_name, h0_result, n20_result, legacy_result = parse_h3_n50_args(args)
    result_dir, result_path, figure_path = h3_n50_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    h0 = validate_h0_predecessor(JSON.parsefile(h0_result))
    n20_reference = bethe_ground_reference(H3_N)
    n20_cross_difference = abs(n20_reference["energy_per_site"] - h0["exact_reference"]["energy_per_site"])
    n20_cross_difference <= 2e-12 || error("Bethe reference does not reproduce the N=20 exact diagonalization")
    exact_reference = bethe_ground_reference(H3_N50)
    exact_energy = exact_reference["energy_per_site"]

    n20_d4 = validate_h3_predecessor(JSON.parsefile(n20_result), H3_N, H3_DEPTHS)
    assert_recorded_input(n20_d4, "h0_result", h0_result)
    h1_n6 = validate_h1_predecessor(JSON.parsefile(recorded_input_path(n20_d4, "h1_result")))
    validate_n20_chain(n20_d4, h0, h1_n6)
    previous_point = last(n20_d4["points"])
    previous_energy = previous_point["energy_per_site"]
    legacy = JSON.parsefile(legacy_result)
    haskey(legacy, "D3-n20") || error("the D=3 comparison does not contain D3-n20")
    d3_legacy = legacy["D3-n20"]
    d3_legacy["status"] == "OPTIMAL" || error("the D=3 comparison is not OPTIMAL")
    isfinite(d3_legacy["E"]) && isfinite(d3_legacy["wall"]) || error("the D=3 comparison contains non-finite values")
    d3_comparison = Dict(
        "source" => relpath(legacy_result, H3_REPO_ROOT),
        "system_scope" => "finite_periodic_chain",
        "system_size" => H3_N50,
        "support_sites" => 40,
        "n_super" => 20,
        "bond_dimension" => 3,
        "status" => d3_legacy["status"],
        "energy_per_site" => d3_legacy["E"],
        "wall_seconds_reported" => d3_legacy["wall"],
        "timing_comparable" => false,
        "timing_note" => "legacy total wall time lacks matched thread and hardware provenance",
    )

    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
    fingerprint == n20_d4["meta"]["coarse_grainer_fingerprint"] || error("D=4 coarse-grainer differs from the N=20 run")

    results = Dict(
        "meta" => Dict(
            "stage" => "H3",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H3_N50,
            "boundary" => "periodic",
            "contract" => "docs/FINITE_SIZE_CONTRACT.md",
            "bond_dimension" => H1_D,
            "block_symmetry" => "U(1)",
            "charge_units" => "2Sz",
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => H3_MOSEK_THREADS,
            "coarse_grainer_fingerprint" => fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
            "inputs" => Dict(
                "h0_result" => relpath(h0_result, H3_REPO_ROOT),
                "n20_result" => relpath(n20_result, H3_REPO_ROOT),
                "d3_result" => relpath(legacy_result, H3_REPO_ROOT),
            ),
        ),
        "reference_validation" => Dict(
            "N20_bethe" => n20_reference,
            "N20_xdiag_energy_per_site" => h0["exact_reference"]["energy_per_site"],
            "absolute_difference" => n20_cross_difference,
            "passed" => true,
        ),
        "exact_reference" => exact_reference,
        "coarse_grainer" => Dict("vumps" => vumps, "largest_forbidden_u1_entry" => largest_forbidden, "map_violations" => map_violations),
        "comparison" => Dict(
            "D4_n10_N20_provenance" => Dict(
                "source" => relpath(n20_result, H3_REPO_ROOT),
                "support_sites" => previous_point["support_sites"],
                "n_super" => previous_point["n_super"],
                "energy_per_site" => previous_point["energy_per_site"],
                "residual_adjusted_lower_bound" => get(previous_point, "residual_adjusted_lower_bound", previous_point["energy_per_site"]),
            ),
            "D3_n20" => d3_comparison,
        ),
        "points" => Dict{String, Any}[],
    )
    write_h1_json(result_path, results)

    for n_super in H3_N50_DEPTHS
        support_sites = 2 * n_super
        support_sites <= H3_N50 || error("support exceeds the finite system size")
        solve = solve_formulation(B, n_super, virtual_charges, super_charges; blocked=true, mosek_threads=H3_MOSEK_THREADS)
        solve["status"] == "OPTIMAL" || error("D=4 n=$n_super did not terminate OPTIMAL")
        energy = solve["energy"]
        energy <= exact_energy || error("raw finite-size ceiling violated at n=$n_super")
        adjusted = solve["residual_adjusted_lower_bound"]
        adjusted < exact_energy || error("residual-adjusted finite-size ceiling violated at n=$n_super")
        energy >= previous_energy - 2e-7 || error("fixed-D monotonicity violated at n=$n_super")
        published = solve["published_lower_bound"]
        point = Dict(
            "system_scope" => "finite_periodic_chain",
            "system_size" => H3_N50,
            "boundary" => "periodic",
            "support_sites" => support_sites,
            "n_super" => n_super,
            "bond_dimension" => H1_D,
            "energy_per_site" => energy,
            "residual_adjusted_lower_bound" => adjusted,
            "published_lower_bound" => published,
            "raw_finite_gap" => exact_energy - energy,
            "finite_gap" => exact_energy - published,
            "solve" => solve,
        )
        push!(results["points"], point)
        previous_energy = energy
        write_h1_json(result_path, results)
        @printf("N=50 D=4 n=%d raw E/N %.9f, published lower bound %.6f, published gap %.3e\n", n_super, energy, published, exact_energy - published)
    end

    d4_n20 = only(point for point in results["points"] if point["n_super"] == 20)
    d4_endpoint = last(results["points"])
    results["gates"] = Dict(
        "all_optimal" => all(point["solve"]["status"] == "OPTIMAL" for point in results["points"]),
        "all_raw_solver_objectives_below_exact" => all(point["energy_per_site"] <= exact_energy for point in results["points"]),
        "all_raw_objectives_within_solver_tolerance" => all(point["energy_per_site"] <= exact_energy + H1_FEAS_TOL for point in results["points"]),
        "all_residual_adjusted_bounds_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "all_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "fixed_D_monotone_from_n10" => true,
        "D4_n20_improvement_over_D3_n20" => d4_n20["energy_per_site"] - d3_comparison["energy_per_site"],
        "D4_endpoint_finite_gap" => d4_endpoint["finite_gap"],
        "timing_comparison_valid" => false,
        "N50_advantage_passed" => d4_n20["energy_per_site"] > d3_comparison["energy_per_site"] + 2e-7,
    )
    write_h1_json(result_path, results)
    plot_h3_n50(results, figure_path)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h3_n50_main(ARGS)
end
