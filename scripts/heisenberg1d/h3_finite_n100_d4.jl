# H3 — finite periodic N=100 U(1)-blocked D=4 ladder.

include(joinpath(@__DIR__, "h3_finite_n50_d4.jl"))

const H3_N100 = 100
const H3_N100_DEPTHS = [30, 40, 50]

function parse_h3_n100_args(args)
    run_name = nothing
    n50_result = nothing
    legacy_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--n50-result=")
            n50_result = input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--d3-result=")
            legacy_result = input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h3-finite-n100-d4")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(n50_result) && error("--n50-result is required")
    isnothing(legacy_result) && error("--d3-result is required")
    isfile(n50_result) || error("N=50 result does not exist: $n50_result")
    isfile(legacy_result) || error("D=3 result does not exist: $legacy_result")
    run_name, n50_result, legacy_result
end

function h3_n100_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "finite_n100_d4.json"), joinpath(result_dir, "figs", "finite_n100_d4.png")
end

function plot_h3_n100(results, figure_path)
    exact_energy = results["exact_reference"]["energy_per_site"]
    d4 = vcat([results["comparison"]["D4_n25_N50_provenance"]], results["points"])
    d3 = results["comparison"]["D3_n50"]
    d4_bounds = [get(point, "residual_adjusted_lower_bound", point["energy_per_site"]) for point in d4]
    d4_gaps = exact_energy .- d4_bounds
    d3_gap = exact_energy - d3["energy_per_site"]
    all_gaps = vcat(d4_gaps, d3_gap)
    all_supports = vcat([point["support_sites"] for point in d4], d3["support_sites"])
    plot_result = plot(
        [point["support_sites"] for point in d4],
        d4_gaps;
        marker=:diamond,
        yscale=:log10,
        xlabel="constraint support (physical sites)",
        ylabel="finite-N gap: exact energy - lower bound",
        title="N=100 periodic Heisenberg chain",
        label="D=4 U(1)",
        legend=:topright,
        size=(800, 500),
        xlims=(0.96 * minimum(all_supports), 1.08 * maximum(all_supports)),
        ylims=(minimum(all_gaps) / 1.2, 1.25 * maximum(all_gaps)),
        right_margin=7mm,
        top_margin=5mm,
    )
    scatter!(
        plot_result,
        [d3["support_sites"]],
        [d3_gap];
        marker=:circle,
        markersize=7,
        label="D=3 dense",
    )
    mkpath(dirname(figure_path))
    savefig(plot_result, figure_path)
end

function h3_n100_main(args)
    run_name, n50_result, legacy_result = parse_h3_n100_args(args)
    result_dir, result_path, figure_path = h3_n100_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    n50_d4 = validate_h3_predecessor(JSON.parsefile(n50_result), H3_N50, H3_N50_DEPTHS)
    assert_recorded_input(n50_d4, "d3_result", legacy_result)
    h0_result = recorded_input_path(n50_d4, "h0_result")
    n20_result = recorded_input_path(n50_d4, "n20_result")
    h0 = validate_h0_predecessor(JSON.parsefile(h0_result))
    n20_d4 = validate_h3_predecessor(JSON.parsefile(n20_result), H3_N, H3_DEPTHS)
    assert_recorded_input(n20_d4, "h0_result", h0_result)
    h1_n6 = validate_h1_predecessor(JSON.parsefile(recorded_input_path(n20_d4, "h1_result")))
    validate_n20_chain(n20_d4, h0, h1_n6)
    exact_reference = bethe_ground_reference(H3_N100)
    exact_energy = exact_reference["energy_per_site"]
    previous_point = last(n50_d4["points"])
    previous_energy = previous_point["energy_per_site"]

    legacy = JSON.parsefile(legacy_result)
    validate_n50_chain(n50_d4, h0, n20_d4, legacy)
    haskey(legacy, "D3-n50") || error("the D=3 comparison does not contain D3-n50")
    d3_legacy = legacy["D3-n50"]
    d3_legacy["status"] == "OPTIMAL" || error("the D=3 comparison is not OPTIMAL")
    isfinite(d3_legacy["E"]) && isfinite(d3_legacy["wall"]) || error("the D=3 comparison contains non-finite values")
    d3_comparison = Dict(
        "source" => relpath(legacy_result, H3_REPO_ROOT),
        "system_scope" => "finite_periodic_chain",
        "system_size" => H3_N100,
        "support_sites" => 100,
        "n_super" => 50,
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
    fingerprint == n50_d4["meta"]["coarse_grainer_fingerprint"] || error("D=4 coarse-grainer differs from the N=50 run")

    results = Dict(
        "meta" => Dict(
            "stage" => "H3",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H3_N100,
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
                "n50_result" => relpath(n50_result, H3_REPO_ROOT),
                "d3_result" => relpath(legacy_result, H3_REPO_ROOT),
            ),
        ),
        "exact_reference" => exact_reference,
        "coarse_grainer" => Dict("vumps" => vumps, "largest_forbidden_u1_entry" => largest_forbidden, "map_violations" => map_violations),
        "comparison" => Dict(
            "D4_n25_N50_provenance" => Dict(
                "source" => relpath(n50_result, H3_REPO_ROOT),
                "support_sites" => previous_point["support_sites"],
                "n_super" => previous_point["n_super"],
                "energy_per_site" => previous_point["energy_per_site"],
                "residual_adjusted_lower_bound" => get(previous_point, "residual_adjusted_lower_bound", previous_point["energy_per_site"]),
            ),
            "D3_n50" => d3_comparison,
        ),
        "points" => Dict{String, Any}[],
    )
    write_h1_json(result_path, results)

    for n_super in H3_N100_DEPTHS
        support_sites = 2 * n_super
        support_sites <= H3_N100 || error("support exceeds the finite system size")
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
            "system_size" => H3_N100,
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
        @printf("N=100 D=4 n=%d raw E/N %.9f, published lower bound %.6f, published gap %.3e\n", n_super, energy, published, exact_energy - published)
    end

    d4_endpoint = last(results["points"])
    results["gates"] = Dict(
        "all_optimal" => all(point["solve"]["status"] == "OPTIMAL" for point in results["points"]),
        "all_raw_solver_objectives_below_exact" => all(point["energy_per_site"] <= exact_energy for point in results["points"]),
        "all_raw_objectives_within_solver_tolerance" => all(point["energy_per_site"] <= exact_energy + H1_FEAS_TOL for point in results["points"]),
        "all_residual_adjusted_bounds_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "all_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "fixed_D_monotone_from_n25" => true,
        "D4_n50_improvement_over_D3_n50" => d4_endpoint["energy_per_site"] - d3_comparison["energy_per_site"],
        "D4_endpoint_finite_gap" => d4_endpoint["finite_gap"],
        "timing_comparison_valid" => false,
        "N100_advantage_passed" => d4_endpoint["energy_per_site"] > d3_comparison["energy_per_site"] + 2e-7,
    )
    write_h1_json(result_path, results)
    plot_h3_n100(results, figure_path)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h3_n100_main(ARGS)
end
