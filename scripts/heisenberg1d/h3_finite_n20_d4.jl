# H3 — extend the finite periodic N=20 ladder with U(1)-blocked D=4 points.

include(joinpath(@__DIR__, "symmetry_probe.jl"))

const H3_N = 20
const H3_DEPTHS = [8, 10]
const H3_MOSEK_THREADS = 16
const H3_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

input_path(spec) = normpath(isabspath(spec) ? spec : joinpath(H3_REPO_ROOT, spec))

function parse_h3_args(args)
    run_name = nothing
    h0_result = nothing
    h1_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--h0-result=")
            h0_result = input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--h1-result=")
            h1_result = input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h3-finite-n20-d4")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(h0_result) && error("--h0-result is required")
    isnothing(h1_result) && error("--h1-result is required")
    isfile(h0_result) || error("H0 result does not exist: $h0_result")
    isfile(h1_result) || error("H1 result does not exist: $h1_result")
    run_name, h0_result, h1_result
end

function h3_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "finite_n20_d4.json"), joinpath(result_dir, "figs", "finite_n20_d4.png")
end

function require_metadata(result, key, expected, label)
    actual = get(get(result, "meta", Dict{String, Any}()), key, nothing)
    actual == expected || error("$label has meta.$key=$actual, expected $expected")
end

include(joinpath(@__DIR__, "predecessor_validation.jl"))

function recorded_input_path(result, key)
    inputs = get(get(result, "meta", Dict{String, Any}()), "inputs", Dict{String, Any}())
    haskey(inputs, key) || error("predecessor does not record input $key")
    path = input_path(inputs[key])
    isfile(path) || error("recorded predecessor input does not exist: $path")
    path
end

function assert_recorded_input(result, key, expected_path)
    recorded = recorded_input_path(result, key)
    realpath(recorded) == realpath(expected_path) || error("predecessor input $key does not match the supplied path")
    recorded
end

function validate_n20_chain(n20, h0, h1)
    n20["exact_reference"] == h0["exact_reference"] || error("N=20 exact-reference snapshot does not match H0")
    comparison = n20["comparison"]
    comparison["D3_ladder"] == h0["compressed"] || error("N=20 D=3 ladder snapshot does not match H0")
    comparison["D3_endpoint_energy_per_site"] == last(h0["compressed"])["energy_per_site"] || error("N=20 D=3 endpoint does not match H0")
    comparison["D4_n6_validated_energy_per_site"] == h1["solves"]["u1_blocks"]["energy"] || error("N=20 prior D=4 energy does not match H1")
    n20["meta"]["coarse_grainer_fingerprint"] == h1["meta"]["coarse_grainer_fingerprint"] || error("N=20 coarse-grainer fingerprint does not match H1")
    n20
end

function validate_h3_point_bound(point, exact_energy, system_size)
    solve = point["solve"]
    adjusted, published = validated_solver_bounds(solve, "H3 nested solve")
    point["energy_per_site"] == solve["energy"] || error("H3 outer energy does not match its solve")
    point["residual_adjusted_lower_bound"] == adjusted || error("H3 outer residual-adjusted bound does not match its solve")
    point["published_lower_bound"] == published || error("H3 outer published bound does not match its solve")
    point["raw_finite_gap"] == exact_energy - solve["energy"] || error("H3 raw finite gap is stale")
    point["finite_gap"] == exact_energy - published || error("H3 published finite gap is stale")
    point["system_scope"] == "finite_periodic_chain" &&
        point["system_size"] == system_size &&
        point["boundary"] == "periodic" &&
        point["bond_dimension"] == H1_D &&
        point["support_sites"] == 2 * point["n_super"] <= system_size &&
        published <= adjusted <= solve["energy"] < exact_energy
end

function validate_h3_predecessor(result, system_size, expected_depths)
    require_metadata(result, "stage", "H3", "H3 predecessor")
    require_metadata(result, "system_scope", "finite_periodic_chain", "H3 predecessor")
    require_metadata(result, "system_size", system_size, "H3 predecessor")
    require_metadata(result, "boundary", "periodic", "H3 predecessor")
    require_metadata(result, "bond_dimension", H1_D, "H3 predecessor")
    require_metadata(result, "numerical_bound_schema_version", NUMERICAL_BOUND_SCHEMA_VERSION, "H3 predecessor")
    exact_energy = result["exact_reference"]["energy_per_site"]
    points = result["points"]
    [point["n_super"] for point in points] == expected_depths || error("H3 predecessor has unexpected hierarchy depths")
    all_optimal = all(point["solve"]["status"] == "OPTIMAL" for point in points)
    all_below = all(validate_h3_point_bound(point, exact_energy, system_size) for point in points)
    energies = [point["energy_per_site"] for point in points]
    previous_energy = if system_size == H3_N
        result["comparison"]["D4_n6_validated_energy_per_site"]
    elseif system_size == H3_N50
        result["comparison"]["D4_n10_N20_provenance"]["energy_per_site"]
    elseif system_size == H3_N100
        result["comparison"]["D4_n25_N50_provenance"]["energy_per_site"]
    else
        error("unsupported H3 predecessor system size")
    end
    monotone = all(diff(vcat(previous_energy, energies)) .>= -2e-7)
    accuracy_passed = if system_size == H3_N
        last(energies) > result["comparison"]["D3_endpoint_energy_per_site"] + 2e-7
    elseif system_size == H3_N50
        only(point["energy_per_site"] for point in points if point["n_super"] == 20) > result["comparison"]["D3_n20"]["energy_per_site"] + 2e-7
    else
        last(energies) > result["comparison"]["D3_n50"]["energy_per_site"] + 2e-7
    end
    gates = result["gates"]
    all_optimal && get(gates, "all_optimal", false) || error("H3 OPTIMAL gate is inconsistent with its points")
    all_below && get(gates, "all_below_exact", false) || error("H3 sandwich gate is inconsistent with its points")
    monotone_gate = system_size == H3_N ? "fixed_D_monotone_from_n6" : system_size == H3_N50 ? "fixed_D_monotone_from_n10" : "fixed_D_monotone_from_n25"
    monotone && get(gates, monotone_gate, false) || error("H3 monotonicity gate is inconsistent with its points")
    all(point["energy_per_site"] < exact_energy for point in points) && get(gates, "all_raw_solver_objectives_below_exact", false) || error("H3 raw-objective gate is inconsistent with its points")
    all(point["residual_adjusted_lower_bound"] < exact_energy for point in points) && get(gates, "all_residual_adjusted_bounds_below_exact", false) || error("H3 residual-adjusted gate is inconsistent with its points")
    gate_key = system_size == H3_N ? "D4_endpoint_beats_D3" : system_size == H3_N50 ? "N50_advantage_passed" : "N100_advantage_passed"
    accuracy_passed && get(gates, gate_key, false) || error("H3 accuracy gate is inconsistent with its points")
    result
end

function plot_h3(results, figure_path)
    d3 = results["comparison"]["D3_ladder"]
    d4 = results["points"]
    exact_energy = results["exact_reference"]["energy_per_site"]
    d3_bounds = [get(point, "residual_adjusted_lower_bound", point["energy_per_site"]) for point in d3]
    d4_bounds = [point["residual_adjusted_lower_bound"] for point in d4]
    plot_result = plot(
        [point["support_sites"] for point in d3],
        exact_energy .- d3_bounds;
        marker=:circle,
        yscale=:log10,
        xlabel="constraint support (physical sites)",
        ylabel="finite-N gap: E_exact(20)/20 - E_lower",
        title="N=20 periodic Heisenberg chain",
        label="D=3 dense",
        legend=:topright,
        size=(800, 500),
    )
    plot!(
        plot_result,
        [point["support_sites"] for point in d4],
        exact_energy .- d4_bounds;
        marker=:diamond,
        label="D=4 U(1)",
    )
    mkpath(dirname(figure_path))
    savefig(plot_result, figure_path)
end

function h3_main(args)
    run_name, h0_result, h1_result = parse_h3_args(args)
    result_dir, result_path, figure_path = h3_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    h0 = validate_h0_predecessor(JSON.parsefile(h0_result))
    h1_n6 = validate_h1_predecessor(JSON.parsefile(h1_result))
    exact_reference = h0["exact_reference"]
    exact_energy = exact_reference["energy_per_site"]
    d3_ladder = h0["compressed"]
    d3_endpoint = last(d3_ladder)
    previous_energy = h1_n6["solves"]["u1_blocks"]["energy"]

    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
    fingerprint == h1_n6["meta"]["coarse_grainer_fingerprint"] || error("D=4 coarse-grainer differs from the validated H1 point")

    results = Dict(
        "meta" => Dict(
            "stage" => "H3",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H3_N,
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
                "h1_result" => relpath(h1_result, H3_REPO_ROOT),
            ),
        ),
        "exact_reference" => exact_reference,
        "coarse_grainer" => Dict("vumps" => vumps, "largest_forbidden_u1_entry" => largest_forbidden, "map_violations" => map_violations),
        "comparison" => Dict(
            "D3_ladder" => d3_ladder,
            "D3_endpoint_energy_per_site" => d3_endpoint["energy_per_site"],
            "D4_n6_validated_energy_per_site" => previous_energy,
        ),
        "points" => Dict{String, Any}[],
    )
    write_h1_json(result_path, results)

    for n_super in H3_DEPTHS
        support_sites = 2 * n_super
        support_sites <= H3_N || error("support exceeds the finite system size")
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
            "system_size" => H3_N,
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
        @printf("D=4 n=%d raw E/N %.9f, published lower bound %.6f, published gap %.3e\n", n_super, energy, published, exact_energy - published)
    end

    d4_endpoint = last(results["points"])
    results["gates"] = Dict(
        "all_optimal" => all(point["solve"]["status"] == "OPTIMAL" for point in results["points"]),
        "all_raw_solver_objectives_below_exact" => all(point["energy_per_site"] <= exact_energy for point in results["points"]),
        "all_raw_objectives_within_solver_tolerance" => all(point["energy_per_site"] <= exact_energy + H1_FEAS_TOL for point in results["points"]),
        "all_residual_adjusted_bounds_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "all_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in results["points"]),
        "fixed_D_monotone_from_n6" => true,
        "endpoint_improvement_over_D3" => d4_endpoint["energy_per_site"] - d3_endpoint["energy_per_site"],
        "D4_endpoint_beats_D3" => d4_endpoint["energy_per_site"] > d3_endpoint["energy_per_site"] + 2e-7,
    )
    write_h1_json(result_path, results)
    plot_h3(results, figure_path)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h3_main(ARGS)
end
