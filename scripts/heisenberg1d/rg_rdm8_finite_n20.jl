# Finite-periodic N=20 certificate for the exact-compatible RG/RDM8 relaxation.

include(joinpath(@__DIR__, "rdm8_common.jl"))

const RG_FINITE_N = 20
const RG_FINITE_N_SUPER = 10
const RG_FINITE_FEAS_TOL = 1e-7

include(joinpath(@__DIR__, "predecessor_validation.jl"))
include(joinpath(@__DIR__, "finite_upper_certificate.jl"))

function parse_rg_rdm8_finite_args(args)
    run_name = nothing
    mosek_threads = RDM8_DEFAULT_MOSEK_THREADS
    h3_result = nothing
    feasibility_tolerance = RG_FINITE_FEAS_TOL
    silent = true
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--h3-result=")
            spec = split(arg, "=", limit=2)[2]
            h3_result = normpath(isabspath(spec) ? spec : joinpath(@__DIR__, "..", "..", spec))
        elseif startswith(arg, "--feasibility-tolerance=")
            feasibility_tolerance = parse(Float64, split(arg, "=", limit=2)[2])
        elseif arg == "--debug"
            silent = false
        else
            error("unknown argument: $arg")
        end
    end
    mosek_threads >= 1 || error("mosek_threads must be positive")
    0 < feasibility_tolerance <= H1_FEAS_TOL || error("feasibility_tolerance must lie in (0, $H1_FEAS_TOL]")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-rg-rdm8-finite-n20")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(h3_result) && error("--h3-result is required")
    isfile(h3_result) || error("H3 N=20 result does not exist: $h3_result")
    run_name, mosek_threads, h3_result, feasibility_tolerance, silent
end

function rg_rdm8_finite_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "finite_n20_rdm8.json"), joinpath(result_dir, "figs", "finite_n20_rdm8.png")
end

function require_h3_metadata(result, key, expected)
    actual = get(get(result, "meta", Dict{String, Any}()), key, nothing)
    actual == expected || error("H3 predecessor has meta.$key=$actual, expected $expected")
end

function load_h3_recorded_input(result, key)
    inputs = get(get(result, "meta", Dict{String, Any}()), "inputs", Dict{String, Any}())
    haskey(inputs, key) || error("H3 predecessor does not record $key")
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    path = normpath(joinpath(repo_root, inputs[key]))
    isfile(path) || error("recorded H3 input does not exist: $path")
    JSON.parsefile(path)
end

validate_h3_n20_endpoint(path::AbstractString) = validate_h3_n20_endpoint(JSON.parsefile(path))

function validate_h3_n20_endpoint(result::AbstractDict)
    require_h3_metadata(result, "stage", "H3")
    require_h3_metadata(result, "system_scope", "finite_periodic_chain")
    require_h3_metadata(result, "system_size", RG_FINITE_N)
    require_h3_metadata(result, "boundary", "periodic")
    require_h3_metadata(result, "bond_dimension", H1_D)
    require_h3_metadata(result, "numerical_bound_schema_version", NUMERICAL_BOUND_SCHEMA_VERSION)

    exact_energy = result["exact_reference"]["energy_per_site"]
    fresh_reference = bethe_ground_reference(RG_FINITE_N)
    isapprox(exact_energy, fresh_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("H3 exact reference does not match an independent Bethe calculation")

    points = result["points"]
    [point["n_super"] for point in points] == [8, 10] || error("H3 predecessor has unexpected hierarchy depths")
    for point in points
        solve = point["solve"]
        adjusted, published = validated_solver_bounds(solve, "H3 N=20 D=4 point")
        point["energy_per_site"] == solve["energy"] || error("H3 point has a stale energy snapshot")
        point["residual_adjusted_lower_bound"] == adjusted || error("H3 point has a stale residual-adjusted snapshot")
        point["published_lower_bound"] == published || error("H3 point has a stale published snapshot")
        point["system_scope"] == "finite_periodic_chain" || error("H3 point has the wrong system scope")
        point["system_size"] == RG_FINITE_N || error("H3 point has the wrong system size")
        point["boundary"] == "periodic" || error("H3 point has the wrong boundary")
        point["bond_dimension"] == H1_D || error("H3 point has the wrong bond dimension")
        point["support_sites"] == 2 * point["n_super"] <= RG_FINITE_N || error("H3 point has invalid support")
    end
    gates = result["gates"]
    all_optimal = all(point["solve"]["status"] == "OPTIMAL" for point in points)
    all_below = all(point["energy_per_site"] < exact_energy for point in points)
    monotone = all(diff([result["comparison"]["D4_n6_validated_energy_per_site"]; [point["energy_per_site"] for point in points]]) .>= -2e-7)
    all_optimal && get(gates, "all_optimal", false) || error("H3 OPTIMAL gate is inconsistent")
    all_below && get(gates, "all_raw_solver_objectives_below_exact", false) || error("H3 exact-reference gate is inconsistent")
    monotone && get(gates, "fixed_D_monotone_from_n6", false) || error("H3 monotonicity gate is inconsistent")
    endpoint_beats_d3 = last(points)["energy_per_site"] > result["comparison"]["D3_endpoint_energy_per_site"] + 2e-7
    endpoint_beats_d3 && get(gates, "D4_endpoint_beats_D3", false) || error("H3 accuracy gate is inconsistent")

    h0 = load_h3_recorded_input(result, "h0_result")
    h1 = load_h3_recorded_input(result, "h1_result")
    validate_h0_predecessor(h0; system_size=RG_FINITE_N)
    validate_h1_predecessor(h1; bond_dimension=H1_D, n_super=6)
    result["exact_reference"] == h0["exact_reference"] || error("H3 exact-reference snapshot does not match H0")
    result["comparison"]["D3_ladder"] == h0["compressed"] || error("H3 D=3 ladder snapshot does not match H0")
    result["comparison"]["D3_endpoint_energy_per_site"] == last(h0["compressed"])["energy_per_site"] || error("H3 D=3 endpoint snapshot is stale")
    result["meta"]["coarse_grainer_fingerprint"] == h1["meta"]["coarse_grainer_fingerprint"] || error("H3 coarse-grainer fingerprint does not match H1")
    result["comparison"]["D4_n6_validated_energy_per_site"] == h1["solves"]["u1_blocks"]["energy"] || error("H3 D=4 n=6 snapshot does not match H1")

    endpoint = last(points)
    endpoint["n_super"] == RG_FINITE_N_SUPER || error("predecessor endpoint has the wrong hierarchy depth")
    endpoint["support_sites"] == RG_FINITE_N || error("predecessor endpoint has the wrong support")
    result
end

function plot_rg_rdm8_finite(results, figure_path)
    comparison = results["comparison"]
    labels = ["compressed", "+ RDM8"]
    gaps = [comparison["baseline_certified_gap"], comparison["hybrid_certified_gap"]]
    gap_plot = bar(
        labels,
        gaps;
        ylabel="finite-N certified gap",
        title="N=20, D=4, n=10",
        legend=false,
        color=[:gray60, :steelblue],
        ylims=(0, 1.12 * maximum(gaps)),
        left_margin=6mm,
        bottom_margin=5mm,
    )
    walls = [
        comparison["baseline"]["solver_wall_seconds"],
        comparison["hybrid"]["solver_wall_seconds"],
    ]
    wall_plot = bar(
        labels,
        walls;
        ylabel="Mosek wall time (s)",
        title="Finite hybrid validation cost",
        legend=false,
        color=[:gray60, :darkorange],
        ylims=(0, 1.08 * maximum(walls)),
        left_margin=5mm,
        bottom_margin=5mm,
    )
    mkpath(dirname(figure_path))
    savefig(plot(gap_plot, wall_plot; layout=(1, 2), size=(1050, 460)), figure_path)
end

function rg_rdm8_finite_main(args)
    run_name, mosek_threads, h3_result_path, feasibility_tolerance, silent = parse_rg_rdm8_finite_args(args)
    result_dir, result_path, figure_path = rg_rdm8_finite_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    h3 = validate_h3_n20_endpoint(h3_result_path)
    exact_reference = h3["exact_reference"]
    upper_certificate = build_finite_upper_certificate(
        RG_FINITE_N,
        joinpath(result_dir, "certificates", "finite-upper.jls"),
    )
    validated_finite_upper_certificate(upper_certificate, "N=20 finite reference")
    certified_upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
    predecessor_fingerprint = h3["meta"]["coarse_grainer_fingerprint"]

    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    B_reference, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    tensor_fingerprint(B_reference) == predecessor_fingerprint || error("source coarse-grainer differs from the H3 N=20 endpoint")
    B = dyadic_coarse_grainer(B_reference, RDM8_DYADIC_BITS)
    compatibility_audit = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = tensor_fingerprint(B)
    rdm8_map = validate_rdm8_map(B, virtual_charges, super_charges)

    GC.gc()
    baseline = solve_formulation(
        B,
        RG_FINITE_N_SUPER,
        virtual_charges,
        super_charges;
        blocked=true,
        mosek_threads,
        feasibility_tolerance,
        strict_certificate=true,
        dual_transcript_path=joinpath(result_dir, "certificates", "baseline.jls"),
        dyadic_bits=RDM8_DYADIC_BITS,
        silent,
    )
    GC.gc()
    hybrid = solve_formulation(
        B,
        RG_FINITE_N_SUPER,
        virtual_charges,
        super_charges;
        blocked=true,
        rdm4=true,
        mosek_threads,
        feasibility_tolerance,
        strict_certificate=true,
        dual_transcript_path=joinpath(result_dir, "certificates", "rdm8.jls"),
        dyadic_bits=RDM8_DYADIC_BITS,
        silent,
    )
    both_optimal = baseline["status"] == "OPTIMAL" && hybrid["status"] == "OPTIMAL"
    if !both_optimal
        failed_results = Dict(
            "meta" => Dict(
                "stage" => "H7",
                "run_name" => run_name,
                "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
                "system_scope" => "finite_periodic_chain",
                "system_size" => RG_FINITE_N,
                "boundary" => "periodic",
                "support_sites" => 2 * RG_FINITE_N_SUPER,
                "n_super" => RG_FINITE_N_SUPER,
                "bond_dimension" => H1_D,
                "strengthening" => "eight_physical_site_rdm",
                "physical_rdm_sites" => 8,
                "solver_feasibility_tolerance" => feasibility_tolerance,
                "mosek_threads" => mosek_threads,
                "coarse_grainer_fingerprint" => fingerprint,
                "dyadic_bits" => RDM8_DYADIC_BITS,
                "inputs" => Dict(
                    "h3_n20_result" => relpath(h3_result_path, normpath(joinpath(@__DIR__, "..", ".."))),
                ),
            ),
            "exact_reference" => exact_reference,
            "finite_energy_upper_certificate" => upper_certificate,
            "attempts" => Dict("baseline" => baseline, "hybrid" => hybrid),
            "gates" => Dict("both_optimal" => false, "stage_passed" => false),
        )
        write_h1_json(result_path, failed_results)
        println("baseline status $(baseline["status"]), hybrid status $(hybrid["status"])")
        println("saved failed attempt $result_path")
        error("both finite validation solves must terminate OPTIMAL")
    end
    validated_solver_bounds(baseline, "N=20 D=4 n=10 compressed baseline")
    validated_solver_bounds(hybrid, "N=20 D=4 n=10 RDM8 hybrid")
    baseline_certificate = validated_dual_certificate(baseline, "N=20 D=4 n=10 compressed baseline")
    hybrid_certificate = validated_dual_certificate(hybrid, "N=20 D=4 n=10 RDM8 hybrid")

    raw_improvement = hybrid["energy"] - baseline["energy"]
    largest_residual = max(solver_residual_scale(baseline), solver_residual_scale(hybrid))
    combined_residual = solver_residual_scale(baseline) + solver_residual_scale(hybrid)
    certified_difference_exact = parse_exact_rational(hybrid_certificate["certified_lower_bound_rational"]) -
        parse_exact_rational(baseline_certificate["certified_lower_bound_rational"])
    certified_difference = Float64(certified_difference_exact)
    baseline_certified = parse_exact_rational(baseline_certificate["certified_lower_bound_rational"])
    hybrid_certified = parse_exact_rational(hybrid_certificate["certified_lower_bound_rational"])
    below_certified_upper = baseline_certified < certified_upper && hybrid_certified < certified_upper
    residual_resolved = raw_improvement > max(RDM8_MIN_IMPROVEMENT, combined_residual)
    physical_compatibility = compatibility_audit["physical_relaxation_compatibility_proven"] &&
        baseline_certificate["physical_relaxation_compatibility_proven"] &&
        hybrid_certificate["physical_relaxation_compatibility_proven"]
    stage_passed = physical_compatibility && below_certified_upper && residual_resolved && certified_difference_exact > 0
    baseline_gap = certified_upper - parse_exact_rational(baseline_certificate["published_certified_lower_bound_rational"])
    hybrid_gap = certified_upper - parse_exact_rational(hybrid_certificate["published_certified_lower_bound_rational"])

    results = Dict(
        "meta" => Dict(
            "stage" => "H7",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => RG_FINITE_N,
            "boundary" => "periodic",
            "support_sites" => 2 * RG_FINITE_N_SUPER,
            "n_super" => RG_FINITE_N_SUPER,
            "bond_dimension" => H1_D,
            "strengthening" => "eight_physical_site_rdm",
            "physical_rdm_sites" => 8,
            "contract" => "docs/FINITE_SIZE_CONTRACT.md + docs/HYBRID_RDM_CONTRACT.md + docs/DUAL_POSTCERTIFICATION_CONTRACT.md",
            "solver_feasibility_tolerance" => feasibility_tolerance,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => fingerprint,
            "source_coarse_grainer_fingerprint" => predecessor_fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
            "dyadic_bits" => RDM8_DYADIC_BITS,
            "inputs" => Dict(
                "h3_n20_result" => relpath(h3_result_path, normpath(joinpath(@__DIR__, "..", ".."))),
            ),
        ),
        "exact_reference" => exact_reference,
        "finite_energy_upper_certificate" => upper_certificate,
        "coarse_grainer" => Dict(
            "vumps" => vumps,
            "largest_forbidden_u1_entry" => largest_forbidden,
            "relative_tensor_distance" => norm(B - B_reference) / norm(B_reference),
            "compatibility_audit" => compatibility_audit,
            "map_violations" => map_violations,
            "rdm8_compression_violation" => rdm8_map["compression_violation"],
            "marginal_identity_residuals" => rdm8_map["marginal_identity_residuals"],
        ),
        "comparison" => Dict(
            "baseline" => baseline,
            "hybrid" => hybrid,
            "raw_improvement" => raw_improvement,
            "certified_lower_bound_difference" => certified_difference,
            "certified_lower_bound_difference_rational" => string(certified_difference_exact),
            "largest_solver_residual_scale" => largest_residual,
            "combined_solver_residual_scale" => combined_residual,
            "baseline_certified_gap" => float_up(baseline_gap),
            "baseline_certified_gap_rational" => string(baseline_gap),
            "hybrid_certified_gap" => float_up(hybrid_gap),
            "hybrid_certified_gap_rational" => string(hybrid_gap),
        ),
        "gates" => Dict(
            "minimum_improvement" => RDM8_MIN_IMPROVEMENT,
            "both_below_certified_upper" => below_certified_upper,
            "sum_residual_resolved_improvement_passed" => residual_resolved,
            "physical_relaxation_compatibility_proven" => physical_compatibility,
            "stage_passed" => stage_passed,
        ),
    )
    write_h1_json(result_path, results)
    plot_rg_rdm8_finite(results, figure_path)

    @printf("baseline %.9f, RDM8 %.9f\n", baseline["energy"], hybrid["energy"])
    @printf("raw improvement %.3e, certified-bound difference %.3e\n", raw_improvement, certified_difference)
    @printf("combined residual %.3e, finite certified gap %.3e, passed=%s\n", combined_residual, results["comparison"]["hybrid_certified_gap"], stage_passed)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    rg_rdm8_finite_main(ARGS)
end
