# Local exact-compatible RDM8 certificate with a dyadic D=4 coarse-grainer.

include(joinpath(@__DIR__, "rdm8_common.jl"))

function parse_rg_rdm8_local_args(args)
    run_name = nothing
    mosek_threads = RDM8_DEFAULT_MOSEK_THREADS
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-rg-rdm8-local")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name, mosek_threads
end

function rg_rdm8_local_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir,
    joinpath(result_dir, "dyadic_d4_n6.json"),
    joinpath(result_dir, "dyadic_d4.jls"),
    joinpath(result_dir, "figs", "dyadic_d4_n6.png")
end

function rg_rdm8_local_main(args)
    run_name, mosek_threads = parse_rg_rdm8_local_args(args)
    result_dir, result_path, tensor_path, figure_path = rg_rdm8_local_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    B_reference, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    B = dyadic_coarse_grainer(B_reference, RDM8_DYADIC_BITS)
    compatibility_audit = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    rdm8_map = validate_rdm8_map(B, virtual_charges, super_charges)

    certificate_dir = joinpath(result_dir, "certificates")
    GC.gc()
    baseline = solve_formulation(
        B,
        RDM8_N_SUPER,
        virtual_charges,
        super_charges;
        blocked=true,
        mosek_threads,
        strict_certificate=true,
        dual_transcript_path=joinpath(certificate_dir, "baseline.jls"),
        dyadic_bits=RDM8_DYADIC_BITS,
    )
    GC.gc()
    hybrid = solve_formulation(
        B,
        RDM8_N_SUPER,
        virtual_charges,
        super_charges;
        blocked=true,
        rdm4=true,
        mosek_threads,
        strict_certificate=true,
        dual_transcript_path=joinpath(certificate_dir, "rdm8.jls"),
        dyadic_bits=RDM8_DYADIC_BITS,
    )
    validated_solver_bounds(baseline, "H7 compressed baseline")
    validated_solver_bounds(hybrid, "H7 RDM8 hybrid")
    baseline_certificate = validated_dual_certificate(baseline, "H7 compressed baseline")
    hybrid_certificate = validated_dual_certificate(hybrid, "H7 RDM8 hybrid")

    raw_improvement = hybrid["energy"] - baseline["energy"]
    combined_residual = solver_residual_scale(baseline) + solver_residual_scale(hybrid)
    exact_difference = parse_exact_rational(hybrid_certificate["certified_lower_bound_rational"]) -
        parse_exact_rational(baseline_certificate["certified_lower_bound_rational"])
    physical_compatibility = compatibility_audit["physical_relaxation_compatibility_proven"] &&
        baseline_certificate["physical_relaxation_compatibility_proven"] &&
        hybrid_certificate["physical_relaxation_compatibility_proven"]
    both_below_bethe = baseline_certificate["certified_lower_bound"] < E0 &&
        hybrid_certificate["certified_lower_bound"] < E0
    residual_resolved = raw_improvement > max(RDM8_MIN_IMPROVEMENT, combined_residual)
    stage_passed = physical_compatibility && both_below_bethe && residual_resolved && exact_difference > 0

    results = Dict(
        "meta" => Dict(
            "stage" => "H7",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 2 * RDM8_N_SUPER,
            "n_super" => RDM8_N_SUPER,
            "bond_dimension" => H1_D,
            "dyadic_bits" => RDM8_DYADIC_BITS,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => tensor_fingerprint(B),
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
        ),
        "coarse_grainer" => Dict(
            "vumps" => vumps,
            "largest_forbidden_u1_entry" => largest_forbidden,
            "relative_tensor_distance" => norm(B - B_reference) / norm(B_reference),
            "map_violations" => map_violations,
            "rdm8_compression_violation" => rdm8_map["compression_violation"],
            "marginal_identity_residuals" => rdm8_map["marginal_identity_residuals"],
            "compatibility_audit" => compatibility_audit,
        ),
        "comparison" => Dict(
            "baseline" => baseline,
            "hybrid" => hybrid,
            "raw_improvement" => raw_improvement,
            "certified_lower_bound_difference" => Float64(exact_difference),
            "certified_lower_bound_difference_rational" => string(exact_difference),
            "combined_solver_residual_scale" => combined_residual,
        ),
        "gates" => Dict(
            "minimum_improvement" => RDM8_MIN_IMPROVEMENT,
            "both_below_bethe" => both_below_bethe,
            "sum_residual_resolved_improvement_passed" => residual_resolved,
            "physical_relaxation_compatibility_proven" => physical_compatibility,
            "stage_passed" => stage_passed,
        ),
    )
    write_h1_json(result_path, results)
    atomic_serialize(tensor_path, (B, Dict(
        "stage" => "H7",
        "run_name" => run_name,
        "dyadic_bits" => RDM8_DYADIC_BITS,
        "coarse_grainer_fingerprint" => tensor_fingerprint(B),
        "physical_relaxation_compatibility_proven" => physical_compatibility,
    )))
    plot_rdm8_comparison(results, figure_path)

    @printf("dyadic distance %.3e\n", results["coarse_grainer"]["relative_tensor_distance"])
    @printf("baseline %.9f, RDM8 %.9f\n", baseline["energy"], hybrid["energy"])
    @printf("raw improvement %.3e, exact replay-endpoint difference %.3e\n", raw_improvement, Float64(exact_difference))
    println("physical compatibility=$physical_compatibility, passed=$stage_passed")
    println("saved $result_path")
    println("saved $tensor_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    rg_rdm8_local_main(ARGS)
end
