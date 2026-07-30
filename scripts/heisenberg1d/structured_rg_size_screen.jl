# Exploratory fixed-parameter structured-NPA + exact-dyadic RG size screen.

const SIZE_SCREEN_PINNED_JULIA = v"1.11.5"

function parse_structured_rg_size_args(args)
    system_size = nothing
    n_super = nothing
    run_name = nothing
    feasibility_tolerance = 1e-7
    mosek_threads = 16
    coefficient_scale = 16384.0
    solver_log = false
    omega_normalization = :none
    for arg in args
        if startswith(arg, "--system-size=")
            system_size = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--n-super=")
            n_super = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--feasibility-tolerance=")
            feasibility_tolerance =
                parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--coefficient-scale=")
            coefficient_scale =
                parse(Float64, split(arg, "=", limit=2)[2])
        elseif arg == "--solver-log"
            solver_log = true
        elseif startswith(arg, "--omega-normalization=")
            omega_normalization =
                Symbol(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    VERSION == SIZE_SCREEN_PINNED_JULIA ||
        error("this screen requires Julia $SIZE_SCREEN_PINNED_JULIA")
    isnothing(system_size) && error("--system-size is required")
    isnothing(n_super) && error("--n-super is required")
    isnothing(run_name) && error("--run-name is required")
    iseven(system_size) || error("system_size must be even")
    system_size >= 20 || error("system_size must be at least 20")
    n_super >= 5 || error("n_super must be at least 5")
    2 * n_super <= system_size ||
        error("the RG support exceeds the requested finite system")
    !isempty(run_name) && basename(run_name) == run_name ||
        error("run_name must be one directory name")
    0 < feasibility_tolerance <= 1e-7 ||
        error("feasibility_tolerance must lie in (0, 1e-7]")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isfinite(coefficient_scale) && coefficient_scale >= 1 ||
        error("coefficient_scale must be finite and at least one")
    isinteger(coefficient_scale) &&
        coefficient_scale <= typemax(Int) &&
        ispow2(Int(coefficient_scale)) ||
        error("coefficient_scale must be an exact power-of-two integer")
    omega_normalization in (:none, :dyadic_trace) ||
        error("omega_normalization must be none or dyadic_trace")
    (
        system_size=system_size,
        n_super=n_super,
        run_name=run_name,
        feasibility_tolerance=feasibility_tolerance,
        mosek_threads=mosek_threads,
        coefficient_scale=coefficient_scale,
        solver_log=solver_log,
        omega_normalization=omega_normalization,
    )
end

const SIZE_SCREEN_CONFIG = parse_structured_rg_size_args(ARGS)

module StructuredRGSizeScreenRuntime

const CONFIG = Main.SIZE_SCREEN_CONFIG
const SYSTEM_SIZE_NEEDLE = "const H10_N = 20"
const COEFFICIENT_SCALE_NEEDLE = "const H11_EQUALITY_SCALE = 16384.0"

function include(path::AbstractString)
    if basename(path) == "h10_qmbcertify_n20.jl"
        source = read(path, String)
        length(findall(SYSTEM_SIZE_NEEDLE, source)) == 1 ||
            error("the pinned N=20 source no longer has the expected size constant")
        instrumented = replace(
            source,
            SYSTEM_SIZE_NEEDLE =>
                "const H10_N = $(CONFIG.system_size)";
            count=1,
        )
        return Base.include_string(@__MODULE__, instrumented, path)
    elseif basename(path) == "h11_full_sohs_certificate.jl"
        source = read(path, String)
        length(findall(COEFFICIENT_SCALE_NEEDLE, source)) == 1 ||
            error("the pinned H11 source no longer has the expected coefficient scale")
        instrumented = replace(
            source,
            COEFFICIENT_SCALE_NEEDLE =>
                "const H11_EQUALITY_SCALE = $(CONFIG.coefficient_scale)";
            count=1,
        )
        return Base.include_string(@__MODULE__, instrumented, path)
    end
    Base.include(@__MODULE__, path)
end

include(joinpath(@__DIR__, "structured_rg_local.jl"))

function size_screen_failure_meta()
    Dict(
        "stage" => "structured_rg_size_screen_failure",
        "run_name" => CONFIG.run_name,
        "julia_version" => string(VERSION),
        "system_size" => CONFIG.system_size,
        "n_super" => CONFIG.n_super,
        "support_sites" => 2 * CONFIG.n_super,
        "formal_certificate" => false,
        "mosek_threads" => CONFIG.mosek_threads,
        "coefficient_scale" => CONFIG.coefficient_scale,
        "solver_log" => CONFIG.solver_log,
        "solver_option_scope" => "structured_npa_and_fusion",
        "omega_normalization" =>
            string(CONFIG.omega_normalization),
        "solver_feasibility_tolerance" =>
            CONFIG.feasibility_tolerance,
    )
end

function replay_residual_diagnostics(residuals, canonical_words)
    length(residuals) == length(canonical_words) ||
        error("replay residuals do not match the canonical words")
    magnitudes = abs.(Float64.(residuals))
    top_indices = partialsortperm(
        magnitudes,
        1:min(10, length(magnitudes));
        rev=true,
    )
    Dict(
        "count" => length(residuals),
        "nonzero_count" => count(!iszero, residuals),
        "absolute_l1" => sum(magnitudes),
        "maximum_absolute_residual" => maximum(magnitudes),
        "absolute_residual_gt_1e_minus_9_count" =>
            count(>(1e-9), magnitudes),
        "absolute_residual_gt_1e_minus_8_count" =>
            count(>(1e-8), magnitudes),
        "absolute_residual_gt_1e_minus_7_count" =>
            count(>(1e-7), magnitudes),
        "largest" => [
            Dict(
                "index" => index,
                "canonical_word" => canonical_words[index],
                "residual" => Float64(residuals[index]),
            )
            for index in top_indices
        ],
    )
end

function replay_psd_repair_diagnostics(groups)
    shifted = [
        let
            minimum_eigenvalue = haskey(
                values,
                "minimum_eigenvalue_lower_bound",
            ) ?
                values["minimum_eigenvalue_lower_bound"] :
                float_down(parse_exact_rational(
                    values[
                        "minimum_eigenvalue_lower_bound_rational"
                    ],
                ))
            diagonal_shift = haskey(values, "diagonal_shift") ?
                values["diagonal_shift"] :
                float_up(parse_exact_rational(
                    values["diagonal_shift_rational"],
                ))
            Dict(
                "group" => group,
                "dimension" => values["dimension"],
                "minimum_eigenvalue_lower_bound" =>
                    minimum_eigenvalue,
                "diagonal_shift" => diagonal_shift,
            )
        end
        for (group, values) in groups
        if parse_exact_rational(
            values["diagonal_shift_rational"],
        ) > 0
    ]
    sort!(
        shifted;
        by=entry -> entry["diagonal_shift"],
        rev=true,
    )
    Dict(
        "block_count" => length(groups),
        "shifted_block_count" => length(shifted),
        "maximum_diagonal_shift" =>
            isempty(shifted) ? 0.0 : shifted[1]["diagonal_shift"],
        "largest_shifts" => shifted[1:min(10, length(shifted))],
    )
end

function replay_summary(evaluated, canonical_words)
    summary = Dict(
        "status" => "DIAGNOSTIC_EXACT_REPLAY",
        "raw_objective" => float_down(evaluated.objective),
        "strict_lower_endpoint" => float_down(evaluated.certified),
        "total_replay_correction" =>
            float_down(evaluated.certified - evaluated.objective),
        "strict_lower_endpoint_rational" =>
            string(evaluated.certified),
        "identity_residual" =>
            Float64(evaluated.residual_bound.identity_residual),
        "identity_residual_rational" =>
            string(evaluated.residual_bound.identity_residual),
        "hamiltonian_residual" =>
            Float64(evaluated.residual_bound.hamiltonian_residual),
        "hamiltonian_residual_rational" =>
            string(evaluated.residual_bound.hamiltonian_residual),
        "other_residual_l1" =>
            float_up(evaluated.residual_bound.other_residual_l1),
        "other_residual_l1_rational" =>
            string(evaluated.residual_bound.other_residual_l1),
        "denominator_rational" =>
            string(evaluated.residual_bound.denominator),
        "legacy_lower_endpoint_rational" =>
            string(evaluated.residual_bound.legacy_lower_bound),
        "energy_aware_lower_endpoint_rational" =>
            string(evaluated.residual_bound.energy_aware_lower_bound),
        "raw_moment_residual_diagnostics" =>
            replay_residual_diagnostics(
                evaluated.raw_moment_residuals,
                canonical_words,
            ),
        "repaired_moment_residual_diagnostics" =>
            replay_residual_diagnostics(
                evaluated.moment_residuals,
                canonical_words,
            ),
        "psd_repair_diagnostics" =>
            replay_psd_repair_diagnostics(evaluated.sohs_groups),
    )
    if hasproperty(evaluated, :rg_groups)
        summary["rg_slack_correction"] =
            float_down(evaluated.rg_correction)
        summary["rg_groups"] = Dict(
            group => Dict(
                "block_count" => values["block_count"],
                "minimum_eigenvalue_lower_bound" =>
                    float_down(parse_exact_rational(
                        values[
                            "minimum_eigenvalue_lower_bound_rational"
                        ],
                    )),
                "trace_upper_bound" =>
                    float_up(parse_exact_rational(
                        values["trace_upper_bound_rational"],
                    )),
                "correction" =>
                    float_down(parse_exact_rational(
                        values["correction_rational"],
                    )),
            )
            for (group, values) in evaluated.rg_groups
        )
    end
    summary
end

function replay_baseline(model, canonical_words)
    structure = h11_model_structure(model, canonical_words)
    transcript = merge(
        structure,
        (
            raw_values=Float64[
                value(variable)
                for variable in all_variables(model)
            ],
        ),
    )
    raw_values = exact_rational.(transcript.raw_values)
    raw_residuals = ExactRational[
        h11_evaluate_linear(equality, raw_values) -
        exact_rational(equality.rhs)
        for equality in transcript.equalities
    ]
    legacy,
    objective,
    _,
    residuals,
    groups = h11_evaluate_transcript(transcript, structure)
    residual_bound = h11_energy_aware_residual_bound(
        objective,
        residuals,
        transcript.canonical_words,
    )
    residual_bound.legacy_lower_bound == legacy ||
        error("baseline replay is internally inconsistent")
    (
        certified=residual_bound.certified_lower_bound,
        objective=objective,
        raw_moment_residuals=raw_residuals,
        moment_residuals=residuals,
        residual_bound=residual_bound,
        sohs_groups=groups,
        structure=structure,
    )
end

function run_size_screen()
    H10_N == CONFIG.system_size ||
        error("the private size instrumentation did not take effect")
    validate_structured_rg_support(
        CONFIG.n_super;
        minimum_depth=5,
    )
    result_dir = joinpath(
        BOUNDGSENERGY_ROOT,
        "results",
        CONFIG.run_name,
    )
    result_path = joinpath(result_dir, "size_screen.json")
    ispath(result_dir) &&
        error("result directory already exists; choose a fresh --run-name")
    mkpath(result_dir)

    Braw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    Bcharged, forbidden = exact_charge_tensor(
        Braw,
        STRUCTURED_RG_VIRTUAL_CHARGES,
        STRUCTURED_RG_SUPER_CHARGES,
    )
    B = dyadic_coarse_grainer(Bcharged, RDM8_DYADIC_BITS)
    compatibility = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)

    rg_only = solve_formulation(
        B,
        CONFIG.n_super,
        STRUCTURED_RG_VIRTUAL_CHARGES,
        STRUCTURED_RG_SUPER_CHARGES;
        blocked=true,
        rdm4=true,
        mosek_threads=CONFIG.mosek_threads,
        feasibility_tolerance=CONFIG.feasibility_tolerance,
    )
    if rg_only["status"] != "OPTIMAL"
        write_h1_json(
            result_path,
            Dict(
                "meta" => size_screen_failure_meta(),
                "failure" => Dict(
                    "phase" => "rg_rdm8",
                    "classification" =>
                        "NON_OPTIMAL_SOLVER_STATUS",
                ),
                "rg_only" => rg_only,
                "gates" => Dict(
                    "all_three_optimal" => false,
                    "strict_improvement_passed" => false,
                ),
            ),
        )
        error("size-screen RG-only parent did not terminate OPTIMAL")
    end
    rg_only["label"] = "rg_rdm8"
    rg_only["system_scope"] = "finite_periodic_chain"
    rg_only["system_size"] = CONFIG.system_size
    rg_only["support_sites"] = 2 * CONFIG.n_super
    rg_only["bond_dimension"] = size(B, 1)
    rg_only["formal_certificate"] = false
    GC.gc()

    build_start = time_ns()
    structured = build_structured_moment_model()
    primal_model = structured.model
    structured_rdm = structured_local_rdm(
        structured.physical_moments,
        structured.canonical_words,
        8,
    )
    local_rdm_exactness = structured_local_rdm_exactness_audit(
        structured_rdm,
        structured.physical_moments,
        structured.canonical_words,
        8,
    )
    hierarchy = add_rdm8_rg_hierarchy!(
        primal_model,
        structured_rdm,
        B,
        CONFIG.n_super,
        STRUCTURED_RG_VIRTUAL_CHARGES,
        STRUCTURED_RG_SUPER_CHARGES,
        omega_normalization=CONFIG.omega_normalization,
    )
    hybrid_dual = dualize(primal_model)
    configure_structured_rg_optimizer!(
        structured.sohs_model,
        CONFIG.mosek_threads,
        CONFIG.feasibility_tolerance,
        solver_log=CONFIG.solver_log,
    )
    configure_structured_rg_optimizer!(
        hybrid_dual,
        CONFIG.mosek_threads,
        CONFIG.feasibility_tolerance,
        solver_log=CONFIG.solver_log,
    )
    build_wall = elapsed_seconds(build_start)

    baseline_wall = @elapsed optimize!(structured.sohs_model)
    baseline = structured_rg_solve_record(
        structured.sohs_model,
        "structured_npa",
        baseline_wall,
    )
    function persist_solver_failure(phase, hybrid=nothing)
        failure_record = Dict(
            "meta" => size_screen_failure_meta(),
            "failure" => Dict(
                "phase" => phase,
                "classification" => "NON_OPTIMAL_SOLVER_STATUS",
            ),
            "model" => Dict(
                "build_wall_seconds" => build_wall,
                "baseline_scalar_variables" =>
                    num_variables(structured.sohs_model),
                "hybrid_scalar_variables" =>
                    num_variables(hybrid_dual),
                "canonical_moment_count" =>
                    length(structured.canonical_words),
                "local_rdm_link_equality_count" =>
                    hierarchy.rho4_link_equality_count,
                "omega4_link_equality_count" =>
                    hierarchy.omega4_link_equality_count,
                "omega_count" => length(hierarchy.omega),
                "omega_trace_scales" =>
                    hierarchy.omega_trace_scales,
                "link_scales" => Dict(
                    "local" => hierarchy.local_link_scale,
                    "omega4" => hierarchy.omega4_link_scale,
                    "recursion" => hierarchy.recursion_scale,
                ),
            ),
            "rg_only" => rg_only,
            "baseline" => baseline,
            "hybrid" => hybrid,
            "gates" => Dict(
                "all_three_optimal" => false,
                "strict_improvement_passed" => false,
            ),
        )
        write_h1_json(result_path, failure_record)
    end
    if termination_status(structured.sohs_model) != MOI.OPTIMAL
        persist_solver_failure("structured_npa")
    end
    require_optimal(
        structured.sohs_model,
        "size-screen structured baseline",
    )
    hybrid_wall = @elapsed optimize!(hybrid_dual)
    hybrid = structured_rg_solve_record(
        hybrid_dual,
        "structured_npa_rg",
        hybrid_wall,
    )
    if termination_status(hybrid_dual) != MOI.OPTIMAL
        persist_solver_failure("structured_npa_rg", hybrid)
    end
    require_optimal(
        hybrid_dual,
        "size-screen structured-NPA + RG",
    )

    baseline_replay_wall = @elapsed baseline_replay =
        replay_baseline(
            structured.sohs_model,
            structured.canonical_words,
        )
    hybrid_structure = structured_rg_model_structure(
        primal_model,
        hybrid_dual,
        structured,
        structured_rdm,
        B,
        CONFIG.n_super,
        link_mode=:rdm8,
        omega_normalization=CONFIG.omega_normalization,
    )
    hybrid_transcript = merge(
        hybrid_structure,
        (
            raw_values=Float64[
                value(variable)
                for variable in all_variables(hybrid_dual)
            ],
        ),
    )
    trace_bounds = structured_rg_trace_bounds(
        B,
        CONFIG.n_super;
        link_mode=:rdm8,
        omega_normalization=CONFIG.omega_normalization,
    )
    hybrid_raw_values =
        exact_rational.(hybrid_transcript.raw_values)
    hybrid_raw_moment_residuals = ExactRational[
        structured_rg_evaluate_affine(
            equality,
            hybrid_raw_values,
        ) / exact_rational(hybrid_transcript.equality_scale)
        for equality in hybrid_transcript.moment_equalities
    ]
    hybrid_replay_wall = @elapsed hybrid_evaluated =
        evaluate_structured_rg_transcript(
            hybrid_transcript,
            hybrid_structure,
            trace_bounds,
        )
    hybrid_replay = merge(
        hybrid_evaluated,
        (raw_moment_residuals=hybrid_raw_moment_residuals,),
    )

    raw_improvement = hybrid["energy"] - baseline["energy"]
    strict_improvement =
        hybrid_replay.certified - baseline_replay.certified
    hybrid_over_rg_only =
        hybrid["energy"] - rg_only["energy"]
    structured_over_rg_only =
        baseline["energy"] - rg_only["energy"]
    exact_reference =
        bethe_ground_reference(CONFIG.system_size)
    exact_energy = exact_reference["energy_per_site"]
    all_dimensions = [
        block.dimension
        for block in vcat(
            baseline_replay.structure.psd_blocks,
            hybrid_structure.psd_variable_blocks,
            hybrid_structure.omega_slacks,
        )
    ]
    record = Dict(
        "meta" => Dict(
            "stage" => "structured_rg_size_screen",
            "run_name" => CONFIG.run_name,
            "julia_version" => string(VERSION),
            "system_size" => CONFIG.system_size,
            "boundary" => "periodic",
            "n_super" => CONFIG.n_super,
            "support_sites" => 2 * CONFIG.n_super,
            "bond_dimension" => size(B, 1),
            "dyadic_bits" => RDM8_DYADIC_BITS,
            "link_mode" => "rdm8",
            "mosek_threads" => CONFIG.mosek_threads,
            "coefficient_scale" => CONFIG.coefficient_scale,
            "solver_log" => CONFIG.solver_log,
            "solver_option_scope" =>
                "structured_npa_and_fusion",
            "omega_normalization" =>
                string(CONFIG.omega_normalization),
            "solver_feasibility_tolerance" =>
                CONFIG.feasibility_tolerance,
            "formal_certificate" => false,
            "scope" =>
                "matched three-route exploratory solve and same-process fusion replay",
            "system_size_instrumentation" =>
                "private-module replacement of the pinned H10_N constant",
            "equality_scales" => Dict(
                "local_link" => hybrid_structure.local_link_scale,
                "omega4_link" => hybrid_structure.omega4_link_scale,
                "recursion" => hybrid_structure.recursion_scale,
            ),
            "structured_configuration" => Dict(
                "order" => structured.configuration.order,
                "extra" => structured.configuration.extra,
                "rdm" => structured.configuration.rdm,
                "pso" => structured.configuration.pso,
                "lso" => structured.configuration.lso,
            ),
            "routes" => [
                "rg_rdm8",
                "structured_npa",
                "structured_npa_rg",
            ],
        ),
        "exact_reference" => exact_reference,
        "coarse_grainer" => Dict(
            "vumps" => vumps,
            "loaded_from_cache" => loaded_from_cache,
            "largest_forbidden_charge_entry" => forbidden,
            "compatibility_audit" => compatibility,
        ),
        "local_rdm_exactness_audit" => Dict(
            string(field) => getproperty(local_rdm_exactness, field)
            for field in propertynames(local_rdm_exactness)
        ),
        "model" => Dict(
            "build_wall_seconds" => build_wall,
            "baseline_scalar_variables" =>
                baseline_replay.structure.variable_count,
            "rg_only_scalar_variables" =>
                rg_only["scalar_variables"],
            "hybrid_scalar_variables" =>
                hybrid_structure.variable_count,
            "baseline_psd_block_count" =>
                length(baseline_replay.structure.psd_blocks),
            "hybrid_sohs_psd_block_count" =>
                length(hybrid_structure.psd_variable_blocks),
            "hybrid_rg_slack_block_count" =>
                length(hybrid_structure.omega_slacks),
            "largest_psd_block" => maximum(all_dimensions),
            "canonical_moment_count" =>
                length(structured.canonical_words),
            "local_rdm_link_equality_count" =>
                hierarchy.rho4_link_equality_count,
            "omega4_link_equality_count" =>
                hierarchy.omega4_link_equality_count,
            "omega_count" => length(hierarchy.omega),
            "omega_trace_scales" =>
                hierarchy.omega_trace_scales,
            "trace_bounds" => Dict(
                group => Dict(
                    "value" => float_up(bound),
                    "rational" => string(bound),
                )
                for (group, bound) in sort!(
                    collect(trace_bounds);
                    by=first,
                )
            ),
        ),
        "rg_only" => rg_only,
        "baseline" => baseline,
        "hybrid" => hybrid,
        "baseline_replay" => replay_summary(
            baseline_replay,
            structured.canonical_words,
        ),
        "hybrid_replay" => replay_summary(
            hybrid_replay,
            structured.canonical_words,
        ),
        "comparison" => Dict(
            "raw_improvement" => raw_improvement,
            "hybrid_over_rg_only" => hybrid_over_rg_only,
            "structured_over_rg_only" =>
                structured_over_rg_only,
            "strict_improvement" => float_down(strict_improvement),
            "strict_improvement_rational" =>
                string(strict_improvement),
            "baseline_replay_wall_seconds" =>
                baseline_replay_wall,
            "hybrid_replay_wall_seconds" =>
                hybrid_replay_wall,
        ),
        "gates" => Dict(
            "both_optimal" => true,
            "all_three_optimal" => true,
            "raw_monotonicity_passed" =>
                hybrid["energy"] +
                STRUCTURED_RG_ENERGY_TOLERANCE >=
                baseline["energy"],
            "fusion_dominates_rg_only" =>
                hybrid["energy"] +
                STRUCTURED_RG_ENERGY_TOLERANCE >=
                rg_only["energy"],
            "all_raw_objectives_below_finite_reference" =>
                all(
                    energy <=
                    exact_energy + STRUCTURED_RG_ENERGY_TOLERANCE
                    for energy in (
                        rg_only["energy"],
                        baseline["energy"],
                        hybrid["energy"],
                    )
                ),
            "strict_improvement_passed" =>
                strict_improvement > 0,
            "local_rdm_exactness_passed" =>
                local_rdm_exactness.status == "PROVEN" &&
                local_rdm_exactness.coefficient_assembly_exact,
            "physical_relaxation_compatibility_proven" =>
                compatibility[
                    "physical_relaxation_compatibility_proven"
                ],
            "omega_normalization_exact" =>
                hierarchy.omega_trace_scales ==
                structured_rg_omega_trace_scales(
                    B,
                    CONFIG.n_super,
                    CONFIG.omega_normalization,
                ),
        ),
    )
    write_h1_json(result_path, record)
    @printf(
        "N=%d n_super=%d RG %.9f structured %.9f fusion %.9f raw fusion gain %.3e strict fusion gain %.3e\n",
        CONFIG.system_size,
        CONFIG.n_super,
        rg_only["energy"],
        baseline["energy"],
        hybrid["energy"],
        raw_improvement,
        Float64(strict_improvement),
    )
    result_path
end

end # module StructuredRGSizeScreenRuntime

if abspath(PROGRAM_FILE) == @__FILE__
    StructuredRGSizeScreenRuntime.run_size_screen()
end
