# Exploratory fixed-parameter structured-NPA + exact-dyadic RG size screen.

const SIZE_SCREEN_PINNED_JULIA = v"1.11.5"

function parse_structured_rg_size_args(args)
    system_size = nothing
    n_super = nothing
    run_name = nothing
    feasibility_tolerance = 1e-7
    mosek_threads = 16
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
    (
        system_size=system_size,
        n_super=n_super,
        run_name=run_name,
        feasibility_tolerance=feasibility_tolerance,
        mosek_threads=mosek_threads,
    )
end

const SIZE_SCREEN_CONFIG = parse_structured_rg_size_args(ARGS)

module StructuredRGSizeScreenRuntime

const CONFIG = Main.SIZE_SCREEN_CONFIG
const SYSTEM_SIZE_NEEDLE = "const H10_N = 20"

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
    end
    Base.include(@__MODULE__, path)
end

include(joinpath(@__DIR__, "structured_rg_local.jl"))

function replay_summary(evaluated)
    Dict(
        "status" => "DIAGNOSTIC_EXACT_REPLAY",
        "raw_objective" => float_down(evaluated.objective),
        "strict_lower_endpoint" => float_down(evaluated.certified),
        "strict_lower_endpoint_rational" =>
            string(evaluated.certified),
        "identity_residual_rational" =>
            string(evaluated.residual_bound.identity_residual),
        "hamiltonian_residual_rational" =>
            string(evaluated.residual_bound.hamiltonian_residual),
        "other_residual_l1_rational" =>
            string(evaluated.residual_bound.other_residual_l1),
        "denominator_rational" =>
            string(evaluated.residual_bound.denominator),
        "legacy_lower_endpoint_rational" =>
            string(evaluated.residual_bound.legacy_lower_bound),
        "energy_aware_lower_endpoint_rational" =>
            string(evaluated.residual_bound.energy_aware_lower_bound),
    )
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
    legacy,
    objective,
    _,
    residuals,
    _ = h11_evaluate_transcript(transcript, structure)
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
        residual_bound=residual_bound,
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
    )
    hybrid_dual = dualize(primal_model)
    configure_structured_rg_optimizer!(
        structured.sohs_model,
        CONFIG.mosek_threads,
        CONFIG.feasibility_tolerance,
    )
    configure_structured_rg_optimizer!(
        hybrid_dual,
        CONFIG.mosek_threads,
        CONFIG.feasibility_tolerance,
    )
    build_wall = elapsed_seconds(build_start)

    baseline_wall = @elapsed optimize!(structured.sohs_model)
    require_optimal(
        structured.sohs_model,
        "size-screen structured baseline",
    )
    baseline = structured_rg_solve_record(
        structured.sohs_model,
        "structured_npa",
        baseline_wall,
    )
    hybrid_wall = @elapsed optimize!(hybrid_dual)
    require_optimal(
        hybrid_dual,
        "size-screen structured-NPA + RG",
    )
    hybrid = structured_rg_solve_record(
        hybrid_dual,
        "structured_npa_rg",
        hybrid_wall,
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
    )
    hybrid_replay_wall = @elapsed hybrid_replay =
        evaluate_structured_rg_transcript(
            hybrid_transcript,
            hybrid_structure,
            trace_bounds,
        )

    raw_improvement = hybrid["energy"] - baseline["energy"]
    strict_improvement =
        hybrid_replay.certified - baseline_replay.certified
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
            "solver_feasibility_tolerance" =>
                CONFIG.feasibility_tolerance,
            "formal_certificate" => false,
            "scope" =>
                "exploratory raw solve and same-process exact replay",
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
        ),
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
        "baseline" => baseline,
        "hybrid" => hybrid,
        "baseline_replay" => replay_summary(baseline_replay),
        "hybrid_replay" => replay_summary(hybrid_replay),
        "comparison" => Dict(
            "raw_improvement" => raw_improvement,
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
            "raw_monotonicity_passed" =>
                hybrid["energy"] +
                STRUCTURED_RG_ENERGY_TOLERANCE >=
                baseline["energy"],
            "strict_improvement_passed" =>
                strict_improvement > 0,
            "local_rdm_exactness_passed" =>
                local_rdm_exactness.status == "PROVEN" &&
                local_rdm_exactness.coefficient_assembly_exact,
            "physical_relaxation_compatibility_proven" =>
                compatibility[
                    "physical_relaxation_compatibility_proven"
                ],
        ),
    )
    write_h1_json(result_path, record)
    @printf(
        "N=%d n_super=%d raw improvement %.3e strict improvement %.3e\n",
        CONFIG.system_size,
        CONFIG.n_super,
        raw_improvement,
        Float64(strict_improvement),
    )
    result_path
end

end # module StructuredRGSizeScreenRuntime

if abspath(PROGRAM_FILE) == @__FILE__
    StructuredRGSizeScreenRuntime.run_size_screen()
end
