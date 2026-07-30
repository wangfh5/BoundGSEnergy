# Local structured-NPA + exact-dyadic RG comparison at N=20.

include(joinpath(@__DIR__, "structured_rg_certificate.jl"))

const STRUCTURED_RG_DEFAULT_N_SUPER = 5

function parse_structured_rg_args(args)
    run_name = nothing
    n_super = STRUCTURED_RG_DEFAULT_N_SUPER
    link_mode = :rho3
    mosek_threads = 16
    feasibility_tolerance = 1e-7
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--n-super=")
            n_super = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--link-mode=")
            link_mode = Symbol(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--feasibility-tolerance=")
            feasibility_tolerance =
                parse(Float64, split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) &&
        (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-structured-rg-n$(n_super)")
    !isempty(run_name) && basename(run_name) == run_name ||
        error("run_name must be one directory name")
    validate_structured_rg_support(n_super; minimum_depth=5)
    link_mode in (:rho3, :rdm8) ||
        error("link_mode must be rho3 or rdm8")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    0 < feasibility_tolerance <= 1e-7 ||
        error("feasibility_tolerance must lie in (0, 1e-7]")
    run_name, n_super, link_mode, mosek_threads, feasibility_tolerance
end

function structured_rg_solve_record(model, label, wall_seconds)
    status = termination_status(model)
    record = Dict{String, Any}(
        "label" => label,
        "status" => string(status),
        "primal_status" => string(primal_status(model)),
        "dual_status" => string(dual_status(model)),
        "wall_seconds" => wall_seconds,
        "solver_wall_seconds" => solve_time(model),
        "scalar_variables" => num_variables(model),
    )
    if status == MOI.OPTIMAL
        record["energy"] = objective_value(model)
        record["objective_bound"] = objective_bound(model)
        record["residuals"] = mosek_residuals(model)
    end
    record
end

function structured_local_rdm_checks(
    local_rdm,
    primal_model,
    values,
    objective,
)
    all(
        owner_model(variable) === primal_model
        for entry in local_rdm
        for (_, variable) in linear_terms(entry)
    ) || error("structured local RDM contains variables outside the primal moment model")
    numeric_local_rdm = evaluate_affine.(local_rdm, Ref(values))
    numeric_rho3, translation_residual = if size(local_rdm) == (64, 64)
        rho2_left = ptrace_left(numeric_local_rdm, 4, 16)
        rho2_right = ptrace_right(numeric_local_rdm, 16, 4)
        numeric_local_rdm, maximum(abs, rho2_left - rho2_right)
    elseif size(local_rdm) == (256, 256)
        rho3_left = ptrace_left(numeric_local_rdm, 4, 64)
        rho3_right = ptrace_right(numeric_local_rdm, 64, 4)
        rho3_left, maximum(abs, rho3_left - rho3_right)
    else
        error("unsupported structured local-RDM dimension")
    end
    rho2_left = ptrace_left(numeric_rho3, 4, 16)
    Dict(
        "trace" => tr(numeric_local_rdm),
        "minimum_eigenvalue" =>
            eigmin(Symmetric(numeric_local_rdm)),
        "symmetry_residual" =>
            maximum(abs, numeric_local_rdm - numeric_local_rdm'),
        "translation_marginal_residual" => translation_residual,
        "local_energy" => tr(HLOC * rho2_left) / 2,
        "objective_energy_residual" =>
            abs(tr(HLOC * rho2_left) / 2 - objective),
    )
end

function configure_structured_rg_optimizer!(
    model,
    mosek_threads,
    feasibility_tolerance,
)
    set_optimizer(model, MosekTools.Optimizer)
    set_silent(model)
    set_attribute(model, "MSK_IPAR_NUM_THREADS", mosek_threads)
    set_attribute(
        model,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS",
        feasibility_tolerance,
    )
    set_attribute(
        model,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS",
        feasibility_tolerance,
    )
    set_attribute(
        model,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",
        feasibility_tolerance,
    )
    model
end

function require_optimal(model, label)
    status = termination_status(model)
    status == MOI.OPTIMAL && return
    diagnostics = [
        "termination_status=$(status)",
        "primal_status=$(primal_status(model))",
        "dual_status=$(dual_status(model))",
        "raw_status=$(raw_status(model))",
    ]
    has_values(model) && push!(diagnostics, "objective=$(objective_value(model))")
    error("$label did not terminate OPTIMAL: " * join(diagnostics, ", "))
end

function plot_structured_rg(record, figure_path)
    labels = ["structured NPA", "structured NPA + RG"]
    energies = [record["baseline"]["energy"], record["hybrid"]["energy"]]
    reference = record["finite_energy_upper_certificate"]["upper_bound"]
    plot = scatter(
        1:2,
        energies;
        xticks=(1:2, labels),
        xlims=(0.5, 2.5),
        ylabel="energy per physical site",
        title="N=20 structured NPA + RG",
        marker=:circle,
        markersize=7,
        legend=false,
        color=[:gray, :steelblue],
        yformatter=value -> @sprintf("%.8f", value),
        bottom_margin=8mm,
        right_margin=5mm,
    )
    hline!(plot, [reference]; color=:black, linestyle=:dash, label=false)
    mkpath(dirname(figure_path))
    savefig(plot, figure_path)
end

function structured_rg_main(args)
    run_name,
    n_super,
    link_mode,
    mosek_threads,
    feasibility_tolerance =
        parse_structured_rg_args(args)
    result_dir = joinpath(BOUNDGSENERGY_ROOT, "results", run_name)
    result_path = joinpath(result_dir, "structured_rg.json")
    figure_path = joinpath(result_dir, "figs", "structured_rg.png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")
    mkpath(result_dir)

    Braw, vumps, coarse_grainer_loaded = load_or_build_coarse_grainer()
    virtual_charges = STRUCTURED_RG_VIRTUAL_CHARGES
    super_charges = STRUCTURED_RG_SUPER_CHARGES
    Bcharged, forbidden = exact_charge_tensor(Braw, virtual_charges, super_charges)
    B = dyadic_coarse_grainer(Bcharged, RDM8_DYADIC_BITS)
    compatibility = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)

    build_start = time_ns()
    structured = build_structured_moment_model()
    primal_model = structured.model
    local_site_count = link_mode == :rho3 ? 6 : 8
    structured_rdm = structured_local_rdm(
        structured.physical_moments,
        structured.canonical_words,
        local_site_count,
    )
    local_rdm_exactness = structured_local_rdm_exactness_audit(
        structured_rdm,
        structured.physical_moments,
        structured.canonical_words,
        local_site_count,
    )
    hierarchy = if link_mode == :rho3
        add_rg_hierarchy!(
            primal_model,
            structured_rdm,
            B,
            n_super,
            virtual_charges,
            super_charges,
        )
    else
        add_rdm8_rg_hierarchy!(
            primal_model,
            structured_rdm,
            B,
            n_super,
            virtual_charges,
            super_charges,
        )
    end
    hybrid_dual = dualize(primal_model)
    configure_structured_rg_optimizer!(
        structured.sohs_model,
        mosek_threads,
        feasibility_tolerance,
    )
    configure_structured_rg_optimizer!(
        hybrid_dual,
        mosek_threads,
        feasibility_tolerance,
    )
    build_wall = elapsed_seconds(build_start)

    baseline_wall = @elapsed optimize!(structured.sohs_model)
    require_optimal(structured.sohs_model, "structured-NPA baseline")
    baseline =
        structured_rg_solve_record(structured.sohs_model, "structured_npa", baseline_wall)
    baseline_moments = structured_moment_values_from_sohs(structured)
    baseline_rdm_checks = structured_local_rdm_checks(
        structured_rdm,
        primal_model,
        baseline_moments,
        objective_value(structured.sohs_model),
    )
    @info "structured local-RDM baseline checks" baseline_rdm_checks
    abs(baseline_rdm_checks["trace"] - 1) <= 1e-7 ||
        error("structured local RDM is not normalized")
    baseline_rdm_checks["minimum_eigenvalue"] >= -1e-7 ||
        error("structured local RDM is not PSD")
    baseline_rdm_checks["translation_marginal_residual"] <= 1e-7 ||
        error("structured local-RDM marginals are not translation invariant")
    baseline_rdm_checks["objective_energy_residual"] <= 1e-7 ||
        error("structured and local-RDM energy conventions disagree")

    hybrid_wall = @elapsed optimize!(hybrid_dual)
    require_optimal(hybrid_dual, "structured-NPA + RG model")
    hybrid =
        structured_rg_solve_record(hybrid_dual, "structured_npa_rg", hybrid_wall)
    hybrid_moments = Dict(
        variable => dualized_primal_value(hybrid_dual, variable)
        for variable in keys(baseline_moments)
    )
    rdm_checks = structured_local_rdm_checks(
        structured_rdm,
        primal_model,
        hybrid_moments,
        objective_value(hybrid_dual),
    )
    abs(rdm_checks["trace"] - 1) <= 1e-7 ||
        error("structured local RDM is not normalized")
    rdm_checks["minimum_eigenvalue"] >= -1e-7 ||
        error("structured local RDM is not PSD")
    rdm_checks["translation_marginal_residual"] <= 1e-7 ||
        error("structured local-RDM marginals are not translation invariant")
    rdm_checks["objective_energy_residual"] <= 1e-7 ||
        error("structured and local-RDM energy conventions disagree")

    upper_certificate = build_finite_upper_certificate(
        H10_N,
        joinpath(result_dir, "certificates", "finite-upper.jls"),
    )
    validated_finite_upper_certificate(upper_certificate, "structured-RG finite upper bound")
    upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
    baseline_certificate = h11_certificate(
        structured.sohs_model,
        (tsupp=structured.canonical_words,),
        joinpath(result_dir, "certificates", "structured-sohs.jls"),
    )
    hybrid_certificate = structured_rg_certificate(
        primal_model,
        hybrid_dual,
        structured,
        structured_rdm,
        B,
        n_super,
        joinpath(result_dir, "certificates", "structured-rg.jls"),
        link_mode=link_mode,
    )
    baseline_certified = parse_exact_rational(
        baseline_certificate["certified_lower_bound_rational"],
    )
    hybrid_certified = parse_exact_rational(
        hybrid_certificate["certified_lower_bound_rational"],
    )
    baseline_certified <= upper ||
        error("structured-NPA certificate crossed the finite upper bound")
    hybrid_certified <= upper ||
        error("structured-RG certificate crossed the finite upper bound")
    raw_improvement = hybrid["energy"] - baseline["energy"]
    hybrid["energy"] + STRUCTURED_RG_ENERGY_TOLERANCE >= baseline["energy"] ||
        error("adding RG constraints lowered the minimization objective")
    hybrid["energy"] <= Float64(upper) + STRUCTURED_RG_ENERGY_TOLERANCE ||
        error("structured-RG lower estimate crossed the finite-energy upper bound")

    record = Dict(
        "meta" => Dict(
            "stage" => "structured_rg_hybrid",
            "run_name" => run_name,
            "julia_version" => string(VERSION),
            "system_size" => H10_N,
            "boundary" => "periodic",
            "n_super" => n_super,
            "support_sites" => 2 * n_super,
            "bond_dimension" => size(B, 1),
            "dyadic_bits" => RDM8_DYADIC_BITS,
            "link_mode" => string(link_mode),
            "local_link_scale" =>
                link_mode == :rdm8 ?
                    STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE :
                    STRUCTURED_RG_EQUALITY_SCALE,
            "omega4_link_scale" =>
                link_mode == :rdm8 ?
                    STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE :
                    STRUCTURED_RG_EQUALITY_SCALE,
            "recursion_scale" =>
                link_mode == :rdm8 ?
                    STRUCTURED_RG_RDM8_RECURSION_SCALE :
                    STRUCTURED_RG_EQUALITY_SCALE,
            "mosek_threads" => mosek_threads,
            "solver_feasibility_tolerance" => feasibility_tolerance,
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
            "loaded_from_cache" => coarse_grainer_loaded,
            "largest_forbidden_charge_entry" => forbidden,
            "compatibility_audit" => compatibility,
            "local_rdm_exactness_audit" => Dict(
                string(field) => getproperty(local_rdm_exactness, field)
                for field in propertynames(local_rdm_exactness)
            ),
        ),
        "model" => Dict(
            "build_wall_seconds" => build_wall,
            "canonical_moment_count" => length(structured.canonical_words),
            "rdm9_block_dimensions" => STRUCTURED_RG_RDM9_DIMENSIONS,
            "local_rdm_link_equality_count" =>
                link_mode == :rho3 ?
                    hierarchy.rho3_link_equality_count :
                    hierarchy.rho4_link_equality_count,
            "omega4_link_equality_count" =>
                link_mode == :rdm8 ?
                    hierarchy.omega4_link_equality_count :
                    0,
            "omega_count" => length(hierarchy.omega),
        ),
        "rho3_checks" => rdm_checks,
        "baseline_rho3_checks" => baseline_rdm_checks,
        "baseline" => baseline,
        "hybrid" => hybrid,
        "baseline_certificate" => baseline_certificate,
        "hybrid_certificate" => hybrid_certificate,
        "comparison" => Dict(
            "raw_improvement" => raw_improvement,
            "hybrid_raw_gap_to_upper" =>
                Float64(upper) - hybrid["energy"],
            "certified_improvement_rational" =>
                string(hybrid_certified - baseline_certified),
            "hybrid_certified_gap_rational" =>
                string(upper - hybrid_certified),
        ),
        "finite_energy_upper_certificate" => upper_certificate,
        "gates" => Dict(
            "both_optimal" => true,
            "local_rdm_exactness_passed" =>
                local_rdm_exactness.status == "PROVEN" &&
                local_rdm_exactness.exact_rational_reconstruction_matches &&
                local_rdm_exactness.coefficient_assembly_exact,
            "rg_monotonicity_passed" =>
                hybrid["energy"] + STRUCTURED_RG_ENERGY_TOLERANCE >= baseline["energy"],
            "finite_ceiling_passed" =>
                hybrid["energy"] <=
                Float64(upper) + STRUCTURED_RG_ENERGY_TOLERANCE,
            "strict_hybrid_certificate_passed" =>
                hybrid_certificate["status"] == "CERTIFIED" &&
                hybrid_certified <= upper,
            "strict_improvement_passed" =>
                hybrid_certified > baseline_certified,
            "target_passed" =>
                upper - hybrid_certified <= 1 // 100000,
        ),
    )
    validated_structured_rg_result(
        record,
        "generated structured-RG result";
        expected_hybrid_structure=structured_rg_model_structure(
            primal_model,
            hybrid_dual,
            structured,
            structured_rdm,
            B,
            n_super,
            virtual_charges,
            super_charges,
            link_mode=link_mode,
        ),
        expected_baseline_structure=h11_model_structure(
            structured.sohs_model,
            structured.canonical_words,
        ),
    )
    write_h1_json(result_path, record)
    plot_structured_rg(record, figure_path)
    @printf(
        "structured NPA %.12f -> structured NPA + RG %.12f (raw improvement %.3e)\n",
        baseline["energy"],
        hybrid["energy"],
        raw_improvement,
    )
    result_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    structured_rg_main(ARGS)
end
