# H11 — full-model exact-rational SOHS certificate for QMBCertify strengthenings.

using Random

include(joinpath(@__DIR__, "h10_qmbcertify_n20.jl"))
include(joinpath(@__DIR__, "finite_upper_certificate.jl"))

const H11_EQUALITY_SCALE = 16384.0
const H11_RANDOM_SEED = 0x48413131
const H11_EXPECTED_QMB_SOURCE_SHA256 = "bc33d5913868e4556651408495b93543e19cf704699ccabc4b5b4ad37f2ab96c"
const H11_EXPECTED_QMB_SOURCE_TREE_SHA256 = "36d3d25dc3fffc716cc057c76c140f173c3592ea03ce49eb9778b3914a25e419"

function h11_source_tree_sha256(package_root)
    sources = String[]
    for (root, _, files) in walkdir(package_root)
        for file in files
            (endswith(file, ".jl") || (root == package_root && file == "Project.toml")) || continue
            push!(sources, joinpath(root, file))
        end
    end
    sort!(sources; by=path -> relpath(path, package_root))
    inventory = join(
        [
            relpath(path, package_root) * "\0" * bytes2hex(sha256(read(path)))
            for path in sources
        ],
        "\n",
    )
    bytes2hex(sha256(inventory))
end

function install_model_returning_gsb()
    source_path = joinpath(dirname(pathof(QMBCertify)), "bound_gsp.jl")
    package_root = dirname(dirname(source_path))
    source = read(source_path, String)
    function_needle = "function GSB("
    signature_needle = "writetofile=false, mosek_setting=mosek_para())"
    return_needle = "return objv,data"
    equality_needle = "    @constraint(model, con, cons==zeros(length(tsupp)))"
    build_only_needle = equality_needle * "\n    if writetofile != false"
    scaled_equality = "    @constraint(model, con, ($(H11_EQUALITY_SCALE) .* cons)==zeros(length(tsupp)))"
    length(findall(function_needle, source)) == 1 || error("unexpected QMBCertify GSB source layout")
    length(findall(signature_needle, source)) == 1 || error("unexpected QMBCertify GSB signature layout")
    length(findall(return_needle, source)) == 1 || error("unexpected QMBCertify GSB return layout")
    length(findall(equality_needle, source)) == 1 || error("unexpected QMBCertify equality source layout")
    instrumented = replace(source, function_needle => "function GSB_with_model("; count=1)
    instrumented = replace(
        instrumented,
        signature_needle => "writetofile=false, mosek_setting=mosek_para(), h11_build_only=false)";
        count=1,
    )
    instrumented = replace(instrumented, return_needle => "return objv,data,model"; count=1)
    instrumented = replace(
        instrumented,
        build_only_needle => scaled_equality * "\n    if h11_build_only\n        return nothing, (tsupp=tsupp,), model\n    end\n    if writetofile != false";
        count=1,
    )
    Base.include_string(QMBCertify, instrumented, source_path * "#model-returning")
    source_sha256 = bytes2hex(sha256(source))
    source_tree_sha256 = h11_source_tree_sha256(package_root)
    source_sha256 == H11_EXPECTED_QMB_SOURCE_SHA256 || error("unexpected QMBCertify bound_gsp.jl source")
    source_tree_sha256 == H11_EXPECTED_QMB_SOURCE_TREE_SHA256 ||
        error("unexpected QMBCertify source tree")
    source_path, source_sha256, source_tree_sha256
end

const H11_QMB_SOURCE_PATH, H11_QMB_SOURCE_SHA256, H11_QMB_SOURCE_TREE_SHA256 =
    install_model_returning_gsb()

function parse_h11_args(args)
    run_name = nothing
    smoke = false
    mosek_threads = 16
    feasibility_tolerance = 1e-8
    upper_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--upper-result=")
            spec = split(arg, "=", limit=2)[2]
            upper_result = normpath(isabspath(spec) ? spec : joinpath(BOUNDGSENERGY_ROOT, spec))
        elseif arg == "--smoke"
            smoke = true
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--feasibility-tolerance=")
            feasibility_tolerance = parse(Float64, split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    mosek_threads >= 1 || error("mosek_threads must be positive")
    0 < feasibility_tolerance <= 1e-8 || error("feasibility_tolerance must lie in (0, 1e-8]")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h11-full-sohs")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    !isnothing(upper_result) && !isfile(upper_result) &&
        error("finite-energy upper result does not exist: $upper_result")
    run_name, smoke, mosek_threads, feasibility_tolerance, upper_result
end

function h11_linear_terms(expression, variable_indices)
    sort!(
        [(variable_indices[variable], Float64(coefficient)) for (coefficient, variable) in linear_terms(expression)];
        by=first,
    )
end

function h11_objective_transcript(model, variable_indices)
    objective = objective_function(model)
    if objective isa VariableRef
        return (constant=0.0, terms=[(variable_indices[objective], 1.0)])
    elseif objective isa AffExpr
        return (
            constant=Float64(objective.constant),
            terms=h11_linear_terms(objective, variable_indices),
        )
    end
    error("unsupported H11 objective type $(typeof(objective))")
end

function h11_model_structure(model, tsupp)
    supported_constraint_types = Set([
        (Vector{AffExpr}, MOI.Zeros),
        (Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle),
    ])
    all(type -> type in supported_constraint_types, list_of_constraint_types(model)) ||
        error("QMBCertify model contains an unsupported constraint type")
    variables = all_variables(model)
    variable_indices = Dict(variable => index for (index, variable) in enumerate(variables))
    equality_constraints = all_constraints(model, Vector{AffExpr}, MOI.Zeros)
    length(equality_constraints) == 1 || error("expected one QMBCertify coefficient-matching constraint")
    equality_object = constraint_object(only(equality_constraints))
    equalities = [
        let normalized = expression / H11_EQUALITY_SCALE
            (
                constant=0.0,
                rhs=-Float64(normalized.constant),
                terms=h11_linear_terms(normalized, variable_indices),
            )
        end
        for expression in equality_object.func
    ]
    length(equalities) == length(tsupp) ||
        error("QMBCertify equality count does not match its canonical word support")

    psd_blocks = NamedTuple[]
    covered_variables = Set{Int}()
    for constraint in all_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
        object = constraint_object(constraint)
        matrix = reshape_vector(object.func, object.shape)
        indices = [variable_indices[matrix[row, column]] for row in axes(matrix, 1), column in axes(matrix, 2)]
        block_variables = Set(indices)
        isempty(intersect(covered_variables, block_variables)) ||
            error("a QMBCertify scalar variable belongs to multiple PSD blocks")
        union!(covered_variables, block_variables)
        base_name = first(split(name(matrix[1, 1]), "[", limit=2))
        push!(psd_blocks, (group=base_name, dimension=size(matrix, 1), variable_indices=indices))
    end
    isempty(psd_blocks) && error("QMBCertify model has no PSD blocks")
    (
        version=3,
        qmbcertify_source_sha256=H11_QMB_SOURCE_SHA256,
        qmbcertify_source_tree_sha256=H11_QMB_SOURCE_TREE_SHA256,
        equality_scale=H11_EQUALITY_SCALE,
        variable_count=length(variables),
        objective_sense=string(objective_sense(model)),
        objective=h11_objective_transcript(model, variable_indices),
        equalities=equalities,
        psd_blocks=psd_blocks,
        canonical_words=[Int.(word) for word in tsupp],
    )
end

function h11_sohs_transcript(model, tsupp)
    structure = h11_model_structure(model, tsupp)
    merge(structure, (raw_values=Float64[value(variable) for variable in all_variables(model)],))
end

function h11_configuration(smoke)
    smoke ?
        (order=3, extra=0, rdm=0, pso=0, lso=false) :
        (order=4, extra=1, rdm=9, pso=3, lso=true)
end

function h11_rebuild_model_structure(configuration)
    Random.seed!(H11_RANDOM_SEED)
    output = redirect_stdout(devnull) do
        QMBCertify.GSB_with_model(
            [[1, 4]],
            [3 / 4],
            H10_N,
            configuration.order;
            lattice="chain",
            rdm=configuration.rdm,
            extra=configuration.extra,
            pso=configuration.pso,
            lso=configuration.lso,
            lol=H10_N,
            three_type=[1, 1],
            Gram=true,
            QUIET=true,
            mosek_setting=QMBCertify.mosek_para(1e-8, 1e-8, 1e-8, 1),
            h11_build_only=true,
        )
    end
    _, data, model = output
    h11_model_structure(model, data.tsupp)
end

const H11_MODEL_STRUCTURE_FIELDS = (
    :version,
    :qmbcertify_source_sha256,
    :qmbcertify_source_tree_sha256,
    :equality_scale,
    :variable_count,
    :objective_sense,
    :objective,
    :equalities,
    :psd_blocks,
    :canonical_words,
)

function h11_validate_model_binding(transcript, expected)
    for field in H11_MODEL_STRUCTURE_FIELDS
        hasproperty(transcript, field) || error("H11 transcript is missing $field")
        getproperty(transcript, field) == getproperty(expected, field) ||
            error("H11 transcript does not match the rebuilt QMBCertify model at $field")
    end
    true
end

function h11_evaluate_linear(expression, values)
    exact_rational(expression.constant) +
        sum(exact_rational(coefficient) * values[index] for (index, coefficient) in expression.terms)
end

function h11_evaluate_transcript(transcript, expected_structure; eig_prec=256)
    h11_validate_model_binding(transcript, expected_structure)
    transcript.version == 3 || error("unsupported H11 transcript version")
    transcript.qmbcertify_source_sha256 == H11_QMB_SOURCE_SHA256 ||
        error("QMBCertify source hash differs from the certificate transcript")
    transcript.qmbcertify_source_tree_sha256 == H11_QMB_SOURCE_TREE_SHA256 ||
        error("QMBCertify source-tree hash differs from the certificate transcript")
    scale = transcript.equality_scale
    scale == H11_EQUALITY_SCALE || error("H11 equality scale differs from the pinned model")
    transcript.objective_sense == string(MOI.MAX_SENSE) ||
        error("H11 objective sense is not maximization")
    length(transcript.raw_values) == transcript.variable_count || error("H11 transcript has stale variable values")
    length(transcript.equalities) == length(transcript.canonical_words) ||
        error("H11 transcript has stale canonical words")
    for word in transcript.canonical_words
        reduced, coefficient = QMBCertify.reduce!(
            UInt16.(word);
            L=H10_N,
            lattice="chain",
            realify=true,
        )
        Int.(reduced) == word && coefficient == 1 ||
            error("H11 transcript contains a noncanonical Pauli word")
    end

    values = exact_rational.(transcript.raw_values)
    groups = Dict{String, Any}()
    for (block_index, block) in enumerate(transcript.psd_blocks)
        matrix = [values[block.variable_indices[row, column]] for row in 1:block.dimension, column in 1:block.dimension]
        matrix == matrix' || error("H11 transcript contains a nonsymmetric PSD block")
        eigenvalue_lower = first(QMBCertify.rigorous_min_eig_bound(Symmetric(matrix); prec=eig_prec))
        diagonal_shift = -min(ExactRational(0), eigenvalue_lower)
        for diagonal in 1:block.dimension
            values[block.variable_indices[diagonal, diagonal]] += diagonal_shift
        end
        groups["$(block.group)#$block_index"] = Dict(
            "dimension" => block.dimension,
            "minimum_eigenvalue_lower_bound" => float_down(eigenvalue_lower),
            "minimum_eigenvalue_lower_bound_rational" => string(eigenvalue_lower),
            "diagonal_shift" => float_up(diagonal_shift),
            "diagonal_shift_rational" => string(diagonal_shift),
        )
    end

    residuals = ExactRational[
        h11_evaluate_linear(equality, values) - exact_rational(equality.rhs)
        for equality in transcript.equalities
    ]
    residual_l1 = sum(abs, residuals)
    objective = h11_evaluate_linear(transcript.objective, values)
    certified_lower_bound = objective - residual_l1
    certified_lower_bound, objective, residual_l1, residuals, groups
end

function h11_certificate(model, data, transcript_path)
    transcript = h11_sohs_transcript(model, data.tsupp)
    atomic_serialize(transcript_path, transcript)
    certified, objective, residual_l1, residuals, groups = h11_evaluate_transcript(
        transcript,
        h11_model_structure(model, data.tsupp),
    )
    published = decimal_floor(certified; digits=9)
    Dict(
        "status" => "CERTIFIED",
        "method" => "full-model exact-dyadic SOHS replay with rigorous PSD repair and Pauli-word l1 residual correction",
        "julia_version" => string(VERSION),
        "operator_norm_argument" => "every canonical Pauli word has operator norm one",
        "qmbcertify_source_sha256" => H11_QMB_SOURCE_SHA256,
        "qmbcertify_source_tree_sha256" => H11_QMB_SOURCE_TREE_SHA256,
        "replay_transcript" => artifact_relative_path(transcript_path),
        "replay_transcript_sha256" => file_sha256(transcript_path),
        "raw_objective" => float_down(objective),
        "raw_objective_rational" => string(objective),
        "identity_residual_l1" => float_up(residual_l1),
        "identity_residual_l1_rational" => string(residual_l1),
        "maximum_identity_residual" => float_up(maximum(abs, residuals)),
        "certified_lower_bound" => float_down(certified),
        "certified_lower_bound_rational" => string(certified),
        "published_certified_lower_bound" => float_down(published),
        "published_certified_lower_bound_rational" => string(published),
        "groups" => groups,
    )
end

function validated_h11_certificate(record, label; expected_structure=nothing)
    record["meta"]["stage"] == "H11" || error("$label stage is stale")
    record["meta"]["system_size"] == H10_N || error("$label system size is stale")
    record["meta"]["boundary"] == "periodic" || error("$label boundary is stale")
    record["meta"]["coefficient_matching_scale"] == H11_EQUALITY_SCALE ||
        error("$label equality scale is stale")
    record["meta"]["julia_version"] == string(PINNED_JULIA_VERSION) ||
        error("$label Julia version is stale")
    record["meta"]["random_seed"] == string(H11_RANDOM_SEED) || error("$label random seed is stale")
    record["meta"]["qmbcertify_source_sha256"] == H11_QMB_SOURCE_SHA256 ||
        error("$label QMBCertify source hash is stale")
    record["meta"]["qmbcertify_source_tree_sha256"] == H11_QMB_SOURCE_TREE_SHA256 ||
        error("$label QMBCertify source-tree hash is stale")
    configuration = h11_configuration(record["meta"]["smoke_test"])
    expected_configuration = Dict(string(key) => value for (key, value) in pairs(configuration))
    record["meta"]["configuration"] == expected_configuration || error("$label configuration is stale")
    record["solve"]["status"] == "OPTIMAL" || error("$label solve did not terminate OPTIMAL")
    record["solve"]["primal_status"] == "FEASIBLE_POINT" ||
        error("$label solve did not return a feasible primal SOHS point")
    certificate = record["certificate"]
    certificate["status"] == "CERTIFIED" || error("$label certificate did not pass")
    certificate["julia_version"] == string(PINNED_JULIA_VERSION) ||
        error("$label certificate Julia version is stale")
    certificate["qmbcertify_source_sha256"] == H11_QMB_SOURCE_SHA256 ||
        error("$label certificate source hash is stale")
    certificate["qmbcertify_source_tree_sha256"] == H11_QMB_SOURCE_TREE_SHA256 ||
        error("$label certificate source-tree hash is stale")
    path = normpath(joinpath(BOUNDGSENERGY_ROOT, certificate["replay_transcript"]))
    isfile(path) || error("$label transcript does not exist")
    file_sha256(path) == certificate["replay_transcript_sha256"] || error("$label transcript hash is stale")
    transcript = deserialize(path)
    rebuilt_structure = isnothing(expected_structure) ?
        h11_rebuild_model_structure(configuration) :
        expected_structure
    certified, objective, residual_l1, residuals, groups = h11_evaluate_transcript(
        transcript,
        rebuilt_structure,
    )
    exact_rational(record["solve"]["numeric_energy_per_site"]) == objective ||
        error("$label numeric objective is stale")
    record["solve"]["scalar_variables"] == transcript.variable_count ||
        error("$label variable count is stale")
    record["solve"]["constraint_count"] == length(transcript.equalities) ||
        error("$label constraint count is stale")
    record["solve"]["largest_psd_block"] == maximum(block.dimension for block in transcript.psd_blocks) ||
        error("$label largest PSD block is stale")
    parse_exact_rational(certificate["raw_objective_rational"]) == objective || error("$label raw objective is stale")
    parse_exact_rational(certificate["identity_residual_l1_rational"]) == residual_l1 ||
        error("$label residual correction is stale")
    parse_exact_rational(certificate["certified_lower_bound_rational"]) == certified ||
        error("$label certified lower bound is stale")
    published = decimal_floor(certified; digits=9)
    parse_exact_rational(certificate["published_certified_lower_bound_rational"]) == published ||
        error("$label published lower bound is stale")
    certificate["published_certified_lower_bound"] == float_down(published) ||
        error("$label published lower-bound summary is stale")
    certificate["maximum_identity_residual"] == float_up(maximum(abs, residuals)) ||
        error("$label maximum residual is stale")
    certificate["groups"] == groups || error("$label PSD repair summary is stale")
    upper_certificate = validated_finite_upper_certificate(
        record["finite_energy_upper_certificate"],
        "$label finite upper bound",
    )
    upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
    gap = upper - certified
    gap >= 0 || error("$label lower bound exceeds its finite upper bound")
    parse_exact_rational(record["comparison"]["certified_gap_rational"]) == gap ||
        error("$label certified gap is stale")
    record["comparison"]["certified_gap"] == float_up(gap) ||
        error("$label certified-gap summary is stale")
    target = 1 // 100000
    record["gates"]["sandwich_passed"] == (gap >= 0) || error("$label sandwich gate is stale")
    record["gates"]["target_passed"] == (gap <= target) || error("$label target gate is stale")
    record["gates"]["stage_passed"] == (gap <= target) || error("$label stage gate is stale")
    certificate
end

function h11_main(args)
    run_name, smoke, mosek_threads, feasibility_tolerance, upper_result_path = parse_h11_args(args)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_path = joinpath(result_dir, "full_sohs_certificate.json")
    log_path = joinpath(result_dir, "qmbcertify.log")
    figure_path = joinpath(result_dir, "figs", "full_sohs_gap.png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")
    mkpath(result_dir)

    configuration = h11_configuration(smoke)
    supp = [[1, 4]]
    coe = [3 / 4]
    output = Ref{Any}()
    Random.seed!(H11_RANDOM_SEED)
    solve_wall = @elapsed open(log_path, "w") do io
        redirect_stdout(io) do
            output[] = QMBCertify.GSB_with_model(
                supp,
                coe,
                H10_N,
                configuration.order;
                lattice="chain",
                rdm=configuration.rdm,
                extra=configuration.extra,
                pso=configuration.pso,
                lso=configuration.lso,
                lol=H10_N,
                three_type=[1, 1],
                Gram=true,
                QUIET=false,
                mosek_setting=QMBCertify.mosek_para(
                    feasibility_tolerance,
                    feasibility_tolerance,
                    feasibility_tolerance,
                    mosek_threads,
                ),
            )
        end
    end
    log_text = read(log_path, String)
    print(log_text)
    optimum, data, model = output[]
    termination_status(model) == MOI.OPTIMAL ||
        error("QMBCertify did not terminate OPTIMAL; see $log_path")
    primal_status(model) == MOI.FEASIBLE_POINT ||
        error("QMBCertify did not return a feasible primal SOHS point; see $log_path")
    certificate_wall = @elapsed certificate = h11_certificate(
        model,
        data,
        joinpath(result_dir, "certificates", "full-sohs.jls"),
    )
    certified_lower = parse_exact_rational(certificate["certified_lower_bound_rational"])
    upper_certificate = if isnothing(upper_result_path)
        build_finite_upper_certificate(
            H10_N,
            joinpath(result_dir, "certificates", "finite-upper.jls"),
        )
    else
        JSON.parsefile(upper_result_path)["finite_energy_upper_certificate"]
    end
    validated_finite_upper_certificate(upper_certificate, "H11 finite upper bound")
    certified_upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
    gap = certified_upper - certified_lower
    gap >= 0 || error("H11 lower bound exceeds the exact Rayleigh upper bound")
    target_gap = 1 // 100000
    largest_psd_block = maximum(group["dimension"] for group in values(certificate["groups"]))

    results = Dict(
        "meta" => Dict(
            "stage" => "H11",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H10_N,
            "boundary" => "periodic",
            "configuration" => Dict(string(key) => value for (key, value) in pairs(configuration)),
            "smoke_test" => smoke,
            "solver_feasibility_tolerance" => feasibility_tolerance,
            "coefficient_matching_scale" => H11_EQUALITY_SCALE,
            "julia_version" => string(VERSION),
            "random_seed" => string(H11_RANDOM_SEED),
            "mosek_threads" => mosek_threads,
            "qmbcertify_source_sha256" => H11_QMB_SOURCE_SHA256,
            "qmbcertify_source_tree_sha256" => H11_QMB_SOURCE_TREE_SHA256,
            "inputs" => isnothing(upper_result_path) ?
                Dict("finite_energy_upper_result" => nothing) :
                Dict("finite_energy_upper_result" => artifact_relative_path(upper_result_path)),
        ),
        "solve" => Dict(
            "status" => "OPTIMAL",
            "primal_status" => "FEASIBLE_POINT",
            "numeric_energy_per_site" => Float64(optimum),
            "solve_wall_seconds" => solve_wall,
            "postcertificate_wall_seconds" => certificate_wall,
            "largest_psd_block" => largest_psd_block,
            "constraint_count" => length(data.tsupp),
            "scalar_variables" => num_variables(model),
            "log" => basename(log_path),
        ),
        "certificate" => certificate,
        "finite_energy_upper_certificate" => upper_certificate,
        "comparison" => Dict(
            "certified_gap" => float_up(gap),
            "certified_gap_rational" => string(gap),
        ),
        "gates" => Dict(
            "target_gap" => Float64(target_gap),
            "sandwich_passed" => gap >= 0,
            "target_passed" => gap <= target_gap,
            "stage_passed" => gap <= target_gap,
        ),
    )
    write_results(result_path, results)
    validated_h11_certificate(results, "fresh H11")

    gap_plot = bar(
        ["certified gap"],
        [float_up(gap)];
        ylabel="energy density",
        title=smoke ? "H11 smoke certificate" : "H11 full strengthened certificate",
        legend=false,
        color=:steelblue,
        ylims=(0, 1.12 * max(float_up(gap), Float64(target_gap))),
    )
    hline!(gap_plot, [Float64(target_gap)]; color=:black, linestyle=:dash, label=false)
    mkpath(dirname(figure_path))
    savefig(gap_plot, figure_path)

    @printf("numeric %.12f, certified %.12f\n", optimum, Float64(certified_lower))
    @printf("identity residual l1 %.3e\n", certificate["identity_residual_l1"])
    @printf("certified gap %.3e, target passed=%s\n", Float64(gap), results["gates"]["target_passed"])
    @printf("solve %.2fs, postcertificate %.2fs\n", solve_wall, certificate_wall)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h11_main(ARGS)
end
