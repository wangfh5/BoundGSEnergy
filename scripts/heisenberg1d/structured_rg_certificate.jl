# Exact mixed SOHS/RG certificate for the structured-NPA + RG formulation.

include(joinpath(@__DIR__, "structured_rg_hybrid.jl"))

const STRUCTURED_RG_ENERGY_TOLERANCE = 2e-7

function structured_rg_constraint_expression(constraint, component)
    object = constraint_object(constraint)
    if object.func isa AffExpr && object.set isa MOI.EqualTo{Float64}
        component in (-1, 1) ||
            error("scalar dual constraint has an invalid component")
        return object.func, object.set.value
    elseif object.func isa Vector{AffExpr} && object.set isa MOI.Zeros
        component in eachindex(object.func) ||
            error("vector dual constraint has an invalid component")
        return object.func[component], 0.0
    end
    error("unsupported structured-RG moment constraint")
end

function structured_rg_affine_transcript(
    expression,
    rhs,
    variable_indices,
)
    (
        rhs=Float64(rhs),
        constant=Float64(expression.constant),
        terms=h11_linear_terms(expression, variable_indices),
    )
end

function structured_rg_psd_variable_blocks(model, variable_indices)
    blocks = NamedTuple[]
    covered = Set{Int}()
    for constraint in
        all_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
        object = constraint_object(constraint)
        dimension = MOI.side_dimension(object.set)
        matrix = reshape_vector(
            object.func,
            SymmetricMatrixShape(dimension),
        )
        indices = [
            variable_indices[matrix[row, column]]
            for row in axes(matrix, 1), column in axes(matrix, 2)
        ]
        block_variables = Set(indices)
        isempty(intersect(covered, block_variables)) ||
            error("a structured-RG dual variable belongs to multiple PSD blocks")
        union!(covered, block_variables)
        push!(
            blocks,
            (
                group=first(split(name(matrix[1, 1]), "[", limit=2)),
                dimension=dimension,
                variable_indices=indices,
            ),
        )
    end
    isempty(blocks) && error("structured-RG dual has no SOHS PSD blocks")
    blocks
end

function structured_rg_omega_slacks(
    primal_model,
    dual_model,
    variable_indices,
)
    slacks = NamedTuple[]
    for primal_constraint in
        all_constraints(primal_model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
        primal_object = constraint_object(primal_constraint)
        primal_block = reshape_vector(primal_object.func, primal_object.shape)
        dual_constraint, _ =
            Dualization._get_dual_constraint(dual_model, primal_block[1, 1])
        dual_constraint isa ConstraintRef ||
            error("RG PSD block is not represented by a dual slack constraint")
        dual_object = constraint_object(dual_constraint)
        dual_object.func isa Vector{AffExpr} ||
            error("RG dual slack is not affine")
        dimension = size(primal_block, 1)
        matrix = reshape_vector(
            dual_object.func,
            SymmetricMatrixShape(dimension),
        )
        entries = NamedTuple[]
        for row in 1:dimension
            for column in row:dimension
                push!(
                    entries,
                    structured_rg_affine_transcript(
                        matrix[row, column],
                        0.0,
                        variable_indices,
                    ),
                )
            end
        end
        push!(
            slacks,
            (
                group=psd_group_name(primal_constraint),
                dimension=dimension,
                upper_entries=entries,
            ),
        )
    end
    isempty(slacks) && error("structured-RG primal has no RG PSD blocks")
    slacks
end

function structured_rg_model_structure(
    primal_model,
    dual_model,
    structured,
    structured_local_rdm,
    B,
    n_super,
    virtual_charges=STRUCTURED_RG_VIRTUAL_CHARGES,
    super_charges=STRUCTURED_RG_SUPER_CHARGES,
    ;
    link_mode=:rho3,
)
    link_mode in (:rho3, :rdm8) ||
        error("unsupported structured-RG link mode")
    site_count = link_mode == :rho3 ? 6 : 8
    objective_sense(dual_model) == MOI.MAX_SENSE ||
        error("structured-RG certificate requires a maximization dual")
    variables = all_variables(dual_model)
    variable_indices =
        Dict(variable => index for (index, variable) in enumerate(variables))
    moment_equalities = [
        let
            constraint, component =
                Dualization._get_dual_constraint(dual_model, moment)
            constraint isa ConstraintRef ||
                error("structured moment is not represented in the hybrid dual")
            expression, rhs =
                structured_rg_constraint_expression(constraint, component)
            structured_rg_affine_transcript(
                expression,
                rhs,
                variable_indices,
            )
        end
        for moment in structured.scaled_moments
    ]
    length(moment_equalities) == length(structured.canonical_words) ||
        error("structured-RG moment equalities do not match canonical words")
    (
        version=link_mode == :rho3 ? 3 : 4,
        julia_version=string(VERSION),
        qmbcertify_source_sha256=H11_QMB_SOURCE_SHA256,
        qmbcertify_source_tree_sha256=H11_QMB_SOURCE_TREE_SHA256,
        system_size=H10_N,
        n_super=n_super,
        local_rdm_exactness_audit=structured_local_rdm_exactness_audit(
            structured_local_rdm,
            structured.physical_moments,
            structured.canonical_words,
            site_count,
        ),
        dyadic_bits=RDM8_DYADIC_BITS,
        equality_scale=H11_EQUALITY_SCALE,
        rg_equality_scale=STRUCTURED_RG_EQUALITY_SCALE,
        local_link_scale=link_mode == :rdm8 ?
            STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE :
            STRUCTURED_RG_EQUALITY_SCALE,
        omega4_link_scale=link_mode == :rdm8 ?
            STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE :
            STRUCTURED_RG_EQUALITY_SCALE,
        recursion_scale=link_mode == :rdm8 ?
            STRUCTURED_RG_RDM8_RECURSION_SCALE :
            STRUCTURED_RG_EQUALITY_SCALE,
        coarse_grainer_fingerprint=bytes2hex(sha1(reinterpret(UInt8, vec(B)))),
        coarse_grainer=copy(B),
        virtual_charges=collect(virtual_charges),
        super_charges=collect(super_charges),
        variable_count=length(variables),
        objective_sense=string(objective_sense(dual_model)),
        objective=h11_objective_transcript(dual_model, variable_indices),
        canonical_words=[Int.(word) for word in structured.canonical_words],
        moment_equalities=moment_equalities,
        psd_variable_blocks=
            structured_rg_psd_variable_blocks(dual_model, variable_indices),
        omega_slacks=structured_rg_omega_slacks(
            primal_model,
            dual_model,
            variable_indices,
        ),
    )
end

const STRUCTURED_RG_MODEL_STRUCTURE_FIELDS = (
    :version,
    :julia_version,
    :qmbcertify_source_sha256,
    :qmbcertify_source_tree_sha256,
    :system_size,
    :n_super,
    :local_rdm_exactness_audit,
    :dyadic_bits,
    :equality_scale,
    :rg_equality_scale,
    :local_link_scale,
    :omega4_link_scale,
    :recursion_scale,
    :coarse_grainer_fingerprint,
    :coarse_grainer,
    :virtual_charges,
    :super_charges,
    :variable_count,
    :objective_sense,
    :objective,
    :canonical_words,
    :moment_equalities,
    :psd_variable_blocks,
    :omega_slacks,
)

function structured_rg_validate_coarse_grainer(structure)
    B = structure.coarse_grainer
    size(B) == (4, 4, 4) ||
        error("structured-RG coarse-grainer has an unexpected shape")
    structure.dyadic_bits == RDM8_DYADIC_BITS ||
        error("structured-RG coarse-grainer uses the wrong dyadic grid")
    bytes2hex(sha1(reinterpret(UInt8, vec(B)))) ==
    structure.coarse_grainer_fingerprint ||
        error("structured-RG coarse-grainer fingerprint is stale")
    exact_B, largest_forbidden = exact_charge_tensor(
        B,
        structure.virtual_charges,
        structure.super_charges,
    )
    largest_forbidden == 0 && exact_B == B ||
        error("structured-RG coarse-grainer violates exact U(1) support")
    all(iszero, values(assert_charge_preserving_maps(
        B,
        structure.virtual_charges,
        structure.super_charges,
    ))) || error("structured-RG maps do not preserve charge exactly")
    audit = dyadic_relaxation_audit(B, structure.dyadic_bits)
    audit["physical_relaxation_compatibility_proven"] ||
        error("structured-RG coarse-grainer failed the exact compatibility audit")
    true
end

function structured_rg_validate_local_rdm_audit(structure)
    audit = structure.local_rdm_exactness_audit
    audit.status == "PROVEN" ||
        error("structured local-RDM exactness audit did not pass")
    expected_site_count = structure.version == 3 ? 6 :
        structure.version == 4 ? 8 :
        error("unsupported structured-RG transcript version")
    audit.site_count == expected_site_count ||
        error("structured local-RDM audit has the wrong support")
    audit.grid_denominator == string(big(2)^audit.site_count) ||
        error("structured local-RDM audit has the wrong dyadic grid")
    audit.physical_moment_scale == string(exact_rational(H11_EQUALITY_SCALE)) ||
        error("structured local-RDM audit has the wrong moment scale")
    audit.coefficient_assembly_exact ||
        error("structured local-RDM coefficient assembly is not exact")
    audit.exact_rational_reconstruction_matches ||
        error("structured local-RDM exact reconstruction does not match")
    parse(BigInt, audit.accumulation_numerator_bound) <
    parse(BigInt, audit.float64_exact_integer_limit) ||
        error("structured local-RDM accumulation exceeds the exact-integer range")
    parse(BigInt, audit.maximum_final_numerator) <
    parse(BigInt, audit.float64_exact_integer_limit) ||
        error("structured local-RDM coefficient exceeds the exact-integer range")
    true
end

function structured_rg_validate_model_binding(transcript, expected)
    for field in STRUCTURED_RG_MODEL_STRUCTURE_FIELDS
        hasproperty(transcript, field) ||
            error("structured-RG transcript is missing $field")
        getproperty(transcript, field) == getproperty(expected, field) ||
            error("structured-RG transcript does not match the rebuilt model at $field")
    end
    structured_rg_validate_coarse_grainer(transcript)
    structured_rg_validate_local_rdm_audit(transcript)
end

function structured_rg_rebuild_model_structure(transcript)
    structured_rg_validate_coarse_grainer(transcript)
    transcript.system_size == H10_N ||
        error("structured-RG transcript has the wrong finite system size")
    validate_structured_rg_support(transcript.n_super)
    structured = build_structured_moment_model()
    site_count = transcript.version == 3 ? 6 :
        transcript.version == 4 ? 8 :
        error("unsupported structured-RG transcript version")
    local_rdm = structured_local_rdm(
        structured.physical_moments,
        structured.canonical_words,
        site_count,
    )
    if transcript.version == 3
        add_rg_hierarchy!(
            structured.model,
            local_rdm,
            transcript.coarse_grainer,
            transcript.n_super,
            transcript.virtual_charges,
            transcript.super_charges,
        )
    else
        add_rdm8_rg_hierarchy!(
            structured.model,
            local_rdm,
            transcript.coarse_grainer,
            transcript.n_super,
            transcript.virtual_charges,
            transcript.super_charges,
            local_link_scale=transcript.local_link_scale,
            omega4_link_scale=transcript.omega4_link_scale,
            recursion_scale=transcript.recursion_scale,
        )
    end
    dual_model = dualize(structured.model)
    structured_rg_model_structure(
        structured.model,
        dual_model,
        structured,
        local_rdm,
        transcript.coarse_grainer,
        transcript.n_super,
        transcript.virtual_charges,
        transcript.super_charges,
        link_mode=transcript.version == 3 ? :rho3 : :rdm8,
    )
end

function structured_rg_evaluate_affine(transcript, values)
    exact_rational(transcript.constant) +
    sum(
        exact_rational(coefficient) * values[index]
        for (index, coefficient) in transcript.terms;
        init=ExactRational(0),
    ) - exact_rational(transcript.rhs)
end

function structured_rg_evaluate_slack(slack, values)
    matrix = Matrix{ExactRational}(
        undef,
        slack.dimension,
        slack.dimension,
    )
    cursor = 1
    for row in 1:slack.dimension, column in row:slack.dimension
        entry = structured_rg_evaluate_affine(
            slack.upper_entries[cursor],
            values,
        )
        cursor += 1
        matrix[row, column] = entry
        matrix[column, row] = entry
    end
    matrix
end

function evaluate_structured_rg_transcript(
    transcript,
    expected_structure,
    trace_bounds;
    eig_prec=256,
)
    structured_rg_validate_model_binding(transcript, expected_structure)
    transcript.version in (3, 4) ||
        error("unsupported structured-RG transcript version")
    transcript.julia_version == string(PINNED_JULIA_VERSION) ||
        error("structured-RG transcript has the wrong Julia version")
    transcript.qmbcertify_source_sha256 == H11_QMB_SOURCE_SHA256 ||
        error("structured-RG transcript has the wrong QMBCertify source")
    transcript.qmbcertify_source_tree_sha256 ==
    H11_QMB_SOURCE_TREE_SHA256 ||
        error("structured-RG transcript has the wrong QMBCertify source tree")
    transcript.objective_sense == string(MOI.MAX_SENSE) ||
        error("structured-RG transcript is not a maximization problem")
    transcript.equality_scale == H11_EQUALITY_SCALE ||
        error("structured-RG moment scale is stale")
    transcript.rg_equality_scale == STRUCTURED_RG_EQUALITY_SCALE ||
        error("structured-RG equality scale is stale")
    expected_scales = transcript.version == 3 ?
        (
            STRUCTURED_RG_EQUALITY_SCALE,
            STRUCTURED_RG_EQUALITY_SCALE,
            STRUCTURED_RG_EQUALITY_SCALE,
        ) :
        (
            STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE,
            STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE,
            STRUCTURED_RG_RDM8_RECURSION_SCALE,
        )
    (
        transcript.local_link_scale,
        transcript.omega4_link_scale,
        transcript.recursion_scale,
    ) == expected_scales ||
        error("structured-RG link scales are stale")
    length(transcript.raw_values) == transcript.variable_count ||
        error("structured-RG transcript has stale variable values")
    length(transcript.moment_equalities) ==
    length(transcript.canonical_words) ||
        error("structured-RG transcript has stale canonical words")

    values = exact_rational.(transcript.raw_values)
    sohs_groups = Dict{String, Any}()
    for (block_index, block) in
        enumerate(transcript.psd_variable_blocks)
        size(block.variable_indices) == (block.dimension, block.dimension) ||
            error("structured-RG transcript has a malformed PSD block")
        matrix = [
            values[block.variable_indices[row, column]]
            for row in 1:block.dimension, column in 1:block.dimension
        ]
        matrix == matrix' ||
            error("structured-RG transcript contains a nonsymmetric SOHS block")
        eigenvalue_lower = first(
            QMBCertify.rigorous_min_eig_bound(
                Symmetric(matrix);
                prec=eig_prec,
            ),
        )
        diagonal_shift = -min(ExactRational(0), eigenvalue_lower)
        for diagonal in 1:block.dimension
            values[block.variable_indices[diagonal, diagonal]] +=
                diagonal_shift
        end
        sohs_groups["$(block.group)#$block_index"] = Dict(
            "dimension" => block.dimension,
            "minimum_eigenvalue_lower_bound_rational" =>
                string(eigenvalue_lower),
            "diagonal_shift_rational" => string(diagonal_shift),
        )
    end

    moment_residuals = ExactRational[
        structured_rg_evaluate_affine(equality, values) /
        exact_rational(transcript.equality_scale)
        for equality in transcript.moment_equalities
    ]
    moment_correction = sum(abs, moment_residuals)

    slack_eigenvalues = Dict{String, Vector{ExactRational}}()
    for slack in transcript.omega_slacks
        matrix = structured_rg_evaluate_slack(slack, values)
        eigenvalue_lower = first(
            QMBCertify.rigorous_min_eig_bound(
                Symmetric(matrix);
                prec=eig_prec,
            ),
        )
        push!(
            get!(slack_eigenvalues, slack.group, ExactRational[]),
            eigenvalue_lower,
        )
    end
    Set(keys(slack_eigenvalues)) == Set(keys(trace_bounds)) ||
        error("structured-RG trace bounds do not match RG slack groups")
    rg_correction = ExactRational(0)
    rg_groups = Dict{String, Any}()
    for group in sort!(collect(keys(slack_eigenvalues)))
        minimum_eigenvalue = minimum(slack_eigenvalues[group])
        correction =
            min(ExactRational(0), minimum_eigenvalue) * trace_bounds[group]
        rg_correction += correction
        rg_groups[group] = Dict(
            "block_count" => length(slack_eigenvalues[group]),
            "minimum_eigenvalue_lower_bound_rational" =>
                string(minimum_eigenvalue),
            "trace_upper_bound_rational" => string(trace_bounds[group]),
            "correction_rational" => string(correction),
        )
    end

    objective =
        h11_evaluate_linear(transcript.objective, values)
    legacy_certified =
        objective - moment_correction + rg_correction
    sharpened = h11_energy_aware_residual_bound(
        objective,
        moment_residuals,
        transcript.canonical_words;
        scalar_correction=rg_correction,
    )
    sharpened.legacy_lower_bound == legacy_certified ||
        error("structured-RG legacy residual bound is inconsistent")
    (
        certified=sharpened.certified_lower_bound,
        legacy_certified=legacy_certified,
        objective=objective,
        moment_correction=moment_correction,
        moment_residuals=moment_residuals,
        rg_correction=rg_correction,
        residual_bound=sharpened,
        sohs_groups=sohs_groups,
        rg_groups=rg_groups,
    )
end

function structured_rg_trace_bounds(B, n_super; link_mode=:rho3)
    link_mode in (:rho3, :rdm8) ||
        error("unsupported structured-RG link mode")
    compressed = compressed_trace_bounds(B, n_super)
    bounds = Dict{String, ExactRational}(
        link_mode == :rho3 ?
            "hybrid_rho3" => compressed["rho3"] :
            "hybrid_rho4" => 1 // 1,
    )
    for depth in 4:n_super
        bounds["hybrid_omega$(depth)"] = compressed["omega$(depth)"]
    end
    bounds
end

function structured_rg_certificate(
    primal_model,
    dual_model,
    structured,
    structured_local_rdm,
    B,
    n_super,
    transcript_path,
    ;
    link_mode=:rho3,
)
    structure = structured_rg_model_structure(
        primal_model,
        dual_model,
        structured,
        structured_local_rdm,
        B,
        n_super,
        link_mode=link_mode,
    )
    transcript = merge(
        structure,
        (raw_values=Float64[value(variable) for variable in all_variables(dual_model)],),
    )
    atomic_serialize(transcript_path, transcript)
    evaluated = evaluate_structured_rg_transcript(
        transcript,
        structure,
        structured_rg_trace_bounds(B, n_super; link_mode=link_mode),
    )
    published = decimal_floor(evaluated.certified; digits=9)
    Dict(
        "status" => "CERTIFIED",
        "method" => link_mode == :rho3 ?
            "exact mixed SOHS/RG replay with an exact six-spin moment link, energy-aware Pauli-word residual, and analytic RG trace corrections" :
            "exact mixed SOHS/RG replay with an exact eight-spin RDM link, dyadic compression, energy-aware Pauli-word residual, and analytic RG trace corrections",
        "physical_relaxation_compatibility_proven" => true,
        "replay_transcript" => artifact_relative_path(transcript_path),
        "replay_transcript_sha256" => file_sha256(transcript_path),
        "solver_objective" => objective_value(dual_model),
        "solver_objective_rational" =>
            string(exact_rational(objective_value(dual_model))),
        "repaired_objective_rational" => string(evaluated.objective),
        "moment_identity_residual_l1_rational" =>
            string(evaluated.moment_correction),
        "maximum_moment_identity_residual_rational" =>
            string(maximum(abs, evaluated.moment_residuals)),
        "rg_slack_correction_rational" =>
            string(evaluated.rg_correction),
        "residual_bound" => Dict(
            "hamiltonian_word" => H11_HAMILTONIAN_WORD,
            "hamiltonian_coefficient_rational" =>
                string(H11_HAMILTONIAN_COEFFICIENT),
            "identity_residual_rational" =>
                string(evaluated.residual_bound.identity_residual),
            "hamiltonian_residual_rational" =>
                string(evaluated.residual_bound.hamiltonian_residual),
            "other_residual_l1_rational" =>
                string(evaluated.residual_bound.other_residual_l1),
            "denominator_rational" =>
                string(evaluated.residual_bound.denominator),
            "legacy_lower_bound_rational" =>
                string(evaluated.legacy_certified),
            "energy_aware_lower_bound_rational" =>
                string(evaluated.residual_bound.energy_aware_lower_bound),
        ),
        "certified_lower_bound" => float_down(evaluated.certified),
        "certified_lower_bound_rational" =>
            string(evaluated.certified),
        "published_certified_lower_bound" => float_down(published),
        "published_certified_lower_bound_rational" =>
            string(published),
        "sohs_groups" => evaluated.sohs_groups,
        "rg_groups" => evaluated.rg_groups,
    )
end

function structured_rg_replay_transcript(certificate, label; expected_structure=nothing)
    certificate["status"] == "CERTIFIED" ||
        error("$label certificate did not pass")
    path = normpath(joinpath(
        BOUNDGSENERGY_ROOT,
        certificate["replay_transcript"],
    ))
    isfile(path) || error("$label transcript does not exist")
    file_sha256(path) == certificate["replay_transcript_sha256"] ||
        error("$label transcript hash is stale")
    transcript = deserialize(path)
    rebuilt = isnothing(expected_structure) ?
        structured_rg_rebuild_model_structure(transcript) :
        expected_structure
    evaluated = evaluate_structured_rg_transcript(
        transcript,
        rebuilt,
        structured_rg_trace_bounds(
            transcript.coarse_grainer,
            transcript.n_super,
            link_mode=transcript.version == 3 ? :rho3 : :rdm8,
        ),
    )
    parse_exact_rational(certificate["solver_objective_rational"]) ==
    exact_rational(certificate["solver_objective"]) ||
        error("$label solver objective is stale")
    parse_exact_rational(certificate["repaired_objective_rational"]) ==
    evaluated.objective ||
        error("$label repaired objective is stale")
    parse_exact_rational(
        certificate["moment_identity_residual_l1_rational"],
    ) == evaluated.moment_correction ||
        error("$label moment correction is stale")
    parse_exact_rational(
        certificate["maximum_moment_identity_residual_rational"],
    ) == maximum(abs, evaluated.moment_residuals) ||
        error("$label maximum moment residual is stale")
    parse_exact_rational(certificate["rg_slack_correction_rational"]) ==
    evaluated.rg_correction ||
        error("$label RG correction is stale")
    residual_bound = certificate["residual_bound"]
    residual_bound["hamiltonian_word"] == H11_HAMILTONIAN_WORD ||
        error("$label Hamiltonian word is stale")
    parse_exact_rational(
        residual_bound["hamiltonian_coefficient_rational"],
    ) == H11_HAMILTONIAN_COEFFICIENT ||
        error("$label Hamiltonian coefficient is stale")
    for (field, value) in (
        "identity_residual_rational" =>
            evaluated.residual_bound.identity_residual,
        "hamiltonian_residual_rational" =>
            evaluated.residual_bound.hamiltonian_residual,
        "other_residual_l1_rational" =>
            evaluated.residual_bound.other_residual_l1,
        "denominator_rational" =>
            evaluated.residual_bound.denominator,
        "legacy_lower_bound_rational" =>
            evaluated.legacy_certified,
        "energy_aware_lower_bound_rational" =>
            evaluated.residual_bound.energy_aware_lower_bound,
    )
        parse_exact_rational(residual_bound[field]) == value ||
            error("$label $field is stale")
    end
    parse_exact_rational(certificate["certified_lower_bound_rational"]) ==
    evaluated.certified ||
        error("$label certified lower bound is stale")
    certificate["certified_lower_bound"] ==
    float_down(evaluated.certified) ||
        error("$label certified lower-bound summary is stale")
    certificate["physical_relaxation_compatibility_proven"] ||
        error("$label physical compatibility summary is stale")
    certificate["sohs_groups"] == evaluated.sohs_groups ||
        error("$label SOHS repair summary is stale")
    certificate["rg_groups"] == evaluated.rg_groups ||
        error("$label RG repair summary is stale")
    published = decimal_floor(evaluated.certified; digits=9)
    parse_exact_rational(
        certificate["published_certified_lower_bound_rational"],
    ) == published ||
        error("$label published lower bound is stale")
    certificate["published_certified_lower_bound"] ==
    float_down(published) ||
        error("$label published lower-bound summary is stale")
    transcript, evaluated, rebuilt
end

function structured_rg_replay_baseline(certificate, label; expected_structure=nothing)
    certificate["status"] == "CERTIFIED" ||
        error("$label baseline certificate did not pass")
    path = normpath(joinpath(
        BOUNDGSENERGY_ROOT,
        certificate["replay_transcript"],
    ))
    isfile(path) || error("$label baseline transcript does not exist")
    file_sha256(path) == certificate["replay_transcript_sha256"] ||
        error("$label baseline transcript hash is stale")
    transcript = deserialize(path)
    rebuilt = isnothing(expected_structure) ?
        h11_rebuild_model_structure(h11_configuration(false)) :
        expected_structure
    legacy_certified, objective, residual_l1, residuals, groups =
        h11_evaluate_transcript(transcript, rebuilt)
    sharpened = h11_energy_aware_residual_bound(
        objective,
        residuals,
        transcript.canonical_words,
    )
    certified = sharpened.certified_lower_bound
    sharpened.legacy_lower_bound == legacy_certified ||
        error("$label baseline legacy residual bound is stale")
    parse_exact_rational(certificate["raw_objective_rational"]) ==
    objective ||
        error("$label baseline objective is stale")
    certificate["raw_objective"] == float_down(objective) ||
        error("$label baseline objective summary is stale")
    parse_exact_rational(certificate["identity_residual_l1_rational"]) ==
    residual_l1 ||
        error("$label baseline residual correction is stale")
    certificate["identity_residual_l1"] == float_up(residual_l1) ||
        error("$label baseline residual summary is stale")
    residual_bound = certificate["residual_bound"]
    residual_bound["hamiltonian_word"] == H11_HAMILTONIAN_WORD ||
        error("$label baseline Hamiltonian word is stale")
    parse_exact_rational(
        residual_bound["hamiltonian_coefficient_rational"],
    ) == H11_HAMILTONIAN_COEFFICIENT ||
        error("$label baseline Hamiltonian coefficient is stale")
    for (field, value) in (
        "identity_residual_rational" => sharpened.identity_residual,
        "hamiltonian_residual_rational" =>
            sharpened.hamiltonian_residual,
        "other_residual_l1_rational" => sharpened.other_residual_l1,
        "denominator_rational" => sharpened.denominator,
        "legacy_lower_bound_rational" => sharpened.legacy_lower_bound,
        "energy_aware_lower_bound_rational" =>
            sharpened.energy_aware_lower_bound,
    )
        parse_exact_rational(residual_bound[field]) == value ||
            error("$label baseline $field is stale")
    end
    parse_exact_rational(certificate["certified_lower_bound_rational"]) ==
    certified ||
        error("$label baseline lower bound is stale")
    certificate["certified_lower_bound"] == float_down(certified) ||
        error("$label baseline lower-bound summary is stale")
    published = decimal_floor(certified; digits=9)
    parse_exact_rational(
        certificate["published_certified_lower_bound_rational"],
    ) == published ||
        error("$label baseline published lower bound is stale")
    certificate["published_certified_lower_bound"] ==
    float_down(published) ||
        error("$label baseline published lower-bound summary is stale")
    certificate["maximum_identity_residual"] ==
    float_up(maximum(abs, residuals)) ||
        error("$label baseline maximum residual is stale")
    certificate["groups"] == groups ||
        error("$label baseline PSD repair summary is stale")
    certified
end

function structured_rg_validate_summary(
    record,
    label,
    upper,
    baseline_certified,
    baseline_objective,
    hybrid_certified,
    hybrid_objective,
)
    baseline_energy = Float64(baseline_objective)
    hybrid_energy = Float64(hybrid_objective)
    upper_energy = Float64(upper)
    exact_rational(record["baseline"]["energy"]) == baseline_objective ||
        error("$label baseline solver objective is stale")
    exact_rational(record["hybrid"]["energy"]) == hybrid_objective ||
        error("$label hybrid solver objective is stale")

    comparison = record["comparison"]
    comparison["raw_improvement"] == hybrid_energy - baseline_energy ||
        error("$label raw comparison is stale")
    comparison["hybrid_raw_gap_to_upper"] ==
    upper_energy - hybrid_energy ||
        error("$label raw gap is stale")
    parse_exact_rational(
        comparison["certified_improvement_rational"],
    ) == hybrid_certified - baseline_certified ||
        error("$label certified comparison is stale")
    parse_exact_rational(
        comparison["hybrid_certified_gap_rational"],
    ) == upper - hybrid_certified ||
        error("$label certified gap is stale")

    gates = record["gates"]
    monotonicity_passed =
        hybrid_energy + STRUCTURED_RG_ENERGY_TOLERANCE >=
        baseline_energy
    finite_ceiling_passed =
        hybrid_energy <= upper_energy + STRUCTURED_RG_ENERGY_TOLERANCE
    gates["rg_monotonicity_passed"] ==
    monotonicity_passed ||
        error("$label RG monotonicity gate is stale")
    gates["finite_ceiling_passed"] ==
    finite_ceiling_passed ||
        error("$label finite ceiling gate is stale")
    gates["strict_improvement_passed"] ==
    (hybrid_certified > baseline_certified) ||
        error("$label strict improvement gate is stale")
    gates["target_passed"] ==
    (upper - hybrid_certified <= 1 // 100000) ||
        error("$label target gate is stale")
    true
end

function validated_structured_rg_result(
    record,
    label;
    expected_hybrid_structure=nothing,
    expected_baseline_structure=nothing,
)
    meta = record["meta"]
    meta["stage"] == "structured_rg_hybrid" ||
        error("$label stage is stale")
    meta["julia_version"] == string(PINNED_JULIA_VERSION) ||
        error("$label Julia version is stale")
    meta["system_size"] == H10_N ||
        error("$label system size is stale")
    meta["boundary"] == "periodic" ||
        error("$label boundary is stale")
    validate_structured_rg_support(meta["n_super"]; minimum_depth=5)
    meta["support_sites"] == 2 * meta["n_super"] ||
        error("$label support size is stale")
    meta["support_sites"] <= meta["system_size"] ||
        error("$label support exceeds the finite system")
    meta["bond_dimension"] == 4 ||
        error("$label bond dimension is stale")
    meta["dyadic_bits"] == RDM8_DYADIC_BITS ||
        error("$label dyadic grid is stale")
    link_mode = Symbol(get(meta, "link_mode", "rho3"))
    link_mode in (:rho3, :rdm8) ||
        error("$label link mode is stale")
    record["baseline"]["status"] == "OPTIMAL" ||
        error("$label baseline did not terminate OPTIMAL")
    record["hybrid"]["status"] == "OPTIMAL" ||
        error("$label hybrid did not terminate OPTIMAL")
    record["baseline"]["primal_status"] == "FEASIBLE_POINT" ||
        error("$label baseline primal status is stale")
    record["hybrid"]["primal_status"] == "FEASIBLE_POINT" ||
        error("$label hybrid primal status is stale")

    upper_certificate = validated_finite_upper_certificate(
        record["finite_energy_upper_certificate"],
        "$label finite upper bound",
    )
    upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
    baseline_certified = structured_rg_replay_baseline(
        record["baseline_certificate"],
        label;
        expected_structure=expected_baseline_structure,
    )
    transcript, evaluated, _ = structured_rg_replay_transcript(
        record["hybrid_certificate"],
        label;
        expected_structure=expected_hybrid_structure,
    )
    baseline_objective = parse_exact_rational(
        record["baseline_certificate"]["raw_objective_rational"],
    )
    hybrid_objective = parse_exact_rational(
        record["hybrid_certificate"]["solver_objective_rational"],
    )
    transcript.system_size == meta["system_size"] ||
        error("$label transcript system size is stale")
    transcript.n_super == meta["n_super"] ||
        error("$label transcript hierarchy depth is stale")
    transcript.version == (link_mode == :rho3 ? 3 : 4) ||
        error("$label transcript link mode is stale")
    for field in (
        "local_link_scale",
        "omega4_link_scale",
        "recursion_scale",
    )
        exact_rational(meta[field]) ==
        exact_rational(getproperty(transcript, Symbol(field))) ||
            error("$label $field is stale")
    end
    size(transcript.coarse_grainer, 1) == meta["bond_dimension"] ||
        error("$label transcript bond dimension is stale")
    record["coarse_grainer"]["local_rdm_exactness_audit"] == Dict(
        string(field) =>
            getproperty(transcript.local_rdm_exactness_audit, field)
        for field in propertynames(transcript.local_rdm_exactness_audit)
    ) || error("$label local-RDM exactness audit is stale")
    baseline_certified <= upper ||
        error("$label baseline crossed the finite upper bound")
    evaluated.certified <= upper ||
        error("$label hybrid crossed the finite upper bound")
    structured_rg_validate_summary(
        record,
        label,
        upper,
        baseline_certified,
        baseline_objective,
        evaluated.certified,
        hybrid_objective,
    )
    record["gates"]["both_optimal"] ||
        error("$label optimality gate is stale")
    record["gates"]["strict_hybrid_certificate_passed"] ||
        error("$label strict certificate gate is stale")
    record["gates"]["local_rdm_exactness_passed"] ||
        error("$label local-RDM exactness gate is stale")
    isapprox(record["rho3_checks"]["trace"], 1; atol=1e-7) ||
        error("$label local-RDM trace is stale")
    record["rho3_checks"]["minimum_eigenvalue"] >= -1e-7 ||
        error("$label local RDM is not numerically PSD")
    record["rho3_checks"]["translation_marginal_residual"] <= 1e-7 ||
        error("$label rho3 translation residual is stale")
    record["rho3_checks"]["objective_energy_residual"] <= 1e-7 ||
        error("$label rho3 energy convention is stale")
    record["hybrid_certificate"]
end
