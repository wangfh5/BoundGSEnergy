using Test

include(joinpath(@__DIR__, "structured_rg_certificate.jl"))

function parse_structured_rg_test_args(args)
    isempty(args) && return nothing
    length(args) == 1 && startswith(only(args), "--result=") ||
        error("expected at most one --result=path")
    path = split(only(args), "=", limit=2)[2]
    normpath(isabspath(path) ? path : joinpath(BOUNDGSENERGY_ROOT, path))
end

@test STRUCTURED_RG_RDM9_DIMENSIONS == [1, 9, 36, 84, 126]
@test sum(STRUCTURED_RG_RDM9_DIMENSIONS) == 256
@test STRUCTURED_RG_EQUALITY_SCALE == 1.0
@test STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE == 1.0
@test STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE == 16.0
@test STRUCTURED_RG_RDM8_RECURSION_SCALE == 1.0
@test STRUCTURED_RG_VIRTUAL_CHARGES == [-1, 1, -1, 1]
@test STRUCTURED_RG_SUPER_CHARGES == [-2, 0, 0, 2]
@test validate_structured_rg_support(10; minimum_depth=5) == 10
@test_throws ErrorException validate_structured_rg_support(4; minimum_depth=5)
@test_throws ErrorException validate_structured_rg_support(11; minimum_depth=5)
stored_basis = sort(vcat(STRUCTURED_RG_RDM9_BLOCKS...))
flipped_basis = sort(513 .- stored_basis)
@test isempty(intersect(stored_basis, flipped_basis))
@test sort(vcat(stored_basis, flipped_basis)) == collect(1:512)

function local_canonical_words(site_count)
    words = Set{Tuple}()
    labels = zeros(Int, site_count)
    for code in 0:(4^site_count - 1)
        remainder = code
        for site in site_count:-1:1
            labels[site] = remainder % 4
            remainder ÷= 4
        end
        word = UInt16[
            3 * (site - 1) + labels[site]
            for site in 1:site_count
            if labels[site] != 0
        ]
        reduced, coefficient =
            QMBCertify.reduce!(word; L=H10_N, lattice="chain")
        iszero(coefficient) || push!(words, Tuple(reduced))
    end
    [UInt16[word...] for word in sort!(collect(words); by=string)]
end

toy_model = Model()
canonical_words = local_canonical_words(6)
scaled_moments =
    @variable(toy_model, [1:length(canonical_words)], base_name="toy_moment")
physical_moments = -H11_EQUALITY_SCALE .* scaled_moments
rho3 = structured_local_rdm(
    physical_moments,
    canonical_words,
    6,
)
@test size(rho3) == (64, 64)
@test all(rho3[row, column] == rho3[column, row] for row in 1:64, column in 1:64)
local_rdm_audit =
    structured_local_rdm_exactness_audit(
        rho3,
        physical_moments,
        canonical_words,
        6,
    )
@test local_rdm_audit.status == "PROVEN"
@test local_rdm_audit.exact_rational_reconstruction_matches
@test local_rdm_audit.coefficient_assembly_exact
@test structured_rg_validate_local_rdm_audit(
    (version=3, local_rdm_exactness_audit=local_rdm_audit),
)
bad_local_rdm_audit =
    merge(local_rdm_audit, (coefficient_assembly_exact=false,))
@test_throws ErrorException structured_rg_validate_local_rdm_audit(
    (version=3, local_rdm_exactness_audit=bad_local_rdm_audit),
)

virtual_charges = STRUCTURED_RG_VIRTUAL_CHARGES
super_charges = STRUCTURED_RG_SUPER_CHARGES
B = zeros(4, 4, 4)
for left in 1:4, site in 1:4, right in 1:4
    virtual_charges[left] + super_charges[site] == virtual_charges[right] ||
        continue
    B[left, site, right] = 1 / 8
end
audit = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)
@test audit["physical_relaxation_compatibility_proven"]
@test structured_rg_omega_trace_scales(B, 10, :none) == ones(7)
@test structured_rg_rdm8_link_scales(:none) == (
    local_link=1.0,
    omega4=16.0,
    recursion=1.0,
)
@test structured_rg_rdm8_link_scales(:dyadic_trace) == (
    local_link=1.0,
    omega4=32.0,
    recursion=8.0,
)
@test_throws ErrorException structured_rg_rdm8_link_scales(
    :unsupported,
)
normalization_B = 4 .* B
normalization_audit =
    dyadic_relaxation_audit(normalization_B, RDM8_DYADIC_BITS)
@test normalization_audit["physical_relaxation_compatibility_proven"]
normalized_scales =
    structured_rg_omega_trace_scales(
        normalization_B,
        10,
        :dyadic_trace,
    )
@test any(>(1), normalized_scales)
@test any(
    normalized_scales[index] / normalized_scales[index + 1] < 1
    for index in 1:(length(normalized_scales) - 1)
)
trace_bounds = compressed_trace_bounds(normalization_B, 10)
raw_rg_trace_bounds =
    structured_rg_trace_bounds(
        normalization_B,
        10;
        link_mode=:rdm8,
    )
normalized_rg_trace_bounds = structured_rg_trace_bounds(
    normalization_B,
    10;
    link_mode=:rdm8,
    omega_normalization=:dyadic_trace,
)
for (index, depth) in enumerate(4:10)
    scale = exact_rational(normalized_scales[index])
    bound = trace_bounds["omega$(depth)"]
    group = "hybrid_omega$(depth)"
    @test ispow2(Int(normalized_scales[index]))
    @test scale >= bound
    @test scale == 1 || scale / 2 < bound
    @test normalized_rg_trace_bounds[group] ==
        raw_rg_trace_bounds[group] / scale
    @test normalized_rg_trace_bounds[group] <= 1
    if scale != 1
        @test normalized_rg_trace_bounds[group] > 1 // 2
    end
end
@test_throws ErrorException structured_rg_omega_trace_scales(
    B,
    10,
    :unsupported,
)
@test_throws ErrorException structured_rg_trace_bounds(
    B,
    10;
    link_mode=:rho3,
    omega_normalization=:dyadic_trace,
)

@testset "safe replay residual polishing" begin
    polishing_transcript = (
        variable_count=9,
        canonical_words=[
            Int[],
            H11_HAMILTONIAN_WORD,
            [2],
            [3],
        ],
        objective=(
            constant=0.0,
            terms=[(1, 1.0)],
        ),
        moment_equalities=[
            (
                constant=0.0,
                rhs=0.0,
                terms=[(2, 1.0)],
            ),
            (
                constant=0.0,
                rhs=0.0,
                terms=[(3, 1.0)],
            ),
            (
                constant=0.0,
                rhs=0.0,
                terms=[
                    (4, 1.0),
                    (5, 1.0),
                    (8, 1.0),
                ],
            ),
            (
                constant=0.0,
                rhs=0.0,
                terms=[(9, 1.0)],
            ),
        ],
        psd_variable_blocks=[
            (
                group="toy_sohs",
                dimension=1,
                variable_indices=reshape([6], 1, 1),
            ),
        ],
        omega_slacks=[
            (
                group="toy_omega",
                dimension=1,
                upper_entries=[
                    (
                        constant=0.0,
                        rhs=0.0,
                        terms=[(7, 1.0)],
                    ),
                ],
            ),
        ],
    )
    polishing_audit =
        structured_rg_polishing_variable_audit(
            polishing_transcript,
        )
    @test polishing_audit.covered_variable_count == 9
    @test polishing_audit.safe_variable_indices == [4, 5, 8, 9]
    uncovered_transcript = merge(
        polishing_transcript,
        (variable_count=10,),
    )
    @test_throws ErrorException structured_rg_polishing_variable_audit(
        uncovered_transcript,
    )

    matrix = sparse([1.0 0.0 1.0; 0.0 1.0 1.0])
    numerical = structured_rg_sparse_minimum_norm(
        matrix,
        [1.0, 2.0],
    )
    @test matrix * numerical.solution ≈ [1.0, 2.0]
    @test norm(numerical.solution) ≈ sqrt(2)
    @test numerical.solution ≈ [0.0, 1.0, 1.0] atol = 1e-12

    original = (
        certified=1 // 4,
        objective=1 // 2,
        sohs_groups=Dict("sohs" => 1),
        rg_groups=Dict("omega" => 2),
        rg_correction=-1 // 16,
        moment_residuals=ExactRational[
            0,
            0,
            1 // 8,
        ],
    )
    improved = merge(
        original,
        (
            certified=3 // 8,
            moment_residuals=ExactRational[
                0,
                0,
                0,
            ],
        ),
    )
    @test structured_rg_accept_polished_replay(
        original,
        improved,
    )
    @test !structured_rg_accept_polished_replay(
        original,
        merge(improved, (certified=original.certified,)),
    )
    @test_throws ErrorException structured_rg_accept_polished_replay(
        original,
        merge(
            improved,
            (objective=original.objective + 1 // 16,),
        ),
    )
end

@testset "structured-RG constraint partition" begin
    partition_model = Model()
    @variable(partition_model, partition_x[1:3])
    equality_one =
        @constraint(partition_model, 2 * partition_x[1] == 0)
    equality_two =
        @constraint(partition_model, 3 * partition_x[2] == 0)
    sohs_constraint = @constraint(
        partition_model,
        [partition_x[3]] in
        MOI.PositiveSemidefiniteConeTriangle(1),
    )
    omega_constraint = @constraint(
        partition_model,
        [partition_x[1] + partition_x[3]] in
        MOI.PositiveSemidefiniteConeTriangle(1),
    )
    moment_occurrences = Any[
        (equality_one, 1),
        (equality_two, 1),
    ]
    partition_audit = structured_rg_constraint_partition_audit(
        partition_model,
        moment_occurrences,
        Any[sohs_constraint],
        Any[omega_constraint],
    )
    @test partition_audit.equality_component_count == 2
    @test partition_audit.sohs_constraint_count == 1
    @test partition_audit.omega_constraint_count == 1
    @test_throws ErrorException structured_rg_constraint_partition_audit(
        partition_model,
        moment_occurrences[1:1],
        Any[sohs_constraint],
        Any[omega_constraint],
    )
    @test_throws ErrorException structured_rg_constraint_partition_audit(
        partition_model,
        moment_occurrences,
        Any[sohs_constraint, sohs_constraint],
        Any[omega_constraint],
    )
    @test_throws ErrorException structured_rg_constraint_partition_audit(
        partition_model,
        moment_occurrences,
        Any[sohs_constraint],
        Any[],
    )
    @constraint(partition_model, partition_x[1] >= -1)
    @test_throws ErrorException structured_rg_constraint_partition_audit(
        partition_model,
        moment_occurrences,
        Any[sohs_constraint],
        Any[omega_constraint],
    )
end

@testset "dyadic omega variable substitution" begin
    unscaled_model = Model()
    unscaled = add_rdm8_rg_hierarchy!(
        unscaled_model,
        zeros(256, 256),
        normalization_B,
        5,
        virtual_charges,
        super_charges,
        omega_normalization=:none,
    )
    normalized_model = Model()
    normalized = add_rdm8_rg_hierarchy!(
        normalized_model,
        zeros(256, 256),
        normalization_B,
        5,
        virtual_charges,
        super_charges,
        omega_normalization=:dyadic_trace,
    )
    @test any(>(1), normalized.omega_trace_scales)
    @test normalized.omega_trace_scales[1] /
        normalized.omega_trace_scales[2] < 1
    @test (
        normalized.local_link_scale,
        normalized.omega4_link_scale,
        normalized.recursion_scale,
    ) == (1.0, 32.0, 8.0)

    unscaled_psd_blocks = all_constraints(
        unscaled_model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    normalized_psd_blocks = all_constraints(
        normalized_model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    rho4_block_count =
        length(charge_blocks(rho4_charges(super_charges)))
    omega_block_count = length(charge_blocks(
        omega_charges(super_charges, virtual_charges),
    ))
    @test length(unscaled_psd_blocks) ==
        rho4_block_count +
        length(unscaled.omega) * omega_block_count
    @test length(normalized_psd_blocks) ==
        length(unscaled_psd_blocks)
    variable_map = Dict{VariableRef, VariableRef}()
    omega_scales = Dict{VariableRef, Float64}()
    for (block_index, (unscaled_block, normalized_block)) in
        enumerate(zip(unscaled_psd_blocks, normalized_psd_blocks))
        unscaled_variables = constraint_object(unscaled_block).func
        normalized_variables = constraint_object(normalized_block).func
        length(unscaled_variables) == length(normalized_variables) ||
            error("normalization changed a PSD block dimension")
        omega_index = cld(
            max(0, block_index - rho4_block_count),
            omega_block_count,
        )
        scale = omega_index == 0 ?
            1.0 :
            normalized.omega_trace_scales[omega_index]
        for (unscaled_variable, normalized_variable) in
            zip(unscaled_variables, normalized_variables)
            variable_map[normalized_variable] = unscaled_variable
            omega_scales[normalized_variable] = scale
        end
    end
    @test Set(keys(variable_map)) ==
        Set(all_variables(normalized_model))
    @test Set(values(variable_map)) ==
        Set(all_variables(unscaled_model))

    unscaled_values = Dict(
        variable => (mod(index, 17) - 8) / 16
        for (index, variable) in
            enumerate(all_variables(unscaled_model))
    )
    normalized_values = Dict(
        variable =>
            unscaled_values[variable_map[variable]] /
            omega_scales[variable]
        for variable in all_variables(normalized_model)
    )
    unscaled_equalities = all_constraints(
        unscaled_model,
        AffExpr,
        MOI.EqualTo{Float64},
    )
    normalized_equalities = all_constraints(
        normalized_model,
        AffExpr,
        MOI.EqualTo{Float64},
    )
    @test length(unscaled_equalities) ==
        length(normalized_equalities)
    evaluate_equality(constraint, values) = JuMP.value(
        variable -> values[variable],
        constraint_object(constraint).func,
    )
    unscaled_residuals = [
        evaluate_equality(constraint, unscaled_values)
        for constraint in unscaled_equalities
    ]
    normalized_residuals = [
        evaluate_equality(constraint, normalized_values)
        for constraint in normalized_equalities
    ]
    local_link_count = unscaled.rho4_link_equality_count
    omega4_link_count = unscaled.omega4_link_equality_count
    omega4_first = local_link_count + 1
    omega4_last = local_link_count + omega4_link_count
    recursion_first = omega4_last + 1
    @test normalized_residuals[1:local_link_count] ==
        unscaled_residuals[1:local_link_count]
    @test normalized_residuals[omega4_first:omega4_last] ==
        unscaled_residuals[omega4_first:omega4_last] .*
        (
            normalized.omega4_link_scale /
            unscaled.omega4_link_scale /
            normalized.omega_trace_scales[1]
        )
    @test recursion_first <= length(unscaled_residuals)
    @test normalized_residuals[recursion_first:end] ==
        unscaled_residuals[recursion_first:end] .*
        (
            normalized.recursion_scale /
            unscaled.recursion_scale /
            normalized.omega_trace_scales[2]
        )
end

hierarchy = add_rg_hierarchy!(
    toy_model,
    rho3,
    B,
    4,
    virtual_charges,
    super_charges,
)
@test size(hierarchy.rho3) == (64, 64)
@test length(hierarchy.omega) == 1
@test size(only(hierarchy.omega)) == (256, 256)
@test hierarchy.rho3_link_equality_count > 0
@test_throws ErrorException add_rg_hierarchy!(
    Model(),
    rho3,
    B,
    11,
    virtual_charges,
    super_charges,
)

function test_structured_rg_summary(
    record,
    upper,
    baseline_certified,
    baseline_objective,
    hybrid_certified,
    hybrid_objective,
)
    validate_summary = (candidate, label) -> structured_rg_validate_summary(
        candidate,
        label,
        upper,
        baseline_certified,
        baseline_objective,
        hybrid_certified,
        hybrid_objective,
    )
    @test validate_summary(record, "valid summary")

    bad_baseline_energy = deepcopy(record)
    bad_baseline_energy["baseline"]["energy"] += 1e-6
    @test_throws ErrorException validate_summary(
        bad_baseline_energy,
        "bad baseline energy",
    )

    bad_raw_improvement = deepcopy(record)
    bad_raw_improvement["comparison"]["raw_improvement"] += 1e-6
    @test_throws ErrorException validate_summary(
        bad_raw_improvement,
        "bad raw improvement",
    )

    bad_raw_gap = deepcopy(record)
    bad_raw_gap["comparison"]["hybrid_raw_gap_to_upper"] += 1e-6
    @test_throws ErrorException validate_summary(
        bad_raw_gap,
        "bad raw gap",
    )

    for gate in (
        "rg_monotonicity_passed",
        "finite_ceiling_passed",
        "strict_improvement_passed",
        "target_passed",
    )
        bad_gate = deepcopy(record)
        bad_gate["gates"][gate] = !bad_gate["gates"][gate]
        @test_throws ErrorException validate_summary(
            bad_gate,
            "bad $gate",
        )
    end
end

summary_upper = exact_rational(-0.44)
summary_baseline_energy = -0.440004
summary_hybrid_energy = -0.440003
summary_baseline_objective = exact_rational(summary_baseline_energy)
summary_hybrid_objective = exact_rational(summary_hybrid_energy)
summary_baseline_certified = exact_rational(-0.440006)
summary_hybrid_certified = exact_rational(-0.440005)
summary_record = Dict(
    "baseline" => Dict("energy" => summary_baseline_energy),
    "hybrid" => Dict("energy" => summary_hybrid_energy),
    "comparison" => Dict(
        "raw_improvement" =>
            summary_hybrid_energy - summary_baseline_energy,
        "hybrid_raw_gap_to_upper" =>
            Float64(summary_upper) - summary_hybrid_energy,
        "certified_improvement_rational" => string(
            summary_hybrid_certified - summary_baseline_certified,
        ),
        "hybrid_certified_gap_rational" => string(
            summary_upper - summary_hybrid_certified,
        ),
    ),
    "gates" => Dict(
        "rg_monotonicity_passed" => true,
        "finite_ceiling_passed" => true,
        "strict_improvement_passed" => true,
        "target_passed" => true,
    ),
)
test_structured_rg_summary(
    summary_record,
    summary_upper,
    summary_baseline_certified,
    summary_baseline_objective,
    summary_hybrid_certified,
    summary_hybrid_objective,
)

result_path = parse_structured_rg_test_args(ARGS)
if !isnothing(result_path)
    isfile(result_path) ||
        error("structured-RG result does not exist: $result_path")
    result = JSON.parsefile(result_path)
    certificate = result["hybrid_certificate"]
    transcript_path = normpath(joinpath(
        BOUNDGSENERGY_ROOT,
        certificate["replay_transcript"],
    ))
    transcript = deserialize(transcript_path)
    expected_hybrid = structured_rg_rebuild_model_structure(transcript)
    expected_baseline =
        h11_rebuild_model_structure(h11_configuration(false))
    validated_structured_rg_result(
        result,
        "generated structured-RG result";
        expected_hybrid_structure=expected_hybrid,
        expected_baseline_structure=expected_baseline,
    )

    @test result["baseline"]["status"] == "OPTIMAL"
    @test result["hybrid"]["status"] == "OPTIMAL"
    @test result["gates"]["strict_hybrid_certificate_passed"]
    @test result["gates"]["local_rdm_exactness_passed"]

    upper = parse_exact_rational(
        result["finite_energy_upper_certificate"]["upper_bound_rational"],
    )
    baseline_certified = parse_exact_rational(
        result["baseline_certificate"]["certified_lower_bound_rational"],
    )
    baseline_objective = parse_exact_rational(
        result["baseline_certificate"]["raw_objective_rational"],
    )
    hybrid_certified = parse_exact_rational(
        result["hybrid_certificate"]["certified_lower_bound_rational"],
    )
    hybrid_objective = parse_exact_rational(
        result["hybrid_certificate"]["solver_objective_rational"],
    )
    test_structured_rg_summary(
        result,
        upper,
        baseline_certified,
        baseline_objective,
        hybrid_certified,
        hybrid_objective,
    )

    bad_status = deepcopy(result)
    bad_status["hybrid"]["status"] = "SLOW_PROGRESS"
    @test_throws ErrorException validated_structured_rg_result(
        bad_status,
        "bad status";
        expected_hybrid_structure=expected_hybrid,
        expected_baseline_structure=expected_baseline,
    )

    bad_hash = deepcopy(result)
    bad_hash["hybrid_certificate"]["replay_transcript_sha256"] =
        repeat("0", 64)
    @test_throws ErrorException validated_structured_rg_result(
        bad_hash,
        "bad hash";
        expected_hybrid_structure=expected_hybrid,
        expected_baseline_structure=expected_baseline,
    )

    bad_audit = deepcopy(result)
    bad_audit["coarse_grainer"]["local_rdm_exactness_audit"][
        "coefficient_assembly_exact"
    ] = false
    @test_throws ErrorException validated_structured_rg_result(
        bad_audit,
        "bad coefficient audit";
        expected_hybrid_structure=expected_hybrid,
        expected_baseline_structure=expected_baseline,
    )

    mktempdir() do directory
        altered_scale_path = joinpath(directory, "altered-scale.jls")
        atomic_serialize(
            altered_scale_path,
            merge(transcript, (omega4_link_scale=8.0,)),
        )
        altered_scale_result = deepcopy(result)
        altered_scale_result["hybrid_certificate"]["replay_transcript"] =
            altered_scale_path
        altered_scale_result["hybrid_certificate"]["replay_transcript_sha256"] =
            file_sha256(altered_scale_path)
        @test_throws ErrorException validated_structured_rg_result(
            altered_scale_result,
            "altered link scale";
            expected_hybrid_structure=expected_hybrid,
            expected_baseline_structure=expected_baseline,
        )

        altered_slacks = copy(transcript.omega_slacks)
        first_slack = altered_slacks[1]
        entries = copy(first_slack.upper_entries)
        entries[1] = merge(entries[1], (constant=entries[1].constant + 1.0,))
        altered_slacks[1] = merge(first_slack, (upper_entries=entries,))
        altered_path = joinpath(directory, "altered-model.jls")
        atomic_serialize(
            altered_path,
            merge(transcript, (omega_slacks=altered_slacks,)),
        )
        altered_result = deepcopy(result)
        altered_result["hybrid_certificate"]["replay_transcript"] =
            altered_path
        altered_result["hybrid_certificate"]["replay_transcript_sha256"] =
            file_sha256(altered_path)
        @test_throws ErrorException validated_structured_rg_result(
            altered_result,
            "altered model";
            expected_hybrid_structure=expected_hybrid,
            expected_baseline_structure=expected_baseline,
        )
    end
end

println("structured-RG formulation tests passed")
