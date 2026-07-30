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
