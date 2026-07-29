using Test

include(joinpath(@__DIR__, "h11_full_sohs_certificate.jl"))

function parse_h11_test_args(args)
    isempty(args) && return nothing
    length(args) == 1 && startswith(only(args), "--result=") ||
        error("expected at most one --result=path")
    path = split(only(args), "=", limit=2)[2]
    normpath(isabspath(path) ? path : joinpath(BOUNDGSENERGY_ROOT, path))
end

@testset "H11 full SOHS model" begin
    @test VERSION == PINNED_JULIA_VERSION
    run_name, smoke, threads, tolerance, upper_result = parse_h11_args([
        "--run-name=test-h11",
        "--mosek-threads=8",
    ])
    @test run_name == "test-h11"
    @test !smoke
    @test threads == 8
    @test tolerance == 1e-8
    @test isnothing(upper_result)

    expected_structure = h11_rebuild_model_structure(h11_configuration(false))
    @test expected_structure.version == 3
    @test expected_structure.objective_sense == "MAX_SENSE"
    @test expected_structure.variable_count == 59_413
    @test length(expected_structure.equalities) == length(expected_structure.canonical_words) == 3_432
    @test length(expected_structure.psd_blocks) == 49
    @test maximum(block.dimension for block in expected_structure.psd_blocks) == 126

    result_path = parse_h11_test_args(ARGS)
    if !isnothing(result_path)
        isfile(result_path) || error("H11 result does not exist: $result_path")
        result = JSON.parsefile(result_path)
        certificate = validated_h11_certificate(
            result,
            "generated H11";
            expected_structure,
        )
        lower = parse_exact_rational(certificate["certified_lower_bound_rational"])
        upper = parse_exact_rational(
            result["finite_energy_upper_certificate"]["upper_bound_rational"],
        )

        @test result["solve"]["status"] == "OPTIMAL"
        @test result["solve"]["primal_status"] == "FEASIBLE_POINT"
        @test certificate["status"] == "CERTIFIED"
        @test certificate["qmbcertify_source_tree_sha256"] == H11_QMB_SOURCE_TREE_SHA256
        @test lower <= upper
        @test upper - lower <= 1 // 100000
        @test result["gates"]["sandwich_passed"]
        @test result["gates"]["target_passed"]
        @test all(
            group["diagonal_shift_rational"] == "0//1"
            for group in values(certificate["groups"])
        )

        transcript_path = joinpath(BOUNDGSENERGY_ROOT, certificate["replay_transcript"])
        transcript = deserialize(transcript_path)
        @test transcript.equality_scale == H11_EQUALITY_SCALE
        @test transcript.objective_sense == "MAX_SENSE"
        @test transcript.variable_count == length(transcript.raw_values)
        @test length(transcript.equalities) == length(transcript.canonical_words) == 3_432
        @test length(transcript.psd_blocks) == 49

        bad_hash = deepcopy(result)
        bad_hash["certificate"]["replay_transcript_sha256"] = repeat("0", 64)
        @test_throws ErrorException validated_h11_certificate(bad_hash, "bad hash")

        bad_status = deepcopy(result)
        bad_status["solve"]["status"] = "SLOW_PROGRESS"
        @test_throws ErrorException validated_h11_certificate(bad_status, "bad status")

        bad_primal_status = deepcopy(result)
        bad_primal_status["solve"]["primal_status"] = "UNKNOWN_RESULT_STATUS"
        @test_throws ErrorException validated_h11_certificate(bad_primal_status, "bad primal status")

        mktempdir() do directory
            raw_values = copy(transcript.raw_values)
            raw_values[1] += 1.0
            tampered_path = joinpath(directory, "tampered.jls")
            atomic_serialize(tampered_path, merge(transcript, (raw_values=raw_values,)))
            bad_values = deepcopy(result)
            bad_values["certificate"]["replay_transcript"] = tampered_path
            bad_values["certificate"]["replay_transcript_sha256"] = file_sha256(tampered_path)
            @test_throws ErrorException validated_h11_certificate(
                bad_values,
                "bad values";
                expected_structure,
            )

            equalities = copy(transcript.equalities)
            equalities[1] = merge(equalities[1], (rhs=equalities[1].rhs + 1.0,))
            altered_model_path = joinpath(directory, "altered-model.jls")
            atomic_serialize(altered_model_path, merge(transcript, (equalities=equalities,)))
            bad_model = deepcopy(result)
            bad_model["certificate"]["replay_transcript"] = altered_model_path
            bad_model["certificate"]["replay_transcript_sha256"] = file_sha256(altered_model_path)
            @test_throws ErrorException validated_h11_certificate(
                bad_model,
                "altered model";
                expected_structure,
            )

            altered_sense_path = joinpath(directory, "altered-sense.jls")
            atomic_serialize(
                altered_sense_path,
                merge(transcript, (objective_sense="MIN_SENSE",)),
            )
            bad_sense = deepcopy(result)
            bad_sense["certificate"]["replay_transcript"] = altered_sense_path
            bad_sense["certificate"]["replay_transcript_sha256"] = file_sha256(altered_sense_path)
            @test_throws ErrorException validated_h11_certificate(
                bad_sense,
                "altered objective sense";
                expected_structure,
            )
        end
    end
end
