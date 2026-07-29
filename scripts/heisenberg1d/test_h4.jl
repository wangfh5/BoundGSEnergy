using Test

include("h4_finite_n200_d5.jl")

function parse_h4_test_args(args)
    paths = Dict{String, String}()
    for arg in args
        startswith(arg, "--") && occursin("=", arg) || error("expected --name=path, got $arg")
        name, path = split(arg[3:end], "=", limit=2)
        paths[name] = h4_input_path(path)
    end
    paths
end

@test expected_previous_depth(50) === nothing
@test expected_previous_depth(75) == 50
@test expected_previous_depth(100) == 75
@test bethe_ground_reference(200)["energy_per_site"] ≈ -0.44316776626120796 atol=2e-12

artifact_paths = parse_h4_test_args(ARGS)
if !isempty(artifact_paths)
    required = Set(["d5-n4", "d5-n6", "d4-n100", "d5-n100"])
    Set(keys(artifact_paths)) == required ||
        error("artifact replay requires --d5-n4, --d5-n6, --d4-n100, and --d5-n100")
    all(isfile, values(artifact_paths)) || error("one or more artifact paths do not exist")

    d5_n4 = validate_d5_n4_probe(artifact_paths["d5-n4"])
    d5_n6 = validate_d5_n6_result(artifact_paths["d5-n6"])
    validate_d4_n100_result(artifact_paths["d4-n100"])
    d5_n100 = validate_h4_n100_result(artifact_paths["d5-n100"])

    @test d5_n4["gates"]["correctness_passed"]
    @test d5_n6["gates"]["residual_resolved_improvement_passed"]
    @test d5_n100["gates"]["N200_entry_passed"]
    @test d5_n100["gates"]["improvement_to_residual_ratio"] > 100

    run_name, n_super, threads, parsed_n100_path, previous_path = parse_h4_n200_args([
        "--run-name=test-h4-n200",
        "--n-super=50",
        "--mosek-threads=8",
        "--n100-result=$(artifact_paths["d5-n100"])",
    ])
    @test run_name == "test-h4-n200"
    @test n_super == 50
    @test threads == 8
    @test realpath(parsed_n100_path) == realpath(artifact_paths["d5-n100"])
    @test isnothing(previous_path)
    @test_throws ErrorException parse_h4_n200_args([
        "--run-name=test-h4-n200-missing-predecessor",
        "--n-super=75",
        "--n100-result=$(artifact_paths["d5-n100"])",
    ])
end

println("H4 regression checks passed")
