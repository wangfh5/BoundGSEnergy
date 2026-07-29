using Test

include("h3_finite_n100_d4.jl")

legacy = Dict{String, Any}(
    "D4-n4" => Dict{String, Any}(
        "status" => "OPTIMAL",
        "E" => -0.45,
        "gap" => 0.01,
        "physical_n" => 8,
        "super_n" => 4,
        "coarse_grainer_fingerprint" => "test",
    ),
    "LTI" => Dict{String, Any}(
        "4" => Dict{String, Any}("status" => "OPTIMAL", "E" => -0.5),
    ),
    "plateaus" => Dict{String, Any}(
        "D4" => Dict{String, Any}("largest_super_n" => 4, "largest_physical_n" => 8),
    ),
)

parsed_legacy = JSON.parse(JSON.json(legacy))
migrate_result_schema!(parsed_legacy)
migrate_result_schema!(parsed_legacy)
point = parsed_legacy["D4-n4"]
@assert !haskey(point, "physical_n")
@assert !haskey(point, "super_n")
@assert point["support_sites"] == 8
@assert point["n_super"] == 4
@assert point["system_scope"] == SYSTEM_SCOPE
@assert isnothing(point["system_size"])
@assert parsed_legacy["LTI"]["4"]["support_sites"] == 4
@assert parsed_legacy["LTI"]["4"]["system_scope"] == SYSTEM_SCOPE
@assert parsed_legacy["plateaus"]["D4"]["largest_support_sites"] == 8
@assert parsed_legacy["meta"]["schema_version"] == RESULT_SCHEMA_VERSION
@assert validate_cached_point(parsed_legacy, "D4-n4", "test", true)

new_point = result_record("OPTIMAL", 1.0, "test"; E=-0.45, gap=0.01, support_sites=8, n_super=4)
@assert new_point["support_sites"] == 8
@assert new_point["n_super"] == 4
@assert new_point["system_scope"] == SYSTEM_SCOPE
@assert isnothing(new_point["system_size"])
roundtrip = JSON.parse(JSON.json(new_point))
@assert haskey(roundtrip, "system_size")
@assert isnothing(roundtrip["system_size"])

conflicting = Dict{String, Any}(
    "D4-n4" => Dict{String, Any}("physical_n" => 10),
)
@test_throws ErrorException migrate_result_schema!(conflicting)

conflicting_n_super = Dict{String, Any}(
    "D4-n4" => Dict{String, Any}("n_super" => 5),
)
@test_throws ErrorException migrate_result_schema!(conflicting_n_super)

missing_n_super = Dict{String, Any}(
    "D4-n4" => Dict{String, Any}("support_sites" => 8),
)
migrate_result_schema!(missing_n_super)
@assert missing_n_super["D4-n4"]["n_super"] == 4

numeric = Dict{String, Any}(
    "meta" => Dict{String, Any}(),
    "solve" => Dict{String, Any}(
        "status" => "OPTIMAL",
        "energy" => -0.45,
        "residuals" => Dict{String, Any}(
            "primal_feasibility" => 2e-8,
            "dual_feasibility" => 1e-9,
            "primal_objective" => -0.45,
            "dual_objective" => -0.44999999,
            "relative_gap" => 1e-8,
        ),
    ),
)
migrate_numerical_bound_schema!(numeric)
migrate_numerical_bound_schema!(numeric)
@assert numeric["meta"]["numerical_bound_schema_version"] == NUMERICAL_BOUND_SCHEMA_VERSION
@assert numeric["solve"]["published_lower_bound"] <= numeric["solve"]["residual_adjusted_lower_bound"] <= numeric["solve"]["energy"]

function parse_artifact_args(args)
    paths = Dict{String, String}()
    for arg in args
        startswith(arg, "--") && occursin("=", arg) || error("expected --name=path, got $arg")
        name, path = split(arg[3:end], "=", limit=2)
        paths[name] = normpath(isabspath(path) ? path : joinpath(@__DIR__, "..", "..", path))
    end
    paths
end

artifact_paths = parse_artifact_args(ARGS)
if !isempty(artifact_paths)
    required = Set(["h0", "h1", "n20", "n50", "n100", "d3"])
    Set(keys(artifact_paths)) == required ||
        error("artifact replay requires --h0, --h1, --n20, --n50, --n100, and --d3")
    all(isfile, values(artifact_paths)) || error("one or more artifact paths do not exist")

    h0_path = artifact_paths["h0"]
    h1_path = artifact_paths["h1"]
    n20_path = artifact_paths["n20"]
    n50_path = artifact_paths["n50"]
    n100_path = artifact_paths["n100"]
    legacy_path = artifact_paths["d3"]

    h0 = validate_h0_predecessor(JSON.parsefile(h0_path))
    h1 = validate_h1_predecessor(JSON.parsefile(h1_path))
    n20 = validate_h3_predecessor(JSON.parsefile(n20_path), H3_N, H3_DEPTHS)
    n50 = validate_h3_predecessor(JSON.parsefile(n50_path), H3_N50, H3_N50_DEPTHS)
    validate_h3_predecessor(JSON.parsefile(n100_path), H3_N100, H3_N100_DEPTHS)
    legacy_result = JSON.parsefile(legacy_path)
    validate_n20_chain(n20, h0, h1)
    validate_n50_chain(n50, h0, n20, legacy_result)

    wrong_system = deepcopy(h0)
    wrong_system["meta"]["system_size"] = H3_N + 2
    @test_throws ErrorException validate_h0_predecessor(wrong_system)

    wrong_status = deepcopy(n20)
    wrong_status["points"][1]["solve"]["status"] = "INFEASIBLE"
    @test_throws ErrorException validate_h3_predecessor(wrong_status, H3_N, H3_DEPTHS)

    wrong_d3_path = deepcopy(n50)
    wrong_d3_path["meta"]["inputs"]["d3_result"] = relpath(h0_path, H3_REPO_ROOT)
    @test_throws ErrorException assert_recorded_input(wrong_d3_path, "d3_result", legacy_path)

    wrong_depths = deepcopy(h0)
    wrong_depths["compressed"][1]["n_super"] = 3
    @test_throws ErrorException validate_h0_predecessor(wrong_depths)

    wrong_dimension = deepcopy(h0)
    wrong_dimension["compressed"][1]["D"] = 2
    @test_throws ErrorException validate_h0_predecessor(wrong_dimension)

    stale_snapshot = deepcopy(n20)
    stale_snapshot["comparison"]["D3_endpoint_energy_per_site"] -= 1e-3
    @test_throws ErrorException validate_n20_chain(stale_snapshot, h0, h1)

    stale_prior_point = deepcopy(n50)
    stale_prior_point["comparison"]["D4_n10_N20_provenance"]["energy_per_site"] -= 1e-3
    @test_throws ErrorException validate_n50_chain(stale_prior_point, h0, n20, legacy_result)

    nonoptimal_legacy = deepcopy(legacy_result)
    nonoptimal_legacy["D3-n20"]["status"] = "INFEASIBLE"
    matched_nonoptimal_n50 = deepcopy(n50)
    matched_nonoptimal_n50["comparison"]["D3_n20"]["status"] = "INFEASIBLE"
    @test_throws ErrorException validate_n50_chain(matched_nonoptimal_n50, h0, n20, nonoptimal_legacy)

    stale_bound = deepcopy(n20)
    stale_bound["points"][1]["residual_adjusted_lower_bound"] -= 1e-6
    @test_throws ErrorException validate_h3_predecessor(stale_bound, H3_N, H3_DEPTHS)

    tolerance_valid_h1 = deepcopy(h1)
    tolerance_valid_h1["charge_check"]["map_violations"]["W2"] = H1_TOL / 2
    validate_h1_predecessor(tolerance_valid_h1)
end

A1 = reshape(Float64.(1:40), 5, 2, 4) ./ 40
A2 = reshape(Float64.(1:40), 4, 2, 5) ./ 40
unequal_bond_tensor = super_tensor((AL=[A1, A2],))
@test size(unequal_bond_tensor) == (5, 4, 5)
for s1 in 1:2, s2 in 1:2
    @test unequal_bond_tensor[:, (s1 - 1) * 2 + s2, :] ≈ A1[:, s1, :] * A2[:, s2, :]
end

println("metadata migration checks passed")
