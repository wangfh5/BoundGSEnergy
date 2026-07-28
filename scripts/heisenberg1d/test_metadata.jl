using Test

include("t2t3_compressed_lti.jl")

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

println("metadata migration checks passed")
