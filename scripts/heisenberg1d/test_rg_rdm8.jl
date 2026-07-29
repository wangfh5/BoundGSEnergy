using Test

include("rg_rdm8_finite_n20.jl")

virtual_charges = [-1, 1, -1, 1]
super_charges = [-2, 0, 0, 2]
B_raw, _ = build_coarse_grainer()
B_reference, _ = exact_charge_tensor(B_raw, virtual_charges, super_charges)
B = dyadic_coarse_grainer(B_reference, RDM8_DYADIC_BITS)
audit = dyadic_relaxation_audit(B, RDM8_DYADIC_BITS)

@test audit["status"] == "PROVEN"
@test audit["W2_exact_rational_contraction"]
@test audit["quadratic_coefficient_assembly_exact"]
@test audit["physical_relaxation_compatibility_proven"]
@test norm(B - B_reference) / norm(B_reference) < 1e-4
@test maximum(values(assert_charge_preserving_maps(B, virtual_charges, super_charges))) == 0
@test maximum(values(verify_rdm4_marginal_identities(W2map(B)))) <= 1e-12
@test_throws ErrorException dyadic_relaxation_audit(B_reference, RDM8_DYADIC_BITS)
@test_throws ErrorException dyadic_relaxation_audit(
    dyadic_coarse_grainer(B_reference, RDM8_DYADIC_BITS + 1),
    RDM8_DYADIC_BITS + 1,
)

function test_certified_rdm8_result(path; finite=false)
    result = JSON.parsefile(path)
    comparison = result["comparison"]
    baseline_certificate = validated_dual_certificate(comparison["baseline"], "RG/RDM8 baseline")
    hybrid_certificate = validated_dual_certificate(comparison["hybrid"], "RG/RDM8 hybrid")
    exact_difference = parse_exact_rational(hybrid_certificate["certified_lower_bound_rational"]) -
        parse_exact_rational(baseline_certificate["certified_lower_bound_rational"])
    @test result["meta"]["stage"] == "H7"
    @test result["meta"]["dyadic_bits"] == RDM8_DYADIC_BITS
    @test comparison["certified_lower_bound_difference_rational"] == string(exact_difference)
    @test comparison["certified_lower_bound_difference"] == Float64(exact_difference)
    @test result["gates"]["physical_relaxation_compatibility_proven"]
    @test result["gates"]["stage_passed"]
    if finite
        upper_certificate = validated_finite_upper_certificate(
            result["finite_energy_upper_certificate"],
            "RG/RDM8 finite upper bound",
        )
        upper = parse_exact_rational(upper_certificate["upper_bound_rational"])
        baseline_gap = upper - parse_exact_rational(baseline_certificate["published_certified_lower_bound_rational"])
        hybrid_gap = upper - parse_exact_rational(hybrid_certificate["published_certified_lower_bound_rational"])
        @test comparison["baseline_certified_gap_rational"] == string(baseline_gap)
        @test comparison["hybrid_certified_gap_rational"] == string(hybrid_gap)
    end
end

function parse_rg_test_args(args)
    paths = Dict{String, String}()
    for arg in args
        startswith(arg, "--") && occursin("=", arg) || error("expected --name=path, got $arg")
        name, path = split(arg[3:end], "=", limit=2)
        paths[name] = normpath(isabspath(path) ? path : joinpath(@__DIR__, "..", "..", path))
    end
    paths
end

artifact_paths = parse_rg_test_args(ARGS)
if !isempty(artifact_paths)
    required = Set(["h3", "local", "finite"])
    Set(keys(artifact_paths)) == required ||
        error("artifact replay requires --h3, --local, and --finite")
    all(isfile, values(artifact_paths)) || error("one or more artifact paths do not exist")

    h3 = validate_h3_n20_endpoint(artifact_paths["h3"])
    @test h3["gates"]["all_below_exact"]

    stale_schema = deepcopy(h3)
    stale_schema["meta"]["numerical_bound_schema_version"] = 0
    @test_throws ErrorException validate_h3_n20_endpoint(stale_schema)

    failed_gate = deepcopy(h3)
    failed_gate["gates"]["all_optimal"] = false
    @test_throws ErrorException validate_h3_n20_endpoint(failed_gate)

    altered_reference = deepcopy(h3)
    altered_reference["exact_reference"]["energy_per_site"] += 1e-4
    @test_throws ErrorException validate_h3_n20_endpoint(altered_reference)

    poorer_endpoint = deepcopy(h3)
    last(poorer_endpoint["points"])["energy_per_site"] =
        poorer_endpoint["comparison"]["D3_endpoint_energy_per_site"]
    @test_throws ErrorException validate_h3_n20_endpoint(poorer_endpoint)

    run_name, threads, predecessor, tolerance, silent = parse_rg_rdm8_finite_args([
        "--run-name=test-rg-rdm8-finite",
        "--mosek-threads=8",
        "--h3-result=$(artifact_paths["h3"])",
    ])
    @test run_name == "test-rg-rdm8-finite"
    @test threads == 8
    @test realpath(predecessor) == realpath(artifact_paths["h3"])
    @test tolerance == RG_FINITE_FEAS_TOL
    @test silent

    test_certified_rdm8_result(artifact_paths["local"])
    test_certified_rdm8_result(artifact_paths["finite"]; finite=true)
end

println("RG/RDM8 compatibility and certificate checks passed")
