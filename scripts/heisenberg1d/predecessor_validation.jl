function require_predecessor_metadata(result, key, expected, label)
    actual = get(get(result, "meta", Dict{String, Any}()), key, nothing)
    actual == expected || error("$label has meta.$key=$actual, expected $expected")
end

function validate_h0_predecessor(h0; system_size=20)
    require_predecessor_metadata(h0, "stage", "H0", "H0 predecessor")
    require_predecessor_metadata(h0, "system_scope", "finite_periodic_chain", "H0 predecessor")
    require_predecessor_metadata(h0, "system_size", system_size, "H0 predecessor")
    require_predecessor_metadata(h0, "boundary", "periodic", "H0 predecessor")
    require_predecessor_metadata(h0, "numerical_bound_schema_version", NUMERICAL_BOUND_SCHEMA_VERSION, "H0 predecessor")
    exact_energy = h0["exact_reference"]["energy_per_site"]
    compressed = h0["compressed"]
    [point["n_super"] for point in compressed] == [4, 6, 8, 10] || error("H0 predecessor has unexpected hierarchy depths")
    for point in compressed
        record = Dict(
            "status" => point["status"],
            "energy" => point["energy_per_site"],
            "residuals" => point["residuals"],
            "residual_adjusted_lower_bound" => point["residual_adjusted_lower_bound"],
            "published_lower_bound" => point["published_lower_bound"],
        )
        validated_solver_bounds(record, "H0 n=$(point["n_super"]) point")
    end
    all_optimal = all(point["status"] == "OPTIMAL" for point in compressed)
    all_below = all(
        point["system_scope"] == "finite_periodic_chain" &&
        point["system_size"] == system_size &&
        point["boundary"] == "periodic" &&
        point["D"] == 3 &&
        point["support_sites"] == 2 * point["n_super"] <= system_size &&
        point["published_lower_bound"] <= point["residual_adjusted_lower_bound"] <= point["energy_per_site"] < exact_energy
        for point in compressed
    )
    monotone = all(diff([point["energy_per_site"] for point in compressed]) .>= -2e-7)
    structured = h0["structured_npa"]
    structured_valid = structured["status"] == "OPTIMAL" && structured["certified_energy_per_site"] <= exact_energy
    verification = h0["verification"]
    all_optimal && get(verification, "all_compressed_optimal", false) || error("H0 OPTIMAL gate is inconsistent with its points")
    all_below && get(verification, "all_bounds_below_exact", false) || error("H0 sandwich gate is inconsistent with its points")
    monotone && get(verification, "compressed_monotone", false) || error("H0 monotonicity gate is inconsistent with its points")
    structured_valid && get(verification, "structured_rational_postcertificate", false) || error("H0 structured certificate gate is inconsistent")
    h0
end

function validate_h1_predecessor(h1; bond_dimension=H1_D, n_super=6)
    require_predecessor_metadata(h1, "stage", "H1", "H1 predecessor")
    require_predecessor_metadata(h1, "system_scope", "infinite_translation_invariant", "H1 predecessor")
    require_predecessor_metadata(h1, "system_size", nothing, "H1 predecessor")
    require_predecessor_metadata(h1, "bond_dimension", bond_dimension, "H1 predecessor")
    require_predecessor_metadata(h1, "n_super", n_super, "H1 predecessor")
    require_predecessor_metadata(h1, "numerical_bound_schema_version", NUMERICAL_BOUND_SCHEMA_VERSION, "H1 predecessor")
    solves = h1["solves"]
    validated_solver_bounds(solves["dense"], "H1 dense solve")
    validated_solver_bounds(solves["u1_blocks"], "H1 U(1)-blocked solve")
    both_optimal = solves["dense"]["status"] == "OPTIMAL" && solves["u1_blocks"]["status"] == "OPTIMAL"
    energy_difference = abs(solves["dense"]["energy"] - solves["u1_blocks"]["energy"])
    charge_check = h1["charge_check"]
    charge_valid = charge_check["exactified_to_zero"] &&
        charge_check["largest_forbidden_entry"] <= H1_TOL &&
        all(abs(value) <= H1_TOL for value in values(charge_check["map_violations"]))
    gates = h1["gates"]
    both_optimal && get(gates, "both_optimal", false) || error("H1 OPTIMAL gate is inconsistent")
    energy_difference <= gates["energy_agreement_atol"] && get(gates, "correctness_passed", false) || error("H1 correctness gate is inconsistent")
    charge_valid || error("H1 charge-preservation checks failed")
    get(gates, "usefulness_passed", false) || error("H1 usefulness gate did not pass")
    h1
end
