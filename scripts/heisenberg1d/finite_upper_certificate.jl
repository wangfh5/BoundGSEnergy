using XDiag

function periodic_heisenberg_operator(system_size)
    hamiltonian = XDiag.OpSum()
    for site in 1:system_size
        hamiltonian += XDiag.Op("SdotS", [site, mod1(site + 1, system_size)])
    end
    hamiltonian
end

function dyadic_parts(value::Float64)
    rational = exact_rational(value)
    denominator_value = denominator(rational)
    ispow2(denominator_value) || error("a Float64 did not convert to a dyadic rational")
    numerator(rational), ndigits(denominator_value; base=2) - 1
end

function exact_rayleigh_per_site(csr, vector, system_size)
    csr.i0 == 1 || error("the exact Rayleigh evaluator requires one-based CSR indices")
    length(vector) == csr.nrows == csr.ncols || error("CSR matrix and vector dimensions do not match")
    vector_parts = dyadic_parts.(vector)
    numerators = first.(vector_parts)
    exponents = last.(vector_parts)
    maximum_vector_exponent = maximum(exponents)

    norm_numerator = BigInt(0)
    for index in eachindex(vector)
        norm_numerator += numerators[index]^2 << (2 * maximum_vector_exponent - 2 * exponents[index])
    end
    norm_numerator > 0 || error("the upper-certificate vector is zero")

    data_parts = Dict(value => dyadic_parts(value) for value in unique(csr.data))
    maximum_data_exponent = maximum(last, values(data_parts))
    energy_exponent = 2 * maximum_vector_exponent + maximum_data_exponent
    energy_numerator = BigInt(0)
    for row in 1:csr.nrows, cursor in csr.rowptr[row]:csr.rowptr[row + 1] - 1
        column = csr.col[cursor]
        matrix_numerator, matrix_exponent = data_parts[csr.data[cursor]]
        shift = energy_exponent - exponents[row] - matrix_exponent - exponents[column]
        energy_numerator += numerators[row] * matrix_numerator * numerators[column] << shift
    end
    energy = energy_numerator // (BigInt(1) << energy_exponent)
    norm_squared = norm_numerator // (BigInt(1) << (2 * maximum_vector_exponent))
    energy / norm_squared / system_size
end

function build_finite_upper_certificate(system_size, transcript_path; precision=1e-12, random_seed=20260729)
    iseven(system_size) || error("the periodic Heisenberg upper certificate requires even system size")
    hamiltonian = periodic_heisenberg_operator(system_size)
    block = XDiag.Spinhalf(system_size, system_size ÷ 2)
    lanczos_value, state = XDiag.eig0(
        hamiltonian,
        block;
        precision,
        random_seed,
    )
    vector = XDiag.vector(state)
    csr = XDiag.csr_matrix(hamiltonian, block)
    rayleigh = exact_rayleigh_per_site(csr, vector, system_size)
    transcript = Dict(
        "version" => 1,
        "system_size" => system_size,
        "sector" => "Sz_total=0",
        "vector" => vector,
        "precision" => precision,
        "random_seed" => random_seed,
    )
    atomic_serialize(transcript_path, transcript)
    Dict(
        "status" => "CERTIFIED_UPPER_BOUND",
        "method" => "exact dyadic Rayleigh quotient of a stored XDiag Lanczos vector",
        "system_size" => system_size,
        "sector" => "Sz_total=0",
        "lanczos_energy_per_site" => lanczos_value / system_size,
        "upper_bound" => float_up(rayleigh),
        "upper_bound_rational" => string(rayleigh),
        "replay_transcript" => artifact_relative_path(transcript_path),
        "replay_transcript_sha256" => file_sha256(transcript_path),
        "vector_sha256" => bytes2hex(sha256(reinterpret(UInt8, vector))),
    )
end

function validated_finite_upper_certificate(certificate, label)
    certificate["status"] == "CERTIFIED_UPPER_BOUND" || error("$label upper certificate did not pass")
    transcript_path = normpath(joinpath(BOUNDGSENERGY_ROOT, certificate["replay_transcript"]))
    isfile(transcript_path) || error("$label upper-certificate transcript does not exist")
    file_sha256(transcript_path) == certificate["replay_transcript_sha256"] ||
        error("$label upper-certificate transcript hash is stale")
    transcript = deserialize(transcript_path)
    transcript["version"] == 1 || error("$label uses an unsupported upper-certificate version")
    transcript["system_size"] == certificate["system_size"] || error("$label system size is stale")
    transcript["sector"] == certificate["sector"] || error("$label sector is stale")
    vector = transcript["vector"]
    bytes2hex(sha256(reinterpret(UInt8, vector))) == certificate["vector_sha256"] ||
        error("$label vector hash is stale")
    hamiltonian = periodic_heisenberg_operator(certificate["system_size"])
    block = XDiag.Spinhalf(certificate["system_size"], certificate["system_size"] ÷ 2)
    rayleigh = exact_rayleigh_per_site(
        XDiag.csr_matrix(hamiltonian, block),
        vector,
        certificate["system_size"],
    )
    rayleigh == parse_exact_rational(certificate["upper_bound_rational"]) ||
        error("$label exact Rayleigh quotient is stale")
    certificate["upper_bound"] == float_up(rayleigh) || error("$label numeric upper bound is stale")
    certificate
end
