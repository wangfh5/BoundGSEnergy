# Shared structured-NPA + exact-dyadic RG formulation.

using Dualization

include(joinpath(@__DIR__, "h11_full_sohs_certificate.jl"))
include(joinpath(@__DIR__, "rdm8_common.jl"))

const STRUCTURED_RG_RDM9_BLOCKS = [
    [index for index in 1:512 if count_ones(index - 1) == weight]
    for weight in 0:4
]
const STRUCTURED_RG_RDM9_DIMENSIONS = length.(STRUCTURED_RG_RDM9_BLOCKS)
const STRUCTURED_RG_EQUALITY_SCALE = 1.0
const STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE = 1.0
const STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE = 16.0
const STRUCTURED_RG_RDM8_RECURSION_SCALE = 1.0
const STRUCTURED_RG_VIRTUAL_CHARGES = [-1, 1, -1, 1]
const STRUCTURED_RG_SUPER_CHARGES = [-2, 0, 0, 2]
const STRUCTURED_RG_OMEGA_NORMALIZATIONS = (:none, :dyadic_trace)

function structured_rg_rdm8_link_scales(omega_normalization)
    omega_normalization in STRUCTURED_RG_OMEGA_NORMALIZATIONS ||
        error("omega_normalization must be none or dyadic_trace")
    omega_normalization == :none ?
        (
            local_link=STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE,
            omega4=STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE,
            recursion=STRUCTURED_RG_RDM8_RECURSION_SCALE,
        ) :
        (
            local_link=STRUCTURED_RG_RDM8_LOCAL_LINK_SCALE,
            omega4=2 * STRUCTURED_RG_RDM8_OMEGA4_LINK_SCALE,
            recursion=8 * STRUCTURED_RG_RDM8_RECURSION_SCALE,
        )
end

function structured_rg_omega_trace_scales(B, n_super, normalization)
    normalization in STRUCTURED_RG_OMEGA_NORMALIZATIONS ||
        error("omega_normalization must be none or dyadic_trace")
    normalization == :none && return ones(n_super - 3)
    bounds = compressed_trace_bounds(B, n_super)
    [
        let
            bound = bounds["omega$(depth)"]
            scale = one(BigInt)
            while scale * denominator(bound) < numerator(bound)
                scale <<= 1
            end
            exact_scale = ExactRational(scale)
            exact_scale >= bound &&
                (exact_scale == 1 || exact_scale / 2 < bound) ||
                error("failed to construct the minimal dyadic trace scale")
            floating_scale = Float64(scale)
            isfinite(floating_scale) ||
                error("dyadic trace scale exceeds Float64 range")
            floating_scale
        end
        for depth in 4:n_super
    ]
end

function validate_structured_rg_support(n_super; minimum_depth=4)
    n_super >= minimum_depth ||
        error("the structured-RG hierarchy requires n_super >= $minimum_depth")
    2 * n_super <= H10_N ||
        error("the structured-RG support exceeds the finite N=$H10_N system")
    n_super
end

function build_structured_moment_model()
    configuration = h11_configuration(false)
    Random.seed!(H11_RANDOM_SEED)
    output = redirect_stdout(devnull) do
        QMBCertify.GSB_with_model(
            [[1, 4]],
            [3 / 4],
            H10_N,
            configuration.order;
            lattice="chain",
            rdm=configuration.rdm,
            extra=configuration.extra,
            pso=configuration.pso,
            lso=configuration.lso,
            lol=H10_N,
            three_type=[1, 1],
            Gram=true,
            QUIET=true,
            mosek_setting=QMBCertify.mosek_para(1e-8, 1e-8, 1e-8, 1),
            h11_build_only=true,
        )
    end
    _, data, sohs_model = output
    model = dualize(sohs_model)
    coefficient_matching = only(all_constraints(sohs_model, Vector{AffExpr}, MOI.Zeros))
    scaled_moments = Dualization._get_dual_variables(model, coefficient_matching)
    length(scaled_moments) == length(data.tsupp) ||
        error("structured moment variables do not match the canonical word support")
    physical_moments = -H11_EQUALITY_SCALE .* scaled_moments
    (
        model=model,
        sohs_model=sohs_model,
        canonical_words=data.tsupp,
        scaled_moments=scaled_moments,
        physical_moments=physical_moments,
        configuration=configuration,
    )
end

function structured_rdm9_constraints(model, sohs_model)
    sohs_constraints =
        all_constraints(sohs_model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
    length(sohs_constraints) >= length(STRUCTURED_RG_RDM9_DIMENSIONS) ||
        error("structured moment model is missing its RDM9 blocks")
    rdm9_sohs_constraints = sohs_constraints[end-4:end]
    dimensions = [
        MOI.side_dimension(constraint_object(constraint).set)
        for constraint in rdm9_sohs_constraints
    ]
    dimensions == STRUCTURED_RG_RDM9_DIMENSIONS ||
        error("structured moment model has an unexpected RDM9 block layout")
    [
        let
            object = constraint_object(constraint)
            block = reshape_vector(object.func, object.shape)
            dual_constraint, _ = Dualization._get_dual_constraint(model, block[1, 1])
            dual_constraint isa ConstraintRef ||
                error("RDM9 SOHS block is not mapped to a structured moment constraint")
            dual_constraint
        end
        for constraint in rdm9_sohs_constraints
    ]
end

function foreach_local_pauli_entry(
    callback,
    canonical_words,
    site_count,
)
    site_count >= 1 || error("site_count must be positive")
    lookup = Dict(Tuple(word) => index for (index, word) in enumerate(canonical_words))
    dimension = 2^site_count
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
        iszero(coefficient) && continue
        index = get(lookup, Tuple(reduced), nothing)
        isnothing(index) &&
            error("structured support is missing a local Pauli word")

        for column_bits in 0:(dimension - 1)
            row_bits = column_bits
            amplitude = one(ComplexF64)
            for site in 1:site_count
                bit_position = site_count - site
                bit = (column_bits >> bit_position) & 1
                label = labels[site]
                if label == 1
                    row_bits ⊻= 1 << bit_position
                elseif label == 2
                    row_bits ⊻= 1 << bit_position
                    amplitude *= bit == 0 ? im : -im
                elseif label == 3
                    amplitude *= bit == 0 ? 1 : -1
                end
            end
            entry_coefficient = amplitude * coefficient / dimension
            abs(imag(entry_coefficient)) <= 1e-12 ||
                error("structured local RDM has a complex coefficient")
            callback(
                index,
                Float64(real(entry_coefficient)),
                row_bits + 1,
                column_bits + 1,
            )
        end
    end
    nothing
end

function structured_local_rdm(
    physical_moments,
    canonical_words,
    site_count,
)
    length(physical_moments) == length(canonical_words) ||
        error("physical moments do not match the canonical word support")
    dimension = 2^site_count
    rho = [AffExpr(0.0) for _ in 1:dimension, _ in 1:dimension]
    foreach_local_pauli_entry(canonical_words, site_count) do index, coefficient, row, column
        add_to_expression!(
            rho[row, column],
            coefficient,
            physical_moments[index],
        )
    end
    rho
end

function structured_local_rdm_exactness_audit(
    rho,
    physical_moments,
    canonical_words,
    site_count,
)
    dimension = 2^site_count
    size(rho) == (dimension, dimension) ||
        error("structured local RDM has the wrong dimension")
    length(physical_moments) == length(canonical_words) ||
        error("physical moments do not match the canonical word support")
    scale = big(2)^site_count
    moment_variables = VariableRef[]
    moment_scales = ExactRational[]
    for moment in physical_moments
        iszero(moment.constant) ||
            error("a physical moment has an unexpected affine constant")
        terms = collect(linear_terms(moment))
        length(terms) == 1 ||
            error("a physical moment is not a scaled model variable")
        coefficient, variable = only(terms)
        push!(moment_variables, variable)
        push!(moment_scales, exact_rational(coefficient))
    end

    addend_count = 0
    maximum_addend_numerator = big(0)
    expected = [
        Dict{VariableRef, ExactRational}()
        for _ in 1:dimension, _ in 1:dimension
    ]
    foreach_local_pauli_entry(canonical_words, site_count) do index, coefficient, row, column
        exact_coefficient =
            exact_rational(coefficient) * moment_scales[index]
        scaled = exact_coefficient * scale
        denominator(scaled) == 1 ||
            error("a structured local-RDM addend is not on the exact dyadic grid")
        maximum_addend_numerator =
            max(maximum_addend_numerator, abs(numerator(scaled)))
        variable = moment_variables[index]
        expected[row, column][variable] =
            get(expected[row, column], variable, ExactRational(0)) +
            exact_coefficient
        addend_count += 1
    end

    maximum_final_numerator = big(0)
    for index in eachindex(rho)
        expression = rho[index]
        scaled_constant = exact_rational(expression.constant) * scale
        denominator(scaled_constant) == 1 ||
            error("a structured local-RDM constant is not exactly dyadic")
        maximum_final_numerator =
            max(maximum_final_numerator, abs(numerator(scaled_constant)))
        actual = Dict(
            variable => exact_rational(coefficient)
            for (coefficient, variable) in linear_terms(expression)
            if !iszero(coefficient)
        )
        filter!(pair -> !iszero(last(pair)), expected[index])
        actual == expected[index] ||
            error("Float64 local-RDM assembly differs from exact-rational reconstruction")
        for (coefficient, _) in linear_terms(expression)
            scaled_coefficient = exact_rational(coefficient) * scale
            denominator(scaled_coefficient) == 1 ||
                error("a structured local-RDM coefficient is not exactly dyadic")
            maximum_final_numerator =
                max(maximum_final_numerator, abs(numerator(scaled_coefficient)))
        end
    end

    accumulation_bound =
        big(4)^site_count * maximum_addend_numerator
    exact_integer_limit = big(2)^53
    accumulation_bound < exact_integer_limit ||
        error("structured local-RDM assembly exceeds the Float64 exact-integer range")
    maximum_final_numerator < exact_integer_limit ||
        error("structured local-RDM coefficient exceeds the Float64 exact-integer range")

    (
        status="PROVEN",
        site_count=site_count,
        grid_denominator=string(scale),
        physical_moment_scale=string(exact_rational(H11_EQUALITY_SCALE)),
        addend_count=addend_count,
        maximum_addend_numerator=string(maximum_addend_numerator),
        accumulation_numerator_bound=string(accumulation_bound),
        maximum_final_numerator=string(maximum_final_numerator),
        float64_exact_integer_limit=string(exact_integer_limit),
        exact_rational_reconstruction_matches=true,
        coefficient_assembly_exact=true,
    )
end

function add_scaled_symmetric_equalities!(
    model,
    residual;
    scale=STRUCTURED_RG_EQUALITY_SCALE,
)
    size(residual, 1) == size(residual, 2) ||
        error("RG marginal residual must be square")
    count = 0
    for row in axes(residual, 1), column in row:last(axes(residual, 2))
        iszero(residual[row, column]) && continue
        @constraint(model, scale * residual[row, column] == 0)
        count += 1
    end
    count
end

function add_rg_hierarchy!(
    model,
    structured_rho3,
    B,
    n_super,
    virtual_charges,
    super_charges,
)
    validate_structured_rg_support(n_super)
    size(structured_rho3) == (64, 64) ||
        error("the structured RDM must cover six physical spins")
    D = size(B, 1)
    D == size(B, 3) || error("the RG coarse-grainer must have equal boundary dimensions")

    rho3 = psd_matrix(
        model,
        64,
        "hybrid_rho3",
        rho3_charges(super_charges),
    )
    rho3_link_count =
        add_scaled_symmetric_equalities!(model, rho3 - structured_rho3)

    boundary_charges = omega_charges(super_charges, virtual_charges)
    W2 = W2map(B)
    omega = Vector{Any}(undef, n_super - 3)
    omega[1] = psd_matrix(
        model,
        16 * D^2,
        "hybrid_omega4",
        boundary_charges,
    )
    right_rho3_map = kron(I4, W2)
    left_rho3_map = kron(W2, I4)
    dw = 4 * D^2
    add_scaled_symmetric_equalities!(
        model,
        ptrace_right(omega[1], dw, 4) -
        right_rho3_map * rho3 * right_rho3_map',
    )
    add_scaled_symmetric_equalities!(
        model,
        ptrace_left(omega[1], 4, dw) -
        left_rho3_map * rho3 * left_rho3_map',
    )

    Iw = Matrix{Float64}(I, 4 * D, 4 * D)
    right_map = kron(Iw, MRmap(B))
    left_map = kron(MLmap(B), Iw)
    domega = 16 * D^2

    for (index, depth) in enumerate(5:n_super)
        next_omega = psd_matrix(model, domega, "hybrid_omega$(depth)", boundary_charges)
        right_residual =
            ptrace_right(next_omega, dw, 4) -
            right_map * omega[index] * right_map'
        left_residual =
            ptrace_left(next_omega, 4, dw) -
            left_map * omega[index] * left_map'
        add_scaled_symmetric_equalities!(model, right_residual)
        add_scaled_symmetric_equalities!(model, left_residual)
        omega[index + 1] = next_omega
    end
    (
        rho3=rho3,
        omega=omega,
        rho3_link_equality_count=rho3_link_count,
    )
end

function add_rdm8_rg_hierarchy!(
    model,
    structured_rho4,
    B,
    n_super,
    virtual_charges,
    super_charges,
    ;
    local_link_scale=nothing,
    omega4_link_scale=nothing,
    recursion_scale=nothing,
    omega_normalization=:none,
)
    validate_structured_rg_support(n_super; minimum_depth=5)
    size(structured_rho4) == (256, 256) ||
        error("the structured RDM must cover eight physical spins")
    D = size(B, 1)
    D == size(B, 3) ||
        error("the RG coarse-grainer must have equal boundary dimensions")
    default_link_scales =
        structured_rg_rdm8_link_scales(omega_normalization)
    isnothing(local_link_scale) &&
        (local_link_scale = default_link_scales.local_link)
    isnothing(omega4_link_scale) &&
        (omega4_link_scale = default_link_scales.omega4)
    isnothing(recursion_scale) &&
        (recursion_scale = default_link_scales.recursion)

    rho4 = psd_matrix(
        model,
        256,
        "hybrid_rho4",
        rho4_charges(super_charges),
    )
    rho4_link_count =
        add_scaled_symmetric_equalities!(
            model,
            rho4 - structured_rho4;
            scale=local_link_scale,
        )

    boundary_charges = omega_charges(super_charges, virtual_charges)
    W2 = W2map(B)
    omega = Vector{Any}(undef, n_super - 3)
    omega_trace_scales = structured_rg_omega_trace_scales(
        B,
        n_super,
        omega_normalization,
    )
    omega[1] = psd_matrix(
        model,
        16 * D^2,
        "hybrid_omega4",
        boundary_charges,
    )
    compression = kron(I4, W2, I4)
    omega4_link_count = add_scaled_symmetric_equalities!(
        model,
        omega[1] -
        compression * rho4 * compression' / omega_trace_scales[1],
        scale=omega4_link_scale,
    )

    Iw = Matrix{Float64}(I, 4 * D, 4 * D)
    right_map = kron(Iw, MRmap(B))
    left_map = kron(MLmap(B), Iw)
    dw = 4 * D^2
    domega = 16 * D^2
    for (index, depth) in enumerate(5:n_super)
        propagation_scale =
            omega_trace_scales[index] /
            omega_trace_scales[index + 1]
        next_omega = psd_matrix(
            model,
            domega,
            "hybrid_omega$(depth)",
            boundary_charges,
        )
        add_scaled_symmetric_equalities!(
            model,
            ptrace_right(next_omega, dw, 4) -
            propagation_scale *
            right_map * omega[index] * right_map',
            scale=recursion_scale,
        )
        add_scaled_symmetric_equalities!(
            model,
            ptrace_left(next_omega, 4, dw) -
            propagation_scale *
            left_map * omega[index] * left_map',
            scale=recursion_scale,
        )
        omega[index + 1] = next_omega
    end
    (
        rho4=rho4,
        omega=omega,
        omega_normalization=omega_normalization,
        omega_trace_scales=omega_trace_scales,
        local_link_scale=local_link_scale,
        omega4_link_scale=omega4_link_scale,
        recursion_scale=recursion_scale,
        rho4_link_equality_count=rho4_link_count,
        omega4_link_equality_count=omega4_link_count,
    )
end

function dualized_primal_value(dual_model, primal_variable)
    constraint, index =
        Dualization._get_dual_constraint(dual_model, primal_variable)
    constraint isa ConstraintRef ||
        error("primal variable is not represented by a dual constraint")
    multiplier = JuMP.dual(constraint)
    value = multiplier isa Number ? multiplier : multiplier[index]
    -Float64(value)
end

function structured_moment_values_from_sohs(structured)
    coefficient_matching =
        only(all_constraints(structured.sohs_model, Vector{AffExpr}, MOI.Zeros))
    multipliers = JuMP.dual(coefficient_matching)
    length(multipliers) == length(structured.physical_moments) ||
        error("SOHS equality multipliers do not match the structured moments")
    scaled_moments = Float64.(multipliers)
    Dict(
        variable => scaled_moments[index]
        for (index, variable) in
            enumerate(Dualization._get_dual_variables(structured.model, coefficient_matching))
    )
end

function evaluate_affine(expression::AffExpr, values)
    Float64(expression.constant) + sum(
        Float64(coefficient) * values[variable]
        for (coefficient, variable) in linear_terms(expression);
        init=0.0,
    )
end
