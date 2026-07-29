# T2/T3 — Kull-Schuch compressed-LTI relaxation (PRX 14, 021008, Sec. II, Eq. 14)
# for the 1D Heisenberg AFM chain, using super-spins (2 physical -> 1 super, d̂=4)
# and the VUMPS uniform MPS as the coarse-graining map.
# Self-checks: E(n,D) <= e0 = 1/4 - ln2; non-decreasing in n; <= E_LTI at matching physical-site support.

ENV["GKSwstype"] = "100"

using Serialization, LinearAlgebra, Printf, JuMP, MosekTools, Plots, JSON, Dates, SHA
using MPSKit, TensorKit   # needed only so deserialize() can resolve the saved types
import QMBCertify

const PINNED_JULIA_VERSION = v"1.11.5"
VERSION == PINNED_JULIA_VERSION ||
    error("BoundGSEnergy requires Julia $PINNED_JULIA_VERSION, found $VERSION")

const E0 = 1/4 - log(2)
const DEFAULT_DS = [2, 3]
const DEFAULT_RUN_SUFFIX = "kullschuch-mve"
const VUMPS_DEFAULT_SEED = 20260729
const CLAIMED_DIGITS = 6
const BOUNDGSENERGY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const COARSE_GRAINER_DIR = normpath(joinpath(@__DIR__, "..", "..", "artifacts", "coarse_grainers", "heisenberg1d"))
const RESULT_SCHEMA_VERSION = 2
const NUMERICAL_BOUND_SCHEMA_VERSION = 1
const SYSTEM_SCOPE = "infinite_translation_invariant"
const MOI = JuMP.MOI
const ExactRational = Rational{BigInt}

parse_int_list(spec::AbstractString) = [parse(Int, strip(x)) for x in split(spec, ",") if !isempty(strip(x))]

function bethe_ground_reference(system_size; tolerance=1e-13, max_iterations=100)
    iseven(system_size) || error("Bethe ground-state reference requires even N")
    down_spins = system_size ÷ 2
    quantum_numbers = collect(-(down_spins - 1) / 2:(down_spins - 1) / 2)
    rapidities = 0.5 .* tan.(pi .* quantum_numbers ./ system_size)

    function residual_jacobian(values)
        residual = zeros(down_spins)
        jacobian = zeros(down_spins, down_spins)
        for row in 1:down_spins
            residual[row] = system_size * 2 * atan(2 * values[row]) - 2 * pi * quantum_numbers[row]
            jacobian[row, row] = system_size * 4 / (1 + 4 * values[row]^2)
            for column in 1:down_spins
                row == column && continue
                difference = values[row] - values[column]
                phase_derivative = 2 / (1 + difference^2)
                residual[row] -= 2 * atan(difference)
                jacobian[row, row] -= phase_derivative
                jacobian[row, column] += phase_derivative
            end
        end
        residual, jacobian
    end

    for iteration in 1:max_iterations
        residual, jacobian = residual_jacobian(rapidities)
        max_residual = maximum(abs, residual)
        if max_residual <= tolerance
            total_energy = system_size / 4 - sum(0.5 / (rapidity^2 + 0.25) for rapidity in rapidities)
            return Dict(
                "method" => "finite periodic Bethe ansatz",
                "system_size" => system_size,
                "down_spins" => down_spins,
                "total_energy" => total_energy,
                "energy_per_site" => total_energy / system_size,
                "max_equation_residual" => max_residual,
                "newton_iterations" => iteration,
                "largest_absolute_rapidity" => maximum(abs, rapidities),
            )
        end
        step = jacobian \ residual
        current_norm = norm(residual)
        scale = 1.0
        while scale > 2.0^-20
            candidate = rapidities - scale * step
            candidate_residual, _ = residual_jacobian(candidate)
            if norm(candidate_residual) < current_norm
                rapidities = candidate
                break
            end
            scale /= 2
        end
        scale > 2.0^-20 || error("Bethe Newton line search failed")
    end
    error("Bethe ground-state reference did not converge")
end

function default_ns_by_D(D)
    if D <= 2
        [4, 6, 8, 10, 15, 20, 30, 50, 80, 120]
    elseif D == 3
        [4, 6, 8, 10, 15, 20, 30, 50]
    elseif D == 4
        [4, 6, 8, 10, 15, 20, 30, 40, 50]
    elseif D == 5
        [4, 6, 8, 10, 15, 20, 30]
    else
        [4, 6, 8, 10, 15, 20]
    end
end

function parse_args(args)
    Ds = copy(DEFAULT_DS)
    ns_override = nothing
    run_name = Dates.format(Dates.today(), "yyyymmdd") * "-" * DEFAULT_RUN_SUFFIX
    resume = true
    silent = true
    for arg in args
        if startswith(arg, "--Ds=")
            Ds = parse_int_list(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--ns=")
            ns_override = parse_int_list(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif arg == "--no-resume"
            resume = false
        elseif arg == "--debug"
            silent = false
        else
            error("unknown argument: $arg")
        end
    end
    sort!(unique!(Ds))
    ns_by_D = Dict(D => (isnothing(ns_override) ? default_ns_by_D(D) : copy(ns_override)) for D in Ds)
    Ds, ns_by_D, run_name, resume, silent
end

function ensure_run_state(rundir)
    mkpath(joinpath(rundir, "figs"))
    json_path = joinpath(rundir, "compressed_lti.json")
    if isfile(json_path)
        return JSON.parsefile(json_path)
    end
    Dict{String, Any}()
end

function write_results(json_path, results)
    mktemp(dirname(json_path)) do tmp, io
        JSON.print(io, results, 2)
        close(io)
        mv(tmp, json_path; force=true)
    end
end

function point_key(D, n)
    "D$(D)-n$(n)"
end

function take_optional!(entry, key)
    haskey(entry, key) || return nothing
    value = entry[key]
    delete!(entry, key)
    value
end

function migrate_result_schema!(results)
    for (key, entry) in results
        m = match(r"^D(\d+)-n(\d+)$", key)
        isnothing(m) && continue
        entry isa AbstractDict || error("cached point $key is not a result record")
        expected_n_super = parse(Int, m.captures[2])
        expected_support = 2 * expected_n_super
        legacy_support = take_optional!(entry, "physical_n")
        isnothing(legacy_support) || legacy_support == expected_support ||
            error("cached point $key has physical_n=$legacy_support, expected $expected_support")
        legacy_n_super = take_optional!(entry, "super_n")
        isnothing(legacy_n_super) || legacy_n_super == expected_n_super ||
            error("cached point $key has super_n=$legacy_n_super, expected $expected_n_super")
        support_sites = get!(entry, "support_sites", expected_support)
        support_sites == expected_support ||
            error("cached point $key has support_sites=$support_sites, expected $expected_support")
        n_super = get!(entry, "n_super", expected_n_super)
        n_super == expected_n_super ||
            error("cached point $key has n_super=$n_super, expected $expected_n_super")
        get(entry, "system_scope", SYSTEM_SCOPE) == SYSTEM_SCOPE ||
            error("cached point $key belongs to a different system scope")
        isnothing(get(entry, "system_size", nothing)) ||
            error("cached point $key has a finite system_size in an infinite-chain run")
        entry["system_scope"] = SYSTEM_SCOPE
        entry["system_size"] = nothing
    end

    lti = get(results, "LTI", Dict{String, Any}())
    for (key, entry) in lti
        entry isa AbstractDict || continue
        support_sites = parse(Int, key)
        get!(entry, "support_sites", support_sites) == support_sites ||
            error("cached LTI point $key has inconsistent support")
        entry["system_scope"] = SYSTEM_SCOPE
        entry["system_size"] = nothing
    end

    for entry in values(get(results, "plateaus", Dict{String, Any}()))
        entry isa AbstractDict || continue
        legacy_support = take_optional!(entry, "largest_physical_n")
        isnothing(legacy_support) || get!(entry, "largest_support_sites", legacy_support) == legacy_support ||
            error("cached plateau has inconsistent support")
    end

    meta = get!(results, "meta", Dict{String, Any}())
    meta["schema_version"] = RESULT_SCHEMA_VERSION
    meta["system_scope"] = SYSTEM_SCOPE
    meta["system_size"] = nothing
    results
end

function summarize_plateaus(results)
    by_D = Dict{Int, Vector{Tuple{Int, Float64}}}()
    for key in keys(results)
        m = match(r"^D(\d+)-n(\d+)$", key)
        isnothing(m) && continue
        entry = results[key]
        entry["status"] == "OPTIMAL" || continue
        haskey(entry, "gap") || continue
        D = parse(Int, m.captures[1])
        n = parse(Int, m.captures[2])
        push!(get!(by_D, D, Tuple{Int, Float64}[]), (n, entry["gap"]))
    end
    summary = Dict{String, Any}()
    for D in sort!(collect(keys(by_D)))
        completed = sort!(by_D[D]; by=first)
        summary["D$D"] = Dict(
            "largest_super_n" => completed[end][1],
            "largest_support_sites" => 2 * completed[end][1],
            "best_gap_seen" => minimum(last.(completed)),
            "last_gap" => completed[end][2],
        )
    end
    summary
end

function merge_vumps_summaries(old_summary, new_summary)
    merged_runs = Dict{String, Any}(String(k) => v for (k, v) in pairs(get(old_summary, "runs", Dict{String, Any}())))
    for (k, v) in pairs(get(new_summary, "runs", Dict{String, Any}()))
        merged_runs[String(k)] = v
    end
    Ds = sort!([parse(Int, only(match(r"^D(\d+)$", key).captures)) for key in keys(merged_runs)])
    Dict{String, Any}(
        "exact_energy" => get(new_summary, "exact_energy", get(old_summary, "exact_energy", E0)),
        "unit_cell" => get(new_summary, "unit_cell", get(old_summary, "unit_cell", 2)),
        "tol" => get(new_summary, "tol", get(old_summary, "tol", nothing)),
        "maxiter" => get(new_summary, "maxiter", get(old_summary, "maxiter", nothing)),
        "Ds_requested" => Ds,
        "runs" => merged_runs,
    )
end

function merge_meta!(results, Ds, ns_by_D, resume, silent, vumps_summary, coarse_fingerprints, run_name)
    meta = get!(results, "meta", Dict{String, Any}())
    old_Ds = [Int(x) for x in get(meta, "Ds", Any[])]
    meta["run_name"] = run_name
    meta["exact_energy"] = E0
    meta["schema_version"] = RESULT_SCHEMA_VERSION
    meta["system_scope"] = SYSTEM_SCOPE
    meta["system_size"] = nothing
    meta["Ds"] = sort!(unique!(vcat(old_Ds, Ds)))
    meta_ns = Dict{String, Any}(string(k) => v for (k, v) in pairs(get(meta, "ns_by_D", Dict{String, Any}())))
    for D in Ds
        meta_ns[string(D)] = ns_by_D[D]
    end
    for key in keys(results)
        m = match(r"^D(\d+)-n(\d+)$", key)
        isnothing(m) && continue
        D = m.captures[1]
        n = parse(Int, m.captures[2])
        ns = [Int(x) for x in get(meta_ns, D, Any[])]
        meta_ns[D] = sort!(unique!(push!(ns, n)))
    end
    meta["ns_by_D"] = meta_ns
    meta["resume"] = resume
    meta["silent"] = silent
    meta["vumps_summary"] = merge_vumps_summaries(get(meta, "vumps_summary", Dict{String, Any}()), vumps_summary)
    meta["coarse_grainer_fingerprints"] = Dict(string(D) => fp for (D, fp) in pairs(coarse_fingerprints))
    delete!(meta, "vumps_summary_path")
    results["meta"] = meta
end

function runs_compatible_for_backfill(old_summary, current_summary, D)
    old_runs = get(old_summary, "runs", Dict{String, Any}())
    current_runs = get(current_summary, "runs", Dict{String, Any}())
    key = "D$D"
    haskey(old_runs, key) || return false
    haskey(current_runs, key) || return false
    old_runs[key] == current_runs[key]
end

function backfill_point_fingerprints!(results, coarse_fingerprints, current_summary)
    old_summary = get(get(results, "meta", Dict{String, Any}()), "vumps_summary", Dict{String, Any}())
    for key in keys(results)
        m = match(r"^D(\d+)-n(\d+)$", key)
        isnothing(m) && continue
        D = parse(Int, m.captures[1])
        haskey(coarse_fingerprints, D) || continue
        entry = results[key]
        if !haskey(entry, "coarse_grainer_fingerprint")
            runs_compatible_for_backfill(old_summary, current_summary, D) ||
                error("cached point $key predates coarse-grainer fingerprints and its saved VUMPS summary does not match the current cache; choose a new --run-name or rerun with --no-resume")
            entry["coarse_grainer_fingerprint"] = coarse_fingerprints[D]
        end
    end
end

function vumps_summary_subset(summary, Ds)
    out = Dict{String, Any}(
        "exact_energy" => get(summary, "exact_energy", E0),
        "unit_cell" => get(summary, "unit_cell", 2),
        "tol" => get(summary, "tol", nothing),
        "maxiter" => get(summary, "maxiter", nothing),
        "Ds_requested" => Ds,
        "runs" => Dict{String, Any}(),
    )
    runs = get(summary, "runs", Dict{String, Any}())
    for D in Ds
        key = "D$D"
        haskey(runs, key) || continue
        out["runs"][key] = runs[key]
    end
    out
end

function coarse_grainer_fingerprint(psi)
    ctx = SHA.SHA1_CTX()
    for tensor in psi.AL
        arr = real.(convert(Array, tensor))
        update!(ctx, codeunits(join(string.(size(arr)), ",")))
        update!(ctx, reinterpret(UInt8, vec(arr)))
    end
    bytes2hex(digest!(ctx))
end

function validate_cached_point(results, key, fingerprint, resume)
    haskey(results, key) || return false
    !resume && return false
    cached = results[key]
    cached_fp = get(cached, "coarse_grainer_fingerprint", nothing)
    cached_fp == fingerprint || begin
        error("cached point $key uses a different coarse-grainer fingerprint; choose a new --run-name or rerun with --no-resume")
    end
    cached["status"] == "OPTIMAL" || return false
    haskey(cached, "E") || return false
    haskey(cached, "gap") || return false
    true
end

function result_record(status, wall, fingerprint; support_sites, n_super, E=nothing, gap=nothing)
    rec = Dict{String, Any}(
        "status" => string(status),
        "wall" => wall,
        "coarse_grainer_fingerprint" => fingerprint,
        "system_scope" => SYSTEM_SCOPE,
        "system_size" => nothing,
        "support_sites" => support_sites,
        "n_super" => n_super,
    )
    isnothing(E) || (rec["E"] = E)
    isnothing(gap) || (rec["gap"] = gap)
    rec
end

fmt_cert(x) = @sprintf("%.*f", CLAIMED_DIGITS, x)

function lti_record(status, support_sites; E=nothing)
    rec = Dict{String, Any}(
        "status" => string(status),
        "system_scope" => SYSTEM_SCOPE,
        "system_size" => nothing,
        "support_sites" => support_sites,
    )
    isnothing(E) || (rec["E"] = E)
    rec
end

function lti_entry_energy(entry)
    entry isa Real && return nothing
    get(entry, "status", nothing) == "OPTIMAL" || return nothing
    get(entry, "E", nothing)
end

function can_skip_lti(entry, resume)
    resume || return false
    entry isa Dict{String, Any} || return false
    get(entry, "status", nothing) == "OPTIMAL" || return false
    haskey(entry, "E")
end

# ---------- helpers ----------
_zero_container(rho::AbstractArray{<:Number}, n, m) = zeros(eltype(rho), n, m)
_zero_container(rho, n, m) = zeros(AffExpr, n, m)

function ptrace_left(rho, da, db)   # trace left factor of (da x db) system, kron/MSF index order
    r4 = reshape(rho, db, da, db, da)   # r4[j,i,l,k] == rho[(i,j),(k,l)], i left, j right
    out = _zero_container(rho, db, db)
    @inbounds for j in 1:db, l in 1:db, i in 1:da
        out[j, l] += r4[j, i, l, i]
    end
    out
end
function ptrace_right(rho, da, db)  # trace right factor
    r4 = reshape(rho, db, da, db, da)
    out = _zero_container(rho, da, da)
    @inbounds for i in 1:da, k in 1:da, j in 1:db
        out[i, k] += r4[j, i, j, k]
    end
    out
end

# spin-1/2 operators
const sx = [0 1; 1 0] ./ 2
const sy = [0 -im; im 0] ./ 2
const sz = [1 0; 0 -1] ./ 2
const h  = real.(kron(sx, sx) + kron(sy, sy) + kron(sz, sz))  # 4x4, S.S on two spins
const I2 = Matrix{Float64}(I, 2, 2)
const I4 = Matrix{Float64}(I, 4, 4)

# between-super-spin term: h on (phys2 of left, phys1 of right), 16x16, MSF order
function build_hb()
    hb = zeros(16, 16)
    for p1 in 1:2, p2 in 1:2, q1 in 1:2, q2 in 1:2, p2p in 1:2, q1p in 1:2
        p  = (p1 - 1) * 2 + p2;   q  = (q1 - 1) * 2 + q2
        pp = (p1 - 1) * 2 + p2p;  qp = (q1p - 1) * 2 + q2
        hb[(p - 1) * 4 + q, (pp - 1) * 4 + qp] = h[(p2 - 1) * 2 + q1, (p2p - 1) * 2 + q1p]
    end
    hb
end
const HB = build_hb()
const HLOC = 0.5 .* (kron(h, I4) .+ kron(I4, h)) .+ HB   # per-super-site, 16x16

# MPS -> super-tensor B (D, d̂=4, D)
function super_tensor(psi)
    A1 = convert(Array, psi.AL[1]); A2 = convert(Array, psi.AL[2])
    @assert maximum(abs, imag.(A1)) < 1e-10 && maximum(abs, imag.(A2)) < 1e-10 "complex MPS"
    A1 = real(A1); A2 = real(A2)
    Douter, d, Dinner = size(A1)
    Dinner2, d2, Douter2 = size(A2)
    @assert Douter == Douter2 && Dinner == Dinner2 && d == d2 "incompatible two-site MPS tensor dimensions"
    B = zeros(Douter, d * d, Douter)
    for α in 1:Douter, s1 in 1:d, s2 in 1:d, γ in 1:Douter, β in 1:Dinner
        B[α, (s1 - 1) * d + s2, γ] += A1[α, s1, β] * A2[β, s2, γ]
    end
    B
end

# W2 map: (βL,βR) <- (sa,sb), size (D^2, 16), MSF order everywhere
function W2map(B)
    D = size(B, 1)
    W = zeros(eltype(B), D * D, 16)
    for βL in 1:D, βR in 1:D, sa in 1:4, sb in 1:4, γ in 1:D
        W[(βL - 1) * D + βR, (sa - 1) * 4 + sb] += B[βL, sa, γ] * B[γ, sb, βR]
    end
    W
end
# M_R: β <- (β',s), size (D, D*4), MSF input order
function MRmap(B)
    D = size(B, 1)
    M = zeros(eltype(B), D, D * 4)
    for β in 1:D, βp in 1:D, s in 1:4
        M[β, (βp - 1) * 4 + s] = B[βp, s, β]
    end
    M
end
# M_L: β <- (s,β'), size (D, 4*D), MSF input order
function MLmap(B)
    D = size(B, 1)
    M = zeros(eltype(B), D, 4 * D)
    for β in 1:D, βp in 1:D, s in 1:4
        M[β, (s - 1) * D + βp] = B[β, s, βp]
    end
    M
end

function dyadic_coarse_grainer(B, bits)
    0 <= bits <= 26 || error("dyadic coarse-grainer bits must be between 0 and 26")
    scale = 2.0^bits
    round.(B .* scale) ./ scale
end

function fixed_scale_integers(values, scale)
    integers = Array{BigInt}(undef, size(values))
    for index in eachindex(values)
        scaled = exact_rational(values[index]) * scale
        denominator(scaled) == 1 || error("coarse-grainer entry is not on the declared dyadic grid")
        integers[index] = numerator(scaled)
    end
    integers
end

function dyadic_relaxation_audit(B, bits)
    scale = big(2)^bits
    integer_B = fixed_scale_integers(B, scale)
    exact_B = integer_B .// scale
    exact_W2 = W2map(exact_B)
    float_W2 = W2map(B)
    all(exact_rational(float_W2[index]) == exact_W2[index] for index in eachindex(float_W2)) ||
        error("Float64 W2 contraction is not exact")

    integer_W2 = fixed_scale_integers(exact_W2, scale^2)
    max_B_numerator = maximum(abs, integer_B)
    max_W2_numerator = maximum(abs, integer_W2)
    exact_integer_limit = big(2)^53
    numerator_bounds = Dict(
        "W2_contraction" => size(B, 1) * max_B_numerator^2,
        "MR_ML_quadratic_coefficient" => 2 * max_B_numerator^2,
        "W2_quadratic_coefficient" => 2 * max_W2_numerator^2,
    )
    maximum(values(numerator_bounds)) < exact_integer_limit ||
        error("dyadic SDP coefficient assembly exceeds the Float64 exact-integer range")

    Dict(
        "status" => "PROVEN",
        "dyadic_bits" => bits,
        "grid_denominator" => string(scale),
        "maximum_B_numerator" => string(max_B_numerator),
        "maximum_W2_numerator" => string(max_W2_numerator),
        "float64_exact_integer_limit" => string(exact_integer_limit),
        "numerator_bounds" => Dict(name => string(value) for (name, value) in numerator_bounds),
        "W2_exact_rational_contraction" => true,
        "quadratic_coefficient_assembly_exact" => true,
        "physical_relaxation_compatibility_proven" => true,
        "scope" => "compressed and RDM8 formulations",
    )
end

function mosek_residuals(model)
    task = unsafe_backend(model).task
    primal_objective = Mosek.getdouinf(task, Mosek.MSK_DINF_INTPNT_PRIMAL_OBJ)
    dual_objective = Mosek.getdouinf(task, Mosek.MSK_DINF_INTPNT_DUAL_OBJ)
    Dict(
        "primal_feasibility" => Mosek.getdouinf(task, Mosek.MSK_DINF_INTPNT_PRIMAL_FEAS),
        "dual_feasibility" => Mosek.getdouinf(task, Mosek.MSK_DINF_INTPNT_DUAL_FEAS),
        "primal_objective" => primal_objective,
        "dual_objective" => dual_objective,
        "relative_gap" => abs(primal_objective - dual_objective) / max(1.0, abs(primal_objective), abs(dual_objective)),
    )
end

exact_rational(value::Real) = rationalize(BigInt, value; tol=0)

function float_down(value::ExactRational)
    rounded = Float64(value)
    isfinite(rounded) || error("exact certificate does not fit in Float64")
    while exact_rational(rounded) > value
        rounded = prevfloat(rounded)
    end
    rounded
end

function float_up(value::ExactRational)
    rounded = Float64(value)
    isfinite(rounded) || error("exact certificate does not fit in Float64")
    while exact_rational(rounded) < value
        rounded = nextfloat(rounded)
    end
    rounded
end

function decimal_floor(value::ExactRational; digits=CLAIMED_DIGITS)
    scale = big(10)^digits
    floor(BigInt, value * scale) // scale
end

function induced_two_norm_squared_upper(matrix)
    row_bound = maximum(sum(exact_rational(abs(matrix[row, column])) for column in axes(matrix, 2)) for row in axes(matrix, 1))
    column_bound = maximum(sum(exact_rational(abs(matrix[row, column])) for row in axes(matrix, 1)) for column in axes(matrix, 2))
    row_bound * column_bound
end

function compressed_trace_bounds(B, n; rdm4=false, rdm5=false)
    bounds = Dict{String, ExactRational}("rho3" => 1 // 1)
    rdm4 && (bounds["rho4"] = 1 // 1)
    rdm5 && (bounds["rho5"] = 1 // 1)
    omega_bound = induced_two_norm_squared_upper(W2map(B))
    propagation_bound = min(
        induced_two_norm_squared_upper(MRmap(B)),
        induced_two_norm_squared_upper(MLmap(B)),
    )
    for depth in 4:n
        bounds["omega$(depth)"] = omega_bound
        omega_bound *= propagation_bound
    end
    bounds
end

function psd_group_name(constraint)
    object = constraint_object(constraint)
    variables = reshape_vector(object.func, object.shape)
    base = first(split(name(variables[1, 1]), "[", limit=2))
    replace(base, r"_q-?\d+$" => "")
end

function indexed_linear_terms(expression, variable_indices)
    sort!(
        [(variable_indices[variable], Float64(coefficient)) for (coefficient, variable) in linear_terms(expression)];
        by=first,
    )
end

function primal_transcript(model)
    variables = all_variables(model)
    variable_indices = Dict(variable => index for (index, variable) in enumerate(variables))
    objective = objective_function(model)
    equalities = [
        begin
            object = constraint_object(constraint)
            (
                rhs=Float64(object.set.value - object.func.constant),
                terms=indexed_linear_terms(object.func, variable_indices),
            )
        end
        for constraint in all_constraints(model, AffExpr, MOI.EqualTo{Float64})
    ]
    psd_blocks = NamedTuple[]
    covered_variables = Set{Int}()
    for constraint in all_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)
        object = constraint_object(constraint)
        matrix = reshape_vector(object.func, object.shape)
        indices = Int[]
        for row in axes(matrix, 1), column in row:last(axes(matrix, 2))
            index = variable_indices[matrix[row, column]]
            push!(indices, index)
            push!(covered_variables, index)
        end
        push!(psd_blocks, (
            group=psd_group_name(constraint),
            dimension=size(matrix, 1),
            variable_indices=indices,
        ))
    end
    covered_variables == Set(eachindex(variables)) ||
        error("strict dual post-certification requires every scalar variable to belong to a PSD block")
    (
        variable_count=length(variables),
        objective=(constant=Float64(objective.constant), terms=indexed_linear_terms(objective, variable_indices)),
        equalities=equalities,
        psd_blocks=psd_blocks,
    )
end

function reconstruct_dual_slacks(primal, multipliers)
    length(multipliers) == length(primal.equalities) ||
        error("dual transcript has the wrong number of equality multipliers")
    coefficients = fill(ExactRational(0), primal.variable_count)
    dual_constant = exact_rational(primal.objective.constant)
    for (index, coefficient) in primal.objective.terms
        coefficients[index] += exact_rational(coefficient)
    end
    for (equality, multiplier_value) in zip(primal.equalities, multipliers)
        multiplier = exact_rational(multiplier_value)
        dual_constant += multiplier * exact_rational(equality.rhs)
        for (index, coefficient) in equality.terms
            coefficients[index] -= multiplier * exact_rational(coefficient)
        end
    end
    grouped_slacks = Dict{String, Vector{Matrix{ExactRational}}}()
    covered_variables = Set{Int}()
    for block in primal.psd_blocks
        dimension = block.dimension
        slack = Matrix{ExactRational}(undef, dimension, dimension)
        cursor = 1
        for row in 1:dimension, column in row:dimension
            index = block.variable_indices[cursor]
            cursor += 1
            push!(covered_variables, index)
            entry = coefficients[index] / (row == column ? 1 : 2)
            slack[row, column] = entry
            slack[column, row] = entry
        end
        cursor == length(block.variable_indices) + 1 || error("dual transcript has a malformed PSD block")
        push!(get!(grouped_slacks, block.group, Matrix{ExactRational}[]), slack)
    end
    covered_variables == Set(1:primal.variable_count) ||
        error("dual transcript does not cover every scalar variable")
    dual_constant, grouped_slacks
end

function evaluate_dual_transcript(primal, multipliers, trace_bounds; eig_prec)
    dual_constant, grouped_slacks = reconstruct_dual_slacks(primal, multipliers)
    Set(keys(grouped_slacks)) == Set(keys(trace_bounds)) || error("PSD trace bounds do not match the assembled dual blocks")

    correction = ExactRational(0)
    groups = Dict{String, Any}()
    for group in sort!(collect(keys(grouped_slacks)))
        eigenvalue_bounds = [
            first(QMBCertify.rigorous_min_eig_bound(Symmetric(slack); prec=eig_prec))
            for slack in grouped_slacks[group]
        ]
        minimum_eigenvalue = minimum(eigenvalue_bounds)
        group_correction = min(ExactRational(0), minimum_eigenvalue) * trace_bounds[group]
        correction += group_correction
        groups[group] = Dict(
            "block_count" => length(eigenvalue_bounds),
            "minimum_eigenvalue_lower_bound" => float_down(minimum_eigenvalue),
            "minimum_eigenvalue_lower_bound_rational" => string(minimum_eigenvalue),
            "trace_upper_bound" => float_down(trace_bounds[group]),
            "trace_upper_bound_rational" => string(trace_bounds[group]),
            "correction" => float_down(group_correction),
            "correction_rational" => string(group_correction),
        )
    end
    dual_constant, correction, groups
end

function atomic_serialize(path, value)
    mkpath(dirname(path))
    mktemp(dirname(path)) do temporary_path, io
        serialize(io, value)
        close(io)
        mv(temporary_path, path; force=true)
    end
end

function artifact_relative_path(path)
    relative = relpath(abspath(path), BOUNDGSENERGY_ROOT)
    (relative == ".." || startswith(relative, "../")) &&
        error("certificate artifacts must be stored inside the repository")
    relative
end

file_sha256(path) = bytes2hex(sha256(read(path)))

function strict_dual_postcertificate(model, trace_bounds; transcript_path, context, eig_prec=256)
    objective_sense(model) == MOI.MIN_SENSE || error("strict dual post-certification requires a minimization problem")
    dyadic_bits = get(context, "dyadic_bits", nothing)
    compatibility_audit = isnothing(dyadic_bits) ? nothing : dyadic_relaxation_audit(context["B"], dyadic_bits)
    primal = primal_transcript(model)
    equality_constraints = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
    multipliers = Float64[JuMP.dual(constraint) for constraint in equality_constraints]
    transcript = Dict(
        "version" => 1,
        "context" => context,
        "primal" => primal,
        "equality_multipliers" => multipliers,
    )
    atomic_serialize(transcript_path, transcript)
    dual_constant, correction, groups = evaluate_dual_transcript(primal, multipliers, trace_bounds; eig_prec)
    certified = dual_constant + correction
    published = decimal_floor(certified)
    mosek_dual_objective = dual_objective_value(model)
    Dict(
        "status" => "CERTIFIED",
        "method" => "exact-rational dual reconstruction with rigorous Arblib eigenvalue enclosures and analytic PSD trace bounds",
        "coefficient_model" => isnothing(compatibility_audit) ?
            "the assembled Float64 JuMP coefficients treated as exact dyadic rationals" :
            "exactly assembled dyadic coarse-graining relaxation",
        "physical_relaxation_compatibility_proven" => !isnothing(compatibility_audit),
        "dyadic_compatibility_audit" => compatibility_audit,
        "eigenvalue_precision_bits" => eig_prec,
        "replay_transcript" => artifact_relative_path(transcript_path),
        "replay_transcript_sha256" => file_sha256(transcript_path),
        "coarse_grainer_fingerprint" => bytes2hex(sha1(reinterpret(UInt8, vec(context["B"])))),
        "dual_constant" => float_down(dual_constant),
        "dual_constant_rational" => string(dual_constant),
        "mosek_dual_objective" => mosek_dual_objective,
        "dual_objective_reconstruction_error" => abs(Float64(dual_constant) - mosek_dual_objective),
        "correction" => float_down(correction),
        "correction_rational" => string(correction),
        "certified_lower_bound" => float_down(certified),
        "certified_lower_bound_rational" => string(certified),
        "published_certified_lower_bound" => float_down(published),
        "published_certified_lower_bound_rational" => string(published),
        "groups" => groups,
    )
end

parse_exact_rational(value) = parse(ExactRational, value)

function replayed_dual_certificate(certificate, label)
    transcript_path = normpath(joinpath(BOUNDGSENERGY_ROOT, certificate["replay_transcript"]))
    isfile(transcript_path) || error("$label replay transcript does not exist")
    file_sha256(transcript_path) == certificate["replay_transcript_sha256"] ||
        error("$label replay transcript hash is stale")
    transcript = deserialize(transcript_path)
    transcript["version"] == 1 || error("$label uses an unsupported dual-transcript version")
    context = transcript["context"]
    B = context["B"]
    bytes2hex(sha1(reinterpret(UInt8, vec(B)))) == certificate["coarse_grainer_fingerprint"] ||
        error("$label coarse-grainer fingerprint is stale")
    exact_B, largest_forbidden = exact_charge_tensor(
        B,
        context["virtual_charges"],
        context["super_charges"],
    )
    largest_forbidden == 0 && exact_B == B ||
        error("$label coarse-grainer does not have exact U(1) support")
    all(iszero, values(assert_charge_preserving_maps(
        B,
        context["virtual_charges"],
        context["super_charges"],
    ))) || error("$label coarse-graining maps do not preserve charge exactly")
    dyadic_bits = get(context, "dyadic_bits", nothing)
    compatibility_audit = isnothing(dyadic_bits) ? nothing : dyadic_relaxation_audit(B, dyadic_bits)
    get(certificate, "physical_relaxation_compatibility_proven", false) == !isnothing(compatibility_audit) ||
        error("$label physical-relaxation compatibility status is stale")
    get(certificate, "dyadic_compatibility_audit", nothing) == compatibility_audit ||
        error("$label dyadic compatibility audit is stale")
    isempty(get(context, "stationarity_matrices", ())) || error("$label replay context contains unsupported stationarity constraints")
    isnothing(get(context, "optimality_coefficients", nothing)) || error("$label replay context contains unsupported optimality constraints")
    isempty(get(context, "symmetry_generators", ())) || error("$label replay context contains unsupported symmetry constraints")
    model, _ = build_formulation_model(
        B,
        context["n"],
        context["virtual_charges"],
        context["super_charges"];
        blocked=context["blocked"],
        rdm4=context["rdm4"],
        rdm5=get(context, "rdm5", false),
        optimizer=nothing,
    )
    primal_transcript(model) == transcript["primal"] ||
        error("$label replay transcript does not match an independently assembled primal SDP")
    trace_bounds = compressed_trace_bounds(
        B,
        context["n"];
        rdm4=context["rdm4"],
        rdm5=get(context, "rdm5", false),
    )
    dual_constant, correction, groups = evaluate_dual_transcript(
        transcript["primal"],
        transcript["equality_multipliers"],
        trace_bounds;
        eig_prec=certificate["eigenvalue_precision_bits"],
    )
    dual_constant, correction, groups
end

function validated_dual_certificate(record, label)
    certificate = get(record, "dual_certificate", nothing)
    certificate isa AbstractDict || error("$label has no strict dual certificate")
    certificate["status"] == "CERTIFIED" || error("$label dual certificate did not pass")
    dual_constant, correction, replayed_groups = replayed_dual_certificate(certificate, label)
    Set(keys(replayed_groups)) == Set(keys(certificate["groups"])) ||
        error("$label dual-certificate groups are stale")
    for (name, replayed) in replayed_groups
        stored = certificate["groups"][name]
        stored["block_count"] == replayed["block_count"] || error("$label $name block count is stale")
        for field in ("minimum_eigenvalue_lower_bound_rational", "trace_upper_bound_rational", "correction_rational")
            parse_exact_rational(stored[field]) == parse_exact_rational(replayed[field]) ||
                error("$label $name $field is stale")
        end
        for field in ("minimum_eigenvalue_lower_bound", "trace_upper_bound", "correction")
            stored[field] == replayed[field] || error("$label $name $field is stale")
        end
        expected_correction = min(
            ExactRational(0),
            parse_exact_rational(stored["minimum_eigenvalue_lower_bound_rational"]),
        ) * parse_exact_rational(stored["trace_upper_bound_rational"])
        expected_correction == parse_exact_rational(stored["correction_rational"]) ||
            error("$label $name correction formula is invalid")
    end
    dual_constant == parse_exact_rational(certificate["dual_constant_rational"]) ||
        error("$label dual constant is stale")
    correction == parse_exact_rational(certificate["correction_rational"]) || error("$label dual-certificate correction is stale")
    certified = dual_constant + correction
    certified == parse_exact_rational(certificate["certified_lower_bound_rational"]) || error("$label certified lower bound is stale")
    published = decimal_floor(certified)
    published == parse_exact_rational(certificate["published_certified_lower_bound_rational"]) || error("$label published certified lower bound is stale")
    certificate["certified_lower_bound"] == float_down(certified) || error("$label numeric certified lower bound is stale")
    certificate["published_certified_lower_bound"] == float_down(published) || error("$label numeric published certified lower bound is stale")
    certificate
end

# ---------- compressed LTI SDP ----------
function solve_compressed(B, n; silent=true, diagnostics=false)
    D = size(B, 1); dw = 4 * D * D            # size of Tr_phys(omega) output
    domega = 16 * D * D                       # omega matrix size
    W2 = W2map(B); MR = MRmap(B); ML = MLmap(B)
    Iw = Matrix{Float64}(I, 4 * D, 4 * D)     # identity on (s, other-bond) for R_A/L_A

    model = Model(Mosek.Optimizer); silent && set_silent(model)
    # residuals stall at ~2e-8 on this degenerate LTI chain; 1e-7 tolerance is
    # honest at our 1e-5 target precision (digits claimed << tolerance)
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-7)
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-7)
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-7)
    @variable(model, rho3[1:64, 1:64], PSD)
    @constraint(model, tr(rho3) == 1)
    rho2L = ptrace_left(rho3, 4, 16); rho2R = ptrace_right(rho3, 16, 4)
    @constraint(model, rho2L .== rho2R)
    @objective(model, Min, tr(HLOC * rho2L) / 2)

    omega = Vector{Any}(undef, n - 3)
    # link rho3 <-> omega4
    IW2r = kron(I4, W2)          # (4D^2 x 64), W2 on super-spins 2,3
    IW2l = kron(W2, I4)          # on super-spins 1,2
    omega[1] = @variable(model, [1:domega, 1:domega], PSD)
    TrR_o4 = ptrace_right(omega[1], dw, 4)
    TrL_o4 = ptrace_left(omega[1], 4, dw)
    @constraint(model, TrR_o4 .== IW2r * rho3 * IW2r')
    @constraint(model, TrL_o4 .== IW2l * rho3 * IW2l')
    # chain m = 5..n
    I_MR = kron(Iw, MR)          # (dw x domega): identity on (s,betaL), MR on (betaR',s)
    I_ML = kron(ML, Iw)          # identity on (betaR) hmm: (s,betaL',betaR) grouping
    for (k, m) in enumerate(5:n)
        omega[k + 1] = @variable(model, [1:domega, 1:domega], PSD)
        @constraint(model, ptrace_right(omega[k + 1], dw, 4) .== I_MR * omega[k] * I_MR')
        @constraint(model, ptrace_left(omega[k + 1], 4, dw) .== I_ML * omega[k] * I_ML')
    end
    optimize!(model)
    status = termination_status(model)
    E = status == OPTIMAL ? objective_value(model) : nothing
    diagnostics || return E, status
    record = Dict{String, Any}(
        "status" => string(status),
        "solver_wall_seconds" => solve_time(model),
    )
    if status == OPTIMAL
        record["energy"] = E
        record["residuals"] = mosek_residuals(model)
    end
    record
end

function solver_residual_scale(record)
    record["status"] == "OPTIMAL" || error("a residual-adjusted bound requires OPTIMAL status")
    residuals = record["residuals"]
    objective_gap = abs(residuals["primal_objective"] - residuals["dual_objective"])
    maximum((residuals["primal_feasibility"], residuals["dual_feasibility"], objective_gap))
end

function residual_adjusted_lower_bound(record)
    residuals = record["residuals"]
    minimum((record["energy"], residuals["primal_objective"], residuals["dual_objective"])) - solver_residual_scale(record)
end

function published_lower_bound(record; digits=CLAIMED_DIGITS)
    scale = 10.0^digits
    floor(residual_adjusted_lower_bound(record) * scale) / scale
end

function validated_solver_bounds(record, label)
    record["status"] == "OPTIMAL" || error("$label is not OPTIMAL")
    residuals = record["residuals"]
    all(isfinite, (
        record["energy"],
        residuals["primal_feasibility"],
        residuals["dual_feasibility"],
        residuals["primal_objective"],
        residuals["dual_objective"],
    )) || error("$label contains non-finite solver values")
    adjusted = residual_adjusted_lower_bound(record)
    published = published_lower_bound(record)
    record["residual_adjusted_lower_bound"] == adjusted || error("$label residual-adjusted bound is stale")
    record["published_lower_bound"] == published || error("$label published bound is stale")
    adjusted, published
end

function migrate_numerical_bound_records!(value)
    if value isa AbstractDict
        if get(value, "status", nothing) == "OPTIMAL" && haskey(value, "energy") && haskey(value, "residuals")
            value["residual_adjusted_lower_bound"] = residual_adjusted_lower_bound(value)
            value["published_lower_bound"] = published_lower_bound(value)
        end
        foreach(migrate_numerical_bound_records!, values(value))
    elseif value isa AbstractVector
        foreach(migrate_numerical_bound_records!, value)
    end
    value
end

function migrate_numerical_bound_schema!(results)
    migrate_numerical_bound_records!(results)
    get!(results, "meta", Dict{String, Any}())["numerical_bound_schema_version"] = NUMERICAL_BOUND_SCHEMA_VERSION
    results
end

# ---------- uncompressed LTI on the physical chain (d=2) ----------
function solve_LTI(n; silent=true)
    model = Model(Mosek.Optimizer); silent && set_silent(model)
    rho = Vector{Any}(undef, n - 1)
    rho[1] = @variable(model, [1:4, 1:4], PSD)
    @constraint(model, tr(rho[1]) == 1)
    for m in 3:n
        d2m = 2^m
        rho[m - 1] = @variable(model, [1:d2m, 1:d2m], PSD)
        @constraint(model, ptrace_left(rho[m - 1], 2, 2^(m - 1)) .== rho[m - 2])
        @constraint(model, ptrace_right(rho[m - 1], 2^(m - 1), 2) .== rho[m - 2])
    end
    # LTI at top level is implied by the chain; objective
    @objective(model, Min, tr(h * rho[1]))
    optimize!(model)
    status = termination_status(model)
    status == OPTIMAL ? objective_value(model) : nothing, status
end

# ---------- main ----------
function plot_results(rundir, results, Ds, lti)
    mkpath(joinpath(rundir, "figs"))
    p = plot(xscale=:log10, yscale=:log10, xlabel="constraint support (physical sites)", ylabel="Delta E",
             title="Compressed-LTI relaxation, Heisenberg chain",
             legend=:bottomleft, size=(800, 500))
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star6, :xcross]
    marker_sizes = [10, 7, 10, 10, 10, 10, 10]
    for (idx, D) in enumerate(Ds)
        ns = Int[]
        for key in keys(results)
            m = match(Regex("^D$(D)-n(\\d+)\$"), key)
            isnothing(m) && continue
            results[key]["status"] == "OPTIMAL" || continue
            push!(ns, parse(Int, m.captures[1]))
        end
        sort!(unique!(ns))
        isempty(ns) && continue
        ys = [results[point_key(D, n)]["gap"] for n in ns]
        plot!(p, 2 .* ns, ys, marker=markers[mod1(idx, length(markers))],
              markersize=marker_sizes[mod1(idx, length(marker_sizes))],
              markercolor=:white, markerstrokewidth=2, label="compressed, D=$D")
    end
    if !isempty(lti)
        lts = sort([parse(Int, key) for key in keys(lti) if !isnothing(lti_entry_energy(lti[key]))])
        isempty(lts) || plot!(p, lts, [E0 - lti_entry_energy(lti[string(n)]) for n in lts], marker=:hexagon, label="exact LTI")
    end
    savefig(p, joinpath(rundir, "figs", "gap_vs_n.png"))
end

function main(args)
    Ds, ns_by_D, run_name, resume, silent = parse_args(args)
    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    summary_path = joinpath(COARSE_GRAINER_DIR, "vumps_summary.json")
    vumps_summary = isfile(summary_path) ? JSON.parsefile(summary_path) : Dict{String, Any}()
    rundir = joinpath("results", run_name)
    mkpath(joinpath(rundir, "figs"))
    results = resume ? ensure_run_state(rundir) : Dict{String, Any}()
    migrate_result_schema!(results)
    json_path = joinpath(rundir, "compressed_lti.json")
    coarse_fingerprints = Dict{Int, String}()
    for D in keys(psis)
        coarse_fingerprints[D] = coarse_grainer_fingerprint(psis[D])
    end
    summary_subset = vumps_summary_subset(vumps_summary, Ds)
    backfill_point_fingerprints!(results, coarse_fingerprints, summary_subset)
    merge_meta!(results, Ds, ns_by_D, resume, silent, summary_subset, coarse_fingerprints, run_name)

    # helper sanity: ptrace on a product state
    a = [1.0, 0]; b = [0.0, 1, 0, 0]
    Pa = a * a'; Pb = b * b'
    @assert ptrace_left(kron(Pa, Pb), 2, 4) ≈ Pb && ptrace_right(kron(Pa, Pb), 2, 4) ≈ Pa
    println("helper sanity OK"); flush(stdout)

    # uncompressed LTI reference (physical chain), n = 4, 6 (n=8 is too heavy uncompressed)
    lti = get!(results, "LTI", Dict{String, Any}())
    for n in [4, 6]
        key = string(n)
        if haskey(lti, key) && can_skip_lti(lti[key], resume)
            println("Skipping cached LTI support n=$n"); flush(stdout)
            continue
        end
        E, st = solve_LTI(n; silent=silent)
        if isnothing(E)
            lti[key] = lti_record(st, n)
            results["LTI"] = lti
            write_results(json_path, results)
            println("LTI support n=$n: status $st -- no certificate, skipping")
            flush(stdout)
            continue
        end
        lti[key] = lti_record(st, n; E=E)
        @printf("LTI support n=%2d: E = %s  (gap %.2e, status %s)\n", n, fmt_cert(E), E0 - E, st)
        flush(stdout)
        @assert E <= E0 + 1e-9
        results["LTI"] = lti
        write_results(json_path, results)
    end

    # compressed ladder
    for D in Ds
        if !haskey(psis, D)
            println("Skipping D=$D: no converged VUMPS coarse-grainer in vumps_mps.jls"); flush(stdout)
            continue
        end
        B = super_tensor(psis[D])
        fingerprint = coarse_fingerprints[D]
        for n in ns_by_D[D]
            key = point_key(D, n)
            if validate_cached_point(results, key, fingerprint, resume)
                println("Skipping cached D=$D super-n=$n"); flush(stdout)
                continue
            end
            t0 = time()
            E, st = solve_compressed(B, n; silent=silent)
            wall = time() - t0
            if !isnothing(E)
                gap = E0 - E
                @printf("D=%d super-n=%3d: E = %s  gap %.3e  wall %.1fs  status %s\n",
                        D, n, fmt_cert(E), gap, wall, st)
                flush(stdout)
                @assert E <= E0 + 1e-6 "bound ABOVE e0 at D=$D n=$n -- index bug"
                if n == 4 && haskey(lti, "6") && !isnothing(lti_entry_energy(lti["6"]))
                    @printf("   check: %s <= LTI(6) = %s ? %s\n",
                            fmt_cert(E), fmt_cert(lti_entry_energy(lti["6"])), E <= lti_entry_energy(lti["6"]) + 1e-9)
                    flush(stdout)
                end
                results[key] = result_record(st, wall, fingerprint; E=E, gap=gap,
                                             support_sites=2 * n, n_super=n)
            else
                println("   WARNING: non-optimal status $st at D=$D n=$n -- no certificate, skipping")
                flush(stdout)
                results[key] = result_record(st, wall, fingerprint; support_sites=2 * n, n_super=n)
            end
            results["plateaus"] = summarize_plateaus(results)
            write_results(json_path, results)
        end
    end
    results["plateaus"] = summarize_plateaus(results)
    write_results(json_path, results)
    plot_results(rundir, results, Ds, lti)
    println("DONE. figure + json saved to $rundir"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
    exit(0)
end
