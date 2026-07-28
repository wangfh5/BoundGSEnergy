# T2/T3 — Kull-Schuch compressed-LTI relaxation (PRX 14, 021008, Sec. II, Eq. 14)
# for the 1D Heisenberg AFM chain, using super-spins (2 physical -> 1 super, d̂=4)
# and the VUMPS uniform MPS as the coarse-graining map.
# Self-checks: E(n,D) <= e0 = 1/4 - ln2; non-decreasing in n; <= E_LTI(2n physical).

ENV["GKSwstype"] = "100"

using Serialization, LinearAlgebra, Printf, JuMP, MosekTools, Plots, JSON, Dates, SHA
using MPSKit, TensorKit   # needed only so deserialize() can resolve the saved types

const E0 = 1/4 - log(2)
const DEFAULT_DS = [2, 3]
const DEFAULT_RUN_SUFFIX = "kullschuch-mve"
const CLAIMED_DIGITS = 6
const COARSE_GRAINER_DIR = normpath(joinpath(@__DIR__, "..", "..", "artifacts", "coarse_grainers", "heisenberg1d"))

parse_int_list(spec::AbstractString) = [parse(Int, strip(x)) for x in split(spec, ",") if !isempty(strip(x))]

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
            "largest_physical_n" => 2 * completed[end][1],
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

function result_record(status, wall, fingerprint; E=nothing, gap=nothing, physical_n=nothing, super_n=nothing)
    rec = Dict{String, Any}(
        "status" => string(status),
        "wall" => wall,
        "coarse_grainer_fingerprint" => fingerprint,
    )
    isnothing(physical_n) || (rec["physical_n"] = physical_n)
    isnothing(super_n) || (rec["super_n"] = super_n)
    isnothing(E) || (rec["E"] = E)
    isnothing(gap) || (rec["gap"] = gap)
    rec
end

fmt_cert(x) = @sprintf("%.*f", CLAIMED_DIGITS, x)

function lti_record(status; E=nothing)
    rec = Dict{String, Any}("status" => string(status))
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
    D, d, _ = size(A1)
    B = zeros(D, d * d, D)
    for α in 1:D, s1 in 1:d, s2 in 1:d, γ in 1:D, β in 1:D
        B[α, (s1 - 1) * d + s2, γ] += A1[α, s1, β] * A2[β, s2, γ]
    end
    B
end

# W2 map: (βL,βR) <- (sa,sb), size (D^2, 16), MSF order everywhere
function W2map(B)
    D = size(B, 1)
    W = zeros(D * D, 16)
    for βL in 1:D, βR in 1:D, sa in 1:4, sb in 1:4, γ in 1:D
        W[(βL - 1) * D + βR, (sa - 1) * 4 + sb] += B[βL, sa, γ] * B[γ, sb, βR]
    end
    W
end
# M_R: β <- (β',s), size (D, D*4), MSF input order
function MRmap(B)
    D = size(B, 1)
    M = zeros(D, D * 4)
    for β in 1:D, βp in 1:D, s in 1:4
        M[β, (βp - 1) * 4 + s] = B[βp, s, β]
    end
    M
end
# M_L: β <- (s,β'), size (D, 4*D), MSF input order
function MLmap(B)
    D = size(B, 1)
    M = zeros(D, 4 * D)
    for β in 1:D, βp in 1:D, s in 1:4
        M[β, (s - 1) * D + βp] = B[β, s, βp]
    end
    M
end

# ---------- compressed LTI SDP ----------
function solve_compressed(B, n; silent=true)
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
    E, status
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
    p = plot(xscale=:log10, yscale=:log10, xlabel="physical spins N", ylabel="Delta E",
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
            println("Skipping cached LTI physical n=$n"); flush(stdout)
            continue
        end
        E, st = solve_LTI(n; silent=silent)
        if isnothing(E)
            lti[key] = lti_record(st)
            results["LTI"] = lti
            write_results(json_path, results)
            println("LTI physical n=$n: status $st -- no certificate, skipping")
            flush(stdout)
            continue
        end
        lti[key] = lti_record(st; E=E)
        @printf("LTI physical n=%2d: E = %s  (gap %.2e, status %s)\n", n, fmt_cert(E), E0 - E, st)
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
                                             physical_n=2 * n, super_n=n)
            else
                println("   WARNING: non-optimal status $st at D=$D n=$n -- no certificate, skipping")
                flush(stdout)
                results[key] = result_record(st, wall, fingerprint; physical_n=2 * n, super_n=n)
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
