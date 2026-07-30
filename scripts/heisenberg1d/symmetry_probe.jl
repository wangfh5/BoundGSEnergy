# H1 — exact U(1) block decomposition of the compressed-LTI SDP.
# A converged SU(2)-symmetric D=4 VUMPS state supplies a real two-site coarse-grainer.
# Its dense tensor is decomposed by the U(1) subgroup using integer charges 2Sz.

ENV["GKSwstype"] = "100"

using LinearAlgebra, Random, SHA, Printf, JSON, Plots, Dates
using Plots.PlotMeasures
using MPSKit, MPSKitModels, TensorKit

BLAS.set_num_threads(1)

isdefined(@__MODULE__, :BOUNDGSENERGY_ROOT) ||
    include(joinpath(@__DIR__, "t2t3_compressed_lti.jl"))

const H1_D = 4
const H1_DEFAULT_N_SUPER = 4
const H1_TOL = 1e-12
const H1_MAXITER = 500
const H1_FEAS_TOL = 1e-7
const H1_COARSE_GRAINER_PATH = normpath(joinpath(@__DIR__, "..", "..", "artifacts", "coarse_grainers", "heisenberg1d", "symmetry_probe_d4.jls"))

elapsed_seconds(t0::UInt64) = (time_ns() - t0) / 1e9

function parse_h1_args(args)
    n_super = H1_DEFAULT_N_SUPER
    mosek_threads = 1
    run_name = nothing
    for arg in args
        if startswith(arg, "--n-super=")
            n_super = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    n_super >= 4 || error("n_super must be at least 4")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h1-u1-block-d4-n$(n_super)")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    n_super, mosek_threads, run_name
end

function result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "symmetry_probe.json"), joinpath(result_dir, "figs", "block_probe.png")
end

function write_h1_json(path, value)
    migrate_numerical_bound_schema!(value)
    mkpath(dirname(path))
    mktemp(dirname(path)) do tmp, io
        JSON.print(io, value, 2)
        close(io)
        mv(tmp, path; force=true)
    end
end

function build_coarse_grainer()
    Random.seed!(20260728)
    H = heisenberg_XXX(Float64, SU2Irrep, InfiniteChain(2); spin=1 // 2)
    physical = Rep[SU₂](1 // 2 => 1)
    integer_bond = Rep[SU₂](0 => 1, 1 => 1)
    half_integer_bond = Rep[SU₂](1 // 2 => 2)
    psi0 = InfiniteMPS(Float64, [physical, physical], [integer_bond, half_integer_bond])
    t0 = time_ns()
    psi, envs, delta = find_groundstate(psi0, H, VUMPS(; tol=H1_TOL, maxiter=H1_MAXITER, verbosity=1))
    wall = elapsed_seconds(t0)
    bnorm = norm(delta)
    energy = real(sum(expectation_value(psi, H, envs)) / 2)
    @assert bnorm <= H1_TOL "VUMPS did not converge: ||B||=$bnorm"
    @assert energy >= E0 - 1e-12 "variational energy fell below the Bethe value"
    super_tensor(psi), Dict(
        "algorithm" => "VUMPS",
        "tensor_symmetry" => "SU(2)",
        "block_symmetry" => "U(1) subgroup",
        "bond_dimension" => H1_D,
        "energy_density" => energy,
        "bethe_gap" => energy - E0,
        "galerkin_norm" => bnorm,
        "tolerance" => H1_TOL,
        "maxiter" => H1_MAXITER,
        "wall_seconds" => wall,
        "converged" => true,
        "julia_version" => string(VERSION),
    )
end

function load_or_build_coarse_grainer()
    if isfile(H1_COARSE_GRAINER_PATH)
        try
            B, metadata = deserialize(H1_COARSE_GRAINER_PATH)
            get(metadata, "julia_version", nothing) == string(VERSION) ||
                error("cached coarse-grainer was built with another Julia version")
            return B, metadata, true
        catch error
            @warn "rebuilding incompatible coarse-grainer cache" exception=(error, catch_backtrace())
        end
    end
    B, metadata = build_coarse_grainer()
    mkpath(dirname(H1_COARSE_GRAINER_PATH))
    mktemp(dirname(H1_COARSE_GRAINER_PATH)) do tmp, io
        serialize(io, (B, metadata))
        close(io)
        mv(tmp, H1_COARSE_GRAINER_PATH; force=true)
    end
    B, metadata, false
end

function exact_charge_tensor(B, virtual_charges, super_charges)
    out = copy(B)
    largest_forbidden = 0.0
    for left in axes(out, 1), site in axes(out, 2), right in axes(out, 3)
        virtual_charges[left] + super_charges[site] == virtual_charges[right] && continue
        largest_forbidden = max(largest_forbidden, abs(out[left, site, right]))
        out[left, site, right] = 0
    end
    @assert largest_forbidden <= 1e-12 "coarse-grainer violates U(1) support by $largest_forbidden"
    out, largest_forbidden
end

function charge_blocks(charges)
    blocks = Dict{Int, Vector{Int}}()
    for (index, charge) in enumerate(charges)
        push!(get!(blocks, charge, Int[]), index)
    end
    sort(collect(blocks); by=first)
end

function block_inventory(charges)
    Dict(string(charge) => length(indices) for (charge, indices) in charge_blocks(charges))
end

function rho3_charges(super_charges)
    charges = zeros(Int, 64)
    for left in 1:4, middle in 1:4, right in 1:4
        index = (left - 1) * 16 + (middle - 1) * 4 + right
        charges[index] = super_charges[left] + super_charges[middle] + super_charges[right]
    end
    charges
end

function rho4_charges(super_charges)
    charges = zeros(Int, 256)
    for first in 1:4, second in 1:4, third in 1:4, fourth in 1:4
        index = (first - 1) * 64 + (second - 1) * 16 + (third - 1) * 4 + fourth
        charges[index] = super_charges[first] + super_charges[second] + super_charges[third] + super_charges[fourth]
    end
    charges
end

function rho5_charges(super_charges)
    charges = zeros(Int, 1024)
    for first in 1:4, second in 1:4, third in 1:4, fourth in 1:4, fifth in 1:4
        index = (first - 1) * 256 + (second - 1) * 64 + (third - 1) * 16 + (fourth - 1) * 4 + fifth
        charges[index] = super_charges[first] + super_charges[second] + super_charges[third] +
            super_charges[fourth] + super_charges[fifth]
    end
    charges
end

function verify_rdm4_marginal_identities(W2)
    local_dimension = 4
    full_dimension = local_dimension^4
    rng = MersenneTwister(20260728)
    factor = randn(rng, full_dimension, 5)
    seed = factor * factor'
    shift = zeros(full_dimension, full_dimension)
    for first in 1:4, second in 1:4, third in 1:4, fourth in 1:4
        source = ((first - 1) * 4 + second - 1) * 16 + (third - 1) * 4 + fourth
        target = ((second - 1) * 4 + third - 1) * 16 + (fourth - 1) * 4 + first
        shift[target, source] = 1
    end
    rho4 = zeros(full_dimension, full_dimension)
    shifted = seed
    for _ in 1:4
        rho4 .+= shifted
        shifted = shift * shifted * shift'
    end
    rho4 ./= tr(rho4)
    rho3_right = ptrace_right(rho4, 64, 4)
    rho3_left = ptrace_left(rho4, 4, 64)
    compression = kron(I4, W2, I4)
    omega4 = compression * rho4 * compression'
    right_map = kron(I4, W2)
    left_map = kron(W2, I4)
    residuals = Dict(
        "rho3_translation" => maximum(abs, rho3_right - rho3_left),
        "trace_right" => maximum(abs, ptrace_right(omega4, 4 * size(W2, 1), 4) - right_map * rho3_right * right_map'),
        "trace_left" => maximum(abs, ptrace_left(omega4, 4, 4 * size(W2, 1)) - left_map * rho3_left * left_map'),
    )
    maximum(values(residuals)) <= 1e-12 || error("the rho4 compression map violates a marginal identity")
    residuals
end

function validate_rdm8_map(B, virtual_charges, super_charges)
    W2 = W2map(B)
    compression_violation = largest_charge_violation(
        kron(I4, W2, I4),
        omega_charges(super_charges, virtual_charges),
        rho4_charges(super_charges),
    )
    compression_violation <= 1e-12 || error("the RDM8 compression map violates U(1) charge support")
    Dict(
        "compression_violation" => compression_violation,
        "marginal_identity_residuals" => verify_rdm4_marginal_identities(W2),
    )
end

function omega_charges(super_charges, virtual_charges)
    D = length(virtual_charges)
    charges = zeros(Int, 16 * D^2)
    for left_site in 1:4, left_bond in 1:D, right_bond in 1:D, right_site in 1:4
        index = (((left_site - 1) * D + left_bond - 1) * D + right_bond - 1) * 4 + right_site
        charges[index] = super_charges[left_site] - virtual_charges[left_bond] + virtual_charges[right_bond] + super_charges[right_site]
    end
    charges
end

function boundary_charges_right(super_charges, virtual_charges)
    D = length(virtual_charges)
    charges = zeros(Int, 4 * D^2)
    for site in 1:4, left_bond in 1:D, right_bond in 1:D
        index = (site - 1) * D^2 + (left_bond - 1) * D + right_bond
        charges[index] = super_charges[site] - virtual_charges[left_bond] + virtual_charges[right_bond]
    end
    charges
end

function boundary_charges_left(super_charges, virtual_charges)
    D = length(virtual_charges)
    charges = zeros(Int, 4 * D^2)
    for left_bond in 1:D, right_bond in 1:D, site in 1:4
        index = ((left_bond - 1) * D + right_bond - 1) * 4 + site
        charges[index] = -virtual_charges[left_bond] + virtual_charges[right_bond] + super_charges[site]
    end
    charges
end

function largest_charge_violation(matrix, row_charges, column_charges)
    @assert size(matrix) == (length(row_charges), length(column_charges))
    largest = 0.0
    for row in axes(matrix, 1), column in axes(matrix, 2)
        row_charges[row] == column_charges[column] && continue
        largest = max(largest, abs(matrix[row, column]))
    end
    largest
end

function assert_charge_preserving_maps(B, virtual_charges, super_charges)
    D = length(virtual_charges)
    super_pair_charges = [super_charges[first] + super_charges[second] for first in 1:4 for second in 1:4]
    virtual_pair_charges = [-virtual_charges[left] + virtual_charges[right] for left in 1:D for right in 1:D]
    mr_input_charges = [virtual_charges[bond] + super_charges[site] for bond in 1:D for site in 1:4]
    ml_input_charges = [virtual_charges[bond] - super_charges[site] for site in 1:4 for bond in 1:D]
    rho_charges = rho3_charges(super_charges)
    omega_charge_labels = omega_charges(super_charges, virtual_charges)
    right_boundary_charges = boundary_charges_right(super_charges, virtual_charges)
    left_boundary_charges = boundary_charges_left(super_charges, virtual_charges)
    W2 = W2map(B)
    MR = MRmap(B)
    ML = MLmap(B)
    Iw = Matrix{Float64}(I, 4 * D, 4 * D)
    violations = Dict(
        "W2" => largest_charge_violation(W2, virtual_pair_charges, super_pair_charges),
        "MR" => largest_charge_violation(MR, virtual_charges, mr_input_charges),
        "ML" => largest_charge_violation(ML, virtual_charges, ml_input_charges),
        "IW2_right" => largest_charge_violation(kron(I4, W2), right_boundary_charges, rho_charges),
        "IW2_left" => largest_charge_violation(kron(W2, I4), left_boundary_charges, rho_charges),
        "I_MR" => largest_charge_violation(kron(Iw, MR), right_boundary_charges, omega_charge_labels),
        "I_ML" => largest_charge_violation(kron(ML, Iw), left_boundary_charges, omega_charge_labels),
        "HLOC" => largest_charge_violation(HLOC, super_pair_charges, super_pair_charges),
    )
    @assert maximum(values(violations)) <= 1e-12 "a compressed map violates U(1) charge support"
    violations
end

function block_psd_matrix(model, charges, name)
    matrix = [AffExpr(0.0) for _row in eachindex(charges), _column in eachindex(charges)]
    for (charge, indices) in charge_blocks(charges)
        dimension = length(indices)
        block = @variable(model, [1:dimension, 1:dimension], PSD, base_name="$(name)_q$(charge)")
        for (block_row, row) in enumerate(indices), (block_column, column) in enumerate(indices)
            add_to_expression!(matrix[row, column], block[block_row, block_column])
        end
    end
    matrix
end

function psd_matrix(model, dimension, name, charges)
    if isnothing(charges)
        return @variable(model, [1:dimension, 1:dimension], PSD, base_name=name)
    end
    block_psd_matrix(model, charges, name)
end

function build_formulation_model(B, n, virtual_charges, super_charges; blocked, rdm4=false, rdm5=false, stationarity_matrices=(), optimality_coefficients=nothing, symmetry_generators=(), optimizer=Mosek.Optimizer, mosek_threads=1, feasibility_tolerance=H1_FEAS_TOL, silent=true)
    rdm5 && !rdm4 && error("the five-super-site RDM extension requires rdm4=true")
    D = size(B, 1)
    dw = 4 * D^2
    domega = 16 * D^2
    rho_charges = blocked ? rho3_charges(super_charges) : nothing
    boundary_charges = blocked ? omega_charges(super_charges, virtual_charges) : nothing
    W2 = W2map(B)
    MR = MRmap(B)
    ML = MLmap(B)
    Iw = Matrix{Float64}(I, 4 * D, 4 * D)

    model = isnothing(optimizer) ? Model() : Model(optimizer)
    if !isnothing(optimizer)
        silent && set_silent(model)
        set_attribute(model, "MSK_IPAR_NUM_THREADS", mosek_threads)
        set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", feasibility_tolerance)
        set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", feasibility_tolerance)
        set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", feasibility_tolerance)
    end

    rho3 = psd_matrix(model, 64, "rho3", rho_charges)
    @constraint(model, tr(rho3) == 1)
    rho2_left = ptrace_left(rho3, 4, 16)
    rho2_right = ptrace_right(rho3, 16, 4)
    @constraint(model, rho2_left .== rho2_right)
    for matrix in stationarity_matrices
        @assert size(matrix) == (64, 64)
        @constraint(model, tr(matrix * rho3) == 0)
    end
    optimality_dimension = 0
    if !isnothing(optimality_coefficients)
        optimality_dimension = size(optimality_coefficients, 1)
        @assert size(optimality_coefficients) == (optimality_dimension, optimality_dimension)
        optimality_matrix = @variable(model, [1:optimality_dimension, 1:optimality_dimension], PSD, base_name="ground_state_optimality")
        for row in 1:optimality_dimension, column in 1:optimality_dimension
            @assert size(optimality_coefficients[row, column]) == (64, 64)
            @constraint(model, optimality_matrix[row, column] == tr(optimality_coefficients[row, column] * rho3))
        end
    end
    for generator in symmetry_generators
        @assert size(generator) == (64, 64)
        @assert norm(generator - generator') <= 1e-12
        commutator = rho3 * generator - generator * rho3
        for row in 1:63, column in row + 1:64
            @constraint(model, commutator[row, column] == 0)
        end
    end
    @objective(model, Min, tr(HLOC * rho2_left) / 2)

    omega = Vector{Any}(undef, n - 3)
    IW2_right = kron(I4, W2)
    IW2_left = kron(W2, I4)
    omega[1] = psd_matrix(model, domega, "omega4", boundary_charges)
    @constraint(model, ptrace_right(omega[1], dw, 4) .== IW2_right * rho3 * IW2_right')
    @constraint(model, ptrace_left(omega[1], 4, dw) .== IW2_left * rho3 * IW2_left')
    if rdm4
        rdm_charges = blocked ? rho4_charges(super_charges) : nothing
        rho4 = psd_matrix(model, 256, "rho4", rdm_charges)
        @constraint(model, ptrace_right(rho4, 64, 4) .== rho3)
        @constraint(model, ptrace_left(rho4, 4, 64) .== rho3)
        compression = kron(I4, W2, I4)
        @constraint(model, omega[1] .== compression * rho4 * compression')
        if rdm5
            rdm5_block_charges = blocked ? rho5_charges(super_charges) : nothing
            rho5 = psd_matrix(model, 1024, "rho5", rdm5_block_charges)
            @constraint(model, ptrace_right(rho5, 256, 4) .== rho4)
            @constraint(model, ptrace_left(rho5, 4, 256) .== rho4)
        end
    end

    I_MR = kron(Iw, MR)
    I_ML = kron(ML, Iw)
    for (index, hierarchy_depth) in enumerate(5:n)
        omega[index + 1] = psd_matrix(model, domega, "omega$(hierarchy_depth)", boundary_charges)
        @constraint(model, ptrace_right(omega[index + 1], dw, 4) .== I_MR * omega[index] * I_MR')
        @constraint(model, ptrace_left(omega[index + 1], 4, dw) .== I_ML * omega[index] * I_ML')
    end
    model, optimality_dimension
end

function solve_formulation(B, n, virtual_charges, super_charges; blocked, rdm4=false, rdm5=false, stationarity_matrices=(), optimality_coefficients=nothing, symmetry_generators=(), mosek_threads=1, feasibility_tolerance=H1_FEAS_TOL, strict_certificate=false, dual_transcript_path=nothing, dyadic_bits=nothing, silent=true)
    strict_certificate && (!isempty(stationarity_matrices) || !isnothing(optimality_coefficients) || !isempty(symmetry_generators)) &&
        error("strict replayable certification currently supports the compressed and RDM8 formulations")
    strict_certificate && isnothing(dual_transcript_path) &&
        error("strict replayable certification requires dual_transcript_path")

    total_start = time_ns()
    model, optimality_dimension = build_formulation_model(
        B,
        n,
        virtual_charges,
        super_charges;
        blocked,
        rdm4,
        rdm5,
        stationarity_matrices,
        optimality_coefficients,
        symmetry_generators,
        mosek_threads,
        feasibility_tolerance,
        silent,
    )
    build_wall = elapsed_seconds(total_start)

    solve_start = time_ns()
    optimize!(model)
    optimize_wall = elapsed_seconds(solve_start)
    total_wall = elapsed_seconds(total_start)
    status = termination_status(model)
    suffix = string(rdm4 ? "_rdm4" : "", rdm5 ? "_rdm5" : "", isempty(stationarity_matrices) ? "" : "_stationarity", isnothing(optimality_coefficients) ? "" : "_psd_optimality", isempty(symmetry_generators) ? "" : "_su2")
    record = Dict{String, Any}(
        "formulation" => string(blocked ? "u1_blocks" : "dense", suffix),
        "rdm4" => rdm4,
        "rdm5" => rdm5,
        "stationarity_constraint_count" => length(stationarity_matrices),
        "optimality_matrix_dimension" => optimality_dimension,
        "symmetry_generator_count" => length(symmetry_generators),
        "status" => string(status),
        "build_wall_seconds" => build_wall,
        "optimize_wall_seconds" => optimize_wall,
        "total_wall_seconds" => total_wall,
        "solver_wall_seconds" => solve_time(model),
        "scalar_variables" => num_variables(model),
        "mosek_threads" => mosek_threads,
        "solver_feasibility_tolerance" => feasibility_tolerance,
    )
    if status == OPTIMAL
        record["energy"] = objective_value(model)
        record["gap_to_bethe"] = E0 - record["energy"]
        record["residuals"] = mosek_residuals(model)
        record["residual_adjusted_lower_bound"] = residual_adjusted_lower_bound(record)
        record["published_lower_bound"] = published_lower_bound(record)
        record["residual_adjustment_is_rigorous_certificate"] = false
        if strict_certificate
            certificate_start = time_ns()
            record["dual_certificate"] = strict_dual_postcertificate(
                model,
                compressed_trace_bounds(B, n; rdm4, rdm5);
                transcript_path=dual_transcript_path,
                context=Dict(
                    "B" => B,
                    "n" => n,
                    "virtual_charges" => collect(virtual_charges),
                    "super_charges" => collect(super_charges),
                    "blocked" => blocked,
                    "rdm4" => rdm4,
                    "rdm5" => rdm5,
                    "dyadic_bits" => dyadic_bits,
                ),
            )
            record["postcertificate_wall_seconds"] = elapsed_seconds(certificate_start)
        end
    end
    record
end

function h2_objective_delta_plot(point, title, threshold)
    values = [
        abs(point["improvement"]),
        max(solver_residual_scale(point["baseline"]), solver_residual_scale(point["hybrid"])),
        threshold,
    ]
    scatter(
        1:3,
        values;
        yscale=:log10,
        xticks=(1:3, ["|raw objective Δ|", "solver residual scale", "acceptance threshold"]),
        ylabel="energy-density scale",
        title,
        label=false,
        markersize=8,
        markercolor=[:steelblue, :darkorange, :black],
        ylims=(minimum(values) / 2, maximum(values) * 2),
        left_margin=6mm,
        bottom_margin=8mm,
    )
end

function warm_up(mosek_threads)
    B = zeros(1, 4, 1)
    B[1, 2, 1] = 1
    dense = solve_formulation(B, 4, [0], [-2, 0, 0, 2]; blocked=false, mosek_threads)
    blocked = solve_formulation(B, 4, [0], [-2, 0, 0, 2]; blocked=true, mosek_threads)
    @assert dense["status"] == "OPTIMAL" && blocked["status"] == "OPTIMAL" "benchmark warm-up failed"
    Dict("dense_status" => dense["status"], "u1_status" => blocked["status"])
end

function plot_probe(results, figure_path, n_super)
    omega_inventory = results["block_inventory"]["omega_u1"]
    charges = sort(parse.(Int, collect(keys(omega_inventory))))
    dimensions = [omega_inventory[string(charge)] for charge in charges]
    size_plot = bar(string.(charges), dimensions; xlabel="total charge 2Sz", ylabel="PSD dimension", title="U(1) omega blocks (dense: 256)", legend=false, color=:steelblue, left_margin=5mm, bottom_margin=4mm)
    times = [results["solves"]["dense"]["total_wall_seconds"], results["solves"]["u1_blocks"]["total_wall_seconds"]]
    time_plot = bar(["dense", "U(1) blocks"], times; ylabel="total wall time (s)", title="D=4, n=$(n_super) after warm-up", legend=false, color=[:gray50, :darkorange], bottom_margin=4mm)
    mkpath(dirname(figure_path))
    savefig(plot(size_plot, time_plot; layout=(1, 2), size=(1000, 460)), figure_path)
end

function main(args)
    n_super, mosek_threads, run_name = parse_h1_args(args)
    result_dir, result_path, figure_path = result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")
    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_coarse_grainer()
    @assert size(B_raw) == (H1_D, 4, H1_D)
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    B_fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))

    rho_charges = rho3_charges(super_charges)
    boundary_charges = omega_charges(super_charges, virtual_charges)
    rho_inventory = block_inventory(rho_charges)
    omega_inventory = block_inventory(boundary_charges)
    largest_rho = maximum(values(rho_inventory))
    largest_omega = maximum(values(omega_inventory))
    @printf("rho3 blocks: %s (dense 64, largest %d)\n", rho_inventory, largest_rho)
    @printf("omega blocks: %s (dense 256, largest %d)\n", omega_inventory, largest_omega)
    flush(stdout)

    warmup = warm_up(mosek_threads)
    GC.gc()
    dense = solve_formulation(B, n_super, virtual_charges, super_charges; blocked=false, mosek_threads)
    GC.gc()
    blocked = solve_formulation(B, n_super, virtual_charges, super_charges; blocked=true, mosek_threads)

    dense_optimal = dense["status"] == "OPTIMAL"
    blocked_optimal = blocked["status"] == "OPTIMAL"
    @assert dense_optimal && blocked_optimal "both formulations must terminate OPTIMAL"
    energy_difference = abs(dense["energy"] - blocked["energy"])
    correctness = energy_difference <= 2e-7 && dense["residual_adjusted_lower_bound"] < E0 && blocked["residual_adjusted_lower_bound"] < E0
    block_reduction = 256 / largest_omega
    wall_nonregression = blocked["total_wall_seconds"] <= dense["total_wall_seconds"]
    @assert correctness "dense/U(1) energy comparison failed"
    @assert block_reduction >= 2 "largest PSD block did not shrink by a factor of two"

    results = Dict(
        "meta" => Dict(
            "stage" => "H1",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 2 * n_super,
            "n_super" => n_super,
            "bond_dimension" => H1_D,
            "energy_precision_digits" => 6,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "charge_units" => "2Sz",
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => B_fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
        ),
        "vumps" => vumps,
        "charge_check" => Dict("largest_forbidden_entry" => largest_forbidden, "exactified_to_zero" => true, "map_violations" => map_violations),
        "block_inventory" => Dict(
            "rho3_dense" => [64],
            "rho3_u1" => rho_inventory,
            "omega_dense" => [256],
            "omega_u1" => omega_inventory,
            "omega_variable_count" => n_super - 3,
            "largest_rho_reduction" => 64 / largest_rho,
            "largest_omega_reduction" => block_reduction,
        ),
        "warmup" => warmup,
        "solves" => Dict("dense" => dense, "u1_blocks" => blocked),
        "gates" => Dict(
            "both_optimal" => dense_optimal && blocked_optimal,
            "energy_difference" => energy_difference,
            "energy_agreement_atol" => 2e-7,
            "correctness_passed" => correctness,
            "largest_block_reduction" => block_reduction,
            "block_reduction_passed" => block_reduction >= 2,
            "wall_nonregression_passed" => wall_nonregression,
            "usefulness_passed" => block_reduction >= 2 && wall_nonregression,
        ),
    )
    write_h1_json(result_path, results)
    plot_probe(results, figure_path, n_super)
    @printf("dense E = %.9f, total wall %.2fs\n", dense["energy"], dense["total_wall_seconds"])
    @printf("U(1)  E = %.9f, total wall %.2fs\n", blocked["energy"], blocked["total_wall_seconds"])
    @printf("|Delta E| = %.3e, omega reduction %.2fx, wall non-regression %s\n", energy_difference, block_reduction, wall_nonregression)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
