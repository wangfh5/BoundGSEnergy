# H2 — linear ground-state stationarity on the central super-site.

include(joinpath(@__DIR__, "symmetry_probe.jl"))

const SO_MOSEK_THREADS = 16
const SO_MIN_IMPROVEMENT = 2e-7

so_tensor_fingerprint(B) = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
matrix_family_rank(matrices) = isempty(matrices) ? 0 : rank(hcat(vec.(matrices)...); atol=1e-10)

function parse_so_args(args)
    run_name = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h2-linear-stationarity-d3")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name
end

function so_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "state_optimality.json"), joinpath(result_dir, "figs", "state_optimality.png")
end

function central_stationarity_family()
    touching_hamiltonian = kron(HB, I4) + kron(I4, HB) + kron(I4, h, I4)
    labels = String[]
    matrices = Matrix{Float64}[]
    for row in 1:3, column in row + 1:4
        generator = zeros(4, 4)
        generator[row, column] = 1
        generator[column, row] = -1
        embedded_generator = kron(I4, generator, I4)
        stationarity = touching_hamiltonian * embedded_generator - embedded_generator * touching_hamiltonian
        @assert norm(stationarity - stationarity') <= 1e-12
        push!(labels, "A$(row)$(column)")
        push!(matrices, stationarity)
    end
    labels, matrices
end

function projected_family(matrices, charges)
    projected = Matrix{Float64}[]
    retained = Int[]
    for (index, matrix) in enumerate(matrices)
        block_matrix = copy(matrix)
        for row in axes(block_matrix, 1), column in axes(block_matrix, 2)
            charges[row] == charges[column] || (block_matrix[row, column] = 0)
        end
        norm(block_matrix) <= 1e-12 && continue
        push!(projected, block_matrix)
        push!(retained, index)
    end
    retained, projected
end

function embed_two_site(term, first_site, second_site, site_count)
    dimension = 2^site_count
    embedded = zeros(dimension, dimension)
    for source_zero in 0:dimension - 1
        bits = [((source_zero >> (site_count - site)) & 1) + 1 for site in 1:site_count]
        local_source = (bits[first_site] - 1) * 2 + bits[second_site]
        for first_output in 1:2, second_output in 1:2
            amplitude = term[(first_output - 1) * 2 + second_output, local_source]
            iszero(amplitude) && continue
            output_bits = copy(bits)
            output_bits[first_site] = first_output
            output_bits[second_site] = second_output
            target_zero = sum((output_bits[site] - 1) << (site_count - site) for site in 1:site_count)
            embedded[target_zero + 1, source_zero + 1] += amplitude
        end
    end
    embedded
end

function periodic_six_site_ground_state()
    site_count = 6
    hamiltonian = sum(embed_two_site(h, site, mod1(site + 1, site_count), site_count) for site in 1:site_count)
    decomposition = eigen(Symmetric(hamiltonian))
    ground_state = decomposition.vectors[:, 1]
    density = ground_state * ground_state'
    decomposition.values[1] / site_count, density
end

function exact_stationarity_check(matrices)
    energy, density = periodic_six_site_ground_state()
    residuals = [abs(tr(matrix * density)) for matrix in matrices]
    @assert maximum(residuals) <= 1e-12 "the exact periodic ground state violates stationarity"
    energy, residuals
end

function plot_state_optimality(results, figure_path)
    point = results["point"]
    labels = ["compressed", "stationarity"]
    walls = [point["baseline"]["solver_wall_seconds"], point["hybrid"]["solver_wall_seconds"]]
    delta_plot = h2_objective_delta_plot(point, "Linear stationarity effect", results["gates"]["minimum_improvement"])
    wall_plot = bar(labels, walls; ylabel="Mosek wall time (s)", title="Linear stationarity cost", legend=false, ylims=(0, 1.08 * maximum(walls)), color=[:gray60, :darkorange], left_margin=5mm, bottom_margin=5mm)
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1050, 460)), figure_path)
end

function so_main(args)
    run_name = parse_so_args(args)
    result_dir, result_path, figure_path = so_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    labels, matrices = central_stationarity_family()
    family_rank = matrix_family_rank(matrices)
    family_rank == length(matrices) || error("stationarity family is linearly dependent")
    exact_energy, exact_residuals = exact_stationarity_check(matrices)
    super_charges = [-2, 0, 0, 2]
    retained, u1_matrices = projected_family(matrices, rho3_charges(super_charges))

    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    B3 = super_tensor(psis[3])
    baseline = solve_formulation(B3, 4, Int[], Int[]; blocked=false, mosek_threads=SO_MOSEK_THREADS)
    hybrid = solve_formulation(B3, 4, Int[], Int[]; blocked=false, stationarity_matrices=matrices, mosek_threads=SO_MOSEK_THREADS)
    both_optimal = baseline["status"] == "OPTIMAL" && hybrid["status"] == "OPTIMAL"
    improvement = both_optimal ? hybrid["energy"] - baseline["energy"] : nothing
    passed = both_optimal && improvement > SO_MIN_IMPROVEMENT && baseline["residual_adjusted_lower_bound"] < E0 && hybrid["residual_adjusted_lower_bound"] < E0

    results = Dict(
        "meta" => Dict(
            "stage" => "H2",
            "run_name" => run_name,
            "strengthening" => "linear_ground_state_stationarity",
            "contract" => "docs/STATE_OPTIMALITY_CONTRACT.md",
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 8,
            "stationarity_support_sites" => 6,
            "n_super" => 4,
            "bond_dimension" => 3,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => SO_MOSEK_THREADS,
            "coarse_grainer_fingerprint" => so_tensor_fingerprint(B3),
        ),
        "stationarity_family" => Dict(
            "labels" => labels,
            "matrix_norms" => Dict(label => norm(matrix) for (label, matrix) in zip(labels, matrices)),
            "dense_count" => length(matrices),
            "dense_rank" => family_rank,
            "u1_retained_labels" => labels[retained],
            "u1_count" => length(u1_matrices),
            "u1_rank" => matrix_family_rank(u1_matrices),
            "exact_periodic_n6_energy_per_site" => exact_energy,
            "exact_periodic_n6_residuals" => Dict(label => residual for (label, residual) in zip(labels, exact_residuals)),
        ),
        "point" => Dict(
            "baseline" => baseline,
            "hybrid" => hybrid,
            "improvement" => improvement,
            "acceptance_passed" => passed,
        ),
        "gates" => Dict(
            "both_optimal" => both_optimal,
            "minimum_improvement" => SO_MIN_IMPROVEMENT,
            "stage_passed" => passed,
        ),
    )
    write_h1_json(result_path, results)
    plot_state_optimality(results, figure_path)
    @printf("published baseline %.6f, stationarity %.6f\n", baseline["energy"], hybrid["energy"])
    @printf("raw objective improvement %.3e, passed=%s\n", improvement, passed)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    so_main(ARGS)
end
