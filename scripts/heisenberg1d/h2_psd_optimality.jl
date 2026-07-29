# H2 — positive-semidefinite local ground-state optimality matrix.

include(joinpath(@__DIR__, "h2_state_optimality.jl"))

const PSO_MOSEK_THREADS = 16
const PSO_MIN_IMPROVEMENT = 2e-7

function parse_pso_args(args)
    run_name = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h2-psd-optimality-d3")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name
end

function pso_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "psd_optimality.json"), joinpath(result_dir, "figs", "psd_optimality.png")
end

function traceless_super_site_basis()
    labels = String[]
    operators = Matrix{Float64}[]
    for index in 1:3
        operator = zeros(4, 4)
        operator[index, index] = inv(sqrt(2))
        operator[4, 4] = -inv(sqrt(2))
        push!(labels, "D$(index)4")
        push!(operators, operator)
    end
    for row in 1:3, column in row + 1:4
        symmetric = zeros(4, 4)
        symmetric[row, column] = inv(sqrt(2))
        symmetric[column, row] = inv(sqrt(2))
        antisymmetric = zeros(4, 4)
        antisymmetric[row, column] = inv(sqrt(2))
        antisymmetric[column, row] = -inv(sqrt(2))
        push!(labels, "S$(row)$(column)")
        push!(operators, symmetric)
        push!(labels, "A$(row)$(column)")
        push!(operators, antisymmetric)
    end
    @assert rank(hcat(vec.(operators)...); atol=1e-10) == 15
    labels, operators
end

function optimality_coefficients(operators)
    touching_hamiltonian = kron(HB, I4) + kron(I4, HB) + kron(I4, h, I4)
    embedded = [kron(I4, operator, I4) for operator in operators]
    commutators = [touching_hamiltonian * operator - operator * touching_hamiltonian for operator in embedded]
    dimension = length(operators)
    coefficients = Matrix{Matrix{Float64}}(undef, dimension, dimension)
    for row in 1:dimension, column in 1:dimension
        coefficients[row, column] = embedded[row]' * commutators[column]
    end
    coefficients
end

function exact_optimality_check(coefficients)
    energy, density = periodic_six_site_ground_state()
    dimension = size(coefficients, 1)
    matrix = [tr(coefficients[row, column] * density) for row in 1:dimension, column in 1:dimension]
    symmetry_residual = maximum(abs, matrix - matrix')
    eigenvalues = eigvals(Symmetric((matrix + matrix') / 2))
    @assert symmetry_residual <= 1e-12 "the exact optimality matrix is not symmetric"
    @assert minimum(eigenvalues) >= -1e-12 "the exact optimality matrix is not positive semidefinite"
    energy, symmetry_residual, eigenvalues
end

function plot_psd_optimality(results, figure_path)
    point = results["point"]
    labels = ["compressed", "PSD optimality"]
    walls = [point["baseline"]["solver_wall_seconds"], point["hybrid"]["solver_wall_seconds"]]
    delta_plot = h2_objective_delta_plot(point, "PSD optimality effect", results["gates"]["minimum_improvement"])
    wall_plot = bar(labels, walls; ylabel="Mosek wall time (s)", title="PSD optimality cost", legend=false, ylims=(0, 1.08 * maximum(walls)), color=[:gray60, :darkorange], left_margin=5mm, bottom_margin=5mm)
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1050, 460)), figure_path)
end

function pso_main(args)
    run_name = parse_pso_args(args)
    result_dir, result_path, figure_path = pso_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    stationarity_labels, stationarity_matrices = central_stationarity_family()
    basis_labels, basis = traceless_super_site_basis()
    coefficients = optimality_coefficients(basis)
    exact_energy, symmetry_residual, exact_eigenvalues = exact_optimality_check(coefficients)

    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    B3 = super_tensor(psis[3])
    baseline = solve_formulation(B3, 4, Int[], Int[]; blocked=false, mosek_threads=PSO_MOSEK_THREADS)
    hybrid = solve_formulation(B3, 4, Int[], Int[]; blocked=false, stationarity_matrices, optimality_coefficients=coefficients, mosek_threads=PSO_MOSEK_THREADS)
    both_optimal = baseline["status"] == "OPTIMAL" && hybrid["status"] == "OPTIMAL"
    improvement = both_optimal ? hybrid["energy"] - baseline["energy"] : nothing
    passed = both_optimal && improvement > PSO_MIN_IMPROVEMENT && baseline["residual_adjusted_lower_bound"] < E0 && hybrid["residual_adjusted_lower_bound"] < E0

    results = Dict(
        "meta" => Dict(
            "stage" => "H2",
            "run_name" => run_name,
            "strengthening" => "positive_semidefinite_ground_state_optimality",
            "contract" => "docs/PSD_STATE_OPTIMALITY_CONTRACT.md",
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 8,
            "optimality_support_sites" => 6,
            "n_super" => 4,
            "bond_dimension" => 3,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => PSO_MOSEK_THREADS,
            "coarse_grainer_fingerprint" => so_tensor_fingerprint(B3),
        ),
        "optimality_family" => Dict(
            "basis_labels" => basis_labels,
            "basis_dimension" => length(basis),
            "basis_rank" => rank(hcat(vec.(basis)...); atol=1e-10),
            "stationarity_labels" => stationarity_labels,
            "stationarity_count" => length(stationarity_matrices),
            "matrix_dimension" => size(coefficients, 1),
            "exact_periodic_n6_energy_per_site" => exact_energy,
            "exact_symmetry_residual" => symmetry_residual,
            "exact_minimum_eigenvalue" => minimum(exact_eigenvalues),
            "exact_maximum_eigenvalue" => maximum(exact_eigenvalues),
        ),
        "point" => Dict(
            "baseline" => baseline,
            "hybrid" => hybrid,
            "improvement" => improvement,
            "acceptance_passed" => passed,
        ),
        "gates" => Dict(
            "both_optimal" => both_optimal,
            "minimum_improvement" => PSO_MIN_IMPROVEMENT,
            "stage_passed" => passed,
        ),
    )
    write_h1_json(result_path, results)
    plot_psd_optimality(results, figure_path)
    @printf("published baseline %.6f, PSD optimality %.6f\n", baseline["energy"], hybrid["energy"])
    @printf("raw objective improvement %.3e, passed=%s\n", improvement, passed)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    pso_main(ARGS)
end
