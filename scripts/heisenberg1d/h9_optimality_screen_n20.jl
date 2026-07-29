# H9 — screen local ground-state optimality constraints at the active finite N=20 depth.

include(joinpath(@__DIR__, "h2_psd_optimality.jl"))

const H9_N = 20
const H9_N_SUPER = 10
const H9_DYADIC_BITS = 13
const H9_MIN_IMPROVEMENT = 2e-7
const H9_DEFAULT_MOSEK_THREADS = 16

function parse_h9_args(args)
    run_name = nothing
    mosek_threads = H9_DEFAULT_MOSEK_THREADS
    reference_result = nothing
    reference_tensor = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--reference-result=")
            spec = split(arg, "=", limit=2)[2]
            reference_result = normpath(isabspath(spec) ? spec : joinpath(BOUNDGSENERGY_ROOT, spec))
        elseif startswith(arg, "--reference-tensor=")
            spec = split(arg, "=", limit=2)[2]
            reference_tensor = normpath(isabspath(spec) ? spec : joinpath(BOUNDGSENERGY_ROOT, spec))
        else
            error("unknown argument: $arg")
        end
    end
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h9-optimality-screen-n20")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(reference_result) && error("--reference-result is required")
    isnothing(reference_tensor) && error("--reference-tensor is required")
    isfile(reference_result) || error("reference result does not exist: $reference_result")
    isfile(reference_tensor) || error("reference tensor does not exist: $reference_tensor")
    run_name, mosek_threads, reference_result, reference_tensor
end

function plot_h9(results, figure_path)
    labels = ["stationarity", "PSD optimality"]
    points = [results["solves"][label] for label in labels]
    improvements = [point["raw_improvement"] for point in points]
    residuals = [point["combined_solver_residual_scale"] for point in points]
    threshold = results["gates"]["minimum_improvement"]
    scales = vcat(abs.(improvements), residuals, threshold)
    delta_plot = scatter(
        1:2,
        abs.(improvements);
        yscale=:log10,
        xticks=(1:2, labels),
        label="|raw objective Δ|",
        ylabel="energy-density scale",
        title="N=20 active-depth screen",
        markersize=8,
        ylims=(minimum(scales) / 2, maximum(scales) * 2),
        left_margin=6mm,
        bottom_margin=7mm,
    )
    scatter!(delta_plot, 1:2, residuals; marker=:diamond, markersize=8, label="sum residual scales")
    hline!(delta_plot, [threshold]; linestyle=:dash, color=:black, label="threshold")
    walls = [point["solve"]["solver_wall_seconds"] for point in points]
    wall_plot = bar(
        labels,
        walls;
        ylabel="Mosek wall time (s)",
        title="Screening cost",
        legend=false,
        color=[:steelblue, :darkorange],
        ylims=(0, 1.08 * maximum(walls)),
        left_margin=5mm,
        bottom_margin=7mm,
    )
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1100, 460)), figure_path)
end

function h9_main(args)
    run_name, mosek_threads, reference_result_path, reference_tensor_path = parse_h9_args(args)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_path = joinpath(result_dir, "optimality_screen_n20.json")
    figure_path = joinpath(result_dir, "figs", "optimality_screen_n20.png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    reference_result = JSON.parsefile(reference_result_path)
    reference_result["gates"]["stage_passed"] || error("H9 finite H7 reference did not pass")
    reference = reference_result["comparison"]["hybrid"]
    B, tensor_metadata = deserialize(reference_tensor_path)
    so_tensor_fingerprint(B) == reference_result["meta"]["coarse_grainer_fingerprint"] ||
        error("H9 reference tensor fingerprint is stale")
    tensor_metadata["dyadic_bits"] == H9_DYADIC_BITS || error("H9 reference uses the wrong dyadic grid")

    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    stationarity_labels, stationarity_family = central_stationarity_family()
    retained, stationarity_matrices = projected_family(stationarity_family, rho3_charges(super_charges))
    basis_labels, basis = traceless_super_site_basis()
    coefficients = optimality_coefficients(basis)

    solve_specs = [
        ("stationarity", () -> solve_formulation(
            B,
            H9_N_SUPER,
            virtual_charges,
            super_charges;
            blocked=true,
            rdm4=true,
            stationarity_matrices,
            mosek_threads,
        )),
        ("PSD optimality", () -> solve_formulation(
            B,
            H9_N_SUPER,
            virtual_charges,
            super_charges;
            blocked=true,
            rdm4=true,
            stationarity_matrices,
            optimality_coefficients=coefficients,
            mosek_threads,
        )),
    ]
    solves = Dict{String, Any}()
    for (label, run_solve) in solve_specs
        GC.gc()
        solve = run_solve()
        validated_solver_bounds(solve, "H9 $label")
        raw_improvement = solve["energy"] - reference["energy"]
        combined_residual = solver_residual_scale(solve) + solver_residual_scale(reference)
        solves[label] = Dict(
            "solve" => solve,
            "raw_improvement" => raw_improvement,
            "combined_solver_residual_scale" => combined_residual,
            "residual_resolved" => raw_improvement > max(H9_MIN_IMPROVEMENT, combined_residual),
        )
        @printf("%s raw improvement %.3e, combined residual %.3e\n", label, raw_improvement, combined_residual)
        flush(stdout)
    end
    candidate_found = any(point["residual_resolved"] for point in values(solves))

    results = Dict(
        "meta" => Dict(
            "stage" => "H9",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H9_N,
            "boundary" => "periodic",
            "support_sites" => 2 * H9_N_SUPER,
            "n_super" => H9_N_SUPER,
            "bond_dimension" => H1_D,
            "dyadic_bits" => H9_DYADIC_BITS,
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => so_tensor_fingerprint(B),
            "inputs" => Dict(
                "h7_finite_result" => relpath(reference_result_path, BOUNDGSENERGY_ROOT),
                "h7_tensor" => relpath(reference_tensor_path, BOUNDGSENERGY_ROOT),
            ),
        ),
        "families" => Dict(
            "stationarity_labels" => stationarity_labels[retained],
            "stationarity_count" => length(stationarity_matrices),
            "optimality_basis_labels" => basis_labels,
            "optimality_matrix_dimension" => length(basis),
        ),
        "reference_rdm8" => reference,
        "solves" => solves,
        "gates" => Dict(
            "minimum_improvement" => H9_MIN_IMPROVEMENT,
            "numerical_candidate_found" => candidate_found,
            "strict_certificate_implemented" => false,
            "stage_passed" => false,
        ),
    )
    write_h1_json(result_path, results)
    plot_h9(results, figure_path)
    println("numerical candidate found=$candidate_found")
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h9_main(ARGS)
end
