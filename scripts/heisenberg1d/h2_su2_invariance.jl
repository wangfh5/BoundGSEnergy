# H2 — direct physical SU(2) invariance of the three-super-site RDM.

include(joinpath(@__DIR__, "h2_state_optimality.jl"))

const SU2_MOSEK_THREADS = 16
const SU2_MIN_IMPROVEMENT = 2e-7

function parse_su2_args(args)
    run_name = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h2-su2-invariance-d3")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name
end

function su2_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "su2_invariance.json"), joinpath(result_dir, "figs", "su2_invariance.png")
end

function embed_one_site(operator, site, site_count)
    factors = [position == site ? operator : I2 for position in 1:site_count]
    foldl(kron, factors)
end

function total_spin_generators()
    site_count = 6
    sum(embed_one_site(sx, site, site_count) for site in 1:site_count), sum(embed_one_site(sz, site, site_count) for site in 1:site_count)
end

function exact_su2_check(generators)
    energy, density = periodic_six_site_ground_state()
    residuals = [maximum(abs, density * generator - generator * density) for generator in generators]
    @assert maximum(residuals) <= 1e-12 "the exact periodic ground state is not SU(2) invariant"
    energy, residuals
end

function plot_su2_invariance(results, figure_path)
    point = results["point"]
    labels = ["compressed", "SU(2) invariant"]
    walls = [point["baseline"]["solver_wall_seconds"], point["hybrid"]["solver_wall_seconds"]]
    delta_plot = h2_objective_delta_plot(point, "Physical SU(2) effect", results["gates"]["minimum_improvement"])
    wall_plot = bar(labels, walls; ylabel="Mosek wall time (s)", title="Physical SU(2) cost", legend=false, ylims=(0, 1.08 * maximum(walls)), color=[:gray60, :darkorange], left_margin=5mm, bottom_margin=5mm)
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1050, 460)), figure_path)
end

function su2_main(args)
    run_name = parse_su2_args(args)
    result_dir, result_path, figure_path = su2_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    generators = collect(total_spin_generators())
    exact_energy, exact_residuals = exact_su2_check(generators)
    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    B3 = super_tensor(psis[3])
    baseline = solve_formulation(B3, 4, Int[], Int[]; blocked=false, mosek_threads=SU2_MOSEK_THREADS)
    hybrid = solve_formulation(B3, 4, Int[], Int[]; blocked=false, symmetry_generators=generators, mosek_threads=SU2_MOSEK_THREADS)
    both_optimal = baseline["status"] == "OPTIMAL" && hybrid["status"] == "OPTIMAL"
    improvement = both_optimal ? hybrid["energy"] - baseline["energy"] : nothing
    passed = both_optimal && improvement > SU2_MIN_IMPROVEMENT && baseline["residual_adjusted_lower_bound"] < E0 && hybrid["residual_adjusted_lower_bound"] < E0

    results = Dict(
        "meta" => Dict(
            "stage" => "H2",
            "run_name" => run_name,
            "strengthening" => "physical_su2_invariance",
            "contract" => "docs/SU2_INVARIANCE_CONTRACT.md",
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 8,
            "symmetry_support_sites" => 6,
            "n_super" => 4,
            "bond_dimension" => 3,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => SU2_MOSEK_THREADS,
            "coarse_grainer_fingerprint" => so_tensor_fingerprint(B3),
        ),
        "symmetry_family" => Dict(
            "generators" => ["total_Sx", "total_Sz"],
            "generator_count" => length(generators),
            "upper_triangle_equation_count" => length(generators) * 64 * 63 ÷ 2,
            "expected_independent_equation_count" => 2004,
            "invariant_real_symmetric_dimension" => 76,
            "exact_periodic_n6_energy_per_site" => exact_energy,
            "exact_periodic_n6_commutator_residuals" => Dict("total_Sx" => exact_residuals[1], "total_Sz" => exact_residuals[2]),
        ),
        "point" => Dict(
            "baseline" => baseline,
            "hybrid" => hybrid,
            "improvement" => improvement,
            "acceptance_passed" => passed,
        ),
        "gates" => Dict(
            "both_optimal" => both_optimal,
            "minimum_improvement" => SU2_MIN_IMPROVEMENT,
            "stage_passed" => passed,
        ),
    )
    write_h1_json(result_path, results)
    plot_su2_invariance(results, figure_path)
    @printf("published baseline %.6f, SU(2) invariant %.6f\n", baseline["energy"], hybrid["energy"])
    @printf("raw objective improvement %.3e, passed=%s\n", improvement, passed)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    su2_main(ARGS)
end
