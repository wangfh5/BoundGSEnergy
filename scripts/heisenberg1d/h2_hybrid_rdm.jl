# H2 — add an eight-physical-site RDM lift to the compressed-LTI relaxation.

include(joinpath(@__DIR__, "symmetry_probe.jl"))

const H2_MOSEK_THREADS = 16

function parse_h2_args(args)
    run_name = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h2-hybrid-rdm")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name
end

function h2_result_paths(run_name)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_dir, joinpath(result_dir, "hybrid_rdm.json"), joinpath(result_dir, "figs", "hybrid_rdm.png")
end

function tensor_fingerprint(B)
    bytes2hex(sha1(reinterpret(UInt8, vec(B))))
end

function solve_pair(B, virtual_charges, super_charges; blocked)
    GC.gc()
    baseline = solve_formulation(B, 4, virtual_charges, super_charges; blocked, mosek_threads=H2_MOSEK_THREADS)
    GC.gc()
    hybrid = solve_formulation(B, 4, virtual_charges, super_charges; blocked, rdm4=true, mosek_threads=H2_MOSEK_THREADS)
    both_optimal = baseline["status"] == "OPTIMAL" && hybrid["status"] == "OPTIMAL"
    improvement = both_optimal ? hybrid["energy"] - baseline["energy"] : nothing
    passed = both_optimal && improvement > 2e-7 && hybrid["energy"] >= baseline["energy"] - 2e-7 && baseline["residual_adjusted_lower_bound"] < E0 && hybrid["residual_adjusted_lower_bound"] < E0
    Dict(
        "baseline" => baseline,
        "hybrid" => hybrid,
        "improvement" => improvement,
        "acceptance_passed" => passed,
    )
end

function plot_hybrid(results, figure_path)
    labels = ["D=3 dense", "D=4 U(1)"]
    x = collect(1:length(labels))
    points = [results["points"]["D3"], results["points"]["D4"]]
    deltas = [abs(point["improvement"]) for point in points]
    residuals = [max(solver_residual_scale(point["baseline"]), solver_residual_scale(point["hybrid"])) for point in points]
    threshold = results["gates"]["minimum_improvement"]
    scales = vcat(deltas, residuals, threshold)
    delta_plot = scatter(x, deltas; marker=:circle, markersize=8, yscale=:log10, label="|raw objective Δ|", ylabel="energy-density scale", title="H2 RDM effect vs numerical scale", xticks=(x, labels), ylims=(minimum(scales) / 2, maximum(scales) * 2), left_margin=6mm, bottom_margin=6mm)
    scatter!(delta_plot, x, residuals; marker=:diamond, markersize=8, label="solver residual scale")
    plot!(delta_plot, x, fill(threshold, length(x)); linestyle=:dash, color=:black, label="2e-7 acceptance threshold")
    baseline_walls = [results["points"]["D3"]["baseline"]["solver_wall_seconds"], results["points"]["D4"]["baseline"]["solver_wall_seconds"]]
    hybrid_walls = [results["points"]["D3"]["hybrid"]["solver_wall_seconds"], results["points"]["D4"]["hybrid"]["solver_wall_seconds"]]
    wall_limit = 1.08 * maximum(vcat(baseline_walls, hybrid_walls))
    wall_plot = bar(x .- 0.18, baseline_walls; bar_width=0.35, label="compressed", ylabel="Mosek wall time (s)", title="H2 solver cost", color=:gray60, ylims=(0, wall_limit), left_margin=5mm, bottom_margin=5mm)
    bar!(wall_plot, x .+ 0.18, hybrid_walls; bar_width=0.35, label="hybrid RDM", color=:darkorange)
    plot!(wall_plot; xticks=(x, labels))
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1100, 460)), figure_path)
end

function main(args)
    run_name = parse_h2_args(args)
    result_dir, result_path, figure_path = h2_result_paths(run_name)
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")
    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    B3 = super_tensor(psis[3])
    B4_raw, vumps4, loaded_from_cache = load_or_build_coarse_grainer()
    virtual_charges = [-1, 1, -1, 1]
    super_charges = [-2, 0, 0, 2]
    B4, largest_forbidden = exact_charge_tensor(B4_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B4, virtual_charges, super_charges)
    rdm8_map = validate_rdm8_map(B4, virtual_charges, super_charges)
    compression_violation = rdm8_map["compression_violation"]
    marginal_identity_residuals = rdm8_map["marginal_identity_residuals"]

    rank3 = rank(W2map(B3); atol=1e-10)
    rank4 = rank(W2map(B4); atol=1e-10)
    @printf("W2 ranks: D=3 %d/16, symmetric D=4 %d/16\n", rank3, rank4)
    flush(stdout)

    D3 = solve_pair(B3, Int[], Int[]; blocked=false)
    D4 = solve_pair(B4, virtual_charges, super_charges; blocked=true)
    rho4_inventory = block_inventory(rho4_charges(super_charges))

    results = Dict(
        "meta" => Dict(
            "stage" => "H2",
            "run_name" => run_name,
            "strengthening" => "four_super_site_rdm",
            "physical_rdm_sites" => 8,
            "contract" => "docs/HYBRID_RDM_CONTRACT.md",
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => SYSTEM_SCOPE,
            "system_size" => nothing,
            "support_sites" => 8,
            "n_super" => 4,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => H2_MOSEK_THREADS,
        ),
        "coarse_grainers" => Dict(
            "D3" => Dict("fingerprint" => tensor_fingerprint(B3), "W2_rank" => rank3, "W2_columns" => 16, "tensor_symmetry" => "none"),
            "D4" => Dict("fingerprint" => tensor_fingerprint(B4), "W2_rank" => rank4, "W2_columns" => 16, "tensor_symmetry" => "SU(2)", "largest_forbidden_u1_entry" => largest_forbidden, "map_violations" => map_violations, "rho4_compression_violation" => compression_violation, "marginal_identity_residuals" => marginal_identity_residuals, "loaded_from_cache" => loaded_from_cache, "vumps" => vumps4),
        ),
        "rdm4_block_inventory" => Dict(
            "dense" => [256],
            "u1" => rho4_inventory,
            "largest_u1_block" => maximum(values(rho4_inventory)),
        ),
        "points" => Dict("D3" => D3, "D4" => D4),
        "gates" => Dict(
            "D3_passed" => D3["acceptance_passed"],
            "D4_passed" => D4["acceptance_passed"],
            "stage_passed" => D3["acceptance_passed"] && D4["acceptance_passed"],
            "minimum_improvement" => 2e-7,
        ),
    )
    write_h1_json(result_path, results)
    plot_hybrid(results, figure_path)
    for D in ("D3", "D4")
        point = results["points"][D]
        @printf("%s baseline %.9f, hybrid %.9f, improvement %.3e, passed=%s\n", D, point["baseline"]["energy"], point["hybrid"]["energy"], point["improvement"], point["acceptance_passed"])
    end
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
