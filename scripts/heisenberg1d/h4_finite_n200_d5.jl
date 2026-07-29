# H4 — cluster-scale finite N=200 D=5 ladder.

include(joinpath(@__DIR__, "h4_finite_n100_d5.jl"))

const H4_N200 = 200
const H4_N200_DEPTHS = (50, 75, 100)

function parse_h4_n200_args(args)
    run_name = nothing
    n_super = 50
    mosek_threads = H4_MOSEK_THREADS
    n100_result = nothing
    previous_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--n-super=")
            n_super = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--n100-result=")
            n100_result = h4_input_path(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--previous-result=")
            previous_result = h4_input_path(split(arg, "=", limit=2)[2])
        else
            error("unknown argument: $arg")
        end
    end
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h4-finite-n200-d5-n$(n_super)")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    n_super in H4_N200_DEPTHS || error("n_super must be 50, 75, or 100")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(n100_result) && error("--n100-result is required")
    isfile(n100_result) || error("D=5 N=100 gate does not exist: $n100_result")
    if n_super == first(H4_N200_DEPTHS)
        isnothing(previous_result) || error("n=50 starts from the N=100 gate and must not take --previous-result")
    else
        isnothing(previous_result) && error("n=$n_super requires --previous-result")
        isfile(previous_result) || error("previous N=200 result does not exist: $previous_result")
    end
    run_name, n_super, mosek_threads, n100_result, previous_result
end

function expected_previous_depth(n_super)
    index = findfirst(==(n_super), H4_N200_DEPTHS)
    isnothing(index) && error("unsupported N=200 depth")
    index == 1 ? nothing : H4_N200_DEPTHS[index - 1]
end

function validate_h4_n200_result(path)
    result = JSON.parsefile(path)
    meta = result["meta"]
    meta["stage"] == "H4" || error("D=5 N=200 result is not an H4 result")
    meta["system_scope"] == "finite_periodic_chain" && meta["system_size"] == H4_N200 && meta["boundary"] == "periodic" || error("D=5 N=200 result metadata is inconsistent")
    meta["bond_dimension"] == H4_D || error("D=5 N=200 result has the wrong bond dimension")
    meta["n_super"] in H4_N200_DEPTHS || error("D=5 N=200 result has an unsupported depth")

    n100_path = h4_recorded_input(result, "n100_result")
    n100 = validate_h4_n100_result(n100_path)
    meta["coarse_grainer_fingerprint"] == n100["meta"]["coarse_grainer_fingerprint"] || error("D=5 N=200 fingerprint does not match the N=100 gate")
    previous_depth = expected_previous_depth(meta["n_super"])
    previous_energy = if isnothing(previous_depth)
        n100["D5_point"]["energy_per_site"]
    else
        previous = validate_h4_n200_result(h4_recorded_input(result, "previous_result"))
        previous["meta"]["n_super"] == previous_depth || error("D=5 N=200 predecessor has the wrong depth")
        previous["meta"]["coarse_grainer_fingerprint"] == meta["coarse_grainer_fingerprint"] || error("D=5 N=200 predecessor fingerprint is inconsistent")
        realpath(h4_recorded_input(previous, "n100_result")) == realpath(n100_path) || error("D=5 N=200 predecessor uses a different N=100 gate")
        previous["D5_point"]["energy_per_site"]
    end

    expected_reference = bethe_ground_reference(H4_N200)
    exact_energy = result["exact_reference"]["energy_per_site"]
    isapprox(exact_energy, expected_reference["energy_per_site"]; atol=2e-12, rtol=0) || error("D=5 N=200 exact reference is inconsistent")
    point = result["D5_point"]
    adjusted, published = validated_solver_bounds(point["solve"], "D=5 N=200 stored solve")
    point["energy_per_site"] == point["solve"]["energy"] || error("D=5 N=200 energy snapshot is stale")
    point["residual_adjusted_lower_bound"] == adjusted || error("D=5 N=200 residual-adjusted snapshot is stale")
    point["published_lower_bound"] == published || error("D=5 N=200 published snapshot is stale")
    point["system_scope"] == "finite_periodic_chain" && point["system_size"] == H4_N200 && point["boundary"] == "periodic" && point["bond_dimension"] == H4_D || error("D=5 N=200 point metadata is inconsistent")
    point["n_super"] == meta["n_super"] && point["support_sites"] == meta["support_sites"] == 2 * meta["n_super"] || error("D=5 N=200 support metadata is inconsistent")
    point["raw_finite_gap"] == exact_energy - point["energy_per_site"] || error("D=5 N=200 raw finite gap is stale")
    point["finite_gap"] == exact_energy - published || error("D=5 N=200 published finite gap is stale")
    published <= adjusted <= point["energy_per_site"] < exact_energy || error("D=5 N=200 point violates the finite ceiling")
    point["energy_per_site"] >= previous_energy - 2e-7 || error("D=5 N=200 hierarchy is not monotone")
    result["previous_energy_per_site"] == previous_energy || error("D=5 N=200 previous-energy snapshot is stale")
    result["gates"]["fixed_D_monotone"] || error("D=5 N=200 monotonicity gate did not pass")
    result["gates"]["target_gap"] == 1e-5 || error("D=5 N=200 target gate has the wrong threshold")
    !result["gates"]["strict_certificate_implemented"] || error("D=5 N=200 unexpectedly claims strict certification")
    result["gates"]["numerical_target_passed"] == (point["finite_gap"] <= result["gates"]["target_gap"]) || error("D=5 N=200 numerical target gate is stale")
    !result["gates"]["target_passed"] || error("D=5 N=200 cannot pass the challenge target without strict certification")
    result
end

function plot_h4_n200(results, figure_path)
    point = results["D5_point"]
    target = results["gates"]["target_gap"]
    gap_plot = scatter(1:2, [point["finite_gap"], target]; yscale=:log10, xticks=(1:2, ["numerical gap", "1e-5 target"]), ylabel="energy-density gap", title="N=200 D=5 numerical screen", label=false, markersize=9, markercolor=[:steelblue, :black], ylims=(min(point["finite_gap"], target) / 2, max(point["finite_gap"], target) * 2), left_margin=8mm, bottom_margin=8mm)
    values = [point["published_lower_bound"], results["exact_reference"]["energy_per_site"]]
    energy_plot = scatter(1:2, values; xticks=(1:2, ["residual-adjusted", "Bethe reference"]), ylabel="energy density", title="N=200 numerical comparison", label=false, markersize=9, markercolor=[:steelblue, :darkorange], ylims=(minimum(values) - 0.1 * abs(diff(values)[1]), maximum(values) + 0.1 * abs(diff(values)[1])), left_margin=8mm, bottom_margin=8mm)
    mkpath(dirname(figure_path))
    savefig(plot(gap_plot, energy_plot; layout=(1, 2), size=(1100, 470)), figure_path)
end

function h4_n200_main(args)
    run_name, n_super, mosek_threads, n100_path, previous_path = parse_h4_n200_args(args)
    result_dir = joinpath(H4_REPO_ROOT, "results", run_name)
    result_path = joinpath(result_dir, "finite_n200_d5_n$(n_super).json")
    figure_path = joinpath(result_dir, "figs", "finite_n200_d5_n$(n_super).png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")

    n100 = validate_h4_n100_result(n100_path)
    previous_depth = expected_previous_depth(n_super)
    previous_energy = if isnothing(previous_depth)
        n100["D5_point"]["energy_per_site"]
    else
        previous = validate_h4_n200_result(previous_path)
        previous["meta"]["n_super"] == previous_depth || error("previous N=200 result has the wrong depth")
        realpath(h4_recorded_input(previous, "n100_result")) == realpath(n100_path) || error("previous N=200 result uses a different N=100 gate")
        previous["D5_point"]["energy_per_site"]
    end

    exact_reference = bethe_ground_reference(H4_N200)
    exact_energy = exact_reference["energy_per_site"]
    virtual_charges = [0, 0, -2, 0, 2]
    super_charges = [-2, 0, 0, 2]
    B_raw, vumps, loaded_from_cache = load_or_build_d5_coarse_grainer()
    B, largest_forbidden = exact_charge_tensor(B_raw, virtual_charges, super_charges)
    map_violations = assert_charge_preserving_maps(B, virtual_charges, super_charges)
    fingerprint = bytes2hex(sha1(reinterpret(UInt8, vec(B))))
    fingerprint == n100["meta"]["coarse_grainer_fingerprint"] || error("cached D=5 coarse-grainer does not match the N=100 gate")

    solve = solve_formulation(B, n_super, virtual_charges, super_charges; blocked=true, mosek_threads)
    adjusted, published = validated_solver_bounds(solve, "D=5 finite N=200 solve")
    energy = solve["energy"]
    energy < exact_energy || error("D=5 raw objective violates the finite N=200 ceiling")
    adjusted < exact_energy || error("D=5 residual-adjusted bound violates the finite N=200 ceiling")
    energy >= previous_energy - 2e-7 || error("D=5 N=200 hierarchy is not monotone")

    point = Dict(
        "system_scope" => "finite_periodic_chain",
        "system_size" => H4_N200,
        "boundary" => "periodic",
        "support_sites" => 2 * n_super,
        "n_super" => n_super,
        "bond_dimension" => H4_D,
        "energy_per_site" => energy,
        "residual_adjusted_lower_bound" => adjusted,
        "published_lower_bound" => published,
        "raw_finite_gap" => exact_energy - energy,
        "finite_gap" => exact_energy - published,
        "solve" => solve,
    )
    inputs = Dict("n100_result" => relpath(n100_path, H4_REPO_ROOT))
    isnothing(previous_path) || (inputs["previous_result"] = relpath(previous_path, H4_REPO_ROOT))
    results = Dict(
        "meta" => Dict(
            "stage" => "H4",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H4_N200,
            "boundary" => "periodic",
            "contract" => "docs/FINITE_SIZE_CONTRACT.md",
            "support_sites" => 2 * n_super,
            "n_super" => n_super,
            "bond_dimension" => H4_D,
            "solver_feasibility_tolerance" => H1_FEAS_TOL,
            "energy_precision_digits" => 6,
            "charge_units" => "2Sz",
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "mosek_threads" => mosek_threads,
            "coarse_grainer_fingerprint" => fingerprint,
            "coarse_grainer_loaded_from_cache" => loaded_from_cache,
            "inputs" => inputs,
        ),
        "exact_reference" => exact_reference,
        "vumps" => vumps,
        "charge_check" => Dict("largest_forbidden_entry" => largest_forbidden, "map_violations" => map_violations),
        "previous_energy_per_site" => previous_energy,
        "D5_point" => point,
        "gates" => Dict(
            "fixed_D_monotone" => true,
            "target_gap" => 1e-5,
            "strict_certificate_implemented" => false,
            "numerical_target_passed" => point["finite_gap"] <= 1e-5,
            "target_passed" => false,
        ),
    )
    write_h1_json(result_path, results)
    validate_h4_n200_result(result_path)
    plot_h4_n200(results, figure_path)
    @printf("N=200 D=5 n=%d raw E %.9f, residual-adjusted %.6f, diagnostic gap %.3e, numerical target %s\n", n_super, energy, published, point["finite_gap"], results["gates"]["numerical_target_passed"])
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h4_n200_main(ARGS)
end
