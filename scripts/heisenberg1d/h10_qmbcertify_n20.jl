# H10 — one-point QMBCertify structured-NPA ladder for the periodic N=20 chain.

using Random

include(joinpath(@__DIR__, "t2t3_compressed_lti.jl"))

const H10_N = 20
const H10_RANDOM_SEED = 0x48413130

function parse_h10_args(args)
    run_name = nothing
    order = 3
    extra = 0
    rdm = 0
    pso = 0
    lso = false
    screen_only = false
    mosek_threads = 16
    upper_result = nothing
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--order=")
            order = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--extra=")
            extra = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--rdm=")
            rdm = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--pso=")
            pso = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--lso=")
            lso = parse(Int, split(arg, "=", limit=2)[2]) == 1
        elseif startswith(arg, "--mosek-threads=")
            mosek_threads = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--upper-result=")
            spec = split(arg, "=", limit=2)[2]
            upper_result = normpath(isabspath(spec) ? spec : joinpath(BOUNDGSENERGY_ROOT, spec))
        elseif arg == "--screen-only"
            screen_only = true
        else
            error("unknown argument: $arg")
        end
    end
    order >= 2 || error("order must be at least 2")
    extra >= 0 || error("extra must be nonnegative")
    rdm in (0, 8, 9, 10) || error("rdm must be 0, 8, 9, or 10")
    pso >= 0 || error("pso must be nonnegative")
    mosek_threads >= 1 || error("mosek_threads must be positive")
    isnothing(run_name) && (run_name = Dates.format(now(), "yyyymmdd-HHMMSS") * "-h10-qmbcertify-n20")
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    isnothing(upper_result) && error("--upper-result is required")
    isfile(upper_result) || error("finite-energy upper result does not exist: $upper_result")
    run_name, order, extra, rdm, pso, lso, mosek_threads, screen_only, upper_result
end

function plot_h10(results, figure_path)
    gap = results["comparison"]["certified_gap"]
    target = results["gates"]["target_gap"]
    plot = bar(
        ["certified gap"],
        [gap];
        ylabel="energy density",
        title="QMBCertify N=20, d=$(results["meta"]["relaxation_order"])",
        legend=false,
        color=:steelblue,
        ylims=(0, 1.12 * max(gap, target)),
    )
    hline!(plot, [target]; color=:black, linestyle=:dash, label=false)
    mkpath(dirname(figure_path))
    savefig(plot, figure_path)
end

function h10_main(args)
    run_name, order, extra, rdm, pso, lso, mosek_threads, screen_only, upper_result_path = parse_h10_args(args)
    result_dir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    result_path = joinpath(result_dir, "qmbcertify_n20.json")
    log_path = joinpath(result_dir, "qmbcertify.log")
    figure_path = joinpath(result_dir, "figs", "qmbcertify_n20.png")
    ispath(result_dir) && error("result directory already exists; choose a fresh --run-name")
    mkpath(result_dir)

    supp = [[1, 4]]
    coe = [3 / 4]
    output = Ref{Any}()
    Random.seed!(H10_RANDOM_SEED)
    solve_wall = @elapsed open(log_path, "w") do io
        redirect_stdout(io) do
            output[] = QMBCertify.GSB(
                supp,
                coe,
                H10_N,
                order;
                lattice="chain",
                rdm,
                extra,
                pso,
                lso,
                lol=H10_N,
                three_type=[1, 1],
                Gram=true,
                QUIET=false,
                mosek_setting=QMBCertify.mosek_para(1e-8, 1e-8, 1e-8, mosek_threads),
            )
        end
    end
    log_text = read(log_path, String)
    print(log_text)
    occursin("termination status:", log_text) &&
        error("QMBCertify did not terminate OPTIMAL; see $log_path")
    occursin("optimum =", log_text) ||
        error("QMBCertify log did not contain a completed optimum")

    optimum, data = output[]
    upper_result = JSON.parsefile(upper_result_path)
    inventory = match(r"SDP size: n = (\d+), m = (\d+)", log_text)
    if screen_only
        reference = upper_result["comparison"]["hybrid"]
        raw_gap = upper_result["finite_energy_upper_certificate"]["upper_bound"] - Float64(optimum)
        results = Dict(
            "meta" => Dict(
                "stage" => "H10",
                "run_name" => run_name,
                "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
                "system_scope" => "finite_periodic_chain",
                "system_size" => H10_N,
                "boundary" => "periodic",
                "method" => "QMBCertify structured NPA numerical screen",
                "relaxation_order" => order,
                "extra" => extra,
                "rdm" => rdm,
                "pso" => pso,
                "lso" => lso,
                "mosek_threads" => mosek_threads,
                "random_seed" => string(H10_RANDOM_SEED),
            ),
            "solve" => Dict(
                "status" => "OPTIMAL",
                "numeric_energy_per_site" => Float64(optimum),
                "solve_wall_seconds" => solve_wall,
                "reported_sdp_basis_size" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[1]),
                "constraint_count" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[2]),
                "log" => basename(log_path),
            ),
            "comparison" => Dict(
                "numeric_gap_to_exact_rayleigh_upper" => raw_gap,
                "numeric_improvement_over_h7_rdm8_objective" => Float64(optimum) - reference["energy"],
            ),
            "gates" => Dict(
                "exact_postcertificate_run" => false,
                "beats_h7_rdm8_raw_objective" => Float64(optimum) > reference["energy"],
                "stage_passed" => false,
            ),
        )
        write_results(result_path, results)
        @printf("numeric %.12f, raw gap %.3e\n", optimum, raw_gap)
        println("beats H7 RDM8 raw objective=$(results["gates"]["beats_h7_rdm8_raw_objective"])")
        println("saved screen $result_path")
        return
    end
    certificate = nothing
    certificate_wall = @elapsed certificate = QMBCertify.certify_qmb(
        data,
        H10_N,
        coe[1],
        optimum;
        tol_gram=1e-15,
        tol_dft=1e-12,
        snn=false,
        check=true,
    )
    certified_lower = certificate.newbound_rat
    certified_upper = parse_exact_rational(
        upper_result["finite_energy_upper_certificate"]["upper_bound_rational"],
    )
    gap = certified_upper - certified_lower
    gap >= 0 || error("QMBCertify lower bound exceeds the exact Rayleigh upper bound")
    target_gap = 1 // 100000

    results = Dict(
        "meta" => Dict(
            "stage" => "H10",
            "run_name" => run_name,
            "model" => "spin-1/2 Heisenberg antiferromagnetic chain",
            "system_scope" => "finite_periodic_chain",
            "system_size" => H10_N,
            "boundary" => "periodic",
            "method" => "QMBCertify structured NPA",
            "relaxation_order" => order,
            "extra" => extra,
            "rdm" => rdm,
            "pso" => pso,
            "lso" => lso,
            "mosek_threads" => mosek_threads,
            "random_seed" => string(H10_RANDOM_SEED),
        ),
        "solve" => Dict(
            "status" => "OPTIMAL",
            "numeric_energy_per_site" => Float64(optimum),
            "solve_wall_seconds" => solve_wall,
            "postcertificate_wall_seconds" => certificate_wall,
            "reported_sdp_basis_size" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[1]),
            "constraint_count" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[2]),
            "log" => basename(log_path),
        ),
        "certificate" => Dict(
            "status" => "CERTIFIED",
            "method" => "QMBCertify rational SOHS projection with rigorous minimum-eigenvalue shift",
            "raw_bound_rational" => string(certificate.oldbound_rat),
            "certified_lower_bound" => float_down(certified_lower),
            "certified_lower_bound_rational" => string(certified_lower),
            "shift" => Float64(certificate.shift),
            "shift_rational" => string(certificate.shift_rat),
            "minimum_eigenvalues" => Float64.(certificate.mineigs),
            "minimum_eigenvalues_rational" => string.(certificate.mineigs_rat),
            "basis_dimensions" => certificate.dims,
        ),
        "finite_energy_upper_certificate" => upper_result["finite_energy_upper_certificate"],
        "comparison" => Dict(
            "certified_gap" => float_up(gap),
            "certified_gap_rational" => string(gap),
        ),
        "gates" => Dict(
            "target_gap" => Float64(target_gap),
            "sandwich_passed" => gap >= 0,
            "target_passed" => gap <= target_gap,
            "stage_passed" => gap <= target_gap,
        ),
    )
    write_results(result_path, results)
    plot_h10(results, figure_path)

    @printf("numeric %.12f, certified %.12f\n", optimum, Float64(certified_lower))
    @printf("certified gap %.3e, target passed=%s\n", Float64(gap), results["gates"]["target_passed"])
    @printf("solve %.2fs, postcertificate %.2fs\n", solve_wall, certificate_wall)
    println("saved $result_path")
    println("saved $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h10_main(ARGS)
end
