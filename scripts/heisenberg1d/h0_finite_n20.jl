# Stage H0: finite periodic N=20 contract check for the spin-1/2 Heisenberg chain.
# Compares exact diagonalization, structured NPA, and the compressed-LTI ladder.

include(joinpath(@__DIR__, "t2t3_compressed_lti.jl"))

using QMBCertify
using XDiag

const H0_N = 20
const H0_D = 3
const H0_SUPER_DEPTHS = [4, 6, 8, 10]
const H0_STRUCTURED_ORDER = 2

function parse_h0_args(args)
    run_name = Dates.format(Dates.now(), "yyyymmdd-HHMMSS") * "-h0-finite-n20"
    for arg in args
        if startswith(arg, "--run-name=")
            run_name = split(arg, "=", limit=2)[2]
        else
            error("unknown argument: $arg")
        end
    end
    !isempty(run_name) && basename(run_name) == run_name || error("run_name must be one directory name")
    run_name
end

function periodic_ed_energy_per_site(N)
    hamiltonian = XDiag.OpSum()
    for site in 1:N
        hamiltonian += XDiag.Op("SdotS", [site, mod1(site + 1, N)])
    end
    XDiag.eigval0(hamiltonian, XDiag.Spinhalf(N, N ÷ 2); precision=1e-12) / N
end

function run_structured_npa(N, rundir)
    supp = [[1, 4]]
    coe = [3 / 4]
    log_path = joinpath(rundir, "qmbcertify.log")
    output = Ref{Any}()

    open(log_path, "w") do io
        redirect_stdout(io) do
            output[] = QMBCertify.GSB(
                supp,
                coe,
                N,
                H0_STRUCTURED_ORDER;
                lattice="chain",
                rdm=0,
                pso=0,
                lso=0,
                extra=0,
                Gram=true,
                QUIET=false,
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
    certificate = QMBCertify.certify_qmb(
        data,
        N,
        coe[1],
        optimum;
        tol_gram=1e-15,
        tol_dft=1e-12,
        snn=false,
        check=true,
    )
    inventory = match(r"SDP size: n = (\d+), m = (\d+)", log_text)

    Dict{String, Any}(
        "method" => "QMBCertify structured NPA",
        "system_scope" => "finite_periodic_chain",
        "system_size" => N,
        "boundary" => "periodic",
        "relaxation_order" => H0_STRUCTURED_ORDER,
        "extra" => 0,
        "rdm" => 0,
        "pso" => 0,
        "lso" => 0,
        "status" => "OPTIMAL",
        "status_source" => "QMBCertify emits a termination-status line only for non-OPTIMAL solves; the captured log contains a completed optimum and no such line",
        "numeric_energy_per_site" => Float64(certificate.oldbound),
        "certified_energy_per_site" => Float64(certificate.newbound),
        "certificate" => "rational_projection_with_rigorous_eigenvalue_shift",
        "certificate_shift" => Float64(certificate.shift),
        "certificate_min_eigenvalues" => Float64.(certificate.mineigs),
        "max_psd_block_dimension" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[1]),
        "constraint_count" => isnothing(inventory) ? nothing : parse(Int, inventory.captures[2]),
        "log" => basename(log_path),
    )
end

function compressed_record(solve, wall, n_super, fingerprint, exact_energy)
    energy = solve["energy"]
    adjusted = residual_adjusted_lower_bound(solve)
    published = published_lower_bound(solve)
    Dict{String, Any}(
        "method" => "Kull-Schuch compressed LTI",
        "system_scope" => "finite_periodic_chain",
        "system_size" => H0_N,
        "boundary" => "periodic",
        "support_sites" => 2 * n_super,
        "n_super" => n_super,
        "D" => H0_D,
        "coarse_grainer_fingerprint" => fingerprint,
        "status" => solve["status"],
        "energy_per_site" => energy,
        "residual_adjusted_lower_bound" => adjusted,
        "published_lower_bound" => published,
        "raw_finite_gap" => exact_energy - energy,
        "published_finite_gap" => exact_energy - published,
        "residuals" => solve["residuals"],
        "solver_wall_seconds" => solve["solver_wall_seconds"],
        "wall_s" => wall,
        "solver_tolerance" => 1e-7,
    )
end

function plot_h0(rundir, exact_energy, compressed, structured)
    supports = [point["support_sites"] for point in compressed]
    compressed_gaps = exact_energy .- [point["residual_adjusted_lower_bound"] for point in compressed]
    structured_gap = exact_energy - structured["certified_energy_per_site"]

    p = plot(
        supports,
        compressed_gaps;
        xscale=:log10,
        yscale=:log10,
        marker=:circle,
        xlabel="constraint support (physical sites)",
        ylabel="finite-N gap: E_exact(20)/20 - E_lower",
        title="N=20 periodic Heisenberg chain: conservative lower bounds",
        label="compressed LTI, D=$H0_D",
        legend=:bottomleft,
        size=(800, 500),
    )
    scatter!(
        p,
        [H0_N],
        [structured_gap];
        marker=:diamond,
        markersize=7,
        label="structured NPA, d=$H0_STRUCTURED_ORDER",
    )
    savefig(p, joinpath(rundir, "figs", "finite_n20_gap.png"))
end

function h0_main(args)
    run_name = parse_h0_args(args)
    rundir = normpath(joinpath(@__DIR__, "..", "..", "results", run_name))
    ispath(rundir) && error("result path already exists; choose a fresh --run-name: $rundir")
    mkpath(joinpath(rundir, "figs"))
    json_path = joinpath(rundir, "finite_n20.json")

    results = Dict{String, Any}(
        "meta" => Dict(
            "stage" => "H0",
            "numerical_bound_schema_version" => NUMERICAL_BOUND_SCHEMA_VERSION,
            "run_name" => run_name,
            "model" => "spin_half_heisenberg_antiferromagnet",
            "J" => 1,
            "system_scope" => "finite_periodic_chain",
            "system_size" => H0_N,
            "boundary" => "periodic",
            "rerun_command" => "julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/h0_finite_n20.jl --run-name=$run_name",
            "contract" => "docs/FINITE_SIZE_CONTRACT.md",
        ),
    )

    println("Computing exact periodic N=$H0_N reference with XDiag...")
    exact_t0 = time()
    exact_energy = periodic_ed_energy_per_site(H0_N)
    results["exact_reference"] = Dict(
        "method" => "XDiag Lanczos exact diagonalization",
        "sector" => "Sz_total=0",
        "energy_per_site" => exact_energy,
        "wall_s" => time() - exact_t0,
    )
    write_results(json_path, results)
    @printf("Exact E0/N = %.12f\n", exact_energy)

    summary_path = joinpath(COARSE_GRAINER_DIR, "vumps_summary.json")
    isfile(summary_path) || error("missing VUMPS summary; run t1_vumps.jl first")
    vumps_summary = JSON.parsefile(summary_path)
    vumps_summary["seed"] == VUMPS_DEFAULT_SEED || error("D=3 VUMPS base seed is stale")
    vumps_summary["julia_version"] == string(PINNED_JULIA_VERSION) ||
        error("D=3 VUMPS Julia version is stale")
    d3_vumps = vumps_summary["runs"]["D$H0_D"]
    d3_vumps["random_seed"] == VUMPS_DEFAULT_SEED + H0_D ||
        error("D=3 VUMPS run seed is stale")
    d3_vumps["saved"] && d3_vumps["converged"] ||
        error("D=3 VUMPS coarse-grainer did not converge")

    psis = deserialize(joinpath(COARSE_GRAINER_DIR, "vumps_mps.jls"))
    haskey(psis, H0_D) || error("missing D=$H0_D coarse-grainer")
    psi = psis[H0_D]
    fingerprint = coarse_grainer_fingerprint(psi)
    results["coarse_grainer"] = Dict(
        "fingerprint" => fingerprint,
        "vumps" => d3_vumps,
    )
    write_results(json_path, results)
    B = super_tensor(psi)
    compressed = Dict{String, Any}[]

    for n_super in H0_SUPER_DEPTHS
        2 * n_super <= H0_N || error("support exceeds finite system size")
        t0 = time()
        solve = solve_compressed(B, n_super; silent=true, diagnostics=true)
        wall = time() - t0
        solve["status"] == "OPTIMAL" ||
            error("compressed point n_super=$n_super did not terminate OPTIMAL: $(solve["status"])")
        energy = solve["energy"]
        energy <= exact_energy ||
            error("finite-size sandwich violated at n_super=$n_super")
        residual_adjusted_lower_bound(solve) < exact_energy ||
            error("residual-adjusted finite-size sandwich violated at n_super=$n_super")
        push!(compressed, compressed_record(solve, wall, n_super, fingerprint, exact_energy))
        results["compressed"] = compressed
        write_results(json_path, results)
        @printf(
            "Compressed D=%d n_super=%2d: raw E/N = %.9f, published lower bound = %.6f, status %s\n",
            H0_D,
            n_super,
            energy,
            last(compressed)["published_lower_bound"],
            solve["status"],
        )
    end

    compressed_energies = [point["energy_per_site"] for point in compressed]
    all(diff(compressed_energies) .>= -2e-7) ||
        error("compressed fixed-D ladder is not monotone")

    println("Running the uncompressed structured-NPA comparison...")
    structured = run_structured_npa(H0_N, rundir)
    structured["certified_energy_per_site"] <= exact_energy ||
        error("structured-NPA finite-size sandwich violated")
    results["structured_npa"] = structured
    results["verification"] = Dict(
        "all_compressed_optimal" => all(point["status"] == "OPTIMAL" for point in compressed),
        "compressed_monotone" => true,
        "all_raw_solver_objectives_below_exact" => all(point["energy_per_site"] <= exact_energy for point in compressed),
        "all_raw_objectives_within_solver_tolerance" => all(point["energy_per_site"] <= exact_energy + point["solver_tolerance"] for point in compressed),
        "all_residual_adjusted_bounds_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in compressed),
        "all_bounds_below_exact" => all(point["residual_adjusted_lower_bound"] < exact_energy for point in compressed) && structured["certified_energy_per_site"] <= exact_energy,
        "structured_rational_postcertificate" => true,
    )
    write_results(json_path, results)
    plot_h0(rundir, exact_energy, compressed, structured)

    println("DONE. finite-size H0 artifacts saved to $rundir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    h0_main(ARGS)
end
