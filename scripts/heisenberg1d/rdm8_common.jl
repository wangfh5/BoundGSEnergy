# Shared constants and plotting for the exact-compatible RDM8 runs.

include(joinpath(@__DIR__, "h2_hybrid_rdm.jl"))

const RDM8_N_SUPER = 6
const RDM8_DYADIC_BITS = 13
const RDM8_MIN_IMPROVEMENT = 2e-7
const RDM8_DEFAULT_MOSEK_THREADS = 16

function plot_rdm8_comparison(results, figure_path)
    comparison = results["comparison"]
    baseline = comparison["baseline"]
    hybrid = comparison["hybrid"]
    scales = [
        abs(comparison["raw_improvement"]),
        comparison["combined_solver_residual_scale"],
        results["gates"]["minimum_improvement"],
    ]
    delta_plot = scatter(
        1:3,
        scales;
        yscale=:log10,
        xticks=(1:3, ["|raw objective Δ|", "sum of residual scales", "acceptance threshold"]),
        ylabel="energy-density scale",
        title="RDM8 effect at D=4, n=6",
        label=false,
        markersize=8,
        markercolor=[:steelblue, :darkorange, :black],
        ylims=(minimum(scales) / 2, maximum(scales) * 2),
        left_margin=6mm,
        bottom_margin=8mm,
    )
    walls = [baseline["solver_wall_seconds"], hybrid["solver_wall_seconds"]]
    wall_plot = bar(
        ["compressed", "+ RDM8"],
        walls;
        ylabel="Mosek wall time (s)",
        title="Local screening cost",
        legend=false,
        color=[:gray60, :darkorange],
        ylims=(0, 1.08 * maximum(walls)),
        left_margin=5mm,
        bottom_margin=5mm,
    )
    mkpath(dirname(figure_path))
    savefig(plot(delta_plot, wall_plot; layout=(1, 2), size=(1100, 460)), figure_path)
end
