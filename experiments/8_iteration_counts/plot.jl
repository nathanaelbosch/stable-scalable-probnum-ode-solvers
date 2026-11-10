# Appendix: average Krylov iterations per linear solve of the MatfreeEK1, and average IEKF
# iterations per accepted step, on the sweeps of Figure 1 and over the tolerance on the
# problems of Figure 5. Rejected steps stop before the update, so we divide by accepted steps.
# Usage: julia --project experiments/8_iteration_counts/plot.jl
include("../../src/plotting.jl")

# The three kinds of linear solves, in colors that no method uses
const SOLVES = [("mean", "Mean update", :circle, "#1B9E77"),
                ("covariance", "Covariance samples", :rect, "#7570B3"),
                ("calibration", "Calibration", :utriangle, "#E7298A")]
const IEKF = ["DiagonalEK1+IEKF", "MatfreeEK1+IEKF"]
fig1(sweep) = read_rows(joinpath(@__DIR__, "data", sweep))
fig5(problem) = read_rows(joinpath(@__DIR__, "..", "5_workprecision", "data", problem))
const COLUMNS = [("a", "Domain growth", "Dimension d", fig1("domain"), r -> r.d),
                 ("b", "Stiffening", "Stiffness λ", fig1("stiffness"), r -> 0.05 / r.dx^2),
                 ("c", "Refinement", "Resolution 1/Δx", fig1("refinement"), r -> 1 / r.dx),
                 ("d", "Burgers", TAU_REL, fig5("burgers"), r -> r.rtol),
                 ("e", "Fisher-KPP", TAU_REL, fig5("fisher_kpp_2d"), r -> r.rtol)]

# One row of five panels with a shared y-axis
function counts_axes(ylabel, yscale)
    fig = Figure()
    axes = map(enumerate(COLUMNS)) do (k, (letter, title, xlabel, rows, x))
        xticks = xlabel === TAU_REL ? ([1e-8, 1e-4, 1e0], ["10⁻⁸", "10⁻⁴", "10⁰"]) : DECADES
        ax = Axis(fig[1, k]; title=panel_title(letter, title), xlabel, xscale=log10, xticks,
                  xreversed=xlabel === TAU_REL, yscale, ylabel=k == 1 ? ylabel : "",
                  (yscale == log10 ? (; yticks=DECADES) : (;))...)
        k > 1 && hideydecorations!(ax; grid=false)
        ax
    end
    return fig, axes
end

set_theme!(paper_theme(size=(397, 115)))

fig, axes = counts_axes("Avg. Krylov iterations\nper linear solve", log10)
for (k, (ax, (_, _, _, rows, x))) in enumerate(zip(axes, COLUMNS))
    mf = sort(filter(r -> r.method == "MatfreeEK1", rows); by=x)
    for (part, label, marker, color) in SOLVES
        niter = getproperty.(mf, Symbol(:niter_, part))
        nsolve = getproperty.(mf, Symbol(:nsolve_, part))
        scatterlines!(ax, x.(mf), niter ./ nsolve; label, marker, markersize=5,
                      color=(color, 0.6), markercolor=color)
    end
    # O(d) guide from the first mean-update count. The x-axis is d in (a) and 1/Δx = N in (c),
    # where d = 2N².
    if k in (1, 3)
        xs = x.(mf)
        y0 = mf[1].niter_mean / mf[1].nsolve_mean
        guide = lines!(ax, xs, y0 .* (xs ./ xs[1]) .^ (k == 1 ? 1 : 2); color=:black,
                       linewidth=0.6, linestyle=:dash, label="O(d)", yautolimits=false)
        translate!(guide, 0, 0, -1)  # behind the data, but last in the legend
    end
end
linkyaxes!(axes...)
Legend(fig[2, 1:5], axes[1]; orientation=:horizontal, tellheight=true)
save(joinpath(@__DIR__, "data", "krylov_counts.pdf"), fig)

fig, axes = counts_axes("Avg. IEKF iterations\nper step", identity)
for (ax, (_, _, _, rows, x)) in zip(axes, COLUMNS)
    plot_methods!(ax, sort(rows; by=x), IEKF, x, r -> r.nnonliniter / r.naccept; markerscale=0.75)
    ylims!(ax, 0, nothing)
end
linkyaxes!(axes...)
Legend(fig[2, 1:5], axes[1]; orientation=:horizontal, tellheight=true)
save(joinpath(@__DIR__, "data", "iekf_counts.pdf"), fig)
