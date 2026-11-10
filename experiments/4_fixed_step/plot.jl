# Figure 4: the solutions of Burgers' and the Fisher-KPP equation (a, b), and the final error
# of fixed-step solves against the step size (c, d).
# Usage: julia --project experiments/4_fixed_step/plot.jl
include("../../src/plotting.jl")
using DelimitedFiles

const METHODS = ["EK1", "ExpEK", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF"]
const PROBLEMS = [("burgers", "Burgers", "1", :balance, (-0.5, 0.5)),
                  ("fisher_kpp", "Fisher-KPP", "2", :viridis, (0, 1))]

set_theme!(paper_theme(size=(397, 97)))
fig = Figure()
axes = []
for (i, (key, title, T, colormap, colorrange)) in enumerate(PROBLEMS)
    dir = joinpath(@__DIR__, "data", key)
    U = readdlm(joinpath(dir, "solution.txt"), ',')
    ax = Axis(fig[i, 1][1, 1]; xlabel=i == 2 ? "Space" : "", ylabel="Time", xlabelpadding=-1,
              xticks=([1, size(U, 1)], ["0", "1"]), yticks=([1, size(U, 2)], ["0", T]),
              xticklabelsvisible=i == 2, title=panel_title("ab"[i], title))
    hm = heatmap!(ax, U; colormap, colorrange)
    Colorbar(fig[i, 1][1, 2], hm; size=3)

    ax = Axis(fig[1:2, 1+i]; xlabel="Step size", ylabel=i == 1 ? "Final error" : "",
              xscale=log10, yscale=log10, title=panel_title("cd"[i], title))
    plot_methods!(ax, read_rows(dir), METHODS, r -> r.dt, r -> r.error)
    push!(axes, ax)
end
linkyaxes!(axes...)
hideydecorations!(axes[2]; grid=false)
Legend(fig[1:2, 4], axes[1])
colsize!(fig.layout, 1, Relative(0.15))
rowgap!(fig.layout, 1, 2)
save(joinpath(@__DIR__, "data", "fixed_step.pdf"), fig)
