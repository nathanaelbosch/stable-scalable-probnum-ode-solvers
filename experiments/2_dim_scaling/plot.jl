# Figure 2: a Lorenz96 trajectory (a), and the runtime of a single step against the
# dimension d (b), with the slopes of O(d) and O(d³) costs.
# Usage: julia --project experiments/2_dim_scaling/plot.jl
include("../../src/plotting.jl")
using DelimitedFiles

const METHODS = ["EK0", "DiagonalEK1", "EK1", "ExpEK", "DiagonalEK1+IEKF", "MatfreeEK1",
                 "MatfreeEK1+IEKF"]
dir = joinpath(@__DIR__, "data")

set_theme!(paper_theme(size=(397, 104)))
fig = Figure()

U = readdlm(joinpath(dir, "trajectory.txt"), ',')
ax = Axis(fig[1, 1]; xlabel="Time", ylabel="State index", title=panel_title("a", "Lorenz96 (d=40)"))
heatmap!(ax, range(0, 30, length=size(U, 2)), 1:size(U, 1), U'; colormap=:viridis)

ax = Axis(fig[1, 2]; xlabel="ODE dimension d", ylabel="Runtime per step [s]", xscale=log10,
          yscale=log10, xticks=DECADES, yticks=(DECADES[1][1:2:end], DECADES[2][1:2:end]),
          title=panel_title("b", "Cost per step"))
rows = read_rows(dir)
ylims!(ax, minimum(r.time for r in rows) / 3, maximum(r.time for r in rows) * 3)
plot_methods!(ax, rows, METHODS, r -> r.d, r -> r.time)
for (k, linestyle, label) in [(1, :dash, "O(d)"), (3, :dot, "O(d³)")]
    d = [1e3, maximum(r.d for r in rows)]
    lines!(ax, d, 1e-1 .* (d ./ 1e3) .^ k; color=:black, linewidth=0.6, linestyle, label)
end
Legend(fig[1, 3], ax)
colsize!(fig.layout, 1, Relative(0.34))
save(joinpath(dir, "dim_scaling.pdf"), fig)
