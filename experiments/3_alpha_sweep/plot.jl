# Figure 3: number of steps against α (a), and work-precision diagrams at α = 0.02 (b)
# and α = 10 (c).
# Usage: julia --project experiments/3_alpha_sweep/plot.jl
include("../../src/plotting.jl")

const METHODS = ["EK0", "DiagonalEK1", "EK1", "ExpEK", "DiagonalEK1+IEKF", "MatfreeEK1",
                 "MatfreeEK1+IEKF"]
const ALPHAS = ["0.02", "0.1", "0.3", "1", "3", "10"]
rows(alpha) = read_rows(joinpath(@__DIR__, "data", "alpha$alpha"))

set_theme!(paper_theme(size=(397, 94)))
fig = Figure()

ax_sweep = Axis(fig[1, 1]; xlabel="Stiffness coefficient α", ylabel="# steps",
                xscale=log10, yscale=log10, xticks=(parse.(Float64, ALPHAS), ALPHAS),
                title=panel_title("a", "Stiffness sweep"))
sweep_rows = [(; r..., alpha=parse(Float64, a)) for a in ALPHAS for r in rows(a) if r.rtol == 1e-3]
plot_methods!(ax_sweep, sweep_rows, METHODS, r -> r.alpha, nsteps)

axes = map(enumerate(["0.02", "10"])) do (i, alpha)
    ax = Axis(fig[1, 1+i]; xlabel="# steps", ylabel=i == 1 ? "Final error" : "",
              xscale=log10, yscale=log10, title=panel_title("bc"[i], "α = $alpha"))
    plot_methods!(ax, rows(alpha), METHODS, nsteps, r -> r.error)
    ax
end
linkyaxes!(axes...)
hideydecorations!(axes[2]; grid=false)
Legend(fig[1, 4], ax_sweep)
save(joinpath(@__DIR__, "data", "alpha_sweep.pdf"), fig)
