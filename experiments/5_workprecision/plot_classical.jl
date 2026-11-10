# Appendix: Figure 5 enlarged, with the classical solvers in the runtime panels only. Their
# steps are not comparable to those of the ODE filters.
# Usage: julia --project experiments/5_workprecision/plot_classical.jl
include("../../src/plotting.jl")

const METHODS = ["DiagonalEK1", "EK1", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF"]
const CLASSICAL_SOLVERS = ["Rodas5P", "RadauIIA5", "KenCarp4+GMRES"]
const PROBLEMS = [("burgers", "Burgers 1D (d=500)", METHODS),
                  ("fisher_kpp_2d", "Fisher-KPP 2D (d=4096)", filter(!=("EK1"), METHODS))]

# At rtol ≥ 1, KenCarp4 takes two steps and has an error of order one, so we drop these.
keep(r) = !(r.method == "KenCarp4+GMRES" && r.rtol >= 1)

set_theme!(paper_theme(size=(397, 250)))
fig = Figure()
for (i, (key, title, methods)) in enumerate(PROBLEMS)
    rows = filter(keep, read_rows(joinpath(@__DIR__, "data", key)))
    kw = (; xscale=log10, yscale=log10, xticks=DECADES,
          yticks=(DECADES[1][3:4:end], DECADES[2][3:4:end]))  # 10⁻¹², 10⁻⁸, 10⁻⁴, ...
    ax_steps = Axis(fig[i, 1]; xlabel="# steps", ylabel="Final error",
                    title=panel_title("ab"[i], title), kw...)
    ax_time = Axis(fig[i, 2]; xlabel="Runtime [s]", kw...)
    linkyaxes!(ax_steps, ax_time)
    hideydecorations!(ax_time; grid=false)
    plot_methods!(ax_steps, rows, methods, nsteps, r -> r.error; markerscale=0.6)
    plot_methods!(ax_time, rows, methods, r -> r.time, r -> r.error; markerscale=0.6)
    plot_methods!(ax_time, rows, CLASSICAL_SOLVERS, r -> r.time, r -> r.error)
    i == 1 && Legend(fig[1:2, 3], ax_time)
end
save(joinpath(@__DIR__, "data", "workprecision_classical.pdf"), fig)
