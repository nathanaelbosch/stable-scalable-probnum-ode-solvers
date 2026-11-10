# Figure 5: final error against the number of steps and the runtime. The classical
# matrix-free KenCarp4+GMRES is a reference in the runtime panels only, since its steps are
# not comparable to those of the ODE filters.
# Usage: julia --project experiments/5_workprecision/plot.jl
include("../../src/plotting.jl")

const METHODS = ["DiagonalEK1", "EK1", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF"]
# The single EK1 solve on Fisher-KPP only gives a runtime for the text
const PROBLEMS = [("burgers", "Burgers 1D (d=500)", METHODS),
                  ("fisher_kpp_2d", "Fisher-KPP 2D (d=4096)", filter(!=("EK1"), METHODS))]

const REFERENCE = "KenCarp4+GMRES"
# At rtol ≥ 1, KenCarp4 takes two steps and has an error of order one, so we drop these.
keep(r) = !(r.method == REFERENCE && r.rtol >= 1)

set_theme!(paper_theme(size=(397, 92)))
fig = Figure()
for (i, (key, title, methods)) in enumerate(PROBLEMS)
    rows = filter(keep, read_rows(joinpath(@__DIR__, "data", key)))
    kw = (; xscale=log10, yscale=log10, xticks=DECADES,
          yticks=(DECADES[1][3:4:end], DECADES[2][3:4:end]))  # 10⁻¹², 10⁻⁸, 10⁻⁴, ...
    ax_steps = Axis(fig[1, 2i-1]; xlabel="# steps", ylabel=i == 1 ? "Final error" : "",
                    title=panel_title("ab"[i], title), kw...)
    ax_time = Axis(fig[1, 2i]; xlabel="Runtime [s]", kw...)
    linkyaxes!(ax_steps, ax_time)
    hideydecorations!(ax_time; grid=false)
    plot_methods!(ax_steps, rows, methods, nsteps, r -> r.error)
    plot_methods!(ax_time, rows, [methods; REFERENCE], r -> r.time, r -> r.error)
    i == 1 && Legend(fig[1, 5], ax_time)
end
save(joinpath(@__DIR__, "data", "workprecision.pdf"), fig)
