# Appendix: calibration on Lorenz96 at fixed steps. χ²/d of the final state against the step
# size, with the median and the range over 10 seeds for the MatfreeEK1 variants at the default
# m = 8.
# Usage: julia --project experiments/6_calibration/plot.jl
include("../../src/plotting.jl")

rows = read_rows(joinpath(@__DIR__, "data"))
dts = sort(unique(r.dt for r in rows))
select(method, dt) = filter(r -> r.method == method && r.dt == dt, rows)

set_theme!(paper_theme(size=(260, 110)))
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Step size Δt", ylabel="χ²/d", xscale=log10, yscale=log10,
          xticks=(dts, ["2·10⁻³", "5·10⁻³", "10⁻²"]), yticks=DECADES)
for (method, key) in [("EK0", "EK0"), ("DiagonalEK1", "DiagonalEK1"), ("EK1", "EK1"),
                      ("DiagonalEK1+IEKF", "DiagonalEK1+IEKF"),
                      ("MatfreeEK1", "MatfreeEK1, m=8"), ("MatfreeEK1+IEKF", "MatfreeEK1+IEKF, m=8")]
    vals = [[r.chi2_d for r in select(key, dt)] for dt in dts]
    style = STYLES[method]
    length(vals[1]) > 1 && band!(ax, dts, minimum.(vals), maximum.(vals); color=(style.color, 0.25))
    scatterlines!(ax, dts, median.(vals); label=method, style..., color=(style.color, 0.6),
                  markercolor=style.color)
end
hlines!(ax, 1; color=:black, linewidth=0.6, linestyle=:dot, label="Perfect calibration")
Legend(fig[1, 2], ax)
save(joinpath(@__DIR__, "data", "calibration.pdf"), fig)
