# Appendix: final error and χ²/d of the MatfreeEK1 against the number of probes m, on
# Lorenz96 at fixed steps. Median and range over 10 seeds, one shade per step size. To the
# right are the EK1 with block-diagonal covariance (the limit m → ∞) and the EK1.
# Usage: julia --project experiments/6_calibration/plot_probes.jl
include("../../src/plotting.jl")

rows = read_rows(joinpath(@__DIR__, "data"))
const MS = [4, 6, 8, 16, 32, 64, 128, 256, 512]
dts = sort(unique(r.dt for r in rows); rev=true)
colors = cgrad(:Purples, length(dts) + 1; categorical=true)[2:end]
select(method, dt) = filter(r -> r.method == method && r.dt == dt, rows)
dtlabel(dt) = Dict(1e-2 => "10⁻²", 5e-3 => "5·10⁻³", 2e-3 => "2·10⁻³")[dt]

set_theme!(paper_theme(size=(397, 105)))
fig = Figure()
all_pairs = []
for (k, (key, ylabel)) in enumerate(((:error, "Final error"), (:chi2_d, "χ²/d")))
    ax, ax_ref = probes_axes(fig[1, k], MS; xlabel="Number of samples m", ylabel,
                             yscale=log10, yticks=DECADES, title=panel_title("ab"[k], ylabel))
    push!(all_pairs, (ax, ax_ref))
    for (i, dt) in enumerate(dts)
        vals = [getproperty.(select("MatfreeEK1, m=$m", dt), key) for m in MS]
        band!(ax, MS, minimum.(vals), maximum.(vals); color=(colors[i], 0.25))
        scatterlines!(ax, MS, median.(vals); color=colors[i], marker=:star5, markersize=6,
                      label="Δt = " * dtlabel(dt))
        scatter!(ax_ref, [1], getproperty.(select("EK1, block-diagonal", dt), key);
                 TRUNCATED_MARKER..., strokecolor=colors[i])
        scatter!(ax_ref, [2], getproperty.(select("EK1", dt), key); EK1_MARKER...,
                 color=colors[i])
    end
end
probes_legend(fig[1, 3], all_pairs[1][1])
finish_probes_axes!(fig, all_pairs)
save(joinpath(@__DIR__, "data", "calibration_probes.pdf"), fig)
