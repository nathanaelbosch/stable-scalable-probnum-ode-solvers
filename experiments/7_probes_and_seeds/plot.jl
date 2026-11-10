# Appendix: final error, number of steps and runtime of the MatfreeEK1 against the number of
# probes m, on the 1D Brusselator. Median and range over 10 seeds, one shade per tolerance and
# one row per stiffness α. To the right are the EK1 with block-diagonal covariance (the limit
# m → ∞) and the EK1.
# Usage: julia --project experiments/7_probes_and_seeds/plot.jl
include("../../src/plotting.jl")

const MS = [4, 6, 8, 16, 32, 64, 128, 256, 512]
const RTOLS = [1e-3, 1e-5]
colors = cgrad(:Oranges, 4; categorical=true)[[2, 4]]

set_theme!(paper_theme(size=(397, 200)))
fig = Figure()
all_pairs = []
for (i, alpha) in enumerate(["1", "10"])
    rows = filter(r -> r.status == "Success", CSV.File(joinpath(@__DIR__, "data", "alpha$alpha.csv")))
    xlabel = i == 2 ? "Number of samples m" : ""
    kw = (; yscale=log10, yticks=DECADES, xlabel)
    panels = [(probes_axes(fig[i, 1], MS; ylabel="Final error",
                           title=panel_title("ad"[i], "α = $alpha, error"), kw...), r -> r.error),
              (probes_axes(fig[i, 2], MS; ylabel="# steps",
                           title=panel_title("be"[i], "α = $alpha, steps"), kw...), nsteps),
              (probes_axes(fig[i, 3], MS; ylabel="Runtime [s]",
                           title=panel_title("cf"[i], "α = $alpha, runtime"), kw...), r -> r.time)]
    for (k, rtol) in enumerate(RTOLS)
        at(method, m) = filter(r -> r.method == method && r.rtol == rtol && r.m == m, rows)
        for ((ax, ax_ref), y) in panels
            vals = [y.(at("MatfreeEK1", m)) for m in MS]
            band!(ax, MS, minimum.(vals), maximum.(vals); color=(colors[k], 0.25))
            scatterlines!(ax, MS, median.(vals); color=colors[k], marker=:star5, markersize=6,
                          label=rich(TAU_REL, " = 10" * superscript(round(Int, log10(rtol)))))
            scatter!(ax_ref, [1], y.(at("EK1, block-diagonal", 0)); TRUNCATED_MARKER...,
                     strokecolor=colors[k])
            scatter!(ax_ref, [2], y.(at("EK1", 0)); EK1_MARKER..., color=colors[k])
        end
    end
    append!(all_pairs, first.(panels))
    i == 1 && probes_legend(fig[1:2, 4], panels[1][1][1])
end
finish_probes_axes!(fig, all_pairs)
save(joinpath(@__DIR__, "data", "probes_and_seeds.pdf"), fig)

