# Figure 1: the FitzHugh-Nagumo solution (a), and the runtime when the domain grows (b), the
# problem gets stiffer (c), and the discretization gets finer (d).
# Usage: julia --project experiments/1_fhn/plot.jl
include("../../src/plotting.jl")
using DelimitedFiles

const METHODS = ["EK0", "DiagonalEK1", "EK1", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF",
                 "KenCarp4+GMRES"]
dir = joinpath(@__DIR__, "data")
rows(sweep) = read_rows(joinpath(dir, sweep))

set_theme!(paper_theme(size=(397, 103)))
fig = Figure()

gl = fig[1, 1] = GridLayout()
for (i, name) in enumerate(["u", "v"])
    U = readdlm(joinpath(dir, "solution_$name.txt"), ',')
    lim = maximum(abs, U)
    ax = Axis(gl[i, 1]; aspect=DataAspect(), ylabel=name, ylabelrotation=0,
              title=i == 1 ? panel_title("a", "PDE") : "")
    hidedecorations!(ax; label=false)
    heatmap!(ax, U; colormap=:RdBu, colorrange=(-lim, lim))
end

kw = (; xscale=log10, yscale=log10, xticks=DECADES, yticks=DECADES)
ax_domain = Axis(fig[1, 2]; xlabel="Dimension d", ylabel="Runtime [s]",
                 title=panel_title("b", "Domain growth"), kw...)
ax_stiff = Axis(fig[1, 3]; xlabel="Stiffness λ", title=panel_title("c", "Stiffening"), kw...)
ax_refine = Axis(fig[1, 4]; xlabel="Resolution 1/Δx", title=panel_title("d", "Refinement"), kw...,
                 yticks=(DECADES[1][2:2:end], DECADES[2][2:2:end]))

# (b) with O(d) and O(d³) slopes. O(d) starts at the fastest method and O(d³) at the EK1.
domain = rows("domain")
plot_methods!(ax_domain, domain, METHODS, r -> r.d, r -> r.time)
d = collect(extrema(r.d for r in domain))
t_fast = minimum(r.time for r in domain)
t_slow = only(r.time for r in domain if r.d == d[1] && r.method == "EK1")
for (k, t, linestyle, label) in [(1, t_fast, :dash, "O(d)"), (3, t_slow, :dot, "O(d³)")]
    lines!(ax_domain, d, t .* (d ./ d[1]) .^ k; color=:black, linewidth=0.6, linestyle,
           label, yautolimits=false)
end
plot_methods!(ax_stiff, rows("stiffness"), METHODS, r -> 0.05 / r.dx^2, r -> r.time)

# (d) with a least-squares line in log-log space through each method's last four points,
# extrapolated to 10⁷ with its 95% confidence band
refinement = rows("refinement")
plot_methods!(ax_refine, refinement, METHODS, r -> r.d, r -> r.time)
for method in METHODS
    rs = sort(filter(r -> r.method == method, refinement); by=r -> r.d)
    lx, ly = log10.([r.d for r in rs[end-3:end]]), log10.([r.time for r in rs[end-3:end]])
    Sxx = sum(abs2, lx .- mean(lx))
    slope = sum((lx .- mean(lx)) .* (ly .- mean(ly))) / Sxx
    fit(x) = mean(ly) + slope * (x - mean(lx))
    sigma = sqrt(sum(abs2, ly .- fit.(lx)) / (length(lx) - 2))
    xs = range(lx[end], 7, length=50)
    halfwidth = 4.303 .* sigma .* sqrt.(1 / length(lx) .+ (xs .- mean(lx)) .^ 2 ./ Sxx)  # t-quantile, 2 dof
    color = STYLES[method].color
    band!(ax_refine, 10 .^ xs, 10 .^ (fit.(xs) .- halfwidth), 10 .^ (fit.(xs) .+ halfwidth);
          color=(color, 0.15), yautolimits=false)
    lines!(ax_refine, 10 .^ xs, 10 .^ fit.(xs); color, linestyle=:dash, linewidth=0.6,
           yautolimits=false)
end
xlims!(ax_refine, 1e2, 1e7)
ylims!(ax_refine, nothing, 1e10)

Legend(fig[1, 5], ax_domain)
save(joinpath(dir, "fhn.pdf"), fig)
save(joinpath(dir, "fhn.png"), fig; px_per_unit=4)  # for the README
