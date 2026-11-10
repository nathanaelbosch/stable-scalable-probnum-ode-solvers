using CairoMakie, TuePlots, CSV, Tables, Statistics

# The legend does not set the height of its row, so the axes keep the full figure height.
# `size` is in points. The NeurIPS text width is 397.
paper_theme(; size) = merge(
    Theme(; size, figure_padding=(2, 6, 2, 2), Axis=(; titlealign=:left),
          ScatterLines=(; strokewidth=0.2),
          Legend=(; framevisible=false, labelsize=7, patchsize=(12, 8), rowgap=0, padding=2,
                  tellheight=false)),
    Theme(TuePlots.SETTINGS[:NEURIPS]))

const STYLES = Dict(
    "EK0" => (color=:gray80, marker=:circle, markersize=5),
    "DiagonalEK1" => (color=:gray60, marker=:diamond, markersize=5),
    "EK1" => (color=:gray40, marker=:rect, markersize=5),
    "ExpEK" => (color=:black, marker=:dtriangle, markersize=5),
    "DiagonalEK1+IEKF" => (color="#1F78B4", marker=:utriangle, markersize=8),
    "MatfreeEK1" => (color="#FDBF6F", marker=:star5, markersize=8),
    "MatfreeEK1+IEKF" => (color="#FF7F00", marker=:star8, linestyle=:dash, markersize=8),
    "Rodas5P" => (color=:gray60, marker=:pentagon, linestyle=:dot, markersize=4),
    "RadauIIA5" => (color=:gray30, marker=:hexagon, linestyle=:dot, markersize=4),
    "KenCarp4+GMRES" => (color=:gray45, marker=:xcross, linestyle=:dot, markersize=4,
                         markercolor=:white, strokecolor=:gray45, strokewidth=0.5),
)

# The rows of the CSV files in `dir`, without those of failed solves
read_rows(dir) = [r for f in readdir(dir; join=true) if endswith(f, ".csv")
                  for r in Tables.rowtable(CSV.File(f)) if get(r, :status, "Success") == "Success"]

nsteps(r) = r.naccept + r.nreject

# Log-axis ticks at the integer powers of 10 only
const SUPERSCRIPTS = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
superscript(i) = map(c -> SUPERSCRIPTS[c], string(i))
const DECADES = (10.0 .^ (-14:14), ["10" * superscript(i) for i in -14:14])

panel_title(letter, title) = rich(rich("$letter.", font=:bold), " ", title)

const TAU_REL = rich("τ", subscript("rel"))

# Plot `y(row)` against `x(row)` for each method in `methods`
function plot_methods!(ax, rows, methods, x, y; markerscale=1)
    for method in methods
        ok = filter(r -> r.method == method, rows)
        isempty(ok) && continue
        style = STYLES[method]
        scatterlines!(ax, x.(ok), y.(ok); label=method,
                      merge(style, (; color=(style.color, 0.6),
                                    markercolor=get(style, :markercolor, style.color),
                                    markersize=markerscale * style.markersize))...)
    end
end

# Axes for plots over the number of samples m. A log-scale axis for m, then a break in the
# x-axis and a narrow axis with two references. The EK1 with block-diagonal covariance (m = ∞)
# is at x = 1 and the EK1 at x = 2. Call `finish_probes_axes!` once the figure is complete.
const PROBES_MARGIN = 1.4  # space before the first and after the last m, as a factor in m
const DEFAULT_LINE = (; color=:gray50, linewidth=0.8, linestyle=:dash)
# Markers of the two references on the narrow axis
const TRUNCATED_MARKER = (; marker=:star5, markersize=7, color=:white, strokewidth=0.8)
const EK1_MARKER = (; marker=:rect, markersize=5, strokewidth=0.2, strokecolor=:black)

function probes_axes(pos, ms; xlabel, title, kwargs...)
    gl = GridLayout(pos)
    ticks = [4, 16, 64, 512]
    ax = Axis(gl[1, 1]; xscale=log10, xticks=(ticks, string.(ticks)), xlabel, title,
              limits=((ms[1] / PROBES_MARGIN, ms[end] * PROBES_MARGIN), nothing),
              rightspinevisible=false, kwargs...)
    ax_ref = Axis(gl[1, 2]; xticks=([1, 2], ["∞", "EK1"]), leftspinevisible=false,
                  xgridvisible=false, width=18, kwargs...)
    linkyaxes!(ax, ax_ref)
    hideydecorations!(ax_ref; grid=false)
    colgap!(gl, 6)
    vlines!(ax, 8; DEFAULT_LINE...)
    return ax, ax_ref
end

# Center the break between the last m and ∞ and draw its marks. This leaves `gap` points
# after the last m, before ∞ and after the EK1. Needs the final layout.
function finish_probes_axes!(fig, pairs)
    gap = 3
    Makie.update_state_before_display!(fig)
    for (ax, ax_ref) in pairs
        r1, r2 = ax.scene.viewport[], ax_ref.scene.viewport[]
        (lo, hi) = ax.limits[][1]
        m1, m2 = lo * PROBES_MARGIN, hi / PROBES_MARGIN  # the first and the last m
        # A margin of ε in log(m) on both sides is `gap` points wide
        ε = gap * log(m2 / m1) / (r1.widths[1] - 2gap)
        xlims!(ax, m1 * exp(-ε), m2 * exp(ε))
        f = gap / r2.widths[1]  # pads of h around [1, 2]: h / (1 + 2h) = f
        h = f / (1 - 2f)
        xlims!(ax_ref, 1 - h, 2 + h)
        xs = (r1.origin[1] + r1.widths[1], r2.origin[1])
        for x in xs, y in (r1.origin[2], r1.origin[2] + r1.widths[2])
            lines!(fig.scene, [Point2f(x - 1.5, y - 2), Point2f(x + 1.5, y + 2)];
                   color=:black, linewidth=0.6)
        end
    end
end

# The legend of the plots of `ax`, followed by the methods and the default m
function probes_legend(pos, ax)
    plots, labels = Makie.get_labeled_plots(ax; merge=false, unique=false)
    extra = [MarkerElement(; marker=:star5, markersize=6, color=:gray30, strokewidth=0.2,
                           strokecolor=:black),
             MarkerElement(; TRUNCATED_MARKER..., strokecolor=:gray30),
             MarkerElement(; EK1_MARKER..., color=:gray30),
             LineElement(; DEFAULT_LINE...)]
    Legend(pos, [plots; extra],
           [labels; ["MatfreeEK1", "Block-diagonal EK1", "EK1", "Default m = 2(q+1)"]])
end
