# diagnostics/plot_fig4_japan.jl — plot from graphs/data/fig4_japan.csv
include(joinpath(@__DIR__, "common.jl"))
using CairoMakie, CSV, DataFrames

const FS = 11
const FIGDIR = joinpath(REPO, "Proposal&Paper", "graphs")
const DATADIR = joinpath(FIGDIR, "data")

const COL = Dict(
    "total" => RGBf(0.10, 0.10, 0.10),
    "from_australia" => RGBf(0.16, 0.47, 0.84),
    "non_australia" => RGBf(0.85, 0.37, 0.00),
)
const TITLES = Dict(
    "total" => "Total imports",
    "from_australia" => "From Australia",
    "non_australia" => "From other origins",
)
const ORDER = ["total", "from_australia", "non_australia"]

d = CSV.read(joinpath(DATADIR, "fig4_japan.csv"), DataFrame)
Hmax = maximum(d.h)
ymin = minimum(vcat(d.lo, d.point))
ymax = maximum(vcat(d.hi, d.point))
pad = 0.06 * (ymax - ymin)
yl = (ymin - pad, ymax + pad)

fig = Figure(size = (780, 320), fontsize = FS, backgroundcolor = :white)
for (k, key) in enumerate(ORDER)
    s = sort(d[d.series .== key, :], :h)
    col = COL[key]
    ax = Axis(fig[1, k];
              xlabel = "months after disruption",
              ylabel = k == 1 ? "response (log points)" : "",
              xticks = 0:2:Hmax)
    ax.xticklabelsize = FS; ax.yticklabelsize = FS
    ax.xlabelsize = FS; ax.ylabelsize = FS
    ax.spinewidth = 0.8
    ax.xgridvisible = false; ax.ygridvisible = false
    ylims!(ax, yl); xlims!(ax, -0.3, Hmax + 0.3)
    hlines!(ax, [0.0]; color = (:gray, 0.75), linestyle = :dash, linewidth = 1.0)

    band!(ax, s.h, s.lo, s.hi; color = (col, 0.22))
    lines!(ax, s.h, s.point; color = col, linewidth = 2.4)
    ok = isfinite.(s.point)
    scatter!(ax, s.h[ok], s.point[ok]; color = col, markersize = 7)
    i0 = findfirst(==(0), s.h)
    if i0 !== nothing && isfinite(s.point[i0])
        scatter!(ax, [0], [s.point[i0]]; color = col, markersize = 14, marker = :star5)
    end
    text!(ax, 0.03, 0.95; text = TITLES[key], space = :relative,
          align = (:left, :top), fontsize = FS, color = :black)
    k > 1 && hideydecorations!(ax; grid = false, label = false, ticklabels = false, ticks = false)
end
colgap!(fig.layout, 10)

pdf = joinpath(FIGDIR, "fig4_japan.pdf")
save(pdf, fig)
cp(pdf, joinpath(REPO, "Proposal&Paper", "figures", "fig4_japan.pdf"); force = true)
mkpath(joinpath(REPO, "Proposal&Paper", "figures", "data"))
cp(joinpath(DATADIR, "fig4_japan.csv"),
   joinpath(REPO, "Proposal&Paper", "figures", "data", "fig4_japan.csv"); force = true)
println("wrote $pdf")
