# diagnostics/make_fig4_japan.jl
# Japan importer IRFs: total / from-Aus / non-Aus, keep-2022 harmonised LP.
include(joinpath(@__DIR__, "harm_lp.jl"))
using CairoMakie

const NBOOT = 10_000
const Hmax = 12
const FS = 11
const FIGDIR = joinpath(REPO, "Proposal&Paper", "graphs")
const DATADIR = joinpath(FIGDIR, "data")
mkpath(DATADIR)

# Original 06_figures palette
const COL_TOTAL = RGBf(0.10, 0.10, 0.10)          # World / total — black
const COL_AUS   = RGBf(0.16, 0.47, 0.84)          # Australia — blue
const COL_OTHER = RGBf(0.85, 0.37, 0.00)          # non-Aus — burnt orange (distinct in greyscale)

P = load_panel()
z_keep = z_baseline(P)
T = length(P.dates)
qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))

jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
function jpull(part)
    sub = jp[jp.partner .== part, :]
    m = Dict(string(r.ds) => Float64(r.tonnes) for r in eachrow(sub))
    [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
end
jw = jpull("World"); ja = jpull("Australia")
jn = [isfinite(jw[i]) && isfinite(ja[i]) ? jw[i] - ja[i] : NaN for i in 1:T]
series = [
    ("total", "Total imports", COL_TOTAL,
     [isfinite(x) && x > 0 ? log(x) : NaN for x in jw]),
    ("from_australia", "From Australia", COL_AUS,
     [isfinite(x) && x > 0 ? log(x) : NaN for x in ja]),
    ("non_australia", "From other origins", COL_OTHER,
     [isfinite(x) && x > 0 ? log(x) : NaN for x in jn]),
]

println("fig4_japan: point paths + MBB $NBOOT ...")
paths = Dict{String, Vector{Float64}}()
bands = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
for (key, _, _, y) in series
    paths[key] = lp_harm_path(y, z_keep, P.rea, P.dates, Hmax)
end

rng = MersenneTwister(BOOT_SEED)
BS = Dict(key => fill(NaN, NBOOT, Hmax + 1) for (key, _, _, _) in series)
ys = Dict(key => y for (key, _, _, y) in series)
for b in 1:NBOOT
    ii = block_idx(T, BLOCK, rng)
    zb = z_keep[ii]; rb = P.rea[ii]; db = P.dates[ii]
    for (key, _, _, _) in series
        BS[key][b, :] = lp_harm_path(ys[key][ii], zb, rb, db, Hmax)
    end
    b % 2000 == 0 && @printf("  draw %d/%d\n", b, NBOOT)
end
for (key, _, _, _) in series
    lo = [qband(BS[key][:, h + 1])[1] for h in 0:Hmax]
    hi = [qband(BS[key][:, h + 1])[2] for h in 0:Hmax]
    bands[key] = (lo, hi)
end

# CSV
out = DataFrame(series = String[], h = Int[], point = Float64[],
                lo = Float64[], hi = Float64[])
for (key, _, _, _) in series
    lo, hi = bands[key]
    for h in 0:Hmax
        push!(out, (key, h, paths[key][h + 1], lo[h + 1], hi[h + 1]))
    end
end
CSV.write(joinpath(DATADIR, "fig4_japan.csv"), out)
println("  wrote graphs/data/fig4_japan.csv")

# Shared y-axis
ymin = minimum(vcat([bands[k][1] for (k, _, _, _) in series]...))
ymax = maximum(vcat([bands[k][2] for (k, _, _, _) in series]...))
for (key, _, _, _) in series
    global ymin = min(ymin, minimum(filter(isfinite, paths[key])))
    global ymax = max(ymax, maximum(filter(isfinite, paths[key])))
end
pad = 0.06 * (ymax - ymin)
yl = (ymin - pad, ymax + pad)

hs = 0:Hmax
fig = Figure(size = (780, 320), fontsize = FS, backgroundcolor = :white)
for (k, (key, title, col, _)) in enumerate(series)
    ax = Axis(fig[1, k];
              xlabel = "months after disruption",
              ylabel = k == 1 ? "response (log points)" : "",
              xticks = 0:2:Hmax)
    ax.xticklabelsize = FS; ax.yticklabelsize = FS
    ax.xlabelsize = FS; ax.ylabelsize = FS
    ax.spinewidth = 0.8
    ax.xgridvisible = false; ax.ygridvisible = false
    ylims!(ax, yl)
    xlims!(ax, -0.3, Hmax + 0.3)
    hlines!(ax, [0.0]; color = (:gray, 0.75), linestyle = :dash, linewidth = 1.0)

    lo, hi = bands[key]
    y = paths[key]
    band!(ax, hs, lo, hi; color = (col, 0.22))
    lines!(ax, hs, y; color = col, linewidth = 2.4)
    ok = isfinite.(y)
    scatter!(ax, collect(hs)[ok], y[ok]; color = col, markersize = 7)
    # star at h = 0
    if isfinite(y[1])
        scatter!(ax, [0], [y[1]]; color = col, markersize = 14, marker = :star5)
    end
    text!(ax, 0.03, 0.95; text = title, space = :relative,
          align = (:left, :top), fontsize = FS, color = :black)
    k > 1 && hideydecorations!(ax; grid = false, label = false, ticklabels = false, ticks = false)
end
colgap!(fig.layout, 10)

pdf = joinpath(FIGDIR, "fig4_japan.pdf")
save(pdf, fig)
# mirror
cp(pdf, joinpath(REPO, "Proposal&Paper", "figures", "fig4_japan.pdf"); force = true)
mkpath(joinpath(REPO, "Proposal&Paper", "figures", "data"))
cp(joinpath(DATADIR, "fig4_japan.csv"),
   joinpath(REPO, "Proposal&Paper", "figures", "data", "fig4_japan.csv"); force = true)
println("  wrote $pdf")
println("script: diagnostics/make_fig4_japan.jl")
