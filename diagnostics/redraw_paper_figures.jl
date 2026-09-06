# diagnostics/redraw_paper_figures.jl — redraw F1–F4 from saved CSVs (no MBB)
include(joinpath(@__DIR__, "common.jl"))
using CairoMakie, CSV, DataFrames, Statistics, Printf

const FIGDIR = joinpath(REPO, "Proposal&Paper", "graphs")
const FS = 11
const COL_AUS  = RGBf(0.16, 0.47, 0.84)
const COL_WEXC = RGBf(0.92, 0.41, 0.20)
const COL_PRICE = RGBf(0.20, 0.20, 0.20)
const COL_NOREA = RGBf(0.55, 0.55, 0.52)
const COL_REA   = RGBf(0.16, 0.47, 0.84)
const COL_ACTUAL = RGBf(0.80, 0.15, 0.15)

function style_axis!(ax)
    ax.xticklabelsize = FS; ax.yticklabelsize = FS
    ax.xlabelsize = FS; ax.ylabelsize = FS
    ax.spinewidth = 0.8
    ax.xgridvisible = false; ax.ygridvisible = false
end

# ---- F1 ----
d1 = CSV.read(joinpath(DIAG_OUT, "fig_data_F1_keep2022.csv"), DataFrame)
hs = d1.h
ymin = min(minimum(d1.aus_lo), minimum(d1.wexc_lo))
ymax = max(maximum(d1.aus_hi), maximum(d1.wexc_hi))
pad = 0.05 * (ymax - ymin)
yl = (ymin - pad, ymax + pad)

fig1 = Figure(size = (480, 520), fontsize = FS, backgroundcolor = :white)
ga = fig1[1, 1] = GridLayout()
ax1 = Axis(ga[1, 1]; ylabel = "Mt per month")
ax2 = Axis(ga[2, 1]; ylabel = "Mt per month", xlabel = "Horizon h (months)")
for ax in (ax1, ax2)
    style_axis!(ax); ylims!(ax, yl); xlims!(ax, -0.3, 12.3)
    hlines!(ax, [0.0]; color = (:gray, 0.85), linestyle = :dash, linewidth = 1.0)
end
hidexdecorations!(ax1; grid = false)
band!(ax1, hs, d1.aus_lo, d1.aus_hi; color = (COL_AUS, 0.22))
lines!(ax1, hs, d1.aus; color = COL_AUS, linewidth = 2.0)
scatter!(ax1, hs, d1.aus; color = COL_AUS, markersize = 6)
text!(ax1, 0.02, 0.95; text = "(a) Australia", space = :relative,
      align = (:left, :top), fontsize = FS)
band!(ax2, hs, d1.wexc_lo, d1.wexc_hi; color = (COL_WEXC, 0.22))
lines!(ax2, hs, d1.wexc; color = COL_WEXC, linewidth = 2.0)
scatter!(ax2, hs, d1.wexc; color = COL_WEXC, markersize = 6)
text!(ax2, 0.02, 0.95; text = "(b) Non-Australian aggregate  Wexc", space = :relative,
      align = (:left, :top), fontsize = FS)
rowgap!(ga, 8)
save(joinpath(FIGDIR, "fig1_supply_and_replacement.pdf"), fig1)

# ---- F2 ----
draws = CSV.read(joinpath(DIAG_OUT, "H_placebo_draws.csv"), DataFrame)
β_actual = -0.11633394529086633
cq = filter(isfinite, Float64.(draws.beta_q))
p_recomp = mean(abs.(cq) .>= abs(β_actual))
@printf("F2 recomputed p = %.6f\n", p_recomp)
if !(0.0024 <= p_recomp <= 0.0038)
    println("WARNING: recomputed randomisation p = $p_recomp (expected ≈ 0.003)")
end
fig2 = Figure(size = (480, 340), fontsize = FS, backgroundcolor = :white)
ax = Axis(fig2[1, 1]; xlabel = "Placebo β (Australian exports, h = 0, log)", ylabel = "Count")
style_axis!(ax)
hist!(ax, cq; bins = 40, color = (COL_AUS, 0.45), strokewidth = 0)
vlines!(ax, [β_actual]; color = COL_ACTUAL, linewidth = 2.2)
text!(ax, 0.03, 0.92; text = @sprintf("actual = %+.4f\nrand. p = %.3f", β_actual, p_recomp),
      space = :relative, align = (:left, :top), fontsize = FS, color = COL_ACTUAL)
save(joinpath(FIGDIR, "fig2_randomisation.pdf"), fig2)

# ---- F3 ----
d3 = CSV.read(joinpath(DIAG_OUT, "fig_data_F3_price_keep2022.csv"), DataFrame)
s = 0.18; Δ = 0.116
ε_style = Dict(0.1 => :solid, 0.2 => :dash, 0.3 => :dot)
fig3 = Figure(size = (480, 340), fontsize = FS, backgroundcolor = :white)
ax = Axis(fig3[1, 1]; xlabel = "Horizon h (months)", ylabel = "Log points")
style_axis!(ax); xlims!(ax, -0.2, 6.2)
hlines!(ax, [0.0]; color = (:gray, 0.85), linestyle = :dash, linewidth = 1.0)
band!(ax, d3.h, d3.lo, d3.hi; color = (COL_PRICE, 0.18))
lines!(ax, d3.h, d3.price; color = COL_PRICE, linewidth = 2.0)
scatter!(ax, d3.h, d3.price; color = COL_PRICE, markersize = 7)
for ε in (0.1, 0.2, 0.3)
    hlines!(ax, [(1 - 0.0) * s * Δ / ε]; color = (RGBf(0.15, 0.45, 0.25), 0.95),
            linestyle = ε_style[ε], linewidth = 1.4)
    hlines!(ax, [(1 - 0.73) * s * Δ / ε]; color = (RGBf(0.15, 0.45, 0.25), 0.40),
            linestyle = ε_style[ε], linewidth = 1.4)
end
text!(ax, 0.98, 0.98;
      text = "solid r = 0; faded r = 0.73\nline style |ε| = 0.1 / 0.2 / 0.3",
      space = :relative, align = (:right, :top), fontsize = FS - 1)
save(joinpath(FIGDIR, "fig3_price_power.pdf"), fig3)

# ---- F4 from fig_data_F4_artifact.txt parsed manually ----
# numbers from previous run
idn_nr = (2.0211, -2.1790, 6.6325)
idn_r  = (-0.5938, -4.6719, 2.1642)
w_nr   = (1.4306, -2.6682, 6.6010)
w_r    = (-0.5257, -6.1381, 5.7784)
lk_nr  = (0.9265, -2.2948, 4.8483)
lk_r   = (-0.1601, -2.6419, 2.6152)

fig4 = Figure(size = (520, 360), fontsize = FS, backgroundcolor = :white)
ax = Axis(fig4[1, 1];
          yticks = (1:3, ["Impact leakage rate", "Wexc impact (Mt)", "Indonesia impact (Mt)"]))
style_axis!(ax)
vlines!(ax, [0.0]; color = (:gray, 0.85), linestyle = :dash, linewidth = 1.0)
specs = [(3, idn_nr, idn_r), (2, w_nr, w_r), (1, lk_nr, lk_r)]
for (y, nore, rea) in specs
    lines!(ax, [nore[2], nore[3]], [y + 0.12, y + 0.12]; color = COL_NOREA, linewidth = 2.0)
    scatter!(ax, [nore[1]], [y + 0.12]; color = COL_NOREA, markersize = 11, marker = :circle)
    lines!(ax, [rea[2], rea[3]], [y - 0.12, y - 0.12]; color = COL_REA, linewidth = 2.0)
    scatter!(ax, [rea[1]], [y - 0.12]; color = COL_REA, markersize = 11, marker = :diamond)
end
# manual legend markers
scatter!(ax, [NaN], [NaN]; color = COL_NOREA, marker = :circle, markersize = 11, label = "No REA")
scatter!(ax, [NaN], [NaN]; color = COL_REA, marker = :diamond, markersize = 11, label = "With REA")
ylims!(ax, 0.4, 3.6)
axislegend(ax; position = :rb, labelsize = FS - 1, framevisible = false)
save(joinpath(FIGDIR, "fig4_specification_artifact.pdf"), fig4)

figdir2 = joinpath(REPO, "Proposal&Paper", "figures")
for f in readdir(FIGDIR)
    endswith(f, ".pdf") && cp(joinpath(FIGDIR, f), joinpath(figdir2, f); force = true)
end
println("redrawn PDFs in $FIGDIR")
