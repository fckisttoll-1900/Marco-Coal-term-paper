# diagnostics/task_C.jl — Leave-one-out on price (and Aus exports)
include(joinpath(@__DIR__, "common.jl"))

P = load_panel()
z0 = apply_drop!(P.z_w, P.ds, DROP_2022)   # baseline preferred
βp0, _, _, _ = lp_fs(P.lp, z0, P.dates, P.lq, P.lp, P.rea, 0)
βq0, _, _, _ = lp_fs(P.lq, z0, P.dates, P.lq, P.lp, P.rea, 0)

# C1: leave-one-block-out (12-month blocks on calendar)
T = length(P.dates)
nblk = cld(T, BLOCK)
c1 = DataFrame(block = Int[], block_start = String[], block_end = String[],
               beta_price = Float64[], delta_price = Float64[],
               n_obs = Int[])
for b in 1:nblk
    t0 = (b - 1) * BLOCK + 1
    t1 = min(b * BLOCK, T)
    keep = trues(T)
    keep[t0:t1] .= false
    y = copy(P.lp); q = copy(P.lq); p = copy(P.lp); r = copy(P.rea); z = copy(z0)
    y[.!keep] .= NaN; q[.!keep] .= NaN; p[.!keep] .= NaN; r[.!keep] .= NaN
    z[.!keep] .= 0.0
    β, _, n, _ = lp_fs(y, z, P.dates, q, p, r, 0)
    push!(c1, (b, string(P.ds[t0]), string(P.ds[t1]), β,
               isfinite(β) && isfinite(βp0) ? β - βp0 : NaN, n))
end
sort!(c1, :delta_price, by = abs, rev = true)
CSV.write(joinpath(DIAG_OUT, "C1_leave_one_block_price.csv"), c1)

# C2: leave-one-event-out on FULL event list (including 2022), price
ev = CSV.read(EVENTS_PATH, DataFrame)
ev.ds = ym.(ev.date)
zw_full = copy(P.z_w)
βp_full, _, _, _ = lp_fs(P.lp, zw_full, P.dates, P.lq, P.lp, P.rea, 0)

c2 = DataFrame(event_month = String[], event_name = String[], severity = String[],
               z_value = Float64[], beta_price = Float64[], delta_price = Float64[])
for r in eachrow(ev)
    m = string(r.ds)
    z = copy(zw_full)
    for i in eachindex(z)
        string(P.ds[i]) == m && (z[i] = 0.0)
    end
    β, _, _, _ = lp_fs(P.lp, z, P.dates, P.lq, P.lp, P.rea, 0)
    zv = 0.0
    for i in eachindex(zw_full)
        if string(P.ds[i]) == m
            zv = zw_full[i]; break
        end
    end
    push!(c2, (m, string(r.event), string(r.severity), zv, β,
               isfinite(β) && isfinite(βp_full) ? β - βp_full : NaN))
end
sort!(c2, :delta_price, by = abs, rev = true)
CSV.write(joinpath(DIAG_OUT, "C2_leave_one_event_price.csv"), c2)

# C3: leave-one-event-out for Aus exports
βq_full, _, _, _ = lp_fs(P.lq, zw_full, P.dates, P.lq, P.lp, P.rea, 0)
c3 = DataFrame(event_month = String[], event_name = String[], severity = String[],
               z_value = Float64[], beta_q = Float64[], delta_q = Float64[])
for r in eachrow(ev)
    m = string(r.ds)
    z = copy(zw_full)
    for i in eachindex(z)
        string(P.ds[i]) == m && (z[i] = 0.0)
    end
    β, _, _, _ = lp_fs(P.lq, z, P.dates, P.lq, P.lp, P.rea, 0)
    zv = 0.0
    for i in eachindex(zw_full)
        if string(P.ds[i]) == m
            zv = zw_full[i]; break
        end
    end
    push!(c3, (m, string(r.event), string(r.severity), zv, β,
               isfinite(β) && isfinite(βq_full) ? β - βq_full : NaN))
end
sort!(c3, :delta_q, by = abs, rev = true)
CSV.write(joinpath(DIAG_OUT, "C3_leave_one_event_qaus.csv"), c3)

open(joinpath(DIAG_OUT, "C_summary.txt"), "w") do io
    @printf(io, "baseline price h=0 (z_w drop2022): %.6g\n", βp0)
    @printf(io, "full z_w (incl 2022) price h=0: %.6g\n", βp_full)
    @printf(io, "baseline q_aus h=0 (z_w drop2022): %.6g\n", βq0)
    @printf(io, "full z_w q_aus h=0: %.6g\n", βq_full)
    @printf(io, "\nC1 top block by |delta|: start=%s beta=%.6g delta=%.6g\n",
            c1.block_start[1], c1.beta_price[1], c1.delta_price[1])
    @printf(io, "C2 top event by |delta| price: %s (%s) beta=%.6g delta=%.6g\n",
            c2.event_month[1], c2.event_name[1], c2.beta_price[1], c2.delta_price[1])
    @printf(io, "C3 top event by |delta| q_aus: %s (%s) beta=%.6g delta=%.6g\n",
            c3.event_month[1], c3.event_name[1], c3.beta_q[1], c3.delta_q[1])
end
println("wrote C1/C2/C3 CSVs and C_summary.txt")
println("script: diagnostics/task_C.jl")
