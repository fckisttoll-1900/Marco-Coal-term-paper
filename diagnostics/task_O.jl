# diagnostics/task_O.jl — Data gaps
include(joinpath(@__DIR__, "common.jl"))

P = load_panel()

function gaps_in(ds_have, span_start, span_end)
    want = [Dates.format(d, "yyyy-mm")
            for d in Date(span_start * "-01"):Month(1):Date(span_end * "-01")]
    have = Set(string.(ds_have))
    return [w for w in want if !(w in have)]
end

# O1: every series entering LPs
series_info = [
    ("q_aus / aus_exports", String.(P.ds[isfinite.(P.lq)]), "panel"),
    ("p_real", String.(P.ds[isfinite.(P.lp)]), "panel"),
    ("rea", String.(P.ds[isfinite.(P.rea)]), "panel"),
]

# add exporters and japan from raw
for (f, lab) in [("aus_exports.csv", "aus_exports_raw"),
                 ("exp_idn.csv", "Indonesia"), ("exp_usa.csv", "USA"),
                 ("exp_col.csv", "Colombia"), ("exp_can.csv", "Canada"),
                 ("imports_jpn.csv", "Japan_imports_world"),
                 ("bilateral_jpn.csv", "Japan_bilateral_all_partners")]
    path = joinpath(RAW, f)
    !isfile(path) && continue
    d = CSV.read(path, DataFrame)
    d.ds = first.(string.(d.date), 7)
    d = d[d.ds .>= SAMPLE_START, :]
    if lab == "Japan_bilateral_all_partners"
        # gaps for World and Australia partners separately
        for part in ("World", "Australia")
            sub = d[d.partner .== part, :]
            push!(series_info, ("Japan_$part", String.(sub.ds), f))
        end
    else
        push!(series_info, (lab, String.(unique(d.ds)), f))
    end
end

span0, span1 = string(minimum(P.ds)), string(maximum(P.ds))
o1 = DataFrame(series = String[], source = String[], n_obs = Int[],
               n_gaps = Int[], missing_months = String[])
for (lab, have, src) in series_info
    g = gaps_in(have, span0, span1)
    push!(o1, (lab, src, length(have), length(g), join(g, ";")))
end
CSV.write(joinpath(DIAG_OUT, "O1_missing_months.csv"), o1)

# O2: 2019 specifically
o2 = DataFrame(series = String[], missing_2019 = String[])
for (lab, have, src) in series_info
    g = [m for m in gaps_in(have, "2019-01", "2019-12")]
    push!(o2, (lab, isempty(g) ? "(none)" : join(g, ";")))
end
CSV.write(joinpath(DIAG_OUT, "O2_missing_2019.csv"), o2)

# Which event/horizon combos lost: if 2019-02 is an event and some lags need 2019 months
ev = CSV.read(EVENTS_PATH, DataFrame)
ev.ds = ym.(ev.date)
open(joinpath(DIAG_OUT, "O2_event_horizon_impact.txt"), "w") do io
    println(io, "2019 gaps by series (see O2_missing_2019.csv).")
    aus_gaps_2019 = [m for m in gaps_in(String.(P.ds[isfinite.(P.lq)]), "2019-01", "2019-12")]
    println(io, "Aus exports missing in 2019: ", isempty(aus_gaps_2019) ? "none" : join(aus_gaps_2019, ", "))
    println(io, "")
    println(io, "Event in 2019: 2019-02 North Queensland floods (minor, z_w=0.5).")
    println(io, "LP at horizon h uses y[t+h], y[t-1], and lags y[t-1]..y[t-6].")
    println(io, "If months around an event are missing, that event-month row drops from")
    println(io, "the regression for every h that needs a missing observation.")
    # check whether 2019-02 enters h=0 under preferred z
    z = apply_drop!(P.z_w, P.ds, DROP_2022)
    rows0 = rows_fs(P.lq, z, P.dates, P.lq, P.lp, P.rea, 0)
    enters = any(string(P.ds[t]) == "2019-02" for t in rows0)
    println(io, "Does 2019-02 enter h=0 Aus/price FS regression? ", enters ? "YES" : "NO")
    for h in 0:12
        rows = rows_fs(P.lq, z, P.dates, P.lq, P.lp, P.rea, h)
        ent = any(string(P.ds[t]) == "2019-02" for t in rows)
        println(io, "  h=$h enters: $ent")
    end
end

# O3: Comtrade collapse assertion
open(joinpath(DIAG_OUT, "O3_comtrade_unique_date.txt"), "w") do io
    println(io, "Pipeline 00_pull_comtrade.jl collapses with:")
    println(io, "  combine(groupby(df, :date), :tonnes => maximum => :tonnes)")
    println(io, "This enforces one row per date AFTER pull. Assertion check on current files:")
    any_fire = false
    for f in readdir(RAW)
        endswith(f, ".csv") || continue
        startswith(f, "exp_") || startswith(f, "aus_") || startswith(f, "imports_") || continue
        d = CSV.read(joinpath(RAW, f), DataFrame)
        !("date" in names(d)) && continue
        d.ds = first.(string.(d.date), 7)
        g = combine(groupby(d, :ds), nrow => :n)
        dups = g[g.n .> 1, :]
        if nrow(dups) > 0
            any_fire = true
            println(io, "FIRE: $f has duplicate dates: ", join(string.(dups.ds), ", "))
        else
            println(io, "OK: $f — one row per date")
        end
    end
    println(io, "")
    println(io, any_fire ? "ASSERTION WOULD FIRE on at least one file." :
                          "ASSERTION does not fire: all checked export/import files unique by date.")
end

println("wrote O1 O2 O3")
println("script: diagnostics/task_O.jl")
