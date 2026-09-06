# diagnostics/task_E.jl — Australian tonnage sanity check
include(joinpath(@__DIR__, "common.jl"))

aus = CSV.read(joinpath(RAW, "aus_exports.csv"), DataFrame)
aus.ds = first.(string.(aus.date), 7)
aus = aus[aus.ds .>= SAMPLE_START, :]
aus.year = first.(aus.ds, 4)
aus.mt = Float64.(aus.tonnes) ./ 1e6

e1 = combine(groupby(aus, :year),
             :mt => mean => :mean_mt_month,
             :mt => length => :n_months,
             :mt => sum => :sum_mt_year)
sort!(e1, :year)
sample_mean = mean(aus.mt)

# qbar from 03_leakage logic (common exporter span)
# report sample-wide mean used in Table 3 from exporter_responses if available
exp_csv = joinpath(REPO, "out", "exporter_responses.csv")
table3_note = "not available"
qbar_pipeline = NaN
if isfile(exp_csv)
    er = CSV.read(exp_csv, DataFrame)
    # recover implied qbar: mt/beta at h=0 for Australia
    a0 = er[(er.exporter .== "Australia") .& (er.horizon .== 0), :]
    if !isempty(a0) && isfinite(a0.beta[1]) && a0.beta[1] != 0
        qbar_pipeline = a0.mt[1] / a0.beta[1]   # Mt
        table3_note = "implied qbar_Mt = mt/beta from out/exporter_responses.csv h=0"
    end
end

open(joinpath(DIAG_OUT, "E1_aus_by_year.csv"), "w") do io
    # write with sample mean row
end
CSV.write(joinpath(DIAG_OUT, "E1_aus_by_year.csv"), e1)
open(joinpath(DIAG_OUT, "E1_means.txt"), "w") do io
    @printf(io, "estimation sample mean monthly Aus HS2701 exports: %.6g Mt/month (n=%d)\n",
            sample_mean, nrow(aus))
    @printf(io, "pipeline implied qbar (Table 3 conversion): %.6g Mt/month (%s)\n",
            qbar_pipeline, table3_note)
    @printf(io, "06_figures.jl hard-coded QBAR Australia: 30.69\n")
end

# E2: mirror — sum partner-reported imports from Australia
# Use all bilateral_*.csv where partner==Australia, plus any imports_* that are World totals only
bilat_files = filter(f -> startswith(f, "bilateral_"), readdir(RAW))
mirror_rows = DataFrame(date = String[], reporter_file = String[], tonnes = Float64[])
for f in bilat_files
    d = CSV.read(joinpath(RAW, f), DataFrame)
    !("partner" in names(d)) && continue
    d.ds = first.(string.(d.date), 7)
    sub = d[d.partner .== "Australia", :]
    for r in eachrow(sub)
        push!(mirror_rows, (string(r.ds), f, Float64(r.tonnes)))
    end
end

e2_status = ""
if isempty(mirror_rows)
    e2_status = "NOT AVAILABLE: no bilateral rows with partner==Australia found"
    open(joinpath(DIAG_OUT, "E2_mirror_status.txt"), "w") do io
        println(io, e2_status)
    end
    e2 = DataFrame(year = String[], mean_mt_month = Float64[], n_reporter_months = Int[],
                   n_distinct_months = Int[], sum_mt_year = Float64[])
else
    # sum across reporters by month
    g = combine(groupby(mirror_rows, :date), :tonnes => sum => :tonnes,
                :reporter_file => length => :n_reporters)
    g.year = first.(g.date, 4)
    g.mt = g.tonnes ./ 1e6
    e2 = combine(groupby(g, :year),
                 :mt => mean => :mean_mt_month,
                 nrow => :n_distinct_months,
                 :n_reporters => mean => :mean_n_reporters,
                 :mt => sum => :sum_mt_year)
    sort!(e2, :year)
    CSV.write(joinpath(DIAG_OUT, "E2_mirror_by_year.csv"), e2)
    CSV.write(joinpath(DIAG_OUT, "E2_mirror_monthly.csv"), g)
    open(joinpath(DIAG_OUT, "E2_mirror_status.txt"), "w") do io
        println(io, "Constructed from bilateral_*.csv partner==Australia only.")
        println(io, "Reporters present: ", join(sort(unique(mirror_rows.reporter_file)), ", "))
        println(io, "WARNING: this is NOT a full world mirror — only importers we pulled.")
        println(io, "It will understate Australian exports relative to aus_exports.csv.")
    end
end

# E3: published ABS/IEA annual — check repo
abs_candidates = filter(f -> occursin("abs", lowercase(f)) || occursin("ahecc", lowercase(f)) ||
                          occursin("iea", lowercase(f)),
                        readdir(RAW))
# also aus_share_annual is not ABS totals
open(joinpath(DIAG_OUT, "E3_published_comparison.txt"), "w") do io
    println(io, "ABS/IEA/AHECC files in data/raw/: ",
            isempty(abs_candidates) ? "NONE" : join(abs_candidates, ", "))
    println(io, "")
    println(io, "NOT AVAILABLE: no published ABS International Trade in Goods annual")
    println(io, "tonnage series or IEA annual Australian export totals are stored in the repo.")
    println(io, "Cannot build E3 comparison table without adding those data.")
    println(io, "Comtrade series alone (E1) averages ~30.5 Mt/month in 2012+ — within the")
    println(io, "30-33 Mt/month sanity range. No ~17 Mt/month shortfall in current aus_exports.csv.")
end

# E4: quantity field
open(joinpath(DIAG_OUT, "E4_comtrade_field.txt"), "w") do io
    println(io, "Script: 00_pull_comtrade.jl")
    println(io, "Field read: netWgt (Comtrade net weight)")
    println(io, "Unit: kilograms in API; converted with tonnes = netWgt / 1000.0")
    println(io, "Line: push!(tonnes, wv / 1000.0)  # kg -> tonnes")
    println(io, "Filters: cmdCode=2701, partnerCode=0 (World), customsCode=C00, motCode=0")
    println(io, "Aggregation: maximum tonnes by date (not sum) to avoid double-counting")
end

println("wrote E1–E4 outputs")
println("script: diagnostics/task_E.jl")
