#=
01_build_panel.jl
-----------------
Assemble the monthly panel from raw downloads in data/raw/ into data/panel.csv.

Expected files in data/raw/ (see README.md for links). Each should be a CSV
with a `date` column in YYYY-MM (or a parseable date) plus the value column:

    aus_exports.csv    date, tonnes          Australian coal exports (ABS/DFAT)
    coal_price.csv     date, usd_per_t       Newcastle thermal (World Bank)
    imports_jpn.csv    date, tonnes          Japan coal imports (MoF/e-Stat)
    imports_kor.csv    date, tonnes          Korea coal imports (Korea Customs)
    imports_ind.csv    date, tonnes          India coal imports (MoC/Comtrade)
    rea.csv            date, index           Kilian global real activity
    uscpi.csv          date, cpi             US CPI, for deflating (FRED CPIAUCSL)

Missing importer files are simply skipped, so you can start with Korea alone
and add Japan and India later.

Run:  julia --project=. 01_build_panel.jl
=#

using CSV, DataFrames, Dates, Statistics

const RAW = joinpath(@__DIR__, "data", "raw")
const OUT = joinpath(@__DIR__, "data")

"Parse a YYYY-MM or full date string into the first of the month."
function tomonth(x)
    s = strip(string(x))
    d = try
        length(s) == 7 ? Date(s * "-01") : Date(first(split(s, r"[ T]")))
    catch
        return missing
    end
    return Date(year(d), month(d), 1)
end

"Read one raw file, rename its value column, return a two-column frame."
function readseries(fname::String, valcol::Symbol)
    path = joinpath(RAW, fname)
    if !isfile(path)
        @warn "missing, skipping" file = fname
        return nothing
    end
    df = CSV.read(path, DataFrame)
    cols  = names(df)
    lower = lowercase.(cols)

    # Date column: sources disagree on the name. FRED exports use
    # "observation_date" or "DATE", Comtrade "period", ours "date".
    di = findfirst(n -> n == "date" || n == "observation_date" ||
                        n == "period" || n == "month" || n == "time" ||
                        occursin("date", n), lower)
    if di === nothing
        collist = join(cols, ", ")
        error("no date column in $fname -- columns are: $collist")
    end
    dcol = cols[di]

    # Value column: the first column that is not the date column.
    vi = findfirst(i -> i != di, eachindex(cols))
    if vi === nothing
        collist = join(cols, ", ")
        error("$fname has only a date column -- columns are: $collist")
    end
    vcol = cols[vi]

    # both entries must be Pairs -- mixing `date = ...` (keyword) with a Pair
    # is a MethodError in DataFrames
    out = DataFrame(:date => tomonth.(df[!, dcol]), valcol => df[!, vcol])
    dropmissing!(out)
    isempty(out) && error("$fname parsed to zero usable rows " *
                          "(check the date format in column '$dcol')")
    # collapse any duplicate months (e.g. if source is finer than monthly)
    return combine(groupby(out, :date), valcol => sum => valcol)
end

# ---- assemble ---------------------------------------------------------
pieces = Dict{Symbol,Union{Nothing,DataFrame}}(
    :q_aus  => readseries("aus_exports.csv", :q_aus),
    :p_usd  => readseries("coal_price.csv",  :p_usd),
    :x_jpn  => readseries("imports_jpn.csv", :x_jpn),
    :x_kor  => readseries("imports_kor.csv", :x_kor),
    :x_ind  => readseries("imports_ind.csv", :x_ind),
    :rea    => readseries("rea.csv",         :rea),
    :cpi    => readseries("uscpi.csv",       :cpi),
)

frames = DataFrame[df for (_, df) in pieces if df !== nothing]
isempty(frames) && error("no raw files found in $RAW -- see README.md")
panel = reduce((a, b) -> outerjoin(a, b, on = :date), frames)
sort!(panel, :date)

# ---- real price -------------------------------------------------------
if "cpi" in names(panel) && "p_usd" in names(panel)
    base = mean(skipmissing(panel.cpi))
    panel.p_real = panel.p_usd .* (base ./ panel.cpi)
else
    @warn "no US CPI supplied -- using nominal price"
    panel.p_real = panel.p_usd
end

# ---- narrative shock series ------------------------------------------
ev = CSV.read(joinpath(OUT, "events.csv"), DataFrame)
ev.date = tomonth.(ev.date)
panel.z       = [d in ev.date ? 1.0 : 0.0 for d in panel.date]
major         = Set(ev.date[ev.severity .== "major"])
panel.z_major = [d in major ? 1.0 : 0.0 for d in panel.date]

CSV.write(joinpath(OUT, "panel.csv"), panel)

println("wrote data/panel.csv")
println("  rows          : ", nrow(panel))
println("  span          : ", minimum(panel.date), " to ", maximum(panel.date))
println("  events (all)  : ", Int(sum(panel.z)))
println("  events (major): ", Int(sum(panel.z_major)))
println("  columns       : ", join(names(panel), ", "))
