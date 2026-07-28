#=
00b_convert_xlsx.jl
-------------------
Convert the two local spreadsheets into the CSV schema 01_build_panel.jl expects
(a `date` column as YYYY-MM plus one value column):

    CMO-Historical-Data-Monthly.xlsx  →  data/raw/coal_price.csv  (date, usd_per_t)
    igrea.xlsx                        →  data/raw/rea.csv         (date, index)

Pink Sheet: sheet "Monthly Prices", column "Coal, Australian". Header rows sit
above the data; dates are like "1960M01". Non-numeric cells (e.g. "…") are
dropped, not interpolated.

Run:  julia --project=. 00b_convert_xlsx.jl
=#

using XLSX, CSV, DataFrames, Dates, Printf

const ROOT = @__DIR__
const RAW  = joinpath(ROOT, "data", "raw")

"Convert Pink Sheet '1960M01' (or a Date) to YYYY-MM."
function to_yyyymm(x)
    if x isa Date
        return Dates.format(x, "yyyy-mm")
    end
    s = strip(string(x))
    m = match(r"^(\d{4})M(\d{2})$", s)
    m !== nothing && return m.captures[1] * "-" * m.captures[2]
    m = match(r"^(\d{4})-(\d{2})", s)
    m !== nothing && return m.captures[1] * "-" * m.captures[2]
    return nothing
end

"Months missing inside the observed span (no interpolation)."
function gap_report(dates::Vector{String})
    isempty(dates) && return String[]
    d0 = Date(dates[1] * "-01")
    d1 = Date(dates[end] * "-01")
    want = [Dates.format(d, "yyyy-mm") for d in d0:Month(1):d1]
    have = Set(dates)
    return [w for w in want if !(w in have)]
end

function print_summary(label::String, path::String, df::DataFrame)
    gaps = gap_report(String.(df.date))
    @printf("%s\n", label)
    @printf("  -> %s: %d rows, %s to %s\n",
            basename(path), nrow(df), df.date[1], df.date[end])
    if !isempty(gaps)
        shown = join(first(gaps, min(length(gaps), 8)), ", ")
        @printf("     GAPS (%d): %s%s\n", length(gaps), shown,
                length(gaps) > 8 ? " ..." : "")
    else
        println("     no gaps in monthly sequence")
    end
end

"""
Pink Sheet → coal_price.csv.

Finds the "Coal, Australian" header in "Monthly Prices", then reads numeric
cells only (drops ellipsis / blank placeholders).
"""
function convert_coal_price()
    src = joinpath(ROOT, "CMO-Historical-Data-Monthly.xlsx")
    isfile(src) || error("missing $src")

    sh = XLSX.readxlsx(src)["Monthly Prices"]
    nrows, ncols = size(sh[:])

    # locate header row / column for "Coal, Australian"
    hdr_row = nothing
    coal_col = nothing
    for r in 1:min(20, nrows), c in 1:ncols
        v = sh[r, c]
        v === missing && continue
        if occursin(r"(?i)^coal,\s*australian", strip(string(v)))
            hdr_row = r
            coal_col = c
            break
        end
    end
    (hdr_row === nothing || coal_col === nothing) &&
        error("could not find 'Coal, Australian' in Monthly Prices")

    dates = String[]
    vals  = Float64[]
    for r in (hdr_row + 1):nrows
        d = to_yyyymm(sh[r, 1])
        d === nothing && continue
        v = sh[r, coal_col]
        v isa Number || continue          # skip "…", blanks, etc.
        push!(dates, d)
        push!(vals, Float64(v))
    end

    df = DataFrame(date = dates, usd_per_t = vals)
    # collapse any accidental duplicate months
    df = combine(groupby(df, :date), :usd_per_t => last => :usd_per_t)
    sort!(df, :date)

    path = joinpath(RAW, "coal_price.csv")
    CSV.write(path, df)
    print_summary("Coal, Australian (World Bank Pink Sheet)", path, df)
    return df
end

"""Dallas Fed IGREA → rea.csv (date, index)."""
function convert_rea()
    src = joinpath(ROOT, "igrea.xlsx")
    isfile(src) || error("missing $src")

    xf = XLSX.readxlsx(src)
    # prefer the named sheet; fall back to the first worksheet
    names_ = XLSX.sheetnames(xf)
    sn = "Kilian Index" in names_ ? "Kilian Index" : names_[1]
    sh = xf[sn]
    nrows = size(sh[:], 1)

    dates = String[]
    vals  = Float64[]
    for r in 1:nrows
        d = to_yyyymm(sh[r, 1])
        d === nothing && continue
        v = sh[r, 2]
        v isa Number || continue
        push!(dates, d)
        push!(vals, Float64(v))
    end

    df = DataFrame(date = dates, index = vals)
    df = combine(groupby(df, :date), :index => last => :index)
    sort!(df, :date)

    path = joinpath(RAW, "rea.csv")
    CSV.write(path, df)
    print_summary("IGREA / Kilian index (Dallas Fed)", path, df)
    return df
end

function main()
    mkpath(RAW)
    println("converting spreadsheets → data/raw/\n")
    convert_coal_price()
    println()
    convert_rea()
    println("\nDone. Still needed for 01_build_panel.jl: Comtrade series " *
            "(00_pull_comtrade.jl) and uscpi.csv.")
end

main()
