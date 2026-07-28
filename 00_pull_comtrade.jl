#=
00_pull_comtrade.jl
-------------------
Pull monthly coal (HS 2701) trade from UN Comtrade into data/raw/.

TWO MODES
---------
No key  : public preview endpoint. Works without registration, but Comtrade
          allows only ONE period per call, so a 2000-2025 pull is ~312 calls
          per country (~10-15 min for all four). No setup, just slow.

With key: authenticated endpoint, 12 periods (one year) per call, so ~26 calls
          per country (~2 min). Free registration:
              https://uncomtrade.org/docs/api-subscription-keys/
          Then, before running:
              export COMTRADE_KEY="your-key"        # shell
          or inside Julia:
              ENV["COMTRADE_KEY"] = "your-key"

The script chunks adaptively: if the server rejects a batch as too many
periods, it halves the batch and retries, down to one period per call. So it
works whatever the current limit happens to be.

Series written:
    Australia  exports  -> aus_exports.csv
    Japan      imports  -> imports_jpn.csv
    Korea Rep. imports  -> imports_kor.csv
    India      imports  -> imports_ind.csv

Setup:
    julia --project=. -e 'using Pkg; Pkg.add(["HTTP","JSON3","DataFrames","CSV"])'

Run (from inside coalEl/):
    julia --project=. 00_pull_comtrade.jl              # 2000-2025
    julia --project=. 00_pull_comtrade.jl 1995 2025    # custom span
=#

using HTTP, JSON3, DataFrames, CSV, Dates, Printf

const BASE_PREVIEW = "https://comtradeapi.un.org/public/v1/preview/C/M/HS"
const BASE_DATA    = "https://comtradeapi.un.org/data/v1/get/C/M/HS"
const CMD          = "2701"
const RAW          = joinpath(@__DIR__, "data", "raw")
const KEY          = get(ENV, "COMTRADE_KEY", "")

# (filename, M49 reporterCode, flowCode, label)
# Australia is the disrupted exporter; the rest are the competitors whose
# response measures leakage. USA, Canada and Mongolia matter specifically
# because they are the metallurgical-coal substitutes for the Bowen Basin.
const TARGETS = [
    ("aus_exports.csv", "36",  "X", "Australia exports"),
    ("imports_jpn.csv", "392", "M", "Japan imports"),
    ("imports_kor.csv", "410", "M", "Korea imports"),
    ("imports_ind.csv", "699", "M", "India imports"),
    ("exp_idn.csv",     "360", "X", "Indonesia exports"),
    ("exp_rus.csv",     "643", "X", "Russia exports"),
    ("exp_zaf.csv",     "710", "X", "South Africa exports"),
    ("exp_col.csv",     "170", "X", "Colombia exports"),
    ("exp_usa.csv",     "842", "X", "USA exports"),
    ("exp_can.csv",     "124", "X", "Canada exports"),
    ("exp_mng.csv",     "496", "X", "Mongolia exports"),
]

# preview allows 1 period per call; authenticated allows 12
const INIT_CHUNK = isempty(KEY) ? 1 : 12
const SLEEP_OK   = isempty(KEY) ? 0.35 : 1.0
const SLEEP_RETRY = 5.0

"Raw GET. Returns (ok::Bool, records_or_message)."
function try_fetch(reporter::String, flow::String, periods::Vector{String})
    # customsCode/motCode/partner2Code MUST be pinned to their TOTAL values.
    # Without them the API returns the total PLUS every breakdown, which
    # double-counts tonnage (symptom: ">12 months" in a year, inflated Mt).
    query = [
        "reporterCode" => reporter,
        "partnerCode"  => "0",          # World
        "partner2Code" => "0",          # World
        "period"       => join(periods, ","),
        "flowCode"     => flow,
        "cmdCode"      => CMD,
        "customsCode"  => "C00",        # TOTAL customs procedure
        "motCode"      => "0",          # TOTAL mode of transport
    ]
    url     = isempty(KEY) ? BASE_PREVIEW : BASE_DATA
    headers = isempty(KEY) ? Pair{String,String}[] :
                             ["Ocp-Apim-Subscription-Key" => KEY]

    r = try
        HTTP.get(url; query = query, headers = headers,
                 status_exception = false, retries = 0)
    catch err
        return (false, "connection error: $err")
    end

    if r.status == 200
        j = JSON3.read(r.body)
        recs = (haskey(j, :data) && j.data !== nothing) ? collect(j.data) : Any[]
        return (true, recs)
    end
    body = String(r.body)
    return (false, "HTTP $(r.status): $(first(body, min(length(body), 160)))")
end

"Is this the server telling us the batch is too large?"
too_many_periods(msg::AbstractString) =
    occursin("Maximum number of periods", msg) || occursin("too many", lowercase(msg))

rate_limited(msg::AbstractString) =
    occursin("429", msg) || occursin("quota", lowercase(msg)) ||
    occursin("rate limit", lowercase(msg))

"""
    fetch_adaptive(reporter, flow, periods, chunk) -> (records, chunk)

Request `periods` in batches of `chunk`. If the server rejects a batch as too
many periods, halve `chunk` and retry. Returns everything collected plus the
chunk size that actually worked, so the caller can reuse it.
"""
function fetch_adaptive(reporter::String, flow::String,
                        periods::Vector{String}, chunk::Int)
    out = Any[]
    i = 1
    while i <= length(periods)
        batch = periods[i:min(i + chunk - 1, end)]
        ok, res = try_fetch(reporter, flow, batch)

        if ok
            append!(out, res)
            i += length(batch)
            sleep(SLEEP_OK)
        elseif too_many_periods(res) && chunk > 1
            chunk = max(1, chunk ÷ 2)
            @info "reducing batch size" chunk
            # do not advance i -- retry the same periods with a smaller chunk
        elseif rate_limited(res)
            @warn "rate limited, backing off" seconds = SLEEP_RETRY
            sleep(SLEEP_RETRY)
        else
            # genuine failure for this batch: skip it and keep going
            @warn "batch failed" periods = (batch[1], batch[end]) msg = res
            i += length(batch)
            sleep(SLEEP_OK)
        end
    end
    return out, chunk
end

"Pull one series across all years; returns a tidy (date, tonnes) frame."
function pull_series(reporter::String, flow::String, y0::Int, y1::Int)
    dates  = String[]
    tonnes = Float64[]
    chunk  = INIT_CHUNK

    for y in y0:y1
        periods = [@sprintf("%d%02d", y, m) for m in 1:12]
        recs, chunk = fetch_adaptive(reporter, flow, periods, chunk)

        n = 0
        for rec in recs
            # belt and braces: drop any breakdown rows that slip through
            cc = string(get(rec, :customsCode, "C00"))
            mc = string(get(rec, :motCode, "0"))
            (cc == "C00" || cc == "") || continue
            (mc == "0"   || mc == "") || continue

            w = get(rec, :netWgt, nothing)
            (w === nothing || w === missing) && continue
            wv = try Float64(w) catch; continue end
            wv > 0 || continue
            p = string(get(rec, :period, ""))
            length(p) == 6 || continue
            push!(dates, p[1:4] * "-" * p[5:6])
            push!(tonnes, wv / 1000.0)          # kg -> tonnes
            n += 1
        end
        @printf("  %d: %2d months\n", y, n)
        flush(stdout)
    end

    isempty(dates) && return DataFrame(date = String[], tonnes = Float64[])
    df = DataFrame(date = dates, tonnes = tonnes)
    # maximum, NOT sum: there should be exactly one TOTAL row per month, and
    # if a breakdown row survives the filter the total is the larger of the two
    df = combine(groupby(df, :date), :tonnes => maximum => :tonnes)
    sort!(df, :date)
    return df
end

"Months missing inside the observed span."
function gap_report(df::DataFrame)
    isempty(df) && return String[]
    d0 = Date(df.date[1] * "-01")
    d1 = Date(df.date[end] * "-01")
    want = [Dates.format(d, "yyyy-mm") for d in d0:Month(1):d1]
    have = Set(df.date)
    return [w for w in want if !(w in have)]
end

function main()
    y0 = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
    y1 = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2025
    mkpath(RAW)

    if isempty(KEY)
        println("using PUBLIC PREVIEW endpoint (no key)")
        println("  Comtrade allows 1 period per call here, so this is slow:")
        @printf("  roughly %d calls per country. A free key makes it ~12x faster:\n",
                (y1 - y0 + 1) * 12)
        println("  https://uncomtrade.org/docs/api-subscription-keys/")
    else
        println("using AUTHENTICATED endpoint (key found)")
    end
    @printf("sample: %d-%d\n", y0, y1)

    for (fname, reporter, flow, label) in TARGETS
        println("\n", label)
        df = pull_series(reporter, flow, y0, y1)

        if isempty(df)
            println("  -> nothing retrieved, skipping $fname")
            continue
        end

        CSV.write(joinpath(RAW, fname), df)
        gaps = gap_report(df)
        @printf("  -> %s: %d months, %s to %s, mean %.1f Mt/month\n",
                fname, nrow(df), df.date[1], df.date[end],
                sum(df.tonnes) / nrow(df) / 1e6)
        if !isempty(gaps)
            shown = join(first(gaps, min(length(gaps), 8)), ", ")
            @printf("     GAPS (%d): %s%s\n", length(gaps), shown,
                    length(gaps) > 8 ? " ..." : "")
        end
    end

    println("\nStill needed: coal_price.csv, rea.csv (both done via 00b), " *
            "uscpi.csv (FRED CPIAUCSL).")
    println("Then: julia --project=. 01_build_panel.jl")
end

main()
