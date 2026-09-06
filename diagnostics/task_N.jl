# diagnostics/task_N.jl — Pink Sheet gas controls + Task D spec 4
include(joinpath(@__DIR__, "harm_lp.jl"))
using XLSX

# Extract Natural gas, Europe and LNG, Japan from Pink Sheet (diagnostics only)
src = joinpath(REPO, "CMO-Historical-Data-Monthly.xlsx")
isfile(src) || error("missing Pink Sheet at $src")
sh = XLSX.readxlsx(src)["Monthly Prices"]
nrows, ncols = size(sh[:])

function find_col(pat)
    for r in 1:min(20, nrows), c in 1:ncols
        v = sh[r, c]
        v === missing && continue
        occursin(pat, strip(string(v))) && return (r, c)
    end
    return (nothing, nothing)
end
function to_yyyymm(x)
    if x isa Date; return Dates.format(x, "yyyy-mm"); end
    s = strip(string(x))
    m = match(r"^(\d{4})M(\d{2})$", s)
    m !== nothing && return m.captures[1] * "-" * m.captures[2]
    return nothing
end
function extract_series(pat, outname)
    hr, hc = find_col(pat)
    (hr === nothing) && error("column not found: $pat")
    dates = String[]; vals = Float64[]
    for r in (hr + 1):nrows
        d = to_yyyymm(sh[r, 1]); d === nothing && continue
        v = sh[r, hc]; v isa Number || continue
        push!(dates, d); push!(vals, Float64(v))
    end
    df = DataFrame(date = dates, usd = vals)
    df = combine(groupby(df, :date), :usd => last => :usd)
    sort!(df, :date)
    CSV.write(joinpath(DIAG_OUT, outname), df)
    @printf("extracted %s -> %s (%d rows)\n", pat, outname, nrow(df))
    return df
end

gas_eu = extract_series(r"(?i)^Natural gas,\s*Europe", "N_gas_europe_nominal.csv")
lng_jp = extract_series(r"(?i)^Liquefied natural gas,\s*Japan|^LNG,\s*Japan", "N_lng_japan_nominal.csv")

P = load_panel()
# deflate like coal: p_real = p_usd * (mean_cpi / cpi)
cpi = Float64.(coalesce.(P.panel.cpi, NaN))
base = mean(filter(isfinite, cpi))
function align_real(df)
    m = Dict(string(r.date) => Float64(r.usd) for r in eachrow(df))
    nom = [haskey(m, string(s)) ? m[string(s)] : NaN for s in P.ds]
    real = similar(nom)
    for i in eachindex(nom)
        real[i] = (isfinite(nom[i]) && isfinite(cpi[i]) && cpi[i] > 0) ?
                  nom[i] * (base / cpi[i]) : NaN
    end
    return nom, real
end
_, gas_eu_r = align_real(gas_eu)
_, lng_jp_r = align_real(lng_jp)
lgas = [isfinite(x) && x > 0 ? log(x) : NaN for x in gas_eu_r]
llng = [isfinite(x) && x > 0 ? log(x) : NaN for x in lng_jp_r]

CSV.write(joinpath(DIAG_OUT, "N_gas_on_panel.csv"),
          DataFrame(date = String.(P.ds),
                    gas_europe_real = gas_eu_r, lng_japan_real = lng_jp_r,
                    log_gas_europe = lgas, log_lng_japan = llng))

# LP with extra gas controls (level + 6 lags), z_w with 2022 RETAINED
function lp_fs_gas(y, z, dates, q, p, rea, gas_controls::Vector, h)
    # gas_controls is vector of series, each gets level+6 lags
    rows = Int[]
    nlags = NLAGS
    T = length(y)
    for t in (nlags + 1):(T - h)
        ok = isfinite(y[t + h]) && isfinite(y[t - 1]) && isfinite(z[t])
        ok = ok && all(isfinite, @view q[(t - nlags):(t - 1)]) &&
                   all(isfinite, @view p[(t - nlags):(t - 1)])
        ok = ok && isfinite(rea[t]) && all(isfinite, @view rea[(t - nlags):(t - 1)])
        for g in gas_controls
            ok = ok && isfinite(g[t]) && all(isfinite, @view g[(t - nlags):(t - 1)])
        end
        ok && push!(rows, t)
    end
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return (NaN, length(rows))
    X0 = design_fs(rows, z, q, p, rea, dates)
    extras = zeros(length(rows), length(gas_controls) * (nlags + 1))
    c = 1
    for g in gas_controls
        extras[:, c] = g[rows]; c += 1
        for l in 1:nlags
            extras[:, c] = g[rows .- l]; c += 1
        end
    end
    X = hcat(X0, extras)
    yy = [y[t + h] - y[t - 1] for t in rows]
    try
        return (olsb(yy, X)[2], length(rows))
    catch
        return (NaN, length(rows))
    end
end

z_keep = copy(P.z_w)   # 2022 RETAINED
post = Float64[ym(d) >= "2022-02" ? 1.0 : 0.0 for d in P.dates]
z_drop = apply_drop!(P.z_w, P.ds, DROP_2022)

function lp_fs_step(y, z, dates, q, p, rea, dummy, h)
    rows = rows_fs(y, z, dates, q, p, rea, h)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return NaN
    X = hcat(design_fs(rows, z, q, p, rea, dates), dummy[rows])
    yy = [y[t + h] - y[t - 1] for t in rows]
    try; return olsb(yy, X)[2]; catch; return NaN; end
end

out = DataFrame(h = Int[],
                spec1_drop2022 = Float64[],
                spec2_keep2022 = Float64[],
                spec3_keep2022_step = Float64[],
                spec4a_gas_europe = Float64[],
                spec4b_lng_japan = Float64[],
                spec4c_both_gas = Float64[],
                n_spec4a = Int[], n_spec4b = Int[], n_spec4c = Int[])

for h in (0, 3, 6, 12)
    b1, _, _, _ = lp_fs(P.lp, z_drop, P.dates, P.lq, P.lp, P.rea, h)
    b2, _, _, _ = lp_fs(P.lp, z_keep, P.dates, P.lq, P.lp, P.rea, h)
    b3 = lp_fs_step(P.lp, z_keep, P.dates, P.lq, P.lp, P.rea, post, h)
    b4a, n4a = lp_fs_gas(P.lp, z_keep, P.dates, P.lq, P.lp, P.rea, [lgas], h)
    b4b, n4b = lp_fs_gas(P.lp, z_keep, P.dates, P.lq, P.lp, P.rea, [llng], h)
    b4c, n4c = lp_fs_gas(P.lp, z_keep, P.dates, P.lq, P.lp, P.rea, [lgas, llng], h)
    push!(out, (h, b1, b2, b3, b4a, b4b, b4c, n4a, n4b, n4c))
end
CSV.write(joinpath(DIAG_OUT, "N_price_gas_controls.csv"), out)

# verdict
open(joinpath(DIAG_OUT, "N_verdict.txt"), "w") do io
    println(io, "Price h=0 under z_w with 2022 RETAINED:")
    r0 = out[out.h .== 0, :]
    @printf(io, "  spec2 (no gas):     %+.4f\n", r0.spec2_keep2022[1])
    @printf(io, "  spec4a gas Europe:  %+.4f\n", r0.spec4a_gas_europe[1])
    @printf(io, "  spec4b LNG Japan:   %+.4f\n", r0.spec4b_lng_japan[1])
    @printf(io, "  spec4c both:        %+.4f\n", r0.spec4c_both_gas[1])
    @printf(io, "  spec1 (drop 2022):  %+.4f\n", r0.spec1_drop2022[1])
    survives = r0.spec4c_both_gas[1] > 0.02   # still clearly positive
    if survives
        println(io, "VERDICT: positive price response SURVIVES gas controls.")
        println(io, "2022 exclusion is NOT defensible on gas-confounding grounds alone;")
        println(io, "paper should report the positive price response when 2022 is retained,")
        println(io, "or justify exclusion on other grounds.")
    else
        println(io, "VERDICT: positive price response does NOT survive gas controls")
        println(io, "(moves to near zero / negative). 2022 exclusion is more defensible;")
        println(io, "preferred near-zero price result can stand.")
    end
end
println("wrote N_* outputs")
println("script: diagnostics/task_N.jl")
