# diagnostics/task_D.jl — Early-2022 flood month handling
include(joinpath(@__DIR__, "common.jl"))

P = load_panel()

# Spec 4 gas price check
gas_candidates = [
    "gas_ttf.csv", "ttf.csv", "henry_hub.csv", "gas_hh.csv",
    "natgas.csv", "natural_gas.csv"
]
found_gas = String[]
for f in readdir(RAW)
    fl = lowercase(f)
    if occursin("ttf", fl) || occursin("henry", fl) || occursin("natgas", fl) ||
       occursin("natural_gas", fl) || occursin("gas", fl)
        push!(found_gas, f)
    end
end
# also check panel columns
panel_gas = [n for n in names(P.panel) if occursin("gas", lowercase(n)) ||
             occursin("ttf", lowercase(n)) || occursin("henry", lowercase(n))]

open(joinpath(DIAG_OUT, "D_gas_availability.txt"), "w") do io
    println(io, "Searched data/raw/ for TTF / Henry Hub / natural gas CSVs.")
    println(io, "Matching filenames: ", isempty(found_gas) ? "NONE" : join(found_gas, ", "))
    println(io, "Matching panel columns: ", isempty(panel_gas) ? "NONE" : join(panel_gas, ", "))
    println(io, "")
    println(io, "Spec 4 NOT ESTIMATED: no monthly gas price series in the repo.")
    println(io, "To run spec 4, add e.g. data/raw/ttf.csv with columns date, usd_per_mwh")
    println(io, "(or Henry Hub USD/mmbtu) covering the estimation sample, then rebuild.")
end

# Specs 1-3
# 1: z_w with 2022 zeroed (baseline)
# 2: z_w with 2022 included
# 3: z_w with 2022 included + post-2022-M2 step dummy

function lp_fs_extradummy(y, z, dates, q, p, rea, h, dummy)
    rows = rows_fs(y, z, dates, q, p, rea, h)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) && return (NaN, length(rows))
    X0 = design_fs(rows, z, q, p, rea, dates)
    X = hcat(X0, dummy[rows])
    yy = [y[t + h] - y[t - 1] for t in rows]
    try
        b = olsb(yy, X)
        return (b[2], length(rows))  # shock is still first column of X
    catch
        return (NaN, length(rows))
    end
end

post = Float64[ym(d) >= "2022-02" ? 1.0 : 0.0 for d in P.dates]

z1 = apply_drop!(P.z_w, P.ds, DROP_2022)
z2 = copy(P.z_w)
z3 = copy(P.z_w)

out = DataFrame(outcome = String[], h = Int[],
                spec1_baseline_drop2022 = Float64[],
                spec2_include2022 = Float64[],
                spec3_include2022_step = Float64[],
                spec4_gas = Float64[],
                note = String[])

for (oname, y) in (("q_aus", P.lq), ("price", P.lp))
    for h in (0, 3)
        b1, _, _, _ = lp_fs(y, z1, P.dates, P.lq, P.lp, P.rea, h)
        b2, _, _, _ = lp_fs(y, z2, P.dates, P.lq, P.lp, P.rea, h)
        b3, _ = lp_fs_extradummy(y, z3, P.dates, P.lq, P.lp, P.rea, h, post)
        push!(out, (oname, h, b1, b2, b3, NaN,
                    "spec4 NaN: gas series not in repo (see D_gas_availability.txt)"))
    end
end

CSV.write(joinpath(DIAG_OUT, "D_2022_specs.csv"), out)
println("wrote D_2022_specs.csv and D_gas_availability.txt")
println("script: diagnostics/task_D.jl")
