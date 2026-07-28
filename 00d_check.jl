using CSV, DataFrames, Printf
d = CSV.read("data/raw/aus_exports.csv", DataFrame)
d.ds = first.(string.(d.date), 7)          # "2010-09"
sort!(d, :ds)

println("raw tonnes, 2010-09 to 2011-06")
for r in eachrow(d)
    "2010-09" <= r.ds <= "2011-06" && @printf("  %s  %15.1f t\n", r.ds, r.tonnes)
end

println("\nannual totals — Australia is ~300 Mt/yr in 2010, ~390 Mt/yr lately")
d.yr = first.(d.ds, 4)
for g in groupby(d, :yr)
    @printf("  %s  %6.1f Mt  (%d months)\n", g.yr[1], sum(g.tonnes)/1e6, nrow(g))
end
