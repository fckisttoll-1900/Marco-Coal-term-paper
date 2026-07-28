using CSV, DataFrames, Printf
d = CSV.read("data/raw/coal_price.csv", DataFrame)
d.ds = first.(string.(d[!,1]), 7)
v = names(d)[2]
println("Newcastle thermal around Cyclone Debbie (landfall 28 Mar 2017)")
for r in eachrow(d)
    "2017-01" <= r.ds <= "2017-08" && @printf("  %s  %8.2f\n", r.ds, r[v])
end
println("\nfor scale, the 2022 spike")
for r in eachrow(d)
    "2022-01" <= r.ds <= "2022-12" && @printf("  %s  %8.2f\n", r.ds, r[v])
end
