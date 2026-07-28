using CSV, DataFrames, Statistics, Printf
p = CSV.read("data/panel.csv", DataFrame)
d = dropmissing(p, [:q_aus])
ds = string.(d.date)
base = mean(d.q_aus) / 1e6
@printf("q_aus: %d months, %s to %s, mean %.1f Mt\n", nrow(d), ds[1], ds[end], base)
for (lab, win) in [
    ("2010-11 Queensland floods", ["2010-10","2010-11","2010-12","2011-01","2011-02","2011-03","2011-04"]),
    ("Cyclone Debbie 2017",       ["2017-01","2017-02","2017-03","2017-04","2017-05","2017-06"])]
    println("\n", lab)
    for w in win
        i = findfirst(x -> startswith(x, w), ds)
        i === nothing && continue
        v = d.q_aus[i] / 1e6
        @printf("  %s  %5.1f Mt  %+6.1f%%\n", w, v, 100*(v - base)/base)
    end
end
