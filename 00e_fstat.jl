using CSV, DataFrames, Statistics, LinearAlgebra, Printf
p = CSV.read("data/panel.csv", DataFrame)
p = p[first.(string.(p.date),7) .>= "2012-01", :]
d = dropmissing(p, [:q_aus, :p_real])
lq, lp, z = log.(d.q_aus), log.(d.p_real), Float64.(d.z)
L = 6; T = length(lq); rows = (L+1):T
X = hcat(ones(length(rows)), z[rows], [lq[t-l] for t in rows, l in 1:L], [lp[t-l] for t in rows, l in 1:L])
for (nm, y) in [("price", lp[rows]), ("q_aus", lq[rows])]
    b = X \ y; e = y - X*b; s2 = sum(e.^2)/(length(y)-size(X,2))
    se = sqrt(s2 * inv(X'X)[2,2])
    @printf("%-6s  beta=%+.4f  se=%.4f  t=%+.2f  F=%.2f\n", nm, b[2], se, b[2]/se, (b[2]/se)^2)
end
