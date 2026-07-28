#=
02_estimate.jl
--------------
Estimate coal demand elasticities from exogenous Australian supply disruptions.

Two estimators, same identification:

  (1) LP-IV / Wald ratio.  Regress the h-step change in log price and in log
      quantity separately on the narrative shock z_t, controlling for lags.
      The elasticity at horizon h is the ratio of the two coefficients.
      This is indirect least squares: eps_h = beta^q_h / beta^p_h.

  (2) Proxy SVAR.  Estimate a reduced-form VAR, then use z_t as an external
      instrument for the supply shock (Mertens-Ravn 2013, Stock-Watson 2018).
      The impact elasticity is the ratio of impact responses.

Bands: wild bootstrap (Rademacher), which is appropriate with few events.

IMPORTANT on interpretation:
  Japan, Korea and Taiwan have negligible domestic coal production, so their
  IMPORT demand elasticity is also their CONSUMPTION elasticity -- these map
  directly onto COALMOD-World demand nodes.
  India (and China) produce domestically, so the estimate is an IMPORT demand
  elasticity only. It is MORE elastic than consumption and must NOT be fed
  into COALMOD as a consumption parameter -- the model already determines the
  domestic/import margin endogenously.

Run:  julia --project=. 02_estimate.jl
=#

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Printf

Random.seed!(20260727)

const H       = 24      # IRF horizon in months
const NLAGS   = 6       # VAR / control lags
const NBOOT   = 2000    # bootstrap replications
const SHOCKVAR = :z     # use :z_major for the major-events-only robustness

# ======================================================================
# helpers
# ======================================================================

"OLS coefficients with an intercept prepended to X."
function olsb(y::Vector{Float64}, X::Matrix{Float64})
    Xc = hcat(ones(length(y)), X)
    return Xc \ y
end

"Stack lags 1..p of the columns of Y. Returns (Xlags, valid_row_range)."
function lagmat(Y::Matrix{Float64}, p::Int)
    T, n = size(Y)
    X = zeros(T - p, n * p)
    for L in 1:p, j in 1:n
        X[:, (L - 1) * n + j] = Y[(p - L + 1):(T - L), j]
    end
    return X, (p + 1):T
end

"Drop rows with any missing across the requested columns; return clean matrix."
function cleanmat(df::DataFrame, cols::Vector{Symbol})
    sub = df[:, cols]
    keep = .!any(ismissing.(Matrix(sub)), dims = 2)[:]
    return Matrix{Float64}(coalesce.(sub[keep, :], 0.0)), keep
end

# ======================================================================
# (1) LP-IV
# ======================================================================
"""
    lpiv(q, p, z, H, nlags)

Reduced-form responses of `q` and `p` to shock `z` at horizons 0..H, and the
implied elasticity path. Returns (eps, bq, bp).
"""
function lpiv(q::Vector{Float64}, p::Vector{Float64}, z::Vector{Float64},
              H::Int, nlags::Int)
    T = length(q)
    Y = hcat(q, p)
    Xl, rows = lagmat(Y, nlags)
    eps = fill(NaN, H + 1); bq = fill(NaN, H + 1); bp = fill(NaN, H + 1)

    for h in 0:H
        idx = [t for t in rows if t + h <= T]
        isempty(idx) && continue
        r = [findfirst(==(t), collect(rows)) for t in idx]
        Xh = hcat(z[idx], Xl[r, :])
        yq = [q[t + h] - q[t - 1] for t in idx]
        yp = [p[t + h] - p[t - 1] for t in idx]
        cq = olsb(yq, Xh)[2]
        cp = olsb(yp, Xh)[2]
        bq[h + 1] = cq; bp[h + 1] = cp
        eps[h + 1] = cq / cp
    end
    return eps, bq, bp
end

# ======================================================================
# (2) Proxy SVAR
# ======================================================================
"""
    proxysvar(Y, z, nlags, H; norm_idx=1)

Reduced-form VAR + external instrument identification. `Y` columns are the
endogenous variables; `norm_idx` is the variable the shock is normalised on
(the supply variable, i.e. Australian exports). Returns (irf, impact).
`irf` is (H+1) x n.
"""
function proxysvar(Y::Matrix{Float64}, z::Vector{Float64}, nlags::Int, H::Int;
                   norm_idx::Int = 1)
    T, n = size(Y)
    Xl, rows = lagmat(Y, nlags)
    X = hcat(ones(length(rows)), Xl)
    Yd = Y[rows, :]
    B = X \ Yd                       # (1 + n*nlags) x n
    U = Yd - X * B                   # reduced-form residuals
    zz = z[rows]

    # impact vector: covariance of each residual with the instrument
    s = [cov(U[:, j], zz) for j in 1:n]
    b = s ./ s[norm_idx]             # normalise: unit impact on supply variable

    # companion matrix
    A = B[2:end, :]'                 # n x (n*nlags)
    C = zeros(n * nlags, n * nlags)
    C[1:n, :] = A
    if nlags > 1
        C[(n + 1):end, 1:(n * (nlags - 1))] = I(n * (nlags - 1))
    end

    irf = zeros(H + 1, n)
    bb = zeros(n * nlags); bb[1:n] = b
    v = copy(bb)
    for h in 0:H
        irf[h + 1, :] = v[1:n]
        v = C * v
    end
    return irf, b
end

# ======================================================================
# bootstrap
# ======================================================================
"Wild bootstrap (Rademacher) for the LP-IV elasticity path."
function boot_lpiv(q, p, z, H, nlags, nboot)
    T = length(q)
    out = fill(NaN, nboot, H + 1)
    for b in 1:nboot
        w = rand([-1.0, 1.0], T)
        # perturb the instrument, preserving its zero/one support pattern
        zb = z .* w
        e, _, _ = lpiv(q, p, zb, H, nlags)
        out[b, :] = e
    end
    return out
end

# ======================================================================
# main
# ======================================================================
panel = CSV.read(joinpath(@__DIR__, "data", "panel.csv"), DataFrame)
sort!(panel, :date)
# 2010 is reported in tonnes not kg, and 2011-01 is partial: unusable
const SAMPLE_START = "2012-01"
panel = panel[first.(string.(panel.date), 7) .>= SAMPLE_START, :]

# log transforms
panel.lq = log.(coalesce.(panel.q_aus, NaN))
panel.lp = log.(coalesce.(panel.p_real, NaN))
for c in [:x_jpn, :x_kor, :x_ind]
    if string(c) in names(panel)
        panel[!, Symbol("l", c)] = log.(coalesce.(panel[!, c], NaN))
    end
end

z = Float64.(coalesce.(panel[!, SHOCKVAR], 0.0))

println("="^66)
println("Coal demand elasticities from exogenous Australian supply shocks")
println("="^66)
@printf("sample     : %s to %s (%d months)\n", panel.date[1], panel.date[end], nrow(panel))
@printf("shock      : %s, %d event months\n", SHOCKVAR, Int(sum(z)))
println()

# ---- first stage: does the shock actually cut Australian exports? -----
lq = panel.lq; lp = panel.lp
ok = .!isnan.(lq) .& .!isnan.(lp)
epsq, bq, bp = lpiv(lq[ok], lp[ok], z[ok], H, NLAGS)
@printf("FIRST STAGE (impact)\n")
@printf("  Australian exports  : %+.4f log points\n", bq[1])
@printf("  coal price          : %+.4f log points\n", bp[1])
println("  (expect exports DOWN and price UP; if not, stop here)")
println()

# ---- naive OLS, for the motivation section ---------------------------
Yn = hcat(lq[ok], lp[ok])
Xn, rn = lagmat(Yn, NLAGS)
naive = olsb(lq[ok][rn], hcat(lp[ok][rn], Xn))[2]
@printf("NAIVE OLS (no instrument): %+.3f   <- simultaneity-biased benchmark\n\n", naive)

# ---- elasticity by importer ------------------------------------------
labels = Dict(:lx_kor => ("Korea",  "consumption (no domestic production)"),
              :lx_jpn => ("Japan",  "consumption (no domestic production)"),
              :lx_ind => ("India",  "IMPORT ONLY - do not feed to COALMOD"))

results = DataFrame(country = String[], horizon = Int[], elasticity = Float64[],
                    lo90 = Float64[], hi90 = Float64[], interpretation = String[])

for (col, (name, interp)) in labels
    string(col) in names(panel) || continue
    lx = panel[!, col]
    ok2 = .!isnan.(lx) .& .!isnan.(lp)
    sum(ok2) < 60 && (@warn "too few observations, skipping" country = name; continue)

    e, _, _ = lpiv(Float64.(lx[ok2]), Float64.(lp[ok2]), z[ok2], H, NLAGS)
    bs = boot_lpiv(Float64.(lx[ok2]), Float64.(lp[ok2]), z[ok2], H, NLAGS, NBOOT)

    println("-"^66)
    @printf("%s  [%s]\n", uppercase(name), interp)
    for h in [0, 3, 6, 12]
        h > H && continue
        col_b = filter(!isnan, bs[:, h + 1])
        lo = isempty(col_b) ? NaN : quantile(col_b, 0.05)
        hi = isempty(col_b) ? NaN : quantile(col_b, 0.95)
        @printf("  h=%2d   eps = %+.3f   90%% band [%+.3f, %+.3f]\n", h, e[h+1], lo, hi)
        push!(results, (name, h, e[h+1], lo, hi, interp))
    end
    println()
end

# ---- proxy SVAR cross-check ------------------------------------------
if "lx_kor" in names(panel)
    lx = panel.lx_kor
    ok3 = .!isnan.(lq) .& .!isnan.(lp) .& .!isnan.(lx)
    Y = hcat(Float64.(lq[ok3]), Float64.(lp[ok3]), Float64.(lx[ok3]))
    irf, b = proxysvar(Y, z[ok3], NLAGS, H; norm_idx = 1)
    println("-"^66)
    println("PROXY SVAR cross-check (Korea)")
    @printf("  impact elasticity = %+.3f  (LP-IV should be similar)\n", b[3] / b[2])
    CSV.write(joinpath(@__DIR__, "out", "irf_svar.csv"),
              DataFrame(h = 0:H, q_aus = irf[:, 1], price = irf[:, 2], x_kor = irf[:, 3]))
    println("  wrote out/irf_svar.csv")
end

CSV.write(joinpath(@__DIR__, "out", "elasticities.csv"), results)
println()
println("wrote out/elasticities.csv")
