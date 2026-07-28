#=
validate_estimator.jl
---------------------
Validate the LP-IV / Wald-ratio estimator on simulated data with a KNOWN
demand elasticity, before trusting it on real data.

Structural model (logs, deviations from trend):

    supply:   q_t = gamma * p_t + s_t          gamma > 0
    demand:   q_t = eps   * p_t + d_t          eps   < 0   <-- target
              s_t = -delta * z_t + u_s_t       z_t = narrative disruption

Equilibrium:  p_t = (d_t - s_t) / (gamma - eps)

The Wald ratio (dq/dz) / (dp/dz) equals eps exactly, so LP-IV should recover
the true value. An uninstrumented OLS of q on p should not: demand shocks move
price too, so OLS traces a mixture of the two schedules.

Run:  julia --project=. validate_estimator.jl
=#

using Random, Statistics, LinearAlgebra, Printf

const TRUE_EPS = -0.25      # true demand elasticity to be recovered
const GAMMA    =  0.40      # supply elasticity
const T        =  300       # months (~2000-2025)
const NLAGS    =  2

"""
    simulate(seed; n_events, delta)

Return (q, q_aus, p, z): importer quantity, Australian export volume, price,
and the binary narrative shock. Events are drawn inside the cyclone season, so
the shock is seasonally clustered exactly as in the real data.
"""
function simulate(seed::Int; n_events::Int = 11, delta::Float64 = 0.12)
    rng = MersenneTwister(seed)

    # persistent demand shifter -- this is what biases OLS
    d = zeros(T)
    for t in 2:T
        d[t] = 0.85 * d[t-1] + 0.05 * randn(rng)
    end

    # non-event supply noise
    u_s = zeros(T)
    for t in 2:T
        u_s[t] = 0.50 * u_s[t-1] + 0.03 * randn(rng)
    end

    # cyclone season: Dec-Apr
    season = [((t - 1) % 12) in (0, 1, 2, 3, 11) for t in 1:T]
    cand   = findall(season)
    ev     = shuffle(rng, cand)[1:min(n_events, length(cand))]
    z      = zeros(T)
    z[ev] .= 1.0

    s     = -delta .* z .+ u_s
    p     = (d .- s) ./ (GAMMA - TRUE_EPS)
    q     = TRUE_EPS .* p .+ d          # quantity demanded (importer volume)
    q_aus = GAMMA .* p .+ s             # Australian export volume (supply side)
    return q, q_aus, p, z
end

"OLS with an intercept prepended; returns the coefficient vector."
function olsb(y::Vector{Float64}, X::Matrix{Float64})
    return hcat(ones(length(y)), X) \ y
end

"Stack lags 1..nlags of q and p. Returns (Xlags, valid_rows)."
function build_lags(q::Vector{Float64}, p::Vector{Float64}, nlags::Int)
    rows = (nlags + 1):length(q)
    Xl = zeros(length(rows), 2 * nlags)
    for (i, t) in enumerate(rows)
        for L in 1:nlags
            Xl[i, 2L - 1] = q[t - L]
            Xl[i, 2L]     = p[t - L]
        end
    end
    return Xl, collect(rows)
end

"Impact Wald ratio: coefficient on z in the q equation over that in the p equation."
function lp_iv_impact(q::Vector{Float64}, p::Vector{Float64}, z::Vector{Float64};
                      nlags::Int = NLAGS)
    Xl, idx = build_lags(q, p, nlags)
    X  = hcat(z[idx], Xl)
    bq = olsb(q[idx], X)[2]      # [1] is the intercept
    bp = olsb(p[idx], X)[2]
    return bq / bp, bp, bq
end

"Uninstrumented regression of q on p -- the simultaneity-biased benchmark."
function naive_ols_elasticity(q::Vector{Float64}, p::Vector{Float64};
                              nlags::Int = NLAGS)
    Xl, idx = build_lags(q, p, nlags)
    X = hcat(p[idx], Xl)
    return olsb(q[idx], X)[2]
end

"Monte Carlo over `reps` draws; returns the vector of LP-IV estimates."
function monte_carlo(; n_events::Int = 11, delta::Float64 = 0.12, reps::Int = 400)
    iv    = Float64[]
    naive = Float64[]
    dpdz  = Float64[]
    for s in 1:reps
        q, q_aus, p, z = simulate(1000 + s; n_events = n_events, delta = delta)
        e, bp, _ = lp_iv_impact(q, p, z)
        isfinite(e) || continue
        push!(iv, e)
        push!(naive, naive_ols_elasticity(q, p))
        push!(dpdz, bp)
    end
    return iv, naive, dpdz
end

# ---- main -------------------------------------------------------------
function main()
    iv, naive, dpdz = monte_carlo()

    println("="^62)
    @printf("True demand elasticity        : %+.3f\n", TRUE_EPS)
    @printf("LP-IV (Wald ratio)  mean      : %+.3f   sd %.3f\n",
            mean(iv), std(iv))
    @printf("Naive OLS           mean      : %+.3f   sd %.3f\n",
            mean(naive), std(naive))
    @printf("OLS bias                      : %+.3f\n", mean(naive) - TRUE_EPS)
    @printf("Mean price response to shock  : %+.4f log points\n", mean(dpdz))
    @printf("IV 5th-95th percentile        : [%+.3f, %+.3f]\n",
            quantile(iv, 0.05), quantile(iv, 0.95))
    println("="^62)

    println("\nPrecision by number of events (delta = 0.12):")
    @printf("%8s %9s %8s %9s %9s\n", "events", "mean", "sd", "5%", "95%")
    for n in (6, 11, 20, 40)
        e, _, _ = monte_carlo(n_events = n, delta = 0.12, reps = 300)
        @printf("%8d %9.3f %8.3f %9.3f %9.3f\n",
                n, mean(e), std(e), quantile(e, 0.05), quantile(e, 0.95))
    end

    println("\nPrecision by disruption size (11 events):")
    @printf("%8s %9s %8s %9s %9s\n", "delta", "mean", "sd", "5%", "95%")
    for dl in (0.05, 0.12, 0.25)
        e, _, _ = monte_carlo(n_events = 11, delta = dl, reps = 300)
        @printf("%8.2f %9.3f %8.3f %9.3f %9.3f\n",
                dl, mean(e), std(e), quantile(e, 0.05), quantile(e, 0.95))
    end

    println("""

    Reading the output
    ------------------
    If LP-IV recovers TRUE_EPS while naive OLS is attenuated or wrong-signed,
    the estimator logic is sound. In the Python run of this design the naive
    OLS came out POSITIVE -- with demand shocks dominating, an uninstrumented
    regression traces the supply curve. That is the likely explanation for the
    near-zero coal elasticity in Liu (2004), and it is the motivating result
    for the paper.

    Note also that disruption SIZE buys more precision than event COUNT:
    documenting tonnage for the few large events is worth more than adding
    marginal ones.
    """)
end

main()
