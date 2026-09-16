# diagnostics/task_headline_rand.jl
# Recompute keep-2022 randomisation p-values for the three headline objects only.
# Spec: z_w retained, 6 own lags, REA+6, month FE, no price lags (lp_harm_h).
#
# Placebo pool = Nov–Apr ∩ valid LP rows for that (outcome, h), excluding true
# event months. Every draw places the full multiset of *effective* treated-month
# weights (those true events that enter the regression sample) onto distinct
# pool months, so all effective treated months contribute in every draw.
include(joinpath(@__DIR__, "exporters_harm.jl"))

const NBOOT = 10_000

P = load_panel()
z_keep = z_baseline(P)   # copy of z_w; 2022 retained
T = length(P.dates)
E = load_exporters_harm(P)

"""
Seasonal placebo pool restricted to valid regression rows for (y, h).

Returns:
  true_eff  — indices of true treated months that enter rows_harm(y,z,rea,h)
  true_vals — severity weights on those effective months (multiset preserved)
  pool      — Nov–Apr valid rows excluding all true treated months
  n_valid   — number of valid estimation rows
"""
function placebo_pool(y, z_ref, h)
    valid = Set(rows_harm(y, z_ref, P.rea, h))
    true_all = findall(>(0), z_ref)
    true_eff = [i for i in true_all if i in valid]
    true_vals = z_ref[true_eff]
    season = Int[i for i in 1:T if (m = month(P.dates[i]); m >= 11 || m <= 4)]
    true_set = Set(true_all)
    pool = [i for i in season if (i in valid) && !(i in true_set)]
    return true_eff, true_vals, pool, length(valid)
end

function rand_p(y, z_ref, h; nboot = NBOOT, seed = BOOT_SEED)
    true_eff, true_vals, pool, n_valid = placebo_pool(y, z_ref, h)
    n_events = length(true_eff)
    length(pool) < n_events && return (NaN, NaN, 0, n_events, length(pool), n_valid)
    β0, _, nobs = lp_harm_h(y, z_ref, P.rea, P.dates, h)
    !isfinite(β0) && return (β0, NaN, 0, n_events, length(pool), n_valid)
    rng = MersenneTwister(seed)
    B = fill(NaN, nboot)
    for b in 1:nboot
        pick = pool[randperm(rng, length(pool))[1:n_events]]
        zpl = zeros(T)
        for (k, i) in enumerate(pick)
            zpl[i] = true_vals[k]
        end
        B[b], _, _ = lp_harm_h(y, zpl, P.rea, P.dates, h)
    end
    c = filter(isfinite, B)
    p = isempty(c) ? NaN : mean(abs.(c) .>= abs(β0))
    return (β0, p, length(c), n_events, length(pool), nobs)
end

tests = [
    ("Aus_exports_log_h0", P.lq, 0),
    ("Price_log_real_h0", P.lp, 0),
    ("Wexc_mt_h0", E.W_exc, 0),
]

rows = DataFrame(test = String[], h = Int[], beta = Float64[],
                 rand_p = Float64[], n_finite = Int[], n_events = Int[],
                 pool = Int[], n_obs = Int[],
                 shock = String[], drop_2022 = Bool[])

println("Headline keep-2022 randomisation (valid-row pool; N=$NBOOT, seed=$BOOT_SEED)...")
@printf("  T=%d  sum(z_w)=%.1f  N_z_gt0=%d\n",
        T, sum(z_keep), count(>(0), z_keep))

for (name, y, h) in tests
    @printf("  %s ...\n", name)
    te, tv, pool0, nv = placebo_pool(y, z_keep, h)
    @printf("    effective treated=%d  weights=%s  valid_rows=%d  pool=%d\n",
            length(te), string(round.(tv; digits = 1)), nv, length(pool0))
    β, p, nf, ne, pool, nobs = rand_p(y, z_keep, h)
    push!(rows, (name, h, β, p, nf, ne, pool, nobs, "z_w", false))
    @printf("    beta=%.6g  rand_p=%.4f  n_events=%d  pool=%d  n_obs=%d\n",
            β, p, ne, pool, nobs)
end

out = joinpath(DIAG_OUT, "headline_rand_keep2022.csv")
CSV.write(out, rows)

open(joinpath(DIAG_OUT, "headline_rand_keep2022.txt"), "w") do io
    println(io, "Keep-2022 headline randomisation p-values (valid-row pool)")
    println(io, "Spec: z_w retained; 6 own lags; REA+6; month FE; no price lags.")
    println(io, "Placebo: Nov–Apr ∩ valid LP rows for (outcome,h), excluding true event months;")
    println(io, "preserve the multiset of effective treated-month severity weights.")
    println(io, "Exact month-of-year composition and multi-month clustering are NOT preserved.")
    @printf(io, "Draws=%d  seed=%d  T=%d\n\n", NBOOT, BOOT_SEED, T)
    for r in eachrow(rows)
        @printf(io, "%s: beta=%.6g  rand_p=%.4f  (n_events=%d, pool=%d, n_obs=%d)\n",
                r.test, r.beta, r.rand_p, r.n_events, r.pool, r.n_obs)
    end
end

println("wrote $out")
println("script: diagnostics/task_headline_rand.jl")
