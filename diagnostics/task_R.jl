# diagnostics/task_R.jl — Extended randomisation inference
include(joinpath(@__DIR__, "exporters_harm.jl"))

const NBOOT = 10_000

P = load_panel()
z_drop = z_preferred(P)
z_keep = copy(P.z_w)
T = length(P.dates)
E = load_exporters_harm(P)

function placebo_pool(z_ref)
    true_idx = findall(>(0), z_ref)
    season = Int[i for i in 1:T if (m = month(P.dates[i]); m >= 11 || m <= 4)]
    pool = [i for i in season if !(i in Set(true_idx))]
    return true_idx, z_ref[true_idx], pool
end

function rand_p(y, z_ref, h; nboot = NBOOT, seed = BOOT_SEED)
    true_idx, true_vals, pool = placebo_pool(z_ref)
    n_events = length(true_idx)
    length(pool) < n_events && return (NaN, NaN, 0, n_events, length(pool))
    β0, _, _ = lp_harm_h(y, z_ref, P.rea, P.dates, h)
    !isfinite(β0) && return (β0, NaN, 0, n_events, length(pool))
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
    return (β0, p, length(c), n_events, length(pool))
end

rows = DataFrame(test = String[], h = Int[], beta = Float64[],
                 rand_p = Float64[], n_finite = Int[], n_events = Int[], pool = Int[])

tests = [
    ("R1_price_keep2022", P.lp, z_keep, 0),
    ("R2_price_drop2022", P.lp, z_drop, 0),
    ("R2_price_drop2022", P.lp, z_drop, 3),
    ("R2_price_drop2022", P.lp, z_drop, 6),
    ("R3_aus_exports_drop2022", P.lq, z_drop, 0),
    ("R3_aus_exports_drop2022", P.lq, z_drop, 3),
    ("R3_aus_exports_drop2022", P.lq, z_drop, 6),
    ("R3_aus_exports_drop2022", P.lq, z_drop, 12),
    ("R4_world_ex_Aus_mt", E.W_exc, z_drop, 0),
    ("R4_world_ex_Aus_mt", E.W_exc, z_drop, 6),
    ("R4_world_ex_Aus_mt", E.W_exc, z_drop, 12),
]

println("R: 10000 placebos × $(length(tests)) tests...")
for (name, y, z, h) in tests
    @printf("  %s h=%d ...\n", name, h)
    β, p, nf, ne, pool = rand_p(y, z, h)
    push!(rows, (name, h, β, p, nf, ne, pool))
end
CSV.write(joinpath(DIAG_OUT, "R_randomisation.csv"), rows)
println("wrote R_randomisation.csv")
println("script: diagnostics/task_R.jl")
