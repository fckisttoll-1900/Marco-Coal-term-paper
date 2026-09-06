# diagnostics/task_B.jl — Bootstrap harmonisation
include(joinpath(@__DIR__, "common.jl"))

const HLIST = [0, 1, 3, 6, 12]

P = load_panel()
# Paper-preferred instrument for price / first-stage objects
z_pref = apply_drop!(P.z_w, P.ds, DROP_2022)

# ---------- wild / DWB helpers for FS LP ----------
function boot_fs_bands(y, z, dates, q, p, rea, h; nboot = NBOOT_B, seed = BOOT_SEED)
    rows = rows_fs(y, z, dates, q, p, rea, h)
    (length(rows) < 40 || sum(abs, z[rows]) < 1.5) &&
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    X = design_fs(rows, z, q, p, rea, dates)
    yy = [y[t + h] - y[t - 1] for t in rows]
    Z = hcat(ones(length(rows)), X)
    bhat = try Z \ yy catch; return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN) end
    e = yy - Z * bhat
    point = bhat[2]
    # OLS SE
    df = length(rows) - size(Z, 2)
    s2 = sum(e.^2) / max(df, 1)
    se_ols = try sqrt(s2 * inv(Z' * Z)[2, 2]) catch; NaN end
    se_nw = nw_se_fs(y, z, dates, q, p, rea, h; bandwidth = h + 1)

    # wild (observation-level)
    rng = MersenneTwister(seed)
    wild = fill(NaN, nboot)
    for b in 1:nboot
        yb = Z * bhat + e .* rand(rng, [-1.0, 1.0], length(rows))
        try wild[b] = (Z \ yb)[2] catch; end
    end
    # dependent wild: same weight within calendar blocks of the *original* t
    rng2 = MersenneTwister(seed + 1)
    dwb = fill(NaN, nboot)
    # map row index -> block id from original calendar position
    blk = [cld(t, BLOCK) for t in rows]
    ublk = unique(blk)
    for b in 1:nboot
        wmap = Dict(u => rand(rng2, [-1.0, 1.0]) for u in ublk)
        w = [wmap[blk[i]] for i in 1:length(rows)]
        yb = Z * bhat + e .* w
        try dwb[b] = (Z \ yb)[2] catch; end
    end
    # moving-block on full series
    rng3 = MersenneTwister(seed + 2)
    T = length(y)
    mbb = fill(NaN, nboot)
    for b in 1:nboot
        ii = block_idx(T, BLOCK, rng3)
        β, _, _, _ = lp_fs(y[ii], z[ii], dates[ii], q[ii], p[ii], rea[ii], h)
        mbb[b] = β
    end
    function qband(v)
        c = filter(isfinite, v)
        isempty(c) && return (NaN, NaN)
        return (quantile(c, 0.05), quantile(c, 0.95))
    end
    wlo, whi = qband(wild)
    dlo, dhi = qband(dwb)
    mlo, mhi = qband(mbb)
    return (point, se_ols, se_nw, mlo, mhi, wlo, whi, dlo, dhi)
end

function push_row!(df, object, h, point, se_ols, se_nw, mlo, mhi, wlo, whi, dlo, dhi; note = "")
    impl = (isfinite(mhi) && isfinite(mlo)) ? (mhi - mlo) / (2 * 1.645) : NaN
    push!(df, (object, h, point, se_ols, se_nw, mlo, mhi, wlo, whi, dlo, dhi, impl, note,
               BOOT_SEED, NBOOT_B, BLOCK))
end

out = DataFrame(
    object = String[], h = Int[], point = Float64[], se_ols = Float64[],
    se_nw = Float64[],
    mbb_lo90 = Float64[], mbb_hi90 = Float64[],
    wild_lo90 = Float64[], wild_hi90 = Float64[],
    dwb_lo90 = Float64[], dwb_hi90 = Float64[],
    implied_SE_from_block = Float64[], note = String[],
    seed = Int[], nboot = Int[], block = Int[]
)

println("B: Aus exports (log) and price under z_w + drop 2022 (FS controls)...")
for h in HLIST
    for (name, y) in (("aus_exports_log_FS", P.lq), ("price_real_log_FS", P.lp))
        pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi =
            boot_fs_bands(y, z_pref, P.dates, P.lq, P.lp, P.rea, h)
        push_row!(out, name, h, pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi;
                  note = "FS: z_w drop2022 REA monthFE")
    end
end

# ---------- exporters as in 03_leakage (binary z, simple LP) ----------
println("B: loading exporter files for simple-LP objects...")
const EXPORTERS = [
    ("aus_exports.csv", "Australia", true),
    ("exp_idn.csv", "Indonesia", false),
    ("exp_usa.csv", "USA", false),
    ("exp_col.csv", "Colombia", false),
    ("exp_can.csv", "Canada", false),
    ("exp_rus.csv", "Russia", false),
    ("exp_zaf.csv", "South Africa", false),
    ("exp_mng.csv", "Mongolia", false),
]

function readexp(fname)
    path = joinpath(RAW, fname)
    !isfile(path) && return nothing
    d = CSV.read(path, DataFrame)
    d.ds = first.(string.(d.date), 7)
    d = d[d.ds .>= SAMPLE_START, :]
    isempty(d) && return nothing
    return d
end

function calendar(a, b)
    [Dates.format(d, "yyyy-mm") for d in Date(a * "-01"):Month(1):Date(b * "-01")]
end

series = Tuple{String,Bool,DataFrame}[]
for (f, lab, isa_) in EXPORTERS
    d = readexp(f)
    d === nothing && continue
    # coverage gates matching 03_leakage
    dens = nrow(d) / length(calendar(d.ds[1], d.ds[end]))
    keep = isa_ || (nrow(d) >= 110 && d.ds[end] >= "2024-01" && dens >= 0.85)
    keep && push!(series, (lab, isa_, d))
end

start = maximum(d.ds[1] for (_, _, d) in series)
stop  = minimum(d.ds[end] for (_, _, d) in series)
dates_e = calendar(max(start, SAMPLE_START), stop)
J = length(series)
Q = fill(NaN, length(dates_e), J)
labels = String[]; isaus = Bool[]
for (j, (lab, isa_, d)) in enumerate(series)
    m = Dict(d.ds .=> Float64.(d.tonnes))
    for (i, s) in enumerate(dates_e)
        haskey(m, s) && (Q[i, j] = m[s])
    end
    push!(labels, lab); push!(isaus, isa_)
end
evset = Set(string.(ym.(CSV.read(EVENTS_PATH, DataFrame).date)))
z_bin_e = Float64[s in evset ? 1.0 : 0.0 for s in dates_e]
dates_ed = [Date(s * "-01") for s in dates_e]
qbar = [mean(filter(isfinite, Q[:, j])) for j in 1:J]
aidx = findfirst(isaus)

function boot_simple_bands(y, z, dates, h; nboot = NBOOT_B, seed = BOOT_SEED)
    rows = rows_simple(y, z, h)
    (length(rows) < 40 || sum(z[rows]) < 2) &&
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    D = monthdummies(dates)
    L = zeros(length(rows), NLAGS)
    for (i, t) in enumerate(rows), l in 1:NLAGS
        L[i, l] = y[t - l]
    end
    X = hcat(z[rows], L, D[rows, :])
    yy = [y[t + h] - y[t - 1] for t in rows]
    Z = hcat(ones(length(rows)), X)
    bhat = try Z \ yy catch; return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN) end
    e = yy - Z * bhat
    point = bhat[2]
    df = length(rows) - size(Z, 2)
    s2 = sum(e.^2) / max(df, 1)
    se_ols = try sqrt(s2 * inv(Z' * Z)[2, 2]) catch; NaN end
    # NW on simple design
    se_nw = begin
        try
            n = length(rows)
            Qinv = inv(Z' * Z / n)
            scores = [Z[i, :] * e[i] for i in 1:n]
            S = zeros(size(Z, 2), size(Z, 2))
            for i in 1:n; S .+= scores[i] * scores[i]'; end
            Lband = h + 1
            for ell in 1:Lband
                w = 1.0 - ell / (Lband + 1)
                Gam = zeros(size(Z, 2), size(Z, 2))
                for i in (ell + 1):n
                    Gam .+= scores[i] * scores[i - ell]'
                end
                S .+= w * (Gam + Gam')
            end
            S ./= n
            sqrt((Qinv * S * Qinv / n)[2, 2])
        catch
            NaN
        end
    end
    rng = MersenneTwister(seed)
    wild = fill(NaN, nboot)
    for b in 1:nboot
        yb = Z * bhat + e .* rand(rng, [-1.0, 1.0], length(rows))
        try wild[b] = (Z \ yb)[2] catch; end
    end
    rng2 = MersenneTwister(seed + 1)
    dwb = fill(NaN, nboot)
    blk = [cld(t, BLOCK) for t in rows]
    ublk = unique(blk)
    for b in 1:nboot
        wmap = Dict(u => rand(rng2, [-1.0, 1.0]) for u in ublk)
        w = [wmap[blk[i]] for i in 1:length(rows)]
        try dwb[b] = (Z \ (Z * bhat + e .* w))[2] catch; end
    end
    rng3 = MersenneTwister(seed + 2)
    T = length(y)
    mbb = fill(NaN, nboot)
    for b in 1:nboot
        ii = block_idx(T, BLOCK, rng3)
        β, _, _, _ = lp_simple(y[ii], z[ii], dates[ii], h)
        mbb[b] = β
    end
    qband(v) = (c = filter(isfinite, v); isempty(c) ? (NaN, NaN) : (quantile(c, 0.05), quantile(c, 0.95)))
    wlo, whi = qband(wild); dlo, dhi = qband(dwb); mlo, mhi = qband(mbb)
    return (point, se_ols, se_nw, mlo, mhi, wlo, whi, dlo, dhi)
end

println("B: exporter log & Mt bands (03_leakage-style binary z)...")
for j in 1:J
    ylog = log.(Q[:, j])
    for h in HLIST
        pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi =
            boot_simple_bands(ylog, z_bin_e, dates_ed, h)
        push_row!(out, "exporter_log_$(labels[j])", h, pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi;
                  note = "03-style binary z; no REA")
        # Mt = beta * qbar / 1e6 — bands scaled the same way as pipeline (fixed qbar)
        sc = qbar[j] / 1e6
        push_row!(out, "exporter_Mt_$(labels[j])", h,
                  pt * sc, se * sc, nw * sc,
                  mlo * sc, mhi * sc, wlo * sc, whi * sc, dlo * sc, dhi * sc;
                  note = "log beta x fixed qbar/1e6")
    end
end

# leakage ratio with MBB recomputed per draw; wild/DWB on ratio via residual boots of system
println("B: leakage ratio bands...")
function leakage_at(Bpath, qbar, aidx, h)
    den = -Bpath[h + 1, aidx] * qbar[aidx]
    num = 0.0
    for j in 1:length(qbar)
        j == aidx && continue
        isnan(Bpath[h + 1, j]) && continue
        num += Bpath[h + 1, j] * qbar[j]
    end
    (isfinite(den) && den > 0.2e6) ? num / den : NaN
end

function betas_matrix(Q, z, dates, Hmax)
    J = size(Q, 2)
    B = fill(NaN, Hmax + 1, J)
    for j in 1:J
        for h in 0:Hmax
            β, _, _, _ = lp_simple(log.(Q[:, j]), z, dates, h)
            B[h + 1, j] = β
        end
    end
    return B
end

Hmax = maximum(HLIST)
B0 = betas_matrix(Q, z_bin_e, dates_ed, Hmax)
for h in HLIST
    pt = leakage_at(B0, qbar, aidx, h)
    # MBB of ratio (as in 03_leakage)
    rng = MersenneTwister(BOOT_SEED)
    T = length(dates_e)
    mbb = fill(NaN, NBOOT_B)
    for b in 1:NBOOT_B
        ii = block_idx(T, BLOCK, rng)
        Bb = betas_matrix(Q[ii, :], z_bin_e[ii], dates_ed[ii], Hmax)
        mbb[b] = leakage_at(Bb, qbar, aidx, h)
    end
    c = filter(isfinite, mbb)
    mlo = isempty(c) ? NaN : quantile(c, 0.05)
    mhi = isempty(c) ? NaN : quantile(c, 0.95)
    # Wild / DWB for leakage: NOT available as a single residual regression —
    # report explicitly
    push_row!(out, "leakage_ratio", h, pt, NaN, NaN, mlo, mhi, NaN, NaN, NaN, NaN;
              note = "MBB recomputes ratio per draw (as 03_leakage). Wild/DWB/OLS SE not available for ratio without a single-equation residualisation — left as NaN.")
end

# Japan series from bilateral
println("B: Japan total and from-Australia...")
jp = CSV.read(joinpath(RAW, "bilateral_jpn.csv"), DataFrame)
jp.ds = first.(string.(jp.date), 7)
jp = jp[jp.ds .>= SAMPLE_START, :]
cal = calendar(minimum(jp.ds), maximum(jp.ds))
z_j = Float64[s in evset ? 1.0 : 0.0 for s in cal]
dates_j = [Date(s * "-01") for s in cal]
for (src, oname) in (("World", "japan_total_log"), ("Australia", "japan_from_aus_log"))
    sub = jp[jp.partner .== src, :]
    m = Dict(sub.ds .=> Float64.(sub.tonnes))
    y = [haskey(m, s) ? log(m[s]) : NaN for s in cal]
    for h in HLIST
        pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi =
            boot_simple_bands(y, z_j, dates_j, h)
        push_row!(out, oname, h, pt, se, nw, mlo, mhi, wlo, whi, dlo, dhi;
                  note = "05-style binary z")
    end
end

CSV.write(joinpath(DIAG_OUT, "B_bootstrap_harmonisation.csv"), out)
println("wrote B_bootstrap_harmonisation.csv  seed=$BOOT_SEED nboot=$NBOOT_B")
println("script: diagnostics/task_B.jl")
