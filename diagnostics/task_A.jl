# diagnostics/task_A.jl — Effective sample and event count
include(joinpath(@__DIR__, "common.jl"))

P = load_panel()
ev = CSV.read(EVENTS_PATH, DataFrame)
ev.ds = ym.(ev.date)

function z_variant(name::String, drop2022::Bool)
    z = name == "z" ? copy(P.z_bin) :
        name == "z_w" ? copy(P.z_w) :
        name == "z_major" ? copy(P.z_maj) :
        error(name)
    drop2022 && (z = apply_drop!(z, P.ds, DROP_2022))
    return z
end

rows_out = DataFrame(
    outcome = String[], instrument = String[], drop_2022 = Bool[],
    h = Int[], N_obs = Int[], N_z_gt0 = Int[], N_z_eq1 = Int[],
    N_z_eq05 = Int[], df = Float64[], n_blocks = Int[]
)

event_months_h0 = DataFrame(
    outcome = String[], instrument = String[], drop_2022 = Bool[],
    month = String[], z_value = Float64[]
)

for outcome in ("q_aus", "price")
    y = outcome == "q_aus" ? P.lq : P.lp
    for inst in ("z", "z_w", "z_major")
        for drop in (false, true)
            z = z_variant(inst, drop)
            for h in 0:18
                rows = rows_fs(y, z, P.dates, P.lq, P.lp, P.rea, h)
                β, se, n, df = lp_fs(y, z, P.dates, P.lq, P.lp, P.rea, h)
                zv = z[rows]
                n_gt0 = count(>(0), zv)
                n_1   = count(==(1.0), zv)
                n_05  = count(x -> abs(x - 0.5) < 1e-9, zv)
                # number of distinct 12-month blocks intersecting estimation rows
                if isempty(rows)
                    nb = 0
                else
                    tmin, tmax = extrema(rows)
                    nb = max(0, cld(tmax - tmin + 1, BLOCK))
                end
                push!(rows_out, (outcome, inst, drop, h, n, n_gt0, n_1, n_05,
                                 isfinite(df) ? df : NaN, nb))
            end
            # h=0 event months that ENTER the regression
            rows0 = rows_fs(y, z, P.dates, P.lq, P.lp, P.rea, 0)
            for t in rows0
                z[t] > 0 || continue
                push!(event_months_h0, (outcome, inst, drop, string(P.ds[t]), z[t]))
            end
        end
    end
end

CSV.write(joinpath(DIAG_OUT, "A_sample_counts.csv"), rows_out)
CSV.write(joinpath(DIAG_OUT, "A_h0_event_months.csv"), event_months_h0)
println("wrote A_sample_counts.csv and A_h0_event_months.csv")
println("script: diagnostics/task_A.jl")
