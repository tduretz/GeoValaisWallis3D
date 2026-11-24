const SPADA_URL(i) =
    "https://seafile.rlp.net/d/a50881f45aa34cdeb3c0/files/?p=%2FMoho_Map_Data-WesternAlps-SpadaETAL2013_Moho$(i).txt&dl=1"

"""
    download_moho_spada(; force = false) -> Vector of 3 paths

Download the Western Alps Moho maps of Spada et al. (2013, GJI) to `data/raw/`:
1 = European Moho, 2 = Adriatic Moho, 3 = Tyrrhenian/Ligurian Moho.
(Same files as used in the GeophysicalModelGenerator Moho tutorial.)
"""
function download_moho_spada(; force = false)
    mkpath(datadir("raw"))
    map(1:3) do i
        file = datadir("raw", "Spada2013_Moho$(i).txt")
        if force || !isfile(file)
            @info "Downloading Spada et al. (2013) Moho $(i) ..."
            Downloads.download(SPADA_URL(i), file)
        end
        file
    end
end

"""
    read_moho_spada(i) -> (E, N, z)

Read Spada et al. (2013) Moho surface `i` (1 European, 2 Adriatic, 3 Tyrrhenian), converted
to LV95 `E`, `N` (km) and Moho elevation `z` (km, negative below sea level).
"""
function read_moho_spada(i::Integer)
    file = download_moho_spada()[i]
    d = readdlm(file, Float64; comments = true, comment_char = '#')
    E, N = lonlat_to_lv95(d[:, 1], d[:, 2])
    return E, N, -d[:, 3]
end

"""
    interpolate_scattered(xp, yp, vp, x, y; radius = 7.5, maxdist = Inf, power = 2)

Inverse-distance-weighted interpolation of scattered points `(xp, yp, vp)` onto the regular
grid `x`, `y` (all km). Points within `radius` are used (falls back to the nearest point);
grid nodes farther than `maxdist` from any data point are set to `NaN`.
"""
function interpolate_scattered(xp, yp, vp, x, y; radius = 7.5, maxdist = Inf, power = 2)
    # keep only points near the grid
    m = max(radius, isfinite(maxdist) ? maxdist : radius)
    keep = (x[1] - m .<= xp .<= x[end] + m) .& (y[1] - m .<= yp .<= y[end] + m)
    xp, yp, vp = xp[keep], yp[keep], vp[keep]
    out = fill(NaN, length(x), length(y))
    isempty(xp) && return out
    Threads.@threads for j in eachindex(y)
        for i in eachindex(x)
            dmin, imin, wsum, vsum = Inf, 0, 0.0, 0.0
            @inbounds for k in eachindex(xp)
                d = hypot(xp[k] - x[i], yp[k] - y[j])
                d < dmin && (dmin = d; imin = k)
                if d <= radius
                    w = 1 / max(d, 1e-6)^power
                    wsum += w
                    vsum += w * vp[k]
                end
            end
            dmin > maxdist && continue
            out[i, j] = wsum > 0 ? vsum / wsum : vp[imin]
        end
    end
    return out
end
