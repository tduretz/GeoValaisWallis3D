# Gravity anomalies (swisstopo) and the gravitational attraction of the topography.
#
# Sources (all from the federal geoportal, https://map.geo.admin.ch):
# - Gravimetric Atlas of Switzerland 1:100 000, measurement points
#   (layer ch.swisstopo.geologie-gravimetrischer_atlas.messpunkte): station height and complete
#   Bouguer anomaly (terrain-corrected, reduction density 2670 kg/m³).
# - Gravimetric Map of Switzerland 1:500 000, Bouguer anomalies
#   (layer ch.swisstopo.geologie-geodaesie-bouguer_anomalien): contour lines.
# - Swiss National Gravity Network LSN2004 (layer ch.swisstopo.landesschwerenetz): absolute
#   gravity g [mGal] and LN02 height of each station.
# swisstopo publishes no free-air anomaly map. Here it is computed (a) exactly at the LSN2004
# stations, FA = g − γ(φ, H), and (b) on the model surface as FA = BA + g_topo, where g_topo is the
# attraction of the topographic masses above sea level (2670 kg/m³), i.e. the complete Bouguer
# correction that was subtracted to get BA. (a) is used to check (b).

const GEOADMIN_IDENTIFY = "https://api3.geo.admin.ch/rest/services/api/MapServer/identify"
const GRAVITY_ATLAS_LAYER = "ch.swisstopo.geologie-gravimetrischer_atlas.messpunkte"
const BOUGUER_500K_LAYER = "ch.swisstopo.geologie-geodaesie-bouguer_anomalien"
const GRAVITY_NETWORK_LAYER = "ch.swisstopo.landesschwerenetz"
const G_NEWTON = 6.674e-11
const RHO_BOUGUER = 2670.0           # reduction density of the Swiss Bouguer anomalies [kg/m³]
const R_EARTH = 6371.0               # [km]

_bbox_lv95(r, margin) = string.(round.(Int, 1e3 .* (r.E[1] - margin, r.N[1] - margin, r.E[2] + margin, r.N[2] + margin)))

# All features of a geo.admin.ch layer in a box, through the identify service (paged, 200 per
# request). The JSON answer {"results": [GeoJSON features]} is turned into a FeatureCollection and
# read with GDAL. Returns a vector of (vertex coordinates [m], properties Dict).
function _geoadmin_features(layer, r, margin)
    feats = Tuple{Vector{NTuple{2,Float64}},Dict{String,Any}}[]
    bbox = join(_bbox_lv95(r, margin), ",")
    off = 0
    while true
        url = "$(GEOADMIN_IDENTIFY)?geometryType=esriGeometryEnvelope&geometry=$(bbox)&sr=2056" *
              "&layers=all:$(layer)&tolerance=0&returnGeometry=true&geometryFormat=geojson&limit=200&offset=$(off)"
        txt = String(take!(Downloads.download(url, IOBuffer())))
        txt = replace(txt, r"^\s*\{\s*\"results\"\s*:" => "{\"type\": \"FeatureCollection\", \"features\":"; count = 1)
        n = 0
        tmp = tempname() * ".geojson"
        write(tmp, txt)
        try
            AG.read(tmp) do ds
                for f in AG.getlayer(ds, 0)
                    g = AG.getgeom(f)
                    pts = NTuple{2,Float64}[]
                    _collect_points!(pts, g)
                    props = Dict{String,Any}(AG.getname(AG.getfielddefn(f, i)) => AG.getfield(f, i) for i in 0:AG.nfield(f)-1)
                    push!(feats, (pts, props))
                    n += 1
                end
            end
        finally
            rm(tmp; force = true)
        end
        n < 200 && break
        off += 200
    end
    return feats
end

function _collect_points!(pts, g)
    t = AG.getgeomtype(g)
    if t == AG.wkbPoint || t == AG.wkbPoint25D
        push!(pts, (AG.getx(g, 0), AG.gety(g, 0)))
    elseif t in (AG.wkbLineString, AG.wkbLineString25D)
        for i in 0:AG.ngeom(g)-1
            push!(pts, (AG.getx(g, i), AG.gety(g, i)))
        end
    else                                             # multi-geometries
        for i in 0:AG.ngeom(g)-1
            _collect_points!(pts, AG.getgeom(g, i))
        end
    end
    return pts
end

_csvcache(name, r, margin) = datadir("raw", "$(name)_$(_region_tag(r, 1.0))_m$(round(Int, margin))km.csv")

function _cached_table(build, file, header)
    if !isfile(file)
        rows = build()
        open(file, "w") do io
            println(io, join(header, ","))
            for row in rows
                println(io, join(row, ","))
            end
        end
        @info "Saved $(file)"
    end
    d = readdlm(file, ',', Any; skipstart = 1)
    return d
end

"""
    fetch_gravity_atlas(r = valais_region(); margin = 10) -> (; E, N, h, bouguer, name)

Measurement points of the **Gravimetric Atlas of Switzerland 1:100 000** (swisstopo), within the
region extended by `margin` km: LV95 coordinates `E`, `N` [km], station height `h` [km] and complete
Bouguer anomaly `bouguer` [mGal] (density 2670 kg/m³). Cached as CSV in `data/raw`.
"""
function fetch_gravity_atlas(r::ValaisRegion = valais_region(); margin = 10)
    mkpath(datadir("raw"))
    d = _cached_table(_csvcache("gravity_atlas_points", r, margin), ["E_m", "N_m", "h_m", "bouguer_mGal", "name"]) do
        @info "Downloading Gravimetric Atlas measurement points from api3.geo.admin.ch ..."
        [(p[1][1]..., p[2]["altitude"], p[2]["bouguerano"], p[2]["stationnam"]) for p in _geoadmin_features(GRAVITY_ATLAS_LAYER, r, margin)]
    end
    E, N, h, ba = Float64.(d[:, 1]) ./ 1e3, Float64.(d[:, 2]) ./ 1e3, Float64.(d[:, 3]) ./ 1e3, Float64.(d[:, 4])
    # a few points carry 0.0 as a placeholder for a missing value (e.g. 060-1081 near Leuk, whose
    # neighbours are at −130 mGal): drop a 0.0 that is > 30 mGal from the median of the points within 5 km
    bad = [ba[i] == 0 && (nb = [ba[j] for j in eachindex(ba) if j != i && hypot(E[j] - E[i], N[j] - N[i]) < 5];
                          !isempty(nb) && abs(median(nb)) > 30) for i in eachindex(ba)]
    any(bad) && @info "Gravimetric Atlas: dropped $(count(bad)) point(s) with a placeholder value 0: $(join(string.(d[bad, 5]), ", "))"
    keep = .!bad
    return (E = E[keep], N = N[keep], h = h[keep], bouguer = ba[keep], name = string.(d[keep, 5]))
end

"""
    fetch_bouguer_500k(r = valais_region(); margin = 10) -> (; E, N, bouguer)

Vertices of the contour lines of the **Gravimetric Map of Switzerland 1:500 000 (Bouguer
anomalies)** (swisstopo), with the contour value [mGal]. Cached as CSV in `data/raw`.
"""
function fetch_bouguer_500k(r::ValaisRegion = valais_region(); margin = 10)
    mkpath(datadir("raw"))
    d = _cached_table(_csvcache("bouguer500k_contours", r, margin), ["E_m", "N_m", "bouguer_mGal"]) do
        @info "Downloading the 1:500 000 Bouguer anomaly contours from api3.geo.admin.ch ..."
        [(p..., f[2]["et_fromatt"]) for f in _geoadmin_features(BOUGUER_500K_LAYER, r, margin) for p in f[1]]
    end
    return (E = Float64.(d[:, 1]) ./ 1e3, N = Float64.(d[:, 2]) ./ 1e3, bouguer = Float64.(d[:, 3]))
end

"""
    fetch_gravity_network(r = valais_region(); margin = 0) -> (; E, N, H, lat, g, name)

Stations of the **Swiss National Gravity Network LSN2004** (swisstopo): LV95 coordinates [km], LN02
height `H` [km], latitude [°] and absolute gravity `g` [mGal]. Cached as CSV in `data/raw`.
"""
function fetch_gravity_network(r::ValaisRegion = valais_region(); margin = 0)
    mkpath(datadir("raw"))
    d = _cached_table(_csvcache("gravity_network_lsn2004", r, margin), ["E_m", "N_m", "H_m", "lat_deg", "g_mGal", "name"]) do
        @info "Downloading National Gravity Network (LSN2004) stations from api3.geo.admin.ch ..."
        [(p[2]["y_lv95"], p[2]["x_lv95"], p[2]["h_ln02"], p[2]["lat_etrs"], p[2]["gravity"],
          replace(string(p[2]["nr_lsn2004"], " ", p[2]["name"]), "," => " "))
         for p in _geoadmin_features(GRAVITY_NETWORK_LAYER, r, margin)]
    end
    return (E = Float64.(d[:, 1]) ./ 1e3, N = Float64.(d[:, 2]) ./ 1e3, H = Float64.(d[:, 3]) ./ 1e3,
            lat = Float64.(d[:, 4]), g = Float64.(d[:, 5]), name = string.(d[:, 6]))
end

# normal_gravity(lat, h) [mGal, h in km]: GRS80, from GeophysicalModelGenerator (≥ 0.10.1).

# Vertical attraction [mGal] at the origin of a rectangular prism x1..x2, y1..y2, z1..z2 (m, z up),
# density ρ (Nagy et al. 2000). Positive for a mass below the origin.
function _prism_gz(x1, x2, y1, y2, z1, z2, ρ)
    s = 0.0
    for (i, x) in enumerate((x1, x2)), (j, y) in enumerate((y1, y2)), (k, z) in enumerate((z1, z2))
        R = sqrt(x^2 + y^2 + z^2)
        μ = isodd(i + j + k) ? -1.0 : 1.0
        t = 0.0
        abs(y + R) > 1e-9 && (t += x * log(y + R))
        abs(x + R) > 1e-9 && (t += y * log(x + R))
        abs(z) > 1e-9 && (t -= z * atan(x * y / (z * R)))
        s += μ * t
    end
    return G_NEWTON * ρ * s * 1e5
end

# Vertical attraction [mGal] of a vertical line mass (column of section A [m²], density ρ) from z1
# to z2 [m, z up], at horizontal distance r [m] and height zs [m].
@inline _column_gz(A, ρ, r2, zs, z1, z2) =
    G_NEWTON * ρ * A * (1 / sqrt(r2 + (zs - z2)^2) - 1 / sqrt(r2 + (zs - z1)^2)) * 1e5

"""
    topographic_attraction(Es, Ns, hs, r; dx = 0.25, near = 10.0, dx_far = 1.0, margin = 100.0,
                           rho = 2670.0) -> g [mGal]

Vertical attraction of the topographic masses between sea level and the swissALTIRegio surface
(density `rho`), at the points `(Es, Ns, hs)` [km]: the complete Bouguer correction (plate + terrain).
Masses closer than `near` km (square window) use the DEM at `dx` km, those farther away the DEM at
`dx_far` km out to `margin` km around the region `r`; Earth curvature is included. The 3 × 3 cells
around each point are exact prisms, the others vertical line masses.
"""
function topographic_attraction(Es, Ns, hs, r::ValaisRegion; dx = 0.25, near = 10.0, dx_far = 1.0,
                                margin = 100.0, rho = RHO_BOUGUER)
    nf = round(Int, dx_far / dx)
    nf * dx ≈ dx_far || error("dx_far must be a multiple of dx")
    nn = round(Int, near / dx_far)
    grow(m) = ValaisRegion((floor(r.E[1] - m), ceil(r.E[2] + m)), (floor(r.N[1] - m), ceil(r.N[2] + m)), r.z)
    rf = grow(margin)
    xf, yf, tf0 = fetch_topography(rf; dx = dx_far)
    rn = grow(near + dx_far)                                    # fine DEM covers every near window
    xn, yn, tn0 = fetch_topography(rn; dx = dx)
    tf = replace(tf0, NaN => 0.0); tn = replace(tn0, NaN => 0.0)   # (new names: no boxed captures)
    # fine cells inside far cell (I, J): (I - 1) * nf + 1 + offset
    ox = round(Int, (rf.E[1] - rn.E[1]) / dx); oy = round(Int, (rf.N[1] - rn.N[1]) / dx)
    Af, An = (dx_far * 1e3)^2, (dx * 1e3)^2
    out = zeros(length(Es))
    Threads.@threads for p in eachindex(Es)
        E, N, zs = Es[p], Ns[p], hs[p] * 1e3
        I0 = clamp(floor(Int, (E - xf[1]) / dx_far + 1.5), 1, length(xf))    # far cell of the point
        J0 = clamp(floor(Int, (N - yf[1]) / dx_far + 1.5), 1, length(yf))
        g = 0.0
        @inbounds for J in eachindex(yf), I in eachindex(xf)
            (abs(I - I0) <= nn && abs(J - J0) <= nn) && continue
            ddx, ddy = (xf[I] - E) * 1e3, (yf[J] - N) * 1e3
            r2 = ddx^2 + ddy^2
            c = r2 / (2R_EARTH * 1e3)                                       # curvature drop [m]
            t = tf[I, J] * 1e3
            g += _column_gz(Af, rho, r2, zs, min(t, 0.0) - c, max(t, 0.0) - c) * (t >= 0 ? 1 : -1)
        end
        i1 = (I0 - nn - 1) * nf + 1 + ox; i2 = (I0 + nn) * nf + ox
        j1 = (J0 - nn - 1) * nf + 1 + oy; j2 = (J0 + nn) * nf + oy
        @inbounds for j in max(j1, 1):min(j2, length(yn)), i in max(i1, 1):min(i2, length(xn))
            ddx, ddy = (xn[i] - E) * 1e3, (yn[j] - N) * 1e3
            t = tn[i, j] * 1e3
            lo, hi = min(t, 0.0), max(t, 0.0)
            sgn = t >= 0 ? 1 : -1
            if abs(ddx) < 1.5dx * 1e3 && abs(ddy) < 1.5dx * 1e3
                h = dx * 500
                g += sgn * _prism_gz(ddx - h, ddx + h, ddy - h, ddy + h, lo - zs, hi - zs, rho)
            else
                g += sgn * _column_gz(An, rho, ddx^2 + ddy^2, zs, lo, hi)
            end
        end
        out[p] = g
    end
    return out
end

"""
    gravity_maps(r = valais_region(); dx = 0.25, radius = 5.0, maxdist = 5.0, kwargs...)
        -> (; x, y, topo, bouguer, bouguer500k, gtopo, freeair, freeair500k, atlas, network)

Gravity anomaly maps on the surface grid of the model (cell centres at spacing `dx` km):
- `bouguer`: complete Bouguer anomaly [mGal], inverse-distance interpolation (within `radius` km)
  of the Gravimetric Atlas points; `NaN` farther than `maxdist` km from any point.
- `bouguer500k`: the same from the vertices of the 1:500 000 contour lines (smoother; also
  covers the Italian side, where the atlas points are sparse).
- `gtopo`: attraction of the topography (`topographic_attraction`).
- `freeair = bouguer + gtopo`, `freeair500k = bouguer500k + gtopo`: free-air anomaly at the
  topographic surface [mGal].
`network` holds the LSN2004 stations with their exact free-air anomaly `fa_obs`, the value of the
map method at the station `fa_map`, and the difference. Other keywords go to `topographic_attraction`.
"""
function gravity_maps(r::ValaisRegion = valais_region(); dx = 0.25, radius = 5.0, maxdist = 5.0, kwargs...)
    x, y, topo = fetch_topography(r; dx = dx)
    atlas = fetch_gravity_atlas(r)
    c500 = fetch_bouguer_500k(r)
    bouguer = interpolate_scattered(atlas.E, atlas.N, atlas.bouguer, x, y; radius = radius, maxdist = maxdist)
    bouguer500k = interpolate_scattered(c500.E, c500.N, c500.bouguer, x, y; radius = radius, maxdist = maxdist)
    X = [xi for xi in x, _ in y]; Y = [yj for _ in x, yj in y]
    @info "Attraction of the topography at $(length(X)) surface points ..."
    gtopo = reshape(topographic_attraction(vec(X), vec(Y), replace(vec(topo), NaN => 0.0), r; dx = dx, kwargs...), size(X))
    freeair = bouguer .+ gtopo
    freeair500k = bouguer500k .+ gtopo
    # check against the exact free-air anomaly at the national gravity network stations
    lsn = fetch_gravity_network(r)
    fa_obs = lsn.g .- normal_gravity.(lsn.lat, lsn.H)
    ba_st = [interpolate_scattered(atlas.E, atlas.N, atlas.bouguer, [e], [n]; radius = radius, maxdist = maxdist)[1]
             for (e, n) in zip(lsn.E, lsn.N)]
    fa_map = ba_st .+ topographic_attraction(lsn.E, lsn.N, lsn.H, r; dx = dx, kwargs...)
    network = merge(lsn, (; fa_obs, fa_map, misfit = fa_map .- fa_obs))
    return (; x, y, topo, bouguer, bouguer500k, gtopo, freeair, freeair500k, atlas, network)
end

"""
    write_gravity_paraview(gm, outdir = datadir("output"); lift = 0.02)

Write the gravity maps of `gravity_maps` for ParaView, to overlay on the model:
- `valais_gravity.vts`: the topographic surface (raised by `lift` km so that it is drawn above
  `valais_surface.vts`) with fields `Bouguer_mGal`, `Bouguer500k_mGal`, `FreeAir_mGal`,
  `FreeAir500k_mGal`, `TopoAttraction_mGal`, `Topography_km`;
- `valais_gravity_atlas.vtu`, `valais_gravity_network.vtu`: the stations (points).
"""
function write_gravity_paraview(gm, outdir = datadir("output"); lift = 0.02)
    mkpath(outdir)
    X, Y, Z = xyz_grid(gm.x, gm.y, 0.0)
    Z[:, :, 1] .= gm.topo .+ lift
    s(a) = reshape(Float64.(a), size(X))
    surf = CartData(X, Y, Z, (Bouguer_mGal = s(gm.bouguer), Bouguer500k_mGal = s(gm.bouguer500k),
                              FreeAir_mGal = s(gm.freeair), FreeAir500k_mGal = s(gm.freeair500k),
                              TopoAttraction_mGal = s(gm.gtopo),
                              Topography_km = s(gm.topo)))
    write_paraview(surf, joinpath(outdir, "valais_gravity"))
    a = gm.atlas
    write_paraview(CartData(a.E, a.N, a.h .+ lift, (Bouguer_mGal = a.bouguer, Height_km = a.h)),
                   joinpath(outdir, "valais_gravity_atlas"); PointsData = true)
    n = gm.network
    write_paraview(CartData(n.E, n.N, n.H .+ lift, (FreeAir_obs_mGal = n.fa_obs, FreeAir_map_mGal = n.fa_map,
                                                    Misfit_mGal = n.misfit, g_mGal = n.g)),
                   joinpath(outdir, "valais_gravity_network"); PointsData = true)
    return outdir
end

"""
    model_gravity(x, y, z, rho; air = rho .< 10, dx = 1.0, dz = 0.25, rho_topo = 2670.0,
                  ztop = nothing, outdir = tempdir())
        -> (; x, y, z_obs, dg, gradx, grady, drho, xg, yg, zg)

Bouguer anomaly [mGal] of a 3D density model, computed with GeophysicalModelGenerator's
`voxel_grav` (point-mass voxels, observation plane at the top of the grid).

`x`, `y`, `z` [km] are the coordinate vectors of the model grid (z may be stretched), `rho`
[kg/m³] the density (`air`: mask of air cells). The density is first averaged onto a regular grid
(`dx`, `dz` km; `voxel_grav` needs one) and turned into the density contrast that a complete
Bouguer anomaly sees:
- between sea level and the surface: `rho − rho_topo` (the topography was removed with `rho_topo`);
- below sea level: `rho` minus the mean of its depth slice (a laterally homogeneous reference Earth;
  this sets the anomaly only up to a constant);
- air: 0.
`voxel_grav` ignores masses above z = 0, so the grid is shifted down by `ztop` (default: top of the
model, rounded to `dz`): every mass is included, and the anomaly is obtained on the horizontal plane
`z_obs = ztop` km a.s.l. (above the topography, as in an airborne survey).
"""
function model_gravity(x, y, z, rho; air = rho .< 10, dx = 1.0, dz = 0.25, rho_topo = RHO_BOUGUER,
                       ztop = nothing, outdir = tempdir())
    # ---- regular grid: x, y nodes every dx (block mean of model columns), z nodes every dz
    xg = collect(x[1]:dx:x[end]); yg = collect(y[1]:dx:y[end])
    ztop = isnothing(ztop) ? z[1] + floor((z[end] - z[1]) / dz) * dz : ztop
    zg = collect(ztop:-dz:z[1])[end:-1:1]
    # model levels as layers: bounds halfway between levels
    zb = vcat(z[1], (z[1:end-1] .+ z[2:end]) ./ 2, z[end])
    dr = zeros(length(xg), length(yg), length(zg))
    # vertical overlap weights between model layers and regular layers [zg - dz/2, zg + dz/2]
    W = zeros(length(zg), length(z))
    for k in eachindex(zg), l in eachindex(z)
        W[k, l] = max(0.0, min(zg[k] + dz / 2, zb[l + 1]) - max(zg[k] - dz / 2, zb[l]))
    end
    W ./= max.(sum(W, dims = 2), eps())
    Threads.@threads for j in eachindex(yg)
        jj = findall(v -> abs(v - yg[j]) <= dx / 2 + 1e-9, y)
        for i in eachindex(xg)
            ii = findall(v -> abs(v - xg[i]) <= dx / 2 + 1e-9, x)
            # contrast in each model cell of the block, then block mean
            col = zeros(length(z))
            for jm in jj, im in ii, l in eachindex(z)
                col[l] += air[im, jm, l] ? 0.0 : (z[l] > 0 ? rho[im, jm, l] - rho_topo : Float64(rho[im, jm, l]))
            end
            col ./= length(ii) * length(jj)
            dr[i, j, :] .= W * col
        end
    end
    # below sea level: subtract the mean of each depth slice (reference Earth)
    for k in eachindex(zg)
        zg[k] <= 0 && (dr[:, :, k] .-= sum(dr[:, :, k]) / length(dr[:, :, k]))
    end
    X, Y, Z = xyz_grid(xg, yg, zg .- ztop)
    out = joinpath(outdir, "voxel_grav_tmp")
    dg, gx, gy = voxel_grav(X, Y, Z, dr; refMod = zeros(length(zg)), lengthUnit = "km", outName = out, printing = false)
    rm(out * ".vts"; force = true)
    return (; x = xg, y = yg, z_obs = ztop, dg, gradx = gx, grady = gy, drho = dr, xg, yg, zg)
end

"""
    upward_continue(x, y, g, h) -> g_h

Upward continuation of a gridded field `g` (regular grid `x`, `y` [km]) by `h` km, by the Poisson
integral in the space domain. The kernel is normalised over the grid, so that the field is not
biased low near the edges.
"""
function upward_continue(x, y, g, h)
    dA = (x[2] - x[1]) * (y[2] - y[1])
    ok = .!isnan.(g)
    out = fill(NaN, size(g))
    Threads.@threads for j in eachindex(y)
        for i in eachindex(x)
            s = 0.0; w = 0.0
            @inbounds for jj in eachindex(y), ii in eachindex(x)
                ok[ii, jj] || continue
                r2 = (x[ii] - x[i])^2 + (y[jj] - y[j])^2 + h^2
                k = h * dA / (2π * r2^1.5)
                s += k * g[ii, jj]; w += k
            end
            out[i, j] = s / w
        end
    end
    return out
end

# ======================================================================================================
# Density inversion of gravity data (GeophysicalModelGenerator ≥ 0.10.1: voxel_grav + voxel_grav_adjoint)
# ======================================================================================================

# bounds of the layers around the nodes v (halfway between nodes; half a spacing beyond the end nodes)
function _node_bounds(v)
    m = (v[1:end-1] .+ v[2:end]) ./ 2
    return vcat(v[1] - (v[2] - v[1]) / 2, m), vcat(m, v[end] + (v[end] - v[end-1]) / 2)
end

"""
    density_on_grid(x, y, z, rho, air, xg, yg, zg; topo = nothing) -> (; rho, rockfrac)

Average a density model (`rho` [kg/m³] on the cell-centre vectors `x`, `y`, `z` [km], `air` mask)
onto the nodes `xg`, `yg`, `zg` of a coarser rectilinear grid, as used by `voxel_grav` (spacing may
vary). Every node gets the mean over the model cells in its dual cell (halfway to the neighbouring
nodes): `rockfrac` is the fraction of rock (not air), `rho` the mean density of the rock (NaN if
there is none).

Nodes outside the horizontal footprint of the model (lateral padding) take the density of the
nearest model columns below sea level and 2670 kg/m³ above. Their rock fraction is set by `topo`
(topography at the nodes, km), which is required if there are such nodes.
"""
function density_on_grid(x, y, z, rho, air, xg, yg, zg; topo = nothing)
    xl, xu = _node_bounds(xg); yl, yu = _node_bounds(yg); zl, zu = _node_bounds(zg)
    zb = vcat(z[1] - (z[2] - z[1]) / 2, (z[1:end-1] .+ z[2:end]) ./ 2, z[end] + (z[end] - z[end-1]) / 2)
    overlap(a, b, c, d) = max(0.0, min(b, d) - max(a, c))
    W = [overlap(zl[k], zu[k], zb[l], zb[l + 1]) for k in eachindex(zg), l in eachindex(z)]
    W ./= max.(sum(W, dims = 2), eps())
    inside(v, lo, hi) = lo <= v <= hi
    dx2, dy2 = (x[end] - x[1]) / (length(x) - 1) / 2, (y[end] - y[1]) / (length(y) - 1) / 2
    R = fill(NaN, length(xg), length(yg), length(zg)); F = zeros(size(R))
    Threads.@threads for j in eachindex(yg)
        yc = clamp(yg[j], y[1], y[end])
        jj = findall(v -> overlap(yl[j], yu[j], v - dy2, v + dy2) > 0, y)
        isempty(jj) && (jj = [argmin(abs.(y .- yc))])
        padj = !inside(yg[j], y[1] - dy2, y[end] + dy2)
        for i in eachindex(xg)
            xc = clamp(xg[i], x[1], x[end])
            ii = findall(v -> overlap(xl[i], xu[i], v - dx2, v + dx2) > 0, x)
            isempty(ii) && (ii = [argmin(abs.(x .- xc))])
            pad = padj || !inside(xg[i], x[1] - dx2, x[end] + dx2)
            # rock fraction and rock mass per model level, summed over the columns of the block
            nrock = zeros(length(z)); mass = zeros(length(z))
            for jm in jj, im in ii, l in eachindex(z)
                air[im, jm, l] && continue
                nrock[l] += 1; mass[l] += rho[im, jm, l]
            end
            n = length(ii) * length(jj)
            for k in eachindex(zg)
                f = sum(W[k, l] * nrock[l] for l in eachindex(z)) / n
                m = sum(W[k, l] * mass[l] for l in eachindex(z)) / n
                if pad
                    isnothing(topo) && error("nodes outside the model need `topo`")
                    F[i, j, k] = clamp((topo[i, j] - zl[k]) / (zu[k] - zl[k]), 0.0, 1.0)
                    R[i, j, k] = zg[k] > 0 ? RHO_BOUGUER : (f > 0 ? m / f : RHO_BOUGUER)
                else
                    F[i, j, k] = f
                    f > 0 && (R[i, j, k] = m / f)
                end
            end
        end
    end
    return (; rho = R, rockfrac = F)
end

"""
    bouguer_contrast(rho, rockfrac, zg; rho_topo = 2670.0) -> (; contrast, ref)

Density contrast seen by a complete Bouguer anomaly, at the nodes of [`density_on_grid`](@ref):
`rockfrac · (rho − rho_topo)` above sea level (the topography was removed with `rho_topo`; air
counts 0), `rho − ref` below, where `ref[k]` is the mean density of depth level `k` (a laterally
homogeneous reference Earth; the anomaly is then defined up to a constant).
"""
function bouguer_contrast(rho, rockfrac, zg; rho_topo = RHO_BOUGUER)
    C = zeros(size(rho)); ref = zeros(length(zg))
    for k in eachindex(zg)
        r = @view rho[:, :, k]
        if zg[k] > 0
            C[:, :, k] .= replace(rockfrac[:, :, k] .* (r .- rho_topo), NaN => 0.0)
        else
            ref[k] = mean(filter(!isnan, r))
            C[:, :, k] .= replace(r .- ref[k], NaN => 0.0)
        end
    end
    return (; contrast = C, ref)
end

"""
    invert_density(xg, yg, zg, contrast0, active, topo, stations, dobs, σ;
                   z0 = 1.0, β = 2.0, chi2 = 1.0, itmax = 300, regional = true,
                   backend = nothing, float_type = Float64)
        -> (; drho, dg0, dg, offset, iterations, misfit, rms_history)

Least-squares density inversion of a Bouguer anomaly with GeophysicalModelGenerator's
`voxel_grav` and `voxel_grav_adjoint`, solved by LSQR (Krylov.jl), as in the GMG documentation
(example "least-squares inversion").

- Grid: nodes `xg`, `yg`, `zg` [km]; `contrast0` is the density contrast of the starting model
  ([`bouguer_contrast`](@ref)), `active` the nodes whose density may change, `topo` [km] the
  surface at the `xg`, `yg` columns (for the depth weighting).
- Data: `stations = (E, N, h)` [km], `dobs` [mGal] the observed complete Bouguer anomaly there,
  `σ` its standard deviation [mGal] (a number or one per station).
- The unknown is the contrast change `drho` at the active nodes, with the depth weighting of Li &
  Oldenburg (1998), `w = √V (depth + z0)^(−β/2)` (`V` the node volume, `depth` below the surface):
  LSQR solves for `w .* drho`, which gives the model of minimum weighted norm. `regional = true`
  adds a plane `a + b·(E − Ē) + c·(N − N̄)` (unpenalised) for the long-wavelength field from outside
  the grid and for the reference Earth, which `refMod` sets only up to a constant.
- The iterations stop at the discrepancy principle, when the χ² misfit per datum reaches `chi2`
  (or after `itmax`): fitting further would fit the noise.
- `backend` / `float_type` are passed to `voxel_grav`, e.g. `MetalBackend()` with `Float32`. A GPU
  call that fails is repeated (up to 3 times).

Returns the contrast change `drho` (size of the grid), the anomaly at the stations of the starting
(`dg0`) and of the inverted model (`dg`, including the plane), the plane `offset` (mGal, one value
per station), the number of `iterations`, and the rms misfit [mGal] before/after (`misfit`) and at
every iteration (`rms_history`).
"""
function invert_density(xg, yg, zg, contrast0, active, topo, stations, dobs, σ; kwargs...)
    X, Y, Z = xyz_grid(xg, yg, zg)
    # volume of the dual cell of every node, and depth below the topography
    h(v) = (b = _node_bounds(v); b[2] .- b[1])
    V = h(xg) .* h(yg)' .* reshape(h(zg), 1, 1, :)
    depth = max.(reshape(topo, size(topo)..., 1) .- reshape(zg, 1, 1, :), 0.0)
    return _invert_density(CartData(X, Y, Z, (z = Z,)), V, depth, contrast0, active, stations, dobs, σ; kwargs...)
end

"""
    invert_density(Grid::CartData, contrast0, active, stations, dobs, σ; z0 = 1.0, β = 2.0, kwargs...)

[`invert_density`](@ref) on a curvilinear grid (e.g. [`curvilinear_density_model`](@ref)), whose top level is the
topography: the depth weighting uses the depth of every node below the top of its column and the volume of its dual
cell (horizontal spacing times the half-distances to the levels above and below). All other arguments and the result
are as for the rectilinear grid; no node lies in the air, so `active` usually is all `true`.
"""
function invert_density(Grid::CartData, contrast0, active, stations, dobs, σ; kwargs...)
    X, Y, Z = Grid.x.val, Grid.y.val, Grid.z.val
    h(v) = (b = _node_bounds(v); b[2] .- b[1])
    hz = similar(Z, Float64)
    for k in axes(Z, 3)
        lo = k > 1 ? (Z[:, :, k] .- Z[:, :, k - 1]) ./ 2 : 0.0
        hi = k < size(Z, 3) ? (Z[:, :, k + 1] .- Z[:, :, k]) ./ 2 : 0.0
        hz[:, :, k] .= lo .+ hi
    end
    V = h(X[:, 1, 1]) .* h(Y[1, :, 1])' .* hz
    depth = Z[:, :, end:end] .- Z
    return _invert_density(Grid, V, depth, contrast0, active, stations, dobs, σ; kwargs...)
end

function _invert_density(Grid, V, depth, contrast0, active, stations, dobs, σ;
                         z0 = 1.0, β = 2.0, chi2 = 1.0, itmax = 300, regional = true,
                         backend = nothing, float_type = Float64)
    X = Grid.x.val
    St = CartData(reshape(collect(Float64, stations[1]), :, 1, 1), reshape(collect(Float64, stations[2]), :, 1, 1),
                  reshape(collect(Float64, stations[3]), :, 1, 1), (h = reshape(collect(Float64, stations[3]), :, 1, 1),))
    nd = length(dobs)
    ref = zeros(size(X, 3))
    kw = (; stations = St, backend, float_type)
    # a GPU call can be discarded by the driver (e.g. Metal "victim of GPU error/recovery", when the GPU
    # resets because of another process): the operator is deterministic, so such a call is repeated
    function retry(f)
        for attempt in 1:3
            try
                return f()
            catch e
                (isnothing(backend) || attempt == 3) && rethrow()
                @warn "voxel_grav on the GPU failed (attempt $(attempt)), repeating" exception = (e, catch_backtrace())
            end
        end
    end
    fwd(ρ) = retry(() -> vec(voxel_grav(Grid, ρ; refMod = ref, rhoTol = 0, kw...)[1]))
    adj(r) = retry(() -> voxel_grav_adjoint(Grid, reshape(r, :, 1, 1); kw...))
    # depth weighting (Li & Oldenburg 1998) with the volume of the dual cell of every node
    w = (sqrt.(V) .* (depth .+ z0) .^ (-β / 2))[active]
    s = σ isa Number ? fill(Float64(σ), nd) : Float64.(σ)
    # regional plane, scaled so that LSQR hardly penalises it
    E, N = collect(Float64, stations[1]), collect(Float64, stations[2])
    P = regional ? hcat(ones(nd), (E .- mean(E)) ./ 50, (N .- mean(N)) ./ 50) : zeros(nd, 0)
    S = 100.0
    na, np = count(active), size(P, 2)
    scatter(m) = (ρ = zeros(size(X)); ρ[active] .= m; ρ)
    A = LinearMap{Float64}(p -> (fwd(scatter(p[1:na] ./ w)) .+ S .* (P * p[na+1:end])) ./ s,
                           r -> vcat(adj(r ./ s)[active] ./ w, S .* (P' * (r ./ s))), nd, na + np)
    dg0 = fwd(contrast0)
    b = (dobs .- dg0) ./ s
    target = sqrt(chi2 * nd)
    @info "density inversion: $(na) unknowns, $(nd) stations; starting rms misfit $(round(sqrt(mean(abs2, dobs .- dg0 .- mean(dobs .- dg0))), digits = 2)) mGal (mean removed)"
    t = time()
    p, stats = lsqr(A, b; itmax = itmax, history = true, atol = 0.0, btol = 0.0, etol = 0.0,
                    callback = ws -> (r = ws.stats.residuals; length(r) > 1 && r[end] <= target))
    @info "LSQR: $(stats.niter) iterations in $(round(time() - t, digits = 1)) s, normalised misfit $(round(stats.residuals[end] / sqrt(nd), digits = 3))"
    drho = scatter(p[1:na] ./ w)
    offset = S .* (P * p[na+1:end])
    dg = dg0 .+ fwd(drho) .+ offset
    rms(r) = sqrt(mean(abs2, r))
    return (; drho, dg0, dg, offset, iterations = stats.niter,
            misfit = (start = rms(dobs .- dg0 .- mean(dobs .- dg0)), final = rms(dobs .- dg)),
            rms_history = stats.residuals ./ sqrt(nd) .* mean(s))
end

"""
    interpolate_nodes(xg, yg, zg, F, x, y, z) -> array of size (length(x), length(y), length(z))

Trilinear interpolation of `F` given at the nodes of the rectilinear grid `xg`, `yg`, `zg` onto the
points of another rectilinear grid (constant extrapolation outside).
"""
function interpolate_nodes(xg, yg, zg, F, x, y, z)
    function wts(v, p)
        i = clamp(searchsortedlast(v, p), 1, length(v) - 1)
        return i, clamp((p - v[i]) / (v[i + 1] - v[i]), 0.0, 1.0)
    end
    out = zeros(eltype(F), length(x), length(y), length(z))
    wx = wts.(Ref(xg), x); wy = wts.(Ref(yg), y); wz = wts.(Ref(zg), z)
    Threads.@threads for k in eachindex(z)
        kk, tz = wz[k]
        for j in eachindex(y), i in eachindex(x)
            ii, tx = wx[i]; jj, ty = wy[j]
            v = zero(eltype(F))
            for (a, fa) in ((0, 1 - tx), (1, tx)), (c, fc) in ((0, 1 - ty), (1, ty)), (e, fe) in ((0, 1 - tz), (1, tz))
                v += fa * fc * fe * F[ii + a, jj + c, kk + e]
            end
            out[i, j, k] = v
        end
    end
    return out
end
