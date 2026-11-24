"""
Phase numbering of the 3D model.
"""
const PHASES = (
    air = 0,
    quaternary = 1,              # unconsolidated deposits (swissBEDROCK TMUD)
    consolidated_sediments = 2,  # Alpine sedimentary cover/nappes (3D-ALPS; Morcles sedimentary units)
    upper_crust = 3,
    lower_crust = 4,
    mantle = 5,                  # below the European Moho (Spada et al. 2013)
    adriatic_mantle = 6,         # Adriatic mantle wedge (Ivrea body): below Adriatic Moho, above European Moho
)
const PHASE_NAMES = Dict(v => k for (k, v) in pairs(PHASES))

"""
    model_interfaces(x, y, topo, tmud) -> NamedTuple of matrices (km / kg m⁻³)

Layer interfaces on the LV95 grid `x`, `y` from the combined data sets:
- `topo`    : surface (swissALTIRegio)
- `bed`     : bedrock = topo - Quaternary thickness (swissBEDROCK; 0 thickness outside CH)
- `base_cs` : base of consolidated sediments = min(bedrock, 3D-ALPS top of upper crust)
- `moho_eu`, `moho_ad` : European / Adriatic Moho (Spada et al. 2013)
- `conrad`  : upper/lower crust boundary, placed at the 3D-ALPS fraction of the crystalline crust
             thickness between `base_cs` and the Spada Moho
- `rho_*`   : 3D-ALPS densities of consolidated sediments, upper/lower crust, mantle
"""
function model_interfaces(x, y, topo, tmud)
    Ee, Ne, Ze = read_moho_spada(1)
    Ea, Na, Za = read_moho_spada(2)
    moho_eu = interpolate_scattered(Ee, Ne, Ze, x, y; radius = 7.5)
    moho_ad = interpolate_scattered(Ea, Na, Za, x, y; radius = 7.5, maxdist = 5.0)
    a = alps3d_on_grid(x, y)
    t = replace(topo, NaN => 0.0)
    q = clamp.(replace(tmud, NaN => 0.0), 0.0, Inf)
    bed = t .- q
    base_cs = max.(min.(bed, a.top_uc), moho_eu .+ 2.0)
    f = clamp.((a.top_uc .- a.top_lc) ./ (a.top_uc .- a.top_lm), 0.05, 0.95)
    conrad = base_cs .- f .* (base_cs .- moho_eu)
    return (; topo = t, bed, quat_thickness = q, base_cs, conrad, moho_eu, moho_ad,
            rho_cs = a.rho_cs, rho_uc = a.rho_uc, rho_lc = a.rho_lc, rho_lm = a.rho_lm)
end

"""
    build_valais_model(r = valais_region(); dx = 0.25, dx_surface = 0.1, z = nothing,
                       fine_band = :auto, dz_fine = 0.025, dz_max_up = 0.25, dz_max_down = 2.0,
                       growth = 1.08, velocities = DEFAULT_VELOCITIES, morcles = true,
                       temperature = true, air_value = NaN, outdir = datadir("output"), write = true)

Build the Valais 3D geological / seismic velocity model on a rectilinear Cartesian LV95 grid (km):
uniform horizontal spacing `dx`, vertical levels refined (`dz_fine`) in the elevation band of the
Quaternary valley fill (see [`refined_levels`](@ref); `fine_band = :auto` derives it from
swissBEDROCK, or pass a tuple `(zlow, zhigh)`; pass `z` to give the levels explicitly).

Data: swissALTIRegio topography, swissBEDROCK Quaternary thickness, swisstopo tectonic map (surface),
3D-ALPS layers & densities, Spada et al. (2013) Moho, 3D-ALPS-TR temperature, Morcles nappe model.

Returns `(; surface, model, z, tecto)`. `model` is a `CartData` with fields `Phases` (see
[`PHASES`](@ref)), `Vp`, `Vs` (km/s), `Rho` (kg/m³), `Temperature` (°C) and `MorclesUnit`
(index into [`MORCLES_UNITS`](@ref), 0 outside). In air cells (`Phases == 0`) Vp, Vs, Rho and
Temperature are set to `air_value` (default `NaN`, which ParaView can render transparent via the
color map's "NaN Opacity"; use e.g. `air_value = 0` for solvers that do not accept NaN).
"""
function build_valais_model(r::ValaisRegion = valais_region(); dx = 0.25, dx_surface = 0.1, z = nothing,
                            fine_band = :auto, dz_fine = 0.025, dz_max_up = 0.25, dz_max_down = 2.0,
                            growth = 1.08, velocities = DEFAULT_VELOCITIES, morcles = true,
                            temperature = true, air_value = NaN, outdir = datadir("output"), write = true)
    # ---------------- surface (high resolution) ----------------
    xs, ys, topo_s = fetch_topography(r; dx = dx_surface)
    _, _, bem_s, tmud_s = fetch_swissbedrock(r; dx = dx_surface)
    tecto = tectonic_map_on_grid(xs, ys)
    I_s = model_interfaces(xs, ys, topo_s, tmud_s)
    Xs, Ys, Zs = xyz_grid(xs, ys, 0.0)
    Zs[:, :, 1] .= topo_s
    s3(a) = reshape(a, size(a)..., 1)
    surface = CartData(Xs, Ys, Zs, (
        Topography = s3(topo_s), BedrockElevation = s3(bem_s), QuaternaryThickness = s3(tmud_s),
        TectoLevel1 = s3(tecto.level1), TectoLevel2 = s3(tecto.level2), TectoUnit = s3(tecto.unit),
        QuaternaryMap = s3(Int.(tecto.quaternary)),
        BaseConsolidatedSediments = s3(I_s.base_cs), UpperLowerCrust = s3(I_s.conrad),
        MohoEuropean = s3(I_s.moho_eu), MohoAdriatic = s3(I_s.moho_ad),
    ))

    # ---------------- 3D grid ----------------
    xg, yg, topo = fetch_topography(r; dx = dx)
    _, _, _, tmud = fetch_swissbedrock(r; dx = dx)
    I = model_interfaces(xg, yg, topo, tmud)
    if z === nothing
        band = fine_band === :auto ? fine_band_from_bedrock(topo, tmud) : fine_band
        z = refined_levels(r.z[1], r.z[2]; fine = band, dz_fine, dz_max_up, dz_max_down, growth)
        @info "Vertical grid: $(length(z)) levels, dz = $(dz_fine) km in [$(round(band[1], digits=3)), $(round(band[2], digits=3))] km"
    end
    nx, ny, nz = length(xg), length(yg), length(z)
    @info "3D grid: $(nx) × $(ny) × $(nz) = $(round(nx*ny*nz/1e6, digits=1)) M cells"

    mm = morcles ? read_morcles() : nothing
    Phases = zeros(Int8, nx, ny, nz)
    Vp = zeros(Float32, nx, ny, nz); Vs = zeros(Float32, nx, ny, nz); Rho = zeros(Float32, nx, ny, nz)
    Morc = zeros(Int8, nx, ny, nz)
    V = velocities
    Threads.@threads for j in 1:ny
        for i in 1:nx, k in 1:nz
            zk = z[k]
            t = I.topo[i, j]
            zk > t && continue                         # air (set to air_value below)
            d = t - zk                                 # depth below surface
            mu = mm === nothing ? Int8(0) : morcles_unit_at(mm, xg[i], yg[j], zk)
            mad = I.moho_ad[i, j]; meu = I.moho_eu[i, j]
            ph, law, rho = if zk > I.bed[i, j]
                PHASES.quaternary, V.quaternary, NaN
            elseif mu > 0
                mu == 1 ? (PHASES.upper_crust, V.upper_crust, I.rho_uc[i, j]) :
                          (PHASES.consolidated_sediments, V.consolidated_sediments, I.rho_cs[i, j])
            elseif zk > I.base_cs[i, j]
                PHASES.consolidated_sediments, V.consolidated_sediments, I.rho_cs[i, j]
            elseif !isnan(mad) && mad > meu && meu < zk <= mad
                PHASES.adriatic_mantle, V.adriatic_mantle, I.rho_lm[i, j]
            elseif zk > I.conrad[i, j]
                PHASES.upper_crust, V.upper_crust, I.rho_uc[i, j]
            elseif zk > meu
                PHASES.lower_crust, V.lower_crust, I.rho_lc[i, j]
            else
                PHASES.mantle, V.mantle, I.rho_lm[i, j]
            end
            vp, vs = velocity(law, d)
            Phases[i, j, k] = ph
            Vp[i, j, k] = vp; Vs[i, j, k] = vs
            Rho[i, j, k] = isnan(law.rho) ? rho : law.rho
            zk <= I.bed[i, j] && (Morc[i, j, k] = mu)
        end
    end
    air = Phases .== PHASES.air
    for A in (Vp, Vs, Rho)
        A[air] .= Float32(air_value)
    end
    fields = (Phases = Phases, Vp = Vp, Vs = Vs, Rho = Rho, MorclesUnit = Morc)
    if temperature
        T = Float32.(temperature_3dalps(xg, yg, z))
        T[air] .= Float32(air_value)
        fields = merge(fields, (Temperature = T,))
    end
    X, Y, Zg = xyz_grid(xg, yg, z)
    model = CartData(X, Y, Zg, fields)

    if write
        mkpath(outdir)
        write_paraview(surface, joinpath(outdir, "valais_surface"))
        write_paraview(model, joinpath(outdir, "valais_model"))
        _write_legend(joinpath(outdir, "legend.csv"), tecto)
        writedlm(joinpath(outdir, "z_levels_km.txt"), z)
    end
    return (; surface, model, z, tecto)
end

function _write_legend(file, t::TectonicMap)
    open(file, "w") do io
        println(io, "field,code,name")
        for (field, names) in (("TectoLevel1", t.names1), ("TectoLevel2", t.names2), ("TectoUnit", t.names_unit))
            for (i, n) in enumerate(names)
                println(io, field, ",", i, ",\"", n, "\"")
            end
        end
        for (n, p) in pairs(PHASES)
            println(io, "Phases,", p, ",", n)
        end
        for (i, n) in enumerate(MORCLES_UNITS)
            println(io, "MorclesUnit,", i, ",", n)
        end
    end
end

"""
    cartdata_to_geodata(d::CartData) -> GeoData

Convert a `CartData` in LV95 km coordinates to a `GeoData` (lon, lat, depth in km), e.g. to
combine it with seismic tomography or other GeophysicalModelGenerator datasets.
"""
function cartdata_to_geodata(d::CartData)
    lon, lat = lv95_to_lonlat(NumValue(d.x), NumValue(d.y))
    return GeoData(lon, lat, NumValue(d.z), d.fields)
end
