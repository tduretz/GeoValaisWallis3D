# Seismic subsoil classes (Baugrundklassen / classes de sol de fondation, SIA 261).
#
# Federal dataset of the Federal Office for the Environment (BAFU/FOEN), layer
# ch.bafu.gefahren-baugrundklassen on map.geo.admin.ch. It contains the cantonal maps: for Valais
# the indicative 1:25 000 maps of the Rhône valley from Brig to Lake Geneva, made by CREALP for the
# canton (https://www.vs.ch/web/sdm/cartes-de-sol-de-fondation). They apply where no spectral
# seismic microzonation exists (Brig-Visp, Sion, Martigny, Chablais: PDF maps on the cantonal page).

const SUBSOIL_CLASSES_STAC = "https://data.geo.admin.ch/api/stac/v0.9/collections/ch.bafu.gefahren-baugrundklassen/items"
const SUBSOIL_CLASSES_URL =
    "https://data.geo.admin.ch/ch.bafu.gefahren-baugrundklassen/gefahren-baugrundklassen/gefahren-baugrundklassen_2056.shp.zip"

"""
    SUBSOIL_CLASSES

Seismic subsoil classes of SIA 261, in the order of their integer codes 1–6 (0 = not mapped):

| code | class | description (SIA 261, abridged) |
|---|---|---|
| 1 | A | rock, or rock-like ground with at most 5 m of soil at the surface |
| 2 | B | very dense gravel and sand, or very stiff clay, several tens of metres thick |
| 3 | C | normally consolidated, uncompacted gravel and sand and/or moraine, > 30 m thick |
| 4 | D | unconsolidated fine sand, silt and clay, > 30 m thick |
| 5 | E | surface alluvial layer of class C or D, 5–30 m thick, over stiffer ground of class A or B |
| 6 | F | structurally sensitive and organic deposits (e.g. peat, lake marl, landslide masses), > 10 m thick |
"""
const SUBSOIL_CLASSES = ["A", "B", "C", "D", "E", "F"]

"""
    download_subsoil_classes(; force = false) -> path of the shapefile

Download the national seismic subsoil class map (BAFU, `ch.bafu.gefahren-baugrundklassen`,
shapefile, LV95, ~80 MB zip) to `data/raw/subsoil_classes/` and unzip it (needs `unzip`). The file
name carries the release date (e.g. `Baugrundklassen_2026_03_04.shp`). `force = true` downloads
the current release again.
"""
function download_subsoil_classes(; force = false)
    dir = datadir("raw", "subsoil_classes"); mkpath(dir)
    zip = joinpath(dir, "gefahren-baugrundklassen_2056.shp.zip")
    if force || !isfile(zip)
        @info "Downloading the seismic subsoil classes (BAFU) from data.geo.admin.ch (~80 MB) ..."
        Downloads.download(SUBSOIL_CLASSES_URL, zip)
        run(`unzip -o -q $zip -d $dir`)
    end
    shp = filter(f -> endswith(f, ".shp"), readdir(dir; join = true))
    isempty(shp) && error("no shapefile in $dir")
    return last(sort(shp))                        # latest release if several were downloaded
end

"""
    subsoil_release(path = download_subsoil_classes()) -> String

Release date of the subsoil class map, from its file name (e.g. "2026-03-04").
"""
function subsoil_release(path = download_subsoil_classes())
    m = match(r"(\d{4})_(\d{2})_(\d{2})", basename(path))
    return isnothing(m) ? "unknown" : join(m.captures, "-")
end

"""
    subsoil_classes_on_grid(x, y; path = download_subsoil_classes()) -> (; class, revision)

Rasterise the seismic subsoil classes onto the regular LV95 grid with cell-centre coordinates `x`,
`y` (km). `class[i, j]` is the code 1–6 of [`SUBSOIL_CLASSES`](@ref) (A–F), 0 where not mapped.
`revision[i, j]` is the SIA 261 edition of the classification: 2003 or 2020 (0 where not mapped).
A cell gets the class of the polygon containing its centre: use a spacing of ≤ 50 m to keep the
small polygons of the 1:25 000 maps.
"""
function subsoil_classes_on_grid(x, y; path = download_subsoil_classes())
    layer = splitext(basename(path))[1]
    cases = join(["WHEN '$(c)' THEN $(i)" for (i, c) in enumerate(SUBSOIL_CLASSES)], " ")
    sql(expr) = ["-sql", "SELECT $(expr) AS code, Geometry FROM \"$(layer)\"", "-dialect", "SQLite", "-a", "code"]
    class = _rasterize(path, nothing, x, y, sql("CASE bgk $(cases) ELSE 0 END"))
    revision = _rasterize(path, nothing, x, y,
                          sql("CASE WHEN Referenz LIKE '%2020%' THEN 2020 WHEN Referenz LIKE '%2003%' THEN 2003 ELSE 0 END"))
    return (; class, revision)
end

const SPECTRAL_MICROZONATION_URL =
    "https://data.geo.admin.ch/ch.bafu.gefahren-spektral/gefahren-spektral/gefahren-spektral_2056.shp.zip"

"""
    download_spectral_microzonation(; force = false) -> path of the shapefile

Perimeters of the spectral seismic microzonation studies (BAFU, `ch.bafu.gefahren-spektral`,
shapefile, LV95, ~0.2 MB) in `data/raw/spectral_microzonation/`. Inside them, the microzonation
(e.g. Sion, Martigny, Brig-Visp, Chablais: PDF maps and design spectra on the cantonal page)
replaces the subsoil classes, so the class map has gaps there.
"""
function download_spectral_microzonation(; force = false)
    dir = datadir("raw", "spectral_microzonation"); mkpath(dir)
    zip = joinpath(dir, "gefahren-spektral_2056.shp.zip")
    if force || !isfile(zip)
        @info "Downloading the spectral microzonation perimeters (BAFU) from data.geo.admin.ch ..."
        Downloads.download(SPECTRAL_MICROZONATION_URL, zip)
        run(`unzip -o -q $zip -d $dir`)
    end
    return last(sort(filter(f -> endswith(f, ".shp"), readdir(dir; join = true))))
end

"""
    spectral_microzonation_on_grid(x, y; path = download_spectral_microzonation()) -> (; zone, names)

Rasterise the spectral microzonation perimeters onto the grid `x`, `y` (km): `zone[i, j]` indexes
`names` (e.g. "Sion VS"), 0 outside every perimeter.
"""
function spectral_microzonation_on_grid(x, y; path = download_spectral_microzonation())
    layer = splitext(basename(path))[1]
    names = Dict{Int,String}()
    AG.read(path) do ds
        for f in AG.getlayer(ds, 0)
            names[Int(AG.getfid(f)) + 1] = something(AG.getfield(f, "spectral_3"), "zone $(AG.getfid(f) + 1)")
        end
    end
    zone = _rasterize(path, nothing, x, y,
                      ["-sql", "SELECT ROWID + 1 AS code, Geometry FROM \"$(layer)\"", "-dialect", "SQLite", "-a", "code"])
    return (; zone, names = [get(names, i, "") for i in 1:maximum(keys(names))])
end

"""
    subsoil_map(r = valais_region(); dx = 0.05)
        -> (; x, y, topo, class, revision, release, ssm, ssm_names)

Seismic subsoil classes and topography on the surface grid of the region (spacing `dx` km), with
the perimeters of the spectral microzonation studies (`ssm`, indexing `ssm_names`; 0 outside).
"""
function subsoil_map(r::ValaisRegion = valais_region(); dx = 0.05)
    path = download_subsoil_classes()
    x, y, topo = fetch_topography(r; dx = dx)
    s = subsoil_classes_on_grid(x, y; path = path)
    m = spectral_microzonation_on_grid(x, y)
    return (; x, y, topo, s.class, s.revision, release = subsoil_release(path), ssm = m.zone, ssm_names = m.names)
end

"""
    write_subsoil_paraview(sm, outdir = datadir("output"); lift = 0.01)

Write the subsoil class map of [`subsoil_map`](@ref) as `valais_subsoil.vts`: the topographic surface
(raised by `lift` km to draw over `valais_surface.vts`) with fields `SubsoilClass` (1–6 = A–F;
`NaN` where not mapped, transparent with NaN Opacity = 0), `SIA261_revision` and
`SpectralMicrozonation` (study perimeter index, `NaN` outside). The legend is written to
`subsoil_legend.csv`.
"""
function write_subsoil_paraview(sm, outdir = datadir("output"); lift = 0.01)
    mkpath(outdir)
    X, Y, Z = xyz_grid(sm.x, sm.y, 0.0)
    Z[:, :, 1] .= replace(sm.topo, NaN => 0.0) .+ lift
    cls = Float64.(sm.class); cls[cls .== 0] .= NaN
    rev = Float64.(sm.revision); rev[rev .== 0] .= NaN
    ssm = Float64.(sm.ssm); ssm[ssm .== 0] .= NaN
    d = CartData(X, Y, Z, (SubsoilClass = reshape(cls, size(X)), SIA261_revision = reshape(rev, size(X)),
                           SpectralMicrozonation = reshape(ssm, size(X))))
    write_paraview(d, joinpath(outdir, "valais_subsoil"))
    open(joinpath(outdir, "subsoil_legend.csv"), "w") do io
        println(io, "code,class,cells,release")
        for (i, c) in enumerate(SUBSOIL_CLASSES)
            println(io, i, ",", c, ",", count(==(i), sm.class), ",", sm.release)
        end
        println(io, "# SpectralMicrozonation field: index,zone,cells")
        for (i, n) in enumerate(sm.ssm_names)
            count(==(i), sm.ssm) > 0 && println(io, "# ", i, ",\"", n, "\",", count(==(i), sm.ssm))
        end
    end
    return outdir
end
