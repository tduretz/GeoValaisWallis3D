const TECTONIC_MAP_URL =
    "https://data.geo.admin.ch/ch.swisstopo.geologie-tektonische_karte/geologie-tektonische_karte/geologie-tektonische_karte_2056_en.gpkg"

"""
    download_tectonic_map(; force = false) -> path

Download the swisstopo *Tectonic Map of Switzerland 1:500 000* (GeoPackage, LV95, English)
to `data/raw/` (~24 MB).
"""
function download_tectonic_map(; force = false)
    mkpath(datadir("raw"))
    file = datadir("raw", "tektonische_karte_2056_en.gpkg")
    if force || !isfile(file)
        @info "Downloading swisstopo Tectonic Map 1:500k ..."
        Downloads.download(TECTONIC_MAP_URL, file)
    end
    return file
end

"""
    TectonicMap

Surface tectonic map rasterised on a grid. Integer codes index into the name tables
(`0` = no unit / outside map):
- `level1` : main domain (Helvetic, Penninic, Lepontic, South Alpine, ...)
- `level2` : subdomain (Lower/Middle/Upper Penninic, Lower/Upper Helvetic, ...)
- `unit`   : most specific unit (nappe / massif / zone, e.g. Siviez-Mischabel Nappe)
- `quaternary` : `true` where the map shows Quaternary cover
"""
struct TectonicMap
    level1::Matrix{Int}
    level2::Matrix{Int}
    unit::Matrix{Int}
    quaternary::Matrix{Bool}
    names1::Vector{String}
    names2::Vector{String}
    names_unit::Vector{String}
end

# Rasterise a layer of `path` on the cell-centred grid (x, y) [km]; returns Int32 matrix (x, y).
function _rasterize(path, layer, x, y, burnopts)
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    te = [x[1] - dx / 2, y[1] - dy / 2, x[end] + dx / 2, y[end] + dy / 2] .* 1e3
    # layer = nothing: the features are selected by an -sql option in `burnopts`
    opts = string.([(isnothing(layer) ? [] : ["-l", layer])..., burnopts..., "-te", te..., "-ts", length(x), length(y),
                    "-ot", "Int32", "-init", 0, "-of", "GTiff"])
    AG.read(path) do ds
        AG.gdalrasterize(ds, opts; dest = "/vsimem/raster_$(something(layer, "sql")).tif") do r
            a = Int.(AG.read(AG.getband(r, 1)))   # (nx, ny), north -> south
            return reverse(a, dims = 2)
        end
    end
end

"""
    tectonic_map_on_grid(x, y; path = download_tectonic_map()) -> TectonicMap

Rasterise the swisstopo tectonic units onto the regular LV95 grid with cell-centre
coordinates `x`, `y` (km).
"""
function tectonic_map_on_grid(x, y; path = download_tectonic_map())
    oid = _rasterize(path, "tecto_units_en", x, y, ["-a", "OBJECTID"])
    quat = _rasterize(path, "quat_surfaces_en", x, y, ["-burn", 1]) .> 0

    # attribute table: OBJECTID -> (Tecto_1, Tecto_2, most specific unit)
    attrs = Dict{Int,NTuple{3,String}}()
    AG.read(path) do ds
        for f in AG.getlayer(ds, "tecto_units_en")
            t = [something(AG.getfield(f, k), "-") for k in ("Tecto_1_EN", "Tecto_2_EN", "Tecto_3_EN", "Tecto_4_EN")]
            specific = something(findlast(!=("-"), t), 1)
            attrs[AG.getfield(f, "OBJECTID")] = (t[1], t[2] == "-" ? t[1] : t[2], t[specific])
        end
    end
    names = [sort(unique(getindex.(values(attrs), k))) for k in 1:3]
    lut = [Dict(n => i for (i, n) in enumerate(names[k])) for k in 1:3]
    codes = [zeros(Int, size(oid)) for _ in 1:3]
    for I in eachindex(oid)
        a = get(attrs, oid[I], nothing)
        isnothing(a) && continue
        for k in 1:3
            codes[k][I] = lut[k][a[k]]
        end
    end
    return TectonicMap(codes..., quat, names...)
end
