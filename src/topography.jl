const SWISSALTIREGIO_URL =
    "https://data.geo.admin.ch/ch.swisstopo.swissaltiregio/swissaltiregio/swissaltiregio_2056_5728.tif"

_region_tag(r, dx) = "$(round(Int, r.E[1]))-$(round(Int, r.E[2]))_$(round(Int, r.N[1]))-$(round(Int, r.N[2]))_$(round(Int, dx * 1000))m"

# GDAL options to cut the region `r` out of an LV95 raster at spacing `dx` (km).
# Every raster of the package uses these, so all grids are aligned
# (cell centres at E[1] + dx/2, E[1] + 3dx/2, ...).
_window_opts(r, dx) = string.(["-projwin", r.E[1] * 1e3, r.N[2] * 1e3, r.E[2] * 1e3, r.N[1] * 1e3,
                               "-tr", dx * 1e3, dx * 1e3, "-r", "average",
                               "-of", "GTiff", "-co", "COMPRESS=DEFLATE"])

"""
    fetch_topography(r::ValaisRegion = valais_region(); dx = 0.25, cache = true)
        -> (x, y, z)

Surface topography from swisstopo **swissALTIRegio** (10 m DEM, LV95, cloud-optimised GeoTIFF).

Only the window covering `r` is read over HTTP (GDAL `/vsicurl/`), resampled
(block average) to a grid spacing `dx` in km. The result is cached as a GeoTIFF in
`data/raw/`. Returns cell-centre coordinates `x` (easting, km), `y` (northing, km)
and elevation `z[i,j]` in km (NaN where no data).
"""
function fetch_topography(r::ValaisRegion = valais_region(); dx = 0.25, cache = true)
    mkpath(datadir("raw"))
    file = datadir("raw", "swissaltiregio_$(_region_tag(r, dx)).tif")
    if !(cache && isfile(file))
        @info "Reading swissALTIRegio window from data.geo.admin.ch (dx = $(dx) km) ..."
        AG.read("/vsicurl/" * SWISSALTIREGIO_URL) do ds
            AG.gdaltranslate(ds, _window_opts(r, dx); dest = file) do _ end
        end
        @info "Saved $(file)"
    end
    x, y, z = _read_lv95_raster(file)
    return x, y, z ./ 1e3
end

# Read band `band` of an LV95 GeoTIFF -> (x [km], y [km], values) with y increasing.
# Values are returned unscaled, with nodata replaced by NaN.
function _read_lv95_raster(file; band = 1)
    AG.read(file) do ds
        gt = AG.getgeotransform(ds)
        b = AG.getband(ds, band)
        z = Float64.(AG.read(b))           # (ncols, nrows), rows from north to south
        nd = AG.getnodatavalue(b)
        isnothing(nd) || (z[z .== nd] .= NaN)
        nx, ny = size(z)
        x = [(gt[1] + (i - 0.5) * gt[2]) / 1e3 for i in 1:nx]
        y = [(gt[4] + (j - 0.5) * gt[6]) / 1e3 for j in 1:ny]
        return x, reverse(y), reverse(z, dims = 2)
    end
end
