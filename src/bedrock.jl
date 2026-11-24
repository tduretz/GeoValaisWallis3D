const SWISSBEDROCK_URL(E, N) =
    "https://data.geo.admin.ch/ch.swisstopo.geologie-swissbedrock/geologie-swissbedrock_$(E)-$(N)/geologie-swissbedrock_$(E)-$(N)_10_2056.tif"

"""
    fetch_swissbedrock(r::ValaisRegion = valais_region(); dx = 0.05, cache = true)
        -> (x, y, bem, tmud)

swisstopo **swissBEDROCK** (Release 1, Aug. 2025; 10 m, LV95) on the region `r` at spacing `dx` (km):
- `bem`  : bedrock elevation model (km a.s.l.)
- `tmud` : thickness of unconsolidated (Quaternary) deposits (km)

swissBEDROCK only contains values where unconsolidated deposits exist (valleys, basins); bare
bedrock is "no data". Here no-data TMUD is therefore treated as **0** (inside a tile) *before*
averaging, so that coarse cells give the correct mean thickness; `bem` is averaged over valid
pixels only (NaN where no deposits). Tiles that do not exist (outside Switzerland) give NaN.

swissBEDROCK comes as 10 km × 10 km tiles of 11 pixel-interleaved bands (~46 MB each). Each tile
is read once, reduced to bands 1–2 block-averaged at `tile_dx` (default 50 m), and cached in
`data/raw/swissbedrock/` (the first call for the whole of Valais transfers ~5 GB and takes
~10 min). Same grid as [`fetch_topography`](@ref) for equal `r` and `dx`.
"""
function fetch_swissbedrock(r::ValaisRegion = valais_region(); dx = 0.05, tile_dx = 0.05, cache = true)
    tiledir = datadir("raw", "swissbedrock")
    mkpath(tiledir)
    file = datadir("raw", "swissbedrock_$(_region_tag(r, dx)).tif")
    if !(cache && isfile(file))
        tiles = [(E, N) for E in 10 * fld(Int(floor(r.E[1])), 10):10:Int(ceil(r.E[2])) - 1
                        for N in 10 * fld(Int(floor(r.N[1])), 10):10:Int(ceil(r.N[2])) - 1]
        @info "swissBEDROCK: $(length(tiles)) candidate tiles (cached in $(tiledir))"
        done = Threads.Atomic{Int}(0)
        files = Vector{Union{Nothing,String}}(nothing, length(tiles))
        Threads.@threads for k in eachindex(tiles)
            E, N = tiles[k]
            f = joinpath(tiledir, "swissbedrock_v2_$(E)-$(N)_$(round(Int, tile_dx * 1000))m.tif")
            missingflag = f * ".missing"
            if isfile(f)
                files[k] = f
            elseif !isfile(missingflag)
                url = SWISSBEDROCK_URL(E, N)
                resp = Downloads.request(url; method = "HEAD", throw = false)
                if resp isa Downloads.Response && resp.status == 404
                    touch(missingflag)       # tile does not exist (outside Switzerland)
                else
                    tmp = f * ".part"
                    AG.read("/vsicurl/" * url) do ds
                        _reduce_bedrock_tile(ds, tmp, tile_dx)
                    end
                    mv(tmp, f; force = true)
                    files[k] = f
                end
            end
            n = Threads.atomic_add!(done, 1) + 1
            n % 10 == 0 && @info "  swissBEDROCK tiles: $(n)/$(length(tiles))"
        end
        avail = filter(!isnothing, files)
        @info "swissBEDROCK: mosaicking $(length(avail)) tiles"
        dss = [AG.read(f) for f in avail]
        try
            AG.gdalbuildvrt(dss, ["-srcnodata", "-9999", "-vrtnodata", "-9999"]) do vrt
                AG.gdaltranslate(vrt, _window_opts(r, dx); dest = file) do _ end
            end
        finally
            foreach(AG.destroy, dss)
        end
    end
    x, y, bem = _read_lv95_raster(file; band = 1)
    _, _, tmud = _read_lv95_raster(file; band = 2)
    return x, y, bem ./ 1e3, tmud ./ 1e3
end

# Block-average bands 1 (BEM, NaN-skipping) and 2 (TMUD, no data -> 0) of a 10 m swissBEDROCK
# tile to `tile_dx` (km) and write a 2-band Float32 GeoTIFF (nodata = -9999).
function _reduce_bedrock_tile(ds, dest, tile_dx)
    gt = AG.getgeotransform(ds)
    f = round(Int, tile_dx * 1e3 / gt[2])
    nd1 = AG.getnodatavalue(AG.getband(ds, 1)); nd2 = AG.getnodatavalue(AG.getband(ds, 2))
    bem = Float64.(AG.read(AG.getband(ds, 1))); tmud = Float64.(AG.read(AG.getband(ds, 2)))
    valid = bem .!= nd1
    tmud[tmud .== nd2] .= 0.0
    nx, ny = size(bem) .÷ f
    B = fill(-9999.0f0, nx, ny); T = fill(-9999.0f0, nx, ny)
    for j in 1:ny, i in 1:nx
        I = (i-1)*f+1:i*f; J = (j-1)*f+1:j*f
        v = @view valid[I, J]
        nv = count(v)
        nv > 0 && (B[i, j] = Float32(sum(bem[I, J][v]) / nv))
        T[i, j] = Float32(sum(@view tmud[I, J]) / f^2)
    end
    AG.create(dest; driver = AG.getdriver("GTiff"), width = nx, height = ny, nbands = 2,
              dtype = Float32, options = ["COMPRESS=DEFLATE"]) do out
        AG.setgeotransform!(out, [gt[1], gt[2] * f, 0.0, gt[4], 0.0, gt[6] * f])
        AG.setproj!(out, AG.getproj(ds))
        for (b, a) in ((1, B), (2, T))
            band = AG.getband(out, b)
            AG.setnodatavalue!(band, -9999.0)
            AG.write!(band, a)
        end
    end
end
