# 3D-ALPS (Spooner et al. 2019) structural/density model and 3D-ALPS-TR (Spooner et al. 2020)
# thermal model of the Alpine lithosphere. GFZ Data Services, CC BY 4.0.
# Both are given on a 20 km x 20 km grid in UTM zone 32N (EPSG:32632).

const ALPS3D_BASE = "https://datapub.gfz.de/download/10.5880.GFZ.4.5.2019.004"
const ALPS3D_TR_BASE = "https://datapub.gfz.de/download/10.5880.GFZ.4.5.2020.007nbif"
const ALPS3D_LAYERS = ("Unconsolidated Sediments", "Consolidated Sediments", "Upper Crust",
                       "Lower Crust", "Lithospheric Mantle")

"""
    download_3dalps(; force = false) -> Dict(layer name => path)

Download the five layer files of the **3D-ALPS** gravity-constrained density model
(Spooner et al. 2019, doi:10.5880/GFZ.4.5.2019.004, CC BY 4.0) to `data/raw/3dalps/`.
Each file holds, per 20 km node (UTM 32N): top (m a.s.l.), thickness (m), density (kg/m³).
"""
function download_3dalps(; force = false)
    dir = datadir("raw", "3dalps")
    mkpath(dir)
    Dict(map(ALPS3D_LAYERS) do l
        file = joinpath(dir, replace(l, " " => "_") * ".txt")
        if force || !isfile(file)
            @info "Downloading 3D-ALPS layer $(l) ..."
            Downloads.download("$(ALPS3D_BASE)/2019-004_Spooner_$(replace(l, " " => "%20")).txt", file)
        end
        l => file
    end)
end

"""
    download_3dalps_tr(; force = false) -> path

Download and unzip the **3D-ALPS-TR** thermal/rheological model (Spooner et al. 2020,
doi:10.5880/GFZ.4.5.2020.007, CC BY 4.0) to `data/raw/3dalps/`. Requires the `unzip` command.
"""
function download_3dalps_tr(; force = false)
    dir = datadir("raw", "3dalps")
    mkpath(dir)
    txt = joinpath(dir, "3D-ALPS-TR_Final.txt")
    if force || !isfile(txt)
        zip = joinpath(dir, "2020-007_Spooner-et-al_3D-ALPS-TR_data.zip")
        isfile(zip) || (@info "Downloading 3D-ALPS-TR ..."; Downloads.download("$(ALPS3D_TR_BASE)/$(basename(zip))", zip))
        run(`unzip -o -q $zip 3D-ALPS-TR_Final.txt -d $dir`)
    end
    return txt
end

# Regular 2D grid (UTM, m) with node values
struct Grid2D
    x::Vector{Float64}
    y::Vector{Float64}
    v::Matrix{Float64}
end

function _grid2d(xs, ys, vs)
    xu = sort(unique(xs)); yu = sort(unique(ys))
    v = fill(NaN, length(xu), length(yu))
    for k in eachindex(xs)
        v[searchsortedfirst(xu, xs[k]), searchsortedfirst(yu, ys[k])] = vs[k]
    end
    return Grid2D(xu, yu, v)
end

# Bilinear interpolation on a regular grid (clamped to the grid boundary).
function _bilinear_weights(g, xq, yq)
    fx = clamp((xq - g.x[1]) / (g.x[2] - g.x[1]), 0, length(g.x) - 1 - 1e-9)
    fy = clamp((yq - g.y[1]) / (g.y[2] - g.y[1]), 0, length(g.y) - 1 - 1e-9)
    i = floor(Int, fx) + 1; j = floor(Int, fy) + 1
    a = fx - (i - 1); b = fy - (j - 1)
    return i, j, a, b
end
function _bilinear(g::Grid2D, xq, yq)
    i, j, a, b = _bilinear_weights(g, xq, yq)
    (1 - a) * (1 - b) * g.v[i, j] + a * (1 - b) * g.v[i+1, j] + (1 - a) * b * g.v[i, j+1] + a * b * g.v[i+1, j+1]
end

# LV95 grid (km) -> UTM 32N coordinates (m), as matrices (nx, ny)
function _lv95grid_to_utm(x, y)
    E = [xi for xi in x, _ in y]; N = [yj for _ in x, yj in y]
    _transform(E .* 1e3, N .* 1e3, AG.importEPSG(2056), AG.importEPSG(32632))
end

"""
    alps3d_on_grid(x, y) -> NamedTuple

Bilinear interpolation of the 3D-ALPS layers onto the LV95 grid `x`, `y` (km). Returns matrices
`top_<layer>` (elevation, km), `thick_<layer>` (km) and `rho_<layer>` (kg/m³) for
`us` (unconsolidated sed.), `cs` (consolidated sed.), `uc` (upper crust), `lc` (lower crust),
`lm` (lithospheric mantle).
"""
function alps3d_on_grid(x, y)
    files = download_3dalps()
    ux, uy = _lv95grid_to_utm(x, y)
    out = Dict{Symbol,Matrix{Float64}}()
    for (l, s) in zip(ALPS3D_LAYERS, (:us, :cs, :uc, :lc, :lm))
        d = readdlm(files[l], '\t', Float64; skipstart = 4)
        for (col, name, scale) in ((3, :top, 1e-3), (4, :thick, 1e-3), (5, :rho, 1.0))
            g = _grid2d(d[:, 1], d[:, 2], d[:, col] .* scale)
            out[Symbol(name, :_, s)] = [_bilinear(g, ux[I], uy[I]) for I in CartesianIndices(ux)]
        end
    end
    return NamedTuple(out)
end

"""
    temperature_3dalps(x, y, z) -> Array{Float64,3}

Temperature (°C) of the 3D-ALPS-TR steady-state thermal model interpolated onto the LV95 grid
`x`, `y`, `z` (km): bilinear between the 20 km model columns, linear in depth along each column
(constant extrapolation above the column top / below its base). `z` is a vector of levels, or the heights of a
curvilinear grid, an array of size `(length(x), length(y), nz)`.
"""
function temperature_3dalps(x, y, z)
    d = readdlm(download_3dalps_tr(), '\t', Float64; comments = true, comment_char = '#', skipstart = 4)
    xu = sort(unique(d[:, 1])); yu = sort(unique(d[:, 2]))
    zc = [Float64[] for _ in xu, _ in yu]; Tc = [Float64[] for _ in xu, _ in yu]
    for k in axes(d, 1)
        i = searchsortedfirst(xu, d[k, 1]); j = searchsortedfirst(yu, d[k, 2])
        push!(zc[i, j], d[k, 3] / 1e3); push!(Tc[i, j], d[k, 4])
    end
    for I in eachindex(zc)               # sort each column by increasing z
        p = sortperm(zc[I]); zc[I] = zc[I][p]; Tc[I] = Tc[I][p]
    end
    g = Grid2D(xu, yu, zeros(length(xu), length(yu)))
    ux, uy = _lv95grid_to_utm(x, y)
    nz = size(z, ndims(z))
    T = Array{Float64}(undef, length(x), length(y), nz)
    colinterp(zz, TT, zq) = begin
        isempty(zz) && return NaN
        k = searchsortedfirst(zz, zq)
        k == 1 && return TT[1]
        k > length(zz) && return TT[end]
        w = (zq - zz[k-1]) / (zz[k] - zz[k-1]); (1 - w) * TT[k-1] + w * TT[k]
    end
    Threads.@threads for J in CartesianIndices(ux)
        i, j, a, b = _bilinear_weights(g, ux[J], uy[J])
        w = ((i, j, (1 - a) * (1 - b)), (i + 1, j, a * (1 - b)), (i, j + 1, (1 - a) * b), (i + 1, j + 1, a * b))
        for k in 1:nz
            zq = z isa AbstractVector ? z[k] : z[J[1], J[2], k]
            T[J[1], J[2], k] = sum(ww * colinterp(zc[ii, jj], Tc[ii, jj], zq) for (ii, jj, ww) in w)
        end
    end
    return T
end
