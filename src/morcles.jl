# 3D geological model of the Nappe de Morcles (Thornton, Mariethoz & Brunner 2018, Sci. Data 5:180238,
# doi:10.1038/sdata.2018.238; data doi:10.6084/m9.figshare.c.4130759, CC0).

const MORCLES_URL = "https://ndownloader.figshare.com/files/11970980"

"""
Morcles formations in stratigraphic order (codes 1…18), as named in the dataset
(basement = Aiguilles Rouges crystalline basement; `a`, `i…` Triassic–Jurassic,
`c…` Cretaceous, `e…` Eocene–Oligocene; see Thornton et al. 2018, Table 1).
"""
const MORCLES_UNITS = ["basement", "a", "i1i", "i1s", "i2", "i3_4", "i5", "i6_8", "i8c2", "c2k", "c3",
                       "c4", "c4_5", "c6_8", "e5", "e6c", "e6k", "e6o1"]

"""
    MorclesModel

Voxel model (50 m) of the Nappe de Morcles converted to LV95 (km). `unit[i,j,k]` is the index
into [`MORCLES_UNITS`](@ref) (0 = outside / above topography). `x`, `y`, `z` are cell centres.
"""
struct MorclesModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    unit::Array{Int8,3}
end

"""
    read_morcles(; resolution = 50) -> MorclesModel

Download (418 MB zip, once) and read the generic voxel export of the Morcles model.
LV03 coordinates are shifted to LV95 (+2 000 km / +1 000 km; the residual LV03/LV95
distortion of < 2 m is negligible at 50 m resolution).
"""
function read_morcles(; resolution = 50)
    dir = datadir("raw", "morcles")
    mkpath(dir)
    vox = joinpath(dir, "$(resolution)m.vox")
    if !isfile(vox)
        zip = joinpath(dir, "Thornton2018_Morcles_data.zip")
        isfile(zip) || (@info "Downloading Morcles 3D model (418 MB) ..."; Downloads.download(MORCLES_URL, zip))
        run(`unzip -o -q -j $zip "*/Generic voxel format/$(resolution)m.vox" -d $dir`)
    end
    lut = Dict(n => Int8(i) for (i, n) in enumerate(MORCLES_UNITS))
    open(vox) do io
        hdr = Dict{String,String}()
        for _ in 1:10
            k, v = split(strip(readline(io)))
            hdr[k] = v
        end
        nx, ny, nz = parse.(Int, (hdr["nx"], hdr["ny"], hdr["nz"]))
        x0, y0, z0, dx, dy, dz = parse.(Float64, (hdr["x0"], hdr["y0"], hdr["z0"], hdr["dx"], hdr["dy"], hdr["dz"]))
        unit = zeros(Int8, nx, ny, nz)
        for k in 1:nz, j in 1:ny, i in 1:nx          # x fastest, z slowest
            unit[i, j, k] = get(lut, strip(readline(io)), Int8(0))
        end
        x = (x0 .+ dx .* (0:nx-1) .+ 2.0e6) ./ 1e3
        y = (y0 .+ dy .* (0:ny-1) .+ 1.0e6) ./ 1e3
        z = (z0 .+ dz .* (0:nz-1)) ./ 1e3
        return MorclesModel(collect(x), collect(y), collect(z), unit)
    end
end

"""
    morcles_unit_at(m::MorclesModel, x, y, z) -> Int8

Nearest-voxel lookup of the Morcles unit at LV95 point (`x`, `y`, `z` in km); 0 outside.
"""
function morcles_unit_at(m::MorclesModel, x, y, z)
    h = m.x[2] - m.x[1]
    i = round(Int, (x - m.x[1]) / h) + 1
    j = round(Int, (y - m.y[1]) / h) + 1
    k = round(Int, (z - m.z[1]) / (m.z[2] - m.z[1])) + 1
    (1 <= i <= length(m.x) && 1 <= j <= length(m.y) && 1 <= k <= length(m.z)) || return Int8(0)
    return m.unit[i, j, k]
end
