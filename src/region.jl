"""
    ValaisRegion(; E=(2550.0, 2680.0), N=(1075.0, 1170.0), z=(-70.0, 5.0))

Model domain in Swiss LV95 coordinates (EPSG:2056), in km:
`E` easting range, `N` northing range, `z` elevation range (negative below sea level).
The default box covers the canton of Valais/Wallis with a small margin, down to 70 km
(below the deepest European Moho in the area, ~55 km).
"""
Base.@kwdef struct ValaisRegion
    E::NTuple{2,Float64} = (2550.0, 2680.0)
    N::NTuple{2,Float64} = (1075.0, 1170.0)
    z::NTuple{2,Float64} = (-70.0, 5.0)
end

"""
    valais_region(; kwargs...) -> ValaisRegion

Convenience constructor for the default Valais/Wallis model domain.
"""
valais_region(; kwargs...) = ValaisRegion(; kwargs...)

function _transform(xs, ys, from, to)
    @assert size(xs) == size(ys)
    outx = similar(xs, Float64)
    outy = similar(ys, Float64)
    AG.createcoordtrans(from, to) do ct
        for i in eachindex(xs)
            g = AG.createpoint(Float64(xs[i]), Float64(ys[i]))
            AG.transform!(g, ct)
            outx[i] = AG.getx(g, 0)
            outy[i] = AG.gety(g, 0)
        end
    end
    return outx, outy
end

"""
    lv95_to_lonlat(E, N) -> (lon, lat)

Convert LV95 coordinates **in km** (scalars or arrays) to WGS84 longitude/latitude.
"""
function lv95_to_lonlat(E::AbstractArray, N::AbstractArray)
    _transform(E .* 1e3, N .* 1e3, AG.importEPSG(2056), AG.importEPSG(4326; order = :trad))
end
lv95_to_lonlat(E::Real, N::Real) = first.(lv95_to_lonlat([E], [N]))

"""
    lonlat_to_lv95(lon, lat) -> (E, N)

Convert WGS84 longitude/latitude (scalars or arrays) to LV95 coordinates **in km**.
"""
function lonlat_to_lv95(lon::AbstractArray, lat::AbstractArray)
    E, N = _transform(lon, lat, AG.importEPSG(4326; order = :trad), AG.importEPSG(2056))
    return E ./ 1e3, N ./ 1e3
end
lonlat_to_lv95(lon::Real, lat::Real) = first.(lonlat_to_lv95([lon], [lat]))
