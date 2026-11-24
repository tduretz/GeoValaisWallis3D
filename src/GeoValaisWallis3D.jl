"""
    GeoValaisWallis3D

Tools to build a 3D geological/geophysical model of the Valais/Wallis region (Swiss Alps)
from the surface down to the Moho, based on GeophysicalModelGenerator.jl.

Data sources
- Topography: swisstopo swissALTIRegio (10 m DEM, LV95), https://map.geo.admin.ch
- Surface tectonics: swisstopo Tectonic Map of Switzerland 1:500 000
- Quaternary deposits: swisstopo swissBEDROCK (bedrock elevation / thickness of unconsolidated deposits)
- Crustal layers & densities: 3D-ALPS (Spooner et al. 2019); temperature: 3D-ALPS-TR (Spooner et al. 2020)
- Moho: Spada et al. (2013, GJI) Western Alps Moho maps (European + Adriatic Moho)
- Nappe de Morcles 3D geological model (Thornton et al. 2018, Sci. Data)

The model is built in Swiss LV95 coordinates (EPSG:2056), expressed in km.
"""
module GeoValaisWallis3D

using Reexport
@reexport using GeophysicalModelGenerator
import ArchGDAL as AG
using DelimitedFiles
using Downloads

export ValaisRegion, valais_region, datadir
export lv95_to_lonlat, lonlat_to_lv95
export fetch_topography, fetch_swissbedrock
export download_tectonic_map, tectonic_map_on_grid, TectonicMap
export download_moho_spada, read_moho_spada, interpolate_scattered
export download_3dalps, download_3dalps_tr, alps3d_on_grid, temperature_3dalps
export read_morcles, MorclesModel, MORCLES_UNITS, morcles_unit_at
export refined_levels, fine_band_from_bedrock
export LayerVelocity, DEFAULT_VELOCITIES
export build_valais_model, model_interfaces, PHASES, cartdata_to_geodata

"""
    datadir(parts...)

Path inside the data directory (downloads go to `data/raw`, outputs to `data/output`). This is the
package's `data/` folder, or `\$GEOVALAIS_DATADIR` if that environment variable is set (e.g. a scratch
file system on a cluster).
"""
datadir(parts...) = joinpath(get(ENV, "GEOVALAIS_DATADIR", joinpath(dirname(@__DIR__), "data")), parts...)

include("region.jl")
include("topography.jl")
include("bedrock.jl")
include("tectonics.jl")
include("moho.jl")
include("alps3d.jl")
include("morcles.jl")
include("grid.jl")
include("velocity.jl")
include("model.jl")

end # module GeoValaisWallis3D
