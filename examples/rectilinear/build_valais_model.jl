# Build the 3D Valais/Wallis geological & seismic velocity model and write ParaView files.
# Run from the repository root with:
#     julia --project=. -t auto examples/rectilinear/build_valais_model.jl
#
# Data are fetched on first run and cached in data/raw/ (first run: ~5 GB transferred, ~15 min):
#   - swisstopo swissALTIRegio DEM (Valais window read over HTTP)
#   - swisstopo swissBEDROCK (bedrock elevation / Quaternary thickness, 10 m tiles)
#   - swisstopo Tectonic Map of Switzerland 1:500k (GeoPackage)
#   - 3D-ALPS (Spooner et al. 2019) layers + densities; 3D-ALPS-TR temperature (Spooner et al. 2020)
#   - Spada et al. (2013) Western Alps Moho maps
#   - Nappe de Morcles 3D model (Thornton et al. 2018), 50 m voxels
# Outputs go to data/output/ (open the .vts files in ParaView).
using GeoValaisWallis3D

r = valais_region()        # LV95 km: E 2550-2680, N 1075-1170, z -70..5

res = build_valais_model(r;
    dx          = 0.25,    # horizontal spacing [km]
    dx_surface  = 0.1,     # spacing of the surface (map) dataset [km]
    dz_fine     = 0.025,   # vertical spacing in the Quaternary band [km]
    fine_band   = :auto,   # elevation band of fine dz, from swissBEDROCK (or e.g. (-0.6, 1.6))
    dz_max_up   = 0.25,    # max dz above the band (mountains, air) [km]
    dz_max_down = 2.0,     # max dz at depth [km]
)

println(res.model)
println("z levels: ", length(res.z), ", dz min/max = ", extrema(diff(res.z)), " km")
