# Build the Valais/Wallis model on a curvilinear mesh (GeophysicalModelGenerator curvilinear_grid) and write ParaView files.
#     julia --project=. -t auto examples/curvilinear/build_valais_model.jl
# The grid levels follow the topography, the base of the consolidated sediments, the Conrad and the European Moho (see
# meshing_surfaces); the Quaternary fill lies in the top `fine_depth` km, meshed with dz[1]. Same data as
# examples/rectilinear/build_valais_model.jl (cached in data/raw/). Output: data/output/valais_curvilinear_model.vts (≈ 57 M nodes,
# several GB: use a larger dx for a quick look) and valais_curvilinear_surfaces.vts.
using GeoValaisWallis3D

res = build_valais_curvilinear_model(valais_region();
    dx         = 0.25,                          # horizontal spacing [km]
    dz         = (0.025, 0.1, 0.5, 0.5, 1.0),   # spacing where each layer is thickest [km]: near-surface, sediments, upper crust, lower crust, mantle
    fine_depth = 1.0,                           # thickness of the near-surface layer [km]
    spacing    = :uniform,                      # or :sbp (nodes of a summation-by-parts grid, SeismicWavePC.jl)
)
println(res.model)
println("levels of the meshed surfaces: ", res.maps.levels, "; cells per layer: ", res.maps.nz)
