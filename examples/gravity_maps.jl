# Bouguer and free-air anomaly maps of the Valais, to overlay on the 3D model in ParaView.
#     julia --project=. -t auto examples/gravity_maps.jl
#
# Data (swisstopo, api3.geo.admin.ch; ~10 s download, cached in data/raw/):
#   - Gravimetric Atlas of Switzerland 1:100 000: measurement points with complete Bouguer anomaly
#   - Gravimetric Map of Switzerland 1:500 000 (Bouguer anomalies): contour lines
#   - National Gravity Network LSN2004: absolute gravity at 177 stations in the region
#   - swissALTIRegio DEM, for the attraction of the topography (out to 100 km around the region)
# swisstopo publishes no free-air map: it is computed as Bouguer + attraction of the topography
# (2670 kg/m³), and checked against the exact free-air anomaly g − γ at the LSN2004 stations.
# Output in data/output/: valais_gravity.vts (surface, 20 m above valais_surface.vts),
# valais_gravity_atlas.vtu and valais_gravity_network.vtu (stations). ~20 s on 12 threads.
using GeoValaisWallis3D, Statistics

gm = gravity_maps(valais_region(); dx = 0.25)
n = gm.network
println("free-air check at $(length(n.E)) LSN2004 stations: map − exact = ",
        round(mean(n.misfit), digits = 1), " ± ", round(std(n.misfit), digits = 1), " mGal")
write_gravity_paraview(gm)
