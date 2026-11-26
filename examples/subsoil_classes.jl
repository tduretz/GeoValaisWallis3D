# Seismic subsoil classes (SIA 261) of the Valais, to overlay on the 3D model in ParaView.
#     julia --project=. examples/subsoil_classes.jl
#
# Data: BAFU national map of seismic subsoil classes (ch.bafu.gefahren-baugrundklassen, ~80 MB, cached
# in data/raw/subsoil_classes). For Valais it contains the cantonal 1:25 000 maps of the Rhône valley
# (CREALP; https://www.vs.ch/web/sdm/cartes-de-sol-de-fondation).
# Output in data/output/: valais_subsoil.vts (surface at 50 m, 10 m above valais_surface.vts),
# subsoil_legend.csv. Field SubsoilClass: 1-6 = A-F, NaN where not mapped; SpectralMicrozonation: perimeters
# of the microzonation studies (Sion, Martigny, Brig-Visp, Chablais, ...) where the classes are replaced.
using GeoValaisWallis3D

sm = subsoil_map(valais_region(); dx = 0.05)
println("subsoil classes, release ", sm.release)
for (i, c) in enumerate(SUBSOIL_CLASSES)
    println("  ", c, ": ", round(count(==(i), sm.class) * 0.05^2, digits = 1), " km²")
end
println("spectral microzonation studies in the region (the class map has gaps there):")
for (i, n) in enumerate(sm.ssm_names)
    k = count(==(i), sm.ssm); k > 0 && println("  ", n, ": ", round(k * 0.05^2, digits = 1), " km²")
end
write_subsoil_paraview(sm)
