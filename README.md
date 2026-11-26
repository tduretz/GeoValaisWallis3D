# GeoValaisWallis3D

3D geological / seismic velocity model of the Valais/Wallis region (Swiss Alps), from the
topographic surface down to below the Moho, built with
[GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

## Data sources

| Component | Dataset | Source / licence |
|---|---|---|
| Topography | swissALTIRegio (10 m DEM) — Valais window read over HTTP | swisstopo, [data.geo.admin.ch](https://data.geo.admin.ch/browser/index.html#/collections/ch.swisstopo.swissaltiregio), OGD |
| Quaternary deposits | swissBEDROCK3D (release 15 Oct 2026): bedrock elevation, thickness of unconsolidated deposits (10 m) | swisstopo, [swissBEDROCK](https://www.swisstopo.admin.ch/en/3d-swissbedrock-en), OGD |
| Surface geology | Tectonic Map of Switzerland 1:500 000 | swisstopo, `ch.swisstopo.geologie-tektonische_karte`, OGD |
| Crustal layers, density | 3D-ALPS (20 km grid): sediments, upper/lower crust, lithospheric mantle | Spooner et al. (2019), [doi:10.5880/GFZ.4.5.2019.004](https://doi.org/10.5880/GFZ.4.5.2019.004), CC BY 4.0 |
| Temperature | 3D-ALPS-TR steady-state thermal model | Spooner et al. (2020), [doi:10.5880/GFZ.4.5.2020.007](https://doi.org/10.5880/GFZ.4.5.2020.007), CC BY 4.0 |
| Moho | European and Adriatic Moho (5 km grid) | Spada et al. (2013, GJI) |
| Local 3D geology | Nappe de Morcles (18 formations, 50 m voxels) | Thornton et al. (2018), [doi:10.6084/m9.figshare.c.4130759](https://doi.org/10.6084/m9.figshare.c.4130759), CC0 |
| Gravity | Gravimetric Atlas of Switzerland 1:100 000 (points, complete Bouguer anomaly, 2670 kg/m³); Gravimetric Map 1:500 000 (Bouguer contours); National Gravity Network LSN2004 (absolute *g*) | swisstopo, map.geo.admin.ch layers `ch.swisstopo.geologie-gravimetrischer_atlas.messpunkte`, `ch.swisstopo.geologie-geodaesie-bouguer_anomalien`, `ch.swisstopo.landesschwerenetz`, OGD |
| Seismic subsoil classes | SIA 261 classes A–F, national map incl. the Valais 1:25 000 maps of the Rhône valley (CREALP for the canton, [vs.ch](https://www.vs.ch/web/sdm/cartes-de-sol-de-fondation)); spectral microzonation perimeters | BAFU, `ch.bafu.gefahren-baugrundklassen` (release 4 Mar 2026), `ch.bafu.gefahren-spektral`, OGD |

The swisstopo **swissAlps3D (SA3D)** model is not yet publicly downloadable (2026); it can be added
when released. Data are downloaded on first use to `data/raw/` (or `$GEOVALAIS_DATADIR/raw` if set, e.g. on a cluster scratch file system; first run ≈ 5 GB transferred,
mostly swissBEDROCK tiles; the cache is small); outputs go to `data/output/`. Both are git-ignored.

## Model

- **Grid**: rectilinear Cartesian grid in Swiss LV95 (EPSG:2056), km (`x` easting, `y` northing,
  `z` elevation). Uniform horizontal spacing `dx`; vertical spacing `dz_fine` (default 25 m) in the
  elevation band of the Quaternary valley fill (derived from swissBEDROCK, ≈ −0.5 … +1.6 km a.s.l.),
  growing geometrically to `dz_max_up` above and `dz_max_down` at depth (`refined_levels`).
- **Phases** (`PHASES`): `0` air, `1` Quaternary deposits, `2` consolidated sediments
  (Alpine sedimentary nappes), `3` upper crust, `4` lower crust, `5` mantle,
  `6` Adriatic mantle wedge (Ivrea body).
  - Quaternary base = topography − swissBEDROCK thickness (0 outside Switzerland).
    With `min_fill = 0.5` (default 0), deposits that fill less than half of the column's top cell are dropped.
  - Base of consolidated sediments = 3D-ALPS top of upper crust (≤ bedrock).
  - Upper/lower crust boundary at the 3D-ALPS fraction of crystalline-crust thickness, scaled to the Spada Moho.
  - Inside the Morcles model box, its formations override the geology below the bedrock surface
    (basement → upper crust, sedimentary units → consolidated sediments; field `MorclesUnit`).
- **Fields**: `Phases`, `Vp`, `Vs` (km/s), `Rho` (kg/m³; 3D-ALPS densities), `Temperature` (°C, 3D-ALPS-TR),
  `MorclesUnit`. Velocities are linear in depth below the surface per phase (`DEFAULT_VELOCITIES`,
  `LayerVelocity`) — **first-order values to be calibrated**, e.g. against Rhône valley Vs profiles.
  Air cells (`Phases = 0`) have Vp = Vs = Rho = Temperature = NaN (`air_value` keyword; use `air_value = 0`
  for solvers that do not accept NaN).
- **Surface dataset**: topography, bedrock elevation, Quaternary thickness, tectonic units
  (3 levels), model interfaces and Moho depths. Legend of all codes in `data/output/legend.csv`.

## Usage

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t auto examples/build_valais_model.jl
```

```julia
using GeoValaisWallis3D            # re-exports GeophysicalModelGenerator
res = build_valais_model(valais_region(); dx = 0.25, dz_fine = 0.025)
res.model      # CartData: Phases, Vp, Vs, Rho, Temperature, MorclesUnit
res.surface    # CartData: surface maps
res.z          # vertical levels (km)
```

A model too large for one process can be built block by block (e.g. one block per MPI rank): the 2D maps once, then any
range of columns, identical to that part of the whole model.

```julia
maps = valais_model_maps(valais_region(); dx = 0.05, dz_fine = 0.025)   # axes x, y, z; interfaces; Morcles model
block = valais_model_block(maps, 1:650, 1:475; air_value = 0)            # Phases, Vp, Vs, Rho, MorclesUnit of these columns
```

Open `data/output/*.vts` in ParaView. To hide the air: colour by e.g. `Vs`, open the Color Map Editor,
click the gear icon (advanced properties) and set **NaN Opacity** to 0 — or apply a Threshold on `Phases` ≥ 1. The `unzip` command is needed for the 3D-ALPS-TR and Morcles archives.

## Gravity anomalies

```bash
julia --project=. -t auto examples/gravity_maps.jl      # ~20 s, after a ~10 s download
```

`gravity_maps` grids the swisstopo Bouguer anomaly on the model surface (cell centres at `dx`), twice:
from the 1:100 000 atlas points (`Bouguer_mGal`, inverse-distance weighting, NaN more than 5 km from a
point) and from the 1:500 000 contour lines (`Bouguer500k_mGal`, smoother, and complete on the Italian
side). swisstopo publishes no free-air map, so the free-air anomaly is computed as **FA = BA + g_topo**.
g_topo is the attraction of the topography between sea level and the swissALTIRegio surface at
2670 kg/m³, i.e. the complete Bouguer correction. It is computed with exact prisms next to each point,
line masses at 250 m out to 10 km, and at 1 km out to 100 km, with Earth curvature. At the 177 LSN2004
stations, where the free-air anomaly is known exactly (g − γ, GRS80, LN02 heights), the map method is
off by −3.6 ± 2.6 mGal (max 10 mGal), on a free-air range of −134…206 mGal. The mean offset probably
comes from the terrain beyond 100 km (the atlas reduction goes to 167 km) and the normal-gravity
formula of the atlas.

`write_gravity_paraview` writes `data/output/valais_gravity.vts`: the surface, raised 20 m to draw
over `valais_surface.vts`, with `Bouguer_mGal`, `Bouguer500k_mGal`, `FreeAir_mGal`,
`FreeAir500k_mGal`, `TopoAttraction_mGal`. It also writes the stations as points:
`valais_gravity_atlas.vtu` and `valais_gravity_network.vtu` (exact and map free-air anomaly, misfit).

### Density inversion

GeophysicalModelGenerator ≥ 0.10.1 provides the adjoint of its gravity code (`voxel_grav_adjoint`), and
`normal_gravity` (GRS80), which this package now uses instead of its own. On top of them:

- `density_on_grid` averages a density model onto the nodes of a coarser grid. Spacing may vary, and the
  grid may be padded laterally beyond the model.
- `bouguer_contrast` turns it into the contrast a complete Bouguer anomaly sees: ρ − 2670 above sea level,
  ρ − the mean of each depth level below it.
- `invert_density` solves for the density change that explains a Bouguer anomaly at arbitrary stations,
  using LSQR (Krylov.jl) as in the GMG documentation. It uses Li & Oldenburg depth weighting, adds a
  regional plane, and stops when the misfit reaches the data error. It runs on a GPU with `backend`.
- `interpolate_nodes` maps the result back onto the model grid.

ValaisRawilWaves (`scripts/07_density_inversion.jl`) applies this to the wave-propagation model and the
atlas points.

## Seismic subsoil classes

```bash
julia --project=. examples/subsoil_classes.jl      # ~20 s, after an 80 MB download
```

`subsoil_map` rasterises the SIA 261 seismic subsoil classes (A–F) at 50 m on the surface grid, together
with the perimeters of the spectral microzonation studies. Inside those perimeters (Sion, Martigny,
Brig-Visp, Chablais, and in Vaud: Plaine du Rhône, Riviera), the class map has gaps, because the
microzonation replaces it; their maps and design spectra are PDFs and spreadsheets on the cantonal page.
`write_subsoil_paraview` writes `data/output/valais_subsoil.vts`: the surface raised 10 m, with
`SubsoilClass` (1–6 = A–F, NaN where not mapped), `SIA261_revision` (2003/2020) and
`SpectralMicrozonation`. The legend is in `subsoil_legend.csv`. In ParaView, colour by `SubsoilClass`
with an indexed colour map and set NaN Opacity to 0. The classes are indicative (1:25 000): they don't
replace a site study.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'                         # offline
GVW_ONLINE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'   # with downloads
```

## References

- Spada, M., Bianchi, I., Kissling, E., Piana Agostinetti, N., Wiemer, S. (2013). *GJI* 194, 1050–1068.
- Spooner, C., Scheck-Wenderoth, M., Götze, H.-J., Ebbing, J., Hetényi, G. (2019). 3D-ALPS. GFZ Data Services.
- Spooner, C., Scheck-Wenderoth, M., Cacace, M., Anikiev, D. (2020). 3D-ALPS-TR. GFZ Data Services.
- Thornton, J. M., Mariethoz, G., Brunner, P. (2018). A 3D geological model of a structurally complex Alpine region as a basis for interdisciplinary research. *Sci. Data* 5, 180238.
- Wagner, M., Kissling, E., Husen, S. (2012). *GJI* 191, 789–802.
- swisstopo: swissALTIRegio; swissBEDROCK; Tectonic Map of Switzerland 1:500 000.
