# GeoValaisWallis3D

3D geological / seismic velocity model of the Valais/Wallis region (Swiss Alps), from the
topographic surface down to below the Moho, built with
[GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

## Data sources

| Component | Dataset | Source / licence |
|---|---|---|
| Topography | swissALTIRegio (10 m DEM) — Valais window read over HTTP | swisstopo, [data.geo.admin.ch](https://data.geo.admin.ch/browser/index.html#/collections/ch.swisstopo.swissaltiregio), OGD |
| Quaternary deposits | swissBEDROCK R1 (Aug. 2025): bedrock elevation, thickness of unconsolidated deposits (10 m) | swisstopo, [swissBEDROCK](https://www.swisstopo.admin.ch/en/3d-swissbedrock-en), OGD |
| Surface geology | Tectonic Map of Switzerland 1:500 000 | swisstopo, `ch.swisstopo.geologie-tektonische_karte`, OGD |
| Crustal layers, density | 3D-ALPS (20 km grid): sediments, upper/lower crust, lithospheric mantle | Spooner et al. (2019), [doi:10.5880/GFZ.4.5.2019.004](https://doi.org/10.5880/GFZ.4.5.2019.004), CC BY 4.0 |
| Temperature | 3D-ALPS-TR steady-state thermal model | Spooner et al. (2020), [doi:10.5880/GFZ.4.5.2020.007](https://doi.org/10.5880/GFZ.4.5.2020.007), CC BY 4.0 |
| Moho | European and Adriatic Moho (5 km grid) | Spada et al. (2013, GJI) |
| Local 3D geology | Nappe de Morcles (18 formations, 50 m voxels) | Thornton et al. (2018), [doi:10.6084/m9.figshare.c.4130759](https://doi.org/10.6084/m9.figshare.c.4130759), CC0 |

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
