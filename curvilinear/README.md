# The model on a curvilinear mesh

GeophysicalModelGenerator.jl 0.10.2 can build curvilinear meshes with `curvilinear_grid` and `deform_grid`. Their
horizontal grid levels follow surfaces, and `voxel_grav` computes the gravity of such a mesh. This folder redoes the
analysis of the package on such a mesh. The rectilinear grid ([`build_valais_model`](../src/model.jl)) and its figures
in [`docs/figures`](../docs/figures) are unchanged.

## The mesh

`valais_curvilinear_maps` / `build_valais_curvilinear_model` ([src/curvilinear.jl](../src/curvilinear.jl)) build a
mesh with 250 m columns. Its levels follow these surfaces, from top to bottom (`meshing_surfaces`):

| surface | data | level | cells in the layer below (spacing where thickest) |
|---|---|---|---|
| topography | swissALTIRegio | 288 | 40 (25 m) |
| base of the near-surface layer | 1 km below the topography (≥ 250 m, ≤ half-way to the next surface) | 248 | 73 (100 m) |
| base of the consolidated sediments | 3D-ALPS / swissBEDROCK | 175 | 77 (500 m) |
| Conrad | 3D-ALPS, scaled to the Spada Moho | 98 | 59 (500 m) |
| European Moho | Spada et al. (2013) | 39 | 38 (1 km) |
| bottom | −70 km | 1 | |

The mesh has 520 × 380 × 288 = 56.9 M nodes, and none of them lies in the air.

These are the interfaces that can be meshed: they span the whole region and never cross.

- **Exact grid levels:** the topography, the Conrad, the European Moho and the bottom are grid levels in every column.
- **Base of the consolidated sediments:** this is a grid level too, except in the 3.7 % of the columns where the
  sediments pinch out (crystalline massifs). There the mesh keeps a minimum layer thickness, and the grid level lies up
  to 323 m below the interface.
- **Not meshable:**
  - The base of the Quaternary is absent in 80 % of the columns. It lies inside the 25 m layer, whose base is 1 km
    deep, while the deposits are at most 1 km thick.
  - The Adriatic Moho is a wedge that crosses the Conrad.

  Both remain boundaries between materials. The material rules are the same as on the rectilinear grid (a shared
  `_model_point`), applied at every node.

`spacing = :sbp` places the points of each layer at the nodes of a summation-by-parts grid, as used by SeismicWavePC.jl.

## Figures

```bash
GVW_BACKEND=metal julia --project=docs -t auto curvilinear/make_figures.jl          # all
julia --project=docs -t auto curvilinear/make_figures.jl mesh layers sections        # some of: overview mesh layers sections subsoil gravity density morcles showcase
```

`overview` and `showcase` render their 3D panels with GLMakie, which needs a display. `gravity` and `density` call GMG
`voxel_grav` on curvilinear grids. Their results are cached in `data/output/curvilinear_*.jls`.

**Overview.** The fence diagram of the model, with the tectonic map on the relief. The cut faces show the phases, every
4th grid level, and the meshed surfaces.

![Overview](figures/overview.jpg)

**The mesh.** Section B–B' with every 4th level and column. The zoom on the Rhône valley at Sion shows every level and
column: the Quaternary fill (yellow) lies in the 25 m levels that follow the topography. On the right, the cell height
in the column at Sion.

![Mesh](figures/mesh.jpg)

**Meshed surfaces.** Maps of the levels of the mesh, with the Quaternary thickness and the Adriatic Moho, which is not
meshed. The red outline marks the columns where the level of the base of the consolidated sediments is pushed down.

![Meshed surfaces](figures/layers.jpg)

**Sections.** Phases, Vs and temperature at the nodes of sections A–A' and B–B'. Each node is drawn as its dual cell,
so the colours follow the mesh exactly.

![Sections](figures/sections.jpg)

**Subsoil classes and the Quaternary on the mesh.** The CREALP / BAFU map, and N–S sections across the Rhône valley at
Martigny, Sion and Visp:
- Vs on every level of the mesh;
- the base of the Quaternary (white line);
- the surface class of the map at each section (squares).

![Subsoil](figures/subsoil.jpg)

**Gravity.** The attraction of the topography (sea level to the surface, 2670 kg/m³), computed by GMG `voxel_grav` on
two curvilinear layers whose top is the topography: 250 m columns out to 10 km, then 1 km columns out to 100 km.

| comparison | mean | std | range |
|---|---|---|---|
| GMG − the package's prisms and line masses, at 1 km stations | 2.3 mGal | 1.3 mGal | −1.9 to 6.1 mGal |

The GMG grids are flat, while the package accounts for the curvature of the Earth.

![Gravity](figures/gravity.jpg)

**Density inversion.** The Bouguer anomaly (8168 atlas points) is inverted on a curvilinear mesh:
- columns every 1 km, with 20 km of padding;
- levels that follow the same surfaces;
- 170 × 135 × 65 = 1.49 M unknowns, and no air nodes, hence no rock fractions.

The density of a node is the mean of the 250 m model over its dual cell, at depths below the local topography. The
reference is 2670 kg/m³ above sea level and the mean density per 500 m depth bin below it.

| inversion | unknowns | LSQR iterations | time on the Apple M2 Max GPU | rms misfit, start → final |
|---|---|---|---|---|
| curvilinear mesh (this folder) | 1.49 M | 66 | 158 s | 17.2 → 1.9 mGal |
| rectilinear grid ([docs/figures/density.jpg](../docs/figures/density.jpg)) | 1.0 M | 62 | 71 s | 18.1 → 1.9 mGal |

The starting misfit is lower here: the padding now holds the model (topography, 3D-ALPS, Moho) instead of uniform
densities. The Ivrea body again needs about +200 kg/m³.

![Density inversion](figures/density.jpg)

**Morcles.** A W–E section of the Morcles model, sampled on a 2D curvilinear mesh: GMG `curvilinear_grid` with a single
`y`, 50 m columns and 25 m levels from the topography. Compared with the 50 m voxels, the mesh follows the relief
instead of making steps.

![Morcles](figures/morcles.jpg)

**Showcase.** One image for posts: the fence diagram, the mesh at Sion and the gravity inversion.

![Showcase](figures/showcase.jpg)
