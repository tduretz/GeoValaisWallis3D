# The model on a rectilinear grid

The first grid of the package: a rectilinear Cartesian grid in LV95 km. The same analysis on a curvilinear
(terrain-following) mesh is in [curvilinear/](../curvilinear/README.md).

## The grid

- **Horizontal:** uniform spacing `dx`.
- **Vertical:** `dz_fine` (default 25 m) in the elevation band of the Quaternary valley fill. The band is derived from
  swissBEDROCK, ≈ −0.5 … +1.6 km a.s.l. Above and below it, the spacing grows geometrically up to `dz_max_up` and
  `dz_max_down` (`refined_levels`). Neighbouring cells differ by at most the factor `growth`.
- **Air:** the levels are horizontal, so the topography is a staircase of cells. Cells above it are air
  (`Phases = 0`), with Vp = Vs = Rho = Temperature = `air_value` (default NaN; use 0 for solvers that do not accept NaN).
- **Thin deposits:** with `min_fill = 0.5` (default 0), deposits that fill less than half of the column's top cell are
  dropped.

## Usage

```bash
julia --project=. -t auto examples/rectilinear/build_valais_model.jl   # writes data/output/valais_model.vts, valais_surface.vts
```

```julia
using GeoValaisWallis3D            # re-exports GeophysicalModelGenerator
res = build_valais_model(valais_region(); dx = 0.25, dz_fine = 0.025)
res.model      # CartData: Phases, Vp, Vs, Rho, Temperature, MorclesUnit
res.surface    # CartData: surface maps
res.z          # vertical levels (km)
```

A model too large for one process can be built block by block, e.g. one block per MPI rank. First build the 2D maps
once; then any range of columns gives the same values as that part of the whole model.

```julia
maps = valais_model_maps(valais_region(); dx = 0.05, dz_fine = 0.025)   # axes x, y, z; interfaces; Morcles model
block = valais_model_block(maps, 1:650, 1:475; air_value = 0)            # Phases, Vp, Vs, Rho, MorclesUnit of these columns
```

## Figures

```bash
julia --project=plotting -e 'using Pkg; Pkg.instantiate()'
julia --project=plotting -t auto rectilinear/make_figures.jl                    # all, ~2 min with cached data
julia --project=plotting -t auto rectilinear/make_figures.jl subsoil sections   # some of: overview layers sections subsoil gravity density morcles
julia --project=plotting -t auto rectilinear/make_figures.jl showcase           # opt-in: GLMakie (needs a display) and the density cache
```

- **`density`:** runs the inversion itself. With `GVW_BACKEND=metal` it runs on the Apple GPU (Metal, Float32);
  otherwise on the CPU (Float64, much slower). The result is cached in `data/output/density_inversion_figure.jls`;
  delete that file to run the inversion again.
- **Output:** the images are written to `rectilinear/figures/*.jpg`. ImageMagick's `magick` converts them to JPEG to
  keep each file around 1 MB.

The figures are described in the [main README](../README.md#the-data-sets-in-pictures):

| | | |
|---|---|---|
| [overview](figures/overview.jpg) | [layers](figures/layers.jpg) | [sections](figures/sections.jpg) |
| [subsoil](figures/subsoil.jpg) | [gravity](figures/gravity.jpg) | [density](figures/density.jpg) |
| [morcles](figures/morcles.jpg) | [showcase](figures/showcase.jpg) | |
