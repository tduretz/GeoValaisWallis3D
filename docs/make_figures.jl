# Figures of the README: the data sets combined in the model, and the model itself.
#     julia --project=docs -t auto docs/make_figures.jl [figure ...]
# Figures: overview layers sections subsoil gravity morcles (default: all). Output: docs/figures/*.jpg (needs ImageMagick).
# Uses the cached data of data/raw/ (downloads what is missing, as the examples do).
using GeoValaisWallis3D, CairoMakie, Statistics

const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)
const WHICH = isempty(ARGS) ? ["overview", "layers", "sections", "subsoil", "gravity", "morcles"] : ARGS
const R = valais_region()
CairoMakie.activate!(type = "png", px_per_unit = 1.25)
set_theme!(fontsize = 13, Axis = (titlesize = 14, titlefont = :regular, xgridvisible = false, ygridvisible = false),
           Axis3 = (titlefont = :regular,))

# places for orientation (LV95, km)
const PLACES = [("Monthey", 2562.0, 1122.5), ("Martigny", 2571.5, 1105.5), ("Sion", 2594.0, 1120.0),
                ("Sierre", 2607.0, 1127.0), ("Visp", 2634.5, 1127.0), ("Brig", 2642.0, 1129.5),
                ("Zermatt", 2624.5, 1097.0), ("Leukerbad", 2614.5, 1136.5)]
# cross-sections of the "sections" figure
const SEC_A = (name = "A", E = 2594.0)          # N–S through Sion
const SEC_B = (name = "B", N = 1120.0)          # W–E along the Rhône valley at Sion

function hillshade(z, dx; az = 315.0, alt = 45.0, exag = 1.0)
    z = replace(z, NaN => 0.0) .* exag
    gx = similar(z); gy = similar(z)
    nx, ny = size(z)
    for j in 1:ny, i in 1:nx
        gx[i, j] = (z[min(i + 1, nx), j] - z[max(i - 1, 1), j]) / ((min(i + 1, nx) - max(i - 1, 1)) * dx)
        gy[i, j] = (z[i, min(j + 1, ny)] - z[i, max(j - 1, 1)]) / ((min(j + 1, ny) - max(j - 1, 1)) * dx)
    end
    a, e = deg2rad(az), deg2rad(alt)
    lx, ly, lz = cos(e) * sin(a), cos(e) * cos(a), sin(e)
    @. (-gx * lx - gy * ly + lz) / sqrt(gx^2 + gy^2 + 1)
end
shade!(ax, x, y, topo; kw...) = heatmap!(ax, x, y, hillshade(topo, x[2] - x[1]; kw...), colormap = :grays, colorrange = (0.2, 1.05))
function places!(ax; color = :black, names = true)
    for (n, e, nn) in PLACES
        scatter!(ax, [e], [nn], color = color, markersize = 6, strokecolor = :white, strokewidth = 1)
        names && text!(ax, e, nn, text = n, offset = (4, 3), fontsize = 12, color = color, glowcolor = (:white, 0.8), glowwidth = 2)
    end
end
mapaxis(f; kw...) = Axis(f; aspect = DataAspect(), xlabel = "E (km, LV95)", ylabel = "N (km)", kw...)
limits_region!(ax) = limits!(ax, R.E..., R.N...)
# JPEG (quality 85): the shaded relief makes PNGs of 3–8 MB; as JPEG they stay below ~1 MB for the README
function savefig(name, fig)
    png = joinpath(tempdir(), name * ".png"); save(png, fig)
    jpg = joinpath(FIGDIR, name * ".jpg")
    run(`magick $png -quality 85 $jpg`)            # ImageMagick
    println("written docs/figures/$name.jpg ($(round(filesize(jpg) / 1e6, digits = 2)) MB)")
end
nanmask(a, m) = (b = Float64.(a); b[m] .= NaN; b)

# --------------------------------------------------------------------------------------------- overview
if "overview" in WHICH
    dx = 0.1
    x, y, topo = fetch_topography(R; dx)
    t = tectonic_map_on_grid(x, y)
    m = read_morcles()
    fig = Figure(size = (1500, 820))
    ax = mapaxis(fig[1, 1], title = "Topography (swissALTIRegio) and tectonic units (swisstopo 1:500 000, level 2)")
    shade!(ax, x, y, topo)
    used = sort(unique(t.level2[t.level2 .> 0]))
    cols = Makie.to_colormap(:tab20)
    code = zeros(size(t.level2)) .* NaN
    for (k, c) in enumerate(used)
        code[t.level2 .== c] .= k
    end
    heatmap!(ax, x, y, code, colormap = cgrad(cols[1:length(used)], length(used), categorical = true),
             colorrange = (0.5, length(used) + 0.5), alpha = 0.55, nan_color = :transparent)
    contour!(ax, x, y, topo, levels = [0.372], color = (:navy, 0.8), linewidth = 0.8)       # Lake Geneva shore
    mx, my = extrema(m.x), extrema(m.y)
    lines!(ax, [mx[1], mx[2], mx[2], mx[1], mx[1]], [my[1], my[1], my[2], my[2], my[1]], color = :black, linewidth = 2)
    text!(ax, mx[2], my[2], text = "Morcles 3D model", align = (:left, :top), fontsize = 11, offset = (3, -3), glowcolor = (:white, 0.8), glowwidth = 2)
    lines!(ax, [SEC_A.E, SEC_A.E], collect(R.N), color = :black, linestyle = :dash)
    lines!(ax, collect(R.E), [SEC_B.N, SEC_B.N], color = :black, linestyle = :dash)
    text!(ax, SEC_A.E, R.N[2], text = "A'", align = (:left, :top), offset = (4, -4)); text!(ax, SEC_A.E, R.N[1], text = "A", align = (:left, :bottom), offset = (4, 4))
    text!(ax, R.E[1], SEC_B.N, text = "B", align = (:left, :bottom), offset = (4, 4)); text!(ax, R.E[2], SEC_B.N, text = "B'", align = (:right, :bottom), offset = (-4, 4))
    places!(ax)
    limits_region!(ax)
    Legend(fig[1, 2], [PolyElement(color = (cols[k], 0.75)) for k in eachindex(used)], t.names2[used], "tectonic unit",
           framevisible = false, labelsize = 11, rowgap = 0)
    colsize!(fig.layout, 1, Relative(0.78))
    savefig("overview", fig)
end

# ------------------------------------------------------------------------------------- interfaces and layers
if "layers" in WHICH
    dx = 0.25
    x, y, topo = fetch_topography(R; dx)
    _, _, _, tmud = fetch_swissbedrock(R; dx)
    I = model_interfaces(x, y, topo, tmud)
    fig = Figure(size = (1500, 760))
    panels = [
        ("Topography (swissALTIRegio)", topo, :batlowW, (0.3, 4.2), "elevation (km)"),
        ("Quaternary thickness (swissBEDROCK)", nanmask(1e3 .* I.quat_thickness, I.quat_thickness .<= 0.002), Reverse(:lapaz), (0, 600), "thickness (m)"),
        ("Base of consolidated sediments (3D-ALPS)", I.base_cs, :batlow, (-12, 3), "elevation (km)"),
        ("Upper / lower crust boundary (3D-ALPS, Spada Moho)", I.conrad, :batlow, (-35, -10), "elevation (km)"),
        ("European Moho (Spada et al. 2013)", I.moho_eu, :batlow, (-60, -28), "elevation (km)"),
        ("Adriatic Moho (Spada et al. 2013)", I.moho_ad, :batlow, (-60, -10), "elevation (km)"),
    ]
    for (k, (title, A, cm, cr, lab)) in enumerate(panels)
        r, c = (k - 1) ÷ 3 + 1, (k - 1) % 3 + 1
        g = fig[r, c] = GridLayout()
        ax = mapaxis(g[1, 1]; title)
        shade!(ax, x, y, topo)
        hm = heatmap!(ax, x, y, A, colormap = cm, colorrange = cr, alpha = k == 1 ? 0.55 : 0.8, nan_color = :transparent)
        k >= 3 && contour!(ax, x, y, A, levels = range(cr..., step = k >= 5 ? 5 : 2.5), color = (:black, 0.5), linewidth = 0.6)
        places!(ax; names = false)
        limits_region!(ax)
        Colorbar(g[1, 2], hm, label = lab, height = Relative(0.7))
        r == 1 || hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
    end
    savefig("layers", fig)
end

# ------------------------------------------------------------------------------------------ cross-sections
if "sections" in WHICH
    maps = valais_model_maps(R; dx = 0.25, dz_fine = 0.025)
    ia = argmin(abs.(maps.x .- SEC_A.E)); jb = argmin(abs.(maps.y .- SEC_B.N))
    A = valais_model_block(maps, ia:ia, eachindex(maps.y); temperature = true)
    B = valais_model_block(maps, eachindex(maps.x), jb:jb; temperature = true)
    sq(a) = dropdims(a, dims = findfirst(==(1), size(a)))
    phcols = [:white, "#e9c46a", "#8fb8de", "#c0504d", "#7f6a93", "#3d5a40", "#2a9d8f"]
    phnames = ["air", "Quaternary", "consolidated sediments", "upper crust", "lower crust", "mantle", "Adriatic mantle"]
    phcm = cgrad(phcols, 7, categorical = true)
    fig = Figure(size = (1500, 1250))
    secs = [(A, maps.y, "A–A' (N–S, E = $(SEC_A.E) km)", "N (km)"), (B, maps.x, "B–B' (W–E, N = $(SEC_B.N) km)", "E (km)")]
    for (c, (S, h, name, xl)) in enumerate(secs)
        ph = Float64.(sq(S.Phases)); air = ph .== 0
        ax1 = Axis(fig[1, c], title = "$name: phases", ylabel = "elevation (km)")
        heatmap!(ax1, h, maps.z, ph, colormap = phcm, colorrange = (-0.5, 6.5))
        ax2 = Axis(fig[2, c], title = "$name: Vs (first-order velocity laws)", ylabel = "elevation (km)")
        h2 = heatmap!(ax2, h, maps.z, nanmask(sq(S.Vs), air), colormap = :lipari, colorrange = (0.3, 4.6), nan_color = :white)
        ax3 = Axis(fig[3, c], title = "$name: temperature (3D-ALPS-TR)", ylabel = "elevation (km)")
        h3 = heatmap!(ax3, h, maps.z, nanmask(sq(S.Temperature), air), colormap = :lajolla, colorrange = (0, 1200), nan_color = :white)
        contour!(ax3, h, maps.z, nanmask(sq(S.Temperature), air), levels = 200:200:1200, color = (:black, 0.5), linewidth = 0.6, labels = true, labelsize = 9)
        ax4 = Axis(fig[4, c], title = "$name: phases near the surface (vertically stretched)", ylabel = "elevation (km)", xlabel = xl)
        heatmap!(ax4, h, maps.z, ph, colormap = phcm, colorrange = (-0.5, 6.5))
        for ax in (ax1, ax2, ax3)
            ylims!(ax, -70, 5); xlims!(ax, extrema(h)...)
        end
        xlims!(ax4, extrema(h)...); ylims!(ax4, -1.0, 4.6)
        c == 2 && (Colorbar(fig[2, 3], h2, label = "Vs (km/s)"); Colorbar(fig[3, 3], h3, label = "T (°C)"))
    end
    Legend(fig[1, 3], [PolyElement(color = phcols[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:7], phnames[2:7], framevisible = false)
    Legend(fig[4, 3], [PolyElement(color = phcols[k], strokecolor = :grey60, strokewidth = 0.5) for k in 2:4], phnames[2:4], framevisible = false)
    savefig("sections", fig)
end

# ------------------------------------------------------------------------------ seismic subsoil classes (CREALP)
if "subsoil" in WHICH
    sm = subsoil_map(R; dx = 0.05)
    clcols = ["#3b5b92", "#5aa0c8", "#9ccf7a", "#f2c14e", "#e8743b", "#b5179e"]   # A (rock) … E (soft), F (special)
    clcm = cgrad(clcols, 6, categorical = true)
    cls = nanmask(sm.class, sm.class .== 0)
    ssm = Float64.(sm.ssm .> 0)
    fig = Figure(size = (1500, 1450))
    ax = mapaxis(fig[1, 1:2], title = "Seismic subsoil classes SIA 261 (BAFU national map; Rhône valley 1:25 000 maps of CREALP / canton of Valais), release $(sm.release)")
    shade!(ax, sm.x, sm.y, sm.topo)
    heatmap!(ax, sm.x, sm.y, cls, colormap = clcm, colorrange = (0.5, 6.5), nan_color = :transparent)
    contour!(ax, sm.x, sm.y, ssm, levels = [0.5], color = :black, linewidth = 1.2)
    places!(ax); limits_region!(ax)
    labels = ["A rock", "B", "C", "D", "E", "F special"]
    Legend(fig[1, 3], [[PolyElement(color = c) for c in clcols]; LineElement(color = :black)],
           [["$(SUBSOIL_CLASSES[k])" * (k == 1 ? " (rock)" : k == 5 ? " (thin soft layer on rock)" : k == 6 ? " (special)" : "") for k in 1:6];
            "spectral microzonation (classes replaced)"], "class", framevisible = false)
    # the subsoil map draped on the topography: Rhône valley from Martigny to Sierre
    ex, ny_ = (2566.0, 2612.0), (1100.0, 1132.0)
    ix = findall(e -> ex[1] <= e <= ex[2], sm.x); iy = findall(n -> ny_[1] <= n <= ny_[2], sm.y)
    xs, ys, zs = sm.x[ix], sm.y[iy], sm.topo[ix, iy]
    hs = clamp.(hillshade(zs, 0.05), 0, 1)
    rgb = map(CartesianIndices(zs)) do I
        g = 0.35 + 0.6hs[I]
        c = cls[ix[I[1]], iy[I[2]]]
        isnan(c) ? RGBf(g, g, g) : (k = Makie.to_color(clcols[Int(c)]); RGBf(k.r * (0.55 + 0.5g), k.g * (0.55 + 0.5g), k.b * (0.55 + 0.5g)))
    end
    ax3 = Axis3(fig[2, 1:3], title = "The same map draped on the topography: Rhône valley, Martigny–Sion–Sierre (vertical exaggeration 1.5, view from the south-west)",
                aspect = (ex[2] - ex[1], ny_[2] - ny_[1], 1.5 * 4.0), azimuth = -0.35π, elevation = 0.17π,
                xlabel = "E (km)", ylabel = "N (km)", zlabel = "z (km)", protrusions = 30)
    surface!(ax3, xs, ys, zs, color = rgb, shading = NoShading)
    for (n, e, nn) in PLACES
        (ex[1] <= e <= ex[2] && ny_[1] <= nn <= ny_[2]) || continue
        i = argmin(abs.(xs .- e)); j = argmin(abs.(ys .- nn))
        text!(ax3, Point3f(e, nn, zs[i, j] + 0.6), text = n, fontsize = 14, align = (:center, :bottom), color = :black, glowcolor = (:white, 0.85), glowwidth = 3)
    end
    zlims!(ax3, 0.3, 4.3)
    rowsize!(fig.layout, 1, Relative(0.5))
    savefig("subsoil", fig)
end

# --------------------------------------------------------------------------------------------- gravity
if "gravity" in WHICH
    gm = gravity_maps(R; dx = 0.25)
    fig = Figure(size = (1500, 900))
    panels = [
        ("Complete Bouguer anomaly, atlas 1:100 000 (points)", gm.bouguer, (-200, 20), Reverse(:RdBu), true),
        ("Complete Bouguer anomaly, map 1:500 000 (contours)", gm.bouguer500k, (-200, 20), Reverse(:RdBu), false),
        ("Attraction of the topography (2670 kg/m³)", gm.gtopo, (0, 300), :lajolla, false),
        ("Free-air anomaly = Bouguer + topography", gm.freeair, (-200, 200), Reverse(:RdBu), false),
    ]
    for (k, (title, A, cr, cm, pts)) in enumerate(panels)
        g = fig[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1] = GridLayout()
        local ax = mapaxis(g[1, 1]; title)
        shade!(ax, gm.x, gm.y, gm.topo)
        # diverging maps centred on zero (Bouguer: range not symmetric, so centre the colormap explicitly)
        crr = cm === :lajolla ? cr : (-1, 1) .* maximum(abs, cr)
        hm = heatmap!(ax, gm.x, gm.y, A, colormap = cm, colorrange = crr, alpha = 0.8, nan_color = :transparent)
        contour!(ax, gm.x, gm.y, A, levels = range(-300, 300, step = 20), color = (:black, 0.35), linewidth = 0.5)
        pts && scatter!(ax, gm.atlas.E, gm.atlas.N, markersize = 1.5, color = (:black, 0.5))
        k == 4 && scatter!(ax, gm.network.E, gm.network.N, marker = :utriangle, markersize = 8, color = :white, strokecolor = :black, strokewidth = 1)
        places!(ax; names = false); limits_region!(ax)
        Colorbar(g[1, 2], hm, label = "mGal", height = Relative(0.7))
    end
    Label(fig[0, 1:2], "swisstopo gravity: atlas points (dots), LSN2004 stations (triangles, free-air check: −3.6 ± 2.6 mGal)", fontsize = 14)
    savefig("gravity", fig)
end

# ------------------------------------------------------------------------------------------------ Morcles
if "morcles" in WHICH
    m = read_morcles()
    nx, ny, nz = size(m.unit)
    top = zeros(Int, nx, ny)                         # surface geology: the highest unit of each column
    for j in 1:ny, i in 1:nx
        k = findlast(!=(0), @view m.unit[i, j, :]); k === nothing || (top[i, j] = m.unit[i, j, k])
    end
    n = length(MORCLES_UNITS)
    cm = cgrad([Makie.to_color("#9e6b5a"); [get(cgrad(:roma), (k - 1) / (n - 2)) for k in 1:n-1]], n, categorical = true)
    fig = Figure(size = (1500, 700))
    ax = mapaxis(fig[1, 1], title = "Nappe de Morcles 3D model (Thornton et al. 2018): surface units")
    heatmap!(ax, m.x, m.y, nanmask(top, top .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
    jc = argmin(abs.(m.y .- (m.y[1] + m.y[end]) / 2))
    lines!(ax, [m.x[1], m.x[end]], [m.y[jc], m.y[jc]], color = :black, linestyle = :dash)
    places!(ax)
    xlims!(ax, extrema(m.x)...); ylims!(ax, extrema(m.y)...)
    ax2 = Axis(fig[1, 2], title = "W–E section at N = $(round(m.y[jc], digits = 1)) km (dashed line)", xlabel = "E (km)", ylabel = "elevation (km)", aspect = DataAspect())
    heatmap!(ax2, m.x, m.z, nanmask(m.unit[:, jc, :], m.unit[:, jc, :] .== 0), colormap = cm, colorrange = (0.5, n + 0.5), nan_color = :white)
    xlims!(ax2, extrema(m.x)...)
    Legend(fig[1, 3], [PolyElement(color = cm[k]) for k in 1:n], MORCLES_UNITS, "formation (old → young)", framevisible = false, rowgap = 0, labelsize = 11)
    colsize!(fig.layout, 2, Relative(0.55))
    savefig("morcles", fig)
end
