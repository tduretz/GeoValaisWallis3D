"""
    refined_levels(zbot, ztop; fine = (-0.6, 1.6), dz_fine = 0.025,
                   dz_max_up = 0.25, dz_max_down = 2.0, growth = 1.08) -> Vector{Float64}

Vertical grid levels (km, increasing) of a rectilinear Cartesian grid, with uniform spacing `dz_fine`
inside the elevation band `fine = (zlow, zhigh)` and spacing growing geometrically by the factor `growth`
away from it, up to `dz_max_up` above (towards `ztop`) and `dz_max_down` below (towards `zbot`).

The default band covers the Quaternary fill of the Rhône valley and its tributaries: the bedrock
lies at ≈ -0.5 km a.s.l. beneath Lake Geneva and the lower Rhône valley, and valley floors
with significant sediment thickness lie below ≈ 1.6 km a.s.l. [`fine_band_from_bedrock`](@ref)
derives the band from swissBEDROCK.
"""
function refined_levels(zbot, ztop; fine = (-0.6, 1.6), dz_fine = 0.025,
                        dz_max_up = 0.25, dz_max_down = 2.0, growth = 1.08)
    zf1 = max(fine[1], zbot); zf2 = min(fine[2], ztop)
    nf = max(1, round(Int, (zf2 - zf1) / dz_fine))
    z = collect(range(zf1, zf2; length = nf + 1))
    # Spacings growing away from the band until the boundary is passed, then all spacings of
    # the section are scaled down (factor ≤ 1) to end exactly on it. Neighbouring cells therefore
    # never differ by more than the factor `growth` and never exceed `dz_max_*`, as required by
    # stretched-grid wave solvers (e.g. SeismicWavePC: ≤ 25 % change per cell).
    function section(L, dzmax)
        L <= 1e-9 && return Float64[]
        d = Float64[]; dz = dz_fine
        while sum(d) < L
            dz = min(dz * growth, dzmax)
            push!(d, dz)
        end
        return d .* (L / sum(d))
    end
    up = section(ztop - zf2, dz_max_up)
    down = section(zf1 - zbot, dz_max_down)
    zup = zf2 .+ cumsum(up)
    zdown = reverse(zf1 .- cumsum(down))
    z = vcat(zdown, z, zup)
    z[1] = max(z[1], zbot); z[end] = min(z[end], ztop)
    return z
end

"""
    fine_band_from_bedrock(topo, tmud; min_thickness = 0.02, q = 0.99, margin = 0.05) -> (zlow, zhigh)

Elevation band (km) containing the unconsolidated deposits: from the deepest bedrock elevation
(`topo - tmud`) to the `q`-quantile of the surface elevation where the sediment thickness
exceeds `min_thickness` (both in km), extended by `margin`.
"""
function fine_band_from_bedrock(topo, tmud; min_thickness = 0.02, q = 0.99, margin = 0.05)
    m = .!isnan.(tmud) .& .!isnan.(topo) .& (tmud .> min_thickness)
    any(m) || error("no unconsolidated deposits found")
    zlow = minimum(topo[m] .- tmud[m]) - margin
    s = sort(topo[m])
    zhigh = s[clamp(round(Int, q * length(s)), 1, length(s))] + margin
    return zlow, zhigh
end
