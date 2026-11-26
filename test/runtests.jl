using GeoValaisWallis3D
using Test, Statistics

@testset "GeoValaisWallis3D.jl" begin
    @testset "region" begin
        r = valais_region()
        @test r.E[1] < r.E[2] && r.N[1] < r.N[2] && r.z[1] < r.z[2]
    end

    @testset "LV95 <-> lon/lat" begin
        # Matterhorn summit: 7.65861 E, 45.97639 N  ~  LV95 2617.05 / 1091.65 km
        E, N = lonlat_to_lv95(7.65861, 45.97639)
        @test isapprox(E, 2617.05; atol = 0.1)
        @test isapprox(N, 1091.65; atol = 0.1)
        lon, lat = lv95_to_lonlat(E, N)
        @test isapprox(lon, 7.65861; atol = 1e-6) && isapprox(lat, 45.97639; atol = 1e-6)
    end

    @testset "scattered interpolation" begin
        xp = [0.0, 10.0, 0.0, 10.0]; yp = [0.0, 0.0, 10.0, 10.0]; vp = [-30.0, -30.0, -30.0, -30.0]
        out = interpolate_scattered(xp, yp, vp, 0.0:5.0:10.0, 0.0:5.0:10.0; radius = 8.0)
        @test all(out .≈ -30.0)
        out = interpolate_scattered(xp, yp, vp, [50.0], [50.0]; maxdist = 5.0)
        @test isnan(out[1])
    end

    @testset "refined vertical levels" begin
        z = refined_levels(-70.0, 5.0; fine = (-0.6, 1.6), dz_fine = 0.025, dz_max_up = 0.25, dz_max_down = 2.0)
        @test issorted(z) && z[1] ≈ -70.0 && z[end] ≈ 5.0
        dz = diff(z)
        inband = (z[1:end-1] .>= -0.6 - 1e-9) .& (z[2:end] .<= 1.6 + 1e-9)
        @test all(isapprox.(dz[inband], 0.025; atol = 1e-6))
        @test maximum(dz) <= 2.0 + 1e-9
        @test all(dz[z[1:end-1] .>= 1.6] .<= 0.25 + 1e-9)
        # band from bedrock: deepest bedrock and highest significant valley fill
        topo = [0.5 1.0; 2.0 3.0]; tmud = [0.9 0.0; 0.05 0.0]
        lo, hi = fine_band_from_bedrock(topo, tmud; margin = 0.0, q = 1.0)
        @test lo ≈ -0.4 && hi ≈ 2.0
    end

    @testset "gravity: normal gravity and attraction of masses" begin
        # GRS80 normal gravity: 978032.68 mGal at the equator, free-air gradient ≈ 0.3086 mGal/m
        @test isapprox(normal_gravity(0.0, 0.0), 978032.677; atol = 1e-3)
        @test isapprox(normal_gravity(46.0, 0.0) - normal_gravity(46.0, 0.1), 30.85; atol = 0.05)
        G = GeoValaisWallis3D
        # a distant prism acts like a vertical line mass
        p = G._prism_gz(4750.0, 5250.0, -250.0, 250.0, -1000.0, 0.0, 2670.0)
        l = G._column_gz(500.0^2, 2670.0, 5000.0^2, 0.0, -1000.0, 0.0)
        @test isapprox(p, l; rtol = 0.01)
        # 1 km thick slab, 50 × 50 km, station on top: close to the Bouguer plate 2πGρt = 111.96 mGal
        s = 0.0
        for i in -100:99, j in -100:99
            x, y = (i + 0.5) * 250, (j + 0.5) * 250
            s += (abs(x) < 400 && abs(y) < 400) ? G._prism_gz(x - 125, x + 125, y - 125, y + 125, -1000.0, 0.0, 2670.0) :
                 G._column_gz(250.0^2, 2670.0, x^2 + y^2, 0.0, -1000.0, 0.0)
        end
        @test 108 < s < 111.96
    end

    @testset "gravity of a density model (GMG voxel_grav)" begin
        # buried sphere on a stretched z grid: anomaly shape vs analytic (the slice-mean reference
        # only adds a near-constant shift)
        x = collect(0:0.5:40.0); y = collect(0:0.5:40.0)
        z = refined_levels(-20.0, 1.0; fine = (-1.0, 0.5), dz_fine = 0.1, dz_max_up = 0.25, dz_max_down = 1.0)
        rho = [zk > 0 ? 1.0 : 2800.0 + (hypot(xi - 20, yi - 20, zk + 6) <= 2.5 ? 300.0 : 0.0) for xi in x, yi in y, zk in z]
        g = model_gravity(x, y, z, rho; dx = 1.0, dz = 0.25)
        d = (6 + g.z_obs) * 1e3; M = 300 * 4 / 3 * π * 2.5e3^3
        ana(r) = 6.674e-11 * M * d / (r^2 + d^2)^1.5 * 1e5
        j = argmin(abs.(g.y .- 20))
        dmod = g.dg[argmin(abs.(g.x .- 20)), j] - g.dg[argmin(abs.(g.x .- 30)), j]
        @test isapprox(dmod, ana(0.0) - ana(10e3); rtol = 0.12)   # coarse voxels, 8-node averaging in voxel_grav
        # upward continuation of a constant stays constant
        @test all(upward_continue(0:1.0:10, 0:1.0:10, fill(5.0, 11, 11), 2.0) .≈ 5.0)
    end

    @testset "density inversion (GMG voxel_grav_adjoint + LSQR)" begin
        # averaging onto a coarser grid keeps a constant density and the rock fraction
        x = collect(0.125:0.25:9.875); y = copy(x); z = collect(-4.875:0.25:1.875)
        air = [zk > 1.0 for _ in x, _ in y, zk in z]
        g = density_on_grid(x, y, z, fill(2700.0, size(air)), air, 1.0:1.0:9.0, 1.0:1.0:9.0, -4.0:1.0:2.0)
        @test all(g.rho[:, :, 1:5] .≈ 2700) && all(g.rockfrac[:, :, 1:5] .≈ 1) && all(g.rockfrac[:, :, 7] .== 0)
        @test all(g.rockfrac[:, :, 6] .≈ 0.5)                                # node at 1 km: half air
        # trilinear interpolation is exact for a linear field
        xg, yg, zg = 0:2.0:10, 0:2.0:10, -6:1.0:0
        F = [xi + 2yi - 3zk for xi in xg, yi in yg, zk in zg]
        @test interpolate_nodes(xg, yg, zg, F, [1.3], [4.1], [-2.6])[1] ≈ 1.3 + 8.2 + 7.8
        # synthetic inversion: dense block at depth, data with 0.1 mGal noise, regional offset of 20 mGal
        xg = collect(-20.0:2:20); yg = copy(xg); zg = collect(-12.0:1:0)
        X, Y, Z = xyz_grid(xg, yg, zg)
        ρt = 300.0 .* ((abs.(X) .< 5) .& (abs.(Y) .< 5) .& (-9 .< Z .< -3))
        E = vec(X[:, :, 1]); N = vec(Y[:, :, 1]); H = fill(0.0, length(E))
        St = CartData(reshape(E, :, 1, 1), reshape(N, :, 1, 1), reshape(H, :, 1, 1), (h = reshape(H, :, 1, 1),))
        Grid = CartData(X, Y, Z, (z = Z,))
        dobs = vec(voxel_grav(Grid, ρt; stations = St, refMod = zeros(length(zg)), rhoTol = 0)[1]) .+ 20 .+ 0.1 .* randn(length(E))
        inv = invert_density(xg, yg, zg, zeros(size(X)), trues(size(X)), zeros(length(xg), length(yg)), (E, N, H), dobs, 0.1; itmax = 100)
        @test inv.misfit.final < 0.12                                         # at the noise level
        @test inv.misfit.start > 10 * inv.misfit.final
        @test isapprox(mean(inv.offset), 20; atol = 3)                       # the regional offset is recovered
        k = argmax(vec(maximum(inv.drho[11, 11:11, :], dims = 1)))            # depth of the maximum below the centre
        @test -10 <= zg[k] <= -2
    end

    @testset "velocity laws" begin
        l = LayerVelocity(Vp0 = 1.6, dVp = 1.6, Vs0 = 0.35, dVs = 0.9, Vpmax = 2.6, Vsmax = 1.1)
        @test all(GeoValaisWallis3D.velocity(l, 0.0) .≈ (1.6, 0.35))
        @test all(GeoValaisWallis3D.velocity(l, 10.0) .≈ (2.6, 1.1))
        for v in DEFAULT_VELOCITIES
            vp, vs = GeoValaisWallis3D.velocity(v, 1.0)
            @test vp / vs > sqrt(2)     # Lamé λ > 0
        end
    end

    @testset "Morcles voxel lookup" begin
        u = zeros(Int8, 3, 3, 2); u[2, 3, 1] = 5
        m = MorclesModel([2570.0, 2570.05, 2570.1], [1112.0, 1112.05, 1112.1], [0.025, 0.075], u)
        @test morcles_unit_at(m, 2570.05, 1112.1, 0.03) == 5
        @test morcles_unit_at(m, 2600.0, 1112.1, 0.03) == 0
    end

    @testset "model blocks" begin
        # synthetic maps: a block of columns must equal that part of the whole model
        x = collect(2590.0:0.5:2594.5); y = collect(1115.0:0.5:1118.5); z = [-60.0:5.0:-5.0; -4.0:1.0:-1.0; -0.5:0.1:3.0]
        topo = [1.0 + 0.3 * sin(i) + 0.2 * cos(j) for i in eachindex(x), j in eachindex(y)]
        f(v) = fill(v, length(x), length(y))
        I = (; topo, bed = topo .- 0.4, quat_thickness = f(0.4), base_cs = f(-3.0), conrad = f(-20.0), moho_eu = f(-40.0),
             moho_ad = [isodd(i) ? NaN : -30.0 for i in eachindex(x), j in eachindex(y)],
             rho_cs = f(2600.0), rho_uc = f(2700.0), rho_lc = f(2900.0), rho_lm = f(3300.0))
        maps = (; x, y, z, interfaces = I, morcles = nothing)
        whole = valais_model_block(maps; air_value = 0)
        @test size(whole.Vp) == (length(x), length(y), length(z))
        @test all(p -> p in whole.Phases, (PHASES.air, PHASES.quaternary, PHASES.consolidated_sediments, PHASES.upper_crust,
                                           PHASES.lower_crust, PHASES.mantle, PHASES.adriatic_mantle))
        @test all(iszero, whole.Vs[whole.Phases .== PHASES.air])
        for (ix, iy) in ((1:4, 1:8), (3:10, 2:5), (10:10, 8:8))
            part = valais_model_block(maps, ix, iy; air_value = 0)
            @test all(k -> part[k] == whole[k][ix, iy, :], keys(whole))
        end
    end

    @testset "thin deposits" begin
        z = [0.0, 0.025, 0.05, 0.075, 0.125, 0.175]          # 25 m cells, then 50 m
        h = GeoValaisWallis3D.top_cell_height
        @test h(z, 0.05) ≈ 0.025 && h(z, 0.07) ≈ 0.025 && h(z, 0.075) ≈ 0.0375 && h(z, 0.13) ≈ 0.05
        @test h(z, -1.0) ≈ 0.025 && h(z, 1.0) ≈ 0.05
        # top cell 25 m (topo 0.05; NaN counts as 0) or 50 m (topo 0.13); half a cell = 12.5 or 25 m
        topo = [0.05 0.05 0.13 0.13 0.05 NaN]
        tmud = [0.010 0.0125 0.020 0.030 NaN 0.010]
        thin = drop_thin_deposits(topo, tmud, z, 0.5)
        @test isequal(thin, [0.0 0.0125 0.0 0.030 NaN 0.0])
        @test isequal(drop_thin_deposits(topo, tmud, z, 0.0), tmud)
        # the top cell is no longer Quaternary where the deposit was dropped
        f(v) = fill(v, 1, 4)
        t = [0.05 0.05 0.13 0.13]
        mk(q) = (; x = [0.0], y = collect(1.0:4.0), z,
                 interfaces = (; topo = t, bed = t .- q, quat_thickness = q, base_cs = f(-1.0), conrad = f(-2.0),
                               moho_eu = f(-3.0), moho_ad = f(NaN), rho_cs = f(2600.0), rho_uc = f(2700.0),
                               rho_lc = f(2900.0), rho_lm = f(3300.0)), morcles = nothing)
        q = tmud[:, 1:4]
        top(P) = [P[1, j, searchsortedlast(z, t[j])] for j in 1:4]
        @test top(valais_model_block(mk(q)).Phases) == fill(PHASES.quaternary, 4)
        @test top(valais_model_block(mk(drop_thin_deposits(t, q, z, 0.5))).Phases) ==
              [PHASES.consolidated_sediments, PHASES.quaternary, PHASES.consolidated_sediments, PHASES.quaternary]
    end

    # Downloads data (swisstopo, GFZ, figshare, Spada et al.) -> opt-in: GVW_ONLINE_TESTS=true
    if get(ENV, "GVW_ONLINE_TESTS", "false") == "true"
        @testset "online: build model" begin
            r = valais_region(E = (2590.0, 2600.0), N = (1115.0, 1125.0), z = (-60.0, 4.0))   # Sion
            res = build_valais_model(r; dx = 0.5, dx_surface = 0.25, dz_fine = 0.05, morcles = false, write = false)
            P = res.model.fields.Phases
            @test PHASES.quaternary in P && PHASES.mantle in P
            @test maximum(res.surface.fields.QuaternaryThickness) > 0.3      # Rhône valley fill near Sion
            @test minimum(res.model.fields.Vs[P .== PHASES.quaternary]) < 0.5
            # a block built on its own = that part of the model
            maps = valais_model_maps(r; dx = 0.5, dz_fine = 0.05, morcles = false)
            part = valais_model_block(maps, 5:12, 3:20)
            @test maps.z == res.z && part.Phases == P[5:12, 3:20, :] && isequal(part.Vs, res.model.fields.Vs[5:12, 3:20, :])
        end
    end
end
