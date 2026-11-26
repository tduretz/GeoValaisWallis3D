using GeoValaisWallis3D
using Test

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
