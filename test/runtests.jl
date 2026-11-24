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

    # Downloads data (swisstopo, GFZ, figshare, Spada et al.) -> opt-in: GVW_ONLINE_TESTS=true
    if get(ENV, "GVW_ONLINE_TESTS", "false") == "true"
        @testset "online: build model" begin
            r = valais_region(E = (2590.0, 2600.0), N = (1115.0, 1125.0), z = (-60.0, 4.0))   # Sion
            res = build_valais_model(r; dx = 0.5, dx_surface = 0.25, dz_fine = 0.05, morcles = false, write = false)
            P = res.model.fields.Phases
            @test PHASES.quaternary in P && PHASES.mantle in P
            @test maximum(res.surface.fields.QuaternaryThickness) > 0.3      # Rhône valley fill near Sion
            @test minimum(res.model.fields.Vs[P .== PHASES.quaternary]) < 0.5
        end
    end
end
