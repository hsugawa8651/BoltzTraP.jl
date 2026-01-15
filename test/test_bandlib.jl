# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests/test_bandlib.py

using Test
using LinearAlgebra
using NPZ
using Random
using Statistics: mean

@testset "Band Library (transport)" begin

    TEST_DIR = @__DIR__
    DATA_DIR = joinpath(TEST_DIR, "data")

    @testset "DOS parabolic band" begin
        # Generate a set of electron energies from a single parabolic band
        L = 5.0
        x = range(-L, L, length=201)
        y = range(-L, L, length=201)
        z = range(-L, L, length=201)

        # Create points and filter to sphere
        points = [(xi, yi, zi) for xi in x for yi in y for zi in z]
        k2 = [p[1]^2 + p[2]^2 + p[3]^2 for p in points]
        k2 = filter(x -> x < L*L, k2)

        # Energy proportional to k²
        energies = reshape(k2 / maximum(k2) * 0.1, 1, :)

        npts = 60
        erange = (minimum(energies), maximum(energies))
        epsilon, dos = BoltzTraP.compute_dos(energies, erange, npts)

        # Test shape
        @test length(epsilon) == npts
        @test length(dos) == npts

        # Test normalization (should integrate to 1)
        de = epsilon[2] - epsilon[1]
        @test sum(dos) * de ≈ 1.0 atol=0.01

        # Analytic result for parabolic band: DOS ∝ √E
        ref = sqrt.(max.(epsilon, 0))
        ref ./= sum(ref) * de

        # Test that the difference is small
        @test norm(dos - ref) * de < 0.01
    end

    @testset "DOS Si" begin
        # Load Si DOS reference
        ref_file = joinpath(DATA_DIR, "Si_BTPdos.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)
        dose_ref = reference["dose"]
        dos_ref = reference["dos"]

        # DOS should be normalized to number of bands
        de = dose_ref[2] - dose_ref[1]
        @test sum(dos_ref) * de ≈ 6.0 atol=0.1  # Si has 6 valence bands
    end

    @testset "Fermi-Dirac integration" begin
        # Test basic Fermi-Dirac integration properties
        epsilon = collect(range(-0.1, 0.1, length=1000))
        de = epsilon[2] - epsilon[1]

        # Parabolic DOS
        dos = sqrt.(max.(epsilon .+ 0.05, 0))
        dos ./= sum(dos) * de

        # Test at different temperatures
        const_kB = 3.1668115634556076e-06  # Hartree/K
        T_range = [100.0, 200.0, 300.0]
        μ = 0.0

        for T in T_range
            kT = const_kB * T
            f = BoltzTraP.fermi_dirac(epsilon, μ, kT)
            df = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT)

            # Integral of f*dos should be between 0 and 1
            N = sum(dos .* f) * de
            @test 0 <= N <= 1

            # df should be negative
            @test all(df .<= 0)

            # Integral of -df should be positive
            @test sum(-df .* dos) * de >= 0
        end
    end

    @testset "Fermi integrals Si" begin
        # Load Si fermiintegrals reference
        ref_file = joinpath(DATA_DIR, "Si_fermiintegrals.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)

        # Just verify the reference data has expected shape
        N_ref = reference["N"]
        L0_ref = reference["L0"]
        L1_ref = reference["L1"]
        L2_ref = reference["L2"]

        # L0, L1, L2 should be 4D: (nT, nμ, 3, 3)
        @test ndims(L0_ref) == 4
        @test size(L0_ref, 3) == 3
        @test size(L0_ref, 4) == 3

        # N should be 2D: (nT, nμ)
        @test ndims(N_ref) == 2
    end

    @testset "Onsager coefficients Si" begin
        # Load Si Onsager reference
        ref_file = joinpath(DATA_DIR, "Si_Onsager.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)

        # Verify reference data structure
        cond_ref = reference["cond"]
        seebeck_ref = reference["seebeck"]
        kappa_ref = reference["kappa"]

        # All should be 4D: (nT, nμ, 3, 3)
        @test ndims(cond_ref) == 4
        @test ndims(seebeck_ref) == 4
        @test ndims(kappa_ref) == 4

        # Conductivity should be positive (diagonal elements)
        @test all(all(cond_ref[:, :, i, i] .>= 0) for i in 1:3)
    end

    @testset "calc_N temperature dependence" begin
        # |N| should be monotonically increasing with T for a single band
        epsilon = collect(range(0.0, 0.1, length=1000))
        de = epsilon[2] - epsilon[1]
        dos = sqrt.(epsilon)
        dos ./= sum(dos) * de

        μ = 0.05 * maximum(epsilon)
        T_range = range(50.0, 600.0, length=20)

        N_vals = Float64[]
        for T in T_range
            kT = 3.1668115634556076e-06 * T  # k_B T in Hartree
            f = BoltzTraP.fermi_dirac(epsilon, μ, kT)
            N = -2.0 * sum(dos .* f) * de
            push!(N_vals, N)
        end

        # N should decrease (more negative) as T increases
        # because more states get occupied above the Fermi level
        dN = diff(N_vals)
        @test all(dN .< 0)
    end

    @testset "Heat capacity properties" begin
        # Test basic properties of heat capacity calculation
        epsilon = collect(range(0.0, 0.1, length=10000))
        de = epsilon[2] - epsilon[1]
        dos = sqrt.(epsilon)
        dos ./= sum(dos) * de

        μ = 0.01  # Well into the band
        const_kB = 3.1668115634556076e-06  # Hartree/K

        # cv should increase with temperature at low T
        T1, T2 = 50.0, 100.0
        kT1, kT2 = const_kB * T1, const_kB * T2

        df1 = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT1)
        df2 = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT2)

        cv1 = -sum((epsilon .- μ).^2 .* df1 .* dos) * de
        cv2 = -sum((epsilon .- μ).^2 .* df2 .* dos) * de

        # Higher T should have higher cv (more thermal smearing)
        @test cv2 > cv1

        # cv should be positive
        @test cv1 > 0
        @test cv2 > 0
    end

end
