# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests/test_fd.py

using Test

# Constants matching Python BoltzTraP2
const BOLTZMANN = 3.1668115634556076e-06  # Hartree/K
const eV_to_Hartree = 1.0 / 27.211386245988

@testset "Fermi-Dirac Module" begin

    @testset "FD at 0K" begin
        # The Fermi-Dirac occupancies at 0 K should be a step function
        # Note: At exactly T=0, our implementation uses a small kT internally
        # to avoid division by zero, so use very small T instead
        e = collect(range(-5.0, 5.0, length=1001))
        kT_small = 1e-10

        # Test for mu = 0
        f = BoltzTraP.fermi_dirac(e, 0.0, kT_small)
        @test all(f[e .< -0.01] .≈ 1.0)  # well below Fermi level
        @test f[findfirst(e .≈ 0.0)] ≈ 0.5 atol=1e-5
        @test all(f[e .> 0.01] .≈ 0.0)  # well above Fermi level

        # Test for mu != 0
        mu = 3.0
        e_shifted = e .+ mu
        f = BoltzTraP.fermi_dirac(e_shifted, mu, kT_small)
        @test all(f[e_shifted .< mu - 0.01] .≈ 1.0)
        @test f[findfirst(e_shifted .≈ mu)] ≈ 0.5 atol=1e-5
        @test all(f[e_shifted .> mu + 0.01] .≈ 0.0)
    end

    @testset "FD at finite T" begin
        # The Fermi-Dirac occupancies should reproduce tabulated values
        kBT = BOLTZMANN * 300.0

        # f(0, 0, kBT) = 0.5
        f = BoltzTraP.fermi_dirac([0.0], 0.0, kBT)
        @test f[1] ≈ 0.5

        # 30 meV in Hartree
        e = 30.0 * 1e-3 * eV_to_Hartree
        f = BoltzTraP.fermi_dirac([e], 0.0, kBT)
        @test f[1] ≈ 0.23858517713617 rtol=1e-5

        f = BoltzTraP.fermi_dirac([-e], 0.0, kBT)
        @test f[1] ≈ 0.76141482286383 rtol=1e-5

        # With shifted chemical potential
        f = BoltzTraP.fermi_dirac([e], 15.0 * 1e-3 * eV_to_Hartree, kBT)
        @test f[1] ≈ 0.358880600979885 rtol=1e-5

        f = BoltzTraP.fermi_dirac([e], -15.0 * 1e-3 * eV_to_Hartree, kBT)
        @test f[1] ≈ 0.149226850069459 rtol=1e-5

        # Limits
        f = BoltzTraP.fermi_dirac([Inf], 0.0, kBT)
        @test f[1] == 0.0

        f = BoltzTraP.fermi_dirac([-Inf], 0.0, kBT)
        @test f[1] == 1.0
    end

    @testset "dFD/dx" begin
        # The Fermi-Dirac derivative should agree with tabulated values
        x = collect(range(-20.0, 20.0, length=1001))

        # Use reduced units: dFD/dx = dFD/de * kT
        kT = 1.0
        df = BoltzTraP.dfermi_dirac_de(x, 0.0, kT)

        # At x=0, dFD/dx = -0.25
        idx_zero = findfirst(x .≈ 0.0)
        @test df[idx_zero] ≈ -0.25 rtol=1e-10

        # Derivative should always be non-positive
        @test all(df .<= 0.0)

        # Maximum magnitude at x=0
        @test all(df[x .!= 0.0] .> -0.25)

        # Approaches 0 at extremes
        @test df[1] ≈ 0.0 atol=1e-8
        @test df[end] ≈ 0.0 atol=1e-8

        # Specific values
        @test BoltzTraP.dfermi_dirac_de([1.0], 0.0, 1.0)[1] ≈ -0.19661193 rtol=1e-5
        @test BoltzTraP.dfermi_dirac_de([-1.0], 0.0, 1.0)[1] ≈ -0.19661193 rtol=1e-5
    end

    @testset "dFD/de scaling" begin
        # The derivative should scale correctly with temperature
        x = collect(range(-20.0, 20.0, length=1001))
        kT = 33.0

        df_x = BoltzTraP.dfermi_dirac_de(x, 0.0, 1.0)  # In reduced units
        df_e = BoltzTraP.dfermi_dirac_de(x .* kT, 0.0, kT)

        @test all(isapprox.(df_x, df_e .* kT, rtol=1e-10))

        # With shifted mu
        mu = 10.0
        df_e = BoltzTraP.dfermi_dirac_de(x .* kT .+ mu, mu, kT)
        @test all(isapprox.(df_x, df_e .* kT, rtol=1e-10))
    end

end
