# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP

@testset "Unit Constants" begin
    # Test BOHR_TO_ANG against CODATA 2018 value
    # Bohr radius a0 = 0.529177210903 Å (exact by definition in SI)
    @testset "Length conversion" begin
        @test isapprox(BoltzTraP.BOHR_TO_ANG, 0.529177210903, rtol=1e-12)
        @test isapprox(BoltzTraP.ANG_TO_BOHR, 1/0.529177210903, rtol=1e-12)
        @test isapprox(BoltzTraP.BOHR_TO_ANG * BoltzTraP.ANG_TO_BOHR, 1.0, rtol=1e-14)
    end

    # Test HA_TO_EV against CODATA 2018 value
    # 1 Hartree = 27.211386245988 eV (CODATA 2018)
    @testset "Energy conversion" begin
        @test isapprox(BoltzTraP.HA_TO_EV, 27.211386245988, rtol=1e-12)
        @test isapprox(BoltzTraP.EV_TO_HA, 1/27.211386245988, rtol=1e-12)
        @test isapprox(BoltzTraP.HA_TO_EV * BoltzTraP.EV_TO_HA, 1.0, rtol=1e-14)
    end

    # Test KB_AU (Boltzmann constant in atomic units)
    # k_B = 3.166811563e-6 Ha/K (CODATA 2018)
    @testset "Boltzmann constant" begin
        @test isapprox(BoltzTraP.KB_AU, 3.166811563e-6, rtol=1e-6)
    end

    # Test Onsager coefficient conversion factors
    @testset "Onsager conversion factors" begin
        @test BoltzTraP.SIGMA_CONV > 0
        @test BoltzTraP.SEEBECK_CONV > 0
        @test BoltzTraP.KAPPA_CONV > 0
    end
end
