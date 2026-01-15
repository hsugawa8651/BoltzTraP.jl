# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Tests for DFTK extension.

These tests only run if DFTK.jl is available.
DFTK is a weak dependency, so it's not required for core functionality.
"""

using Test
using BoltzTraP

# Check if DFTK is available
const HAS_DFTK = try
    using DFTK
    true
catch
    false
end

@testset "DFTK Extension" begin
    if !HAS_DFTK
        @test_skip "DFTK not available"
    else
        @testset "load_dftk is defined" begin
            # Verify the function is exported and callable
            @test isdefined(BoltzTraP, :load_dftk)
            @test hasmethod(BoltzTraP.load_dftk, Tuple{Any})
        end

        @testset "load_dftk returns DFTData{1}" begin
            # Create a minimal DFTK calculation for testing
            # This is a very small Si calculation just for type checking

            a = 10.26  # Bohr
            lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
            Si = ElementPsp(:Si; psp=load_psp("hgh/lda/Si-q4"))
            atoms = [Si, Si]
            positions = [ones(3)/8, -ones(3)/8]

            # Very coarse settings for speed
            model = model_LDA(lattice, atoms, positions)
            basis = PlaneWaveBasis(model; Ecut=5, kgrid=[2, 2, 2])

            # Run minimal SCF
            scfres = self_consistent_field(basis; tol=1e-2, maxiter=10)

            # Test that load_dftk returns DFTData{1}
            data = load_dftk(scfres)

            @test data isa DFTData{1}
            @test isconcretetype(typeof(data))
            @test typeof(data) == DFTData{1}

            # Verify fields
            @test size(data.lattice) == (3, 3)
            @test size(data.kpoints, 1) == 3
            @test size(data.ebands, 3) == 1  # nspin = 1
            @test data.nelect > 0
            @test length(data.species) == 2
        end
    end
end
