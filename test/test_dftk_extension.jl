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

        @testset "load_dftk non-magnetic (Si) returns DFTData{1}" begin
            # Create a minimal DFTK calculation for testing
            # This is a very small Si calculation just for type checking

            a = 10.26  # Bohr
            lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
            Si = ElementPsp(:Si; psp = load_psp("hgh/lda/Si-q4"))
            atoms = [Si, Si]
            positions = [ones(3)/8, -ones(3)/8]

            # Very coarse settings for speed
            model = model_LDA(lattice, atoms, positions)
            basis = PlaneWaveBasis(model; Ecut = 5, kgrid = [2, 2, 2])

            # Run minimal SCF
            scfres = self_consistent_field(basis; tol = 1e-2, maxiter = 10)

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

            # Non-magnetic: is_magnetic should return false
            @test !BoltzTraP.is_magnetic(data)
        end

        @testset "load_dftk collinear (Fe) returns DFTData{2}" begin
            # Create a minimal Fe collinear calculation
            # Fe BCC primitive cell

            a = 5.42  # Bohr (Fe BCC lattice constant)
            lattice = a / 2 * [[-1 1 1]; [1 -1 1]; [1 1 -1]]  # BCC primitive

            Fe = ElementPsp(:Fe; psp=load_psp("hgh/lda/Fe-q8"))
            atoms = [Fe]
            positions = [[0.0, 0.0, 0.0]]

            # Collinear magnetic
            magnetic_moments = [4.0]

            model = model_LDA(lattice, atoms, positions;
                              temperature=0.01,
                              magnetic_moments)

            # Very coarse settings for speed
            basis = PlaneWaveBasis(model; Ecut=10, kgrid=[2, 2, 2])

            # Run minimal SCF
            ρ0 = guess_density(basis, magnetic_moments)
            scfres = self_consistent_field(basis; ρ=ρ0, tol=1e-2, maxiter=50)

            # Test that load_dftk returns DFTData{2}
            data = load_dftk(scfres)

            @test data isa DFTData{2}
            @test isconcretetype(typeof(data))
            @test typeof(data) == DFTData{2}

            # Verify fields
            @test size(data.lattice) == (3, 3)
            @test size(data.kpoints, 1) == 3
            @test size(data.ebands, 3) == 2  # nspin = 2
            @test data.nelect > 0
            @test length(data.species) == 1
            @test data.species[1] == "Fe"

            # Collinear magnetic: is_magnetic should return true
            @test BoltzTraP.is_magnetic(data)

            # magmom should be set
            @test !isnothing(data.magmom)
            @test length(data.magmom) == 1  # 1 atom
        end
    end
end
