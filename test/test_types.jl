# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Tests for DFTData{NSpin} type.
"""

using Test
using BoltzTraP

@testset "DFTData Type" begin

    # Test data
    lattice = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]
    positions = [0.0 0.25; 0.0 0.25; 0.0 0.25]  # 3×2
    species = ["Si", "Si"]
    kpoints = [0.0 0.5; 0.0 0.5; 0.0 0.5]  # 3×2
    weights = [0.5, 0.5]
    ebands_nspin1 = reshape([-0.1, 0.1, -0.15, 0.05], 2, 2, 1)  # 2 bands, 2 kpts, 1 spin
    ebands_nspin2 = cat(ebands_nspin1, ebands_nspin1 .+ 0.01; dims=3)  # 2 spins
    occupations_nspin1 = reshape([1.0, 0.0, 1.0, 0.0], 2, 2, 1)
    occupations_nspin2 = cat(occupations_nspin1, occupations_nspin1; dims=3)
    fermi = 0.0
    nelect = 4.0

    @testset "Type construction" begin
        # Non-spin-polarized (NSpin=1)
        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )
        @test data1 isa DFTData{1}

        # Spin-polarized (NSpin=2)
        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
        )
        @test data2 isa DFTData{2}
    end

    @testset "NSpin inference from ebands" begin
        # NSpin should be inferred from size(ebands, 3)
        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )
        @test nspin(data1) == 1

        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
        )
        @test nspin(data2) == 2
    end

    @testset "Accessor functions" begin
        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )

        @test nspin(data1) == 1
        @test is_magnetic(data1) == false

        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
        )

        @test nspin(data2) == 2
        @test is_magnetic(data2) == true
    end

    @testset "Type aliases" begin
        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )
        @test data1 isa NonMagneticData

        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
        )
        @test data2 isa SpinPolarizedData
    end

    @testset "Optional magmom field" begin
        # Without magmom
        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )
        @test isnothing(data1.magmom)

        # With magmom
        magmom = [1.0, -1.0]
        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
            magmom = magmom,
        )
        @test data2.magmom == magmom
    end

    @testset "Field access" begin
        data = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )

        @test data.lattice == lattice
        @test data.positions == positions
        @test data.species == species
        @test data.kpoints == kpoints
        @test data.weights == weights
        @test data.ebands == ebands_nspin1
        @test data.occupations == occupations_nspin1
        @test data.fermi == fermi
        @test data.nelect == nelect
    end

    @testset "Type concreteness (type stability)" begin
        # Verify that DFTData always has concrete type parameter
        # This ensures loaders don't return abstract DFTData without annotation

        data1 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin1,
            occupations = occupations_nspin1,
            fermi = fermi,
            nelect = nelect,
        )
        # Type parameter must be concrete Int, not abstract
        @test typeof(data1) == DFTData{1}
        @test isconcretetype(typeof(data1))

        data2 = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands_nspin2,
            occupations = occupations_nspin2,
            fermi = fermi,
            nelect = nelect,
        )
        @test typeof(data2) == DFTData{2}
        @test isconcretetype(typeof(data2))
    end

end
