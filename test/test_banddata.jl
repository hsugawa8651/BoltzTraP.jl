# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - BandData tests for v0.3 magnetic support

@testset "BandData" begin
    # Create minimal test data
    lattice = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]  # Simple cubic, Bohr
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'  # 2 atoms, 3×2
    species = [:Si, :Si]
    kpoints = [0.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.5 0.0]'  # 3 k-points, 3×3
    fermi = 0.2
    nelect = 8.0

    @testset "Unpolarized construction" begin
        # Unpolarized: nspin=1, no magmom
        ebands = rand(4, 3)  # 4 bands, 3 k-points (2D for unpolarized)

        data = BandData(;
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            ebands = ebands,
            fermi = fermi,
            nelect = nelect,
        )

        @test data isa BandData{Unpolarized}
        @test BoltzTraP.get_dosweight(data) == 2.0
        @test BoltzTraP.get_nspinor(data) == 1
        @test isnothing(data.magmom)
    end

    @testset "Collinear construction" begin
        # Collinear: nspin=2, scalar magmom per atom
        ebands = rand(4, 3, 2)  # 4 bands, 3 k-points, 2 spins (3D)
        magmom = [1.0, -1.0]  # 2 atoms, scalar moments

        data = BandData(;
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            ebands = ebands,
            fermi = fermi,
            nelect = nelect,
            magmom = magmom,
        )

        @test data isa BandData{Collinear}
        @test BoltzTraP.get_dosweight(data) == 1.0
        @test BoltzTraP.get_nspinor(data) == 2
        @test data.magmom == [1.0, -1.0]
    end

    @testset "NonCollinear construction" begin
        # NonCollinear: vector magmom per atom
        ebands = rand(4, 3)  # 4 bands, 3 k-points
        magmom = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]  # 2 atoms, 3D vectors

        data = BandData(;
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            ebands = ebands,
            fermi = fermi,
            nelect = nelect,
            magmom = magmom,
        )

        @test data isa BandData{NonCollinear}
        @test BoltzTraP.get_dosweight(data) == 1.0
        @test BoltzTraP.get_nspinor(data) == 1
        @test length(data.magmom) == 2
        @test data.magmom[1] == [0.0, 0.0, 1.0]
    end

    @testset "Error for invalid input" begin
        # Invalid type should throw
        @test_throws ArgumentError BandData(123)
    end

    @testset "Accessors" begin
        ebands = rand(4, 3)
        data = BandData(;
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            ebands = ebands,
            fermi = fermi,
            nelect = nelect,
        )

        @test size(data.ebands) == (4, 3)
        @test size(data.lattice) == (3, 3)
        @test size(data.positions) == (3, 2)
        @test data.fermi == fermi
        @test data.nelect == nelect
    end
end
