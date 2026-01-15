# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Tests for DFTData integration with run_interpolate/run_integrate.
"""

using Test
using LinearAlgebra
using BoltzTraP

@testset "DFTData Integration" begin

    # Create minimal test data (similar to Si)
    lattice = [
        0.0 2.715 2.715
        2.715 0.0 2.715
        2.715 2.715 0.0
    ] * BoltzTraP.ANG_TO_BOHR  # Convert to Bohr

    positions = [
        0.0 0.25
        0.0 0.25
        0.0 0.25
    ]  # 3×2 fractional

    species = ["Si", "Si"]

    # Simple 2×2 k-grid (4 k-points)
    kpoints = [
        0.0 0.5 0.0 0.5
        0.0 0.0 0.5 0.5
        0.0 0.0 0.0 0.0
    ]  # 3×4

    weights = [0.25, 0.25, 0.25, 0.25]

    # 2 bands, 4 k-points, 1 spin
    # Band energies in Hartree (simple parabolic-like)
    ebands = reshape([
        # k1    k2    k3    k4
        -0.3, -0.2, -0.2, -0.1,  # band 1 (valence-like)
        0.1, 0.2, 0.2, 0.3,   # band 2 (conduction-like)
    ], 2, 4, 1)

    occupations = reshape([
        2.0, 2.0, 2.0, 2.0,  # band 1 fully occupied
        0.0, 0.0, 0.0, 0.0,  # band 2 empty
    ], 2, 4, 1)

    fermi = 0.0  # Ha
    nelect = 8.0

    @testset "run_interpolate accepts DFTData{1}" begin
        data = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands,
            occupations = occupations,
            fermi = fermi,
            nelect = nelect,
        )

        @test data isa DFTData{1}

        # Test that method exists and is callable
        # (test data is too simple for full interpolation, so we check method dispatch)
        @test hasmethod(run_interpolate, Tuple{DFTData{1}})
    end

    @testset "DFTData and NamedTuple produce same error with insufficient data" begin
        # With insufficient k-points, both should throw the same error
        # This verifies DFTData method delegates correctly to NamedTuple method
        data_dft = DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands,
            occupations = occupations,
            fermi = fermi,
            nelect = nelect,
        )

        data_nt = (
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands,
            occupations = occupations,
            fermi = fermi,
            nelect = nelect,
        )

        # Both should throw SingularException (insufficient k-points for fitting)
        @test_throws LinearAlgebra.SingularException run_interpolate(data_dft; kpoints=10)
        @test_throws LinearAlgebra.SingularException run_interpolate(data_nt; kpoints=10)
    end

end
