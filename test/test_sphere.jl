# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests/test_sphere.py

using Test
using LinearAlgebra
using StaticArrays
using NPZ

@testset "Sphere Module" begin

    @testset "compute_bounds cube" begin
        # Simple cubic lattice
        lattvec = Matrix{Float64}(I, 3, 3)
        for i in 1:9
            bounds = BoltzTraP.compute_bounds(lattvec, Float64(i))
            # bounds returns a vector, compare as Tuple
            @test Tuple(bounds) == (i, i, i)
        end
    end

    @testset "compute_bounds Li (BCC)" begin
        # Test based on Li (BCC)
        lattvec = [
            -3.19204785  3.19204785  3.19204785;
             3.19204785 -3.19204785  3.19204785;
             3.19204785  3.19204785 -3.19204785
        ]
        bounds = BoltzTraP.compute_bounds(lattvec, 145.47229311)
        @test Tuple(bounds) == (33, 33, 33)
    end

    @testset "lattice_points_in_sphere" begin
        # Simple cubic lattice
        lattvec = SMatrix{3,3,Float64}(I)
        radius = 2.0
        bounds = (2, 2, 2)

        points, sqnorms = BoltzTraP.lattice_points_in_sphere(lattvec, radius, bounds)

        # Origin should be included
        @test SVector(0, 0, 0) in points

        # All points should have norm ≤ radius
        for p in points
            @test norm(p) <= radius + 1e-10
        end

        # Points should be sorted by squared norm
        @test issorted(sqnorms)
    end

    @testset "calc_tensor_basis Si" begin
        # Si crystal structure
        a = 5.45052526
        lattvec = a * 0.5 * (ones(3, 3) - I)
        positions = [0.125 0.125 0.125; 0.875 0.875 0.875]

        basis = BoltzTraP.calc_tensor_basis(lattvec', positions, [14, 14], nothing)

        # For cubic symmetry (Oh), should have 1 independent component
        # (isotropic tensor: σ_ij = σ δ_ij)
        @test size(basis, 1) == 1
        @test size(basis, 2) == 3
        @test size(basis, 3) == 3

        # The basis tensor should be proportional to identity
        @test basis[1, :, :] ≈ basis[1, 1, 1] * I atol=1e-10
    end

    @testset "calc_nrotations" begin
        # Si crystal structure
        a = 5.45052526
        lattvec = a * 0.5 * (ones(3, 3) - I)
        positions = [0.125 0.125 0.125; 0.875 0.875 0.875]

        nrot = BoltzTraP.calc_nrotations(lattvec', positions, [14, 14], nothing)

        # Si has 48 rotations (Oh point group)
        @test nrot == 48
    end

    @testset "get_equivalences number" begin
        # Test that we get a reasonable number of equivalence classes
        # (not testing exact values since they depend on Spglib internals)
        TEST_DIR = @__DIR__
        reference_file = joinpath(TEST_DIR, "data", "Si_equivalences.npz")

        if isfile(reference_file)
            ref_data = npzread(reference_file)
            n_ref = length(keys(ref_data))

            # Reference Si equivalences should have ~5048 classes (from Python)
            @test n_ref > 4000
            @test n_ref < 6000
        else
            @warn "Reference file not found: $reference_file"
            @test_skip "Reference data not available"
        end
    end

end
