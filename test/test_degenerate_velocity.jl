# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Unit tests for the degeneracy-aware velocity-tensor construction used by the
# Wannier transport path: at a degenerate multiplet the band velocities are
# averaged over approach directions so the transport tensor does not depend on
# the arbitrary eigenbasis a LAPACK build picks inside the multiplet. These
# exercise the building blocks directly on synthetic input, independent of any
# interpolator.

using Test
using BoltzTraP
using LinearAlgebra
using Random

const _fib = BoltzTraP._fibonacci_directions
const _groups = BoltzTraP._degenerate_groups
const _degvv! = BoltzTraP._degenerate_vv!

# Random Hermitian L×L matrix.
function _rand_herm(rng, L)
    M = randn(rng, ComplexF64, L, L)
    return (M + M') / 2
end

# Run `_degenerate_vv!` on a single-k synthetic block of size L; return the
# (L, 3, 3) per-band velocity tensor for that block.
function _block_vv(Ablocks, L; ndirs = 200)
    A = zeros(ComplexF64, L, L, 3, 1)
    for a = 1:3
        A[:, :, a, 1] = Ablocks[a]
    end
    vvband = zeros(L, 3, 3, 1)
    _degvv!(vvband, A, collect(1:L), 1, _fib(ndirs))
    return vvband
end

@testset "degenerate velocity construction" begin
    @testset "Fibonacci directions are unit vectors" begin
        d = _fib(200)
        @test length(d) == 200
        @test all(x -> isapprox(sqrt(x[1]^2 + x[2]^2 + x[3]^2), 1.0; atol = 1e-12), d)
    end

    @testset "degenerate grouping" begin
        thr = 1e-6
        # two tight pairs (gaps below thr) + one isolated band
        groups = _groups([0.0, 5e-7, 1.0, 1.0 + 4e-7, 2.0], thr)
        @test sort(length.(groups)) == [1, 2, 2]
        # all gaps above thr → all singletons
        @test all(==(1), length.(_groups([0.0, 1.0, 2.0], thr)))
    end

    @testset "singleton reduces to v⊗v" begin
        v = (0.3, -0.7, 1.1)
        A = (
            fill(ComplexF64(v[1]), 1, 1),
            fill(ComplexF64(v[2]), 1, 1),
            fill(ComplexF64(v[3]), 1, 1),
        )
        vv = _block_vv(A, 1)
        for i = 1:3, j = 1:3
            @test vv[1, i, j, 1] ≈ v[i] * v[j] rtol = 1e-12
        end
    end

    @testset "gauge invariance under a block unitary" begin
        rng = MersenneTwister(20260530)
        for L in (2, 3)
            A = ntuple(_ -> _rand_herm(rng, L), 3)
            vv1 = _block_vv(A, L)
            Q = Matrix(qr(randn(rng, ComplexF64, L, L)).Q)   # unitary
            Arot = ntuple(a -> Q' * A[a] * Q, 3)
            vv2 = _block_vv(Arot, L)
            # the eigenvector rotation cancels exactly, so the per-band tensor
            # (not just the block total) is invariant
            @test vv1 ≈ vv2 rtol = 1e-9
        end
    end

    @testset "determinism" begin
        rng = MersenneTwister(7)
        A = ntuple(_ -> _rand_herm(rng, 3), 3)
        @test _block_vv(A, 3) == _block_vv(A, 3)
    end
end
