# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

# Helper subtype for abstract-type fallback test (issue #54).
# Defined at top level because Julia struct declarations cannot be inside
# function/begin/let blocks.
struct _DummySubtypeInterp <: BoltzTraP.AbstractInterpolator end

@testset "AbstractBZSampling and interface hygiene" begin
    @testset "type hierarchy" begin
        @test isa(UniformMesh(), BoltzTraP.AbstractBZSampling)
        @test isa(UniformMesh((4, 4, 4)), BoltzTraP.AbstractBZSampling)
    end

    @testset "UniformMesh constructors" begin
        # Default constructor: nk unspecified (interpolator decides at runtime)
        @test UniformMesh().nk === nothing
        # Explicit nk: user-specified mesh
        @test UniformMesh((4, 4, 4)).nk == (4, 4, 4)
        @test UniformMesh((8, 8, 8)).nk == (8, 8, 8)
    end

    @testset "AbstractInterpolator abstract-type fallback (subtype without impl, #54)" begin
        # interpolate_bands: dispatch falls to abstract-type ::AbstractInterpolator method,
        # which throws ErrorException with a message naming the actual subtype
        # (rather than the misleading "expects a concrete subtype" from ::Any fallback).
        @test_throws ErrorException BoltzTraP.interpolate_bands(
            _DummySubtypeInterp(),
            zeros(0, 3),
        )

        # interpolate_velocities: same abstract-type behavior
        @test_throws ErrorException BoltzTraP.interpolate_velocities(
            _DummySubtypeInterp(),
            zeros(0, 3),
        )

        # interpolate(::AbstractInterpolator, kpoints) at src/interpolation.jl:64
        # is the existing default delegate that calls interpolate_bands and
        # interpolate_velocities. For _DummySubtypeInterp it cascades through the
        # new abstract-type fallback and throws ErrorException.
        @test_throws ErrorException BoltzTraP.interpolate(
            _DummySubtypeInterp(),
            zeros(0, 3),
        )
    end
end
