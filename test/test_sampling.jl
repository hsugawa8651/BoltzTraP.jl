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

@testset "TransportSystem" begin
    # Synthetic InterpolationResult helpers (avoid heavy run_interpolate path)
    function _synth_interp(;
        fermi::Float64 = 0.5,
        nelect::Float64 = 8.0,
        dosweight::Float64 = 2.0,
        spintype::String = "Unpolarized",
    )
        coeffs = ComplexF64[1.0+0im 2.0+0im; 3.0+0im 4.0+0im]
        equivs = [reshape([0, 0, 0], 1, 3), reshape([1, 0, 0], 1, 3)]
        lattvec = 5.0 * Matrix{Float64}(LinearAlgebra.I, 3, 3)
        metadata = Dict{String,Any}(
            "fermi" => fermi,
            "nelect" => nelect,
            "dosweight" => dosweight,
            "spintype" => spintype,
        )
        return BoltzTraP.InterpolationResult(coeffs, equivs, lattvec, nothing, metadata)
    end

    @testset "construction from InterpolationResult (non-magnetic Si-like)" begin
        ir = _synth_interp(;
            fermi = 0.5,
            nelect = 8.0,
            dosweight = 2.0,
            spintype = "Unpolarized",
        )
        sys = BoltzTraP.TransportSystem(ir)

        @test sys.lattice ≈ ir.lattvec
        @test sys.fermi == 0.5
        @test sys.nelect == 8.0
        @test sys.dosweight == 2.0
        @test sys.spintype isa Unpolarized
    end

    @testset "construction (collinear magnetic Fe-like)" begin
        ir = _synth_interp(;
            fermi = 0.3,
            nelect = 26.0,
            dosweight = 1.0,
            spintype = "Collinear",
        )
        sys = BoltzTraP.TransportSystem(ir)

        @test sys.fermi == 0.3
        @test sys.nelect == 26.0
        @test sys.dosweight == 1.0
        @test sys.spintype isa Collinear
    end

    @testset "unknown spintype throws" begin
        ir = _synth_interp(; spintype = "AlienSpin")
        @test_throws ErrorException BoltzTraP.TransportSystem(ir)
    end
end
