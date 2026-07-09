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

@testset "FourierInterpolator from InterpolationResult" begin
    coeffs = ComplexF64[1.0+0im 2.0+0im; 3.0+0im 4.0+0im]
    equivs = [reshape([0, 0, 0], 1, 3), reshape([1, 0, 0], 1, 3)]
    lattvec = 5.0 * Matrix{Float64}(LinearAlgebra.I, 3, 3)
    ir = BoltzTraP.InterpolationResult(coeffs, equivs, lattvec, nothing, Dict{String,Any}())

    fi_manual = BoltzTraP.FourierInterpolator(
        ir.coeffs,
        ir.equivalences,
        SMatrix{3,3,Float64}(ir.lattvec),
    )
    fi_conv = BoltzTraP.FourierInterpolator(ir)

    @test fi_conv isa BoltzTraP.FourierInterpolator
    @test fi_conv.coeffs == fi_manual.coeffs
    @test fi_conv.equivalences == fi_manual.equivalences
    @test fi_conv.lattvec == fi_manual.lattvec
end

# Helper: dummy AbstractBZSampling subtype for abstract-type fallback test
struct _DummyBZSampling <: BoltzTraP.AbstractBZSampling end

@testset "solve_transport dispatch" begin
    @testset "rejects wrong argument types (::Any fallback)" begin
        # Build minimal valid interp + sys for the "other two args correct" cases
        ir = let
            coeffs = ComplexF64[1.0+0im 2.0+0im; 3.0+0im 4.0+0im]
            equivs = [reshape([0, 0, 0], 1, 3), reshape([1, 0, 0], 1, 3)]
            lattvec = 5.0 * Matrix{Float64}(LinearAlgebra.I, 3, 3)
            metadata = Dict{String,Any}(
                "fermi" => 0.5,
                "nelect" => 8.0,
                "dosweight" => 2.0,
                "spintype" => "Unpolarized",
            )
            BoltzTraP.InterpolationResult(coeffs, equivs, lattvec, nothing, metadata)
        end
        fi = BoltzTraP.FourierInterpolator(
            ir.coeffs,
            ir.equivalences,
            SMatrix{3,3,Float64}(ir.lattvec),
        )
        sys = BoltzTraP.TransportSystem(ir)

        # Wrong sampling type
        @test_throws ArgumentError solve_transport("not_sampling", fi, sys)
        # Wrong interpolator type
        @test_throws ArgumentError solve_transport(UniformMesh(), "not_interp", sys)
        # Wrong sys type
        @test_throws ArgumentError solve_transport(UniformMesh(), fi, "not_sys")
    end

    @testset "rejects subtype without concrete method (::AbstractBZSampling fallback)" begin
        ir = let
            coeffs = ComplexF64[1.0+0im 2.0+0im; 3.0+0im 4.0+0im]
            equivs = [reshape([0, 0, 0], 1, 3), reshape([1, 0, 0], 1, 3)]
            lattvec = 5.0 * Matrix{Float64}(LinearAlgebra.I, 3, 3)
            metadata = Dict{String,Any}(
                "fermi" => 0.5,
                "nelect" => 8.0,
                "dosweight" => 2.0,
                "spintype" => "Unpolarized",
            )
            BoltzTraP.InterpolationResult(coeffs, equivs, lattvec, nothing, metadata)
        end
        fi = BoltzTraP.FourierInterpolator(
            ir.coeffs,
            ir.equivalences,
            SMatrix{3,3,Float64}(ir.lattvec),
        )
        sys = BoltzTraP.TransportSystem(ir)

        # _DummyBZSampling has no concrete solve_transport method, falls to
        # the abstract-type fallback that throws ErrorException.
        @test_throws ErrorException solve_transport(_DummyBZSampling(), fi, sys)
    end

    @testset "concrete UniformMesh + FourierInterpolator on Si VASP" begin
        datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
        if isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            sys = BoltzTraP.TransportSystem(ir)
            fi = BoltzTraP.FourierInterpolator(
                ir.coeffs,
                ir.equivalences,
                SMatrix{3,3,Float64}(ir.lattvec),
            )
            result = solve_transport(UniformMesh(), fi, sys; temperatures = [300.0])

            @test result isa TransportResult
            @test size(result.sigma, 1) == 3
            @test size(result.sigma, 2) == 3
            @test size(result.sigma, 3) == 1   # 1 temperature
            @test size(result.sigma, 4) == length(result.mu_values)
            @test !any(isnan, result.sigma)
            @test !any(isnan, result.seebeck)
            @test !any(isnan, result.kappa)
        else
            @warn "Skipping solve_transport Si VASP integration test: data not found at $datadir"
        end
    end
end

@testset "run_integrate (backward-compat wrapper)" begin
    @testset "auto μ grid wrapper preserves source metadata" begin
        datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
        if isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            result = BoltzTraP.run_integrate(ir; temperatures = [300.0])

            @test result isa TransportResult
            @test size(result.sigma, 3) == 1
            @test !any(isnan, result.sigma)
            # source must come from IR metadata (wrapper restores), not the
            # solve_transport delegate's hardcoded "solve_transport".
            @test result.metadata["source"] != "solve_transport"
            # ... and it must be the value the interpolation stage recorded,
            # not the "unknown" fallback.
            @test result.metadata["source"] == ir.metadata["source"]
            @test haskey(result.metadata, "fermi_Ha")
            @test haskey(result.metadata, "vuc_ang3")
        else
            @warn "Skipping run_integrate auto-μ wrapper test: data not found at $datadir"
        end
    end

    @testset "explicit μ grid wrapper (positional mur)" begin
        datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
        if isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            mur_test = collect(0.0:0.05:0.5)
            result = BoltzTraP.run_integrate(ir, mur_test; temperatures = [300.0])

            @test result isa TransportResult
            @test length(result.mu_values) == length(mur_test)
            @test !any(isnan, result.sigma)
            @test result.metadata["source"] != "solve_transport"
            @test result.metadata["source"] == ir.metadata["source"]
        else
            @warn "Skipping run_integrate explicit-μ wrapper test: data not found at $datadir"
        end
    end
end

@testset "run_integrate (interpolator + system)" begin
    datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
    have_data = isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))

    @testset "canonical form reproduces the wrapper" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            wrapped = BoltzTraP.run_integrate(ir; temperatures = [300.0])
            canonical = BoltzTraP.run_integrate(
                BoltzTraP.FourierInterpolator(ir),
                BoltzTraP.TransportSystem(ir);
                temperatures = [300.0],
            )

            @test canonical.sigma == wrapped.sigma
            @test canonical.seebeck == wrapped.seebeck
            @test canonical.kappa == wrapped.kappa
            @test canonical.mu_values == wrapped.mu_values

            # Only the wrapper knows the interpolation metadata, so only the
            # wrapper can restore the source. The canonical form leaves the
            # value stamped by solve_transport.
            @test wrapped.metadata["source"] == ir.metadata["source"]
            @test canonical.metadata["source"] == "solve_transport"
        else
            @warn "Skipping canonical run_integrate test: data not found at $datadir"
        end
    end

    @testset "canonical form with explicit μ grid" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            mur_test = collect(0.0:0.05:0.5)
            wrapped = BoltzTraP.run_integrate(ir, mur_test; temperatures = [300.0])
            canonical = BoltzTraP.run_integrate(
                BoltzTraP.FourierInterpolator(ir),
                BoltzTraP.TransportSystem(ir),
                mur_test;
                temperatures = [300.0],
            )

            @test canonical.sigma == wrapped.sigma
            @test canonical.mu_values == wrapped.mu_values
        else
            @warn "Skipping canonical run_integrate mur test: data not found at $datadir"
        end
    end

    @testset "sampling keyword is forwarded to solve_transport" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            fi = BoltzTraP.FourierInterpolator(ir)
            sys = BoltzTraP.TransportSystem(ir)
            default_mesh = BoltzTraP.run_integrate(fi, sys; temperatures = [300.0])
            explicit_mesh = BoltzTraP.run_integrate(
                fi,
                sys;
                sampling = UniformMesh(),
                temperatures = [300.0],
            )
            @test explicit_mesh.sigma == default_mesh.sigma
            @test default_mesh.metadata["sampling"] == "UniformMesh"
        else
            @warn "Skipping sampling keyword test: data not found at $datadir"
        end
    end
end

@testset "integration grid on the Fourier path" begin
    datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
    have_data = isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))

    @testset "an explicit mesh size warns and is ignored" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            fi = BoltzTraP.FourierInterpolator(ir)
            sys = BoltzTraP.TransportSystem(ir)

            derived = solve_transport(UniformMesh(), fi, sys; temperatures = [300.0])
            requested = @test_logs (:warn,) match_mode = :any solve_transport(
                UniformMesh((8, 8, 8)),
                fi,
                sys;
                temperatures = [300.0],
            )
            # The warning is the only effect: the grid comes from the basis.
            @test requested.sigma == derived.sigma
            @test requested.metadata["sampling_nk"] == derived.metadata["sampling_nk"]

            # The default mesh must stay silent.
            @test_logs solve_transport(UniformMesh(), fi, sys; temperatures = [300.0])
        else
            @warn "Skipping explicit mesh warning test: data not found at $datadir"
        end
    end

    @testset "sampling_nk records the grid that was used" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            dims = BoltzTraP.determine_fft_dims(ir.equivalences)
            result = BoltzTraP.run_integrate(ir; temperatures = [300.0])

            @test result.metadata["sampling_nk"] == dims
            @test prod(result.metadata["sampling_nk"]) == result.metadata["npts_fft"]

            mktempdir() do dir
                path = joinpath(dir, "sampling_nk_roundtrip.jld2")
                save_integrate(path, result)
                @test load_integrate(path).metadata["sampling_nk"] == dims
            end
        else
            @warn "Skipping sampling_nk test: data not found at $datadir"
        end
    end
end

@testset "provenance metadata" begin
    datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")
    have_data = isdir(datadir) && isfile(joinpath(datadir, "vasprun.xml"))

    @testset "solve_transport stamps sampling/interpolator + preserves legacy keys" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            sys = BoltzTraP.TransportSystem(ir)
            fi = BoltzTraP.FourierInterpolator(
                ir.coeffs,
                ir.equivalences,
                SMatrix{3,3,Float64}(ir.lattvec),
            )
            result = solve_transport(UniformMesh(), fi, sys; temperatures = [300.0])

            # New provenance keys.
            @test result.metadata["sampling"] == "UniformMesh"
            @test result.metadata["interpolator"] == "Fourier"

            # Legacy metadata keys must remain (regression: stamping is additive).
            for key in ("nelect", "dosweight", "fermi_Ha", "vuc_ang3", "nbands", "npts_dos")
                @test haskey(result.metadata, key)
            end
        else
            @warn "Skipping provenance metadata (solve_transport): data not found at $datadir"
        end
    end

    @testset "run_integrate wrapper preserves provenance + source override" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            result = BoltzTraP.run_integrate(ir; temperatures = [300.0])

            # Wrapper passes through new keys from solve_transport.
            @test result.metadata["sampling"] == "UniformMesh"
            @test result.metadata["interpolator"] == "Fourier"

            # Wrapper's source override does not clobber the provenance keys.
            @test result.metadata["source"] != "solve_transport"
        else
            @warn "Skipping provenance metadata (wrapper): data not found at $datadir"
        end
    end

    @testset "JLD2 round-trip preserves provenance keys" begin
        if have_data
            ir = BoltzTraP.run_interpolate(datadir; kpoints = 200, verbose = false)
            result = BoltzTraP.run_integrate(ir; temperatures = [300.0])

            mktempdir() do dir
                path = joinpath(dir, "provenance_roundtrip.jld2")
                save_integrate(path, result)
                loaded = load_integrate(path)
                @test loaded.metadata["sampling"] == "UniformMesh"
                @test loaded.metadata["interpolator"] == "Fourier"
            end
        else
            @warn "Skipping provenance metadata (JLD2 round-trip): data not found at $datadir"
        end
    end
end
