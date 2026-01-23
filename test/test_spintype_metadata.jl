# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP
using BoltzTraP: FourierInterpolator
using StaticArrays
using JLD2

@testset "spintype metadata" begin

    # Create minimal test data
    coeffs = rand(ComplexF64, 4, 100)
    equivalences = [rand(-5:5, 10, 3) for _ = 1:100]
    lattvec = SMatrix{3,3}(rand(3, 3))

    @testset "InterpolationResult contains spintype" begin
        interp = FourierInterpolator(coeffs, equivalences, lattvec)
        result = InterpolationResult(interp)

        @test haskey(result.metadata, "spintype")
        @test result.metadata["spintype"] == "Unpolarized"
    end

    @testset "InterpolationResult JLD2 round-trip" begin
        interp = FourierInterpolator(coeffs, equivalences, lattvec)
        result = InterpolationResult(interp)

        mktempdir() do dir
            path = joinpath(dir, "test_interp.jld2")
            save_interpolation(path, result)
            loaded = load_interpolation(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

    @testset "InterpolationResult BT2 round-trip" begin
        interp = FourierInterpolator(coeffs, equivalences, lattvec)
        result = InterpolationResult(interp)

        mktempdir() do dir
            path = joinpath(dir, "test_interp.bt2")
            save_interpolation(path, result)
            loaded = load_interpolation(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

    @testset "TransportResult contains spintype" begin
        # Create minimal TransportResult
        temperatures = [300.0]
        mu_values = collect(-1.0:0.1:1.0)
        nT, nμ = length(temperatures), length(mu_values)

        sigma = zeros(3, 3, nT, nμ)
        seebeck = zeros(3, 3, nT, nμ)
        kappa = zeros(3, 3, nT, nμ)
        metadata = Dict{String,Any}(
            "source" => "test",
            "nelect" => 8.0,
            "dosweight" => 2.0,
            "spintype" => "Unpolarized",
        )

        result = TransportResult(
            temperatures,
            mu_values,
            sigma,
            seebeck,
            kappa,
            nothing,
            metadata,
        )

        @test haskey(result.metadata, "spintype")
        @test result.metadata["spintype"] == "Unpolarized"
    end

    @testset "TransportResult JLD2 round-trip" begin
        temperatures = [300.0]
        mu_values = collect(-1.0:0.1:1.0)
        nT, nμ = length(temperatures), length(mu_values)

        sigma = zeros(3, 3, nT, nμ)
        seebeck = zeros(3, 3, nT, nμ)
        kappa = zeros(3, 3, nT, nμ)
        metadata = Dict{String,Any}("source" => "test", "spintype" => "Unpolarized")

        result = TransportResult(
            temperatures,
            mu_values,
            sigma,
            seebeck,
            kappa,
            nothing,
            metadata,
        )

        mktempdir() do dir
            path = joinpath(dir, "test_transport.jld2")
            save_integrate(path, result)
            loaded = load_integrate(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

    # ========================================================================
    # v0.1 Backward Compatibility Tests (Sprint 7)
    # ========================================================================

    @testset "v0.1 backward compat: InterpolationResult JLD2 without spintype" begin
        # Simulate v0.1 file by saving without spintype in metadata
        mktempdir() do dir
            path = joinpath(dir, "v01_interp.jld2")

            # Create v0.1-style data (no spintype in metadata)
            coeffs_v01 = rand(ComplexF64, 4, 100)
            equivs_v01 = [rand(-5:5, 10, 3) for _ in 1:100]
            lattvec_v01 = rand(3, 3)
            metadata_v01 = Dict{String,Any}(
                "version" => "0.1.0",
                "created" => "2026-01-01",
                "nbands" => 4,
                "neq" => 100,
                # Note: no "spintype" key - simulating v0.1 file
            )

            # Save directly using JLD2 (bypassing our save function that adds spintype)
            jldsave(path;
                coeffs=coeffs_v01,
                equivalences=equivs_v01,
                lattvec=lattvec_v01,
                atoms=nothing,
                metadata=metadata_v01
            )

            # Load using our function - should add spintype default
            loaded = load_interpolation(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

    @testset "v0.1 backward compat: TransportResult JLD2 without spintype" begin
        mktempdir() do dir
            path = joinpath(dir, "v01_transport.jld2")

            # Create v0.1-style data (no spintype in metadata)
            temperatures = [300.0]
            mu_values = collect(-1.0:0.1:1.0)
            nT, nμ = length(temperatures), length(mu_values)

            sigma_v01 = zeros(3, 3, nT, nμ)
            seebeck_v01 = zeros(3, 3, nT, nμ)
            kappa_v01 = zeros(3, 3, nT, nμ)
            metadata_v01 = Dict{String,Any}(
                "source" => "test",
                # Note: no "spintype" key - simulating v0.1 file
            )

            # Save directly using JLD2
            jldsave(path;
                temperatures=temperatures,
                mu_values=mu_values,
                sigma=sigma_v01,
                seebeck=seebeck_v01,
                kappa=kappa_v01,
                dos=nothing,
                metadata=metadata_v01
            )

            # Load using our function - should add spintype default
            loaded = load_integrate(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

    # ========================================================================
    # Collinear Spintype Tests (Sprint 7)
    # ========================================================================

    @testset "Collinear spintype: InterpolationResult JLD2 round-trip" begin
        interp = FourierInterpolator(coeffs, equivalences, lattvec)
        # Create result with Collinear spintype
        metadata_collinear = Dict{String,Any}(
            "spintype" => "Collinear",
            "dosweight" => 1.0,  # spin-polarized
        )
        result = InterpolationResult(interp; metadata=metadata_collinear)

        @test result.metadata["spintype"] == "Collinear"

        mktempdir() do dir
            path = joinpath(dir, "collinear_interp.jld2")
            save_interpolation(path, result)
            loaded = load_interpolation(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Collinear"
            @test loaded.metadata["dosweight"] == 1.0
        end
    end

    @testset "Collinear spintype: InterpolationResult BT2 round-trip" begin
        interp = FourierInterpolator(coeffs, equivalences, lattvec)
        metadata_collinear = Dict{String,Any}(
            "spintype" => "Collinear",
            "dosweight" => 1.0,
        )
        result = InterpolationResult(interp; metadata=metadata_collinear)

        mktempdir() do dir
            path = joinpath(dir, "collinear_interp.bt2")
            save_interpolation(path, result)
            loaded = load_interpolation(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Collinear"
            @test loaded.metadata["dosweight"] == 1.0
        end
    end

    @testset "Collinear spintype: TransportResult JLD2 round-trip" begin
        temperatures = [300.0]
        mu_values = collect(-1.0:0.1:1.0)
        nT, nμ = length(temperatures), length(mu_values)

        sigma = zeros(3, 3, nT, nμ)
        seebeck = zeros(3, 3, nT, nμ)
        kappa = zeros(3, 3, nT, nμ)
        metadata = Dict{String,Any}(
            "source" => "test",
            "spintype" => "Collinear",
            "dosweight" => 1.0,
        )

        result = TransportResult(temperatures, mu_values, sigma, seebeck, kappa, nothing, metadata)

        mktempdir() do dir
            path = joinpath(dir, "collinear_transport.jld2")
            save_integrate(path, result)
            loaded = load_integrate(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Collinear"
            @test loaded.metadata["dosweight"] == 1.0
        end
    end

end
