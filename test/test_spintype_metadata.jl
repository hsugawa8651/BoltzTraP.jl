# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP
using BoltzTraP: FourierInterpolator
using StaticArrays

@testset "spintype metadata" begin

    # Create minimal test data
    coeffs = rand(ComplexF64, 4, 100)
    equivalences = [rand(-5:5, 10, 3) for _ in 1:100]
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

        result = TransportResult(temperatures, mu_values, sigma, seebeck, kappa, nothing, metadata)

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
        metadata = Dict{String,Any}(
            "source" => "test",
            "spintype" => "Unpolarized",
        )

        result = TransportResult(temperatures, mu_values, sigma, seebeck, kappa, nothing, metadata)

        mktempdir() do dir
            path = joinpath(dir, "test_transport.jld2")
            save_integrate(path, result)
            loaded = load_integrate(path)

            @test haskey(loaded.metadata, "spintype")
            @test loaded.metadata["spintype"] == "Unpolarized"
        end
    end

end
