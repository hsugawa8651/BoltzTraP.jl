# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - SpinType tests for v0.3 magnetic support

@testset "SpinType" begin
    @testset "Type hierarchy" begin
        @test Unpolarized <: SpinType
        @test Collinear <: SpinType
        @test NonCollinear <: SpinType
    end

    @testset "dosweight" begin
        @test BoltzTraP.dosweight(Unpolarized) == 2.0
        @test BoltzTraP.dosweight(Collinear) == 1.0
        @test BoltzTraP.dosweight(NonCollinear) == 1.0
    end

    @testset "nspinor" begin
        @test BoltzTraP.nspinor(Unpolarized) == 1
        @test BoltzTraP.nspinor(Collinear) == 2
        @test BoltzTraP.nspinor(NonCollinear) == 1
    end
end
