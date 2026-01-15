# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Tests for DFT format auto-detection and load_dft.

Verifies that:
- detected_format correctly identifies each DFT format
- load_dft auto-loads data without explicit format specification
"""

using Test
using BoltzTraP

const BOLTZTRAP_DATA = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data")

# Check if NCDatasets is available for ABINIT tests
const HAS_NCDATASETS = try
    using NCDatasets
    true
catch
    false
end

@testset "Format Detection" begin

    @testset "detected_format" begin
        @testset "VASP" begin
            vasp_dir = joinpath(BOLTZTRAP_DATA, "Si.vasp")
            if isdir(vasp_dir)
                @test BoltzTraP.detected_format(vasp_dir) == "VASP"
            else
                @test_skip "Test data not found: Si.vasp"
            end
        end

        @testset "QE" begin
            qe_dir = joinpath(BOLTZTRAP_DATA, "Si.ESPRESSO", "out")
            if isdir(qe_dir)
                @test BoltzTraP.detected_format(qe_dir) == "QE"
            else
                @test_skip "Test data not found: Si.ESPRESSO/out"
            end
        end

        @testset "Wien2k" begin
            wien2k_dir = joinpath(BOLTZTRAP_DATA, "Si")
            if isdir(wien2k_dir)
                @test BoltzTraP.detected_format(wien2k_dir) == "Wien2k"
            else
                @test_skip "Test data not found: Si"
            end
        end

        @testset "GENE" begin
            gene_dir = joinpath(BOLTZTRAP_DATA, "Si.GENE")
            if isdir(gene_dir)
                @test BoltzTraP.detected_format(gene_dir) == "GENE"
            else
                @test_skip "Test data not found: Si.GENE"
            end
        end

        @testset "ABINIT" begin
            abinit_dir = joinpath(BOLTZTRAP_DATA, "Si.abinit")
            if isdir(abinit_dir)
                @test BoltzTraP.detected_format(abinit_dir) == "ABINIT"
            else
                @test_skip "Test data not found: Si.abinit"
            end
        end

        @testset "Unknown directory" begin
            @test isnothing(BoltzTraP.detected_format("/nonexistent"))

            # Empty temp directory
            mktempdir() do tmpdir
                @test isnothing(BoltzTraP.detected_format(tmpdir))
            end
        end
    end

    @testset "load_dft auto-detection" begin
        @testset "VASP via load_dft" begin
            vasp_dir = joinpath(BOLTZTRAP_DATA, "Si.vasp")
            if isdir(vasp_dir)
                data = BoltzTraP.load_dft(vasp_dir)
                @test data isa DFTData
                @test "Si" in data.species
            else
                @test_skip "Test data not found: Si.vasp"
            end
        end

        @testset "QE via load_dft" begin
            qe_dir = joinpath(BOLTZTRAP_DATA, "Si.ESPRESSO", "out")
            if isdir(qe_dir)
                data = BoltzTraP.load_dft(qe_dir)
                @test data isa DFTData
                @test "Si" in data.species
            else
                @test_skip "Test data not found: Si.ESPRESSO/out"
            end
        end

        @testset "Wien2k via load_dft" begin
            wien2k_dir = joinpath(BOLTZTRAP_DATA, "Si")
            if isdir(wien2k_dir)
                data = BoltzTraP.load_dft(wien2k_dir)
                @test data isa DFTData
                @test "Si" in data.species
            else
                @test_skip "Test data not found: Si"
            end
        end

        @testset "GENE via load_dft" begin
            gene_dir = joinpath(BOLTZTRAP_DATA, "Si.GENE")
            if isdir(gene_dir)
                data = BoltzTraP.load_dft(gene_dir)
                @test data isa DFTData
                @test "Si" in data.species
            else
                @test_skip "Test data not found: Si.GENE"
            end
        end

        @testset "ABINIT via load_dft" begin
            abinit_dir = joinpath(BOLTZTRAP_DATA, "Si.abinit")
            if !isdir(abinit_dir)
                @test_skip "Test data not found: Si.abinit"
            elseif !HAS_NCDATASETS
                @test_skip "NCDatasets not available"
            else
                data = BoltzTraP.load_dft(abinit_dir)
                @test data isa DFTData
                @test "Si" in data.species
            end
        end

        @testset "Error on unknown format" begin
            mktempdir() do tmpdir
                @test_throws ErrorException BoltzTraP.load_dft(tmpdir)
            end
        end

        @testset "Error on nonexistent directory" begin
            @test_throws ErrorException BoltzTraP.load_dft("/nonexistent/path")
        end
    end

    @testset "load_dft returns same as explicit loader" begin
        @testset "VASP: load_dft == load_vasp" begin
            vasp_dir = joinpath(BOLTZTRAP_DATA, "Si.vasp")
            if isdir(vasp_dir)
                data_auto = BoltzTraP.load_dft(vasp_dir)
                data_explicit = load_vasp(vasp_dir)

                @test data_auto.lattice == data_explicit.lattice
                @test data_auto.kpoints == data_explicit.kpoints
                @test data_auto.ebands == data_explicit.ebands
                @test data_auto.fermi == data_explicit.fermi
            else
                @test_skip "Test data not found: Si.vasp"
            end
        end

        @testset "QE: load_dft == load_qe" begin
            qe_dir = joinpath(BOLTZTRAP_DATA, "Si.ESPRESSO", "out")
            if isdir(qe_dir)
                data_auto = BoltzTraP.load_dft(qe_dir)
                data_explicit = load_qe(qe_dir)

                @test data_auto.lattice == data_explicit.lattice
                @test data_auto.kpoints == data_explicit.kpoints
                @test data_auto.ebands == data_explicit.ebands
                @test data_auto.fermi == data_explicit.fermi
            else
                @test_skip "Test data not found: Si.ESPRESSO/out"
            end
        end
    end
end
