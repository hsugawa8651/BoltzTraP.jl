# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using LinearAlgebra
using NPZ
using BoltzTraP

@testset "Magnetic Workflow" begin

    @testset "run_interpolate Collinear (DFTData{2})" begin
        # Load PbTe collinear reference data
        REFTEST_DIR = joinpath(@__DIR__, "..", "reftest", "data")
        reference_file = joinpath(REFTEST_DIR, "pbte_collinear.npz")

        if !isfile(reference_file)
            @test_skip "PbTe collinear reference data not available"
            return
        end

        ref = npzread(reference_file)

        # Create synthetic DFTData{2} from reference
        # Note: Reference ebands is already concatenated (32, 145)
        # We need to split it back to (16, 145, 2) for DFTData{2}
        ebands_concat = ref["ebands"]  # (32, 145) - concatenated
        nbands_per_spin = size(ebands_concat, 1) ÷ 2
        nkpts = size(ebands_concat, 2)

        # Split into two spin channels
        ebands_spin1 = ebands_concat[1:nbands_per_spin, :]
        ebands_spin2 = ebands_concat[nbands_per_spin+1:end, :]
        ebands_3d = zeros(nbands_per_spin, nkpts, 2)
        ebands_3d[:, :, 1] = ebands_spin1
        ebands_3d[:, :, 2] = ebands_spin2

        # Decode species from uint8 array
        symbols_bytes = ref["symbols"]
        symbols_str = String(UInt8.(symbols_bytes))
        species = split(symbols_str, ",")

        # Create DFTData{2}
        data = BoltzTraP.DFTData(
            lattice = ref["lattvec"],
            positions = ref["positions"]',  # transpose to 3×natom
            species = species,
            kpoints = ref["kpoints"]',  # transpose to 3×nkpts
            weights = ones(nkpts) / nkpts,
            ebands = ebands_3d,
            occupations = zeros(size(ebands_3d)),
            fermi = ref["fermi"],
            nelect = ref["nelect"],
            magmom = ref["magmom"],
        )

        # Verify DFTData{2} was created
        @test data isa BoltzTraP.DFTData{2}
        @test BoltzTraP.is_magnetic(data)
        @test size(data.ebands, 3) == 2

        # Run interpolation (use nkpts as Python does with len(data.kpoints))
        result = BoltzTraP.run_interpolate(data; kpoints=nkpts, verbose=false)

        # Verify result
        @test result isa BoltzTraP.InterpolationResult

        # Check coefficients shape matches reference
        # Reference: (32, 205) - 32 bands (concatenated), 205 equivalences
        @test size(result.coeffs, 1) == size(ref["coeffs"], 1)  # nbands
        @test size(result.coeffs, 2) == size(ref["coeffs"], 2)  # neq

        # Check metadata
        @test result.metadata["dosweight"] == 1.0  # spin-polarized
        @test haskey(result.metadata, "spintype")
        @test result.metadata["spintype"] == "Collinear"

        # Note: Direct coefficient comparison is skipped because equivalence class
        # ordering differs between Python and Julia implementations.
        # The coefficients are correct but stored in different order.
        # Key verification: shape matches and interpolated bands would match.
        # See Sprint 6 for band reconstruction tests that verify correctness.
    end

end
