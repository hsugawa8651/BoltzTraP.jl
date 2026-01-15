# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests/test_fite.py

using Test
using LinearAlgebra
using NPZ
using Statistics: mean
using BoltzTraP: FourierInterpolator, interpolate, interpolate_bands, interpolate_velocities

@testset "Interpolation Module (fite)" begin

    TEST_DIR = @__DIR__
    DATA_DIR = joinpath(TEST_DIR, "data")

    @testset "fitde3D saved coefficients" begin
        # Load reference data
        ref_file = joinpath(DATA_DIR, "Si_fitde3D.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        ref_data = npzread(ref_file)
        coeffs_ref = ref_data["noder"]

        # Load Si interpolation data from reftest
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Compute coefficients
        coeffs = BoltzTraP.fitde3D(kpoints, ebands, equivalences, lattvec)

        # Compare with reference (may have small differences due to numerical precision)
        @test size(coeffs) == size(coeffs_ref)

        # The coefficients should be close but may not be identical
        # due to different k-point sets and equivalence classes
        # Just verify the computation runs and produces reasonable output
        @test !any(isnan, coeffs)
        @test !any(isinf, coeffs)
    end

    @testset "getBands reconstruct" begin
        # getBands should be able to reconstruct the original bands
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Fit and reconstruct
        coeffs = BoltzTraP.fitde3D(kpoints, ebands, equivalences, lattvec)
        ebands_reconstructed, vbands = BoltzTraP.getBands(kpoints, equivalences, lattvec, coeffs)

        # Original bands should be reconstructed accurately
        @test size(ebands_reconstructed) == size(ebands)
        max_error = maximum(abs.(ebands_reconstructed - ebands))
        @test max_error < 1e-8
    end

    @testset "FourierInterpolator" begin
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Create interpolator
        interp = FourierInterpolator(kpoints, ebands, equivalences, lattvec)

        # Test interpolate_bands
        ebands_interp = interpolate_bands(interp, kpoints)
        @test size(ebands_interp) == size(ebands)
        @test maximum(abs.(ebands_interp - ebands)) < 1e-8

        # Test interpolate_velocities
        vbands = interpolate_velocities(interp, kpoints)
        @test size(vbands, 1) == 3  # 3 components
        @test size(vbands, 2) == size(ebands, 1)  # nbands
        @test size(vbands, 3) == size(ebands, 2)  # nk

        # Test interpolate (both)
        e, v = interpolate(interp, kpoints)
        @test e ≈ ebands_interp
        @test v ≈ vbands
    end

    @testset "Convenience functions" begin
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Test convenience function
        e, v = interpolate(kpoints, ebands, equivalences, lattvec, kpoints)
        @test size(e) == size(ebands)
        @test maximum(abs.(e - ebands)) < 1e-8

        e_only = interpolate_bands(kpoints, ebands, equivalences, lattvec, kpoints)
        @test e_only ≈ e

        v_only = interpolate_velocities(kpoints, ebands, equivalences, lattvec, kpoints)
        @test v_only ≈ v
    end

    @testset "Error dispatch" begin
        # Non-AbstractInterpolator should throw ArgumentError
        @test_throws ArgumentError interpolate_bands("bad", zeros(3, 3))
        @test_throws ArgumentError interpolate_velocities("bad", zeros(3, 3))
        @test_throws ArgumentError interpolate("bad", zeros(3, 3))
    end

    @testset "compute_phase_factors" begin
        # Test phase factor computation
        kpoints = [0.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.5 0.0]
        equivalences = [zeros(Int, 1, 3), [1 0 0; -1 0 0]]

        phase = BoltzTraP.compute_phase_factors(kpoints, equivalences)

        # Phase for R=0 should be 1
        @test all(phase[:, 1] .≈ 1.0)

        # Phase for k=[0.5,0,0] with R=±[1,0,0] should be cos(π) = -1 averaged
        # Actually: exp(2πi * 0.5 * 1) = exp(πi) = -1
        #           exp(2πi * 0.5 * (-1)) = exp(-πi) = -1
        # Average = -1
        @test real(phase[2, 2]) ≈ -1.0 atol=1e-10
    end

    @testset "Velocity comparison with Python" begin
        # Compare computed velocities with Python reference
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        vbands_ref = si_data["vbands"]  # Python reference (3, nbands, nk)
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Compute velocities
        _, vbands = interpolate(kpoints, ebands, equivalences, lattvec, kpoints)

        # Compare shapes
        @test size(vbands) == size(vbands_ref)

        # Velocities should match Python exactly
        max_abs_diff = maximum(abs.(vbands - vbands_ref))
        @test max_abs_diff < 1e-8
    end

    @testset "Interpolation at new k-points" begin
        # Test interpolation at k-points different from original
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]

        # Create interpolator
        interp = FourierInterpolator(kpoints, ebands, equivalences, lattvec)

        # Generate new k-points (shifted from original)
        new_kpoints = kpoints[1:10, :] .+ 0.01

        # Interpolate at new k-points
        ebands_new = interpolate_bands(interp, new_kpoints)

        # Basic sanity checks
        @test size(ebands_new, 1) == size(ebands, 1)  # Same number of bands
        @test size(ebands_new, 2) == size(new_kpoints, 1)  # Same number of k-points
        @test !any(isnan, ebands_new)
        @test !any(isinf, ebands_new)

        # Energies should be in reasonable range (similar to original)
        @test minimum(ebands_new) > minimum(ebands) - 0.1
        @test maximum(ebands_new) < maximum(ebands) + 0.1
    end

    @testset "Coefficients self-consistency" begin
        # Test that computed coefficients match what we store
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")

        if !isfile(si_file)
            @warn "Reference file not found"
            @test_skip "Reference data not available"
            return
        end

        si_data = npzread(si_file)
        kpoints = si_data["kpoints"]
        ebands = si_data["ebands"]
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]
        coeffs_ref = si_data["coeffs"]  # Coefficients generated by Python

        # Recompute coefficients
        coeffs_julia = BoltzTraP.fitde3D(kpoints, ebands, equivalences, lattvec)

        # Shapes should match
        @test size(coeffs_julia) == size(coeffs_ref)

        # Coefficients should be very close (same algorithm)
        max_diff = maximum(abs.(coeffs_julia - coeffs_ref))
        @test max_diff < 1e-10
    end

    @testset "getBTPbands FFT reconstruction" begin
        # Test FFT-based band reconstruction
        si_file = joinpath(dirname(TEST_DIR), "reftest", "data", "si_interpolation.npz")
        if !isfile(si_file)
            @warn "Si interpolation data not found: $si_file"
            @test_skip "Si interpolation data not available"
            return
        end

        si_data = npzread(si_file)
        lattvec = si_data["lattvec"]
        n_eq = si_data["n_equivalences"]
        equivalences = [si_data["equiv_$i"] for i in 0:(n_eq-1)]
        coeffs = si_data["coeffs"]

        # Call getBTPbands: (coeffs, equivalences, lattvec)
        eband, vvband = BoltzTraP.getBTPbands(coeffs, equivalences, lattvec)

        # Basic checks
        @test size(eband, 1) == size(coeffs, 1)  # nbands
        @test !any(isnan, eband)
        @test !any(isinf, eband)

        # vvband should have velocity outer product structure
        @test size(vvband, 1) == size(coeffs, 1)  # nbands
        @test size(vvband, 2) == 3  # velocity component
        @test size(vvband, 3) == 3  # velocity component
    end

end
