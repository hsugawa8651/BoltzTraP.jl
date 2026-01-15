# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Functional test: DFT format auto-detection
#
# Usage:
#   julia --project=. ftest/test_format_detection.jl
#
# This script tests that DFT format auto-detection works correctly.
# Run with JULIA_DEBUG=BoltzTraP to see debug output.

using Logging

# Enable debug logging
ENV["JULIA_DEBUG"] = "BoltzTraP"
global_logger(ConsoleLogger(stderr, Logging.Debug))

using BoltzTraP

# Test data paths (relative to project root)
const BOLTZTRAP2_DATA = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data")

function main()
    println("=" ^ 60)
    println("Functional Test: DFT Format Auto-Detection")
    println("=" ^ 60)
    println()

    # Test 1: Detect VASP format
    println("-" ^ 60)
    println("Test 1: Detect VASP format")
    println("-" ^ 60)
    vasp_dir = joinpath(BOLTZTRAP2_DATA, "Si.vasp")
    if isdir(vasp_dir)
        detected = BoltzTraP.detected_format(vasp_dir)
        println("  Directory: $vasp_dir")
        println("  Detected: $detected")
        @assert detected == "VASP" "Expected VASP, got $detected"
        println("  ✓ PASS")
    else
        println("  ⚠ SKIP (test data not found: $vasp_dir)")
    end
    println()

    # Test 2: Load VASP with auto-detection
    println("-" ^ 60)
    println("Test 2: Load VASP with auto-detection (load_dft)")
    println("-" ^ 60)
    if isdir(vasp_dir)
        @info "Calling load_dft for VASP"
        data = BoltzTraP.load_dft(vasp_dir)
        println("  Loaded successfully")
        println("  → kpoints: $(size(data.kpoints, 2))")
        println("  → bands: $(size(data.ebands, 1))")
        println("  → fermi: $(data.fermi) Ha")
        @assert size(data.kpoints, 2) > 0 "Expected kpoints > 0"
        @assert size(data.ebands, 1) > 0 "Expected bands > 0"
        println("  ✓ PASS")
    else
        println("  ⚠ SKIP")
    end
    println()

    # Test 3: Load with explicit format (vasp)
    println("-" ^ 60)
    println("Test 3: Load with explicit format (load_vasp)")
    println("-" ^ 60)
    if isdir(vasp_dir)
        @info "Calling load_vasp explicitly"
        data = BoltzTraP.load_vasp(vasp_dir)
        println("  Loaded successfully")
        println("  → kpoints: $(size(data.kpoints, 2))")
        println("  ✓ PASS")
    else
        println("  ⚠ SKIP")
    end
    println()

    # Test 4: Non-existent directory
    println("-" ^ 60)
    println("Test 4: Non-existent directory")
    println("-" ^ 60)
    fake_dir = "/nonexistent/directory"
    detected = BoltzTraP.detected_format(fake_dir)
    println("  Directory: $fake_dir")
    println("  Detected: $detected")
    @assert isnothing(detected) "Expected nothing for non-existent directory"
    println("  ✓ PASS")
    println()

    # Test 5: Empty directory (no format detected)
    println("-" ^ 60)
    println("Test 5: Empty directory")
    println("-" ^ 60)
    empty_dir = mktempdir()
    detected = BoltzTraP.detected_format(empty_dir)
    println("  Directory: $empty_dir")
    println("  Detected: $detected")
    @assert isnothing(detected) "Expected nothing for empty directory"
    rm(empty_dir)
    println("  ✓ PASS")
    println()

    # Test 6: load_dft error on unknown format
    println("-" ^ 60)
    println("Test 6: load_dft error on unknown format")
    println("-" ^ 60)
    unknown_dir = mktempdir()
    touch(joinpath(unknown_dir, "random.txt"))
    try
        BoltzTraP.load_dft(unknown_dir)
        @assert false "Expected error for unknown format"
    catch e
        println("  Directory: $unknown_dir")
        println("  Error (expected): $(typeof(e))")
        println("  ✓ PASS")
    finally
        rm(unknown_dir; recursive=true)
    end
    println()

    println("=" ^ 60)
    println("All tests passed!")
    println("=" ^ 60)
end

main()
