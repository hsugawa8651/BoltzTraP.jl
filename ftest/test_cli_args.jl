# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Functional test: CLI argument propagation
#
# Usage:
#   julia --project=. ftest/test_cli_args.jl [vasp_directory]
#
# This script tests that CLI options are correctly passed to workflow functions.
# Run with JULIA_DEBUG=BoltzTraP to see debug output.

using Logging

# Enable debug logging for BoltzTraP module
ENV["JULIA_DEBUG"] = "BoltzTraP"
global_logger(ConsoleLogger(stderr, Logging.Debug))

using BoltzTraP

# Default test data path
const DEFAULT_VASP_DIR = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data", "Si.vasp")

function main()
    vasp_dir = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_VASP_DIR

    if !isdir(vasp_dir)
        @error "VASP directory not found" vasp_dir
        println("\nUsage: julia --project=. ftest/test_cli_args.jl [vasp_directory]")
        println("\nExpected directory structure:")
        println("  vasp_directory/")
        println("    vasprun.xml")
        exit(1)
    end

    println("=" ^ 60)
    println("Functional Test: CLI Argument Propagation")
    println("=" ^ 60)
    println()

    # Test 1: Default arguments
    println("-" ^ 60)
    println("Test 1: Default arguments (kpoints=5000)")
    println("-" ^ 60)
    @info "Calling run_interpolate with defaults"
    result1 = BoltzTraP.run_interpolate(vasp_dir; verbose=true)
    println("  → nkpt_target in metadata: $(result1.metadata["nkpt_target"])")
    @assert result1.metadata["nkpt_target"] == 5000 "Expected nkpt_target=5000"
    println("  ✓ PASS")
    println()

    # Test 2: Custom kpoints
    println("-" ^ 60)
    println("Test 2: Custom kpoints=1000")
    println("-" ^ 60)
    @info "Calling run_interpolate with kpoints=1000"
    result2 = BoltzTraP.run_interpolate(vasp_dir; kpoints=1000, verbose=true)
    println("  → nkpt_target in metadata: $(result2.metadata["nkpt_target"])")
    @assert result2.metadata["nkpt_target"] == 1000 "Expected nkpt_target=1000"
    println("  ✓ PASS")
    println()

    # Test 3: Energy range
    println("-" ^ 60)
    println("Test 3: Energy range emin=-0.5, emax=0.5")
    println("-" ^ 60)
    @info "Calling run_interpolate with energy range"
    result3 = BoltzTraP.run_interpolate(vasp_dir; kpoints=500, emin=-0.5, emax=0.5, verbose=true)
    println("  → emin in metadata: $(result3.metadata["emin"])")
    println("  → emax in metadata: $(result3.metadata["emax"])")
    @assert result3.metadata["emin"] == -0.5 "Expected emin=-0.5"
    @assert result3.metadata["emax"] == 0.5 "Expected emax=0.5"
    println("  ✓ PASS")
    println()

    # Test 4: Multiplier instead of kpoints
    println("-" ^ 60)
    println("Test 4: Multiplier=2")
    println("-" ^ 60)
    @info "Calling run_interpolate with multiplier=2"
    result4 = BoltzTraP.run_interpolate(vasp_dir; multiplier=2, verbose=true)
    expected_nkpt = result4.metadata["nkpt_original"] * 2
    println("  → nkpt_original: $(result4.metadata["nkpt_original"])")
    println("  → nkpt_target in metadata: $(result4.metadata["nkpt_target"])")
    @assert result4.metadata["nkpt_target"] == expected_nkpt "Expected nkpt_target=$expected_nkpt"
    println("  ✓ PASS")
    println()

    # Test 5: Absolute energy mode
    println("-" ^ 60)
    println("Test 5: Absolute energy mode")
    println("-" ^ 60)
    @info "Calling run_interpolate with absolute=true"
    result5 = BoltzTraP.run_interpolate(vasp_dir; kpoints=500, emin=-1.0, emax=1.0, absolute=true, verbose=true)
    println("  → absolute in metadata: $(result5.metadata["absolute"])")
    @assert result5.metadata["absolute"] == true "Expected absolute=true"
    println("  ✓ PASS")
    println()

    # Test 6: run_integrate with temperatures
    println("-" ^ 60)
    println("Test 6: run_integrate with temperatures")
    println("-" ^ 60)
    @info "Calling run_integrate with temperatures=[100, 200, 300]"
    transport = BoltzTraP.run_integrate(result1; temperatures=[100.0, 200.0, 300.0], verbose=true)
    println("  → temperatures: $(transport.temperatures)")
    @assert transport.temperatures == [100.0, 200.0, 300.0] "Expected temperatures=[100, 200, 300]"
    println("  ✓ PASS")
    println()

    # Test 7: run_integrate with custom bins
    println("-" ^ 60)
    println("Test 7: run_integrate with bins=200")
    println("-" ^ 60)
    @info "Calling run_integrate with bins=200"
    transport2 = BoltzTraP.run_integrate(result1; temperatures=[300.0], bins=200, verbose=true)
    println("  → npts_dos in metadata: $(transport2.metadata["npts_dos"])")
    @assert transport2.metadata["npts_dos"] == 200 "Expected npts_dos=200"
    println("  ✓ PASS")
    println()

    println("=" ^ 60)
    println("All tests passed!")
    println("=" ^ 60)
end

main()
