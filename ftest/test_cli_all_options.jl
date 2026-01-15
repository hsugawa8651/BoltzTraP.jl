# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Functional test: All CLI options
#
# Usage:
#   julia --project=. ftest/test_cli_all_options.jl [vasp_directory]
#
# This script tests all CLI options for both `interpolate` and `integrate` commands.
# Run with JULIA_DEBUG=BoltzTraP to see debug output.

using Logging

# Enable debug logging for BoltzTraP module
ENV["JULIA_DEBUG"] = "BoltzTraP"
global_logger(ConsoleLogger(stderr, Logging.Debug))

using BoltzTraP

# Default test data path
const DEFAULT_VASP_DIR = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data", "Si.vasp")

function test_interpolate_options(vasp_dir::String)
    println()
    println("=" ^ 70)
    println("Testing: interpolate command options")
    println("=" ^ 70)

    # ------------------------------------------------------------------
    # Test: --format auto (default)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: --format auto (default)")
    println("-" ^ 70)
    @info "Testing format=auto"
    detected = BoltzTraP.detected_format(vasp_dir)
    println("  Detected format: $detected")
    data = BoltzTraP.load_dft(vasp_dir)
    result = BoltzTraP.run_interpolate(data; source=vasp_dir, kpoints=500, verbose=true)
    @assert result.metadata["source"] == vasp_dir
    println("  ✓ PASS: format auto-detection works")

    # ------------------------------------------------------------------
    # Test: --format vasp (explicit)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: --format vasp (explicit)")
    println("-" ^ 70)
    @info "Testing format=vasp"
    data_vasp = BoltzTraP.load_vasp(vasp_dir)
    result_vasp = BoltzTraP.run_interpolate(data_vasp; source=vasp_dir, kpoints=500, verbose=true)
    println("  ✓ PASS: explicit VASP loading works")

    # ------------------------------------------------------------------
    # Test: -k, --kpoints
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -k, --kpoints <n>")
    println("-" ^ 70)
    for kpts in [500, 1000, 2000]
        @info "Testing kpoints=$kpts"
        result = BoltzTraP.run_interpolate(data; source=vasp_dir, kpoints=kpts, verbose=false)
        actual = result.metadata["nkpt_target"]
        println("  kpoints=$kpts → nkpt_target=$actual")
        @assert actual == kpts "Expected nkpt_target=$kpts, got $actual"
    end
    println("  ✓ PASS: kpoints option works")

    # ------------------------------------------------------------------
    # Test: -m, --multiplier
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -m, --multiplier <n>")
    println("-" ^ 70)
    nkpt_orig = size(data.kpoints, 2)
    for mult in [2, 3]
        @info "Testing multiplier=$mult"
        result = BoltzTraP.run_interpolate(data; source=vasp_dir, multiplier=mult, verbose=false)
        expected = nkpt_orig * mult
        actual = result.metadata["nkpt_target"]
        println("  multiplier=$mult, nkpt_original=$nkpt_orig → nkpt_target=$actual (expected $expected)")
        @assert actual == expected "Expected nkpt_target=$expected, got $actual"
    end
    println("  ✓ PASS: multiplier option works")

    # ------------------------------------------------------------------
    # Test: --emin, --emax
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: --emin <e>, --emax <e>")
    println("-" ^ 70)
    test_cases = [
        (emin=-Inf, emax=Inf, desc="no filter"),
        (emin=-0.5, emax=0.5, desc="narrow range"),
        (emin=-1.0, emax=1.0, desc="wide range"),
    ]
    for tc in test_cases
        @info "Testing emin=$(tc.emin), emax=$(tc.emax)"
        result = BoltzTraP.run_interpolate(
            data; source=vasp_dir, kpoints=500, emin=tc.emin, emax=tc.emax, verbose=false
        )
        println("  $(tc.desc): emin=$(result.metadata["emin"]), emax=$(result.metadata["emax"])")
        @assert result.metadata["emin"] == tc.emin
        @assert result.metadata["emax"] == tc.emax
    end
    println("  ✓ PASS: emin/emax options work")

    # ------------------------------------------------------------------
    # Test: --absolute
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: --absolute")
    println("-" ^ 70)
    for absolute in [false, true]
        @info "Testing absolute=$absolute"
        result = BoltzTraP.run_interpolate(
            data; source=vasp_dir, kpoints=500, emin=-0.5, emax=0.5, absolute=absolute, verbose=false
        )
        println("  absolute=$absolute → metadata[\"absolute\"]=$(result.metadata["absolute"])")
        @assert result.metadata["absolute"] == absolute
    end
    println("  ✓ PASS: absolute option works")

    # ------------------------------------------------------------------
    # Test: -o, --output (file creation)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -o, --output <file>")
    println("-" ^ 70)
    tmpdir = mktempdir()
    outfile = joinpath(tmpdir, "test_output.jld2")
    @info "Testing output=$outfile"
    result = BoltzTraP.run_interpolate(data; source=vasp_dir, kpoints=500, output=outfile, verbose=false)
    @assert isfile(outfile) "Output file not created"
    println("  Output file created: $outfile")
    # Verify can be loaded back
    loaded = BoltzTraP.load_interpolation(outfile)
    @assert length(loaded.equivalences) == length(result.equivalences)
    println("  Output file can be reloaded")
    rm(tmpdir; recursive=true)
    println("  ✓ PASS: output option works")

    # ------------------------------------------------------------------
    # Test: -v, --verbose (no assertion, just check no error)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -v, --verbose")
    println("-" ^ 70)
    @info "Testing verbose=true"
    result = BoltzTraP.run_interpolate(data; source=vasp_dir, kpoints=500, verbose=true)
    println("  ✓ PASS: verbose option works")

    return result  # Return for integrate tests
end

function test_integrate_options(interp_result)
    println()
    println("=" ^ 70)
    println("Testing: integrate command options")
    println("=" ^ 70)

    # ------------------------------------------------------------------
    # Test: -t, --temperature (single value)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -t, --temperature (single value)")
    println("-" ^ 70)
    @info "Testing temperature=\"300\""
    temps = BoltzTraP._parse_temperatures("300")
    @assert temps == [300.0]
    result = BoltzTraP.run_integrate(interp_result; temperatures=temps, verbose=false)
    @assert result.temperatures == [300.0]
    println("  temperature=\"300\" → $temps")
    println("  ✓ PASS")

    # ------------------------------------------------------------------
    # Test: -t, --temperature (range format)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -t, --temperature (range: start:stop:step)")
    println("-" ^ 70)
    @info "Testing temperature=\"100:500:100\""
    temps = BoltzTraP._parse_temperatures("100:500:100")
    expected = [100.0, 200.0, 300.0, 400.0, 500.0]
    @assert temps == expected "Expected $expected, got $temps"
    result = BoltzTraP.run_integrate(interp_result; temperatures=temps, verbose=false)
    @assert result.temperatures == expected
    println("  temperature=\"100:500:100\" → $temps")
    println("  ✓ PASS")

    # ------------------------------------------------------------------
    # Test: -t, --temperature (range with default step)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -t, --temperature (range: start:stop, default step=100)")
    println("-" ^ 70)
    @info "Testing temperature=\"200:400\""
    temps = BoltzTraP._parse_temperatures("200:400")
    expected = [200.0, 300.0, 400.0]
    @assert temps == expected "Expected $expected, got $temps"
    println("  temperature=\"200:400\" → $temps")
    println("  ✓ PASS")

    # ------------------------------------------------------------------
    # Test: -t, --temperature (list format)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -t, --temperature (list: T1,T2,T3)")
    println("-" ^ 70)
    @info "Testing temperature=\"100,250,400\""
    temps = BoltzTraP._parse_temperatures("100,250,400")
    expected = [100.0, 250.0, 400.0]
    @assert temps == expected "Expected $expected, got $temps"
    result = BoltzTraP.run_integrate(interp_result; temperatures=temps, verbose=false)
    @assert result.temperatures == expected
    println("  temperature=\"100,250,400\" → $temps")
    println("  ✓ PASS")

    # ------------------------------------------------------------------
    # Test: -b, --bins
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -b, --bins <n>")
    println("-" ^ 70)
    for bins in [100, 200, 500]
        @info "Testing bins=$bins"
        result = BoltzTraP.run_integrate(interp_result; temperatures=[300.0], bins=bins, verbose=false)
        actual = result.metadata["npts_dos"]
        println("  bins=$bins → npts_dos=$actual")
        @assert actual == bins "Expected npts_dos=$bins, got $actual"
    end
    println("  ✓ PASS: bins option works")

    # ------------------------------------------------------------------
    # Test: -o, --output (file creation)
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -o, --output <file>")
    println("-" ^ 70)
    tmpdir = mktempdir()
    outfile = joinpath(tmpdir, "test_transport.jld2")
    @info "Testing output=$outfile"
    result = BoltzTraP.run_integrate(interp_result; temperatures=[300.0], output=outfile, verbose=false)
    @assert isfile(outfile) "Output file not created"
    println("  Output file created: $outfile")
    # Verify can be loaded back
    loaded = BoltzTraP.load_integrate(outfile)
    @assert length(loaded.temperatures) == length(result.temperatures)
    println("  Output file can be reloaded")
    rm(tmpdir; recursive=true)
    println("  ✓ PASS: output option works")

    # ------------------------------------------------------------------
    # Test: -v, --verbose
    # ------------------------------------------------------------------
    println()
    println("-" ^ 70)
    println("Test: -v, --verbose")
    println("-" ^ 70)
    @info "Testing verbose=true"
    result = BoltzTraP.run_integrate(interp_result; temperatures=[300.0], verbose=true)
    println("  ✓ PASS: verbose option works")
end

function main()
    vasp_dir = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_VASP_DIR

    if !isdir(vasp_dir)
        @error "VASP directory not found" vasp_dir
        println("\nUsage: julia --project=. ftest/test_cli_all_options.jl [vasp_directory]")
        println("\nExpected directory structure:")
        println("  vasp_directory/")
        println("    vasprun.xml")
        exit(1)
    end

    println("=" ^ 70)
    println("Functional Test: All CLI Options")
    println("=" ^ 70)
    println("VASP directory: $vasp_dir")

    # Test interpolate options
    interp_result = test_interpolate_options(vasp_dir)

    # Test integrate options
    test_integrate_options(interp_result)

    println()
    println("=" ^ 70)
    println("ALL TESTS PASSED!")
    println("=" ^ 70)
end

main()
