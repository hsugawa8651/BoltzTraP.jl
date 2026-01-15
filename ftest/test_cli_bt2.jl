# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Functional test: CLI commands with .bt2 files
#
# Usage:
#   julia --project=. ftest/test_cli_bt2.jl [bt2_file]
#
# This script tests CLI commands (describe, plotbands) with Python-generated .bt2 files.

using BoltzTraP

# Default test data path
const DEFAULT_BT2_FILE = joinpath(
    @__DIR__, "..", "..", "BoltzTraP2-public", "data", "Si.vasp", "interpolation.bt2"
)

function main()
    bt2_file = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_BT2_FILE

    if !isfile(bt2_file)
        @error "bt2 file not found" bt2_file
        println("\nUsage: julia --project=. ftest/test_cli_bt2.jl [bt2_file]")
        println("\nExpected: Python BoltzTraP2 generated .bt2 file")
        exit(1)
    end

    println("=" ^ 60)
    println("Functional Test: CLI bt2 Commands")
    println("=" ^ 60)
    println("bt2 file: $bt2_file")
    println()

    # Test 1: describe command with .bt2 file
    println("-" ^ 60)
    println("Test 1: describe command with .bt2 file")
    println("-" ^ 60)
    try
        BoltzTraP.describe(bt2_file)
        println("  ✓ PASS: describe command succeeded")
    catch e
        println("  ✗ FAIL: $(e)")
        exit(1)
    end
    println()

    # Test 2: load_interpolation with .bt2 file
    println("-" ^ 60)
    println("Test 2: load_interpolation with .bt2 file")
    println("-" ^ 60)
    result = load_interpolation(bt2_file)
    @assert !isempty(result.coeffs) "coeffs should not be empty"
    @assert !isempty(result.equivalences) "equivalences should not be empty"
    println("  → coeffs shape: $(size(result.coeffs))")
    println("  → equivalences count: $(length(result.equivalences))")
    println("  → has spacegroup: $(haskey(result.metadata, "spacegroup_number"))")
    println("  ✓ PASS")
    println()

    # Test 3: plotbands command with manual k-path (no spacegroup required)
    println("-" ^ 60)
    println("Test 3: plotbands with manual --kpath")
    println("-" ^ 60)
    output_file = tempname() * ".png"
    kpath_str = "G:0,0,0;X:0.5,0,0.5;L:0.5,0.5,0.5;W:0.5,0.25,0.75|G-X-W-L-G"
    try
        BoltzTraP.plotbands(bt2_file; npoints=30, kpath=kpath_str, output=output_file)
        if isfile(output_file)
            fsize = filesize(output_file)
            println("  → output file: $output_file ($fsize bytes)")
            rm(output_file)
            println("  ✓ PASS: plotbands with manual kpath succeeded")
        else
            println("  ✗ FAIL: output file not created")
            exit(1)
        end
    catch e
        println("  ✗ FAIL: $(e)")
        exit(1)
    end
    println()

    # Test 4: plotbands error without kpath (when no spacegroup)
    println("-" ^ 60)
    println("Test 4: plotbands error without kpath (when no spacegroup)")
    println("-" ^ 60)
    if !haskey(result.metadata, "spacegroup_number")
        try
            BoltzTraP.plotbands(bt2_file; npoints=10, output=tempname() * ".png")
            println("  ✗ FAIL: should have raised error")
            exit(1)
        catch e
            if occursin("Space group not found", string(e))
                println("  → Got expected error: Space group not found")
                println("  ✓ PASS: error message is correct")
            else
                println("  ✗ FAIL: unexpected error: $(e)")
                exit(1)
            end
        end
    else
        println("  → SKIP: file has spacegroup metadata")
    end
    println()

    # Test 5: plot_bands API with manual kpath NamedTuple
    println("-" ^ 60)
    println("Test 5: plot_bands API with manual kpath NamedTuple")
    println("-" ^ 60)
    output_file = tempname() * ".png"
    kpath_tuple = (
        points = Dict(
            "G" => [0.0, 0.0, 0.0],
            "X" => [0.5, 0.0, 0.5],
            "L" => [0.5, 0.5, 0.5],
        ),
        paths = [["G", "X", "L", "G"]],
    )
    try
        p = plot_bands(result; npoints=30, kpath=kpath_tuple, output=output_file)
        if isfile(output_file)
            fsize = filesize(output_file)
            println("  → output file: $output_file ($fsize bytes)")
            rm(output_file)
            println("  ✓ PASS")
        else
            println("  ✗ FAIL: output file not created")
            exit(1)
        end
    catch e
        println("  ✗ FAIL: $(e)")
        exit(1)
    end
    println()

    # Test 6: _parse_kpath_string validation
    println("-" ^ 60)
    println("Test 6: _parse_kpath_string validation")
    println("-" ^ 60)
    # Valid format
    kpath = BoltzTraP._parse_kpath_string("A:0,0,0;B:1,0,0|A-B-A")
    @assert haskey(kpath.points, "A")
    @assert haskey(kpath.points, "B")
    @assert kpath.paths[1] == ["A", "B", "A"]
    println("  → valid format parsed correctly")

    # Invalid format (missing |)
    try
        BoltzTraP._parse_kpath_string("A:0,0,0;B:1,0,0")
        println("  ✗ FAIL: should have raised error for missing |")
        exit(1)
    catch e
        if occursin("Invalid kpath format", string(e))
            println("  → missing | detected correctly")
        else
            println("  ✗ FAIL: unexpected error: $(e)")
            exit(1)
        end
    end

    # Invalid format (unknown label)
    try
        BoltzTraP._parse_kpath_string("A:0,0,0|A-C-A")
        println("  ✗ FAIL: should have raised error for unknown label")
        exit(1)
    catch e
        if occursin("Unknown point label", string(e))
            println("  → unknown label detected correctly")
        else
            println("  ✗ FAIL: unexpected error: $(e)")
            exit(1)
        end
    end
    println("  ✓ PASS")
    println()

    println("=" ^ 60)
    println("All tests passed!")
    println("=" ^ 60)
end

main()
