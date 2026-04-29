# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Validation script for calc_cv against Python BoltzTraP2.

Usage:
    julia --project validation/compare_cv.jl [--verbose]

This script compares the Julia calc_cv implementation against
Python BoltzTraP2 reference data.
"""

using NPZ
using BoltzTraP
using Printf

function main()
    verbose = "--verbose" in ARGS || "-v" in ARGS

    println("=" ^ 60)
    println("calc_cv Validation: Julia vs Python BoltzTraP2")
    println("=" ^ 60)

    # Reference data path
    ref_path = joinpath(@__DIR__, "..", "reftest", "data", "cv_reference.npz")

    if !isfile(ref_path)
        println("\nError: Reference data not found at: $ref_path")
        println("\nTo generate reference data, run:")
        println("  cd BoltzTraP2-public")
        println("  python ../BoltzTraP.jl/reftest/generate_cv_reference.py")
        exit(1)
    end

    # Load reference data
    ref = npzread(ref_path)

    println("\n[Input Data]")
    println("  Material: Si (from BoltzTraP2 bundled data)")
    @printf(
        "  epsilon range: %.4f to %.4f Ha\n",
        minimum(ref["epsilon"]),
        maximum(ref["epsilon"])
    )
    println("  DOS points: $(length(ref["dos"]))")
    println("  mu points: $(length(ref["mu_range"]))")
    println("  T range: $(ref["T_range"]) K")
    println("  dosweight: $(ref["dosweight"])")

    # Julia calculation
    cv_julia = BoltzTraP.calc_cv(
        ref["epsilon"],
        ref["dos"],
        ref["mu_range"],
        ref["T_range"];
        dosweight = ref["dosweight"],
    )
    cv_python = ref["cv"]

    println("\n[Results Comparison]")
    println("  Shape: Julia $(size(cv_julia)) vs Python $(size(cv_python))")

    # Find Fermi level index
    fermi_idx = argmin(abs.(ref["mu_range"] .- ref["fermi"]))

    println("\n  cv values at mu = Fermi level (index $fermi_idx):")
    println("  " * "-" ^ 54)
    @printf(
        "  %6s  %15s  %15s  %12s\n",
        "T [K]",
        "Julia [J/K]",
        "Python [J/K]",
        "Rel. Diff"
    )
    println("  " * "-" ^ 54)
    for (i, T) in enumerate(ref["T_range"])
        jl = cv_julia[i, fermi_idx]
        py = cv_python[i, fermi_idx]
        rel_diff = abs(jl - py) / max(abs(py), 1e-30)
        @printf("  %6.0f  %15.6e  %15.6e  %12.2e\n", T, jl, py, rel_diff)
    end

    if verbose
        println("\n  Full comparison (all mu values at T = 300 K):")
        println("  " * "-" ^ 54)
        @printf(
            "  %6s  %15s  %15s  %12s\n",
            "mu idx",
            "Julia [J/K]",
            "Python [J/K]",
            "Rel. Diff"
        )
        println("  " * "-" ^ 54)
        T_idx = findfirst(==(300.0), ref["T_range"])
        if !isnothing(T_idx)
            for i = 1:length(ref["mu_range"])
                jl = cv_julia[T_idx, i]
                py = cv_python[T_idx, i]
                rel_diff = abs(jl - py) / max(abs(py), 1e-30)
                @printf("  %6d  %15.6e  %15.6e  %12.2e\n", i, jl, py, rel_diff)
            end
        end
    end

    # Statistics
    println("\n[Statistics]")
    max_abs = maximum(abs.(cv_julia - cv_python))
    max_rel = maximum(abs.(cv_julia - cv_python) ./ max.(abs.(cv_python), 1e-30))
    mean_rel =
        sum(abs.(cv_julia - cv_python) ./ max.(abs.(cv_python), 1e-30)) / length(cv_python)

    @printf("  Max absolute difference: %.3e J/K\n", max_abs)
    @printf("  Max relative difference: %.3e\n", max_rel)
    @printf("  Mean relative difference: %.3e\n", mean_rel)

    # Pass/fail criteria
    tol = 1e-5
    passed = max_rel < tol

    println("\n[Validation Result]")
    println("  Tolerance: rtol = $tol")
    if passed
        println("  Status: PASS")
    else
        println("  Status: FAIL")
    end

    println("\n" * "=" ^ 60)

    return passed ? 0 : 1
end

exit(main())
