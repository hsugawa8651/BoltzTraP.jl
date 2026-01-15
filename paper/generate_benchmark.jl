# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Generate Figure 2 for JOSS paper:
# Scaling benchmark comparing Python BoltzTraP2 vs Julia BoltzTraP.jl
#
# Usage:
#   julia --project paper/generate_figure2.jl
#
# Requires:
#   benchmarks/scaling_results_python.json
#   benchmarks/scaling_results_julia.json
#
# Output:
#   paper/figure2.png

using JSON
using Plots
using LaTeXStrings

gr()

function load_benchmark_results()
    julia_file = joinpath(@__DIR__, "..", "benchmarks", "scaling_results_julia.json")
    python_file = joinpath(@__DIR__, "..", "benchmarks", "scaling_results_python.json")

    if !isfile(julia_file)
        error("Julia benchmark results not found: $julia_file\n" *
              "Run: julia --project benchmarks/scaling_benchmark.jl")
    end
    if !isfile(python_file)
        error("Python benchmark results not found: $python_file\n" *
              "Run: python benchmarks/scaling_benchmark.py")
    end

    julia_data = JSON.parsefile(julia_file)
    python_data = JSON.parsefile(python_file)

    return julia_data, python_data
end

function create_benchmark_figure(julia_data, python_data, output_path)
    nk_values = julia_data["nk_values"]

    # Extract median times
    python_times = [python_data["python"]["$nk"]["median"] for nk in nk_values]
    julia_times = [julia_data["optimized"]["$nk"]["median"] for nk in nk_values]

    # Extract Julia allocations (MiB)
    julia_allocs = [julia_data["optimized"]["$nk"]["alloc_mb"] for nk in nk_values]

    # Convert to seconds for cleaner display
    python_times_s = python_times ./ 1000
    julia_times_s = julia_times ./ 1000

    # Create figure with dual y-axis
    fig = plot(
        size=(700, 400),
        dpi=300,
        left_margin=5Plots.mm,
        right_margin=12Plots.mm,
        bottom_margin=5Plots.mm,
        fontfamily="Computer Modern"
    )

    # Left axis: Time
    scatter!(
        fig,
        nk_values, python_times_s,
        label="Python BoltzTraP2 (time)",
        markersize=6,
        markershape=:circle,
        markercolor=:blue,
        markerstrokewidth=0
    )
    plot!(
        fig,
        nk_values, python_times_s,
        label="",
        linewidth=1.5,
        linestyle=:dash,
        color=:blue,
        alpha=0.5
    )

    scatter!(
        fig,
        nk_values, julia_times_s,
        label="Julia BoltzTraP.jl (time)",
        markersize=6,
        markershape=:star5,
        markercolor=:red,
        markerstrokewidth=0
    )
    plot!(
        fig,
        nk_values, julia_times_s,
        label="",
        linewidth=1.5,
        linestyle=:solid,
        color=:red,
        alpha=0.5
    )

    plot!(
        fig,
        xlabel=L"Number of $k$-points",
        ylabel="Time (s)",
        title="getBands Performance Scaling",
        legend=:topleft,
        grid=true,
        gridalpha=0.3,
        xlims=(0, maximum(nk_values) * 1.1),
        ylims=(0, maximum(python_times_s) * 1.1)
    )

    # Right axis: Memory allocation (Julia only)
    fig2 = twinx(fig)
    scatter!(
        fig2,
        nk_values, julia_allocs,
        label="Julia allocation (MiB)",
        markersize=5,
        markershape=:diamond,
        markercolor=:green,
        markerstrokewidth=0
    )
    plot!(
        fig2,
        nk_values, julia_allocs,
        label="",
        linewidth=1.5,
        linestyle=:dot,
        color=:green,
        alpha=0.5
    )
    plot!(
        fig2,
        ylabel="Allocation (MiB)",
        legend=:topright,
        ylims=(0, maximum(julia_allocs) * 1.2)
    )

    savefig(fig, output_path)
    println("Figure saved: $output_path")

    # Print speedup summary
    println("\nSpeedup Summary (vs Python):")
    println("  nk     | Python (ms) | Julia (ms) | Speedup | Alloc (MiB)")
    println("  -------|-------------|------------|---------|------------")
    for (i, nk) in enumerate(nk_values)
        speedup = python_times[i] / julia_times[i]
        println("  $(lpad(nk, 5)) | $(lpad(round(python_times[i], digits=0), 11)) | $(lpad(round(julia_times[i], digits=0), 10)) | $(lpad(round(speedup, digits=1), 7))x | $(lpad(round(julia_allocs[i], digits=1), 10))")
    end

    return fig
end

function main()
    println("=" ^ 60)
    println("Generating JOSS Figure 2: Benchmark Scaling")
    println("=" ^ 60)
    println()

    julia_data, python_data = load_benchmark_results()

    output_path = joinpath(@__DIR__, "figure2.png")
    create_benchmark_figure(julia_data, python_data, output_path)

    println("\n" * "=" ^ 60)
    println("Done!")
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
