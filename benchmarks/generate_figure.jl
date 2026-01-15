# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Generate benchmark comparison figure from JSON results.
#
# Usage:
#   julia --project benchmarks/generate_figure.jl
#
# Input:
#   benchmarks/scaling_results_python.json
#   benchmarks/scaling_results_julia.json
#
# Output:
#   benchmarks/benchmark.png
#   benchmarks/benchmark_breakdown.png

using JSON
using PythonPlot

function generate_figure(;
    python_file::String = "benchmarks/scaling_results_python.json",
    julia_file::String = "benchmarks/scaling_results_julia.json",
    output_file::String = "benchmarks/benchmark.png",
)
    println("Loading benchmark results...")

    # Load results
    python_data = JSON.parsefile(python_file)
    julia_data = JSON.parsefile(julia_file)

    # Extract data
    kpoints = [r["kpoints"] for r in julia_data["results"]]
    python_times = [r["total_ms"] / 1000 for r in python_data["results"]]
    julia_times = [r["total_ms"] / 1000 for r in julia_data["results"]]
    julia_memory = [r["peak_memory_MiB"] for r in julia_data["results"]]

    # Breakdown data
    python_interp = [r["interpolate_ms"] / 1000 for r in python_data["results"]]
    python_integ = [r["integrate_ms"] / 1000 for r in python_data["results"]]
    julia_interp = [r["interpolate_ms"] / 1000 for r in julia_data["results"]]
    julia_integ = [r["integrate_ms"] / 1000 for r in julia_data["results"]]

    println("  Python: $(python_data["language"]) $(python_data["version"])")
    println("  Julia: $(julia_data["language"]) $(julia_data["version"])")
    println("  K-points: $kpoints")
    println()

    plt = PythonPlot.pyplot

    # ========================================
    # Figure 1: End-to-End Performance Scaling
    # ========================================
    fig1, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(kpoints, python_times, "o--", color="blue", markersize=10, linewidth=2,
             label="BoltzTraP2 (time)")
    ax1.plot(kpoints, julia_times, "*-", color="red", markersize=12, linewidth=2,
             label="BoltzTraP.jl (time)")

    ax1.set_xlabel("Number of k-points", fontsize=14)
    ax1.set_ylabel("Time (s)", fontsize=14)
    ax1.set_title("End-to-End Performance Scaling", fontsize=16)
    ax1.set_xlim(0, maximum(kpoints) * 1.1)
    ax1.set_ylim(0, maximum(python_times) * 1.15)
    ax1.grid(true, alpha=0.3)
    ax1.tick_params(direction="in", which="both", labelsize=12)

    # Secondary y-axis for memory
    ax1b = ax1.twinx()
    ax1b.plot(kpoints, julia_memory, "D:", color="green", markersize=8, linewidth=2,
              label="BoltzTraP.jl allocation (MiB)")
    ax1b.set_ylabel("Allocation (MiB)", fontsize=14)
    ax1b.set_ylim(0, maximum(julia_memory) * 1.15)
    ax1b.tick_params(direction="in", labelsize=12)

    # Combined legend - add dummy line to ax1 for allocation
    ax1.plot([], [], "D:", color="green", markersize=8, linewidth=2,
             label="BoltzTraP.jl allocation (MiB)")
    ax1.legend(loc="upper left", fontsize=12)

    fig1.tight_layout()
    fig1.savefig(output_file, dpi=150)
    println("Saved: $output_file")
    plt.close(fig1)

    # ========================================
    # Figure 2: Performance Breakdown
    # ========================================
    fig2, ax2 = plt.subplots(figsize=(10, 6))

    bar_width = 120
    x_py = kpoints .- bar_width/2
    x_jl = kpoints .+ bar_width/2

    # Python bars (stacked)
    ax2.bar(x_py, python_interp, bar_width, color="royalblue", alpha=0.8)
    ax2.bar(x_py, python_integ, bar_width, bottom=python_interp, color="blue")

    # Julia bars (stacked)
    ax2.bar(x_jl, julia_interp, bar_width, color="lightsalmon", alpha=0.8)
    ax2.bar(x_jl, julia_integ, bar_width, bottom=julia_interp, color="red")

    ax2.set_xlabel("Number of k-points", fontsize=16)
    ax2.set_ylabel("Time (s)", fontsize=16)
    ax2.set_title("Performance Breakdown", fontsize=18)
    ax2.set_xticks(kpoints)
    ax2.set_xlim(0, maximum(kpoints) + 500)
    ax2.set_ylim(0, maximum(python_times) * 1.25)
    ax2.grid(true, axis="y", alpha=0.3)
    ax2.tick_params(direction="in", which="both", labelsize=14)

    # Add speedup annotations above each bar pair
    for (i, k) in enumerate(kpoints)
        speedup = python_times[i] / julia_times[i]
        y_pos = python_times[i] + 0.25
        # Shift the last annotation to the left to avoid overlap with allocation marker
        x_offset = (i == length(kpoints)) ? -250 : 0
        ax2.text(k + x_offset, y_pos, "$(round(speedup, digits=1))x", ha="center", va="bottom",
                 fontsize=14, fontweight="bold", color="black")
    end

    # Secondary y-axis for memory
    ax2b = ax2.twinx()
    ax2b.plot(kpoints, julia_memory, "D:", color="green", markersize=8, linewidth=2,
              label="allocation (MiB)")
    ax2b.set_ylabel("Allocation (MiB)", fontsize=16)
    ax2b.set_ylim(0, maximum(julia_memory) * 1.15)
    ax2b.tick_params(direction="in", labelsize=14)

    # Custom 2x2 legend for colors + allocation
    # Create proxy artists for legend
    mpatches = plt.matplotlib.patches

    # Header row (empty + BoltzTraP2 + BoltzTraP.jl)
    # integrate row
    patch_py_integ = mpatches.Patch(color="blue", label="")
    patch_jl_integ = mpatches.Patch(color="red", label="")
    # interpolate row
    patch_py_interp = mpatches.Patch(color="royalblue", alpha=0.8, label="")
    patch_jl_interp = mpatches.Patch(color="lightsalmon", alpha=0.8, label="")

    # Create table-like legend using annotations
    ax2.annotate("", xy=(0.02, 0.98), xycoords="axes fraction")

    # Use text annotations for the legend table
    legend_x = 0.02
    legend_y = 0.96
    row_height = 0.07
    col_width = 0.16

    # Headers
    ax2.text(legend_x + col_width, legend_y, "BoltzTraP2", transform=ax2.transAxes,
             fontsize=12, fontweight="bold", ha="center", va="top", color="blue")
    ax2.text(legend_x + 2*col_width, legend_y, "BoltzTraP.jl", transform=ax2.transAxes,
             fontsize=12, fontweight="bold", ha="center", va="top", color="red")

    # Row labels and color boxes
    ax2.text(legend_x, legend_y - row_height, "integrate", transform=ax2.transAxes,
             fontsize=12, ha="left", va="top")
    ax2.text(legend_x, legend_y - 2*row_height, "interpolate", transform=ax2.transAxes,
             fontsize=12, ha="left", va="top")

    # Color boxes (using small rectangles)
    box_width = 0.06
    box_height = 0.03
    # integrate row
    ax2.add_patch(plt.matplotlib.patches.Rectangle(
        (legend_x + col_width - box_width/2, legend_y - row_height - box_height*0.7),
        box_width, box_height, transform=ax2.transAxes, color="blue", clip_on=false))
    ax2.add_patch(plt.matplotlib.patches.Rectangle(
        (legend_x + 2*col_width - box_width/2, legend_y - row_height - box_height*0.7),
        box_width, box_height, transform=ax2.transAxes, color="red", clip_on=false))
    # interpolate row
    ax2.add_patch(plt.matplotlib.patches.Rectangle(
        (legend_x + col_width - box_width/2, legend_y - 2*row_height - box_height*0.7),
        box_width, box_height, transform=ax2.transAxes, color="royalblue", alpha=0.8, clip_on=false))
    ax2.add_patch(plt.matplotlib.patches.Rectangle(
        (legend_x + 2*col_width - box_width/2, legend_y - 2*row_height - box_height*0.7),
        box_width, box_height, transform=ax2.transAxes, color="lightsalmon", alpha=0.8, clip_on=false))

    # Allocation legend (separate, below the table)
    ax2.plot([], [], "D:", color="green", markersize=8, linewidth=2, label="BoltzTraP.jl allocation (MiB)")
    ax2.legend(loc="upper left", fontsize=12, bbox_to_anchor=(0.0, 0.74))

    fig2.tight_layout()
    breakdown_file = replace(output_file, ".png" => "_breakdown.png")
    fig2.savefig(breakdown_file, dpi=150)
    println("Saved: $breakdown_file")
    plt.close(fig2)

    # Print summary
    println()
    println("=" ^ 55)
    println("Performance Summary")
    println("=" ^ 55)
    println()
    println(lpad("kpoints", 10), " | ", lpad("Python (s)", 12), " | ", lpad("Julia (s)", 12), " | ", lpad("Speedup", 10))
    println("-" ^ 55)
    for (i, k) in enumerate(kpoints)
        speedup = python_times[i] / julia_times[i]
        println(
            lpad(k, 10), " | ",
            lpad(round(python_times[i], digits=2), 12), " | ",
            lpad(round(julia_times[i], digits=2), 12), " | ",
            lpad("$(round(speedup, digits=2))x", 10)
        )
    end
    println()
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_figure()
end
