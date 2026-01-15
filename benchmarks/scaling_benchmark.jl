# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# End-to-end benchmark: measure interpolate + integrate performance.
# Compares Julia BoltzTraP.jl with Python BoltzTraP2.
#
# Usage:
#   julia --project benchmarks/scaling_benchmark.jl
#
# Output:
#   benchmarks/scaling_results_julia.json

using BoltzTraP
using JSON
using Dates

"""
    run_single_benchmark(data, nkpt; temperature=300.0, bins=500)

Run a single interpolate + integrate cycle and return elapsed times.

# Returns
Tuple of (total_ms, interpolate_ms, integrate_ms)
"""
function run_single_benchmark(data, nkpt; temperature=300.0, bins=500)
    # Interpolate
    t_interp = @elapsed begin
        interp = run_interpolate(data; kpoints=nkpt, verbose=false)
    end
    interpolate_ms = t_interp * 1000

    # Integrate
    t_integ = @elapsed begin
        result = run_integrate(interp; temperatures=[temperature], bins=bins, verbose=false)
    end
    integrate_ms = t_integ * 1000

    total_ms = interpolate_ms + integrate_ms
    return total_ms, interpolate_ms, integrate_ms
end

"""
    measure_memory(data, nkpt; temperature=300.0, bins=500)

Measure memory allocation for a single run using @allocated.

# Returns
Peak memory in MiB (approximated by total allocations)
"""
function measure_memory(data, nkpt; temperature=300.0, bins=500)
    # Force GC before measurement
    GC.gc()

    bytes = @allocated begin
        interp = run_interpolate(data; kpoints=nkpt, verbose=false)
        result = run_integrate(interp; temperatures=[temperature], bins=bins, verbose=false)
    end

    return bytes / (1024 * 1024)  # Convert to MiB
end

"""
    run_scaling_benchmark(; datadir, nkpt_values, n_runs, temperature, bins)

Run the full scaling benchmark.
"""
function run_scaling_benchmark(;
    datadir::String = "benchmarks/data/Si.vasp",
    nkpt_values::Vector{Int} = [1000, 2000, 4000],
    n_runs::Int = 5,
    temperature::Float64 = 300.0,
    bins::Int = 500,
)
    println("=" ^ 60)
    println("Julia BoltzTraP.jl End-to-End Benchmark")
    println("=" ^ 60)
    println()

    # Load data (I/O excluded from measurement)
    println("Loading VASP data: $datadir")
    data = load_vasp(datadir)
    println("  Bands: $(size(data.ebands, 1))")
    println("  K-points (DFT): $(size(data.kpoints, 2))")
    println("  Fermi: $(round(data.fermi, digits=6)) Ha")
    println()

    # Warm up (important for JIT compilation)
    println("Warming up (JIT compilation)...")
    run_single_benchmark(data, 500; temperature, bins)
    println()

    results = Dict{String,Any}(
        "language" => "Julia",
        "version" => string(VERSION),
        "package" => "BoltzTraP.jl",
        "timestamp" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "parameters" => Dict{String,Any}(
            "temperature_K" => temperature,
            "dos_bins" => bins,
            "n_runs" => n_runs,
        ),
        "results" => []
    )

    println("-" ^ 60)
    println("Benchmarking (T=$(temperature)K, DOS bins=$bins)")
    println("-" ^ 60)
    println()

    for nkpt in nkpt_values
        println("nkpt=$nkpt:")

        # Time measurements
        times_total = Float64[]
        times_interp = Float64[]
        times_integ = Float64[]

        for run in 1:n_runs
            total, interp, integ = run_single_benchmark(data, nkpt; temperature, bins)
            push!(times_total, total)
            push!(times_interp, interp)
            push!(times_integ, integ)
            println("  Run $run: $(round(total, digits=1)) ms (interp: $(round(interp, digits=1)), integ: $(round(integ, digits=1)))")
        end

        # Memory measurement (separate run)
        peak_mem = measure_memory(data, nkpt; temperature, bins)

        median_total = median(times_total)
        median_interp = median(times_interp)
        median_integ = median(times_integ)
        std_total = std(times_total)

        println("  Median: $(round(median_total, digits=1)) ms (interp: $(round(median_interp, digits=1)), integ: $(round(median_integ, digits=1)))")
        println("  Allocations: $(round(peak_mem, digits=1)) MiB")
        println()

        push!(results["results"], Dict{String,Any}(
            "kpoints" => nkpt,
            "total_ms" => median_total,
            "interpolate_ms" => median_interp,
            "integrate_ms" => median_integ,
            "std_ms" => std_total,
            "peak_memory_MiB" => peak_mem,
        ))
    end

    # Save results
    script_dir = dirname(@__FILE__)
    output_file = joinpath(script_dir, "scaling_results_julia.json")
    open(output_file, "w") do f
        JSON.print(f, results, 2)
    end
    println("Results saved: $output_file")
    println()

    # Summary table
    println("=" ^ 60)
    println("Summary")
    println("=" ^ 60)
    println()
    println(lpad("kpoints", 8), " | ", lpad("Total (ms)", 12), " | ", lpad("Interp (ms)", 12), " | ", lpad("Integ (ms)", 12), " | ", lpad("Memory (MiB)", 12))
    println("-" ^ 70)
    for r in results["results"]
        println(
            lpad(r["kpoints"], 8), " | ",
            lpad(round(r["total_ms"], digits=1), 12), " | ",
            lpad(round(r["interpolate_ms"], digits=1), 12), " | ",
            lpad(round(r["integrate_ms"], digits=1), 12), " | ",
            lpad(round(r["peak_memory_MiB"], digits=1), 12)
        )
    end
    println()

    return results
end

# Use Statistics for median and std
using Statistics

# Run benchmark
if abspath(PROGRAM_FILE) == @__FILE__
    run_scaling_benchmark()
end
