# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si VASP basic workflow
#
# Basic workflow: interpolation and transport coefficient calculation
# for silicon using VASP output (vasprun.xml + POSCAR).
#
# Usage:
#   julia --project=. examples/110_si_vasp.jl
#
# CLI equivalent:
#   boltztrap interpolate -k 5000 -o si_interp.jld2 benchmarks/data/Si.vasp
#   boltztrap integrate si_interp.jld2 -t 200:100:500 -o si_transport.jld2

using BoltzTraP

# Path to Si VASP data (directory containing vasprun.xml and POSCAR)
datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")

# Step 1: Interpolate band structure
println("Step 1: Interpolating band structure...")
interp = run_interpolate(datadir; kpoints=5000, verbose=true)

println("  Equivalence classes: $(length(interp.equivalences))")
println("  Bands: $(size(interp.coeffs, 1))")

# Step 2: Compute transport coefficients
println("\nStep 2: Computing transport coefficients...")
temperatures = [200.0, 300.0, 400.0, 500.0]
transport = run_integrate(interp; temperatures=temperatures, verbose=true)

# Step 3: Print results
println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
println("  Tensor shape (σ): $(size(transport.sigma))")

# Step 4: Save results (optional)
output_file = joinpath(@__DIR__, "si_transport.jld2")
save_integrate(output_file, transport)
println("\nSaved to $output_file")
