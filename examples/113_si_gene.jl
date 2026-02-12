# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si GENE (BoltzTraP1) workflow
#
# Interpolation and transport calculation for silicon
# using GENE/generic format (case.structure, case.energy).
#
# Usage:
#   julia --project=. examples/113_si_gene.jl

using BoltzTraP

# Path to GENE data
datadir = joinpath(@__DIR__, "data", "Si.GENE")

# Step 1: Load GENE data explicitly
println("Loading GENE data...")
data = load_gene(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="GENE", kpoints=5000, verbose=true)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
