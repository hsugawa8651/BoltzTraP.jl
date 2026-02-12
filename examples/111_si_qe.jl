# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si Quantum ESPRESSO workflow
#
# Interpolation and transport calculation for silicon
# using Quantum ESPRESSO output (data-file-schema.xml).
#
# Usage:
#   julia --project=. examples/111_si_qe.jl

using BoltzTraP

# Path to QE data (directory containing silicon.xml)
datadir = joinpath(@__DIR__, "data", "Si.ESPRESSO", "out")

# Step 1: Load QE data explicitly
println("Loading QE data...")
data = load_qe(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="QE", kpoints=5000, verbose=true)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
