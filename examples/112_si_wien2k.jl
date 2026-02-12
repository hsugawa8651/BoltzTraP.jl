# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si Wien2k workflow
#
# Interpolation and transport calculation for silicon
# using Wien2k output (case.struct, case.energy, case.scf).
#
# Data: Wien2k Si data is NOT bundled due to size (15 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Si
#       Place the downloaded "Si" directory as examples/data/Si.wien2k/
#
# Usage:
#   julia --project=. examples/112_si_wien2k.jl

using BoltzTraP

# Path to Wien2k data
datadir = joinpath(@__DIR__, "data", "Si.wien2k")

if !isdir(datadir)
    error("""Wien2k Si data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Si
    Place as: examples/data/Si.wien2k/""")
end

# Step 1: Load Wien2k data explicitly
println("Loading Wien2k data...")
data = load_wien2k(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="Wien2k", kpoints=5000, verbose=true)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
