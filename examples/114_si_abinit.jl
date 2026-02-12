# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si ABINIT workflow
#
# Interpolation and transport calculation for silicon
# using ABINIT output (GSR.nc NetCDF file).
#
# Note: Requires NCDatasets.jl package for reading NetCDF files.
#   using Pkg; Pkg.add("NCDatasets")
#
# Usage:
#   julia --project=. examples/114_si_abinit.jl

using NCDatasets  # Required before loading BoltzTraP for ABINIT support
using BoltzTraP

# Path to ABINIT data (directory containing *_GSR.nc)
datadir = joinpath(@__DIR__, "data", "Si.abinit")

# Step 1: Load ABINIT data explicitly
println("Loading ABINIT data...")
data = load_abinit(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="ABINIT", kpoints=5000, verbose=true)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
