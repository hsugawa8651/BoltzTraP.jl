# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si VASP band structure plot
#
# Band structure interpolation and visualization for silicon
# using VASP output (vasprun.xml + POSCAR).
#
# DFT conditions:
#   Code: VASP (PAW_PBE Si 05Jan2001)
#   XC functional: PBE-GGA
#   k-grid: 17×17×17 (165 IBZ k-points)
#   Lattice constant: a = 5.467 Å (10.33 bohr)
#   Source: BoltzTraP2 distribution (Si.vasp)
#
# Usage:
#   julia --project=. examples/710_si_vasp_plotbands.jl

using BoltzTraP

# Path to Si VASP data
datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 1: Interpolate band structure
println("Interpolating band structure...")
interp = run_interpolate(datadir; kpoints=5000, verbose=true)

# Step 2: Plot band structure along high-symmetry path
println("\nPlotting band structure...")
output_png = joinpath(@__DIR__, stem * "_bands.png")
plot_bands(interp; emin=-5.0, emax=5.0, output=output_png)
println("Saved plot to $output_png")
