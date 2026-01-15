# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# End-to-end test: DFTK → BoltzTraP.jl transport calculation
#
# Usage:
#   julia --project ftest/test_dftk_e2e.jl

println("=" ^ 60)
println("DFTK → BoltzTraP.jl End-to-End Test")
println("=" ^ 60)
println()

using DFTK
using BoltzTraP

# Note: DFTK will print timing summary at the end

println("Step 1: Setting up Si crystal structure...")
# Silicon in diamond structure
a = 10.26  # Lattice constant in Bohr
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]

Si = ElementPsp(:Si; psp=load_psp("hgh/lda/si-q4"))
atoms = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

println("  Lattice constant: $a Bohr")
println("  Atoms: Si (2 atoms, diamond structure)")
println()

println("Step 2: Running DFTK SCF calculation...")
println("  (This may take a minute...)")
# Use moderate k-grid for testing (production would use 10x10x10 or more)
model = model_LDA(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut=15, kgrid=[6, 6, 6])

println("  Ecut: 15 Ha")
println("  k-grid: 6×6×6")
println("  Number of k-points: $(length(basis.kpoints))")

scfres = self_consistent_field(basis; tol=1e-6)
println("  SCF converged!")
println("  Fermi energy: $(scfres.εF) Ha")
println()

println("Step 3: Extracting data with load_dftk...")
data = load_dftk(scfres)
println("  Lattice shape: $(size(data.lattice))")
println("  K-points: $(size(data.kpoints, 2))")
println("  Bands: $(size(data.ebands, 1))")
println("  Fermi: $(data.fermi) Ha")
println("  n_elect: $(data.nelect)")
println()

println("Step 4: Running BoltzTraP interpolation...")
# Use lower kpoints for faster testing
interp = run_interpolate(data; source="DFTK", kpoints=1000, verbose=false)
println("  Equivalences: $(length(interp.equivalences))")
println("  Coefficients shape: $(size(interp.coeffs))")
println()

println("Step 5: Running BoltzTraP integration...")
transport = run_integrate(interp; temperatures=[300.0], verbose=false)
println("  Temperatures: $(transport.temperatures) K")
println("  μ range: $(minimum(transport.mu_values)) to $(maximum(transport.mu_values)) Ha")
println()

println("Step 6: Checking results...")
# Check conductivity tensor shape
sigma = transport.sigma
println("  σ tensor shape: $(size(sigma))")

# Get conductivity at Fermi level (μ closest to 0)
fermi_idx = argmin(abs.(transport.mu_values .- data.fermi))
sigma_at_fermi = sigma[:, :, 1, fermi_idx]
println("  σ at Fermi level (T=300K):")
println("    σ_xx = $(sigma_at_fermi[1,1]) S/(m·s)")
println("    σ_yy = $(sigma_at_fermi[2,2]) S/(m·s)")
println("    σ_zz = $(sigma_at_fermi[3,3]) S/(m·s)")

# Check Seebeck coefficient
seebeck = transport.seebeck
S_at_fermi = seebeck[:, :, 1, fermi_idx]
println("  S at Fermi level (T=300K):")
println("    S_xx = $(S_at_fermi[1,1]) V/K")

# For silicon (semiconductor), conductivity should be near zero at Fermi level
# Seebeck should be undefined (0/0) but numerically small
println()

println("=" ^ 60)
println("End-to-end test completed successfully!")
println("=" ^ 60)
