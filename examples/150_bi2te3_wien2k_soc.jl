# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Bi2Te3 Wien2k SOC workflow
#
# Interpolation and transport calculation for bismuth telluride (Bi2Te3)
# using Wien2k output with spin-orbit coupling (case.struct, case.energyso, case.scf).
#
# DFT conditions:
#   Code: Wien2k (FP-LAPW + SOC)
#   Spin-orbit coupling: yes (case.energyso)
#   dosweight: 1 (SOC includes both spin channels)
#   k-grid: 4960 IBZ k-points
#   Lattice: a = 4.386 Angstrom (8.288 bohr), c = 30.50 Angstrom (57.63 bohr)
#   Space group: 166 R-3m (rhombohedral)
#   Source: BoltzTraP2 distribution (Bi2Te3)
#
# Data: Bi2Te3 Wien2k data is NOT bundled due to size (33 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Bi2Te3
#       Place the downloaded directory as examples/data/Bi2Te3/
#
# Usage:
#   julia --project=. examples/150_bi2te3_wien2k_soc.jl

using BoltzTraP

# Path to Wien2k data
datadir = joinpath(@__DIR__, "data", "Bi2Te3")

if !isdir(datadir)
    error("""Bi2Te3 Wien2k data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Bi2Te3
    Place as: examples/data/Bi2Te3/""")
end

# Step 1: Load Wien2k data explicitly
println("Loading Wien2k data (SOC)...")
data = load_wien2k(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="Wien2k", kpoints=5000, dosweight=1.0, verbose=true)

# Step 3: Compute transport
# Bi2Te3 has 174 bands spanning a huge energy range (~80 eV).
# The default auto mu range would be too sparse for plotting near Fermi level.
# Use explicit mu range: Fermi ± 0.1 Ha (≈ ±2.7 eV), 200 points.
println("\nComputing transport coefficients...")
fermi_Ha = interp.metadata["fermi"]
mur = collect(range(fermi_Ha - 0.1, fermi_Ha + 0.1, length=200))
transport = run_integrate(interp, mur; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")

# Step 4: Plot transport coefficients (S, sigma, kappa vs mu at T=300K)
using Plots
fermi_dft_eV = transport.metadata["fermi_dft_eV"]
mu = transport.mu_values .- fermi_dft_eV
iT = findfirst(t -> isapprox(t, 300.0, atol=1.0), transport.temperatures)
p1 = plot(mu, transport.seebeck[1,1,iT,:] .* 1e6;
    title="T = 300 K", ylabel="S_xx (μV/K)", ylims=(-300, 300),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p2 = plot(mu, transport.sigma[1,1,iT,:];
    ylabel="σ_xx/τ (S/m)", yscale=:log10, ylims=(1e21, 1e25),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p3 = plot(mu, transport.kappa[1,1,iT,:];
    ylabel="κ_xx/τ (W/m/K)", yscale=:log10, ylims=(1e16, 1e20),
    xlabel="μ - E_F (eV)", legend=false, linewidth=2)
p = plot(p1, p2, p3; layout=(3, 1), size=(700, 900), xlims=(-0.3, 0.3), left_margin=5Plots.mm)
output_png = joinpath(@__DIR__, "150_transport_Bi2Te3_Wien2k_SOC_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
