# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: CoSb3 Wien2k workflow
#
# Interpolation and transport calculation for cobalt antimonide (CoSb3)
# using Wien2k output. CoSb3 is a skutterudite thermoelectric material
# with 16 atoms per unit cell.
#
# DFT conditions:
#   Code: Wien2k (FP-LAPW)
#   XC: GGA (from BoltzTraP2 distribution)
#   Spin-polarization: No (non-magnetic, dosweight=2)
#   k-grid: 1030 IBZ k-points
#   Space group: 204 Im-3 (cubic skutterudite)
#   Source: BoltzTraP2 distribution (CoSb3)
#
# Data: CoSb3 Wien2k data is NOT bundled due to size (18 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/CoSb3
#       Place the downloaded directory as examples/data/CoSb3/
#
# Usage:
#   julia --project=. examples/160_cosb3_wien2k.jl

using BoltzTraP

# Path to Wien2k data
datadir = joinpath(@__DIR__, "data", "CoSb3")

if !isdir(datadir)
    error("""CoSb3 Wien2k data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/CoSb3
    Place as: examples/data/CoSb3/""")
end

# Step 1: Load Wien2k data
println("Loading Wien2k data...")
data = load_wien2k(datadir)

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="Wien2k", kpoints=5000, verbose=true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

# Step 3: Compute transport
# CoSb3 has 174 bands spanning a huge energy range.
# Use explicit mu range: Fermi ± 0.02 Ha (≈ ±0.54 eV), 200 points.
println("\nComputing transport coefficients...")
fermi_Ha = interp.metadata["fermi"]
mur = collect(range(fermi_Ha - 0.02, fermi_Ha + 0.02, length=200))
transport = run_integrate(interp, mur; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")

# Step 4: Plot transport coefficients (S, sigma, kappa vs mu at T=300K)
using Plots
fermi_dft_eV = transport.metadata["fermi_dft_eV"]
mu = transport.mu_values .- fermi_dft_eV
iT = findfirst(t -> isapprox(t, 300.0, atol=1.0), transport.temperatures)
p1 = plot(mu, transport.seebeck[1,1,iT,:] .* 1e6;
    title="T = 300 K", ylabel="S_xx (μV/K)", ylims=(-500, 500),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p2 = plot(mu, transport.sigma[1,1,iT,:];
    ylabel="σ_xx/τ (S/m)", yscale=:log10, ylims=(1e15, 1e21),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p3 = plot(mu, transport.kappa[1,1,iT,:];
    ylabel="κ_xx/τ (W/m/K)", yscale=:log10, ylims=(1e11, 1e16),
    xlabel="μ - E_F (eV)", legend=false, linewidth=2)
p = plot(p1, p2, p3; layout=(3, 1), size=(700, 900), xlims=(-0.5, 0.5), left_margin=5Plots.mm)
output_png = joinpath(@__DIR__, stem * "_transport_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
