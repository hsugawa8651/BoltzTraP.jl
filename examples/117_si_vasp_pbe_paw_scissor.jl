# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si VASP PBE-PAW scissor operator (Eg=1.17 eV)
#
# Demonstrates the scissor correction for silicon.
# DFT (LDA/GGA) underestimates the Si band gap (~0.5 eV);
# the scissor operator shifts conduction bands to match the
# experimental gap (1.17 eV).
#
# DFT conditions (same as 110_si_vasp_pbe_paw.jl):
#   Code: VASP (PAW_PBE Si 05Jan2001)
#   XC functional: PBE-GGA
#   k-grid: 17×17×17 (165 IBZ k-points)
#   Lattice constant: a = 5.467 Å (10.33 bohr)
#   Source: BoltzTraP2 distribution (Si.vasp)
#
# Usage:
#   julia --project=. examples/117_si_vasp_pbe_paw_scissor.jl

using BoltzTraP

# Path to Si VASP data
datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")

# Step 1: Interpolate band structure
println("Step 1: Interpolating band structure...")
interp = run_interpolate(datadir; kpoints=5000, verbose=true)

# Step 2: Transport WITHOUT scissor correction
println("\nStep 2: Transport without scissor correction...")
transport_nosci = run_integrate(interp; temperatures=[300.0], verbose=true)

# Step 3: Transport WITH scissor correction (target gap = 1.17 eV)
println("\nStep 3: Transport with scissor correction (gap = 1.17 eV)...")
transport_sci = run_integrate(interp; temperatures=[300.0], scissor=1.17, verbose=true)

# Step 4: Compare results
println("\nComparison at T = 300 K:")
println("  Without scissor:")
println("    Chemical potential range: $(round(transport_nosci.mu_values[1], digits=4)) .. $(round(transport_nosci.mu_values[end], digits=4)) eV")
println("    Number of mu points: $(length(transport_nosci.mu_values))")
println("  With scissor (1.17 eV):")
println("    Chemical potential range: $(round(transport_sci.mu_values[1], digits=4)) .. $(round(transport_sci.mu_values[end], digits=4)) eV")
println("    Number of mu points: $(length(transport_sci.mu_values))")

# Step 5: Plot transport comparison (S, σ, κ vs μ at T=300K)
using Plots
fermi_dft_eV = transport_nosci.metadata["fermi_dft_eV"]
mu_nosci = transport_nosci.mu_values .- fermi_dft_eV
mu_sci = transport_sci.mu_values .- transport_sci.metadata["fermi_dft_eV"]
iT = 1  # only T=300K
p1 = plot(mu_nosci, transport_nosci.seebeck[1,1,iT,:] .* 1e6;
    title="T = 300 K", ylabel="S_xx (μV/K)", ylims=(-1100, 1100),
    xlabel="", xformatter=_->"", label="no scissor", linewidth=2, color=:blue)
plot!(p1, mu_sci, transport_sci.seebeck[1,1,iT,:] .* 1e6;
    label="scissor (1.17 eV)", linewidth=2, color=:red)
p2 = plot(mu_nosci, transport_nosci.sigma[1,1,iT,:];
    ylabel="σ_xx/τ (S/m)", yscale=:log10, ylims=(1e14, 1e21),
    xlabel="", xformatter=_->"", legend=false, linewidth=2, color=:blue)
plot!(p2, mu_sci, transport_sci.sigma[1,1,iT,:];
    linewidth=2, color=:red)
p3 = plot(mu_nosci, transport_nosci.kappa[1,1,iT,:];
    ylabel="κ_xx/τ (W/m/K)", yscale=:log10, ylims=(1e10, 1e16),
    xlabel="μ - E_F (eV)", legend=false, linewidth=2, color=:blue)
plot!(p3, mu_sci, transport_sci.kappa[1,1,iT,:];
    linewidth=2, color=:red)
p = plot(p1, p2, p3; layout=(3, 1), size=(700, 900), xlims=(-1.0, 1.5), left_margin=5Plots.mm)
output_png = joinpath(@__DIR__, "117_transport_Si_VASP_PBE-PAW_scissor_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
