# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Li GENE PBE-PAW (BoltzTraP1) workflow
#
# Interpolation and transport calculation for lithium (Li, BCC metal)
# using GENE/generic format (case.structure, case.energy).
#
# DFT conditions (converted from VASP data):
#   Code: VASP → GENE format (PAW_PBE Li_sv 10Sep2004)
#   XC functional: PBE-GGA
#   k-grid: 21×21×21 (455 IBZ k-points)
#   Lattice constant: a = 3.420 Angstrom (6.46 bohr), BCC Im-3m
#   Source: same VASP calculation as 130_li_vasp_pbe_paw.jl
#
# Note: GENE format does not store spin information.
#       This data was converted from a spin-polarized VASP calculation,
#       but is loaded as non-magnetic (DFTData{1}, dosweight=2).
#
# Usage:
#   julia --project=. examples/133_li_gene_pbe_paw.jl

using BoltzTraP

# Path to GENE data
datadir = joinpath(@__DIR__, "data", "Li.GENE")

# Step 1: Load GENE data explicitly
println("Loading GENE data...")
data = load_gene(datadir)

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="GENE", kpoints=5000, verbose=true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)
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
    title="T = 300 K", ylabel="S_xx (μV/K)", ylims=(-50, 50),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p2 = plot(mu, transport.sigma[1,1,iT,:];
    ylabel="σ_xx/τ (S/m)", yscale=:log10, ylims=(1e18, 1e22),
    xlabel="", xformatter=_->"", legend=false, linewidth=2)
p3 = plot(mu, transport.kappa[1,1,iT,:];
    ylabel="κ_xx/τ (W/m/K)", yscale=:log10, ylims=(1e13, 1e17),
    xlabel="μ - E_F (eV)", legend=false, linewidth=2)
p = plot(p1, p2, p3; layout=(3, 1), size=(700, 900), xlims=(-0.5, 0.5), left_margin=5Plots.mm)
output_png = joinpath(@__DIR__, stem * "_transport_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
