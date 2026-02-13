# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Li VASP PBE-PAW workflow
#
# Interpolation and transport calculation for lithium (Li, BCC metal)
# using VASP output (vasprun.xml + POSCAR).
#
# DFT conditions:
#   Code: VASP (PAW_PBE Li_sv 10Sep2004)
#   XC functional: PBE-GGA
#   Spin: ISPIN=2 (spin-polarized)
#   k-grid: 21x21x21 (455 IBZ k-points)
#   Cutoff: default (PREC=High)
#   Lattice constant: a = 3.420 Angstrom (6.46 bohr), BCC Im-3m
#   Source: BoltzTraP2 distribution (Li.vasp)
#
# Data: Li VASP data is NOT bundled due to size (9.8 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Li.vasp
#       Place the downloaded directory as examples/data/Li.vasp/
#
# Usage:
#   julia --project=. examples/130_li_vasp_pbe_paw.jl

using BoltzTraP

# Path to Li VASP data
datadir = joinpath(@__DIR__, "data", "Li.vasp")

if !isdir(datadir)
    error("""Li VASP data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Li.vasp
    Place as: examples/data/Li.vasp/""")
end

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 1: Interpolate band structure
println("Step 1: Interpolating band structure...")
interp = run_interpolate(datadir; kpoints=5000, verbose=true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

println("  Equivalence classes: $(length(interp.equivalences))")
println("  Bands: $(size(interp.coeffs, 1))")

# Step 2: Compute transport coefficients
println("\nStep 2: Computing transport coefficients...")
temperatures = [200.0, 300.0, 400.0, 500.0]
transport = run_integrate(interp; temperatures=temperatures, verbose=true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

# Step 3: Print results
println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
println("  Tensor shape (sigma): $(size(transport.sigma))")

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
