# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si VASP PBE-PAW basic workflow
#
# Basic workflow: interpolation and transport coefficient calculation
# for silicon using VASP output (vasprun.xml + POSCAR).
#
# DFT conditions:
#   Code: VASP (PAW_PBE Si 05Jan2001)
#   XC functional: PBE-GGA
#   k-grid: 17×17×17 (165 IBZ k-points)
#   Cutoff: ENMAX = 306.7 eV
#   Lattice constant: a = 5.467 Å (10.33 bohr)
#   Source: BoltzTraP2 distribution (Si.vasp)
#
# Usage:
#   julia --project=. examples/110_si_vasp_pbe_paw.jl
#
# CLI equivalent:
#   boltztrap interpolate -k 5000 -o si_interp.jld2 benchmarks/data/Si.vasp
#   boltztrap integrate si_interp.jld2 -t 200:100:500 -o si_transport.jld2

using BoltzTraP

# Path to Si VASP data (directory containing vasprun.xml and POSCAR)
datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 1: Interpolate band structure
println("Step 1: Interpolating band structure...")
interp = run_interpolate(datadir; kpoints = 5000, verbose = true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

println("  Equivalence classes: $(length(interp.equivalences))")
println("  Bands: $(size(interp.coeffs, 1))")

# Step 2: Compute transport coefficients
println("\nStep 2: Computing transport coefficients...")
temperatures = [200.0, 300.0, 400.0, 500.0]
transport = run_integrate(interp; temperatures = temperatures, verbose = true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

# Step 3: Print results
println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
println("  Tensor shape (σ): $(size(transport.sigma))")

# Step 4: Plot transport coefficients (S, σ, κ vs μ at T=300K)
using Plots
fermi_dft_eV = transport.metadata["fermi_dft_eV"]
mu = transport.mu_values .- fermi_dft_eV
iT = findfirst(t -> isapprox(t, 300.0, atol = 1.0), transport.temperatures)
p1 = plot(
    mu,
    transport.seebeck[1, 1, iT, :] .* 1e6;
    title = "T = 300 K",
    ylabel = "S_xx (μV/K)",
    ylims = (-1100, 1100),
    xlabel = "",
    xformatter = _->"",
    legend = false,
    linewidth = 2,
)
p2 = plot(
    mu,
    transport.sigma[1, 1, iT, :];
    ylabel = "σ_xx/τ (S/m)",
    yscale = :log10,
    ylims = (1e14, 1e21),
    xlabel = "",
    xformatter = _->"",
    legend = false,
    linewidth = 2,
)
p3 = plot(
    mu,
    transport.kappa[1, 1, iT, :];
    ylabel = "κ_xx/τ (W/m/K)",
    yscale = :log10,
    ylims = (1e10, 1e16),
    xlabel = "μ - E_F (eV)",
    legend = false,
    linewidth = 2,
)
p = plot(
    p1,
    p2,
    p3;
    layout = (3, 1),
    size = (700, 900),
    xlims = (-0.5, 0.5),
    left_margin = 5Plots.mm,
)
output_png = joinpath(@__DIR__, stem * "_transport_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
