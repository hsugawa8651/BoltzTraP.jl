# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si ABINIT PBE-PAW workflow
#
# Interpolation and transport calculation for silicon
# using ABINIT output (GSR.nc NetCDF file).
#
# DFT conditions:
#   Code: ABINIT (PAW)
#   XC functional: PBE-GGA (ixc = -101130)
#   k-grid: ~16×16×16 (408 IBZ k-points)
#   Cutoff: ecut = 20 Ha
#   Lattice constant: a = 5.469 Å (10.33 bohr)
#   Source: independent ABINIT calculation
#
# Note: Requires NCDatasets.jl package for reading NetCDF files.
#   using Pkg; Pkg.add("NCDatasets")
#
# Usage:
#   julia --project=. examples/114_si_abinit_pbe_paw.jl

using NCDatasets  # Required before loading BoltzTraP for ABINIT support
using BoltzTraP

# Path to ABINIT data (directory containing *_GSR.nc)
datadir = joinpath(@__DIR__, "data", "Si.abinit")

# Step 1: Load ABINIT data explicitly
println("Loading ABINIT data...")
data = load_abinit(datadir)

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source = "ABINIT", kpoints = 5000, verbose = true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport =
    run_integrate(interp; temperatures = [200.0, 300.0, 400.0, 500.0], verbose = true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")

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
