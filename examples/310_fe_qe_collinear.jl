# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Fe QE collinear magnetic workflow
#
# Interpolation and transport calculation for iron (Fe, BCC metal)
# using Quantum ESPRESSO output with collinear spin polarization.
#
# DFT conditions:
#   Code: Quantum ESPRESSO (collinear, nspin=2)
#   XC functional: PBE-GGA
#   k-grid: 47 IBZ k-points
#   Lattice: a = 5.217 bohr, BCC Im-3m
#   Source: BoltzTraP2 distribution (Fe.ESPRESSO.collinear)
#
# Usage:
#   julia --project=. examples/310_fe_qe_collinear.jl

using BoltzTraP
using BoltzTraP: HA_TO_EV

# Path to Fe QE collinear data
datadir = joinpath(@__DIR__, "data", "Fe.ESPRESSO.collinear")

if !isdir(datadir)
    error("Fe QE data not found at $datadir")
end

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 1: Load QE data and interpolate band structure
println("Loading QE data (collinear)...")
data = BoltzTraP.load_qe(datadir)

println("Interpolating band structure...")
interp = run_interpolate(data; source = "QE", kpoints = 5000, verbose = true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

println("  Equivalence classes: $(length(interp.equivalences))")
println("  Bands: $(size(interp.coeffs, 1))")

# Step 2: Compute transport coefficients
# Fe has many bands spanning a wide energy range.
# Use explicit mu range: Fermi ± 2 eV, 200 points for fine resolution.
println("\nComputing transport coefficients...")
fermi_ha = interp.metadata["fermi"]
mur = collect(range(fermi_ha - 2.0 / HA_TO_EV, fermi_ha + 2.0 / HA_TO_EV, length = 200))
temperatures = [200.0, 300.0, 400.0, 1000.0]
transport = run_integrate(interp, mur; temperatures = temperatures, verbose = true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

println("\nResults:")
println("  Temperatures: $(transport.temperatures) K")
println("  Chemical potential points: $(length(transport.mu_values))")
println("  Tensor shape (sigma): $(size(transport.sigma))")

# Step 3: Plot transport coefficients (S, sigma, kappa vs mu at T=300K)
using Plots
fermi_dft_eV = transport.metadata["fermi_dft_eV"]
mu = transport.mu_values .- fermi_dft_eV
iT = findfirst(t -> isapprox(t, 300.0, atol = 1.0), transport.temperatures)
p1 = plot(
    mu,
    transport.seebeck[1, 1, iT, :] .* 1e6;
    title = "T = 300 K",
    ylabel = "S_xx (μV/K)",
    ylims = (-500, 500),
    xlabel = "",
    xformatter = _ -> "",
    legend = false,
    linewidth = 2,
)
p2 = plot(
    mu,
    transport.sigma[1, 1, iT, :];
    ylabel = "σ_xx/τ (S/m)",
    yscale = :log10,
    ylims = (1e18, 1e22),
    xlabel = "",
    xformatter = _ -> "",
    legend = false,
    linewidth = 2,
)
p3 = plot(
    mu,
    transport.kappa[1, 1, iT, :];
    ylabel = "κ_xx/τ (W/m/K)",
    yscale = :log10,
    ylims = (1e13, 1e17),
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
