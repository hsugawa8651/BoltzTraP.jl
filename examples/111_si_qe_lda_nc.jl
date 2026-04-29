# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si QE LDA-NC workflow
#
# Interpolation and transport calculation for silicon
# using Quantum ESPRESSO output (data-file-schema.xml).
#
# DFT conditions:
#   Code: Quantum ESPRESSO (Si.pz-vbc.UPF, norm-conserving)
#   XC functional: LDA (Perdew-Zunger)
#   k-grid: 11×11×11 (56 IBZ k-points)
#   Cutoff: ecutwfc = 50 Ry, ecutrho = 400 Ry
#   Lattice constant: a = 5.402 Å (10.21 bohr)
#   Source: QE tutorial (Pranab Das, CC-BY 4.0)
#
# Usage:
#   julia --project=. examples/111_si_qe_lda_nc.jl

using BoltzTraP

# Path to QE data (directory containing silicon.xml)
datadir = joinpath(@__DIR__, "data", "Si.ESPRESSO", "out")

# Step 1: Load QE data explicitly
println("Loading QE data...")
data = load_qe(datadir)

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source = "QE", kpoints = 5000, verbose = true)
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
