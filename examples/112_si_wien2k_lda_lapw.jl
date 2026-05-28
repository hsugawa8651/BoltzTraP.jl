# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si Wien2k LDA-LAPW workflow
#
# Interpolation and transport calculation for silicon
# using Wien2k output (case.struct, case.energy, case.scf).
#
# DFT conditions:
#   Code: Wien2k (FP-LAPW)
#   XC functional: LDA
#   k-grid: 25×25×25 (455 IBZ k-points)
#   Lattice constant: a = 5.451 Å (10.30 bohr)
#   Source: BoltzTraP2 distribution (Si)
#
# Data: Wien2k Si data is NOT bundled due to size (15 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Si
#       Place the downloaded "Si" directory as examples/data/Si.wien2k/
#
# Usage:
#   julia --project=. examples/112_si_wien2k_lda_lapw.jl

using BoltzTraP

# Path to Wien2k data
datadir = joinpath(@__DIR__, "data", "Si.wien2k")

if !isdir(datadir)
    error("""Wien2k Si data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Si
    Place as: examples/data/Si.wien2k/""")
end

# Step 1: Load Wien2k data explicitly
println("Loading Wien2k data...")
data = load_wien2k(datadir)

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source = "Wien2k", kpoints = 5000, verbose = true)
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
    xformatter = _ -> "",
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
    xformatter = _ -> "",
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
