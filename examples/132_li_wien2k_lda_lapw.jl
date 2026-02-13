# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Li Wien2k LDA-LAPW workflow
#
# Interpolation and transport calculation for lithium (Li, BCC metal)
# using Wien2k output (case.struct, case.energy, case.scf).
#
# DFT conditions:
#   Code: Wien2k (FP-LAPW)
#   XC functional: LDA
#   Spin: non-magnetic (dosweight=2)
#   k-grid: 220 IBZ k-points
#   Lattice constant: a = 3.378 Angstrom (6.384 bohr), BCC Im-3m
#   Source: BoltzTraP2 distribution (Li.W2K)
#
# Data: Li Wien2k data is NOT bundled due to size (2.1 MB).
#       Download from BoltzTraP2:
#       https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Li.W2K
#       Place the downloaded directory as examples/data/Li.W2K/
#
# Usage:
#   julia --project=. examples/132_li_wien2k_lda_lapw.jl

using BoltzTraP

# Path to Wien2k data
datadir = joinpath(@__DIR__, "data", "Li.W2K")

if !isdir(datadir)
    error("""Li Wien2k data not found at $datadir
    Download from: https://gitlab.com/sousaw/BoltzTraP2/-/tree/master/data/Li.W2K
    Place as: examples/data/Li.W2K/""")
end

# Step 1: Load Wien2k data explicitly
println("Loading Wien2k data...")
data = load_wien2k(datadir)

# Step 2: Interpolate
println("Interpolating band structure...")
interp = run_interpolate(data; source="Wien2k", kpoints=5000, verbose=true)

# Step 3: Compute transport
println("\nComputing transport coefficients...")
transport = run_integrate(interp; temperatures=[200.0, 300.0, 400.0, 500.0], verbose=true)

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
output_png = joinpath(@__DIR__, "132_transport_Li_Wien2k_LDA-LAPW_300K.png")
savefig(p, output_png)
println("Saved plot to $output_png")
