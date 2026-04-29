# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: Si VASP PBE-PAW electronic heat capacity
#
# Computes the electronic contribution to the heat capacity C_v
# for silicon at multiple temperatures.
#
# DFT conditions (same as 110_si_vasp_pbe_paw.jl):
#   Code: VASP (PAW_PBE Si 05Jan2001)
#   XC functional: PBE-GGA
#   k-grid: 17×17×17 (165 IBZ k-points)
#   Lattice constant: a = 5.467 Å (10.33 bohr)
#   Source: BoltzTraP2 distribution (Si.vasp)
#
# Usage:
#   julia --project=. examples/118_si_vasp_pbe_paw_cv.jl

using BoltzTraP

# Path to Si VASP data
datadir = joinpath(@__DIR__, "..", "benchmarks", "data", "Si.vasp")

# Output file stem
stem = splitext(basename(@__FILE__))[1]

# Step 1: Interpolate band structure
println("Step 1: Interpolating band structure...")
interp = run_interpolate(datadir; kpoints = 10000, verbose = true)
save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

# Step 2: Compute transport (DOS is needed for C_v)
println("\nStep 2: Computing transport coefficients...")
temperatures = [100.0, 200.0, 300.0, 400.0, 500.0]
transport = run_integrate(interp; temperatures = temperatures, verbose = true)
save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

# Step 3: Compute electronic heat capacity
println("\nStep 3: Computing electronic heat capacity...")
cv = calc_cv(transport)  # Returns (nT, nmu) matrix in J/K

# Step 4: Print results
println("\nResults: C_v matrix shape = $(size(cv)) (nT × nμ)")

# At the Fermi level (in the gap for Si → C_v ≈ 0)
fermi_eV = transport.metadata["fermi_eV"]
mu_idx = argmin(abs.(transport.mu_values .- fermi_eV))
println("\nC_v at Fermi level (mu ≈ $(round(transport.mu_values[mu_idx], digits=2)) eV):")
println("  Si is a semiconductor → C_v ≈ 0 in the gap")

# Maximum C_v across all mu (occurs near band edges)
for (i, T) in enumerate(temperatures)
    cv_max = maximum(cv[i, :])
    mu_max = transport.mu_values[argmax(cv[i, :])]
    println(
        "  T = $(Int(T)) K:  max C_v = $(round(cv_max, sigdigits=4)) J/K  at mu = $(round(mu_max, digits=2)) eV",
    )
end

# Step 5: Plot C_v vs μ at multiple temperatures
using Plots
fermi_dft_eV = transport.metadata["fermi_dft_eV"]
mu = transport.mu_values .- fermi_dft_eV
p = plot(;
    xlabel = "μ - E_F (eV)",
    ylabel = "C_v (10⁻²⁴ J/K)",
    xlims = (-7, -4),
    legend = :topright,
    linewidth = 2,
    size = (700, 500),
    left_margin = 5Plots.mm,
)
for (i, T) in enumerate(temperatures)
    plot!(p, mu, cv[i, :] .* 1e24; label = "$(Int(T)) K", linewidth = 2)
end
output_png = joinpath(@__DIR__, stem * ".png")
savefig(p, output_png)
println("Saved plot to $output_png")
