# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - v0.2 Performance benchmark

using BoltzTraP
using BenchmarkTools
using NPZ

println("=== v0.2 Performance Benchmark ===\n")

# Si data load from reftest data
data_path = joinpath(@__DIR__, "..", "reftest", "data", "si_interpolation.npz")
println("Loading Si data from: $data_path")
si_data = npzread(data_path)
println("Data loaded successfully\n")

# Extract data for benchmarks
lattvec = si_data["lattvec"]
n_eq = si_data["n_equivalences"]
equivalences = [si_data["equiv_$i"] for i = 0:(n_eq-1)]
coeffs = si_data["coeffs"]

# 1. getBTPbands benchmark
println("1. getBTPbands (vvband symmetry optimization)")
b1 = @benchmark BoltzTraP.getBTPbands($coeffs, $equivalences, $lattvec)
display(b1)
println()

# Get eband and vvband for subsequent benchmarks
eband, vvband = BoltzTraP.getBTPbands(coeffs, equivalences, lattvec)

# 2. BTPDOS benchmark
println("2. BTPDOS (@view optimization)")
b2 = @benchmark BoltzTraP.BTPDOS($eband, $vvband)
display(b2)
println()

# 3. calc_onsager_coefficients benchmark
println("3. calc_onsager_coefficients (@view + multiply optimization)")
epsilon, dos, vvdos = BoltzTraP.BTPDOS(eband, vvband)
T_range = [300.0]
μ_range = collect(-0.5:0.01:0.5)
vuc = 100.0

# Run fermi_integrals first
N, L0, L1, L2 = BoltzTraP.fermi_integrals(epsilon, dos, vvdos, μ_range, T_range)

b3 = @benchmark BoltzTraP.calc_onsager_coefficients($L0, $L1, $L2, $T_range, $vuc)
display(b3)
println()

println("=== Benchmark Complete ===")
