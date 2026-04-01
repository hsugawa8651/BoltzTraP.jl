# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Plotting data types

"""
    BandPlotData

Precomputed band structure data for plotting, independent of plotting library.

# Fields
- `kdist::Vector{Float64}` — cumulative k-path distances
- `ebands_eV::Matrix{Float64}` — band energies in eV relative to Fermi (nbands × nk)
- `unique_labels::Vector{String}` — high-symmetry point labels (merged with "|")
- `unique_positions::Vector{Float64}` — high-symmetry point positions
- `emin::Float64` — display energy minimum (eV)
- `emax::Float64` — display energy maximum (eV)
- `fermi_Ha::Float64` — Fermi energy in Hartree
- `title::String` — plot title
"""
struct BandPlotData
    kdist::Vector{Float64}
    ebands_eV::Matrix{Float64}
    unique_labels::Vector{String}
    unique_positions::Vector{Float64}
    emin::Float64
    emax::Float64
    fermi_Ha::Float64
    title::String
end

"""
    TransportPlotData

Precomputed transport property data for plotting, independent of plotting library.

# Fields
- `xdata::Vector{Float64}` — x-axis data (mu-Ef in eV, or temperature in K)
- `ydata::Vector{Float64}` — y-axis data (already unit-converted)
- `xlabel::String` — x-axis label
- `ylabel::String` — y-axis label
- `title::String` — plot title
"""
struct TransportPlotData
    xdata::Vector{Float64}
    ydata::Vector{Float64}
    xlabel::String
    ylabel::String
    title::String
end
