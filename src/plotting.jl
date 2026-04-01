# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Plotting stubs
#
# Implementations moved to ext/BoltzTraPPlotsExt.jl

"""
    plot_bands(result::InterpolationResult; kwargs...)
    plot_bands(file::String; kwargs...)

Plot band structure. Requires `using Plots`.
"""
function plot_bands end

"""
    plot_transport(result::TransportResult; kwargs...)
    plot_transport(file::String; kwargs...)

Plot transport properties. Requires `using Plots`.
"""
function plot_transport end

"""
    savefig_publication(data, path; kwargs...)

Save publication-quality figure. Requires `using PythonPlot`.
"""
function savefig_publication end
