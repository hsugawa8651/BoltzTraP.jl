# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - RecipesBase extension

module BoltzTraPRecipesBaseExt

using RecipesBase
using BoltzTraP: BandPlotData, TransportPlotData

@recipe function f(bpd::BandPlotData)
    xlabel --> ""
    ylabel --> "Energy (eV)"
    title --> bpd.title
    legend --> false
    ylims --> (bpd.emin, bpd.emax)
    xticks --> (bpd.unique_positions, bpd.unique_labels)
    framestyle --> :box

    # Band lines
    nbands = size(bpd.ebands_eV, 1)
    for b = 1:nbands
        @series begin
            seriestype := :path
            linecolor --> :black
            linewidth --> 1.5
            label := ""
            bpd.kdist, bpd.ebands_eV[b, :]
        end
    end

    # Fermi level
    @series begin
        seriestype := :hline
        primary := false
        linecolor := :red
        linewidth := 0.5
        linestyle := :dash
        label := ""
        [0.0]
    end

    # Vertical lines at high-symmetry points
    for pos in bpd.unique_positions
        @series begin
            seriestype := :vline
            primary := false
            linecolor := :gray
            linewidth := 0.5
            linestyle := :dot
            label := ""
            [pos]
        end
    end
end

@recipe function f(tpd::TransportPlotData)
    xlabel --> tpd.xlabel
    ylabel --> tpd.ylabel
    title --> tpd.title
    legend --> false
    linewidth --> 2

    tpd.xdata, tpd.ydata
end

end # module
