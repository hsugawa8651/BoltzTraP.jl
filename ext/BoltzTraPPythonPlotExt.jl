# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - PythonPlot extension
#
# Convention: do NOT use `PythonPlot.subplots()` — it returns a raw
# `PythonCall.Py` whose `PythonPlot.close(fig)` triggers a `MethodError`.
# Always create figures via `PythonPlot.figure(figsize=...)` and add axes
# with `fig.add_axes([l, b, w, h])` (or `fig.add_subplot(...)`).

module BoltzTraPPythonPlotExt

using BoltzTraP
import BoltzTraP: BandPlotData, TransportPlotData, savefig_publication
import PythonPlot

const CM_PER_INCH = 2.54

# ── Helper: layout ──

function _layout_axes(axis_width_cm, axis_height_cm, n;
        margin_left_cm=1.5, margin_right_cm=0.3,
        margin_bottom_cm=1.0, margin_top_cm=0.8,
        hgap_cm=1.8, vgap_cm=1.5,
        nrows=1, ncols=1)
    widths = axis_width_cm isa AbstractVector ? axis_width_cm : fill(axis_width_cm, ncols)
    heights = axis_height_cm isa AbstractVector ? axis_height_cm : fill(axis_height_cm, nrows)

    fig_w_cm = margin_left_cm + sum(widths) + hgap_cm * (ncols - 1) + margin_right_cm
    fig_h_cm = margin_bottom_cm + sum(heights) + vgap_cm * (nrows - 1) + margin_top_cm

    fig_w = fig_w_cm / CM_PER_INCH
    fig_h = fig_h_cm / CM_PER_INCH

    positions = Vector{NTuple{4,Float64}}()
    for row in 1:nrows
        for col in 1:ncols
            left = (margin_left_cm + sum(widths[1:col-1]) + hgap_cm * (col - 1)) / fig_w_cm
            bottom = (margin_bottom_cm + sum(heights[row+1:end]) + vgap_cm * (nrows - row)) / fig_h_cm
            w = widths[col] / fig_w_cm
            h = heights[row] / fig_h_cm
            push!(positions, (left, bottom, w, h))
        end
    end
    return fig_w, fig_h, positions
end

# ── Helper: render BandPlotData onto an Axes ──

function _plot_bands_on_ax!(ax, bpd::BandPlotData;
        color="black", linewidth=1.0, linestyle="-",
        title::AbstractString="")
    nbands = size(bpd.ebands_eV, 1)
    for i in 1:nbands
        ax.plot(bpd.kdist, bpd.ebands_eV[i, :];
            color, linewidth, linestyle)
    end

    # Fermi level
    ax.axhline(0.0; color="gray", linewidth=0.5, linestyle="--")

    # High-symmetry vertical lines
    for pos in bpd.unique_positions
        ax.axvline(pos; color="gray", linewidth=0.5)
    end

    ax.set_xticks(bpd.unique_positions)
    ax.set_xticklabels(bpd.unique_labels)
    ax.set_ylabel("E - E_F (eV)")
    ax.set_xlim(bpd.kdist[1], bpd.kdist[end])
    ax.set_ylim(bpd.emin, bpd.emax)

    eff_title = isempty(title) ? bpd.title : title
    isempty(eff_title) || ax.set_title(eff_title)
end

# ── Helper: render TransportPlotData onto an Axes ──

function _plot_transport_on_ax!(ax, tpd::TransportPlotData;
        color="black", linewidth=1.5, linestyle="-",
        title::AbstractString="")
    ax.plot(tpd.xdata, tpd.ydata; color, linewidth, linestyle)
    ax.set_xlabel(tpd.xlabel)
    ax.set_ylabel(tpd.ylabel)
    eff_title = isempty(title) ? tpd.title : title
    isempty(eff_title) || ax.set_title(eff_title)
end

# ── Default styles for overlay ──

const _OVERLAY_STYLES = [
    (color="black", linewidth=1.5, linestyle="-"),
    (color="blue", linewidth=1.2, linestyle="--"),
    (color="red", linewidth=1.2, linestyle="-."),
    (color="green", linewidth=1.2, linestyle=":"),
    (color="purple", linewidth=1.2, linestyle="--"),
]

# ── Single BandPlotData ──

function BoltzTraP.savefig_publication(
        bpd::BandPlotData, path::AbstractString;
        axis_width_cm=8.0, axis_height_cm=6.0,
        ylims=nothing,
        kwargs...)

    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
    fig = PythonPlot.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes(collect(positions[1]))

    _plot_bands_on_ax!(ax, bpd; kwargs...)
    isnothing(ylims) || ax.set_ylim(ylims...)

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

# ── Vector{BandPlotData} (subplot / overlay) ──

function BoltzTraP.savefig_publication(
        bpds::AbstractVector{<:BandPlotData}, path::AbstractString;
        axis_width_cm=8.0, axis_height_cm=6.0,
        layout=(1, length(bpds)),
        overlay::Bool=false,
        ylims=nothing,
        title::AbstractString="",
        kwargs...)

    if overlay
        fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
        fig = PythonPlot.figure(figsize=(fig_w, fig_h))
        ax = fig.add_axes(collect(positions[1]))
        for (i, bpd) in enumerate(bpds)
            style = _OVERLAY_STYLES[mod1(i, length(_OVERLAY_STYLES))]
            _plot_bands_on_ax!(ax, bpd;
                color=style.color, linewidth=style.linewidth, linestyle=style.linestyle)
        end
        isnothing(ylims) || ax.set_ylim(ylims...)
        isempty(title) || ax.set_title(title)
    else
        nrows, ncols = layout
        fig_w, fig_h, positions = _layout_axes(
            axis_width_cm, axis_height_cm, length(bpds); nrows, ncols)
        fig = PythonPlot.figure(figsize=(fig_w, fig_h))
        for (i, bpd) in enumerate(bpds)
            i > length(positions) && break
            ax = fig.add_axes(collect(positions[i]))
            _plot_bands_on_ax!(ax, bpd; kwargs...)
            isnothing(ylims) || ax.set_ylim(ylims...)
        end
    end

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

# ── Single TransportPlotData ──

function BoltzTraP.savefig_publication(
        tpd::TransportPlotData, path::AbstractString;
        axis_width_cm=8.0, axis_height_cm=6.0,
        ylims=nothing,
        kwargs...)

    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
    fig = PythonPlot.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes(collect(positions[1]))

    _plot_transport_on_ax!(ax, tpd; kwargs...)
    isnothing(ylims) || ax.set_ylim(ylims...)

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

# ── Vector{TransportPlotData} (subplot only; overlay deferred to v0.3.3+
# because it requires unit compatibility checks across panels) ──

function BoltzTraP.savefig_publication(
        tpds::AbstractVector{<:TransportPlotData}, path::AbstractString;
        axis_width_cm=8.0, axis_height_cm=6.0,
        layout=(1, length(tpds)),
        ylims=nothing,
        kwargs...)

    nrows, ncols = layout
    fig_w, fig_h, positions = _layout_axes(
        axis_width_cm, axis_height_cm, length(tpds); nrows, ncols)
    fig = PythonPlot.figure(figsize=(fig_w, fig_h))
    for (i, tpd) in enumerate(tpds)
        i > length(positions) && break
        ax = fig.add_axes(collect(positions[i]))
        _plot_transport_on_ax!(ax, tpd; kwargs...)
        isnothing(ylims) || ax.set_ylim(ylims...)
    end

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

end # module
