# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Generate Fe collinear transport plot at T=1000K.
# Compares Julia BoltzTraP.jl with Python BoltzTraP2.
#
# Usage:
#   julia --project validation/plot_fe_transport.jl

using NPZ
using BoltzTraP
using Plots
using LaTeXStrings
using Printf
using Statistics

gr()

const HA_TO_EV = 27.211386245988
const ANG_TO_BOHR = 1.8897261246257702

function generate_fe_transport_plot()
    println("=" ^ 60)
    println("Fe Collinear Transport @ T=1000K")
    println("=" ^ 60)

    # Load reference data
    ref = npzread("reftest/data/fe_collinear_e2e.npz")

    # Data dimensions
    ebands_concat = ref["ebands"]  # (24, 47)
    nbands_per_spin = size(ebands_concat, 1) ÷ 2  # 12
    nkpts = size(ebands_concat, 2)  # 47

    ebands_up = ebands_concat[1:nbands_per_spin, :]
    ebands_down = ebands_concat[nbands_per_spin+1:end, :]
    ebands_3d = zeros(nbands_per_spin, nkpts, 2)
    ebands_3d[:, :, 1] = ebands_up
    ebands_3d[:, :, 2] = ebands_down

    println("Bands per spin: $nbands_per_spin")
    println("K-points: $nkpts")

    # Unit conversion
    lattvec = ref["lattvec"] * ANG_TO_BOHR
    fermi = ref["fermi"]
    nelect = ref["nelect"]
    mur = ref["mur"]
    temperatures = [1000.0]

    # === Collinear (Total) calculation ===
    println("\n--- Collinear (Total) ---")
    data = BoltzTraP.DFTData(
        lattice = lattvec,
        positions = zeros(3, 1),
        species = ["Fe"],
        kpoints = ref["kpoints"]',
        weights = ones(nkpts) / nkpts,
        ebands = ebands_3d,
        occupations = zeros(nbands_per_spin, nkpts, 2),
        fermi = fermi,
        nelect = nelect,
        magmom = ref["magmom"],
    )
    interp = BoltzTraP.run_interpolate(data; kpoints=nkpts, verbose=false)
    transport = BoltzTraP.run_integrate(interp, mur; temperatures=temperatures, verbose=false)

    # Mu relative to Fermi (in eV)
    fermi_eV = fermi * HA_TO_EV
    mu_rel = transport.mu_values .- fermi_eV

    # Extract transport data
    sigma_jl = vec(transport.sigma[1, 1, 1, :])
    S_jl = vec(transport.seebeck[1, 1, 1, :]) .* 1e6
    kappa_jl = vec(transport.kappa[1, 1, 1, :])

    # Python reference
    iT_1000K = findfirst(t -> t ≈ 1000.0, ref["Tr"])
    sigma_py = ref["sigma"][iT_1000K, :, 1, 1]
    S_py = ref["S"][iT_1000K, :, 1, 1] .* 1e6
    kappa_py = ref["kappa"][iT_1000K, :, 1, 1]

    # Common plot settings
    xlims = (-2.5, 2.5)
    idx = 1:5:length(mu_rel)  # Subsample for markers

    # Compute y-axis limits with integer powers
    function int_pow_limits(data; margin=0.5)
        pos_data = data[data .> 0]
        if isempty(pos_data)
            return (1e15, 1e25)
        end
        min_pow = floor(Int, log10(minimum(pos_data)))
        max_pow = ceil(Int, log10(maximum(pos_data)))
        return (10.0^min_pow, 10.0^max_pow)
    end

    sigma_ylims = int_pow_limits(vcat(sigma_jl, sigma_py))
    kappa_ylims = int_pow_limits(vcat(kappa_jl, kappa_py))

    println("\n--- Generating Plot ---")

    # === Total Plot (Python vs Julia comparison) ===
    p1 = plot(
        mu_rel[idx], S_py[idx],
        seriestype=:scatter, marker=:circle, markersize=4, color=:red,
        label="Python BoltzTraP2",
        ylabel=L"$S_{xx}$ ($\mu$V/K)",
        title="Fe Collinear @ 1000 K",
        legend=:topright, legendfontsize=7,
        xformatter=_->"",
        xlims=xlims,
    )
    plot!(p1, mu_rel[idx], S_jl[idx],
          seriestype=:scatter, marker=:xcross, markersize=4, color=:blue,
          label="Julia BoltzTraP.jl")

    p2 = plot(
        mu_rel[idx], sigma_py[idx],
        seriestype=:scatter, marker=:circle, markersize=4, color=:red,
        label="",
        ylabel=L"$\sigma_{xx}/\tau$ (1/$\Omega$ms)",
        yscale=:log10, ylims=sigma_ylims,
        xformatter=_->"",
        xlims=xlims,
    )
    plot!(p2, mu_rel[idx], sigma_jl[idx],
          seriestype=:scatter, marker=:xcross, markersize=4, color=:blue,
          label="")

    p3 = plot(
        mu_rel[idx], kappa_py[idx],
        seriestype=:scatter, marker=:circle, markersize=4, color=:red,
        label="",
        ylabel=L"$\kappa_{xx}/\tau$ (W/mKs)",
        xlabel=L"$\mu - E_F$ (eV)",
        yscale=:log10, ylims=kappa_ylims,
        xlims=xlims,
    )
    plot!(p3, mu_rel[idx], kappa_jl[idx],
          seriestype=:scatter, marker=:xcross, markersize=4, color=:blue,
          label="")

    fig = plot(p1, p2, p3, layout=(3, 1), size=(400, 500),
               left_margin=3Plots.mm, right_margin=2Plots.mm,
               bottom_margin=2Plots.mm, top_margin=2Plots.mm)

    output_path = joinpath(@__DIR__, "transport_Fe_1000K_total.png")
    savefig(fig, output_path)
    println("Saved: $output_path")

    # Agreement statistics
    println("\n--- Agreement Statistics (T=1000K) ---")

    function compute_stats(jl, py, name)
        valid = abs.(py) .> 1e-10
        if !any(valid)
            println("  $name: No valid points")
            return
        end
        rel_err = abs.(jl[valid] .- py[valid]) ./ abs.(py[valid])
        max_err = maximum(rel_err) * 100
        mean_err = mean(rel_err) * 100
        correlation = cor(jl, py)
        @printf("  %s: max rel error = %.4f%%, mean = %.4f%%, correlation = %.6f\n",
                name, max_err, mean_err, correlation)
    end

    compute_stats(sigma_jl, sigma_py, "σ_xx")
    compute_stats(S_jl, S_py, "S_xx")
    compute_stats(kappa_jl, kappa_py, "κ_xx")

    println("\n✓ Julia BoltzTraP.jl matches Python BoltzTraP2")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_fe_transport_plot()
end
