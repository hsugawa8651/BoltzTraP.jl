# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Example: publication-quality transport figures
#
# Render publication-quality PDFs for Si transport coefficients using
# `build_transport_plot_data` + `savefig_publication`. Reuses the
# transport result already produced by 110_si_vasp_pbe_paw.jl.
#
# DFT conditions: see 110_si_vasp_pbe_paw.jl.
#
# Usage:
#   julia --project=. examples/110_si_vasp_pbe_paw.jl       # produce *_transport.jld2
#   julia --project=. examples/119_si_vasp_pbe_paw_publication.jl
#
# Output:
#   examples/119_publication_Si_VASP_seebeck_300K.pdf
#   examples/119_publication_Si_VASP_sigma_300K.pdf
#   examples/119_publication_Si_VASP_kappa_300K.pdf
#   examples/119_publication_Si_VASP_combined_300K.pdf

using BoltzTraP
using PythonPlot

transport_jld2 = joinpath(@__DIR__, "110_si_vasp_pbe_paw_transport.jld2")
transport = load_integrate(transport_jld2)

stem = "119_publication_Si_VASP_"

# Single-quantity figures at T = 300 K
for (q, tag) in (("seebeck", "seebeck"), ("sigma", "sigma"), ("kappa", "kappa"))
    tpd = build_transport_plot_data(
        transport;
        quantity = q,
        component = "xx",
        abscissa = "mu",
        temperature = 300.0,
    )
    out = joinpath(@__DIR__, stem * tag * "_300K.pdf")
    savefig_publication(tpd, out; axis_width_cm = 8.0, axis_height_cm = 5.0)
    println("Saved: $out")
end

# Combined 1×3 subplot showing S, σ, κ vs μ at T = 300 K
tpds = [
    build_transport_plot_data(transport;
        quantity = q, component = "xx", abscissa = "mu", temperature = 300.0)
    for q in ("seebeck", "sigma", "kappa")
]
out_combined = joinpath(@__DIR__, stem * "combined_300K.pdf")
savefig_publication(
    tpds, out_combined;
    axis_width_cm = 6.0, axis_height_cm = 4.5,
    layout = (1, 3),
)
println("Saved: $out_combined")
