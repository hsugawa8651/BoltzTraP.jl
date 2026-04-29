# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Tests for plotting data types

@testset "BandPlotData" begin
    # Convention: ebands_eV is (nbands × nk); see src/plottypes.jl docstring.
    bpd = BandPlotData(
        [0.0, 0.1, 0.2],                # kdist (3 kpts)
        [1.0 1.1 1.2; 2.0 2.1 2.2],     # ebands_eV (2 bands × 3 kpts)
        ["Γ", "X"],                     # labels
        [0.0, 0.2],                     # positions
        -5.0, 5.0,                      # emin, emax
        0.1,                            # fermi_Ha
        "Test",                         # title
    )
    @test bpd isa BandPlotData
    @test size(bpd.ebands_eV) == (2, 3)
    @test size(bpd.ebands_eV, 2) == length(bpd.kdist)  # nk consistency
    @test bpd.unique_labels == ["Γ", "X"]
    @test bpd.emin == -5.0
    @test bpd.title == "Test"
end

@testset "TransportPlotData" begin
    tpd = TransportPlotData(
        [0.0, 0.1, 0.2],     # xdata
        [10.0, 20.0, 30.0],  # ydata
        "μ - E_F (eV)",      # xlabel
        "S (μV/K)",          # ylabel
        "Seebeck xx",        # title
    )
    @test tpd isa TransportPlotData
    @test length(tpd.xdata) == 3
    @test tpd.xlabel == "μ - E_F (eV)"
end

# ── Phase B' regression: build_transport_plot_data builder ──
#
# Lightweight TransportResult fixture used by both the builder unit tests
# and the plot_transport (Plots ext) regression tests further below.
const _PT_NT, _PT_NMU = 2, 5
let
    sigma   = zeros(3, 3, _PT_NT, _PT_NMU)
    seebeck = zeros(3, 3, _PT_NT, _PT_NMU)
    kappa   = zeros(3, 3, _PT_NT, _PT_NMU)
    for d in 1:3
        sigma[d, d, :, :] .= 1.0
        kappa[d, d, :, :] .= 0.1
    end
    # Seebeck values in V/K — builder must scale by 1e6 to μV/K.
    seebeck[1, 1, 1, :] .= [1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6]
    seebeck[1, 1, 2, :] .= [6.0e-6, 7.0e-6, 8.0e-6, 9.0e-6, 1.0e-5]
    md = Dict{String,Any}("fermi_Ha" => 0.0)
    global _PT_RESULT = TransportResult(
        [100.0, 300.0],
        collect(LinRange(-0.5, 0.5, _PT_NMU)),
        sigma,
        seebeck,
        kappa,
        nothing,
        md,
    )
end

@testset "build_transport_plot_data: happy paths" begin
    # 4 quantity × 2 abscissa × 2 component coverage matrix.
    quantities = ("seebeck", "sigma", "kappa", "S")
    abscissae = ("mu", "T")
    components = ("xx", "zz")
    for q in quantities, a in abscissae, c in components
        tpd = build_transport_plot_data(
            _PT_RESULT;
            quantity = q,
            component = c,
            abscissa = a,
            temperature = 300.0,
        )
        @test tpd isa TransportPlotData
        @test length(tpd.xdata) == length(tpd.ydata)
    end
end

@testset "build_transport_plot_data: labels use display symbols" begin
    tpd = build_transport_plot_data(
        _PT_RESULT;
        quantity = "seebeck",
        component = "xx",
        abscissa = "mu",
        temperature = 300.0,
    )
    @test tpd.xlabel == "μ - E_F (eV)"
    @test tpd.ylabel == "S_xx (μV/K)"
    @test tpd.title == "T = 300 K"

    @test build_transport_plot_data(_PT_RESULT;
        quantity = "S", component = "xx", abscissa = "mu").ylabel == "S_xx (μV/K)"
    @test build_transport_plot_data(_PT_RESULT;
        quantity = "sigma", component = "xx").ylabel == "σ_xx (S/m)"
    @test build_transport_plot_data(_PT_RESULT;
        quantity = "kappa", component = "yy").ylabel == "κ_yy (W/m/K)"

    tpd_T = build_transport_plot_data(
        _PT_RESULT;
        quantity = "sigma",
        component = "zz",
        abscissa = "T",
        mu_index = 3,
    )
    @test tpd_T.xlabel == "Temperature (K)"
    @test tpd_T.xdata == [100.0, 300.0]
end

@testset "build_transport_plot_data: rejects unknown inputs" begin
    @test_throws ArgumentError build_transport_plot_data(_PT_RESULT; abscissa = "bogus")
    @test_throws ArgumentError build_transport_plot_data(_PT_RESULT; quantity = "banana")
    @test_throws ArgumentError build_transport_plot_data(_PT_RESULT; component = "qq")
    @test_throws ArgumentError build_transport_plot_data(_PT_RESULT; temperature = 999.0)
end

@testset "build_transport_plot_data: Seebeck V/K → μV/K conversion" begin
    tpd = build_transport_plot_data(
        _PT_RESULT;
        quantity = "seebeck",
        component = "xx",
        abscissa = "mu",
        temperature = 300.0,
    )
    @test tpd.ydata ≈ [6.0, 7.0, 8.0, 9.0, 10.0]
end

# ── Phase B' regression: plot_transport thin wrapper (Plots ext) ──

# BoltzTraPPlotsExt is a multi-trigger extension (["Plots", "Brillouin"]),
# so both packages must be loaded for the ext methods to materialize.
using Plots
using Brillouin

@testset "plot_transport thin wrapper preserves API" begin
    p = plot_transport(
        _PT_RESULT;
        quantity = "seebeck",
        component = "xx",
        abscissa = "mu",
        temperature = 300.0,
    )
    @test p isa Plots.Plot

    # Unknown abscissa is rejected at the public API level.
    @test_throws ArgumentError plot_transport(_PT_RESULT; abscissa = "bogus")
end

@testset "plot_transport ylabel uses display symbols" begin
    # The ylabel is set by the recipe from tpd.ylabel; ensure the
    # display symbols ("S_xx", "σ_xx", "κ_xx") propagate end-to-end.
    p_S = plot_transport(_PT_RESULT; quantity = "seebeck", component = "xx")
    p_sig = plot_transport(_PT_RESULT; quantity = "sigma", component = "xx")
    p_kap = plot_transport(_PT_RESULT; quantity = "kappa", component = "xx")
    @test occursin("S_xx", p_S.subplots[1].attr[:yaxis][:guide])
    @test occursin("σ_xx", p_sig.subplots[1].attr[:yaxis][:guide])
    @test occursin("κ_xx", p_kap.subplots[1].attr[:yaxis][:guide])
end
