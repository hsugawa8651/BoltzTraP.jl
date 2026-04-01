# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Tests for plotting data types

@testset "BandPlotData" begin
    bpd = BandPlotData(
        [0.0, 0.1, 0.2],              # kdist
        [1.0 2.0; 1.1 2.1; 1.2 2.2],  # ebands_eV (3 kpts × 2 bands)
        ["Γ", "X"],                     # labels
        [0.0, 0.2],                     # positions
        -5.0, 5.0,                      # emin, emax
        0.1,                            # fermi_Ha
        "Test"                          # title
    )
    @test bpd isa BandPlotData
    @test size(bpd.ebands_eV) == (3, 2)
    @test bpd.unique_labels == ["Γ", "X"]
    @test bpd.emin == -5.0
    @test bpd.title == "Test"
end

@testset "TransportPlotData" begin
    tpd = TransportPlotData(
        [0.0, 0.1, 0.2],     # xdata
        [10.0, 20.0, 30.0],  # ydata
        "μ - E_F (eV)",      # xlabel
        "S (μV/K)",           # ylabel
        "Seebeck xx"          # title
    )
    @test tpd isa TransportPlotData
    @test length(tpd.xdata) == 3
    @test tpd.xlabel == "μ - E_F (eV)"
end
