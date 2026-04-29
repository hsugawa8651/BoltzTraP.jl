# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - PythonPlot extension tests

using Test
using BoltzTraP
using BoltzTraP: BandPlotData, TransportPlotData
using PythonPlot

@testset "savefig_publication (PythonPlot)" begin
    bpd = BandPlotData(
        [0.0, 1.0, 2.0, 3.0],                    # kdist
        Float64[-2.0 -1.5 -1.0 -1.5;             # ebands_eV (3 bands × 4 kpts)
                 0.0  0.5  1.0  0.5;
                 2.0  2.5  3.0  2.5],
        ["Γ", "X", "M", "Γ"],                    # unique_labels
        [0.0, 1.0, 2.0, 3.0],                    # unique_positions
        -3.0, 4.0,                               # emin, emax
        0.0,                                     # fermi_Ha
        "Si bands (test fixture)",               # title
    )

    tpd = TransportPlotData(
        collect(-1.0:0.1:1.0),                   # xdata (μ-Ef in eV)
        sin.(collect(-1.0:0.1:1.0)) .* 100.0,    # ydata
        "μ - E_F (eV)",
        "S_xx (μV/K)",
        "Si Seebeck (test fixture)",
    )

    @testset "BandPlotData single" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands.pdf")
            ret = savefig_publication(bpd, path)
            @test isfile(path)
            @test ret == path
            @test filesize(path) > 0
        end
    end

    @testset "axis_width_cm / axis_height_cm custom" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands_custom.pdf")
            savefig_publication(bpd, path; axis_width_cm = 10.0, axis_height_cm = 5.0)
            @test isfile(path)
        end
    end

    @testset "title / ylims kwargs" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands_titled.pdf")
            savefig_publication(bpd, path; title = "Test", ylims = (-2.0, 3.0))
            @test isfile(path)
        end
    end

    @testset "subplot layout=(1,2)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands_1x2.pdf")
            savefig_publication([bpd, bpd], path; layout = (1, 2))
            @test isfile(path)
        end
    end

    @testset "subplot layout=(2,1)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands_2x1.pdf")
            savefig_publication([bpd, bpd], path; layout = (2, 1))
            @test isfile(path)
        end
    end

    @testset "overlay mode" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands_overlay.pdf")
            savefig_publication([bpd, bpd], path; overlay = true)
            @test isfile(path)
        end
    end

    @testset "PNG output" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bands.png")
            savefig_publication(bpd, path)
            @test isfile(path)
        end
    end

    @testset "TransportPlotData single" begin
        mktempdir() do tmp
            path = joinpath(tmp, "transport.pdf")
            ret = savefig_publication(tpd, path)
            @test isfile(path)
            @test ret == path
        end
    end

    @testset "TransportPlotData subplot=(1,2)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "transport_1x2.pdf")
            savefig_publication([tpd, tpd], path; layout = (1, 2))
            @test isfile(path)
        end
    end
end
