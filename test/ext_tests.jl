# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - extension tests

using BoltzTraP
using Test

@testset "BoltzTraP.jl (ext)" begin
    include("test_publication.jl")
    if get(ENV, "TEST_WANNIER", "false") == "true"
        include("test_wannier_ext.jl")
        si_wannier_reftest_data =
            joinpath(@__DIR__, "..", "reftest", "data", "si_wannier_transport.jld2")
        if isfile(si_wannier_reftest_data)
            include("test_wannier_reftest.jl")
        else
            @warn "Si Wannier reftest golden not found — skipping regression" path =
                si_wannier_reftest_data
        end
    end
end
