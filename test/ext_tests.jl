# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - extension tests

using BoltzTraP
using Test

@testset "BoltzTraP.jl (ext)" begin
    include("test_publication.jl")
end
