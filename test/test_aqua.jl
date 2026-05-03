# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Aqua

@testset "Aqua quality assurance" begin
    Aqua.test_all(BoltzTraP)
end
