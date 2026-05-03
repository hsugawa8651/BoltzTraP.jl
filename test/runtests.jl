# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests

using Test

const TEST_GROUP = isempty(ARGS) ? "all" : lowercase(ARGS[1])
@assert TEST_GROUP in ("all", "core", "ext") "TEST_GROUP must be all/core/ext, got: $TEST_GROUP"

if TEST_GROUP in ("all", "core")
    include("core_tests.jl")
end
if TEST_GROUP in ("all", "ext")
    include("ext_tests.jl")
end
