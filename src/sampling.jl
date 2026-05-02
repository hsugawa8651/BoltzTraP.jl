# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

#=
    AbstractBZSampling

Abstract type for Brillouin-zone sampling strategies used by transport
calculations. Concrete subtypes (e.g. `UniformMesh`) carry the parameters
that determine how `solve_transport` integrates over the BZ.

This is the dispatch axis for BZ sampling strategies; alternative methods
(adaptive grids, AutoBZ-style integration, wavelet meshes) can be added
later as additional concrete subtypes.
=#
abstract type AbstractBZSampling end

#=
    UniformMesh(nk=nothing)
    UniformMesh((nk1, nk2, nk3))

Uniform Monkhorst-Pack-style mesh. The `nk` field selects the mesh density:

  - `nk === nothing`: the interpolator decides the mesh size (matches the
    legacy `run_integrate` behavior, which derives the mesh from the
    interpolator's k-point density).
  - `nk::Tuple{Int,Int,Int}`: explicit user-specified mesh dimensions.
=#
struct UniformMesh <: AbstractBZSampling
    nk::Union{Nothing,Tuple{Int,Int,Int}}
end

UniformMesh() = UniformMesh(nothing)
