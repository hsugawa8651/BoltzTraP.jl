# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
    AbstractBZSampling

Abstract type for Brillouin-zone sampling strategies used by transport
calculations. Concrete subtypes (e.g. `UniformMesh`) carry the parameters
that determine how `solve_transport` integrates over the BZ.

This is the dispatch axis for BZ sampling strategies; alternative methods
(adaptive grids, AutoBZ-style integration, wavelet meshes) can be added
later as additional concrete subtypes.
"""
abstract type AbstractBZSampling end

"""
    UniformMesh(nk=nothing)
    UniformMesh((nk1, nk2, nk3))

Uniform Monkhorst-Pack-style mesh. The `nk` field selects the mesh density:

  - `nk === nothing`: the interpolator decides the mesh size (matches the
    legacy `run_integrate` behavior, which derives the mesh from the
    interpolator's k-point density).
  - `nk::Tuple{Int,Int,Int}`: explicit user-specified mesh dimensions.
"""
struct UniformMesh <: AbstractBZSampling
    nk::Union{Nothing,Tuple{Int,Int,Int}}
end

UniformMesh() = UniformMesh(nothing)

"""
    TransportSystem(lattice, fermi, nelect, dosweight, spintype)
    TransportSystem(interp::InterpolationResult)

System-level constants needed for transport calculations: lattice, Fermi
energy, electron count, DOS weight, and spin type. These quantities are
conventionally spread across `InterpolationResult.metadata` and the
underlying band data; this struct unifies them as a single, dispatched
data carrier for the new `solve_transport` API.

Fields:
- `lattice::SMatrix{3,3,Float64,9}`: 3×3 lattice vectors (columns) in Bohr
- `fermi::Float64`: Fermi energy in Hartree
- `nelect::Float64`: number of electrons per unit cell
- `dosweight::Float64`: DOS weight (2.0 = non-magnetic, 1.0 = spin-polarized)
- `spintype::SpinType`: `Unpolarized()` / `Collinear()` / `NonCollinear()`

The `TransportSystem(::InterpolationResult)` constructor is defined in
`io/serialization.jl` (after `InterpolationResult` itself is defined).
"""
struct TransportSystem
    lattice::SMatrix{3,3,Float64,9}
    fermi::Float64
    nelect::Float64
    dosweight::Float64
    spintype::SpinType
end
