# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
Wannier.jl extension for BoltzTraP.jl.

Provide a `WannierInterpolator` constructor that consumes a Wannier90
prefix (`.win`, `.chk`, etc.) via `Wannier.read_w90_interp` and stores
the resulting `Wannier.InterpModel` for band/velocity interpolation.
"""
module BoltzTraPWannierExt

using BoltzTraP
using BoltzTraP: WannierInterpolator, AbstractInterpolator
using BoltzTraP: SpinType, Unpolarized
using BoltzTraP: BOHR_TO_ANG, HA_TO_EV
using StaticArrays
using Wannier

"""
    WannierInterpolator(prefix::String; spintype=Unpolarized()) -> WannierInterpolator

Construct a `WannierInterpolator` from a Wannier90 file prefix by reading
`<prefix>.win` plus the unitary matrices from `<prefix>.chk` (or
`<prefix>.chk.fmt`) via `Wannier.read_w90_interp`, then delegating to
`WannierInterpolator(::Wannier.InterpModel)`.
"""
function BoltzTraP.WannierInterpolator(prefix::String; spintype = Unpolarized())
    return BoltzTraP.WannierInterpolator(Wannier.read_w90_interp(prefix); spintype)
end

"""
    WannierInterpolator(imodel::Wannier.InterpModel; spintype=Unpolarized()) -> WannierInterpolator

Construct a `WannierInterpolator` from an existing `Wannier.InterpModel`.

This entry point lets callers feed an InterpModel produced by any
in-process Wannierization workflow (`Wannier.read_w90` followed by
`disentangle` and `Wannier.InterpModel`, custom initial projections,
etc.) without first writing a `.chk` binary to disk.

The lattice is converted from the Wannier.jl Ångström convention to
Bohr to match the project's atomic-units policy (Hartree for energy,
Bohr for length).
"""
function BoltzTraP.WannierInterpolator(
    imodel::Wannier.InterpModel;
    spintype = Unpolarized(),
)
    lattice_ang = SMatrix{3,3,Float64}(imodel.lattice)
    lattice_bohr = lattice_ang ./ BOHR_TO_ANG
    return WannierInterpolator{typeof(spintype)}(imodel, lattice_bohr, spintype)
end

"""
    interpolate_bands(interp::WannierInterpolator, kpoints) -> ebands

Interpolate band energies at fractional `kpoints` (`nk × 3`, BoltzTraP
convention) using the stored `Wannier.InterpModel`. Wannier.jl expects
`(3, nk)` and returns energies in eV; this method transposes the input
and converts the output from eV to Hartree to match the BoltzTraP
atomic-units policy (Hartree for energy, Bohr for length).

Return shape: `(n_wann, nk)` in Hartree.
"""
function BoltzTraP.interpolate_bands(interp::WannierInterpolator, kpoints)
    kpts_3xN = Matrix{Float64}(transpose(kpoints))
    E_eV = Wannier.interpolate(interp.model, kpts_3xN)
    return E_eV ./ HA_TO_EV
end

"""
    interpolate_velocities(interp::WannierInterpolator, kpoints) -> vbands

Interpolate group velocities at fractional `kpoints` (`nk × 3`, BoltzTraP
convention) using the stored `Wannier.InterpModel`. `Wannier.velocity`
returns shape `(n_wann, nk, 3)` in eV·Å; this method permutes to the
BoltzTraP shape `(3, n_wann, nk)` and converts to Hartree·Bohr to match
the project's atomic-units policy.

Wannier.jl's `velocity` requires the MDRS R-vectors, which are produced
by default by `Wannier.read_w90_interp` (the `WannierInterpolator(::String)`
constructor); a non-MDRS model would fail here.

Return shape: `(3, n_wann, nk)` in Hartree·Bohr.
"""
function BoltzTraP.interpolate_velocities(interp::WannierInterpolator, kpoints)
    kpts_3xN = Matrix{Float64}(transpose(kpoints))
    Rvecs = interp.model.kRvectors.Rvectors
    H = interp.model.H
    _E, v_eVA = Wannier.velocity(Rvecs, H, kpts_3xN)
    v_BT = permutedims(v_eVA, (3, 1, 2))
    return v_BT ./ (HA_TO_EV * BOHR_TO_ANG)
end

"""
    velocity_matrices(interp::WannierInterpolator, kpoints) -> A

Interpolate the full band-gauge covariant velocity matrices
``\\bar{H}_\\alpha = U^\\dagger (\\partial_\\alpha H^W) U`` at fractional `kpoints`
(`nk × 3`, BoltzTraP convention). The real diagonal is the band group velocity
(equal to [`interpolate_velocities`](@ref)); the off-diagonal entries within a
degenerate multiplet are the intra-manifold velocity matrix elements required for a
gauge-invariant transport tensor at degenerate k-points.

The covariant part is read from `Wannier._get_D` (its `Haᴴ` return); the Berry
connection matrix `D` that `_get_D` also forms is discarded — `Haᴴ` itself contains
no `1/ΔE` and is finite at degeneracies. Like `interpolate_velocities`, this
requires the MDRS R-vectors produced by default by `Wannier.read_w90_interp` (the
`WannierInterpolator(::String)` constructor); a non-MDRS model would fail in
`_get_D`. `_get_D` is a Wannier internal — its diagonal is regression-guarded
against `interpolate_velocities` (see `test/test_wannier_ext.jl`).

Return shape: `(n_wann, n_wann, 3, nk)` in Hartree·Bohr.
"""
function BoltzTraP.velocity_matrices(interp::WannierInterpolator, kpoints)
    kpts_3xN = Matrix{Float64}(transpose(kpoints))
    Rvecs = interp.model.kRvectors.Rvectors
    H = interp.model.H
    _E, _U, HaH, _D = Wannier._get_D(Rvecs, H, kpts_3xN; use_degen_pert = false)
    # HaH: (n_wann, n_wann, nk, 3) in eV·Å, band gauge.
    A = permutedims(HaH, (1, 2, 4, 3)) ./ (HA_TO_EV * BOHR_TO_ANG)
    return A
end

end # module BoltzTraPWannierExt
