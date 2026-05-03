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

Construct a `WannierInterpolator` from a Wannier90 file prefix.

`Wannier.read_w90_interp` reads `<prefix>.win` plus the unitary matrices
from `<prefix>.chk` (or `<prefix>.chk.fmt`) and returns an `InterpModel`
ready for k-space interpolation. The lattice is converted from the
Wannier.jl Ångström convention to Bohr so downstream BoltzTraP transport
code stays in atomic units per the project's units policy.
"""
function BoltzTraP.WannierInterpolator(prefix::String; spintype = Unpolarized())
    imodel = Wannier.read_w90_interp(prefix)
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
    v_eVA = Wannier.velocity(Rvecs, H, kpts_3xN)
    v_BT = permutedims(v_eVA, (3, 1, 2))
    return v_BT ./ (HA_TO_EV * BOHR_TO_ANG)
end

end # module BoltzTraPWannierExt
