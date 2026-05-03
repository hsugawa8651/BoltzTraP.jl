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
using BoltzTraP: BOHR_TO_ANG
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

end # module BoltzTraPWannierExt
