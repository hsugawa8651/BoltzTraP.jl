# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# CoSb3 conventional cubic cell — single-source geometry for reftest
# baselines.
#
# Used as the lattice / position source by:
#
#   - `docs/src/wannier.md` Validation section (Wannier path documented
#     golden), and
#   - the Fourier path reference baseline.
#
# Keeping a single copy of the geometry in this file prevents drift
# between the two regression baselines.
#
# Source: BoltzTraP2-public/data/CoSb3/CoSb3.struct (Wien2k .struct
# format, Bohr by default). Space group 204 (Im-3, body-centered cubic
# skutterudite). The conventional 16-atom cell is used directly (4 Co +
# 12 Sb) to match the BT2 Wien2k reference data layout.
#
# A development copy with DFTK pseudopotential loading is kept at
# `boltztrap-jl-dev/wannier-smoke/cosb3-dftk/geometry.jl`; that copy
# `include`s this file's positions / lattice manually today and must be
# updated in lock-step if these constants change.

using StaticArrays
using LinearAlgebra: I

# Cubic conventional cell, `a` is in Bohr (Wien2k .struct default unit).
# Literature value (Pizzi 2014 reference context) is a ≈ 9.04 Å.
const a_bohr = 17.080296                  # Bohr
const a_ang  = a_bohr * 0.529177210903    # ≈ 9.0385 Å

# Lattice sanity check — catches Bohr/Å swap errors before launching
# multi-hour DFT runs. The repository carries a recorded incident
# (2026-05-15) where a previous revision mislabeled the unit and gave
# sigma about forty times below the Pizzi 2014 reference; this assert
# blocks that mode of failure at include time.
@assert 0.5 < a_ang < 30.0 "Lattice a = $a_ang Å is outside physical range — unit error?"
@assert abs(a_ang - 9.04) / 9.04 < 0.10 "Lattice a = $a_ang Å diverges >10% from CoSb3 literature value 9.04 Å — unit error?"

# Lattice vectors (column-wise, Bohr) for the conventional cubic cell.
const cosb3_lattice = a_bohr * SMatrix{3,3,Float64}(I)

# 4 Co atoms in conventional cell (Wyckoff 8c-like positions).
const cosb3_co_positions = [
    [0.25, 0.25, 0.25],
    [0.75, 0.75, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
]

# 12 Sb atoms in conventional cell (Wyckoff 24g-like positions).
const cosb3_sb_positions = [
    [0.00000, 0.33537, 0.15788],
    [0.00000, 0.66463, 0.84212],
    [0.00000, 0.66463, 0.15788],
    [0.00000, 0.33537, 0.84212],
    [0.15788, 0.00000, 0.33537],
    [0.84212, 0.00000, 0.66463],
    [0.84212, 0.00000, 0.33537],
    [0.15788, 0.00000, 0.66463],
    [0.33537, 0.15788, 0.00000],
    [0.66463, 0.84212, 0.00000],
    [0.33537, 0.84212, 0.00000],
    [0.66463, 0.15788, 0.00000],
]

# Combined fractional positions and matching species symbols, in the
# order baseline reference data expects (4 Co followed by 12 Sb).
const cosb3_positions = vcat(cosb3_co_positions, cosb3_sb_positions)
const cosb3_species   = vcat(fill(:Co, length(cosb3_co_positions)),
                             fill(:Sb, length(cosb3_sb_positions)))
