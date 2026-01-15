#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl reftest
"""
Generate reference data for calc_cv (electronic heat capacity).

Uses BoltzTraP2 bundled Si data to generate reference values.

Usage:
    cd BoltzTraP2-public
    python ../BoltzTraP.jl/reftest/generate_cv_reference.py
"""
import numpy as np
import os
import sys

from common import get_material_path, OUTPUT_DIR

import BoltzTraP2.dft
import BoltzTraP2.sphere
import BoltzTraP2.fite
import BoltzTraP2.bandlib


def main():
    """Generate calc_cv reference data."""
    # Load Si data (bundled with BoltzTraP2)
    DATA_DIR = str(get_material_path("Si.vasp"))
    print(f"Loading data from: {DATA_DIR}")
    data = BoltzTraP2.dft.DFTData(DATA_DIR)

    print(f"  ebands shape: {data.ebands.shape}")
    print(f"  fermi: {data.fermi}")
    print(f"  nelect: {data.nelect}")
    print(f"  dosweight: {data.dosweight}")

    # Get equivalences for interpolation
    equivalences = BoltzTraP2.sphere.get_equivalences(data.atoms, data.magmom, len(data.kpoints))
    print(f"  n_equivalences: {len(equivalences)}")

    # Interpolation
    coeffs = BoltzTraP2.fite.fitde3D(data, equivalences)
    print(f"  coeffs shape: {coeffs.shape}")

    # Reconstruct bands on dense grid
    eband, vvband, _ = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, data.get_lattvec())
    print(f"  eband shape: {eband.shape}")

    # Calculate DOS
    epsilon, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband, npts=10000)
    print(f"  epsilon shape: {epsilon.shape}")
    print(f"  dos shape: {dos.shape}")

    # Define mu and T ranges for calc_cv
    mu_range = np.linspace(data.fermi - 0.1, data.fermi + 0.1, 21)  # 21 points around Fermi
    T_range = np.array([100.0, 200.0, 300.0, 400.0, 500.0])  # 5 temperatures

    print(f"  mu_range: {mu_range.shape}, min={mu_range.min():.6f}, max={mu_range.max():.6f}")
    print(f"  T_range: {T_range}")

    # Calculate heat capacity using BoltzTraP2
    cv = BoltzTraP2.bandlib.calc_cv(epsilon, dos, mu_range, T_range, dosweight=data.dosweight)
    print(f"  cv shape: {cv.shape}")  # (nT, nmu)
    print(f"  cv min: {cv.min():.6e}, max: {cv.max():.6e}")

    # Save reference data
    OUTPUT = os.path.join(OUTPUT_DIR, "cv_reference.npz")
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

    np.savez(OUTPUT,
        epsilon=epsilon,
        dos=dos,
        mu_range=mu_range,
        T_range=T_range,
        dosweight=data.dosweight,
        fermi=data.fermi,
        cv=cv,
    )
    print(f"Saved: {OUTPUT}")
    print("Done!")


if __name__ == '__main__':
    main()
