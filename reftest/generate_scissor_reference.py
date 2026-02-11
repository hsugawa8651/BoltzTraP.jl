#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Generate reference data for apply_scissor tests

"""
Generate reference data for apply_scissor function.

Usage:
    python generate_scissor_reference.py

Output:
    reftest/data/scissor_reference.npz
"""

import sys
from pathlib import Path

import numpy as np

# Add BoltzTraP2 to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "BoltzTraP2-public"))

import BoltzTraP2.dft as dft
import BoltzTraP2.fite as fite
import BoltzTraP2.sphere as sphere
import BoltzTraP2.bandlib as bandlib


def main():
    # Constants
    eV = 1 / 27.211386245988  # Ha per eV

    # Output path
    output_dir = Path(__file__).parent / "data"
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "scissor_reference.npz"

    print("=== Generating scissor reference data ===")

    # 1.2.1 Load Si data
    data_path = Path(__file__).parent.parent.parent / "BoltzTraP2-public" / "data" / "Si.vasp"
    data = dft.DFTData(str(data_path))
    print(f"Loaded: {data_path}")
    print(f"  nelect: {data.nelect}")
    print(f"  dosweight: {data.dosweight}")

    # 1.2.2 Interpolation
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, 5000)
    coeffs = fite.fitde3D(data, equivalences)
    eband_before, vvband, _ = fite.getBTPbands(equivalences, coeffs, data.lattvec)
    print(f"  eband shape: {eband_before.shape}")

    # 1.2.3 Calculate DOS
    epsilon, dos, vvdos, _ = bandlib.BTPDOS(eband_before, vvband)
    print(f"  epsilon shape: {epsilon.shape}")
    print(f"  dos shape: {dos.shape}")

    # Check current gap
    fermi = bandlib.solve_for_mu(
        epsilon, dos, data.nelect, 0.0, data.dosweight, refine=False, try_center=False
    )
    pos = np.abs(epsilon - fermi).argmin()
    assert dos[pos] == 0, "No gap at Fermi level"

    # Find current gap
    for i in range(pos, -1, -1):
        if dos[i] != 0:
            vbm = epsilon[i]
            break
    for i in range(pos, len(dos)):
        if dos[i] != 0:
            cbm = epsilon[i]
            break
    current_gap = cbm - vbm
    print(f"  Current gap: {current_gap:.6f} Ha = {current_gap / eV:.4f} eV")

    # 1.2.4 Apply scissor (target 1.1 eV)
    desired_gap_eV = 1.1
    desired_gap = desired_gap_eV * eV
    print(f"  Desired gap: {desired_gap:.6f} Ha = {desired_gap_eV} eV")

    eband_after = bandlib.apply_scissor(
        epsilon, dos, data.nelect, eband_before, desired_gap, data.dosweight
    )

    # Verify new gap
    epsilon2, dos2, _, _ = bandlib.BTPDOS(eband_after, vvband)
    fermi2 = bandlib.solve_for_mu(
        epsilon2, dos2, data.nelect, 0.0, data.dosweight, refine=False, try_center=False
    )
    pos2 = np.abs(epsilon2 - fermi2).argmin()
    for i in range(pos2, -1, -1):
        if dos2[i] != 0:
            vbm2 = epsilon2[i]
            break
    for i in range(pos2, len(dos2)):
        if dos2[i] != 0:
            cbm2 = epsilon2[i]
            break
    actual_gap_after = cbm2 - vbm2
    print(f"  Actual gap after: {actual_gap_after:.6f} Ha = {actual_gap_after / eV:.4f} eV")

    # 1.2.5 Save reference data
    np.savez(
        output_file,
        # Input data
        epsilon=epsilon,
        dos=dos,
        N0=data.nelect,
        eband_before=eband_before,
        desired_gap=desired_gap,
        dosweight=data.dosweight,
        # Output data
        eband_after=eband_after,
        # For verification
        current_gap=current_gap,
        actual_gap_after=actual_gap_after,
        fermi=fermi,
    )

    print(f"\nSaved: {output_file}")
    print(f"  Keys: {list(np.load(output_file).keys())}")


if __name__ == "__main__":
    main()
