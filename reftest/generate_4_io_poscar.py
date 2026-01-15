#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data for POSCAR reading test."""

import numpy as np
import os

# Use ASE for POSCAR reading (same as BoltzTraP2)
import ase.io

from common import get_material_path, OUTPUT_DIR


def generate_poscar_reference():
    """Generate reference data from POSCAR files."""
    print("=" * 60)
    print("POSCAR Reading Reference Data Generation")
    print("=" * 60)

    test_cases = [
        ("Si.vasp", "si_poscar"),
        ("Li.vasp", "li_poscar"),
        ("PbTe.vasp.unpolarized", "pbte_poscar"),
    ]

    for dirname, output_name in test_cases:
        poscar_path = get_material_path(dirname) / "POSCAR"
        if not os.path.exists(poscar_path):
            print(f"Skipping {dirname}: POSCAR not found")
            continue

        print(f"\nProcessing {dirname}...")

        # Read with ASE (same as BoltzTraP2)
        atoms = ase.io.read(poscar_path)

        # Extract data in BoltzTraP2 convention
        # lattvec: column vectors (3x3), atoms.get_cell() returns row vectors
        lattvec = atoms.get_cell().T
        positions = atoms.get_scaled_positions()  # fractional coordinates
        types = atoms.numbers  # atomic numbers
        symbols = atoms.get_chemical_symbols()

        print(f"  Lattice shape: {lattvec.shape}")
        print(f"  Positions shape: {positions.shape}")
        print(f"  Types: {types}")
        print(f"  Symbols: {symbols}")

        # Save reference data
        # Note: symbols saved as bytes to avoid NPZ.jl Unicode issues
        np.savez(
            OUTPUT_DIR / f"{output_name}.npz",
            lattvec=lattvec,
            positions=positions,
            types=types,
            n_atoms=len(atoms),
        )
        print(f"  Saved: {output_name}.npz")


if __name__ == "__main__":
    generate_poscar_reference()
    print("\nDone!")
