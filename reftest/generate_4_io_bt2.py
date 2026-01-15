#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Generate reference data for bt2 file validation.

This script reads a Python BoltzTraP2 generated .bt2 file and saves
the key data as a .npz file for comparison with Julia implementation.
"""

import os
import numpy as np

from BoltzTraP2 import serialization

from common import get_material_path, OUTPUT_DIR


def generate_si_bt2_reference():
    """Generate reference for Si.vasp bt2 file."""
    bt2_path = get_material_path("Si.vasp") / "interpolation.bt2"

    if not os.path.exists(bt2_path):
        print(f"Error: bt2 file not found: {bt2_path}")
        return False

    print(f"Loading: {bt2_path}")
    data, equivalences, coeffs, metadata = serialization.load_calculation(bt2_path)

    print(f"  Coefficients shape: {coeffs.shape}")
    print(f"  Equivalences count: {len(equivalences)}")
    print(f"  Format version: {metadata.get('format_version', 'unknown')}")

    # Get lattice vectors from ASE Atoms
    lattvec = data.atoms.get_cell()[:]

    # Save reference data
    # Note: NPZ.jl doesn't support Unicode strings, so encode as bytes
    output_path = OUTPUT_DIR / 'si_bt2_reference.npz'
    np.savez(
        output_path,
        # Metadata for provenance check (as bytes for NPZ.jl compatibility)
        format_version=np.array(list(metadata.get("format_version", "").encode('utf-8')), dtype=np.uint8),
        program_version=np.array(list(metadata.get("program_version", "").encode('utf-8')), dtype=np.uint8),
        source=np.array(list(metadata.get("source", "").encode('utf-8')), dtype=np.uint8),
        ctime=np.array(list(metadata.get("ctime", "").encode('utf-8')), dtype=np.uint8),
        # Coefficients
        coeffs_real=coeffs.real,
        coeffs_imag=coeffs.imag,
        coeffs_shape=np.array(coeffs.shape),
        # Equivalences
        n_equivalences=len(equivalences),
        equiv_0=equivalences[0] if len(equivalences) > 0 else np.array([]),
        # Lattice vectors (row convention in Python)
        lattvec=lattvec,
        # DFT data
        fermi=data.fermi,
        nelect=data.nelect,
        dosweight=data.dosweight,
        n_bands=data.ebands.shape[0],
        n_kpoints=data.ebands.shape[1],
    )
    print(f"Generated: {output_path}")
    return True


if __name__ == "__main__":
    import sys
    success = generate_si_bt2_reference()
    sys.exit(0 if success else 1)
