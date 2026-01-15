#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Generate reference data for DFT loader validation.

This script loads DFT data using Python BoltzTraP2 loaders and saves
the key data as .npz files for comparison with Julia implementations.

Test Materials Policy:
- Reftest: All materials (including magnetic) for comprehensive validation
- Docs: Only non-magnetic examples (v0.1 scope)
"""

import os

import numpy as np

import BoltzTraP2.dft as dft
from BoltzTraP2.units import Angstrom  # 1 Angstrom in Bohr

from common import get_material_path, OUTPUT_DIR


def save_loader_reference(loader_data, output_name, source_format):
    """Save loader output as reference data.

    Parameters
    ----------
    loader_data : DFTData-like object
        Object with atoms, kpoints, ebands, fermi, dosweight, nelect attributes
    output_name : str
        Output filename without extension (e.g., 'vasp_si')
    source_format : str
        Source format name (e.g., 'VASP', 'Wien2k')
    """
    # Get lattice vectors from ASE Atoms (row convention -> transpose for column convention)
    # Convert from Angstrom to Bohr to match Julia's internal units
    lattvec = loader_data.atoms.get_cell()[:].T * Angstrom  # 3x3, columns are vectors, in Bohr

    # Get atomic positions and symbols
    positions = loader_data.atoms.get_scaled_positions()  # (natom, 3) fractional
    symbols = loader_data.atoms.get_chemical_symbols()  # list of strings

    # Encode symbols as bytes for NPZ.jl compatibility
    # Join with comma and encode as uint8 array
    symbols_str = ",".join(symbols)
    symbols_encoded = np.array(list(symbols_str.encode('utf-8')), dtype=np.uint8)

    # Get band structure data
    kpoints = loader_data.kpoints  # (nkpts, 3)
    ebands = loader_data.ebands  # (nbands, nkpts) or (nbands, nkpts, nspin)
    fermi = loader_data.fermi
    dosweight = loader_data.dosweight
    nelect = getattr(loader_data, 'nelect', -1.0)

    # Determine nspin from dosweight
    nspin = 1 if dosweight == 2.0 else 2

    # Handle magmom if present
    magmom = getattr(loader_data, 'magmom', None)
    has_magmom = magmom is not None

    output_path = OUTPUT_DIR / f"{output_name}.npz"

    save_dict = {
        # Structure
        'lattvec': lattvec,
        'positions': positions,
        'symbols': symbols_encoded,
        # K-points and bands
        'kpoints': kpoints,
        'ebands': ebands,
        'fermi': np.array(fermi),
        'dosweight': np.array(dosweight),
        'nelect': np.array(nelect),
        'nspin': np.array(nspin),
        # Metadata
        'source_format': np.array(list(source_format.encode('utf-8')), dtype=np.uint8),
    }

    if has_magmom:
        save_dict['magmom'] = np.asarray(magmom)

    np.savez(output_path, **save_dict)
    print(f"  Generated: {output_name}.npz")
    print(f"    lattvec shape: {lattvec.shape}")
    print(f"    positions shape: {positions.shape}")
    print(f"    kpoints shape: {kpoints.shape}")
    print(f"    ebands shape: {ebands.shape}")
    print(f"    fermi: {fermi:.6f} Ha")
    print(f"    dosweight: {dosweight}, nspin: {nspin}")
    print(f"    nelect: {nelect}")
    if has_magmom:
        print(f"    magmom shape: {np.asarray(magmom).shape}")


def generate_vasp_reference():
    """Generate reference data for VASP loader.

    Test materials:
    - Si: Non-magnetic semiconductor (docs example)
    - PbTe: Non-magnetic thermoelectric (docs example)
    """
    print("\n=== VASP Loader Reference ===")

    test_cases = [
        # (directory_name, output_name, description)
        ("Si.vasp", "vasp_si", "Si - non-magnetic semiconductor"),
        ("PbTe.vasp.unpolarized", "vasp_pbte", "PbTe - non-magnetic thermoelectric"),
    ]

    for dirname, output_name, description in test_cases:
        path = get_material_path(dirname)
        if not os.path.exists(path):
            print(f"  Skipping {dirname}: not found")
            continue

        print(f"\n  Loading {dirname} ({description})...")
        try:
            data = dft.DFTData(str(path))
            save_loader_reference(data, output_name, "VASP")
        except Exception as e:
            print(f"  Error loading {dirname}: {e}")


def generate_qe_reference():
    """Generate reference data for Quantum ESPRESSO loader.

    Test materials:
    - Si: Non-magnetic semiconductor (docs example)
    - Fe: Collinear magnetic metal (reftest only)
    - CrI3: Antiferromagnetic (reftest only)

    Note: QE loader expects directory containing *.xml file,
    so we use the 'out/' subdirectory.
    """
    print("\n=== QE Loader Reference ===")

    test_cases = [
        # (directory_name, output_name, description)
        ("Si.ESPRESSO/out", "qe_si", "Si - non-magnetic semiconductor"),
        ("Fe.ESPRESSO.collinear/out", "qe_fe", "Fe - collinear magnetic"),
        ("CrI3.ESPRESSO.antiferro/out", "qe_cri3", "CrI3 - antiferromagnetic"),
        ("nitinol.ESPRESSO/out", "qe_nitinol", "NiTi - monoclinic structure"),
    ]

    for dirname, output_name, description in test_cases:
        path = get_material_path(dirname)
        if not os.path.exists(path):
            print(f"  Skipping {dirname}: not found")
            continue

        print(f"\n  Loading {dirname} ({description})...")
        try:
            data = dft.DFTData(str(path))
            save_loader_reference(data, output_name, "QE")
        except Exception as e:
            print(f"  Error loading {dirname}: {e}")


def generate_wien2k_reference():
    """Generate reference data for Wien2k loader.

    Test materials:
    - Si: Non-magnetic semiconductor (docs example)
    - CoSb3: Skutterudite (reftest only)
    - Bi2Te3: Topological insulator with SOC (reftest only)
    """
    print("\n=== Wien2k Loader Reference ===")

    test_cases = [
        ("Si", "wien2k_si", "Si - non-magnetic semiconductor"),
        ("CoSb3", "wien2k_cosb3", "CoSb3 - skutterudite"),
        ("Bi2Te3", "wien2k_bi2te3", "Bi2Te3 - topological insulator (SOC)"),
    ]

    for dirname, output_name, description in test_cases:
        path = get_material_path(dirname)
        if not os.path.exists(path):
            print(f"  Skipping {dirname}: not found")
            continue

        print(f"\n  Loading {dirname} ({description})...")
        try:
            data = dft.DFTData(str(path))
            save_loader_reference(data, output_name, "Wien2k")
        except Exception as e:
            print(f"  Error loading {dirname}: {e}")


def generate_gene_reference():
    """Generate reference data for GENE/Generic loader.

    Test materials:
    - Si: Non-magnetic semiconductor (converted from VASP)
    """
    print("\n=== GENE Loader Reference ===")

    test_cases = [
        ("Si.GENE", "gene_si", "Si - non-magnetic semiconductor"),
    ]

    for dirname, output_name, description in test_cases:
        path = get_material_path(dirname)
        if not os.path.exists(path):
            print(f"  Skipping {dirname}: not found")
            continue

        print(f"\n  Loading {dirname} ({description})...")
        try:
            data = dft.DFTData(str(path))
            save_loader_reference(data, output_name, "GENE")
        except Exception as e:
            print(f"  Error loading {dirname}: {e}")


def generate_abinit_reference():
    """Generate reference data for ABINIT loader.

    Test materials:
    - Si.abinit: Non-magnetic semiconductor (docs example)
    """
    print("\n=== ABINIT Loader Reference ===")

    test_cases = [
        ("Si.abinit", "abinit_si", "Si - non-magnetic semiconductor"),
    ]

    for dirname, output_name, description in test_cases:
        path = get_material_path(dirname)
        if not os.path.exists(path):
            print(f"  Skipping {dirname}: not found")
            continue

        print(f"\n  Loading {dirname} ({description})...")
        try:
            data = dft.DFTData(str(path))
            save_loader_reference(data, output_name, "ABINIT")
        except Exception as e:
            print(f"  Error loading {dirname}: {e}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate loader reference data")
    parser.add_argument(
        '--data-dir', type=str, default=None,
        help='Path to BoltzTraP2 data directory (passed to common.py)'
    )
    parser.add_argument(
        'formats',
        nargs='*',
        default=['vasp'],
        choices=['vasp', 'qe', 'wien2k', 'gene', 'abinit', 'all'],
        help='Formats to generate (default: vasp)'
    )
    args = parser.parse_args()

    formats = args.formats
    if 'all' in formats:
        formats = ['vasp', 'qe', 'wien2k', 'gene', 'abinit']

    generators = {
        'vasp': generate_vasp_reference,
        'qe': generate_qe_reference,
        'wien2k': generate_wien2k_reference,
        'gene': generate_gene_reference,
        'abinit': generate_abinit_reference,
    }

    for fmt in formats:
        generators[fmt]()

    print("\nDone!")
