#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Convert DFT data to GENE format.

This script converts any DFT data readable by BoltzTraP2 (VASP, QE, Wien2k,
ABINIT) to GENE format (.structure and .energy files).

Usage:
    python convert_dft_to_gene.py <input_dir> [--output <output_dir>] [--name <name>]

Examples:
    python convert_dft_to_gene.py /path/to/Si.vasp
    python convert_dft_to_gene.py /path/to/Si.vasp --output /path/to/Si.GENE
    python convert_dft_to_gene.py /path/to/Si.vasp --name Si
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# BoltzTraP2 units and loader
from BoltzTraP2.units import Angstrom, Ha  # 1 Angstrom in Bohr, 1 Hartree in atomic units
import BoltzTraP2.dft as dft

# 1 Rydberg = 0.5 Hartree (in atomic units, Ha = 1.0)
Rydberg = 0.5


def write_gene_structure(filepath, lattvec, positions, species, title="Structure"):
    """Write GENE .structure file.

    Parameters
    ----------
    filepath : Path
        Output file path
    lattvec : ndarray
        Lattice vectors (3x3), columns are vectors, in Bohr
    positions : ndarray
        Cartesian coordinates (natom, 3), in Bohr
    species : list of str
        Element symbols
    title : str
        Title line for the file
    """
    with open(filepath, 'w') as f:
        # Line 1: Title
        f.write(f"{title}\n")

        # Lines 2-4: Lattice vectors (row-major, Bohr)
        # lattvec has columns as vectors, so we write rows
        for i in range(3):
            f.write(f"{lattvec[0, i]:20.14f} {lattvec[1, i]:20.14f} {lattvec[2, i]:20.14f}\n")

        # Line 5: Number of atoms
        f.write(f"{len(species)}\n")

        # Lines 6+: Element + Cartesian coordinates (Bohr)
        for i, elem in enumerate(species):
            f.write(f"{elem:2s} {positions[i, 0]:20.14f} {positions[i, 1]:20.14f} {positions[i, 2]:20.14f}\n")


def write_gene_energy(filepath, kpoints, ebands, fermi, nspin=1, title="Band energies"):
    """Write GENE .energy file.

    Parameters
    ----------
    filepath : Path
        Output file path
    kpoints : ndarray
        K-points (nkpts, 3), fractional coordinates
    ebands : ndarray
        Band energies (nbands, nkpts) or (nbands, nkpts, nspin), in Hartree
    fermi : float
        Fermi level in Hartree
    nspin : int
        Number of spin channels (1 or 2)
    title : str
        Title line for the file
    """
    nkpts = kpoints.shape[0]

    # Handle different ebands shapes
    if ebands.ndim == 2:
        nbands = ebands.shape[0]
        ebands_3d = ebands[:, :, np.newaxis]  # (nbands, nkpts, 1)
    else:
        nbands = ebands.shape[0]
        ebands_3d = ebands

    # Convert Hartree to Rydberg for GENE format
    fermi_ry = fermi / Rydberg

    with open(filepath, 'w') as f:
        # Line 1: Title
        f.write(f"{title}\n")

        # Line 2: nk nspin efermi(Ry)
        f.write(f"{nkpts} {nspin} {fermi_ry:20.14f}\n")

        # For each spin channel
        for ispin in range(nspin):
            # For each k-point
            for ik in range(nkpts):
                # Header: kx ky kz nband
                kx, ky, kz = kpoints[ik]
                f.write(f"{kx:20.14f} {ky:20.14f} {kz:20.14f} {nbands}\n")

                # Band energies in Rydberg
                for ib in range(nbands):
                    energy_ry = ebands_3d[ib, ik, ispin] / Rydberg
                    f.write(f"{energy_ry:20.14f}\n")


def convert_dft_to_gene(input_dir, output_dir=None, name=None, verify=True):
    """Convert DFT data to GENE format.

    Parameters
    ----------
    input_dir : Path or str
        Input directory containing DFT data (VASP, QE, Wien2k, ABINIT)
    output_dir : Path or str, optional
        Output directory for GENE files. Default: <input_dir>/../<name>.GENE
    name : str, optional
        Base name for output files. Default: derived from input directory
    verify : bool
        If True, verify output by loading with BoltzTraP2

    Returns
    -------
    Path
        Path to output directory
    """
    input_dir = Path(input_dir)

    # Derive name from input directory
    if name is None:
        name = input_dir.name
        # Handle subdirectories like Si.ESPRESSO/out
        if name == 'out':
            name = input_dir.parent.name
        # Remove common suffixes like .vasp, .vasp.unpolarized, .ESPRESSO, etc.
        suffixes_to_remove = [
            '.vasp.unpolarized', '.vasp.polarized',
            '.vasp', '.VASP',
            '.ESPRESSO', '.espresso',
            '.abinit', '.ABINIT',
            '.W2K', '.w2k',
        ]
        for suffix in suffixes_to_remove:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break

    # Default output directory
    if output_dir is None:
        output_dir = input_dir.parent / f"{name}.GENE"
    else:
        output_dir = Path(output_dir)

    print("=" * 60)
    print(f"Converting DFT to GENE format")
    print(f"  Input:  {input_dir}")
    print(f"  Output: {output_dir}")
    print(f"  Name:   {name}")
    print("=" * 60)

    # Load DFT data
    print(f"\nLoading DFT data...")
    data = dft.DFTData(str(input_dir))

    # Get lattice vectors (convert from Angstrom to Bohr)
    lattvec = data.atoms.get_cell()[:].T * Angstrom  # (3, 3), columns are vectors, Bohr

    # Get atomic positions in Cartesian coordinates (Bohr)
    positions = data.atoms.get_positions() * Angstrom  # (natom, 3), Bohr

    # Get species
    species = data.atoms.get_chemical_symbols()

    # Get k-points and bands
    kpoints = data.kpoints  # (nkpts, 3), fractional
    ebands = data.ebands    # (nbands, nkpts), Hartree
    fermi = data.fermi      # Hartree

    # Determine nspin from dosweight
    nspin = 1 if data.dosweight == 2.0 else 2

    print(f"  Species: {species}")
    print(f"  K-points: {kpoints.shape[0]}")
    print(f"  Bands: {ebands.shape[0]}")
    print(f"  Fermi level: {fermi:.6f} Ha")
    print(f"  dosweight: {data.dosweight}, nspin: {nspin}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write structure file
    structure_file = output_dir / f"{name}.structure"
    print(f"\nWriting {structure_file}...")
    write_gene_structure(
        structure_file, lattvec, positions, species,
        title=f"{name} structure"
    )

    # Write energy file
    energy_file = output_dir / f"{name}.energy"
    print(f"Writing {energy_file}...")
    write_gene_energy(
        energy_file, kpoints, ebands, fermi, nspin,
        title=f"{name} band energies"
    )

    print(f"\nGenerated GENE files in {output_dir}")

    # Verify by loading with BoltzTraP2
    if verify:
        print("\nVerifying with BoltzTraP2 loader...")
        try:
            verify_data = dft.DFTData(str(output_dir))
            print(f"  ✓ Load successful")

            # Compare with original
            fermi_diff = abs(verify_data.fermi - data.fermi)
            ebands_diff = np.max(np.abs(verify_data.ebands - data.ebands))
            print(f"  Fermi diff: {fermi_diff:.2e} Ha")
            print(f"  Max ebands diff: {ebands_diff:.2e} Ha")

            if fermi_diff < 1e-10 and ebands_diff < 1e-10:
                print("  ✓ Verification passed!")
            else:
                print("  ⚠ Verification shows differences")

        except Exception as e:
            print(f"  ✗ Verification failed: {e}")

    print("=" * 60)
    return output_dir


def main():
    parser = argparse.ArgumentParser(
        description="Convert DFT data to GENE format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_dft_to_gene.py /path/to/Si.vasp
  python convert_dft_to_gene.py /path/to/Si.ESPRESSO/out --name Si
  python convert_dft_to_gene.py /path/to/data/Si --output ./Si.GENE

Supported input formats:
  - VASP (vasprun.xml)
  - Quantum ESPRESSO (*.xml)
  - Wien2k (case.struct, case.energy)
  - ABINIT (*_GSR.nc)
"""
    )
    parser.add_argument(
        'input_dir',
        type=str,
        help='Input directory containing DFT data'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output directory for GENE files (default: <input>/../<name>.GENE)'
    )
    parser.add_argument(
        '--name', '-n',
        type=str,
        default=None,
        help='Base name for output files (default: derived from input directory)'
    )
    parser.add_argument(
        '--no-verify',
        action='store_true',
        help='Skip verification step'
    )

    args = parser.parse_args()

    if not Path(args.input_dir).exists():
        print(f"Error: Input directory not found: {args.input_dir}")
        sys.exit(1)

    convert_dft_to_gene(
        args.input_dir,
        output_dir=args.output,
        name=args.name,
        verify=not args.no_verify
    )


if __name__ == "__main__":
    main()
