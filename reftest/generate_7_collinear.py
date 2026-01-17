#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Generate reference data for collinear magnetic calculations.

This script loads PbTe.vasp.collinear from BoltzTraP2 and generates
reference data for validating the Julia collinear implementation.

Usage:
    python generate_7_collinear.py
    python generate_7_collinear.py --data-dir /path/to/BoltzTraP2/data
"""

import numpy as np

import BoltzTraP2.dft as dft
import BoltzTraP2.sphere as sphere
import BoltzTraP2.fite as fite
import BoltzTraP2.bandlib as bandlib
from BoltzTraP2.units import Angstrom, eV

from common import get_material_path, OUTPUT_DIR


def generate_pbte_collinear():
    """Generate reference data for PbTe collinear calculation."""
    print("=" * 60)
    print("Generating PbTe collinear reference data")
    print("=" * 60)

    # Load data
    material_path = get_material_path("PbTe.vasp.collinear")
    print(f"Loading from: {material_path}")

    data = dft.DFTData(str(material_path))

    # Print basic info
    print(f"\n--- Basic Info ---")
    print(f"magmom: {data.magmom}")
    print(f"magmom shape: {data.magmom.shape}")
    print(f"dosweight: {data.dosweight}")
    print(f"ebands shape: {data.ebands.shape}")
    print(f"fermi: {data.fermi} Ha = {data.fermi / eV} eV")
    print(f"nelect: {data.nelect}")

    # Get lattice vectors (convert to column convention, Bohr)
    lattvec = data.atoms.get_cell()[:].T * Angstrom  # 3x3, columns are vectors

    # Get atomic positions and types
    positions = data.atoms.get_scaled_positions()  # (natom, 3) fractional
    symbols = data.atoms.get_chemical_symbols()
    symbols_str = ",".join(symbols)
    symbols_encoded = np.array(list(symbols_str.encode('utf-8')), dtype=np.uint8)

    # Get atomic numbers for types
    types = np.array([atom.number for atom in data.atoms])

    print(f"\n--- Structure ---")
    print(f"lattvec shape: {lattvec.shape}")
    print(f"positions shape: {positions.shape}")
    print(f"symbols: {symbols}")
    print(f"types: {types}")

    # Symmetry analysis with magmom
    print(f"\n--- Symmetry Analysis ---")
    nrotations = sphere.calc_nrotations(data.atoms, data.magmom)
    print(f"nrotations: {nrotations}")

    # Compute equivalences
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, len(data.kpoints))
    neq = len(equivalences)
    print(f"n_equivalences: {neq}")

    # Interpolation
    print(f"\n--- Interpolation ---")
    coeffs = fite.fitde3D(data, equivalences)
    print(f"coeffs shape: {coeffs.shape}")

    # Reconstruct bands on dense grid (uses equivalences to determine grid)
    eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
    print(f"eband shape: {eband.shape}")
    print(f"vvband shape: {vvband.shape}")

    # Transport integration
    print(f"\n--- Transport Integration ---")
    temps = np.array([300.0])
    epsilon, dos, vvdos, cdos = bandlib.BTPDOS(
        eband, vvband, erange=None, npts=10000
    )
    print(f"epsilon shape: {epsilon.shape}")
    print(f"dos shape: {dos.shape}")

    mur = np.array([data.fermi])

    N, L0, L1, L2, L11 = bandlib.fermiintegrals(
        epsilon, dos, vvdos, mur=mur, Tr=temps, dosweight=data.dosweight
    )

    # Unit cell volume (in Bohr^3)
    vuc = np.abs(np.linalg.det(lattvec))
    print(f"vuc: {vuc} Bohr^3")

    sigma, seebeck, kappa, hall = bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, temps, vuc
    )
    print(f"sigma shape: {sigma.shape}")
    print(f"seebeck shape: {seebeck.shape}")
    print(f"kappa shape: {kappa.shape}")

    # Save reference data
    output_file = OUTPUT_DIR / "pbte_collinear.npz"
    print(f"\n--- Saving to {output_file} ---")

    np.savez(
        output_file,
        # Structure
        lattvec=lattvec,
        positions=positions,
        symbols=symbols_encoded,
        types=types,
        # K-points
        kpoints=data.kpoints,
        # Band structure
        ebands=data.ebands,
        fermi=data.fermi,
        nelect=data.nelect,
        # Magnetic specific
        magmom=data.magmom,
        dosweight=data.dosweight,
        # Symmetry
        nrotations=nrotations,
        neq=neq,
        # Interpolation
        coeffs=coeffs,
        # Transport (at T=300K, mu=fermi)
        temps=temps,
        sigma=sigma,
        seebeck=seebeck,
        kappa=kappa,
    )

    print("Done!")
    print(f"\nSaved keys: {list(np.load(output_file).keys())}")


if __name__ == "__main__":
    generate_pbte_collinear()
