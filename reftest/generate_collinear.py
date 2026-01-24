#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data for collinear magnetic calculations.

This script generates reference data for validating collinear magnetic
calculations using test materials from BoltzTraP2.

Usage:
    python generate_collinear.py                    # Generate all materials
    python generate_collinear.py pbte               # Generate PbTe only
    python generate_collinear.py fe                 # Generate Fe only
    python generate_collinear.py --data-dir /path   # Specify data directory
"""

import argparse
import numpy as np
import os

import BoltzTraP2.dft as dft
import BoltzTraP2.fite as fite
import BoltzTraP2.sphere as sphere
import BoltzTraP2.bandlib as bandlib

from common import get_material_path, OUTPUT_DIR


# Material configurations
MATERIALS = {
    'pbte': {
        'path': 'PbTe.vasp.collinear',
        'subdir': None,  # No subdirectory
        'temperatures': np.array([300.0]),
        'npts_dos': 10000,
        'mu_npoints': 1,  # Single point at Fermi
        'mu_margin': 0.0,  # Not used when npoints=1
        'output_file': 'pbte_collinear.npz',
        'description': 'PbTe collinear (VASP)',
    },
    'fe': {
        'path': 'Fe.ESPRESSO.collinear',
        'subdir': 'out',  # QE uses 'out' subdirectory
        'temperatures': np.array([200.0, 300.0, 400.0, 1000.0]),
        'npts_dos': 500,
        'mu_npoints': 200,  # Range around Fermi
        'mu_margin': 0.1,  # ±0.1 Ha
        'output_file': 'fe_collinear_e2e.npz',
        'description': 'Fe collinear (QE)',
    },
}


def generate_collinear(material_name):
    """Generate reference data for a collinear magnetic material.

    Args:
        material_name: Key in MATERIALS dict ('pbte' or 'fe')
    """
    config = MATERIALS[material_name]

    print("=" * 60)
    print(f"Generating {config['description']} reference data")
    print("=" * 60)

    # Load DFT data
    material_path = get_material_path(config['path'])
    if config['subdir']:
        data_path = str(material_path) + "/" + config['subdir']
    else:
        data_path = str(material_path)

    data = dft.DFTData(data_path)

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi:.6f} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")
    print(f"  magmom: {data.magmom}")

    # Get lattice vectors (column convention, Angstrom from ASE)
    # Important: Keep in Angstrom for consistency with other reftest scripts.
    # Julia test must convert Å→Bohr when loading.
    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))
    print(f"  Unit cell volume: {vuc:.4f} Å³")

    # Get symmetry info
    nrotations = sphere.calc_nrotations(data.atoms, data.magmom)
    print(f"  Rotations: {nrotations}")

    # Get equivalences
    nkpt_target = len(data.kpoints)
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)
    print(f"  Equivalences: {len(equivalences)}")

    # Fit coefficients
    coeffs = fite.fitde3D(data, equivalences)
    print(f"  Coefficients: {coeffs.shape}")

    # Test reconstruction at original k-points
    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )
    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    # Build save dictionary with common fields
    save_dict = {
        # Input data
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'magmom': data.magmom,
        'vuc': vuc,
        # Symmetry and equivalence info
        'nrotations': nrotations,
        'n_equivalences': len(equivalences),
        'neq': len(equivalences),  # Alias for compatibility
        # Coefficients
        'coeffs': coeffs,
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        # Reconstruction test
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    # Add structure data (useful for some tests)
    positions = data.atoms.get_scaled_positions()
    symbols = data.atoms.get_chemical_symbols()
    symbols_str = ",".join(symbols)
    symbols_encoded = np.array(list(symbols_str.encode('utf-8')), dtype=np.uint8)
    types = np.array([atom.number for atom in data.atoms])

    save_dict.update({
        'positions': positions,
        'symbols': symbols_encoded,
        'types': types,
    })

    # Transport calculation
    print(f"\n  Computing transport...")
    Tr = config['temperatures']

    # Reconstruct bands via FFT
    eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
    print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

    # Compute DOS
    npts = config['npts_dos']
    epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=npts)
    print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

    # Define mu range
    if config['mu_npoints'] == 1:
        mur = np.array([data.fermi])
    else:
        mu_margin = config['mu_margin']
        mur = np.linspace(data.fermi - mu_margin, data.fermi + mu_margin,
                          config['mu_npoints'])
    print(f"    Tr: {Tr}, mur: {len(mur)} points")

    # Compute Fermi integrals
    N, L0, L1, L2, *_ = bandlib.fermiintegrals(
        epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
    )
    print(f"    L0: {L0.shape}")

    # Compute Onsager coefficients
    sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, vuc
    )
    print(f"    sigma: {sigma.shape}")

    # Solve for intrinsic mu at each temperature
    mu0 = np.zeros(len(Tr))
    for iT, T in enumerate(Tr):
        mu0[iT] = bandlib.solve_for_mu(
            epsilon, dos, data.nelect, T, data.dosweight,
            refine=True, try_center=True
        )
    print(f"    mu0: {mu0}")

    # Add transport data to save dictionary
    save_dict.update({
        # getBTPbands output
        'eband': eband,
        'vvband': vvband,
        # BTPDOS output
        'epsilon': epsilon,
        'dos': dos,
        'vvdos': vvdos,
        # Fermi integrals
        'Tr': Tr,
        'mur': mur,
        'N': N,
        'L0': L0,
        'L1': L1,
        'L2': L2,
        # Onsager coefficients (use consistent naming)
        'sigma': sigma,
        'S': S,
        'seebeck': S,  # Alias for compatibility
        'kappa': kappa,
        # solve_for_mu
        'mu0': mu0,
        # Config used
        'temps': Tr,  # Alias for compatibility
        'npts_dos': npts,
    })

    # Save reference data
    output_path = os.path.join(OUTPUT_DIR, config['output_file'])
    np.savez(output_path, **save_dict)

    print(f"\n  Saved: {output_path}")
    print(f"  Keys: {sorted(save_dict.keys())}")
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Generate collinear magnetic reference data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python generate_collinear.py              # Generate all materials
    python generate_collinear.py pbte         # Generate PbTe only
    python generate_collinear.py fe           # Generate Fe only
    python generate_collinear.py pbte fe      # Generate both explicitly
        """
    )
    parser.add_argument('materials', nargs='*', default=[],
                        help=f'Materials to generate: {", ".join(MATERIALS.keys())} (default: all)')
    parser.add_argument('--data-dir', type=str, default=None,
                        help='Path to BoltzTraP2 data directory')

    args = parser.parse_args()

    # Handle empty list (generate all)
    materials = args.materials if args.materials else list(MATERIALS.keys())

    # Validate materials
    for m in materials:
        if m not in MATERIALS:
            parser.error(f"Unknown material: {m}. Choose from: {', '.join(MATERIALS.keys())}")

    for material in materials:
        generate_collinear(material)
        print()

    print("Done!")


if __name__ == '__main__':
    main()
