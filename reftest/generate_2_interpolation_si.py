#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data for Phase 3 (Interpolation).

Run this script to create reference data for compute_phase_factors,
fitde3D, and getBTPbands comparison.
"""

import numpy as np

import BoltzTraP2.dft
import BoltzTraP2.sphere
import BoltzTraP2.fite

from common import get_material_path, OUTPUT_DIR


def generate_si_interpolation():
    """Generate interpolation test data for Si."""
    print("Loading Si DFT data...")
    data = BoltzTraP2.dft.DFTData(str(get_material_path("Si")))

    print(f"  kpoints: {data.kpoints.shape}")
    print(f"  ebands: {data.ebands.shape}")
    print(f"  lattvec: {data.get_lattvec().shape}")

    # Compute equivalences with enhancement factor = 2
    print("Computing equivalence classes...")
    equivalences = BoltzTraP2.sphere.get_equivalences(
        data.atoms, data.magmom, 2 * len(data.kpoints)
    )
    print(f"  Number of equivalence classes: {len(equivalences)}")

    # Compute phase factors
    print("Computing phase factors...")
    tpii = 2j * np.pi
    nk = len(data.kpoints)
    neq = len(equivalences)
    phase = np.zeros((nk, neq), dtype=complex)

    for j, equiv in enumerate(equivalences):
        nstar = len(equiv)
        phase0 = np.exp(tpii * data.kpoints @ equiv.T)
        phase[:, j] = np.sum(phase0, axis=1) / nstar

    print(f"  Phase matrix shape: {phase.shape}")

    # Fit coefficients
    print("Fitting coefficients (fitde3D)...")
    coeffs = BoltzTraP2.fite.fitde3D(data, equivalences)
    print(f"  Coefficients shape: {coeffs.shape}")

    # Reconstruct bands using getBands
    print("Reconstructing bands (getBands)...")
    ebands_reconstructed, vbands = BoltzTraP2.fite.getBands(
        data.kpoints, equivalences, data.get_lattvec(), coeffs
    )
    print(f"  Reconstructed ebands shape: {ebands_reconstructed.shape}")
    print(f"  Velocities shape: {vbands.shape}")

    # Verify reconstruction quality
    max_error = np.max(np.abs(data.ebands - ebands_reconstructed))
    print(f"  Max reconstruction error: {max_error:.2e}")

    # Save all data
    print("Saving to si_interpolation.npz...")
    save_dict = {
        # Input data
        'lattvec': data.get_lattvec(),
        'kpoints': data.kpoints,
        'ebands': data.ebands,  # (nbands, nk)
        'n_equivalences': len(equivalences),
        # Intermediate
        'phase': phase,  # (nk, neq)
        # Output
        'coeffs': coeffs,  # (nbands, neq)
        'ebands_reconstructed': ebands_reconstructed,
        'vbands': vbands,  # (3, nbands, nk)
    }

    # Add equivalences
    for i, equiv in enumerate(equivalences):
        save_dict[f'equiv_{i}'] = equiv

    np.savez(OUTPUT_DIR / "si_interpolation.npz", **save_dict)
    print("Done!")


def generate_simple_test_case():
    """Generate a simple synthetic test case for debugging."""
    print("\nGenerating simple test case...")

    # Simple cubic lattice
    lattvec = np.eye(3) * 5.0  # 5 Angstrom

    # Simple 2x2x2 k-point grid
    kgrid = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                kgrid.append([i/2, j/2, k/2])
    kpoints = np.array(kgrid)

    # Simple parabolic band (free electron)
    # E(k) = k^2 (in units where hbar^2/2m = 1)
    ebands = np.zeros((1, len(kpoints)))
    for i, k in enumerate(kpoints):
        # Wrap to [-0.5, 0.5)
        kw = k - np.floor(k + 0.5)
        ebands[0, i] = np.sum(kw**2)

    # Simple equivalences (just origin for this test)
    # For phase factor test: use known R vectors
    equiv_0 = np.array([[0, 0, 0]])  # Origin
    equiv_1 = np.array([[1, 0, 0], [-1, 0, 0]])  # ±x
    equiv_2 = np.array([[0, 1, 0], [0, -1, 0]])  # ±y
    equiv_3 = np.array([[0, 0, 1], [0, 0, -1]])  # ±z
    equivalences = [equiv_0, equiv_1, equiv_2, equiv_3]

    # Compute phase factors manually
    tpii = 2j * np.pi
    nk = len(kpoints)
    neq = len(equivalences)
    phase = np.zeros((nk, neq), dtype=complex)

    for j, equiv in enumerate(equivalences):
        nstar = len(equiv)
        for ik, k in enumerate(kpoints):
            s = 0.0
            for R in equiv:
                s += np.exp(tpii * np.dot(k, R))
            phase[ik, j] = s / nstar

    print(f"  kpoints: {kpoints.shape}")
    print(f"  ebands: {ebands.shape}")
    print(f"  phase: {phase.shape}")

    # Save
    save_dict = {
        'lattvec': lattvec,
        'kpoints': kpoints,
        'ebands': ebands,
        'n_equivalences': neq,
        'phase': phase,
    }
    for i, equiv in enumerate(equivalences):
        save_dict[f'equiv_{i}'] = equiv

    np.savez(OUTPUT_DIR / "simple_interpolation.npz", **save_dict)
    print("Saved to simple_interpolation.npz")


if __name__ == "__main__":

    # Generate simple test case first (for debugging)
    generate_simple_test_case()

    # Generate Si interpolation data
    generate_si_interpolation()
