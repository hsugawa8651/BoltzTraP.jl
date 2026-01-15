#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data for Phase 4 (Transport).

Run this script to create reference data for BTPDOS, fermiintegrals,
and calc_Onsager_coefficients comparison.
"""

import numpy as np
import os

import BoltzTraP2.dft
import BoltzTraP2.sphere
import BoltzTraP2.fite
import BoltzTraP2.bandlib
from BoltzTraP2.units import *

from common import get_material_path, OUTPUT_DIR


def generate_simple_transport_test():
    """Generate a simple synthetic test case for transport debugging."""
    print("\n" + "=" * 60)
    print("Simple Transport Test Case")
    print("=" * 60)

    # Simple parabolic band: E(k) = k^2
    # DOS ~ sqrt(E) for 3D parabolic band
    print("\n1. Generating simple parabolic band...")

    # Generate energies on a grid
    ngrid = 10
    eband = np.zeros((1, ngrid**3))  # 1 band
    vvband = np.zeros((1, 3, 3, ngrid**3))

    idx = 0
    for i in range(ngrid):
        for j in range(ngrid):
            for k in range(ngrid):
                # k-point in [-0.5, 0.5)
                kx = (i / ngrid) - 0.5
                ky = (j / ngrid) - 0.5
                kz = (k / ngrid) - 0.5

                # Energy: E = k^2 (in arbitrary units)
                eband[0, idx] = kx**2 + ky**2 + kz**2

                # Velocity: v = dE/dk = 2k
                vx, vy, vz = 2*kx, 2*ky, 2*kz

                # v⊗v
                vvband[0, 0, 0, idx] = vx * vx
                vvband[0, 0, 1, idx] = vx * vy
                vvband[0, 0, 2, idx] = vx * vz
                vvband[0, 1, 0, idx] = vy * vx
                vvband[0, 1, 1, idx] = vy * vy
                vvband[0, 1, 2, idx] = vy * vz
                vvband[0, 2, 0, idx] = vz * vx
                vvband[0, 2, 1, idx] = vz * vy
                vvband[0, 2, 2, idx] = vz * vz

                idx += 1

    print(f"   eband shape: {eband.shape}")
    print(f"   vvband shape: {vvband.shape}")
    print(f"   Energy range: [{eband.min():.4f}, {eband.max():.4f}]")

    # Compute BTPDOS
    print("\n2. Computing BTPDOS...")
    epsilon, dos, vvdos, *_ = BoltzTraP2.bandlib.BTPDOS(eband, vvband, npts=100)
    print(f"   epsilon: {epsilon.shape}")
    print(f"   dos: {dos.shape}")
    print(f"   vvdos: {vvdos.shape}")

    # Simple parameters
    Tr = np.array([300.0])  # Single temperature
    mur = np.linspace(0.0, 0.5, 20)  # Chemical potential range
    vuc = 1.0  # Unit volume
    dosweight = 2.0

    # Fermi integrals
    print("\n3. Computing Fermi integrals...")
    N, L0, L1, L2, *_ = BoltzTraP2.bandlib.fermiintegrals(
        epsilon, dos, vvdos, mur, Tr, dosweight=dosweight
    )
    print(f"   L0 shape: {L0.shape}")

    # Onsager coefficients
    print("\n4. Computing Onsager coefficients...")
    sigma, S, kappa, *_ = BoltzTraP2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, vuc
    )
    print(f"   sigma shape: {sigma.shape}")

    # Save
    print("\n5. Saving to simple_transport.npz...")
    save_dict = {
        'eband': eband,
        'vvband': vvband,
        'epsilon': epsilon,
        'dos': dos,
        'vvdos': vvdos,
        'Tr': Tr,
        'mur': mur,
        'vuc': vuc,
        'dosweight': dosweight,
        'N': N,
        'L0': L0,
        'L1': L1,
        'L2': L2,
        'sigma': sigma,
        'S': S,
        'kappa': kappa,
    }
    np.savez(OUTPUT_DIR / "simple_transport.npz", **save_dict)
    print("Done!")

    return save_dict


def generate_calc_N_solve_mu_reference():
    """Generate reference data for calc_N and solve_for_mu functions.

    Uses Si end-to-end data which contains DOS information.
    """
    print("\n" + "=" * 60)
    print("calc_N and solve_for_mu Reference Data")
    print("=" * 60)

    # Load Si end-to-end data with DOS
    print("\n1. Loading Si end-to-end data...")
    si_data = np.load(OUTPUT_DIR / "si_end2end.npz")
    epsilon = si_data['epsilon']
    dos = si_data['dos']
    fermi = float(si_data['fermi'])
    dosweight = float(si_data['dosweight'])

    print(f"   epsilon shape: {epsilon.shape}")
    print(f"   dos shape: {dos.shape}")
    print(f"   Fermi level: {fermi:.6f} Ha")
    print(f"   dosweight: {dosweight}")

    # Test calc_N at various chemical potentials and temperatures
    print("\n2. Computing calc_N at various (mu, T) points...")
    mu_test = np.linspace(fermi - 0.05, fermi + 0.05, 21)  # 21 points around Fermi
    T_test = np.array([0.0, 100.0, 200.0, 300.0, 400.0, 500.0])  # Various temperatures

    calc_N_results = np.zeros((len(T_test), len(mu_test)))
    for iT, T in enumerate(T_test):
        for imu, mu in enumerate(mu_test):
            calc_N_results[iT, imu] = BoltzTraP2.bandlib.calc_N(
                epsilon, dos, mu, T, dosweight
            )
        print(f"   T={T:5.0f}K: N range [{calc_N_results[iT].min():.4f}, {calc_N_results[iT].max():.4f}]")

    # Test solve_for_mu at various N0 values and temperatures
    print("\n3. Computing solve_for_mu for various (N0, T) points...")

    # Get N0 at Fermi level (number of valence electrons)
    N0_fermi = BoltzTraP2.bandlib.calc_N(epsilon, dos, fermi, 0.0, dosweight)
    print(f"   N0 at Fermi level (T=0): {-N0_fermi:.4f}")

    # Test with different N0 values (doping levels)
    N0_test = np.array([
        -N0_fermi - 0.5,  # hole doped
        -N0_fermi - 0.2,
        -N0_fermi,        # intrinsic
        -N0_fermi + 0.2,
        -N0_fermi + 0.5,  # electron doped
    ])

    # solve_for_mu results for different refinement options
    solve_mu_basic = np.zeros((len(T_test), len(N0_test)))
    solve_mu_refine = np.zeros((len(T_test), len(N0_test)))
    solve_mu_center = np.zeros((len(T_test), len(N0_test)))

    for iT, T in enumerate(T_test):
        for iN, N0 in enumerate(N0_test):
            try:
                solve_mu_basic[iT, iN] = BoltzTraP2.bandlib.solve_for_mu(
                    epsilon, dos, N0, T, dosweight, refine=False, try_center=False
                )
            except Exception as e:
                solve_mu_basic[iT, iN] = np.nan
                print(f"   Warning: solve_for_mu failed for T={T}, N0={N0}: {e}")

            try:
                solve_mu_refine[iT, iN] = BoltzTraP2.bandlib.solve_for_mu(
                    epsilon, dos, N0, T, dosweight, refine=True, try_center=False
                )
            except Exception as e:
                solve_mu_refine[iT, iN] = np.nan

            try:
                solve_mu_center[iT, iN] = BoltzTraP2.bandlib.solve_for_mu(
                    epsilon, dos, N0, T, dosweight, refine=True, try_center=True
                )
            except Exception as e:
                solve_mu_center[iT, iN] = np.nan

        print(f"   T={T:5.0f}K: mu_refine = {solve_mu_refine[iT, 2]:.6f} Ha (N0={N0_test[2]:.4f})")

    # Save results
    print("\n4. Saving to calc_N_solve_mu.npz...")
    save_dict = {
        'epsilon': epsilon,
        'dos': dos,
        'fermi': fermi,
        'dosweight': dosweight,
        # calc_N test data
        'mu_test': mu_test,
        'T_test': T_test,
        'calc_N_results': calc_N_results,
        # solve_for_mu test data
        'N0_test': N0_test,
        'solve_mu_basic': solve_mu_basic,
        'solve_mu_refine': solve_mu_refine,
        'solve_mu_center': solve_mu_center,
    }
    np.savez(OUTPUT_DIR / "calc_N_solve_mu.npz", **save_dict)
    print("Done!")

    return save_dict


if __name__ == "__main__":

    # Generate simple test case first
    generate_simple_transport_test()

    # Generate calc_N and solve_for_mu reference (requires si_end2end.npz)
    generate_calc_N_solve_mu_reference()
