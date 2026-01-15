#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data for end-to-end interpolation test.

Compares Python BoltzTraP2 with Julia BoltzTraP.jl.
"""

import numpy as np
import os

import BoltzTraP2.dft as dft
import BoltzTraP2.fite as fite
import BoltzTraP2.sphere as sphere
import BoltzTraP2.bandlib as bandlib

from common import get_material_path, OUTPUT_DIR

# Default transport parameters
DEFAULT_TEMPERATURES = np.array([200.0, 300.0, 400.0])  # Kelvin
DEFAULT_NPTS_DOS = 500  # DOS histogram bins


def generate_si_interpolate(include_transport=True):
    """Generate reference for Si interpolation (end-to-end).

    Args:
        include_transport: If True, also compute transport coefficients.
    """
    print("Generating Si interpolation reference...")

    # Load DFT data
    data_path = get_material_path("Si.vasp")
    data = dft.DFTData(str(data_path))

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")

    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))  # Unit cell volume in Bohr^3
    print(f"  Unit cell volume: {vuc:.4f} Bohr^3")

    # Get equivalences with target ~5000 k-points
    nkpt_target = 5000
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)

    print(f"  Equivalences: {len(equivalences)}")

    # Fit coefficients (using all bands, no energy filter for simplicity)
    coeffs = fite.fitde3D(data, equivalences)

    print(f"  Coefficients: {coeffs.shape}")

    # Test reconstruction at original k-points
    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )

    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    # Build save dictionary
    save_dict = {
        # Input data
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'vuc': vuc,
        # Equivalence info
        'n_equivalences': len(equivalences),
        'equiv_sizes': np.array([len(eq) for eq in equivalences]),
        # Representative points from each equivalence class (first point)
        'equiv_reps': np.array([eq[0] for eq in equivalences]),
        # Output coefficients
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        # Reconstruction test
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    # Transport calculation (E2E integrate)
    if include_transport:
        print("\n  Computing transport (E2E integrate)...")

        # 1. Reconstruct bands via FFT (getBTPbands)
        print("    getBTPbands...")
        eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
        print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

        # 2. Compute DOS and transport DOS (BTPDOS)
        print("    BTPDOS...")
        epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=DEFAULT_NPTS_DOS)
        print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

        # 3. Define temperature and chemical potential ranges
        Tr = DEFAULT_TEMPERATURES
        # Chemical potential range around Fermi level
        mu_margin = 0.1  # Ha
        mu_min = data.fermi - mu_margin
        mu_max = data.fermi + mu_margin
        mur = np.linspace(mu_min, mu_max, 200)
        print(f"    Tr: {Tr}, mur: {len(mur)} points")

        # 4. Compute Fermi integrals
        print("    fermiintegrals...")
        N, L0, L1, L2, *_ = bandlib.fermiintegrals(
            epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
        )
        print(f"    L0: {L0.shape}")

        # 5. Compute Onsager coefficients
        print("    calc_Onsager_coefficients...")
        sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
            L0, L1, L2, mur, Tr, vuc
        )
        print(f"    sigma: {sigma.shape}")

        # 6. Solve for intrinsic mu at each temperature
        print("    solve_for_mu...")
        mu0 = np.zeros(len(Tr))
        for iT, T in enumerate(Tr):
            mu0[iT] = bandlib.solve_for_mu(
                epsilon, dos, data.nelect, T, data.dosweight,
                refine=True, try_center=True
            )
            print(f"      T={T}K: mu0={mu0[iT]:.6f} Ha")

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
            # Onsager coefficients
            'sigma': sigma,
            'S': S,
            'kappa': kappa,
            # solve_for_mu
            'mu0': mu0,
        })

    # Save reference data
    output_path = os.path.join(OUTPUT_DIR, 'si_end2end.npz')
    np.savez(output_path, **save_dict)

    print(f"  Saved: {output_path}")
    return output_path


def generate_pbte_interpolate(include_transport=True):
    """Generate reference for PbTe interpolation (thermoelectric material)."""
    print("\nGenerating PbTe interpolation reference...")

    data_path = get_material_path("PbTe.vasp.unpolarized")
    data = dft.DFTData(str(data_path))

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")

    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))
    print(f"  Unit cell volume: {vuc:.4f} Bohr^3")

    nkpt_target = 5000
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)

    print(f"  Equivalences: {len(equivalences)}")

    coeffs = fite.fitde3D(data, equivalences)

    print(f"  Coefficients: {coeffs.shape}")

    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )

    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    save_dict = {
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'vuc': vuc,
        'n_equivalences': len(equivalences),
        'equiv_sizes': np.array([len(eq) for eq in equivalences]),
        'equiv_reps': np.array([eq[0] for eq in equivalences]),
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    if include_transport:
        print("\n  Computing transport (E2E integrate)...")

        print("    getBTPbands...")
        eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
        print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

        print("    BTPDOS...")
        epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=DEFAULT_NPTS_DOS)
        print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

        Tr = DEFAULT_TEMPERATURES
        mu_margin = 0.1
        mu_min = data.fermi - mu_margin
        mu_max = data.fermi + mu_margin
        mur = np.linspace(mu_min, mu_max, 200)
        print(f"    Tr: {Tr}, mur: {len(mur)} points")

        print("    fermiintegrals...")
        N, L0, L1, L2, *_ = bandlib.fermiintegrals(
            epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
        )
        print(f"    L0: {L0.shape}")

        print("    calc_Onsager_coefficients...")
        sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
            L0, L1, L2, mur, Tr, vuc
        )
        print(f"    sigma: {sigma.shape}")

        print("    solve_for_mu...")
        mu0 = np.zeros(len(Tr))
        for iT, T in enumerate(Tr):
            mu0[iT] = bandlib.solve_for_mu(
                epsilon, dos, data.nelect, T, data.dosweight,
                refine=True, try_center=True
            )
            print(f"      T={T}K: mu0={mu0[iT]:.6f} Ha")

        save_dict.update({
            'eband': eband,
            'vvband': vvband,
            'epsilon': epsilon,
            'dos': dos,
            'vvdos': vvdos,
            'Tr': Tr,
            'mur': mur,
            'N': N,
            'L0': L0,
            'L1': L1,
            'L2': L2,
            'sigma': sigma,
            'S': S,
            'kappa': kappa,
            'mu0': mu0,
        })

    output_path = os.path.join(OUTPUT_DIR, 'pbte_end2end.npz')
    np.savez(output_path, **save_dict)

    print(f"  Saved: {output_path}")
    return output_path


def generate_qe_si_interpolate(include_transport=True):
    """Generate reference for QE Si interpolation (end-to-end).

    Args:
        include_transport: If True, also compute transport coefficients.
    """
    print("\nGenerating QE Si interpolation reference...")

    # Load DFT data from QE output
    data_path = get_material_path("Si.ESPRESSO/out")
    data = dft.DFTData(str(data_path))

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")

    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))  # Unit cell volume in Bohr^3
    print(f"  Unit cell volume: {vuc:.4f} Bohr^3")

    # Get equivalences with target ~5000 k-points
    nkpt_target = 5000
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)

    print(f"  Equivalences: {len(equivalences)}")

    # Fit coefficients (using all bands, no energy filter for simplicity)
    coeffs = fite.fitde3D(data, equivalences)

    print(f"  Coefficients: {coeffs.shape}")

    # Test reconstruction at original k-points
    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )

    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    # Build save dictionary
    save_dict = {
        # Input data
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'vuc': vuc,
        # Equivalence info
        'n_equivalences': len(equivalences),
        'equiv_sizes': np.array([len(eq) for eq in equivalences]),
        # Representative points from each equivalence class (first point)
        'equiv_reps': np.array([eq[0] for eq in equivalences]),
        # Output coefficients
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        # Reconstruction test
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    # Transport calculation (E2E integrate)
    if include_transport:
        print("\n  Computing transport (E2E integrate)...")

        # 1. Reconstruct bands via FFT (getBTPbands)
        print("    getBTPbands...")
        eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
        print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

        # 2. Compute DOS and transport DOS (BTPDOS)
        print("    BTPDOS...")
        epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=DEFAULT_NPTS_DOS)
        print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

        # 3. Define temperature and chemical potential ranges
        Tr = DEFAULT_TEMPERATURES
        # Chemical potential range around Fermi level
        mu_margin = 0.1  # Ha
        mu_min = data.fermi - mu_margin
        mu_max = data.fermi + mu_margin
        mur = np.linspace(mu_min, mu_max, 200)
        print(f"    Tr: {Tr}, mur: {len(mur)} points")

        # 4. Compute Fermi integrals
        print("    fermiintegrals...")
        N, L0, L1, L2, *_ = bandlib.fermiintegrals(
            epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
        )
        print(f"    L0: {L0.shape}")

        # 5. Compute Onsager coefficients
        print("    calc_Onsager_coefficients...")
        sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
            L0, L1, L2, mur, Tr, vuc
        )
        print(f"    sigma: {sigma.shape}")

        # 6. Solve for intrinsic mu at each temperature
        print("    solve_for_mu...")
        mu0 = np.zeros(len(Tr))
        for iT, T in enumerate(Tr):
            mu0[iT] = bandlib.solve_for_mu(
                epsilon, dos, data.nelect, T, data.dosweight,
                refine=True, try_center=True
            )
            print(f"      T={T}K: mu0={mu0[iT]:.6f} Ha")

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
            # Onsager coefficients
            'sigma': sigma,
            'S': S,
            'kappa': kappa,
            # solve_for_mu
            'mu0': mu0,
        })

    # Save reference data
    output_path = os.path.join(OUTPUT_DIR, 'qe_si_end2end.npz')
    np.savez(output_path, **save_dict)

    print(f"  Saved: {output_path}")
    return output_path


def generate_wien2k_si_interpolate(include_transport=True):
    """Generate reference for Wien2k Si interpolation (end-to-end).

    Args:
        include_transport: If True, also compute transport coefficients.
    """
    print("\nGenerating Wien2k Si interpolation reference...")

    # Load DFT data from Wien2k output
    data_path = get_material_path("Si")
    data = dft.DFTData(str(data_path))

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")

    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))  # Unit cell volume in Bohr^3
    print(f"  Unit cell volume: {vuc:.4f} Bohr^3")

    # Get equivalences with target ~5000 k-points
    nkpt_target = 5000
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)

    print(f"  Equivalences: {len(equivalences)}")

    # Fit coefficients (using all bands, no energy filter for simplicity)
    coeffs = fite.fitde3D(data, equivalences)

    print(f"  Coefficients: {coeffs.shape}")

    # Test reconstruction at original k-points
    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )

    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    # Build save dictionary
    save_dict = {
        # Input data
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'vuc': vuc,
        # Equivalence info
        'n_equivalences': len(equivalences),
        'equiv_sizes': np.array([len(eq) for eq in equivalences]),
        # Representative points from each equivalence class (first point)
        'equiv_reps': np.array([eq[0] for eq in equivalences]),
        # Output coefficients
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        # Reconstruction test
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    # Transport calculation (E2E integrate)
    if include_transport:
        print("\n  Computing transport (E2E integrate)...")

        # 1. Reconstruct bands via FFT (getBTPbands)
        print("    getBTPbands...")
        eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
        print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

        # 2. Compute DOS and transport DOS (BTPDOS)
        print("    BTPDOS...")
        epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=DEFAULT_NPTS_DOS)
        print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

        # 3. Define temperature and chemical potential ranges
        Tr = DEFAULT_TEMPERATURES
        # Chemical potential range around Fermi level
        mu_margin = 0.1  # Ha
        mu_min = data.fermi - mu_margin
        mu_max = data.fermi + mu_margin
        mur = np.linspace(mu_min, mu_max, 200)
        print(f"    Tr: {Tr}, mur: {len(mur)} points")

        # 4. Compute Fermi integrals
        print("    fermiintegrals...")
        N, L0, L1, L2, *_ = bandlib.fermiintegrals(
            epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
        )
        print(f"    L0: {L0.shape}")

        # 5. Compute Onsager coefficients
        print("    calc_Onsager_coefficients...")
        sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
            L0, L1, L2, mur, Tr, vuc
        )
        print(f"    sigma: {sigma.shape}")

        # 6. Solve for intrinsic mu at each temperature
        print("    solve_for_mu...")
        mu0 = np.zeros(len(Tr))
        for iT, T in enumerate(Tr):
            mu0[iT] = bandlib.solve_for_mu(
                epsilon, dos, data.nelect, T, data.dosweight,
                refine=True, try_center=True
            )
            print(f"      T={T}K: mu0={mu0[iT]:.6f} Ha")

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
            # Onsager coefficients
            'sigma': sigma,
            'S': S,
            'kappa': kappa,
            # solve_for_mu
            'mu0': mu0,
        })

    # Save reference data
    output_path = os.path.join(OUTPUT_DIR, 'wien2k_si_end2end.npz')
    np.savez(output_path, **save_dict)

    print(f"  Saved: {output_path}")
    return output_path


def generate_abinit_si_interpolate(include_transport=True):
    """Generate reference for ABINIT Si interpolation (end-to-end).

    Args:
        include_transport: If True, also compute transport coefficients.
    """
    print("\nGenerating ABINIT Si interpolation reference...")

    # Load DFT data from ABINIT output
    data_path = get_material_path("Si.abinit")
    data = dft.DFTData(str(data_path))

    print(f"  Loaded: {data_path}")
    print(f"  K-points: {data.kpoints.shape}")
    print(f"  Bands: {data.ebands.shape}")
    print(f"  Fermi: {data.fermi} Ha")
    print(f"  dosweight: {data.dosweight}")
    print(f"  nelect: {data.nelect}")

    lattvec = data.atoms.get_cell().T
    vuc = np.abs(np.linalg.det(lattvec))  # Unit cell volume in Bohr^3
    print(f"  Unit cell volume: {vuc:.4f} Bohr^3")

    # Get equivalences with target ~5000 k-points
    nkpt_target = 5000
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt_target)

    print(f"  Equivalences: {len(equivalences)}")

    # Fit coefficients (using all bands, no energy filter for simplicity)
    coeffs = fite.fitde3D(data, equivalences)

    print(f"  Coefficients: {coeffs.shape}")

    # Test reconstruction at original k-points
    ebands_reconstructed, _ = fite.getBands(
        data.kpoints, equivalences, lattvec, coeffs
    )

    max_error = np.max(np.abs(ebands_reconstructed - data.ebands))
    print(f"  Reconstruction error: {max_error:.2e}")

    # Build save dictionary
    save_dict = {
        # Input data
        'lattvec': lattvec,
        'kpoints': data.kpoints,
        'ebands': data.ebands,
        'fermi': data.fermi,
        'nelect': data.nelect,
        'dosweight': data.dosweight,
        'vuc': vuc,
        # Equivalence info
        'n_equivalences': len(equivalences),
        'equiv_sizes': np.array([len(eq) for eq in equivalences]),
        # Representative points from each equivalence class (first point)
        'equiv_reps': np.array([eq[0] for eq in equivalences]),
        # Output coefficients
        'coeffs_real': coeffs.real,
        'coeffs_imag': coeffs.imag,
        # Reconstruction test
        'ebands_reconstructed': ebands_reconstructed,
        'max_reconstruction_error': max_error,
    }

    # Transport calculation (E2E integrate)
    if include_transport:
        print("\n  Computing transport (E2E integrate)...")

        # 1. Reconstruct bands via FFT (getBTPbands)
        print("    getBTPbands...")
        eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
        print(f"    eband: {eband.shape}, vvband: {vvband.shape}")

        # 2. Compute DOS and transport DOS (BTPDOS)
        print("    BTPDOS...")
        epsilon, dos, vvdos, *_ = bandlib.BTPDOS(eband, vvband, npts=DEFAULT_NPTS_DOS)
        print(f"    epsilon: {epsilon.shape}, dos: {dos.shape}")

        # 3. Define temperature and chemical potential ranges
        Tr = DEFAULT_TEMPERATURES
        # Chemical potential range around Fermi level
        mu_margin = 0.1  # Ha
        mu_min = data.fermi - mu_margin
        mu_max = data.fermi + mu_margin
        mur = np.linspace(mu_min, mu_max, 200)
        print(f"    Tr: {Tr}, mur: {len(mur)} points")

        # 4. Compute Fermi integrals
        print("    fermiintegrals...")
        N, L0, L1, L2, *_ = bandlib.fermiintegrals(
            epsilon, dos, vvdos, mur, Tr, dosweight=data.dosweight
        )
        print(f"    L0: {L0.shape}")

        # 5. Compute Onsager coefficients
        print("    calc_Onsager_coefficients...")
        sigma, S, kappa, *_ = bandlib.calc_Onsager_coefficients(
            L0, L1, L2, mur, Tr, vuc
        )
        print(f"    sigma: {sigma.shape}")

        # 6. Solve for intrinsic mu at each temperature
        print("    solve_for_mu...")
        mu0 = np.zeros(len(Tr))
        for iT, T in enumerate(Tr):
            mu0[iT] = bandlib.solve_for_mu(
                epsilon, dos, data.nelect, T, data.dosweight,
                refine=True, try_center=True
            )
            print(f"      T={T}K: mu0={mu0[iT]:.6f} Ha")

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
            # Onsager coefficients
            'sigma': sigma,
            'S': S,
            'kappa': kappa,
            # solve_for_mu
            'mu0': mu0,
        })

    # Save reference data
    output_path = os.path.join(OUTPUT_DIR, 'abinit_si_end2end.npz')
    np.savez(output_path, **save_dict)

    print(f"  Saved: {output_path}")
    return output_path


if __name__ == '__main__':

    generate_si_interpolate()
    generate_pbte_interpolate()
    generate_qe_si_interpolate()
    generate_wien2k_si_interpolate()
    generate_abinit_si_interpolate()

    print("\nDone!")
