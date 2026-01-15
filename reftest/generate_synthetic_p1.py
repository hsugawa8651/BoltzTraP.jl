#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Generate synthetic triclinic P1 test data in GENE format.

This script creates a synthetic low-symmetry (P1) structure with
fake band energies based on a free-electron model for testing
BoltzTraP.jl loaders with triclinic cells.
"""

import os
import numpy as np

from common import get_boltztrap2_data_dir, OUTPUT_DIR


def triclinic_cell():
    """Create a triclinic P1 cell (all angles different from 90 degrees).

    Returns lattice vectors in Bohr (column convention after transpose).
    """
    # Triclinic cell parameters (in Bohr)
    a, b, c = 5.1, 4.8, 6.2
    alpha, beta, gamma = 85.0, 80.0, 75.0  # degrees

    # Convert to radians
    alpha_r = np.radians(alpha)
    beta_r = np.radians(beta)
    gamma_r = np.radians(gamma)

    # Build lattice vectors (standard convention)
    # a1 along x
    a1 = np.array([a, 0.0, 0.0])

    # a2 in xy plane
    a2 = np.array([b * np.cos(gamma_r), b * np.sin(gamma_r), 0.0])

    # a3 general
    c_x = c * np.cos(beta_r)
    c_y = c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r)
    c_z = np.sqrt(c**2 - c_x**2 - c_y**2)
    a3 = np.array([c_x, c_y, c_z])

    return np.array([a1, a2, a3])  # (3, 3) row convention


def monkhorst_pack(n1, n2, n3):
    """Generate Monkhorst-Pack k-point grid.

    Returns k-points in fractional coordinates.
    """
    kpoints = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                kx = (2 * i - n1 + 1) / (2 * n1)
                ky = (2 * j - n2 + 1) / (2 * n2)
                kz = (2 * k - n3 + 1) / (2 * n3)
                kpoints.append([kx, ky, kz])
    return np.array(kpoints)


def free_electron_bands(kpoints, lattvec, nbands=8, offset=0.0):
    """Generate free-electron-like band energies.

    E(k) = hbar^2 |k + G|^2 / (2m)

    In atomic units: E = |k + G|^2 / 2

    Parameters
    ----------
    kpoints : array (nk, 3)
        K-points in fractional coordinates
    lattvec : array (3, 3)
        Lattice vectors (row convention)
    nbands : int
        Number of bands to generate
    offset : float
        Energy offset (e.g., for simulating band gap)

    Returns
    -------
    ebands : array (nbands, nk)
        Band energies in Hartree
    """
    # Reciprocal lattice vectors
    recip = 2 * np.pi * np.linalg.inv(lattvec).T

    nk = len(kpoints)
    ebands = np.zeros((nbands, nk))

    # G-vectors to consider (sorted by |G|^2)
    G_list = []
    for h in range(-2, 3):
        for k in range(-2, 3):
            for l in range(-2, 3):
                G = h * recip[0] + k * recip[1] + l * recip[2]
                G_list.append((np.dot(G, G), G))

    G_list.sort(key=lambda x: x[0])
    G_list = G_list[:nbands * 2]  # Take enough G-vectors

    for ik, kfrac in enumerate(kpoints):
        # k in Cartesian
        k_cart = kfrac[0] * recip[0] + kfrac[1] * recip[1] + kfrac[2] * recip[2]

        # Calculate |k + G|^2 / 2 for each G
        energies = []
        for _, G in G_list:
            k_plus_G = k_cart + G
            E = 0.5 * np.dot(k_plus_G, k_plus_G)
            energies.append(E)

        # Sort and take lowest nbands
        energies.sort()
        ebands[:, ik] = energies[:nbands]

    # Add offset to simulate valence/conduction bands
    ebands += offset

    return ebands


def write_gene_structure(filename, lattvec, positions, species):
    """Write GENE .structure file.

    Parameters
    ----------
    filename : str
        Output file path
    lattvec : array (3, 3)
        Lattice vectors in Bohr (row convention)
    positions : array (natom, 3)
        Cartesian positions in Bohr
    species : list of str
        Element symbols
    """
    with open(filename, 'w') as f:
        f.write("Synthetic triclinic P1 structure for BoltzTraP.jl testing\n")
        for i in range(3):
            f.write(f"  {lattvec[i, 0]:18.12e}  {lattvec[i, 1]:18.12e}  {lattvec[i, 2]:18.12e}\n")
        f.write(f"  {len(species)}\n")
        for s, pos in zip(species, positions):
            f.write(f"{s} {pos[0]:.12e} {pos[1]:.12e} {pos[2]:.12e}\n")


def write_gene_energy(filename, kpoints, ebands, fermi, nspin=1):
    """Write GENE .energy file.

    Parameters
    ----------
    filename : str
        Output file path
    kpoints : array (nk, 3)
        K-points in fractional coordinates
    ebands : array (nbands, nk)
        Band energies in Hartree
    fermi : float
        Fermi energy in Hartree (will be converted to Rydberg)
    nspin : int
        Number of spin channels
    """
    nk = len(kpoints)
    nbands = ebands.shape[0]

    # Convert Fermi energy to Rydberg
    fermi_ry = fermi * 2.0

    with open(filename, 'w') as f:
        f.write("Synthetic triclinic P1 bands for BoltzTraP.jl testing\n")
        f.write(f"  {nk}  {nspin}  {fermi_ry:18.12e}\n")

        for ik in range(nk):
            kpt = kpoints[ik]
            f.write(f"  {kpt[0]:18.12e}  {kpt[1]:18.12e}  {kpt[2]:18.12e}  {nbands}\n")

            for ib in range(nbands):
                # Energy in Rydberg, velocities set to zero
                E_ry = ebands[ib, ik] * 2.0
                f.write(f"  {E_ry:18.12e}  0.0  0.0  0.0\n")


def generate_triclinic_p1():
    """Generate synthetic triclinic P1 test data."""
    print("Generating synthetic triclinic P1 data...")

    # Create output directory in BoltzTraP2 data
    gene_output_dir = get_boltztrap2_data_dir() / "Triclinic.GENE.synthetic"
    os.makedirs(gene_output_dir, exist_ok=True)

    # Triclinic cell
    lattvec = triclinic_cell()
    print(f"  Lattice vectors (Bohr):\n{lattvec}")

    # Check angles
    def angle(v1, v2):
        cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(np.clip(cos_a, -1, 1)))

    alpha = angle(lattvec[1], lattvec[2])
    beta = angle(lattvec[0], lattvec[2])
    gamma = angle(lattvec[0], lattvec[1])
    print(f"  Angles: α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")

    # Two atoms at general positions (to break all symmetry)
    frac_positions = np.array([
        [0.12, 0.23, 0.34],
        [0.67, 0.78, 0.56]
    ])

    # Convert to Cartesian
    cart_positions = frac_positions @ lattvec
    species = ["C", "N"]  # Different elements to break inversion

    print(f"  Atoms: {species}")
    print(f"  Fractional positions:\n{frac_positions}")

    # K-points (4x4x4 grid = 64 points)
    kpoints = monkhorst_pack(4, 4, 4)
    print(f"  K-points: {len(kpoints)}")

    # Free-electron bands (8 bands)
    nbands = 8
    ebands = free_electron_bands(kpoints, lattvec, nbands=nbands, offset=-0.5)

    # Set Fermi level at band 4 (4 valence electrons)
    fermi = np.mean([ebands[3, :].max(), ebands[4, :].min()])
    nelect = 4.0

    print(f"  Bands: {nbands}")
    print(f"  Fermi energy: {fermi:.6f} Ha")
    print(f"  Band gap estimate: {ebands[4, :].min() - ebands[3, :].max():.4f} Ha")

    # Write files
    struct_file = gene_output_dir / "Triclinic.structure"
    energy_file = gene_output_dir / "Triclinic.energy"

    write_gene_structure(struct_file, lattvec, cart_positions, species)
    write_gene_energy(energy_file, kpoints, ebands, fermi)

    print(f"  Written: {struct_file}")
    print(f"  Written: {energy_file}")

    # Also save NPZ reference for Julia comparison
    npz_file = OUTPUT_DIR / 'gene_triclinic_p1.npz'

    np.savez(npz_file,
        lattvec=lattvec.T,  # Column convention
        positions=frac_positions,
        symbols=np.array(list(",".join(species).encode('utf-8')), dtype=np.uint8),
        kpoints=kpoints,
        ebands=ebands,
        fermi=np.array(fermi),
        nelect=np.array(nelect),
        dosweight=np.array(2.0),
        nspin=np.array(1),
        source_format=np.array(list("GENE".encode('utf-8')), dtype=np.uint8),
    )
    print(f"  Written: {npz_file}")

    return gene_output_dir


if __name__ == '__main__':
    generate_triclinic_p1()
    print("\nDone!")
