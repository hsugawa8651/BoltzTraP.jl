#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""Generate synthetic low-symmetry test data for all DFT formats.

Creates monoclinic and triclinic test structures for:
- GENE (generic format)
- VASP (vasprun.xml)
- Wien2k (case.struct, case.energy, case.scf)
- QE (XML output)
- ABINIT (NetCDF GSR.nc)

These synthetic structures use free-electron band models for testing
BoltzTraP.jl's ability to handle low-symmetry crystals.
"""

import os
import numpy as np
from xml.etree import ElementTree as ET
from xml.dom import minidom

from common import get_boltztrap2_data_dir, OUTPUT_DIR

# Constants
BOHR_TO_ANG = 0.529177210903
ANG_TO_BOHR = 1.0 / BOHR_TO_ANG
HA_TO_RY = 2.0
RY_TO_HA = 0.5
HA_TO_EV = 27.211386245988


# =============================================================================
# Structure Definitions
# =============================================================================

def monoclinic_cell():
    """Create monoclinic cell (β ≠ 90°, others = 90°).

    Returns lattice vectors in Bohr (row convention).
    """
    a, b, c = 5.0, 6.0, 7.0  # Bohr
    beta = 100.0  # degrees (only β ≠ 90°)

    beta_r = np.radians(beta)

    # Standard monoclinic: a along x, b along y, c in xz plane
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([0.0, b, 0.0])
    a3 = np.array([c * np.cos(beta_r), 0.0, c * np.sin(beta_r)])

    return np.array([a1, a2, a3])


def triclinic_cell():
    """Create triclinic cell (all angles ≠ 90°).

    Returns lattice vectors in Bohr (row convention).
    """
    a, b, c = 5.1, 4.8, 6.2  # Bohr
    alpha, beta, gamma = 85.0, 80.0, 75.0  # degrees

    alpha_r = np.radians(alpha)
    beta_r = np.radians(beta)
    gamma_r = np.radians(gamma)

    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([b * np.cos(gamma_r), b * np.sin(gamma_r), 0.0])

    c_x = c * np.cos(beta_r)
    c_y = c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r)
    c_z = np.sqrt(c**2 - c_x**2 - c_y**2)
    a3 = np.array([c_x, c_y, c_z])

    return np.array([a1, a2, a3])


def get_atoms_general_position():
    """Return atoms at general positions (breaks all symmetry).

    Returns (fractional_positions, species, atomic_numbers)
    """
    frac_pos = np.array([
        [0.12, 0.23, 0.34],
        [0.67, 0.78, 0.56]
    ])
    species = ["C", "N"]
    atomic_numbers = [6, 7]
    return frac_pos, species, atomic_numbers


def monkhorst_pack(n1, n2, n3):
    """Generate Monkhorst-Pack k-point grid."""
    kpoints = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                kx = (2 * i - n1 + 1) / (2 * n1)
                ky = (2 * j - n2 + 1) / (2 * n2)
                kz = (2 * k - n3 + 1) / (2 * n3)
                kpoints.append([kx, ky, kz])
    return np.array(kpoints)


def free_electron_bands(kpoints, lattvec, nbands=8, offset=-0.3):
    """Generate free-electron-like band energies.

    E(k) = |k + G|^2 / 2 (atomic units)
    """
    recip = 2 * np.pi * np.linalg.inv(lattvec).T
    nk = len(kpoints)

    G_list = []
    for h in range(-2, 3):
        for l in range(-2, 3):
            for m in range(-2, 3):
                G = h * recip[0] + l * recip[1] + m * recip[2]
                G_list.append((np.dot(G, G), G))

    G_list.sort(key=lambda x: x[0])

    ebands = np.zeros((nbands, nk))
    for ik, kfrac in enumerate(kpoints):
        k_cart = kfrac @ recip
        energies = sorted([0.5 * np.dot(k_cart + G, k_cart + G) for _, G in G_list[:nbands*3]])
        ebands[:, ik] = energies[:nbands]

    return ebands + offset


# =============================================================================
# GENE Format Writer
# =============================================================================

def write_gene(output_dir, lattvec, frac_pos, species, kpoints, ebands, fermi):
    """Write GENE format files."""
    os.makedirs(output_dir, exist_ok=True)

    # case.structure
    case_name = os.path.basename(output_dir).split('.')[0]
    struct_file = os.path.join(output_dir, f"{case_name}.structure")
    energy_file = os.path.join(output_dir, f"{case_name}.energy")

    cart_pos = frac_pos @ lattvec

    with open(struct_file, 'w') as f:
        f.write("Synthetic structure for BoltzTraP.jl testing\n")
        for i in range(3):
            f.write(f"  {lattvec[i, 0]:18.12e}  {lattvec[i, 1]:18.12e}  {lattvec[i, 2]:18.12e}\n")
        f.write(f"  {len(species)}\n")
        for s, pos in zip(species, cart_pos):
            f.write(f"{s} {pos[0]:.12e} {pos[1]:.12e} {pos[2]:.12e}\n")

    # case.energy
    nk, nbands = len(kpoints), ebands.shape[0]
    fermi_ry = fermi * HA_TO_RY

    with open(energy_file, 'w') as f:
        f.write("Synthetic bands for BoltzTraP.jl testing\n")
        f.write(f"  {nk}  1  {fermi_ry:18.12e}\n")
        for ik in range(nk):
            f.write(f"  {kpoints[ik, 0]:18.12e}  {kpoints[ik, 1]:18.12e}  {kpoints[ik, 2]:18.12e}  {nbands}\n")
            for ib in range(nbands):
                E_ry = ebands[ib, ik] * HA_TO_RY
                f.write(f"  {E_ry:18.12e}  0.0  0.0  0.0\n")

    print(f"  Written: {struct_file}")
    print(f"  Written: {energy_file}")


# =============================================================================
# Wien2k Format Writer
# =============================================================================

def write_wien2k(output_dir, lattvec, frac_pos, species, kpoints, ebands, fermi):
    """Write Wien2k format files (case.struct, case.energy, case.scf)."""
    os.makedirs(output_dir, exist_ok=True)

    case_name = os.path.basename(output_dir)
    struct_file = os.path.join(output_dir, f"{case_name}.struct")
    energy_file = os.path.join(output_dir, f"{case_name}.energy")
    scf_file = os.path.join(output_dir, f"{case_name}.scf")

    # Calculate cell parameters from lattice vectors
    a = np.linalg.norm(lattvec[0])
    b = np.linalg.norm(lattvec[1])
    c = np.linalg.norm(lattvec[2])

    cos_alpha = np.dot(lattvec[1], lattvec[2]) / (b * c)
    cos_beta = np.dot(lattvec[0], lattvec[2]) / (a * c)
    cos_gamma = np.dot(lattvec[0], lattvec[1]) / (a * b)

    alpha = np.degrees(np.arccos(np.clip(cos_alpha, -1, 1)))
    beta = np.degrees(np.arccos(np.clip(cos_beta, -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(cos_gamma, -1, 1)))

    # case.struct
    with open(struct_file, 'w') as f:
        f.write("Synthetic structure\n")
        f.write(f"P   LATTICE,NONEQUIV.ATOMS:{len(species):3d}\n")
        f.write("MODE OF CALC=RELA unit=bohr\n")
        f.write(f"{a:10.6f}{b:10.6f}{c:10.6f}{alpha:10.6f}{beta:10.6f}{gamma:10.6f}\n")

        for iatom, (s, pos) in enumerate(zip(species, frac_pos)):
            f.write(f"ATOM{iatom+1:4d}: X={pos[0]:10.8f} Y={pos[1]:10.8f} Z={pos[2]:10.8f}\n")
            f.write(f"          MULT= 1          ISPLIT= 8\n")
            f.write(f"{s:2s}         NPT=  781  R0=0.00010000 RMT=    2.0000   Z:{iatom+6:5.1f}\n")
            f.write(f"LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000\n")
            f.write(f"                     0.0000000 1.0000000 0.0000000\n")
            f.write(f"                     0.0000000 0.0000000 1.0000000\n")

        f.write("   1      NUMBER OF SYMMETRY OPERATIONS\n")
        f.write(" 1 0 0 0.00000000\n")
        f.write(" 0 1 0 0.00000000\n")
        f.write(" 0 0 1 0.00000000\n")
        f.write("       1\n")

    # case.energy
    nk, nbands = len(kpoints), ebands.shape[0]

    with open(energy_file, 'w') as f:
        for ik in range(nk):
            kpt = kpoints[ik]
            # Wien2k energy format: kx ky kz  kpt_idx  atomidx  nbands  weight
            # Fixed-width Fortran format: kx/ky/kz each 19 chars (incl. sign/space)
            # Real Wien2k format: columns [1:19], [20:38], [39:57] for k-coords
            # nband at [74:79], so total before nband = 73 chars
            f.write(f"{kpt[0]:19.12E}{kpt[1]:19.12E}{kpt[2]:19.12E}{ik+1:10d}{1:6d}{nbands:6d}  1.0   \n")
            for ib in range(nbands):
                E_ry = ebands[ib, ik] * HA_TO_RY
                f.write(f"{ib+1:12d}  {E_ry:19.12f}\n")

    # case.scf (minimal, just need Fermi energy)
    fermi_ry = fermi * HA_TO_RY
    with open(scf_file, 'w') as f:
        f.write(":VOL  : UNIT CELL VOLUME =     200.00000\n")
        f.write(f":FER  : F E R M I - ENERGY(TETRAH.M.)=   {fermi_ry:.10f}\n")

    print(f"  Written: {struct_file}")
    print(f"  Written: {energy_file}")
    print(f"  Written: {scf_file}")


# =============================================================================
# VASP Format Writer
# =============================================================================

def write_vasp(output_dir, lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi):
    """Write VASP format files (vasprun.xml, POSCAR)."""
    os.makedirs(output_dir, exist_ok=True)

    vasprun_file = os.path.join(output_dir, "vasprun.xml")
    poscar_file = os.path.join(output_dir, "POSCAR")

    nk, nbands = len(kpoints), ebands.shape[0]
    natoms = len(species)

    # Convert to Angstrom for VASP
    lattvec_ang = lattvec * BOHR_TO_ANG

    # Build vasprun.xml
    root = ET.Element("modeling")

    # Generator info
    gen = ET.SubElement(root, "generator")
    ET.SubElement(gen, "i", name="program").text = "synthetic"
    ET.SubElement(gen, "i", name="version").text = "1.0"

    # Parameters section (NELECT)
    parameters = ET.SubElement(root, "parameters")
    ET.SubElement(parameters, "i", name="NELECT").text = " 8.0"

    # Atom info
    atominfo = ET.SubElement(root, "atominfo")
    ET.SubElement(atominfo, "atoms").text = str(natoms)
    ET.SubElement(atominfo, "types").text = str(len(set(species)))

    # atoms array with proper structure
    array = ET.SubElement(atominfo, "array", name="atoms")
    dimension = ET.SubElement(array, "dimension", dim="1")
    dimension.text = "ion"
    field1 = ET.SubElement(array, "field", type="string")
    field1.text = "element"
    field2 = ET.SubElement(array, "field", type="int")
    field2.text = "atomtype"
    set_elem = ET.SubElement(array, "set")
    unique_species = sorted(set(species))
    for s in species:
        rc = ET.SubElement(set_elem, "rc")
        c1 = ET.SubElement(rc, "c")
        c1.text = s
        c2 = ET.SubElement(rc, "c")
        c2.text = str(unique_species.index(s) + 1)

    # atomtypes array
    array = ET.SubElement(atominfo, "array", name="atomtypes")
    dimension = ET.SubElement(array, "dimension", dim="1")
    dimension.text = "type"
    field1 = ET.SubElement(array, "field", type="int")
    field1.text = "atomspertype"
    field2 = ET.SubElement(array, "field", type="string")
    field2.text = "element"
    field3 = ET.SubElement(array, "field")
    field3.text = "mass"
    field4 = ET.SubElement(array, "field", type="string")
    field4.text = "pseudopotential"
    set_elem = ET.SubElement(array, "set")
    for s in unique_species:
        rc = ET.SubElement(set_elem, "rc")
        c1 = ET.SubElement(rc, "c")
        c1.text = str(species.count(s))
        c2 = ET.SubElement(rc, "c")
        c2.text = s
        c3 = ET.SubElement(rc, "c")
        c3.text = "12.0"
        c4 = ET.SubElement(rc, "c")
        c4.text = "PAW"

    # Structure
    structure = ET.SubElement(root, "structure", name="finalpos")
    crystal = ET.SubElement(structure, "crystal")

    varray = ET.SubElement(crystal, "varray", name="basis")
    for i in range(3):
        v = ET.SubElement(varray, "v")
        v.text = f" {lattvec_ang[i, 0]:16.12f} {lattvec_ang[i, 1]:16.12f} {lattvec_ang[i, 2]:16.12f}"

    volume = abs(np.linalg.det(lattvec_ang))
    ET.SubElement(crystal, "i", name="volume").text = f" {volume:.12f}"

    varray = ET.SubElement(structure, "varray", name="positions")
    for pos in frac_pos:
        v = ET.SubElement(varray, "v")
        v.text = f" {pos[0]:16.12f} {pos[1]:16.12f} {pos[2]:16.12f}"

    # Calculation
    calc = ET.SubElement(root, "calculation")

    # DOS - Fermi energy
    dos = ET.SubElement(calc, "dos")
    ET.SubElement(dos, "i", name="efermi").text = f" {fermi * HA_TO_EV:.8f}"

    # Eigenvalues
    eigenvalues = ET.SubElement(calc, "eigenvalues")
    array = ET.SubElement(eigenvalues, "array")
    ET.SubElement(array, "dimension", dim="1").text = "band"
    ET.SubElement(array, "dimension", dim="2").text = "kpoint"
    ET.SubElement(array, "dimension", dim="3").text = "spin"
    ET.SubElement(array, "field").text = "eigene"
    ET.SubElement(array, "field").text = "occ"

    set_elem = ET.SubElement(array, "set")
    set_spin = ET.SubElement(set_elem, "set", comment="spin 1")

    for ik in range(nk):
        set_kpt = ET.SubElement(set_spin, "set", comment=f"kpoint {ik+1}")
        for ib in range(nbands):
            r = ET.SubElement(set_kpt, "r")
            E_ev = ebands[ib, ik] * HA_TO_EV
            occ = 2.0 if E_ev < fermi * HA_TO_EV else 0.0
            r.text = f" {E_ev:16.8f} {occ:8.4f}"

    # K-points
    kpoints_elem = ET.SubElement(calc, "kpoints")
    varray = ET.SubElement(kpoints_elem, "varray", name="kpointlist")
    for kpt in kpoints:
        v = ET.SubElement(varray, "v")
        v.text = f" {kpt[0]:16.12f} {kpt[1]:16.12f} {kpt[2]:16.12f}"

    varray = ET.SubElement(kpoints_elem, "varray", name="weights")
    weight = 1.0 / nk
    for _ in range(nk):
        v = ET.SubElement(varray, "v")
        v.text = f" {weight:.12f}"

    # Write XML
    xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
    # Remove extra blank lines
    xml_str = '\n'.join([line for line in xml_str.split('\n') if line.strip()])

    with open(vasprun_file, 'w') as f:
        f.write(xml_str)

    # Write POSCAR
    with open(poscar_file, 'w') as f:
        f.write("Synthetic structure\n")
        f.write("1.0\n")
        for i in range(3):
            f.write(f"  {lattvec_ang[i, 0]:16.12f}  {lattvec_ang[i, 1]:16.12f}  {lattvec_ang[i, 2]:16.12f}\n")
        f.write("  ".join(sorted(set(species))) + "\n")
        f.write("  ".join([str(species.count(s)) for s in sorted(set(species))]) + "\n")
        f.write("Direct\n")
        for pos in frac_pos:
            f.write(f"  {pos[0]:16.12f}  {pos[1]:16.12f}  {pos[2]:16.12f}\n")

    print(f"  Written: {vasprun_file}")
    print(f"  Written: {poscar_file}")


# =============================================================================
# QE Format Writer
# =============================================================================

def write_qe(output_dir, lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi):
    """Write Quantum ESPRESSO XML output."""
    os.makedirs(output_dir, exist_ok=True)

    xml_file = os.path.join(output_dir, "synthetic.xml")

    nk, nbands = len(kpoints), ebands.shape[0]
    natoms = len(species)

    # Build QE XML
    root = ET.Element("qes:espresso", {
        "xmlns:qes": "http://www.quantum-espresso.org/ns/qes/qes-1.0",
        "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance"
    })

    # Output section
    output = ET.SubElement(root, "output")

    # Atomic structure
    atomic_structure = ET.SubElement(output, "atomic_structure", nat=str(natoms), alat="1.0")

    cell = ET.SubElement(atomic_structure, "cell")
    for i, name in enumerate(["a1", "a2", "a3"]):
        v = ET.SubElement(cell, name)
        v.text = f" {lattvec[i, 0]:18.12e} {lattvec[i, 1]:18.12e} {lattvec[i, 2]:18.12e}"

    atomic_positions = ET.SubElement(atomic_structure, "atomic_positions")
    for s, pos in zip(species, frac_pos):
        cart = pos @ lattvec
        atom = ET.SubElement(atomic_positions, "atom", name=s)
        atom.text = f" {cart[0]:18.12e} {cart[1]:18.12e} {cart[2]:18.12e}"

    # Band structure
    band_structure = ET.SubElement(output, "band_structure")
    ET.SubElement(band_structure, "nbnd").text = str(nbands)
    ET.SubElement(band_structure, "nks").text = str(nk)
    ET.SubElement(band_structure, "nelec").text = "8.0"
    ET.SubElement(band_structure, "fermi_energy").text = f" {fermi:18.12e}"

    # Each k-point gets its own <ks_energies> block (matching real QE format)
    for ik in range(nk):
        ks_energies = ET.SubElement(band_structure, "ks_energies")

        kpoint = ET.SubElement(ks_energies, "k_point", weight=f"{1.0/nk:.12f}")
        kpoint.text = f" {kpoints[ik, 0]:18.12e} {kpoints[ik, 1]:18.12e} {kpoints[ik, 2]:18.12e}"

        eigenvalues = ET.SubElement(ks_energies, "eigenvalues", size=str(nbands))
        eigenvalues.text = " ".join([f"{ebands[ib, ik]:18.12e}" for ib in range(nbands)])

        occupations = ET.SubElement(ks_energies, "occupations", size=str(nbands))
        occupations.text = " ".join(["2.0" if ebands[ib, ik] < fermi else "0.0" for ib in range(nbands)])

    # Write XML
    xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
    xml_str = '\n'.join([line for line in xml_str.split('\n') if line.strip()])

    with open(xml_file, 'w') as f:
        f.write(xml_str)

    print(f"  Written: {xml_file}")


# =============================================================================
# ABINIT Format Writer (NetCDF)
# =============================================================================

def write_abinit(output_dir, lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi):
    """Write ABINIT NetCDF GSR file."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        print("  Skipped ABINIT: netCDF4 not installed (pip install netCDF4)")
        return

    os.makedirs(output_dir, exist_ok=True)

    nc_file = os.path.join(output_dir, "synthetic_GSR.nc")

    nk, nbands = len(kpoints), ebands.shape[0]
    natoms = len(species)
    unique_species = list(sorted(set(species)))
    ntypat = len(unique_species)

    with Dataset(nc_file, 'w', format='NETCDF4') as ds:
        # Dimensions - matching actual ABINIT GSR.nc structure
        ds.createDimension('three', 3)
        ds.createDimension('number_of_atoms', natoms)
        ds.createDimension('number_of_atom_species', ntypat)
        ds.createDimension('number_of_kpoints', nk)
        ds.createDimension('max_number_of_states', nbands)
        ds.createDimension('number_of_spinor_components', 1)
        ds.createDimension('number_of_spins', 1)
        ds.createDimension('symbol_length', 80)

        # Lattice vectors (3, 3) in Bohr - rows are vectors
        pv = ds.createVariable('primitive_vectors', 'f8', ('three', 'three'))
        pv[:] = lattvec

        # Atomic positions (natom, 3) fractional - matching original ABINIT format
        pos = ds.createVariable('reduced_atom_positions', 'f8', ('number_of_atoms', 'three'))
        pos[:] = frac_pos

        # Species index (1-based)
        sp_idx = ds.createVariable('atom_species', 'i4', ('number_of_atoms',))
        sp_idx[:] = np.array([unique_species.index(s) + 1 for s in species], dtype=np.int32)

        # Species names (ntypat, symbol_length) - character array matching original format
        sp_names = ds.createVariable('atom_species_names', 'S1', ('number_of_atom_species', 'symbol_length'))
        for i, s in enumerate(unique_species):
            # Pad with spaces
            name_chars = list(s.ljust(80))
            sp_names[i, :] = np.array(name_chars, dtype='S1')

        # K-points (nk, 3) fractional - matching original ABINIT format
        kpts = ds.createVariable('reduced_coordinates_of_kpoints', 'f8', ('number_of_kpoints', 'three'))
        kpts[:] = kpoints

        # Eigenvalues (nspin, nk, nbands) in Hartree - matching original ABINIT format
        eig = ds.createVariable('eigenvalues', 'f8', ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'))
        eig[0, :, :] = ebands.T  # Transpose to (nk, nbands)

        # Fermi energy in Hartree (scalar)
        fermie = ds.createVariable('fermie', 'f8')
        fermie[...] = fermi

        # Number of electrons (scalar)
        nelect = ds.createVariable('nelect', 'f8')
        nelect[...] = 8.0

    print(f"  Written: {nc_file}")


# =============================================================================
# Main Generator
# =============================================================================

def generate_all_formats(cell_func, name_suffix):
    """Generate test data in all formats for a given cell type."""
    print(f"\n=== Generating {name_suffix} data for all formats ===")

    lattvec = cell_func()
    frac_pos, species, atomic_nums = get_atoms_general_position()
    kpoints = monkhorst_pack(4, 4, 4)
    ebands = free_electron_bands(kpoints, lattvec, nbands=8)
    fermi = np.mean([ebands[3, :].max(), ebands[4, :].min()])

    # Print cell info
    a = np.linalg.norm(lattvec[0])
    b = np.linalg.norm(lattvec[1])
    c = np.linalg.norm(lattvec[2])
    cos_alpha = np.dot(lattvec[1], lattvec[2]) / (b * c)
    cos_beta = np.dot(lattvec[0], lattvec[2]) / (a * c)
    cos_gamma = np.dot(lattvec[0], lattvec[1]) / (a * b)
    alpha = np.degrees(np.arccos(np.clip(cos_alpha, -1, 1)))
    beta = np.degrees(np.arccos(np.clip(cos_beta, -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(cos_gamma, -1, 1)))

    print(f"  Cell: a={a:.2f}, b={b:.2f}, c={c:.2f} Bohr")
    print(f"  Angles: α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")
    print(f"  K-points: {len(kpoints)}, Bands: {ebands.shape[0]}")
    print(f"  Fermi: {fermi:.6f} Ha")

    # Generate each format
    base_dir = get_boltztrap2_data_dir()

    print("\n  GENE format:")
    write_gene(base_dir / f"Synthetic.GENE.{name_suffix}",
               lattvec, frac_pos, species, kpoints, ebands, fermi)

    print("\n  Wien2k format:")
    write_wien2k(base_dir / f"Synthetic.W2K.{name_suffix}",
                 lattvec, frac_pos, species, kpoints, ebands, fermi)

    print("\n  VASP format:")
    write_vasp(base_dir / f"Synthetic.vasp.{name_suffix}",
               lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi)

    print("\n  QE format:")
    write_qe(base_dir / f"Synthetic.ESPRESSO.{name_suffix}" / "out",
             lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi)

    print("\n  ABINIT format:")
    write_abinit(base_dir / f"Synthetic.abinit.{name_suffix}",
                 lattvec, frac_pos, species, atomic_nums, kpoints, ebands, fermi)

    # Save NPZ reference
    npz_file = OUTPUT_DIR / f"synthetic_{name_suffix}.npz"
    np.savez(npz_file,
        lattvec=lattvec.T,  # Column convention
        positions=frac_pos,
        symbols=np.array(list(",".join(species).encode('utf-8')), dtype=np.uint8),
        kpoints=kpoints,
        ebands=ebands,
        fermi=np.array(fermi),
        nelect=np.array(8.0),
        dosweight=np.array(2.0),
        nspin=np.array(1),
    )
    print(f"\n  Reference: {npz_file}")


def main():
    """Generate all synthetic test data."""
    # Monoclinic
    generate_all_formats(monoclinic_cell, "monoclinic")

    # Triclinic
    generate_all_formats(triclinic_cell, "triclinic")

    print("\n" + "=" * 60)
    print("Done! Generated synthetic data for all formats.")
    print("=" * 60)


if __name__ == '__main__':
    main()
