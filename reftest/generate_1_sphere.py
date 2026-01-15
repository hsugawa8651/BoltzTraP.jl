#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Generate reference data from Python BoltzTraP2.

Run this script to create .npz files containing reference outputs
for comparison with Julia implementation.
"""

import numpy as np
import ase
from pathlib import Path

import BoltzTraP2.sphere as sphere

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)


def si_diamond():
    """Si diamond structure test case (48 rotations, Oh symmetry)."""
    # FCC lattice for Si
    a = 5.45052526  # Angstrom
    lattvec = a * 0.5 * (np.ones((3, 3)) - np.eye(3))

    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,  # ASE uses row vectors, but get_cell().T gives column
        scaled_positions=np.array([
            [0.125, 0.125, 0.125],
            [0.875, 0.875, 0.875]
        ]),
        pbc=True,
    )

    # Parameters
    radius = 127.8
    symprec = 1e-5

    # Compute bounds
    bounds = sphere.compute_bounds(lattvec, radius)

    # Compute rotations
    nrotations = sphere.calc_nrotations(atoms, None, symprec)

    # Compute tensor basis
    tensor_basis = sphere.calc_tensor_basis(atoms, None, symprec)

    # Compute equivalences
    equivalences = sphere.calc_sphere_quotient_set(
        atoms, None, radius, bounds, symprec
    )

    # Save to npz
    np.savez(
        DATA_DIR / "si_diamond.npz",
        # Input data
        lattvec=atoms.get_cell().T,  # Column vectors
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        radius=radius,
        symprec=symprec,
        # Output: sphere module
        bounds=bounds,
        nrotations=nrotations,
        tensor_basis=tensor_basis,
        n_equivalences=len(equivalences),
        # Equivalences saved as separate arrays (variable length)
        **{f"equiv_{i}": eq for i, eq in enumerate(equivalences)}
    )
    print(f"Generated: si_diamond.npz ({nrotations} rotations, {len(equivalences)} equiv classes)")


def simple_cubic():
    """Simple cubic single atom (48 rotations, Oh symmetry)."""
    a = 5.0  # Angstrom
    lattvec = np.eye(3) * a

    atoms = ase.Atoms(
        "X",
        cell=lattvec.T,
        scaled_positions=[[0, 0, 0]],
        pbc=True,
    )

    radius = 50.0
    symprec = 1e-5
    bounds = sphere.compute_bounds(lattvec, radius)
    nrotations = sphere.calc_nrotations(atoms, None, symprec)
    tensor_basis = sphere.calc_tensor_basis(atoms, None, symprec)

    np.savez(
        DATA_DIR / "simple_cubic.npz",
        lattvec=atoms.get_cell().T,
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        radius=radius,
        symprec=symprec,
        bounds=bounds,
        nrotations=nrotations,
        tensor_basis=tensor_basis,
    )
    print(f"Generated: simple_cubic.npz ({nrotations} rotations)")


def fe_bcc_magnetic():
    """Fe BCC with collinear magnetic moments (reduced symmetry)."""
    a = 2.87  # Angstrom
    # BCC lattice vectors
    lattvec = a * 0.5 * np.array([
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1]
    ], dtype=float).T  # Column vectors

    atoms = ase.Atoms(
        "Fe",
        cell=lattvec.T,
        scaled_positions=[[0, 0, 0]],
        pbc=True,
    )

    # Magnetic moment (collinear)
    magmom = np.array([2.2])

    radius = 50.0
    symprec = 1e-5
    bounds = sphere.compute_bounds(lattvec, radius)
    nrotations = sphere.calc_nrotations(atoms, magmom, symprec)
    tensor_basis = sphere.calc_tensor_basis(atoms, magmom, symprec)

    np.savez(
        DATA_DIR / "fe_bcc_magnetic.npz",
        lattvec=atoms.get_cell().T,
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        magmom=magmom,
        radius=radius,
        symprec=symprec,
        bounds=bounds,
        nrotations=nrotations,
        tensor_basis=tensor_basis,
    )
    print(f"Generated: fe_bcc_magnetic.npz ({nrotations} rotations)")


def unit_cube():
    """Unit cube (a=1) for simple boundary testing."""
    lattvec = np.eye(3)

    atoms = ase.Atoms(
        "X",
        cell=lattvec.T,
        scaled_positions=[[0, 0, 0]],
        pbc=True,
    )

    # Test various radii
    radii = [1.0, 2.5, 5.0, 10.0]
    symprec = 1e-5

    results = {}
    for r in radii:
        bounds = sphere.compute_bounds(lattvec, r)
        results[f"bounds_r{r}"] = bounds

    nrotations = sphere.calc_nrotations(atoms, None, symprec)

    np.savez(
        DATA_DIR / "unit_cube.npz",
        lattvec=lattvec,
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        radii=radii,
        symprec=symprec,
        nrotations=nrotations,
        **results
    )
    print(f"Generated: unit_cube.npz ({nrotations} rotations)")


def triclinic_p1():
    """Triclinic P1 structure (2 rotations, minimal symmetry).

    Note: 2 rotations = identity + time-reversal (k -> -k) for non-magnetic.
    This is the minimum for band structure calculations.
    """
    # Triclinic lattice (all angles != 90 degrees)
    lattvec = np.array([
        [5.1, 0.0, 0.0],
        [0.7, 4.8, 0.0],
        [0.3, 0.4, 6.2]
    ])

    # Two different atoms at general positions
    atoms = ase.Atoms(
        ["C", "N"],
        cell=lattvec.T,
        scaled_positions=np.array([
            [0.123, 0.234, 0.345],
            [0.567, 0.678, 0.789],
        ]),
        pbc=True,
    )

    radius = 50.0
    symprec = 1e-5
    bounds = sphere.compute_bounds(lattvec, radius)
    nrotations = sphere.calc_nrotations(atoms, None, symprec)
    tensor_basis = sphere.calc_tensor_basis(atoms, None, symprec)

    # Compute equivalences
    equivalences = sphere.calc_sphere_quotient_set(
        atoms, None, radius, bounds, symprec
    )

    np.savez(
        DATA_DIR / "triclinic_p1.npz",
        lattvec=atoms.get_cell().T,
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        radius=radius,
        symprec=symprec,
        bounds=bounds,
        nrotations=nrotations,
        tensor_basis=tensor_basis,
        n_equivalences=len(equivalences),
        **{f"equiv_{i}": eq for i, eq in enumerate(equivalences)}
    )
    print(f"Generated: triclinic_p1.npz ({nrotations} rotations, {len(equivalences)} equiv classes)")


def monoclinic():
    """Monoclinic structure (4 rotations, low symmetry)."""
    # Monoclinic lattice (beta != 90 degrees)
    lattvec = np.array([
        [5.0, 0.0, 0.0],
        [0.0, 6.0, 0.0],
        [1.5, 0.0, 7.0]
    ])

    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,
        scaled_positions=np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ]),
        pbc=True,
    )

    radius = 50.0
    symprec = 1e-5
    bounds = sphere.compute_bounds(lattvec, radius)
    nrotations = sphere.calc_nrotations(atoms, None, symprec)
    tensor_basis = sphere.calc_tensor_basis(atoms, None, symprec)

    # Compute equivalences
    equivalences = sphere.calc_sphere_quotient_set(
        atoms, None, radius, bounds, symprec
    )

    np.savez(
        DATA_DIR / "monoclinic.npz",
        lattvec=atoms.get_cell().T,
        positions=atoms.get_scaled_positions(),
        types=atoms.numbers,
        radius=radius,
        symprec=symprec,
        bounds=bounds,
        nrotations=nrotations,
        tensor_basis=tensor_basis,
        n_equivalences=len(equivalences),
        **{f"equiv_{i}": eq for i, eq in enumerate(equivalences)}
    )
    print(f"Generated: monoclinic.npz ({nrotations} rotations, {len(equivalences)} equiv classes)")


if __name__ == "__main__":
    print("Generating reference data from Python BoltzTraP2...")
    print()

    # High symmetry (48 rotations)
    unit_cube()
    simple_cubic()
    si_diamond()

    # Magnetic (reduced by spin)
    fe_bcc_magnetic()

    # Low symmetry
    monoclinic()     # 4 rotations
    triclinic_p1()   # 2 rotations (minimal)

    print()
    print("Done. Reference data saved to:", DATA_DIR)
