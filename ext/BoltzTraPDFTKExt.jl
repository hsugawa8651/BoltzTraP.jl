# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
DFTK.jl extension for BoltzTraP.jl.

Provide `load_dftk` function to extract band structure data from DFTK SCF results.
No unit conversion is needed since DFTK uses atomic units (Hartree, Bohr) internally.
"""
module BoltzTraPDFTKExt

using BoltzTraP
using DFTK

"""
    load_dftk(scfres) -> DFTData{1} or DFTData{2}

Extract band structure data from DFTK SCF result.

DFTK uses atomic units internally, so no unit conversion is needed.
Returns `DFTData{1}` for non-magnetic or `DFTData{2}` for collinear magnetic calculations.

# Arguments

  - `scfres`: SCF result from DFTK `self_consistent_field`

# Returns

`DFTData{N}` with fields:

  - `lattice`: Lattice vectors (3×3) in Bohr (columns are vectors)
  - `positions`: Atomic positions (3×natoms) in fractional coordinates
  - `species`: Atom species names
  - `kpoints`: K-points (3×nkpts) in fractional coordinates
  - `weights`: K-point weights (nkpts,)
  - `ebands`: Band energies in Hartree
      + Non-magnetic: (nbands, nkpts, 1)
      + Collinear: (nbands, nkpts, 2) where [:,:,1]=up, [:,:,2]=down
  - `occupations`: Occupations (same shape as ebands)
  - `fermi`: Fermi energy in Hartree
  - `nelect`: Number of electrons
  - `magmom`: Magnetic moments per atom (collinear only)

# Example (Non-magnetic)

```julia
using DFTK
using BoltzTraP

model = model_LDA(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut = 30, kgrid = [8, 8, 8])
scfres = self_consistent_field(basis)

data = load_dftk(scfres)
@assert data isa DFTData{1}
```

# Example (Collinear magnetic)

```julia
using DFTK
using BoltzTraP

magnetic_moments = [4.0]  # Initial moment for Fe
model = model_LDA(lattice, atoms, positions; magnetic_moments)
basis = PlaneWaveBasis(model; Ecut = 30, kgrid = [8, 8, 8])
ρ0 = guess_density(basis, magnetic_moments)
scfres = self_consistent_field(basis; ρ = ρ0)

data = load_dftk(scfres)
@assert data isa DFTData{2}
@assert BoltzTraP.is_magnetic(data)
```

# Notes

  - Dense k-grid is required for accurate interpolation (10×10×10 or more recommended)
  - Non-collinear calculations are not supported    # Check spin polarization type
"""
function BoltzTraP.load_dftk(scfres)
    basis = scfres.basis
    model = basis.model

    # Check spin polarization type
    spin_pol = model.spin_polarization

    if spin_pol == :none
        return _load_dftk_nonmagnetic(scfres)
    elseif spin_pol == :collinear
        return _load_dftk_collinear(scfres)
    else
        error(
            "Unsupported spin_polarization: $spin_pol. " *
            "Only :none and :collinear are supported.",
        )
    end
end

#=
    _load_dftk_nonmagnetic(scfres) -> DFTData{1}

Load non-magnetic DFTK SCF result into DFTData{1}.
=#
function _load_dftk_nonmagnetic(scfres)
    basis = scfres.basis
    model = basis.model

    # Lattice vectors (already in Bohr, stored as columns)
    lattice = Matrix{Float64}(model.lattice)

    # K-points: fractional coordinates (3 × nkpts)
    nkpts = length(basis.kpoints)
    kpoints = zeros(3, nkpts)
    for (ik, kpt) in enumerate(basis.kpoints)
        kpoints[:, ik] = kpt.coordinate
    end

    # K-point weights (from symmetry reduction)
    weights = Float64.(basis.kweights)

    # Eigenvalues: (nbands, nkpts, nspin=1)
    nbands = length(scfres.eigenvalues[1])
    ebands = zeros(nbands, nkpts, 1)
    for ik = 1:nkpts
        ebands[:, ik, 1] = scfres.eigenvalues[ik]
    end

    # Occupations: same structure as eigenvalues
    occupations = zeros(nbands, nkpts, 1)
    for ik = 1:nkpts
        occupations[:, ik, 1] = scfres.occupation[ik]
    end

    # Fermi energy (already in Hartree)
    fermi = scfres.εF

    # Number of electrons
    nelect = Float64(model.n_electrons)

    # Atomic positions and species
    positions, species = _extract_atoms(model)

    # Return DFTData{1} (non-magnetic)
    return BoltzTraP.DFTData(;
        lattice = lattice,
        positions = positions,
        species = species,
        kpoints = kpoints,
        weights = weights,
        ebands = ebands,
        occupations = occupations,
        fermi = fermi,
        nelect = nelect,
    )
end

#=
    _load_dftk_collinear(scfres) -> DFTData{2}

Load collinear magnetic DFTK SCF result into DFTData{2}.

In DFTK collinear calculations, basis.kpoints contains each k-point twice:
first all spin-up k-points, then all spin-down k-points.
=#
function _load_dftk_collinear(scfres)
    basis = scfres.basis
    model = basis.model

    # Lattice vectors (already in Bohr, stored as columns)
    lattice = Matrix{Float64}(model.lattice)

    # For collinear: kpoints list has [up_1, up_2, ..., up_n, down_1, down_2, ..., down_n]
    # Each k-point coordinate appears twice (once for each spin)
    nkpts_total = length(basis.kpoints)
    nkpts = nkpts_total ÷ 2  # Actual number of unique k-points

    # K-points: fractional coordinates (3 × nkpts) - take only spin-up half
    kpoints = zeros(3, nkpts)
    for ik = 1:nkpts
        kpoints[:, ik] = basis.kpoints[ik].coordinate
    end

    # K-point weights - take only spin-up half (they should be the same)
    weights = Float64.(basis.kweights[1:nkpts])

    # Eigenvalues: (nbands, nkpts, 2)
    # First half of scfres.eigenvalues is spin-up, second half is spin-down
    nbands = length(scfres.eigenvalues[1])
    ebands = zeros(nbands, nkpts, 2)
    for ik = 1:nkpts
        ebands[:, ik, 1] = scfres.eigenvalues[ik]          # spin-up
        ebands[:, ik, 2] = scfres.eigenvalues[ik+nkpts]  # spin-down
    end

    # Occupations: same structure
    occupations = zeros(nbands, nkpts, 2)
    for ik = 1:nkpts
        occupations[:, ik, 1] = scfres.occupation[ik]          # spin-up
        occupations[:, ik, 2] = scfres.occupation[ik+nkpts]  # spin-down
    end

    # Fermi energy (already in Hartree)
    fermi = scfres.εF

    # Number of electrons
    nelect = Float64(model.n_electrons)

    # Atomic positions and species
    positions, species = _extract_atoms(model)

    # Compute magnetic moments from occupation difference
    natoms = length(model.positions)
    magmom = _compute_magmom(scfres, nkpts, natoms)

    # Return DFTData{2} (collinear magnetic)
    return BoltzTraP.DFTData(;
        lattice = lattice,
        positions = positions,
        species = species,
        kpoints = kpoints,
        weights = weights,
        ebands = ebands,
        occupations = occupations,
        fermi = fermi,
        nelect = nelect,
        magmom = magmom,
    )
end

#=
    _extract_atoms(model) -> (positions, species)

Extract atomic positions and species from DFTK model.
=#
function _extract_atoms(model)
    natoms = length(model.positions)
    positions = zeros(3, natoms)
    for (ia, pos) in enumerate(model.positions)
        positions[:, ia] = pos
    end

    species = String[]
    for atom in model.atoms
        sym = string(atom.species)
        push!(species, sym)
    end

    return positions, species
end

#=
    _compute_magmom(scfres, nkpts, natoms) -> Vector{Float64}

Compute magnetic moments per atom from SCF result.
Returns total spin moment distributed equally among atoms as approximation.
=#
function _compute_magmom(scfres, nkpts, natoms)
    basis = scfres.basis

    # Compute total magnetic moment from occupation difference
    # For each k-point pair (up, down), sum occupation differences
    total_up = 0.0
    total_down = 0.0
    for ik = 1:nkpts
        total_up += sum(scfres.occupation[ik]) * basis.kweights[ik]
        total_down += sum(scfres.occupation[ik+nkpts]) * basis.kweights[ik+nkpts]
    end

    # Total spin moment (in units of electrons)
    total_moment = total_up - total_down

    # Distribute equally among atoms (simple approximation)
    # For accurate per-atom moments, would need Mulliken or Bader analysis
    magmom = fill(total_moment / natoms, natoms)

    return magmom
end

end # module BoltzTraPDFTKExt
