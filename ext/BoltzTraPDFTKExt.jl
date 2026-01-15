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
    load_dftk(scfres) -> DFTData{1}

Extract band structure data from DFTK SCF result.

DFTK uses atomic units internally, so no unit conversion is needed.
Return `DFTData{1}` (non-magnetic) in the same format as other loaders.

# Arguments
- `scfres`: SCF result from DFTK `self_consistent_field`

# Returns
`DFTData{1}` with fields:
- `lattice`: Lattice vectors (3×3) in Bohr (columns are vectors)
- `positions`: Atomic positions (3×natoms) in fractional coordinates
- `species`: Atom species names
- `kpoints`: K-points (3×nkpts) in fractional coordinates
- `weights`: K-point weights (nkpts,)
- `ebands`: Band energies (nbands, nkpts, 1) in Hartree
- `occupations`: Occupations (nbands, nkpts, 1)
- `fermi`: Fermi energy in Hartree
- `nelect`: Number of electrons

# Example
```julia
using DFTK
using BoltzTraP

# Run DFT calculation
model = model_LDA(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut=30, kgrid=[8, 8, 8])
scfres = self_consistent_field(basis)

# Extract data for BoltzTraP
data = load_dftk(scfres)
@assert data isa DFTData{1}

# Run interpolation
result = run_interpolate(data)
```

# Notes
- Only non-spin-polarized calculations are supported
- Dense k-grid is required for accurate interpolation (10x10x10 or more recommended)
"""
function BoltzTraP.load_dftk(scfres)
    basis = scfres.basis
    model = basis.model

    # Validate: only non-spin-polarized supported
    if model.spin_polarization != :none
        error("Spin-polarized calculations not supported. " *
              "Use model.spin_polarization = :none")
    end

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
    # scfres.eigenvalues is Vector{Vector{Float64}} - eigenvalues per k-point
    nbands = length(scfres.eigenvalues[1])
    ebands = zeros(nbands, nkpts, 1)
    for ik in 1:nkpts
        ebands[:, ik, 1] = scfres.eigenvalues[ik]
    end

    # Occupations: same structure as eigenvalues
    occupations = zeros(nbands, nkpts, 1)
    for ik in 1:nkpts
        occupations[:, ik, 1] = scfres.occupation[ik]
    end

    # Fermi energy (already in Hartree)
    fermi = scfres.εF

    # Number of electrons
    nelect = Float64(model.n_electrons)

    # Atomic positions (fractional coordinates, 3 × natoms)
    natoms = length(model.positions)
    positions = zeros(3, natoms)
    for (ia, pos) in enumerate(model.positions)
        positions[:, ia] = pos
    end

    # Species names (extract element symbol from atom type)
    # DFTK atoms have a .species field (AtomsBase.ChemicalSpecies)
    species = String[]
    for atom in model.atoms
        sym = string(atom.species)
        push!(species, sym)
    end

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

end # module BoltzTraPDFTKExt
