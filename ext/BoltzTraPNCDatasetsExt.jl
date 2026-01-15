# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

#=
NCDatasets.jl extension for BoltzTraP.jl.

Provide `load_abinit` function to read ABINIT NetCDF output files (*_GSR.nc).
ABINIT uses atomic units (Hartree, Bohr) internally, so no unit conversion is needed.
=#
module BoltzTraPNCDatasetsExt

using BoltzTraP
using BoltzTraP: DFTData
using NCDatasets

#=
    load_abinit(directory::String) -> DFTData{NSpin}

Load DFT data from ABINIT NetCDF Ground State Results file.

# Arguments
- `directory`: Path to directory containing `*_GSR.nc` file

# Returns
- `DFTData{NSpin}`: DFT calculation data

# File Format (actual structure from Si.abinit)
- `primitive_vectors`: (3, 3) Lattice vectors in Bohr
- `reduced_atom_positions`: (3, natom) Fractional coordinates
- `atom_species`: (natom,) Species index per atom (1-based)
- `atom_species_names`: (80, ntypat) Element names as Char array
- `reduced_coordinates_of_kpoints`: (3, nkpt) K-points in fractional
- `eigenvalues`: (nband, nkpt, nspin) Band energies in Hartree
- `fermie`: Fermi energy in Hartree
- `nelect`: Number of electrons
=#
function BoltzTraP.load_abinit(directory::String)
    # Find GSR.nc file
    gsr_file = _find_gsr_file(directory)

    NCDataset(gsr_file, "r") do ds
        # Lattice vectors: (3, 3) in Bohr - column vectors
        # NCDatasets.jl handles column-major conversion automatically
        lattice = Matrix{Float64}(ds["primitive_vectors"][:, :])

        # Atomic positions: (3, natom) fractional - already in correct format
        positions = Matrix{Float64}(ds["reduced_atom_positions"][:, :])

        # Species: atom_species is 1-based index into atom_species_names
        atom_species_idx = ds["atom_species"][:]
        atom_species_names_raw = ds["atom_species_names"][:, :]

        # Parse species names from null-padded Char array
        species = String[]
        for idx in atom_species_idx
            # Get the column for this species type
            chars = atom_species_names_raw[:, idx]
            # Convert to string, removing null padding
            name = _chars_to_string(chars)
            push!(species, name)
        end

        # K-points: (3, nkpt) fractional - already in correct format
        kpoints = Matrix{Float64}(ds["reduced_coordinates_of_kpoints"][:, :])

        # K-point weights
        nkpts = size(kpoints, 2)
        if haskey(ds, "kpoint_weights")
            weights = Vector{Float64}(ds["kpoint_weights"][:])
        else
            weights = ones(nkpts) / nkpts
        end

        # Eigenvalues: (nband, nkpt, nspin) in Hartree - already in correct format
        ebands = Array{Float64}(ds["eigenvalues"][:, :, :])

        # Fermi energy in Hartree
        fermi = Float64(ds["fermie"][])

        # Number of electrons
        nelect = Float64(ds["nelect"][])

        # Occupations: read from file or use placeholder
        if haskey(ds, "occupations")
            occupations = Array{Float64}(ds["occupations"][:, :, :])
        else
            occupations = zeros(size(ebands))
        end

        # Create DFTData
        return DFTData(
            lattice = lattice,
            positions = positions,
            species = species,
            kpoints = kpoints,
            weights = weights,
            ebands = ebands,
            occupations = occupations,
            fermi = fermi,
            nelect = nelect,
            magmom = nothing
        )
    end
end

#=
    _find_gsr_file(directory::String) -> String

Find the GSR.nc file in the given directory.
=#
function _find_gsr_file(directory::String)
    !isdir(directory) && error("Directory not found: $directory")

    for f in readdir(directory)
        if endswith(f, "_GSR.nc") && isfile(joinpath(directory, f))
            return joinpath(directory, f)
        end
    end

    error("No *_GSR.nc file found in: $directory")
end

#=
    _chars_to_string(chars::AbstractVector{Char}) -> String

Convert a null-padded Char array to a trimmed String.
=#
function _chars_to_string(chars::AbstractVector{Char})
    # Find first null character or end
    len = length(chars)
    for i in 1:len
        if chars[i] == '\0'
            len = i - 1
            break
        end
    end
    return strip(String(chars[1:len]))
end

end # module BoltzTraPNCDatasetsExt
