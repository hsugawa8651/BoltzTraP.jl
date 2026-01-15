# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Type definitions for BoltzTraP.jl.
"""

"""
    DFTData{NSpin}

Container for DFT calculation results with spin parameter.

# Type Parameter
- `NSpin::Int`: Number of spin channels
  - `1`: Non-spin-polarized (non-magnetic)
  - `2`: Spin-polarized (collinear magnetic, ISPIN=2)

# Fields
- `lattice::Matrix{Float64}`: Lattice vectors (3×3) in Bohr, columns are vectors
- `positions::Matrix{Float64}`: Atomic positions (3×natom) in fractional coordinates
- `species::Vector{String}`: Element symbols
- `kpoints::Matrix{Float64}`: K-points (3×nkpts) in fractional coordinates
- `weights::Vector{Float64}`: K-point weights (nkpts,)
- `ebands::Array{Float64,3}`: Band energies (nbands, nkpts, nspin) in Hartree
- `occupations::Array{Float64,3}`: Occupations (nbands, nkpts, nspin)
- `fermi::Float64`: Fermi energy in Hartree
- `nelect::Float64`: Number of electrons
- `magmom::Union{Nothing,Vector{Float64},Matrix{Float64}}`: Magnetic moments (optional)

# Notes
- Non-collinear (spinor) calculations are not supported in v0.1
- When `NSpin=2`, bands from both spin channels are stored separately
"""
struct DFTData{NSpin}
    lattice::Matrix{Float64}
    positions::Matrix{Float64}
    species::Vector{String}
    kpoints::Matrix{Float64}
    weights::Vector{Float64}
    ebands::Array{Float64,3}
    occupations::Array{Float64,3}
    fermi::Float64
    nelect::Float64
    magmom::Union{Nothing,Vector{Float64},Matrix{Float64}}
end

# Type aliases for clarity
const NonMagneticData = DFTData{1}
const SpinPolarizedData = DFTData{2}

#=
    DFTData(; lattice, positions, species, kpoints, weights,
              ebands, occupations, fermi, nelect, magmom=nothing) -> DFTData{NSpin}

Construct DFTData with NSpin inferred from ebands shape.
=#
function DFTData(;
    lattice::AbstractMatrix,
    positions::AbstractMatrix,
    species::AbstractVector{<:AbstractString},
    kpoints::AbstractMatrix,
    weights::AbstractVector,
    ebands::AbstractArray{<:Real,3},
    occupations::AbstractArray{<:Real,3},
    fermi::Real,
    nelect::Real,
    magmom::Union{Nothing,AbstractVector,AbstractMatrix} = nothing,
)
    nspin = size(ebands, 3)
    DFTData{nspin}(
        Matrix{Float64}(lattice),
        Matrix{Float64}(positions),
        Vector{String}(species),
        Matrix{Float64}(kpoints),
        Vector{Float64}(weights),
        Array{Float64,3}(ebands),
        Array{Float64,3}(occupations),
        Float64(fermi),
        Float64(nelect),
        isnothing(magmom) ? nothing :
        (magmom isa AbstractMatrix ? Matrix{Float64}(magmom) : Vector{Float64}(magmom)),
    )
end

"""
    nspin(data::DFTData{N}) -> Int

Return the number of spin channels.
"""
nspin(::DFTData{N}) where {N} = N

"""
    is_magnetic(data::DFTData) -> Bool

Return `true` if data is from spin-polarized calculation.
"""
is_magnetic(data::DFTData) = nspin(data) > 1
