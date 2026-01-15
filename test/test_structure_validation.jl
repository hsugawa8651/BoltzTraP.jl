# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Crystal structure validation tests.

Tests that all loaders correctly read crystal structure data,
and that the same material produces consistent results across formats.

Note:
- Different DFT codes use different unit cell representations (primitive vs conventional)
- This affects lattice type detection and atom counts
- Volume per atom is used for cross-loader comparison
"""

using Test
using BoltzTraP
using LinearAlgebra

const BOLTZTRAP_DATA = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data")

# Check if NCDatasets is available for ABINIT tests
const HAS_NCDATASETS = try
    using NCDatasets
    true
catch
    false
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    validate_lattice(lattice)

Validate lattice vectors.

# Checks
- Shape is (3, 3)
- Vectors are non-zero
- Volume is positive
"""
function validate_lattice(lattice)
    @test size(lattice) == (3, 3)

    # Each lattice vector should have non-zero length
    for i in 1:3
        @test norm(lattice[:, i]) > 1e-10
    end

    # Volume should be positive
    volume = abs(det(lattice))
    @test volume > 1e-10

    return volume
end

"""
    validate_positions(positions, natoms)

Validate atomic positions.

# Checks
- Shape is (3, natoms)
- Fractional coordinates are in valid range
"""
function validate_positions(positions, natoms)
    @test size(positions) == (3, natoms)

    # Fractional coordinates should be approximately in [0, 1)
    # Allow small tolerance for floating point errors
    @test all(-1e-6 .<= positions .< 1.0 + 1e-6)
end

"""
    validate_species(species, natoms)

Validate species list.

# Checks
- Length matches natoms
- All species are non-empty strings
"""
function validate_species(species, natoms)
    @test length(species) == natoms
    @test all(s -> !isempty(s) && s isa AbstractString, species)
end

"""
    get_lattice_angles(lattice)

Get lattice angles (α, β, γ) in degrees.
"""
function get_lattice_angles(lattice)
    a1, a2, a3 = lattice[:, 1], lattice[:, 2], lattice[:, 3]
    len1, len2, len3 = norm(a1), norm(a2), norm(a3)

    cos_alpha = dot(a2, a3) / (len2 * len3)
    cos_beta = dot(a1, a3) / (len1 * len3)
    cos_gamma = dot(a1, a2) / (len1 * len2)

    alpha = acosd(clamp(cos_alpha, -1, 1))
    beta = acosd(clamp(cos_beta, -1, 1))
    gamma = acosd(clamp(cos_gamma, -1, 1))

    return (alpha=alpha, beta=beta, gamma=gamma)
end

"""
    validate_structure(data; min_atoms=1)

Comprehensive structure validation.

# Returns
NamedTuple with natoms, volume, volume_per_atom, species_counts
"""
function validate_structure(data; min_atoms=1)
    natoms = size(data.positions, 2)
    @test natoms >= min_atoms

    volume = validate_lattice(data.lattice)
    validate_positions(data.positions, natoms)
    validate_species(data.species, natoms)

    # Compute species counts
    species_counts = Dict{String,Int}()
    for s in data.species
        species_counts[s] = get(species_counts, s, 0) + 1
    end

    volume_per_atom = volume / natoms
    angles = get_lattice_angles(data.lattice)

    return (natoms=natoms, volume=volume, volume_per_atom=volume_per_atom,
            species_counts=species_counts, angles=angles)
end

"""
    compare_structures_relaxed(data1, data2; vol_rtol=0.05)

Compare structures from two loaders for the same material (relaxed comparison).

# Checks
- Same species composition (element counts match)
- Similar volume per atom (within tolerance)
"""
function compare_structures_relaxed(data1, data2; vol_rtol=0.05)
    info1 = validate_structure(data1)
    info2 = validate_structure(data2)

    # Same species counts (order may differ, atom count may differ due to cell choice)
    @test info1.species_counts == info2.species_counts

    # Similar volume per atom
    @test isapprox(info1.volume_per_atom, info2.volume_per_atom; rtol=vol_rtol)

    return (info1=info1, info2=info2)
end

# =============================================================================
# Per-Loader Structure Tests
# =============================================================================

@testset "Structure Validation" begin

    @testset "VASP structures" begin
        test_cases = [
            # (dirname, min_atoms, expected_species)
            ("Si.vasp", 2, Dict("Si" => 2)),
            ("Li.vasp", 1, Dict("Li" => 1)),  # primitive BCC cell
            ("PbTe.vasp.unpolarized", 2, Dict("Pb" => 1, "Te" => 1)),
        ]

        for (dirname, min_atoms, expected_species) in test_cases
            @testset "$dirname" begin
                data_path = joinpath(BOLTZTRAP_DATA, dirname)
                if !isdir(data_path)
                    @test_skip "Test data not found: $dirname"
                    continue
                end

                data = load_vasp(data_path)
                info = validate_structure(data; min_atoms=min_atoms)

                # Check species composition
                @test info.species_counts == expected_species
            end
        end
    end

    @testset "QE structures" begin
        test_cases = [
            ("Si.ESPRESSO/out", 2, Dict("Si" => 2)),
            ("Fe.ESPRESSO.collinear/out", 1, Dict("Fe" => 1)),
            ("nitinol.ESPRESSO/out", 4, Dict("Ni" => 2, "Ti" => 2)),  # Monoclinic NiTi
        ]

        for (dirname, min_atoms, expected_species) in test_cases
            @testset "$dirname" begin
                data_path = joinpath(BOLTZTRAP_DATA, dirname)
                if !isdir(data_path)
                    @test_skip "Test data not found: $dirname"
                    continue
                end

                data = load_qe(data_path)
                info = validate_structure(data; min_atoms=min_atoms)

                @test info.species_counts == expected_species
            end
        end
    end

    @testset "Wien2k structures" begin
        test_cases = [
            ("Si", 2, Dict("Si" => 2)),
            ("Li.W2K", 1, Dict("Li" => 1)),
            ("CoSb3", 16, Dict("Co" => 4, "Sb" => 12)),
            ("Bi2Te3", 5, Dict("Bi" => 2, "Te" => 3)),
        ]

        for (dirname, min_atoms, expected_species) in test_cases
            @testset "$dirname" begin
                data_path = joinpath(BOLTZTRAP_DATA, dirname)
                if !isdir(data_path)
                    @test_skip "Test data not found: $dirname"
                    continue
                end

                data = load_wien2k(data_path)
                info = validate_structure(data; min_atoms=min_atoms)

                @test info.species_counts == expected_species
            end
        end
    end

    @testset "ABINIT structures" begin
        # Note: ABINIT requires NCDatasets.jl extension
        @testset "Si.abinit" begin
            data_path = joinpath(BOLTZTRAP_DATA, "Si.abinit")
            if !isdir(data_path)
                @test_skip "Test data not found: Si.abinit"
            elseif !HAS_NCDATASETS
                @test_skip "load_abinit not available (NCDatasets extension not loaded)"
            else
                data = load_abinit(data_path)
                info = validate_structure(data; min_atoms=2)

                @test info.species_counts == Dict("Si" => 2)
            end
        end
    end

    @testset "GENE structures" begin
        @testset "Li.GENE.fromvasp (full load)" begin
            data_path = joinpath(BOLTZTRAP_DATA, "Li.GENE.fromvasp")
            if !isdir(data_path)
                @test_skip "Test data not found"
            else
                data = load_gene(data_path)
                info = validate_structure(data; min_atoms=1)

                # Li should be present
                @test haskey(info.species_counts, "Li")
                @test info.species_counts["Li"] >= 1
            end
        end

    end

    # =========================================================================
    # Cross-Loader Consistency Tests
    # =========================================================================

    @testset "Cross-loader consistency" begin

        @testset "Li: VASP vs Wien2k vs GENE" begin
            vasp_path = joinpath(BOLTZTRAP_DATA, "Li.vasp")
            wien2k_path = joinpath(BOLTZTRAP_DATA, "Li.W2K")
            gene_path = joinpath(BOLTZTRAP_DATA, "Li.GENE.fromvasp")

            skip_test = !isdir(vasp_path) || !isdir(wien2k_path) || !isdir(gene_path)
            if skip_test
                @test_skip "Some test data not found"
            else
                vasp_data = load_vasp(vasp_path)
                wien2k_data = load_wien2k(wien2k_path)
                gene_data = load_gene(gene_path)

                info_vasp = validate_structure(vasp_data)
                info_wien2k = validate_structure(wien2k_data)
                info_gene = validate_structure(gene_data)

                # All should have Li
                @test haskey(info_vasp.species_counts, "Li")
                @test haskey(info_wien2k.species_counts, "Li")
                @test haskey(info_gene.species_counts, "Li")

                # Volume per atom should be similar (within 5%)
                @test isapprox(info_vasp.volume_per_atom, info_wien2k.volume_per_atom; rtol=0.05)
                @test isapprox(info_vasp.volume_per_atom, info_gene.volume_per_atom; rtol=0.05)
            end
        end

        @testset "Si: VASP vs Wien2k vs QE vs ABINIT" begin
            vasp_path = joinpath(BOLTZTRAP_DATA, "Si.vasp")
            wien2k_path = joinpath(BOLTZTRAP_DATA, "Si")
            qe_path = joinpath(BOLTZTRAP_DATA, "Si.ESPRESSO", "out")
            abinit_path = joinpath(BOLTZTRAP_DATA, "Si.abinit")

            skip_test = !isdir(vasp_path) || !isdir(wien2k_path) || !isdir(qe_path)
            if skip_test
                @test_skip "Some test data not found"
            else
                vasp_data = load_vasp(vasp_path)
                wien2k_data = load_wien2k(wien2k_path)
                qe_data = load_qe(qe_path)

                info_vasp = validate_structure(vasp_data)
                info_wien2k = validate_structure(wien2k_data)
                info_qe = validate_structure(qe_data)

                # All should have Si
                @test haskey(info_vasp.species_counts, "Si")
                @test haskey(info_wien2k.species_counts, "Si")
                @test haskey(info_qe.species_counts, "Si")

                # Volume per atom should be similar (within 5%)
                @test isapprox(info_vasp.volume_per_atom, info_wien2k.volume_per_atom; rtol=0.05)
                @test isapprox(info_vasp.volume_per_atom, info_qe.volume_per_atom; rtol=0.05)

                # Include ABINIT if NCDatasets extension is loaded
                if isdir(abinit_path) && HAS_NCDATASETS
                    abinit_data = load_abinit(abinit_path)
                    info_abinit = validate_structure(abinit_data)

                    @test haskey(info_abinit.species_counts, "Si")
                    @test isapprox(info_vasp.volume_per_atom, info_abinit.volume_per_atom; rtol=0.05)
                end
            end
        end

    end
end
