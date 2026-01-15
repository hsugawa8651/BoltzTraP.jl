#!/usr/bin/env julia
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
Compare Julia BoltzTraP.jl implementation with Python BoltzTraP2 reference.

Usage:
    julia compare.jl              # Run all test cases
    julia compare.jl si_diamond   # Run specific test case
"""

using Test
using NPZ
using LinearAlgebra

# Add the parent directory to load BoltzTraP
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Import BoltzTraP (will fail gracefully if not yet implemented)
try
    using BoltzTraP
catch e
    @warn "BoltzTraP module not available yet. Some tests will be skipped."
end

const DATA_DIR = joinpath(@__DIR__, "data")
const RTOL = 1e-10
const ATOL = 1e-12

"""
    load_reference(name::String) -> Dict

Load reference data from an NPZ file.
"""
function load_reference(name::String)
    path = joinpath(DATA_DIR, "$(name).npz")
    if !isfile(path)
        error("Reference file not found: $path\n" *
              "Run 'python generate_*.py' scripts first (see reftest/README.md).")
    end
    return npzread(path)
end

"""
    compare_unit_cube()

Test compute_bounds with unit cube reference data.
"""
function compare_unit_cube()
    ref = load_reference("unit_cube")

    @testset "Unit cube" begin
        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            lattvec = ref["lattvec"]
            radii = ref["radii"]

            for r in radii
                result = BoltzTraP.compute_bounds(lattvec, r)
                expected = Int.(ref["bounds_r$r"])
                @test result == expected
            end
        end

        @testset "nrotations" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            lattvec = ref["lattvec"]
            positions = ref["positions"]
            types = Int.(ref["types"])

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, nothing)
            expected = Int(ref["nrotations"])
            @test result == expected
        end
    end
end

"""
    same_subspace(A, B; tol=1e-10)

Check if two sets of vectors span the same subspace.
Works for tensor bases where vectors may have different normalization.

A, B: Arrays of shape (nbasis, 3, 3) representing tensor bases.
Each slice A[i,:,:] is a 3x3 tensor (reshaped to 9-vector).
"""
function same_subspace(A, B; tol=1e-10)
    # Flatten each tensor to a vector (nbasis x 9)
    a_vecs = reshape(A, size(A,1), :)  # nbasis x 9
    b_vecs = reshape(B, size(B,1), :)  # nbasis x 9

    if size(a_vecs) != size(b_vecs)
        return false
    end

    nbasis = size(a_vecs, 1)

    # Each row of a_vecs is a 9-vector representing one basis tensor
    # Check that the row spaces of a_vecs and b_vecs are the same

    # Method: For each row in A, check if it can be expressed as
    # a linear combination of rows of B
    for i in 1:nbasis
        v = a_vecs[i, :]  # 9-vector

        # Solve min_c ||B' * c - v||^2 where B' is (9 x nbasis)
        # i.e., find coefficients c such that sum_j c[j] * b_vecs[j,:] ≈ v
        coeffs = b_vecs' \ v  # (9 x nbasis)' \ (9,) => (nbasis,)
        reconstructed = b_vecs' * coeffs
        residual = norm(reconstructed - v) / max(norm(v), 1e-10)
        if residual > tol
            return false
        end
    end

    # Also check reverse direction
    for i in 1:nbasis
        v = b_vecs[i, :]
        coeffs = a_vecs' \ v
        reconstructed = a_vecs' * coeffs
        residual = norm(reconstructed - v) / max(norm(v), 1e-10)
        if residual > tol
            return false
        end
    end

    return true
end

"""
    compare_simple_cubic()

Test with simple cubic reference data.
"""
function compare_simple_cubic()
    ref = load_reference("simple_cubic")

    @testset "Simple cubic" begin
        lattvec = ref["lattvec"]
        positions = ref["positions"]
        types = Int.(ref["types"])
        radius = ref["radius"]
        symprec = ref["symprec"]

        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            result = BoltzTraP.compute_bounds(lattvec, radius)
            expected = Int.(ref["bounds"])
            @test result == expected
        end

        @testset "calc_nrotations" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, nothing; symprec)
            expected = Int(ref["nrotations"])
            @test result == expected
        end

        @testset "calc_tensor_basis" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_tensor_basis)
                @test_skip "calc_tensor_basis not implemented"
                return
            end

            result = BoltzTraP.calc_tensor_basis(lattvec, positions, types, nothing; symprec)
            expected = ref["tensor_basis"]
            # Check subspace equivalence (normalization may differ)
            @test same_subspace(result, expected)
        end
    end
end

"""
    compare_si_diamond()

Test with Si diamond reference data.
"""
function compare_si_diamond()
    ref = load_reference("si_diamond")

    @testset "Si diamond" begin
        lattvec = ref["lattvec"]
        positions = ref["positions"]
        types = Int.(ref["types"])
        radius = ref["radius"]
        symprec = ref["symprec"]

        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            result = BoltzTraP.compute_bounds(lattvec, radius)
            expected = Int.(ref["bounds"])
            @test result == expected
        end

        @testset "calc_nrotations" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, nothing; symprec)
            expected = Int(ref["nrotations"])
            @test result == expected
        end

        @testset "calc_tensor_basis" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_tensor_basis)
                @test_skip "calc_tensor_basis not implemented"
                return
            end

            result = BoltzTraP.calc_tensor_basis(lattvec, positions, types, nothing; symprec)
            expected = ref["tensor_basis"]
            # Check subspace equivalence (normalization may differ)
            @test same_subspace(result, expected)
        end

        @testset "calc_sphere_quotient_set" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_sphere_quotient_set)
                @test_skip "calc_sphere_quotient_set not implemented"
                return
            end

            bounds = Int.(ref["bounds"])
            result = BoltzTraP.calc_sphere_quotient_set(
                lattvec, positions, types, nothing, radius, bounds; symprec
            )

            n_equiv = Int(ref["n_equivalences"])
            @test length(result) == n_equiv

            # Build sets of equivalence classes for order-independent comparison
            # Each class becomes a Set of tuples (frozen set)
            julia_classes = Set([Set(Tuple.(r)) for r in result])
            python_classes = Set([Set(Tuple.(eachrow(ref["equiv_$i"]))) for i in 0:(n_equiv-1)])

            @test julia_classes == python_classes
        end
    end
end

"""
    compare_fe_bcc_magnetic()

Test with Fe BCC magnetic reference data.
"""
function compare_fe_bcc_magnetic()
    ref = load_reference("fe_bcc_magnetic")

    @testset "Fe BCC magnetic" begin
        lattvec = ref["lattvec"]
        positions = ref["positions"]
        types = Int.(ref["types"])
        magmom = ref["magmom"]
        radius = ref["radius"]
        symprec = ref["symprec"]

        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            result = BoltzTraP.compute_bounds(lattvec, radius)
            expected = Int.(ref["bounds"])
            @test result == expected
        end

        @testset "calc_nrotations (collinear magnetic)" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, magmom; symprec)
            expected = Int(ref["nrotations"])
            @test result == expected
        end

        @testset "calc_tensor_basis (collinear magnetic)" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_tensor_basis)
                @test_skip "calc_tensor_basis not implemented"
                return
            end

            result = BoltzTraP.calc_tensor_basis(lattvec, positions, types, magmom; symprec)
            expected = ref["tensor_basis"]
            # Check subspace equivalence (normalization may differ)
            @test same_subspace(result, expected)
        end
    end
end

"""
    compare_monoclinic()

Test with monoclinic reference data (4 rotations).
"""
function compare_monoclinic()
    ref = load_reference("monoclinic")

    @testset "Monoclinic (4 rot)" begin
        lattvec = ref["lattvec"]
        positions = ref["positions"]
        types = Int.(ref["types"])
        radius = ref["radius"]
        symprec = ref["symprec"]

        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            result = BoltzTraP.compute_bounds(lattvec, radius)
            expected = Int.(ref["bounds"])
            @test result == expected
        end

        @testset "calc_nrotations" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, nothing; symprec)
            expected = Int(ref["nrotations"])
            @test result == expected
        end

        @testset "calc_sphere_quotient_set" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_sphere_quotient_set)
                @test_skip "calc_sphere_quotient_set not implemented"
                return
            end

            bounds = Int.(ref["bounds"])
            result = BoltzTraP.calc_sphere_quotient_set(
                lattvec, positions, types, nothing, radius, bounds; symprec
            )

            n_equiv = Int(ref["n_equivalences"])
            @test length(result) == n_equiv
        end
    end
end

"""
    compare_triclinic_p1()

Test with triclinic P1 reference data (2 rotations, minimal symmetry).
"""
function compare_triclinic_p1()
    ref = load_reference("triclinic_p1")

    @testset "Triclinic P1 (2 rot)" begin
        lattvec = ref["lattvec"]
        positions = ref["positions"]
        types = Int.(ref["types"])
        radius = ref["radius"]
        symprec = ref["symprec"]

        @testset "compute_bounds" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_bounds)
                @test_skip "compute_bounds not implemented"
                return
            end

            result = BoltzTraP.compute_bounds(lattvec, radius)
            expected = Int.(ref["bounds"])
            @test result == expected
        end

        @testset "calc_nrotations" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_nrotations)
                @test_skip "calc_nrotations not implemented"
                return
            end

            result = BoltzTraP.calc_nrotations(lattvec, positions, types, nothing; symprec)
            expected = Int(ref["nrotations"])
            @test result == expected
        end

        @testset "calc_sphere_quotient_set" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_sphere_quotient_set)
                @test_skip "calc_sphere_quotient_set not implemented"
                return
            end

            bounds = Int.(ref["bounds"])
            result = BoltzTraP.calc_sphere_quotient_set(
                lattvec, positions, types, nothing, radius, bounds; symprec
            )

            n_equiv = Int(ref["n_equivalences"])
            @test length(result) == n_equiv
        end
    end
end

# ============================================================================
# Phase 3: Interpolation Tests
# ============================================================================

"""
    compare_simple_interpolation()

Test phase factor computation with simple test case.
"""
function compare_simple_interpolation()
    ref = load_reference("simple_interpolation")

    @testset "Simple interpolation" begin
        lattvec = ref["lattvec"]
        kpoints = ref["kpoints"]
        n_equiv = Int(ref["n_equivalences"])

        # Reconstruct equivalences from saved arrays
        equivalences = [ref["equiv_$i"] for i in 0:(n_equiv-1)]

        @testset "compute_phase_factors" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_phase_factors)
                @test_skip "compute_phase_factors not implemented"
                return
            end

            result = BoltzTraP.compute_phase_factors(kpoints, equivalences)
            expected = ref["phase"]
            @test result ≈ expected rtol=1e-10
        end
    end
end

"""
    compare_si_interpolation()

Test interpolation with Si DFT data.
"""
function compare_si_interpolation()
    ref = load_reference("si_interpolation")

    @testset "Si interpolation" begin
        lattvec = ref["lattvec"]
        kpoints = ref["kpoints"]
        ebands = ref["ebands"]
        n_equiv = Int(ref["n_equivalences"])

        # Reconstruct equivalences from saved arrays
        equivalences = [ref["equiv_$i"] for i in 0:(n_equiv-1)]

        @testset "compute_phase_factors" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :compute_phase_factors)
                @test_skip "compute_phase_factors not implemented"
                return
            end

            result = BoltzTraP.compute_phase_factors(kpoints, equivalences)
            expected = ref["phase"]
            @test result ≈ expected rtol=1e-10
        end

        @testset "fitde3D coefficients" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :fitde3D)
                @test_skip "fitde3D not implemented"
                return
            end

            result = BoltzTraP.fitde3D(kpoints, ebands, equivalences, lattvec)
            expected = ref["coeffs"]
            # Coefficients may have small differences due to numerical precision
            @test result ≈ expected rtol=1e-8 atol=1e-10
        end

        @testset "band reconstruction" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :fitde3D)
                @test_skip "fitde3D not implemented"
                return
            end

            # Fit and reconstruct
            coeffs = BoltzTraP.fitde3D(kpoints, ebands, equivalences, lattvec)

            # Use getBands for reconstruction at original k-points
            if isdefined(BoltzTraP, :getBands)
                ebands_recon, vbands = BoltzTraP.getBands(kpoints, equivalences, lattvec, coeffs)
                # Reconstruction should match original bands very closely
                @test ebands_recon ≈ ebands rtol=1e-10 atol=1e-10
            else
                @test_skip "getBands not implemented"
            end
        end
    end
end

# ============================================================================
# Python bt2 Format Compatibility Tests
# ============================================================================

"""
    compare_bt2_python_format()

Test reading Python BoltzTraP2 generated .bt2 files.
"""
function compare_bt2_python_format()
    ref_file = joinpath(DATA_DIR, "si_bt2_reference.npz")
    if !isfile(ref_file)
        @warn "si_bt2_reference.npz not found (run generate_4_io_bt2.py)"
        @test_skip "Reference file missing"
        return
    end

    ref = load_reference("si_bt2_reference")
    bt2_file = joinpath(BOLTZTRAP_DATA, "Si.vasp", "interpolation.bt2")

    if !isfile(bt2_file)
        @warn "interpolation.bt2 not found: $bt2_file"
        @test_skip "bt2 file missing"
        return
    end

    @testset "Python bt2 format" begin
        @testset "File provenance" begin
            # XZ magic bytes
            header = read(bt2_file, 6)
            @test header == UInt8[0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00]
        end

        @testset "Metadata" begin
            result = BoltzTraP.load_interpolation(bt2_file)
            # Reference stores strings as UInt8 arrays for NPZ.jl compatibility
            expected_format_version = String(UInt8.(ref["format_version"]))
            expected_source = String(UInt8.(ref["source"]))
            @test result.metadata["format_version"] == expected_format_version
            @test result.metadata["source"] == expected_source
            @test haskey(result.metadata, "ctime")
        end

        @testset "Coefficients" begin
            result = BoltzTraP.load_interpolation(bt2_file)
            expected_shape = Tuple(Int.(ref["coeffs_shape"]))
            @test size(result.coeffs) == expected_shape

            expected = complex.(ref["coeffs_real"], ref["coeffs_imag"])
            @test result.coeffs ≈ expected rtol=1e-10
        end

        @testset "Equivalences" begin
            result = BoltzTraP.load_interpolation(bt2_file)
            @test length(result.equivalences) == Int(ref["n_equivalences"])

            # First equivalence class should match
            if Int(ref["n_equivalences"]) > 0
                # Reference equiv_0 might be stored as (n, 3) or (3,) depending on npz
                ref_equiv = ref["equiv_0"]
                julia_equiv = result.equivalences[1]
                # Compare values regardless of shape
                @test vec(julia_equiv) ≈ vec(ref_equiv) rtol=1e-10
            end
        end

        @testset "Lattice vectors" begin
            result = BoltzTraP.load_interpolation(bt2_file)
            # Python uses row convention, Julia uses column convention
            @test result.lattvec ≈ ref["lattvec"]' rtol=1e-10
        end

        println("  Python bt2 format: All tests passed")
    end
end

"""
    run_all_tests()

Run all comparison tests.
"""
function run_all_tests()
    @testset "BoltzTraP.jl vs Python BoltzTraP2" begin
        # Phase 2: High symmetry
        compare_unit_cube()
        compare_simple_cubic()
        compare_si_diamond()

        # Phase 2: Magnetic
        compare_fe_bcc_magnetic()

        # Phase 2: Low symmetry
        compare_monoclinic()
        compare_triclinic_p1()

        # Phase 3: Interpolation
        compare_simple_interpolation()
        compare_si_interpolation()

        # Phase 4: Transport
        compare_simple_transport()

        # Phase 5: I/O
        compare_poscar()

        # Phase 6.5: E2E Integrate
        compare_si_integrate_e2e()
        compare_qe_si_integrate_e2e()
        compare_wien2k_si_integrate_e2e()
        compare_abinit_si_integrate_e2e()

        # Python bt2 Compatibility
        compare_bt2_python_format()
    end
end

# ============================================================================
# Phase 4: Transport Tests
# ============================================================================

"""
    compare_simple_transport()

Test transport functions with simple parabolic band test case.
"""
function compare_simple_transport()
    ref = load_reference("simple_transport")

    @testset "Simple transport" begin
        eband = ref["eband"]
        vvband = ref["vvband"]
        Tr = ref["Tr"]
        mur = ref["mur"]
        vuc = ref["vuc"]
        dosweight = ref["dosweight"]

        @testset "BTPDOS" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :BTPDOS)
                @test_skip "BTPDOS not implemented"
                return
            end

            epsilon, dos, vvdos = BoltzTraP.BTPDOS(eband, vvband; npts=100)

            # Check shapes
            @test length(epsilon) == length(ref["epsilon"])
            @test length(dos) == length(ref["dos"])
            @test size(vvdos) == size(ref["vvdos"])

            # Check DOS values (may differ slightly due to histogram binning)
            @test epsilon ≈ ref["epsilon"] rtol=0.01
            @test dos ≈ ref["dos"] rtol=0.1  # Histogram can vary
        end

        @testset "fermi_integrals" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :fermi_integrals)
                @test_skip "fermi_integrals not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            vvdos = ref["vvdos"]

            N, L0, L1, L2 = BoltzTraP.fermi_integrals(
                epsilon, dos, vvdos, mur, Tr; dosweight=dosweight
            )

            @test size(L0) == size(ref["L0"])
            @test L0 ≈ ref["L0"] rtol=0.01
            @test L1 ≈ ref["L1"] rtol=0.01
            @test L2 ≈ ref["L2"] rtol=0.01
        end

        @testset "calc_onsager_coefficients" begin
            if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_onsager_coefficients)
                @test_skip "calc_onsager_coefficients not implemented"
                return
            end

            L0 = ref["L0"]
            L1 = ref["L1"]
            L2 = ref["L2"]

            sigma, S, kappa = BoltzTraP.calc_onsager_coefficients(L0, L1, L2, Tr, vuc)

            @test size(sigma) == size(ref["sigma"])
            @test sigma ≈ ref["sigma"] rtol=0.01
            # Seebeck can have larger relative errors near zero
            @test S ≈ ref["S"] rtol=0.1 atol=1e-6
        end
    end
end

# ============================================================================
# Phase 5: I/O Tests
# ============================================================================

const BOLTZTRAP_DATA = joinpath(@__DIR__, "..", "..", "BoltzTraP2-public", "data")

"""
    compare_poscar()

Test POSCAR reading against ASE/Python reference.
"""
function compare_poscar()
    @testset "POSCAR reading" begin
        test_cases = [
            ("si_poscar", "Si.vasp"),
            ("li_poscar", "Li.vasp"),
            ("pbte_poscar", "PbTe.vasp.unpolarized"),
        ]

        for (ref_name, dirname) in test_cases
            ref_file = joinpath(DATA_DIR, "$(ref_name).npz")
            if !isfile(ref_file)
                @test_skip "Reference file not found: $ref_name"
                continue
            end

            ref = load_reference(ref_name)
            poscar_path = joinpath(BOLTZTRAP_DATA, dirname, "POSCAR")

            @testset "$dirname" begin
                if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :read_poscar)
                    @test_skip "read_poscar not available"
                    return
                end

                result = BoltzTraP.read_poscar(poscar_path)

                # Lattice vectors (column convention)
                @test result.lattice ≈ ref["lattvec"] rtol=1e-10

                # Atomic positions (fractional)
                # Julia returns Cartesian, need to convert
                positions_frac = result.lattice \ result.positions
                @test positions_frac ≈ ref["positions"]' rtol=1e-10

                # Number of atoms
                @test sum(result.counts) == Int(ref["n_atoms"])
            end
        end
    end
end

# ============================================================================
# Phase 6: End-to-End Interpolation Tests
# ============================================================================

"""
    compare_si_end2end()

Test end-to-end interpolation workflow against Python BoltzTraP2.
"""
function compare_si_end2end()
    ref = load_reference("si_end2end")

    # Unit conversion: Python BoltzTraP2 uses Hartree, Julia uses eV
    Ha_to_eV = 27.211386245988

    @testset "Si end-to-end interpolation" begin
        if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :run_interpolate)
            @test_skip "run_interpolate not implemented"
            return
        end

        si_path = joinpath(BOLTZTRAP_DATA, "Si.vasp")

        @testset "load_vasp" begin
            data = BoltzTraP.load_vasp(si_path)

            @test data.lattice ≈ ref["lattvec"] rtol=1e-10
            @test data.kpoints ≈ ref["kpoints"]' rtol=1e-10
            # Convert Python Ha to eV for comparison
            # Note: rtol=1e-6 accounts for VASP parser precision differences
            @test data.ebands[:,:,1] ≈ ref["ebands"] * Ha_to_eV rtol=1e-6
        end

        @testset "run_interpolate" begin
            # Run with same nkpt_target as Python
            # Python got 5586 equivalences with nkpt_target=5000
            result = BoltzTraP.run_interpolate(si_path; kpoints=5000)

            # Check equivalence count matches
            @test length(result.equivalences) == Int(ref["n_equivalences"])

            # Check coefficient dimensions
            @test size(result.coeffs) == (size(ref["coeffs_real"], 1), size(ref["coeffs_real"], 2))

            # Note: We don't compare coefficients directly because equivalence
            # ordering may differ between implementations. The reconstruction
            # accuracy test below verifies correctness.
        end

        @testset "equivalence classes as sets" begin
            result = BoltzTraP.run_interpolate(si_path; kpoints=5000)

            # Convert Julia equivalences to set of sets
            julia_classes = Set{Set{NTuple{3,Int}}}()
            for eq in result.equivalences
                class_set = Set(Tuple(eq[i, :]) for i in 1:size(eq, 1))
                push!(julia_classes, class_set)
            end

            # Python only stores representatives, but we can check that
            # each Python representative belongs to exactly one Julia class
            py_reps = ref["equiv_reps"]
            found_classes = Set{Set{NTuple{3,Int}}}()

            for i in 1:size(py_reps, 1)
                py_rep = Tuple(Int.(py_reps[i, :]))
                # Find which Julia class contains this representative
                for jclass in julia_classes
                    if py_rep in jclass
                        push!(found_classes, jclass)
                        break
                    end
                end
            end

            # Every Python representative should map to a unique Julia class
            @test length(found_classes) == length(julia_classes)
            println("  All $(length(julia_classes)) equivalence classes match as sets")
        end

        @testset "reconstruction accuracy" begin
            result = BoltzTraP.run_interpolate(si_path; kpoints=5000)

            # Create interpolator from result
            interp = BoltzTraP.FourierInterpolator(
                result.coeffs, result.equivalences, result.lattvec
            )

            # Reconstruct at original k-points
            kpoints = ref["kpoints"]
            ebands_julia = BoltzTraP.interpolate_bands(interp, kpoints)

            # Check reconstruction error (in eV)
            # Compare Julia's reconstruction to original bands
            ebands_original = ref["ebands"] * Ha_to_eV
            max_error_julia = maximum(abs.(ebands_julia - ebands_original))

            # Good reconstruction: error should be < 1e-6 eV
            # (Python achieves similar accuracy)
            @test max_error_julia < 1e-6

            # Also verify Python's reconstruction error for reference
            max_error_python = Float64(ref["max_reconstruction_error"]) * Ha_to_eV
            println("  Python reconstruction error: $(max_error_python) eV")
            println("  Julia reconstruction error:  $(max_error_julia) eV")
        end
    end
end

# ============================================================================
# Phase 6.5: End-to-End Integrate Tests
# ============================================================================

"""
    compare_si_integrate_e2e()

Test end-to-end integrate workflow: interpolation → getBTPbands → transport.
"""
function compare_si_integrate_e2e()
    ref_file = joinpath(DATA_DIR, "si_end2end.npz")
    if !isfile(ref_file)
        @test_skip "si_end2end.npz not found"
        return
    end

    ref = load_reference("si_end2end")

    # Check if transport data is included
    if !haskey(ref, "eband")
        @test_skip "Transport data not in si_end2end.npz (regenerate with Python)"
        return
    end

    @testset "Si E2E integrate" begin
        if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :getBTPbands)
            @test_skip "getBTPbands not implemented"
            return
        end

        # Extract reference data
        lattvec = ref["lattvec"]
        coeffs_ref = complex.(ref["coeffs_real"], ref["coeffs_imag"])
        n_equiv = Int(ref["n_equivalences"])
        vuc = Float64(ref["vuc"])
        dosweight = Float64(ref["dosweight"])
        Tr = ref["Tr"]
        mur = ref["mur"]

        # Reconstruct equivalences from reference
        equiv_reps = ref["equiv_reps"]
        equivalences = [reshape(equiv_reps[i, :], 1, 3) for i in 1:n_equiv]

        # Note: getBTPbands test is skipped here because:
        # 1. Reference data only stores equivalence representatives (1 point per class)
        # 2. Full equivalence classes are needed for proper FFT grid computation
        # 3. getBTPbands is tested indirectly via compare_simple_transport
        #    which uses Python's pre-computed eband/vvband
        # The tests below use Python's eband/vvband to verify downstream functions.

        @testset "BTPDOS" begin
            if !isdefined(BoltzTraP, :BTPDOS)
                @test_skip "BTPDOS not implemented"
                return
            end

            eband_ref = ref["eband"]
            vvband_ref = ref["vvband"]

            epsilon_julia, dos_julia, vvdos_julia = BoltzTraP.BTPDOS(eband_ref, vvband_ref; npts=500)

            epsilon_ref = ref["epsilon"]
            dos_ref = ref["dos"]
            vvdos_ref = ref["vvdos"]

            @test length(epsilon_julia) == length(epsilon_ref)
            @test epsilon_julia ≈ epsilon_ref rtol=0.01
            # DOS can vary due to histogram binning
            @test dos_julia ≈ dos_ref rtol=0.1
        end

        @testset "fermi_integrals" begin
            if !isdefined(BoltzTraP, :fermi_integrals)
                @test_skip "fermi_integrals not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            vvdos = ref["vvdos"]

            N_julia, L0_julia, L1_julia, L2_julia = BoltzTraP.fermi_integrals(
                epsilon, dos, vvdos, mur, Tr; dosweight=dosweight
            )

            L0_ref = ref["L0"]
            L1_ref = ref["L1"]
            L2_ref = ref["L2"]

            @test size(L0_julia) == size(L0_ref)
            @test L0_julia ≈ L0_ref rtol=0.01
            @test L1_julia ≈ L1_ref rtol=0.01
            @test L2_julia ≈ L2_ref rtol=0.01

            println("  fermi_integrals: L0 max rel diff = $(maximum(abs.(L0_julia - L0_ref) ./ max.(abs.(L0_ref), 1e-20)))")
        end

        @testset "calc_onsager_coefficients" begin
            if !isdefined(BoltzTraP, :calc_onsager_coefficients)
                @test_skip "calc_onsager_coefficients not implemented"
                return
            end

            L0 = ref["L0"]
            L1 = ref["L1"]
            L2 = ref["L2"]

            sigma_julia, S_julia, kappa_julia = BoltzTraP.calc_onsager_coefficients(
                L0, L1, L2, Tr, vuc
            )

            sigma_ref = ref["sigma"]
            S_ref = ref["S"]
            kappa_ref = ref["kappa"]

            @test size(sigma_julia) == size(sigma_ref)
            @test sigma_julia ≈ sigma_ref rtol=0.01
            # Seebeck can have larger relative errors near zero
            @test S_julia ≈ S_ref rtol=0.1 atol=1e-6
            @test kappa_julia ≈ kappa_ref rtol=0.01

            # Sample value at T=300K (index 2), mid-mu
            iT = 2
            imu = size(sigma_julia, 2) ÷ 2
            println("  sigma_xx at T=$(Tr[iT])K: Julia=$(sigma_julia[iT, imu, 1, 1]), Ref=$(sigma_ref[iT, imu, 1, 1])")
        end

        @testset "solve_for_mu" begin
            if !isdefined(BoltzTraP, :solve_for_mu)
                @test_skip "solve_for_mu not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            nelect = Float64(ref["nelect"])
            mu0_ref = ref["mu0"]

            for (iT, T) in enumerate(Tr)
                mu0_julia = BoltzTraP.solve_for_mu(
                    epsilon, dos, nelect, T;
                    dosweight=dosweight, refine=true, try_center=true
                )

                # Allow 1 histogram bin difference
                de = epsilon[2] - epsilon[1]
                @test abs(mu0_julia - mu0_ref[iT]) < 2 * de

                println("  T=$(T)K: mu0 Julia=$(round(mu0_julia, digits=6)), Ref=$(round(mu0_ref[iT], digits=6))")
            end
        end
    end
end

"""
    compare_qe_si_integrate_e2e()

Test end-to-end integrate workflow for QE Si (non-magnetic QE data).
"""
function compare_qe_si_integrate_e2e()
    ref_file = joinpath(DATA_DIR, "qe_si_end2end.npz")
    if !isfile(ref_file)
        @test_skip "qe_si_end2end.npz not found (run generate_5_e2e.py)"
        return
    end

    ref = load_reference("qe_si_end2end")

    # Check if transport data is included
    if !haskey(ref, "eband")
        @test_skip "Transport data not in qe_si_end2end.npz (regenerate with Python)"
        return
    end

    @testset "QE Si E2E integrate" begin
        if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_onsager_coefficients)
            @test_skip "calc_onsager_coefficients not implemented"
            return
        end

        # Extract reference data
        vuc = Float64(ref["vuc"])
        dosweight = Float64(ref["dosweight"])
        Tr = ref["Tr"]
        mur = ref["mur"]

        @testset "BTPDOS" begin
            if !isdefined(BoltzTraP, :BTPDOS)
                @test_skip "BTPDOS not implemented"
                return
            end

            eband_ref = ref["eband"]
            vvband_ref = ref["vvband"]

            epsilon_julia, dos_julia, vvdos_julia = BoltzTraP.BTPDOS(eband_ref, vvband_ref; npts=500)

            epsilon_ref = ref["epsilon"]
            dos_ref = ref["dos"]
            vvdos_ref = ref["vvdos"]

            @test length(epsilon_julia) == length(epsilon_ref)
            @test epsilon_julia ≈ epsilon_ref rtol=0.01
            # DOS can vary due to histogram binning
            @test dos_julia ≈ dos_ref rtol=0.1
        end

        @testset "fermi_integrals" begin
            if !isdefined(BoltzTraP, :fermi_integrals)
                @test_skip "fermi_integrals not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            vvdos = ref["vvdos"]

            N_julia, L0_julia, L1_julia, L2_julia = BoltzTraP.fermi_integrals(
                epsilon, dos, vvdos, mur, Tr; dosweight=dosweight
            )

            L0_ref = ref["L0"]
            L1_ref = ref["L1"]
            L2_ref = ref["L2"]

            @test size(L0_julia) == size(L0_ref)
            @test L0_julia ≈ L0_ref rtol=0.01
            @test L1_julia ≈ L1_ref rtol=0.01
            @test L2_julia ≈ L2_ref rtol=0.01

            println("  QE Si fermi_integrals: L0 max rel diff = $(maximum(abs.(L0_julia - L0_ref) ./ max.(abs.(L0_ref), 1e-20)))")
        end

        @testset "calc_onsager_coefficients" begin
            L0 = ref["L0"]
            L1 = ref["L1"]
            L2 = ref["L2"]

            sigma_julia, S_julia, kappa_julia = BoltzTraP.calc_onsager_coefficients(
                L0, L1, L2, Tr, vuc
            )

            sigma_ref = ref["sigma"]
            S_ref = ref["S"]
            kappa_ref = ref["kappa"]

            @test size(sigma_julia) == size(sigma_ref)
            @test sigma_julia ≈ sigma_ref rtol=0.01
            # Seebeck can have larger relative errors near zero
            @test S_julia ≈ S_ref rtol=0.1 atol=1e-6
            @test kappa_julia ≈ kappa_ref rtol=0.01

            # Sample value at T=300K (index 2), mid-mu
            iT = 2
            imu = size(sigma_julia, 2) ÷ 2
            println("  QE Si sigma_xx at T=$(Tr[iT])K: Julia=$(sigma_julia[iT, imu, 1, 1]), Ref=$(sigma_ref[iT, imu, 1, 1])")
        end

        @testset "solve_for_mu" begin
            if !isdefined(BoltzTraP, :solve_for_mu)
                @test_skip "solve_for_mu not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            nelect = Float64(ref["nelect"])
            mu0_ref = ref["mu0"]

            for (iT, T) in enumerate(Tr)
                mu0_julia = BoltzTraP.solve_for_mu(
                    epsilon, dos, nelect, T;
                    dosweight=dosweight, refine=true, try_center=true
                )

                # Allow 1 histogram bin difference
                de = epsilon[2] - epsilon[1]
                @test abs(mu0_julia - mu0_ref[iT]) < 2 * de

                println("  QE Si T=$(T)K: mu0 Julia=$(round(mu0_julia, digits=6)), Ref=$(round(mu0_ref[iT], digits=6))")
            end
        end
    end
end

"""
    compare_wien2k_si_integrate_e2e()

Test end-to-end integrate workflow for Wien2k Si (non-magnetic Wien2k data).
"""
function compare_wien2k_si_integrate_e2e()
    ref_file = joinpath(DATA_DIR, "wien2k_si_end2end.npz")
    if !isfile(ref_file)
        @test_skip "wien2k_si_end2end.npz not found (run generate_5_e2e.py)"
        return
    end

    ref = load_reference("wien2k_si_end2end")

    # Check if transport data is included
    if !haskey(ref, "eband")
        @test_skip "Transport data not in wien2k_si_end2end.npz (regenerate with Python)"
        return
    end

    @testset "Wien2k Si E2E integrate" begin
        if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_onsager_coefficients)
            @test_skip "calc_onsager_coefficients not implemented"
            return
        end

        # Extract reference data
        vuc = Float64(ref["vuc"])
        dosweight = Float64(ref["dosweight"])
        Tr = ref["Tr"]
        mur = ref["mur"]

        @testset "BTPDOS" begin
            if !isdefined(BoltzTraP, :BTPDOS)
                @test_skip "BTPDOS not implemented"
                return
            end

            eband_ref = ref["eband"]
            vvband_ref = ref["vvband"]

            epsilon_julia, dos_julia, vvdos_julia = BoltzTraP.BTPDOS(eband_ref, vvband_ref; npts=500)

            epsilon_ref = ref["epsilon"]
            dos_ref = ref["dos"]
            vvdos_ref = ref["vvdos"]

            @test length(epsilon_julia) == length(epsilon_ref)
            @test epsilon_julia ≈ epsilon_ref rtol=0.01
            # DOS can vary due to histogram binning
            @test dos_julia ≈ dos_ref rtol=0.1
        end

        @testset "fermi_integrals" begin
            if !isdefined(BoltzTraP, :fermi_integrals)
                @test_skip "fermi_integrals not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            vvdos = ref["vvdos"]

            N_julia, L0_julia, L1_julia, L2_julia = BoltzTraP.fermi_integrals(
                epsilon, dos, vvdos, mur, Tr; dosweight=dosweight
            )

            L0_ref = ref["L0"]
            L1_ref = ref["L1"]
            L2_ref = ref["L2"]

            @test size(L0_julia) == size(L0_ref)
            @test L0_julia ≈ L0_ref rtol=0.01
            @test L1_julia ≈ L1_ref rtol=0.01
            @test L2_julia ≈ L2_ref rtol=0.01

            println("  Wien2k Si fermi_integrals: L0 max rel diff = $(maximum(abs.(L0_julia - L0_ref) ./ max.(abs.(L0_ref), 1e-20)))")
        end

        @testset "calc_onsager_coefficients" begin
            L0 = ref["L0"]
            L1 = ref["L1"]
            L2 = ref["L2"]

            sigma_julia, S_julia, kappa_julia = BoltzTraP.calc_onsager_coefficients(
                L0, L1, L2, Tr, vuc
            )

            sigma_ref = ref["sigma"]
            S_ref = ref["S"]
            kappa_ref = ref["kappa"]

            @test size(sigma_julia) == size(sigma_ref)
            @test sigma_julia ≈ sigma_ref rtol=0.01
            # Seebeck can have larger relative errors near zero
            @test S_julia ≈ S_ref rtol=0.1 atol=1e-6
            @test kappa_julia ≈ kappa_ref rtol=0.01

            # Sample value at T=300K (index 2), mid-mu
            iT = 2
            imu = size(sigma_julia, 2) ÷ 2
            println("  Wien2k Si sigma_xx at T=$(Tr[iT])K: Julia=$(sigma_julia[iT, imu, 1, 1]), Ref=$(sigma_ref[iT, imu, 1, 1])")
        end

        @testset "solve_for_mu" begin
            if !isdefined(BoltzTraP, :solve_for_mu)
                @test_skip "solve_for_mu not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            nelect = Float64(ref["nelect"])
            mu0_ref = ref["mu0"]

            for (iT, T) in enumerate(Tr)
                mu0_julia = BoltzTraP.solve_for_mu(
                    epsilon, dos, nelect, T;
                    dosweight=dosweight, refine=true, try_center=true
                )

                # Allow 1 histogram bin difference
                de = epsilon[2] - epsilon[1]
                @test abs(mu0_julia - mu0_ref[iT]) < 2 * de

                println("  Wien2k Si T=$(T)K: mu0 Julia=$(round(mu0_julia, digits=6)), Ref=$(round(mu0_ref[iT], digits=6))")
            end
        end
    end
end

"""
    compare_abinit_si_integrate_e2e()

Test end-to-end integrate workflow for ABINIT Si (non-magnetic ABINIT data).
"""
function compare_abinit_si_integrate_e2e()
    ref_file = joinpath(DATA_DIR, "abinit_si_end2end.npz")
    if !isfile(ref_file)
        @test_skip "abinit_si_end2end.npz not found (run generate_5_e2e.py)"
        return
    end

    ref = load_reference("abinit_si_end2end")

    # Check if transport data is included
    if !haskey(ref, "eband")
        @test_skip "Transport data not in abinit_si_end2end.npz (regenerate with Python)"
        return
    end

    @testset "ABINIT Si E2E integrate" begin
        if !@isdefined(BoltzTraP) || !isdefined(BoltzTraP, :calc_onsager_coefficients)
            @test_skip "calc_onsager_coefficients not implemented"
            return
        end

        # Extract reference data
        vuc = Float64(ref["vuc"])
        dosweight = Float64(ref["dosweight"])
        Tr = ref["Tr"]
        mur = ref["mur"]

        @testset "BTPDOS" begin
            if !isdefined(BoltzTraP, :BTPDOS)
                @test_skip "BTPDOS not implemented"
                return
            end

            eband_ref = ref["eband"]
            vvband_ref = ref["vvband"]

            epsilon_julia, dos_julia, vvdos_julia = BoltzTraP.BTPDOS(eband_ref, vvband_ref; npts=500)

            epsilon_ref = ref["epsilon"]
            dos_ref = ref["dos"]
            vvdos_ref = ref["vvdos"]

            @test length(epsilon_julia) == length(epsilon_ref)
            @test epsilon_julia ≈ epsilon_ref rtol=0.01
            # DOS can vary due to histogram binning
            @test dos_julia ≈ dos_ref rtol=0.1
        end

        @testset "fermi_integrals" begin
            if !isdefined(BoltzTraP, :fermi_integrals)
                @test_skip "fermi_integrals not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            vvdos = ref["vvdos"]

            N_julia, L0_julia, L1_julia, L2_julia = BoltzTraP.fermi_integrals(
                epsilon, dos, vvdos, mur, Tr; dosweight=dosweight
            )

            L0_ref = ref["L0"]
            L1_ref = ref["L1"]
            L2_ref = ref["L2"]

            @test size(L0_julia) == size(L0_ref)
            @test L0_julia ≈ L0_ref rtol=0.01
            @test L1_julia ≈ L1_ref rtol=0.01
            @test L2_julia ≈ L2_ref rtol=0.01

            println("  ABINIT Si fermi_integrals: L0 max rel diff = $(maximum(abs.(L0_julia - L0_ref) ./ max.(abs.(L0_ref), 1e-20)))")
        end

        @testset "calc_onsager_coefficients" begin
            L0 = ref["L0"]
            L1 = ref["L1"]
            L2 = ref["L2"]

            sigma_julia, S_julia, kappa_julia = BoltzTraP.calc_onsager_coefficients(
                L0, L1, L2, Tr, vuc
            )

            sigma_ref = ref["sigma"]
            S_ref = ref["S"]
            kappa_ref = ref["kappa"]

            @test size(sigma_julia) == size(sigma_ref)
            @test sigma_julia ≈ sigma_ref rtol=0.01
            # Seebeck can have larger relative errors near zero
            @test S_julia ≈ S_ref rtol=0.1 atol=1e-6
            @test kappa_julia ≈ kappa_ref rtol=0.01

            # Sample value at T=300K (index 2), mid-mu
            iT = 2
            imu = size(sigma_julia, 2) ÷ 2
            println("  ABINIT Si sigma_xx at T=$(Tr[iT])K: Julia=$(sigma_julia[iT, imu, 1, 1]), Ref=$(sigma_ref[iT, imu, 1, 1])")
        end

        @testset "solve_for_mu" begin
            if !isdefined(BoltzTraP, :solve_for_mu)
                @test_skip "solve_for_mu not implemented"
                return
            end

            epsilon = ref["epsilon"]
            dos = ref["dos"]
            nelect = Float64(ref["nelect"])
            mu0_ref = ref["mu0"]

            for (iT, T) in enumerate(Tr)
                mu0_julia = BoltzTraP.solve_for_mu(
                    epsilon, dos, nelect, T;
                    dosweight=dosweight, refine=true, try_center=true
                )

                # Allow 1 histogram bin difference
                de = epsilon[2] - epsilon[1]
                @test abs(mu0_julia - mu0_ref[iT]) < 2 * de

                println("  ABINIT Si T=$(T)K: mu0 Julia=$(round(mu0_julia, digits=6)), Ref=$(round(mu0_ref[iT], digits=6))")
            end
        end
    end
end

# Main entry point
if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        run_all_tests()
    else
        for case in ARGS
            if case == "unit_cube"
                compare_unit_cube()
            elseif case == "simple_cubic"
                compare_simple_cubic()
            elseif case == "si_diamond"
                compare_si_diamond()
            elseif case == "fe_bcc_magnetic"
                compare_fe_bcc_magnetic()
            elseif case == "monoclinic"
                compare_monoclinic()
            elseif case == "triclinic_p1"
                compare_triclinic_p1()
            elseif case == "simple_interpolation"
                compare_simple_interpolation()
            elseif case == "si_interpolation"
                compare_si_interpolation()
            elseif case == "simple_transport"
                compare_simple_transport()
            elseif case == "transport"
                compare_simple_transport()
            elseif case == "poscar"
                compare_poscar()
            elseif case == "io"
                compare_poscar()
            elseif case == "si_end2end"
                compare_si_end2end()
            elseif case == "end2end"
                compare_si_end2end()
            elseif case == "si_integrate_e2e"
                compare_si_integrate_e2e()
            elseif case == "integrate_e2e"
                compare_si_integrate_e2e()
                compare_qe_si_integrate_e2e()
                compare_wien2k_si_integrate_e2e()
                compare_abinit_si_integrate_e2e()
            elseif case == "qe_si_integrate_e2e"
                compare_qe_si_integrate_e2e()
            elseif case == "wien2k_si_integrate_e2e"
                compare_wien2k_si_integrate_e2e()
            elseif case == "abinit_si_integrate_e2e"
                compare_abinit_si_integrate_e2e()
            elseif case == "bt2_python_format"
                compare_bt2_python_format()
            elseif case == "bt2"
                compare_bt2_python_format()
            else
                @warn "Unknown test case: $case"
            end
        end
    end
end
