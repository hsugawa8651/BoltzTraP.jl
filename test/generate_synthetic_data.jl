# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Pure-Julia synthetic DFT data generator for loader testing.
# Users can verify loader functionality without Python or external tools.
#
# Usage:
#   julia --project=. test/generate_synthetic_data.jl
#   julia --project=. test/generate_synthetic_data.jl --verify

using LinearAlgebra
using Printf

# Conditionally load BoltzTraP for verification
const BOLTZTRAP_AVAILABLE = try
    @eval using BoltzTraP
    true
catch
    false
end

# Conditionally load NCDatasets for ABINIT format
const NCDATASETS_AVAILABLE = try
    @eval using NCDatasets
    true
catch
    false
end

# Constants
const BOHR_TO_ANG = 0.529177210903
const HA_TO_EV = 27.211386245988
const HA_TO_RY = 2.0

"""
Build lattice vectors from cell parameters (a, b, c, α, β, γ).

Returns 3×3 matrix where columns are lattice vectors in Bohr.
"""
function build_lattice(a, b, c, alpha_deg, beta_deg, gamma_deg)
    α = deg2rad(alpha_deg)
    β = deg2rad(beta_deg)
    γ = deg2rad(gamma_deg)

    # Standard crystallographic convention
    ax = a
    bx = b * cos(γ)
    by = b * sin(γ)
    cx = c * cos(β)
    cy = c * (cos(α) - cos(β) * cos(γ)) / sin(γ)
    cz = sqrt(c^2 - cx^2 - cy^2)

    return [ax 0.0 0.0; bx by 0.0; cx cy cz]'  # Columns are vectors
end

"""
Generate Monkhorst-Pack k-point grid.
"""
function monkhorst_pack(n1, n2, n3)
    kpoints = Matrix{Float64}(undef, 3, n1 * n2 * n3)
    idx = 1
    for i in 0:(n1-1), j in 0:(n2-1), k in 0:(n3-1)
        kpoints[1, idx] = (2i - n1 + 1) / (2n1)
        kpoints[2, idx] = (2j - n2 + 1) / (2n2)
        kpoints[3, idx] = (2k - n3 + 1) / (2n3)
        idx += 1
    end
    return kpoints
end

"""
Calculate free-electron band energies: E(k) = |k + G|² / 2 (atomic units).
"""
function free_electron_bands(kpoints, rec_lattice, nbands)
    nk = size(kpoints, 2)
    ebands = zeros(nbands, nk)

    # Generate G-vectors (reciprocal lattice translations)
    gvecs = Vector{Vector{Int}}()
    for i in -2:2, j in -2:2, k in -2:2
        push!(gvecs, [i, j, k])
    end

    for ik in 1:nk
        k_frac = kpoints[:, ik]
        k_cart = rec_lattice * k_frac

        # Calculate |k + G|² for all G-vectors
        energies = Float64[]
        for g in gvecs
            kg = k_cart + rec_lattice * g
            push!(energies, dot(kg, kg) / 2)
        end
        sort!(energies)

        # Take lowest nbands
        ebands[:, ik] = energies[1:nbands]
    end

    return ebands
end

# ============================================================================
# GENE Format Writer
# ============================================================================

function write_gene(dir, lattice, positions, species, kpoints, ebands, fermi)
    mkpath(dir)
    natoms = size(positions, 2)
    nk, nbands = size(kpoints, 2), size(ebands, 1)

    # .structure file
    # GENE format: Element + Cartesian coordinates (Bohr)
    open(joinpath(dir, "Synthetic.structure"), "w") do f
        println(f, "Synthetic low-symmetry structure")
        for i in 1:3
            @printf(f, "%20.14f %20.14f %20.14f\n", lattice[1, i], lattice[2, i], lattice[3, i])
        end
        println(f, natoms)
        # Convert fractional to Cartesian coordinates
        for i in 1:natoms
            cart = lattice * positions[:, i]
            @printf(f, "%s %20.14f %20.14f %20.14f\n",
                species[i], cart[1], cart[2], cart[3])
        end
    end

    # .energy file (Rydberg units)
    # GENE format: Title, nk nspin efermi(Ry), then k-point blocks
    open(joinpath(dir, "Synthetic.energy"), "w") do f
        println(f, "Synthetic GENE energy file")
        @printf(f, "%6d %6d %20.10e\n", nk, 1, fermi * HA_TO_RY)  # nk, nspin, efermi
        for ik in 1:nk
            @printf(f, "%20.10e %20.10e %20.10e %6d\n",
                kpoints[1, ik], kpoints[2, ik], kpoints[3, ik], nbands)
            for ib in 1:nbands
                @printf(f, "%20.10e\n", ebands[ib, ik] * HA_TO_RY)
            end
        end
    end
end

# ============================================================================
# VASP Format Writer
# ============================================================================

function write_vasp(dir, lattice, positions, species, kpoints, ebands, fermi, nelect)
    mkpath(dir)
    lattice_ang = lattice * BOHR_TO_ANG
    nk, nbands = size(kpoints, 2), size(ebands, 1)
    natoms = size(positions, 2)

    open(joinpath(dir, "vasprun.xml"), "w") do f
        println(f, """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
<generator>
 <i name="program" type="string">synthetic</i>
</generator>
<parameters>
 <i name="NELECT">      $(Float64(nelect))</i>
</parameters>
<atominfo>
 <atoms>$(natoms)</atoms>
 <types>$(length(unique(species)))</types>
 <array name="atoms">
  <set>""")
        for sp in species
            println(f, "   <rc><c>   $(sp)</c><c>     1</c></rc>")
        end
        println(f, """  </set>
 </array>
</atominfo>
<structure name="finalpos">
 <crystal>
  <varray name="basis">""")
        for i in 1:3
            @printf(f, "   <v>%20.14f %20.14f %20.14f</v>\n",
                lattice_ang[1, i], lattice_ang[2, i], lattice_ang[3, i])
        end
        println(f, """  </varray>
  <varray name="rec_basis">""")
        rec = 2π * inv(lattice_ang)'
        for i in 1:3
            @printf(f, "   <v>%20.14f %20.14f %20.14f</v>\n", rec[1, i], rec[2, i], rec[3, i])
        end
        println(f, """  </varray>
 </crystal>
 <varray name="positions">""")
        for i in 1:natoms
            @printf(f, "  <v>%20.14f %20.14f %20.14f</v>\n",
                positions[1, i], positions[2, i], positions[3, i])
        end
        println(f, """ </varray>
</structure>
<kpoints>
 <varray name="kpointlist">""")
        for ik in 1:nk
            @printf(f, "  <v>%20.14f %20.14f %20.14f</v>\n",
                kpoints[1, ik], kpoints[2, ik], kpoints[3, ik])
        end
        println(f, """ </varray>
 <varray name="weights">""")
        w = 1.0 / nk
        for _ in 1:nk
            @printf(f, "  <v>%20.14f</v>\n", w)
        end
        println(f, """ </varray>
</kpoints>
<calculation>
 <dos>
  <i name="efermi">$(fermi * HA_TO_EV)</i>
 </dos>
 <eigenvalues>
  <array>
   <set>
    <set comment="spin 1">""")
        for ik in 1:nk
            println(f, "     <set comment=\"kpoint $(ik)\">")
            for ib in 1:nbands
                occ = ebands[ib, ik] < fermi ? 2.0 : 0.0
                @printf(f, "      <r>%20.14f %20.14f</r>\n", ebands[ib, ik] * HA_TO_EV, occ)
            end
            println(f, "     </set>")
        end
        println(f, """    </set>
   </set>
  </array>
 </eigenvalues>
</calculation>
</modeling>""")
    end
end

# ============================================================================
# Quantum ESPRESSO Format Writer
# ============================================================================

function write_qe(dir, lattice, positions, species, kpoints, ebands, fermi)
    outdir = joinpath(dir, "out")
    mkpath(outdir)
    lattice_ang = lattice * BOHR_TO_ANG
    nk, nbands = size(kpoints, 2), size(ebands, 1)
    natoms = size(positions, 2)

    open(joinpath(outdir, "synthetic.xml"), "w") do f
        println(f, """<?xml version="1.0" encoding="UTF-8"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
<output>
<atomic_structure nat="$(natoms)">
<cell>
<a1>$(@sprintf("%20.14f %20.14f %20.14f", lattice_ang[1,1], lattice_ang[2,1], lattice_ang[3,1]))</a1>
<a2>$(@sprintf("%20.14f %20.14f %20.14f", lattice_ang[1,2], lattice_ang[2,2], lattice_ang[3,2]))</a2>
<a3>$(@sprintf("%20.14f %20.14f %20.14f", lattice_ang[1,3], lattice_ang[2,3], lattice_ang[3,3]))</a3>
</cell>
<atomic_positions>""")
        for i in 1:natoms
            cart = lattice_ang * positions[:, i]
            println(f, "<atom name=\"$(species[i])\" index=\"$(i)\">$(@sprintf("%20.14f %20.14f %20.14f", cart...))</atom>")
        end
        println(f, """</atomic_positions>
</atomic_structure>
<band_structure>
<fermi_energy>$(fermi * HA_TO_EV / HA_TO_EV)</fermi_energy>
<nks>$(nk)</nks>
<nbnd>$(nbands)</nbnd>""")
        w = 1.0 / nk
        for ik in 1:nk
            println(f, "<ks_energies>")
            @printf(f, "<k_point weight=\"%.12f\">%20.14f %20.14f %20.14f</k_point>\n",
                w, kpoints[1, ik], kpoints[2, ik], kpoints[3, ik])
            print(f, "<eigenvalues size=\"$(nbands)\">")
            for ib in 1:nbands
                @printf(f, " %.14e", ebands[ib, ik])  # Hartree
            end
            println(f, "</eigenvalues>")
            print(f, "<occupations size=\"$(nbands)\">")
            for ib in 1:nbands
                occ = ebands[ib, ik] < fermi ? 2.0 : 0.0
                @printf(f, " %.14e", occ)
            end
            println(f, "</occupations>")
            println(f, "</ks_energies>")
        end
        println(f, """</band_structure>
</output>
</qes:espresso>""")
    end
end

# ============================================================================
# Wien2k Format Writer
# ============================================================================

function write_wien2k(dir, lattice, positions, species, kpoints, ebands, fermi)
    mkpath(dir)
    nk, nbands = size(kpoints, 2), size(ebands, 1)
    natoms = size(positions, 2)

    # Calculate cell parameters
    a_vec, b_vec, c_vec = eachcol(lattice)
    a, b, c = norm(a_vec), norm(b_vec), norm(c_vec)
    alpha = acosd(dot(b_vec, c_vec) / (b * c))
    beta = acosd(dot(a_vec, c_vec) / (a * c))
    gamma = acosd(dot(a_vec, b_vec) / (a * b))

    # .struct file
    # Wien2k fixed-width format: X at cols 13-22, Y at 26-35, Z at 39-48
    open(joinpath(dir, "Synthetic.struct"), "w") do f
        println(f, "Synthetic low-symmetry structure")
        @printf(f, "P   LATTICE,NONEQUIV.ATOMS:%3d\n", natoms)
        println(f, "MODE OF CALC=RELA unit=bohr")
        @printf(f, "%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n", a, b, c, alpha, beta, gamma)
        for i in 1:natoms
            # Position line: "ATOM  -N: X=xxxxxxxxxx Y=xxxxxxxxxx Z=xxxxxxxxxx"
            # Columns: X at 13-22, Y at 26-35, Z at 39-48 (1-indexed)
            @printf(f, "ATOM  -%d: X=%10.8f Y=%10.8f Z=%10.8f\n",
                i, positions[1, i], positions[2, i], positions[3, i])
            println(f, "          MULT= 1          ISPLIT= 8")
            @printf(f, "%-2s         NPT=  781  R0=0.00010000 RMT=   2.00000   Z:%5.1f\n",
                species[i], i == 1 ? 6.0 : 7.0)
            println(f, "LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000")
            println(f, "                     0.0000000 1.0000000 0.0000000")
            println(f, "                     0.0000000 0.0000000 1.0000000")
        end
        println(f, "   1      NUMBER OF SYMMETRY OPERATIONS")
        println(f, " 1 0 0 0.00000000")
        println(f, " 0 1 0 0.00000000")
        println(f, " 0 0 1 0.00000000")
        println(f, "       1")
    end

    # .energy file (Rydberg units)
    open(joinpath(dir, "Synthetic.energy"), "w") do f
        for ik in 1:nk
            @printf(f, "%19.12E%19.12E%19.12E%10d%6d%6d  1.0   \n",
                kpoints[1, ik], kpoints[2, ik], kpoints[3, ik], ik, 1, nbands)
            for ib in 1:nbands
                @printf(f, "%4d%19.12E\n", ib, ebands[ib, ik] * HA_TO_RY)
            end
        end
    end

    # .scf file (Fermi energy)
    # Format: ":FER  : F E R M I - ENERGY(...)=   value" with value at col 39+
    open(joinpath(dir, "Synthetic.scf"), "w") do f
        @printf(f, ":FER  : F E R M I - ENERGY(TETRAH.M.)=   %.10f\n", fermi * HA_TO_RY)
    end
end

# ============================================================================
# ABINIT Format Writer (requires NCDatasets.jl)
# ============================================================================

function write_abinit(dir, lattice, positions, species, kpoints, ebands, fermi, nelect)
    if !NCDATASETS_AVAILABLE
        return false
    end

    mkpath(dir)
    nk, nbands = size(kpoints, 2), size(ebands, 1)
    natoms = size(positions, 2)

    # Count unique species types
    unique_species = unique(species)
    ntypat = length(unique_species)

    gsr_file = joinpath(dir, "synthetic_GSR.nc")

    NCDatasets.NCDataset(gsr_file, "c") do ds
        # Dimensions
        NCDatasets.defDim(ds, "number_of_cartesian_directions", 3)
        NCDatasets.defDim(ds, "number_of_reduced_dimensions", 3)
        NCDatasets.defDim(ds, "three", 3)
        NCDatasets.defDim(ds, "number_of_atoms", natoms)
        NCDatasets.defDim(ds, "number_of_kpoints", nk)
        NCDatasets.defDim(ds, "number_of_spins", 1)
        NCDatasets.defDim(ds, "max_number_of_states", nbands)
        NCDatasets.defDim(ds, "number_of_atom_species", ntypat)
        NCDatasets.defDim(ds, "symbol_length", 80)

        # Lattice vectors (Bohr) - (3, 3)
        lat_var = NCDatasets.defVar(ds, "primitive_vectors", Float64,
            ("number_of_cartesian_directions", "number_of_reduced_dimensions"))
        lat_var[:, :] = lattice

        # Atomic positions (fractional) - (3, natom) per ABINIT format
        pos_var = NCDatasets.defVar(ds, "reduced_atom_positions", Float64,
            ("three", "number_of_atoms"))
        pos_var[:, :] = positions

        # Species index per atom (1-based)
        species_idx = [findfirst(==(sp), unique_species) for sp in species]
        species_var = NCDatasets.defVar(ds, "atom_species", Int32, ("number_of_atoms",))
        species_var[:] = Int32.(species_idx)

        # Species names as Char array (80, ntypat)
        names_var = NCDatasets.defVar(ds, "atom_species_names", Char,
            ("symbol_length", "number_of_atom_species"))
        for (i, sp) in enumerate(unique_species)
            chars = fill('\0', 80)
            for (j, c) in enumerate(sp)
                chars[j] = c
            end
            names_var[:, i] = chars
        end

        # K-points (fractional) - (3, nk) per ABINIT format
        kpt_var = NCDatasets.defVar(ds, "reduced_coordinates_of_kpoints", Float64,
            ("three", "number_of_kpoints"))
        kpt_var[:, :] = kpoints

        # K-point weights
        wt_var = NCDatasets.defVar(ds, "kpoint_weights", Float64, ("number_of_kpoints",))
        wt_var[:] = ones(nk) / nk

        # Eigenvalues (Hartree) - (nbands, nk, nspin) per ABINIT format
        eig_var = NCDatasets.defVar(ds, "eigenvalues", Float64,
            ("max_number_of_states", "number_of_kpoints", "number_of_spins"))
        eig_var[:, :, 1] = ebands

        # Fermi energy (scalar)
        fermi_var = NCDatasets.defVar(ds, "fermie", Float64, ())
        fermi_var[] = fermi

        # Number of electrons (scalar)
        nelect_var = NCDatasets.defVar(ds, "nelect", Float64, ())
        nelect_var[] = Float64(nelect)
    end

    return true
end

# ============================================================================
# Main Generator
# ============================================================================

function generate_all(output_base::String)
    # Cell definitions: (name, a, b, c, α, β, γ) in Bohr and degrees
    cells = [
        ("monoclinic", 5.1, 4.8, 6.2, 90.0, 100.0, 90.0),
        ("triclinic", 5.1, 4.8, 6.2, 85.0, 80.0, 75.0),
    ]

    # Atoms at general positions (breaks all symmetry)
    frac_positions = [0.12 0.67; 0.23 0.78; 0.34 0.56]
    species = ["C", "N"]

    # Band parameters
    nk_each = 4
    nbands = 8
    nelect = 8
    fermi = 0.3  # Hartree

    for (name, a, b, c, α, β, γ) in cells
        println("Generating $name structure...")

        lattice = build_lattice(a, b, c, α, β, γ)
        rec_lattice = 2π * inv(lattice)'
        kpoints = monkhorst_pack(nk_each, nk_each, nk_each)
        ebands = free_electron_bands(kpoints, rec_lattice, nbands)

        # Write all formats
        write_gene(joinpath(output_base, "Synthetic.GENE.$name"),
            lattice, frac_positions, species, kpoints, ebands, fermi)
        println("  - GENE: done")

        write_vasp(joinpath(output_base, "Synthetic.vasp.$name"),
            lattice, frac_positions, species, kpoints, ebands, fermi, nelect)
        println("  - VASP: done")

        write_qe(joinpath(output_base, "Synthetic.ESPRESSO.$name"),
            lattice, frac_positions, species, kpoints, ebands, fermi)
        println("  - QE: done")

        write_wien2k(joinpath(output_base, "Synthetic.W2K.$name"),
            lattice, frac_positions, species, kpoints, ebands, fermi)
        println("  - Wien2k: done")

        # ABINIT requires NCDatasets.jl
        if write_abinit(joinpath(output_base, "Synthetic.abinit.$name"),
                lattice, frac_positions, species, kpoints, ebands, fermi, nelect)
            println("  - ABINIT: done")
        else
            println("  - ABINIT: skipped (NCDatasets.jl not available)")
        end
    end

    println("\nGenerated files in: $output_base")
end

# ============================================================================
# Verification
# ============================================================================

function verify_loaders(data_base::String)
    if !BOLTZTRAP_AVAILABLE
        println("BoltzTraP not available. Run with: julia --project=. test/generate_synthetic_data.jl --verify")
        return false
    end

    println("Verifying loaders...")
    println("=" ^ 60)

    cells = [
        ("monoclinic", (α=90.0, β=100.0, γ=90.0)),
        ("triclinic", (α=85.0, β=80.0, γ=75.0)),
    ]

    function calc_angles(lattice)
        a, b, c = norm.(eachcol(lattice))
        α = acosd(dot(lattice[:, 2], lattice[:, 3]) / (b * c))
        β = acosd(dot(lattice[:, 1], lattice[:, 3]) / (a * c))
        γ = acosd(dot(lattice[:, 1], lattice[:, 2]) / (a * b))
        return (α=α, β=β, γ=γ)
    end

    loaders = [
        ("GENE", "Synthetic.GENE", BoltzTraP.load_gene, false),
        ("VASP", "Synthetic.vasp", BoltzTraP.load_vasp, false),
        ("QE", "Synthetic.ESPRESSO", BoltzTraP.load_qe, true),
        ("Wien2k", "Synthetic.W2K", BoltzTraP.load_wien2k, false),
        ("ABINIT", "Synthetic.abinit", BoltzTraP.load_abinit, false),
    ]

    all_pass = true

    for (cell_name, expected_angles) in cells
        println("\n[$cell_name] Expected: α=$(expected_angles.α)°, β=$(expected_angles.β)°, γ=$(expected_angles.γ)°")
        println("-" ^ 60)

        for (loader_name, prefix, loader_fn, use_out) in loaders
            dir = joinpath(data_base, "$(prefix).$(cell_name)")
            if use_out
                dir = joinpath(dir, "out")
            end

            if !isdir(dir)
                println("  $loader_name: SKIP (directory not found)")
                continue
            end

            try
                data = loader_fn(dir)
                angles = calc_angles(data.lattice)

                α_ok = isapprox(angles.α, expected_angles.α, atol=0.5)
                β_ok = isapprox(angles.β, expected_angles.β, atol=0.5)
                γ_ok = isapprox(angles.γ, expected_angles.γ, atol=0.5)

                if α_ok && β_ok && γ_ok
                    @printf("  %-8s: PASS (α=%.1f°, β=%.1f°, γ=%.1f°)\n",
                        loader_name, angles.α, angles.β, angles.γ)
                else
                    @printf("  %-8s: FAIL (α=%.1f°, β=%.1f°, γ=%.1f°)\n",
                        loader_name, angles.α, angles.β, angles.γ)
                    all_pass = false
                end
            catch e
                println("  $loader_name: ERROR - $e")
                all_pass = false
            end
        end
    end

    println("\n" * "=" ^ 60)
    println(all_pass ? "All tests PASSED" : "Some tests FAILED")
    return all_pass
end

# ============================================================================
# Entry Point
# ============================================================================

function main()
    # Default output directory
    script_dir = @__DIR__
    project_dir = dirname(script_dir)
    output_base = joinpath(project_dir, "test", "synthetic_data")

    if "--verify" in ARGS
        # Verify existing data
        if !isdir(output_base)
            println("No synthetic data found. Run without --verify first.")
            return
        end
        verify_loaders(output_base)
    else
        # Generate synthetic data
        generate_all(output_base)
        println("\nTo verify loaders, run:")
        println("  julia --project=. test/generate_synthetic_data.jl --verify")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
