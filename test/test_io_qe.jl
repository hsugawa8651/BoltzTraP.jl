# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP: read_qe_bands_gnu, read_qe_kpoints, read_qe_xml, read_qe_output, load_qe
using BoltzTraP: _generate_qe_tags

@testset "Quantum ESPRESSO I/O" begin
    @testset "read_qe_bands_gnu" begin
        # Create a minimal bands.dat.gnu file
        # Format: k_distance energy
        # Bands separated by blank lines
        bands_gnu_content = """
  0.0000000   -5.1234
  0.2500000   -4.5678
  0.5000000   -3.9012

  0.0000000    6.7890
  0.2500000    5.6789
  0.5000000    4.5678

  0.0000000    6.7890
  0.2500000    7.8901
  0.5000000    8.9012
"""
        mktempdir() do tmpdir
            bands_file = joinpath(tmpdir, "bands.dat.gnu")
            write(bands_file, bands_gnu_content)

            data = read_qe_bands_gnu(bands_file)

            # Check k-distances
            @test length(data.kdist) == 3
            @test data.kdist[1] ≈ 0.0
            @test data.kdist[2] ≈ 0.25
            @test data.kdist[3] ≈ 0.5

            # Check bands shape
            @test size(data.ebands) == (3, 3)  # nbands, nkpts

            # Check band energies
            @test data.ebands[1, 1] ≈ -5.1234
            @test data.ebands[2, 1] ≈ 6.789
            @test data.ebands[3, 2] ≈ 7.8901
        end
    end

    @testset "read_qe_kpoints" begin
        # Create a minimal bands.dat file
        # Format: nbnd nks on first line, then k-points and energies
        bands_dat_content = """           4           3
    0.000000    0.000000    0.000000
   -5.1234    6.7890    6.7890    6.7890
    0.500000    0.500000    0.000000
   -3.4567    4.5678    5.6789    7.8901
    0.500000    0.000000    0.500000
   -2.3456    3.4567    6.7890    8.9012
"""
        mktempdir() do tmpdir
            bands_file = joinpath(tmpdir, "bands.dat")
            write(bands_file, bands_dat_content)

            data = read_qe_kpoints(bands_file)

            # Check k-points
            @test size(data.kpoints) == (3, 3)
            @test data.kpoints[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.kpoints[:, 2] ≈ [0.5, 0.5, 0.0]
            @test data.kpoints[:, 3] ≈ [0.5, 0.0, 0.5]

            # Check bands
            @test size(data.ebands) == (4, 3)
            @test data.ebands[1, 1] ≈ -5.1234
            @test data.ebands[2, 2] ≈ 4.5678
            @test data.ebands[4, 3] ≈ 8.9012
        end

        # Test with multi-line energy format
        bands_dat_multiline = """           8           2
    0.000000    0.000000    0.000000
   -5.1234    6.7890    6.7890    6.7890
    7.8901    8.9012    9.0123   10.1234
    0.500000    0.500000    0.000000
   -3.4567    4.5678    5.6789    7.8901
    8.9012    9.0123   10.1234   11.2345
"""
        mktempdir() do tmpdir
            bands_file = joinpath(tmpdir, "bands.dat")
            write(bands_file, bands_dat_multiline)

            data = read_qe_kpoints(bands_file)

            @test size(data.ebands) == (8, 2)
            @test data.ebands[5, 1] ≈ 7.8901
            @test data.ebands[8, 2] ≈ 11.2345
        end
    end

    @testset "read_qe_xml" begin
        # Create a minimal data-file-schema.xml
        qe_xml_content = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
 <output>
  <atomic_structure>
   <cell>
    <a1>5.13 5.13 0.00</a1>
    <a2>0.00 5.13 5.13</a2>
    <a3>5.13 0.00 5.13</a3>
   </cell>
   <atomic_positions>
    <atom name="Si" index="1">0.0 0.0 0.0</atom>
    <atom name="Si" index="2">0.25 0.25 0.25</atom>
   </atomic_positions>
  </atomic_structure>
  <band_structure>
   <fermi_energy>0.2345</fermi_energy>
   <nks>2</nks>
   <nbnd>4</nbnd>
   <ks_energies>
    <k_point weight="0.5">0.0 0.0 0.0</k_point>
    <eigenvalues size="4">-0.1890  0.2500  0.2500  0.2500</eigenvalues>
   </ks_energies>
   <ks_energies>
    <k_point weight="0.5">0.5 0.5 0.0</k_point>
    <eigenvalues size="4">-0.1270  0.1680  0.2090  0.2900</eigenvalues>
   </ks_energies>
  </band_structure>
 </output>
</qes:espresso>
"""
        mktempdir() do tmpdir
            xml_file = joinpath(tmpdir, "data-file-schema.xml")
            write(xml_file, qe_xml_content)

            data = read_qe_xml(xml_file)

            # Check Fermi energy
            @test data.fermi ≈ 0.2345

            # Check lattice vectors (in Bohr)
            @test size(data.lattice) == (3, 3)
            @test data.lattice[:, 1] ≈ [5.13, 5.13, 0.0]
            @test data.lattice[:, 2] ≈ [0.0, 5.13, 5.13]
            @test data.lattice[:, 3] ≈ [5.13, 0.0, 5.13]

            # Check k-points
            @test size(data.kpoints) == (3, 2)
            @test data.kpoints[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.kpoints[:, 2] ≈ [0.5, 0.5, 0.0]

            # Check weights
            @test length(data.weights) == 2
            @test data.weights[1] ≈ 0.5
            @test data.weights[2] ≈ 0.5

            # Check eigenvalues (in Hartree)
            @test size(data.ebands) == (4, 2)
            @test data.ebands[1, 1] ≈ -0.189
            @test data.ebands[2, 2] ≈ 0.168

            # Check species
            @test length(data.species) == 2
            @test all(s -> s == "Si", data.species)

            # Check positions
            @test size(data.positions) == (3, 2)
        end
    end

    @testset "read_qe_output" begin
        # Create a minimal pw.x output file
        qe_output_content = """
     Program PWSCF v.7.0 starts on 10Jan2024 at 12: 0: 0

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.2600  a.u.

     number of k points=     3

     k = 0.0000 0.0000 0.0000 (   283 PWs)   bands (ev):

    -5.8094   6.2539   6.2539   6.2539

     k = 0.5000 0.5000 0.0000 (   290 PWs)   bands (ev):

    -3.4567   4.5678   5.6789   7.8901

     k = 0.5000 0.0000 0.5000 (   295 PWs)   bands (ev):

    -2.3456   3.4567   6.7890   8.9012

     the Fermi energy is     5.6789 ev

     Writing output data file ./pwscf.save/

     End of band structure calculation

"""
        mktempdir() do tmpdir
            output_file = joinpath(tmpdir, "pw.out")
            write(output_file, qe_output_content)

            data = read_qe_output(output_file)

            # Check k-points
            @test size(data.kpoints) == (3, 3)
            @test data.kpoints[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.kpoints[:, 2] ≈ [0.5, 0.5, 0.0]
            @test data.kpoints[:, 3] ≈ [0.5, 0.0, 0.5]

            # Check bands
            @test size(data.ebands) == (4, 3)
            @test data.ebands[1, 1] ≈ -5.8094
            @test data.ebands[2, 1] ≈ 6.2539
            @test data.ebands[1, 2] ≈ -3.4567
            @test data.ebands[4, 3] ≈ 8.9012
        end

        # Test with different output format (occupation numbers visible)
        qe_output_occ = """
     k = 0.0000 0.0000 0.0000 (   283 PWs)   bands (ev):

    -5.8094   6.2539   6.2539   6.2539   8.1234   9.2345

     highest occupied level (ev):     6.2539

"""
        mktempdir() do tmpdir
            output_file = joinpath(tmpdir, "pw.out")
            write(output_file, qe_output_occ)

            data = read_qe_output(output_file)

            @test size(data.kpoints) == (3, 1)
            @test size(data.ebands) == (6, 1)
            @test data.ebands[6, 1] ≈ 9.2345
        end
    end

    @testset "_generate_qe_tags" begin
        # Test tags generation following Python BoltzTraP2 behavior
        # Same extended species get same tag (0-based)
        @test _generate_qe_tags(["Fe1", "Fe1", "Fe2", "Fe2"]) == [0, 0, 1, 1]
        @test _generate_qe_tags(["Si", "Si"]) == [0, 0]
        @test _generate_qe_tags(["Cr1", "Cr2", "Cr1", "Cr2"]) == [0, 1, 0, 1]
        @test _generate_qe_tags(["Fe_up", "Fe_down"]) == [0, 1]
    end

    @testset "load_qe with extended species (magmom generation)" begin
        # Synthetic QE XML with extended species (e.g., Fe1, Fe2 for antiferromagnetic)
        qe_xml_extended = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
 <output>
  <atomic_structure>
   <cell>
    <a1>5.40 0.00 0.00</a1>
    <a2>0.00 5.40 0.00</a2>
    <a3>0.00 0.00 5.40</a3>
   </cell>
   <atomic_positions>
    <atom name="Fe1" index="1">0.0 0.0 0.0</atom>
    <atom name="Fe1" index="2">2.70 2.70 0.0</atom>
    <atom name="Fe2" index="3">2.70 0.0 2.70</atom>
    <atom name="Fe2" index="4">0.0 2.70 2.70</atom>
   </atomic_positions>
  </atomic_structure>
  <band_structure>
   <lsda>false</lsda>
   <fermi_energy>0.2500</fermi_energy>
   <nelec>32.0</nelec>
   <nks>1</nks>
   <nbnd>4</nbnd>
   <ks_energies>
    <k_point weight="1.0">0.0 0.0 0.0</k_point>
    <eigenvalues size="4">-0.10 0.20 0.20 0.30</eigenvalues>
    <occupations size="4">1.0 1.0 1.0 1.0</occupations>
   </ks_energies>
  </band_structure>
 </output>
</qes:espresso>
"""
        mktempdir() do tmpdir
            xml_file = joinpath(tmpdir, "data-file-schema.xml")
            write(xml_file, qe_xml_extended)

            # Test read_qe_xml returns extended_species
            raw = read_qe_xml(xml_file)
            @test raw.extended_species == ["Fe1", "Fe1", "Fe2", "Fe2"]
            @test raw.species == ["Fe", "Fe", "Fe", "Fe"]

            # Test load_qe generates dummy magmom
            data = load_qe(tmpdir)
            @test !isnothing(data.magmom)
            @test length(data.magmom) == 4
            # magmom = tags + 1.0 = [0, 0, 1, 1] + 1.0 = [1.0, 1.0, 2.0, 2.0]
            @test data.magmom ≈ [1.0, 1.0, 2.0, 2.0]
        end

        # Test with simple species (no extended names)
        qe_xml_simple = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">
 <output>
  <atomic_structure>
   <cell>
    <a1>5.13 5.13 0.00</a1>
    <a2>0.00 5.13 5.13</a2>
    <a3>5.13 0.00 5.13</a3>
   </cell>
   <atomic_positions>
    <atom name="Si" index="1">0.0 0.0 0.0</atom>
    <atom name="Si" index="2">2.565 2.565 2.565</atom>
   </atomic_positions>
  </atomic_structure>
  <band_structure>
   <lsda>false</lsda>
   <fermi_energy>0.2345</fermi_energy>
   <nelec>8.0</nelec>
   <nks>1</nks>
   <nbnd>4</nbnd>
   <ks_energies>
    <k_point weight="1.0">0.0 0.0 0.0</k_point>
    <eigenvalues size="4">-0.189 0.25 0.25 0.25</eigenvalues>
    <occupations size="4">1.0 1.0 1.0 1.0</occupations>
   </ks_energies>
  </band_structure>
 </output>
</qes:espresso>
"""
        mktempdir() do tmpdir
            xml_file = joinpath(tmpdir, "data-file-schema.xml")
            write(xml_file, qe_xml_simple)

            data = load_qe(tmpdir)
            # For simple species, tags = [0, 0], magmom = [1.0, 1.0]
            @test !isnothing(data.magmom)
            @test data.magmom ≈ [1.0, 1.0]
        end
    end

end
