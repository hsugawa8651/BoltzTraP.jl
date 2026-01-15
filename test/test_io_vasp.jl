# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP: read_poscar, read_vasprun, read_eigenval

@testset "VASP I/O" begin

    @testset "read_poscar" begin
        # Create a temporary POSCAR file (Si diamond)
        poscar_content = """
Si2 diamond
5.43
  0.5  0.5  0.0
  0.0  0.5  0.5
  0.5  0.0  0.5
Si
2
Direct
  0.00  0.00  0.00
  0.25  0.25  0.25
"""
        mktempdir() do tmpdir
            poscar_file = joinpath(tmpdir, "POSCAR")
            write(poscar_file, poscar_content)

            data = read_poscar(poscar_file)

            @test data.comment == "Si2 diamond"
            @test data.species == ["Si"]
            @test data.counts == [2]
            @test size(data.lattice) == (3, 3)
            @test size(data.positions) == (3, 2)

            # Check lattice vectors (scaled by 5.43)
            @test data.lattice[:, 1] ≈ [0.5, 0.5, 0.0] * 5.43
            @test data.lattice[:, 2] ≈ [0.0, 0.5, 0.5] * 5.43
            @test data.lattice[:, 3] ≈ [0.5, 0.0, 0.5] * 5.43
        end

        # Test VASP 4 format (no element names)
        poscar_v4 = """
Si2 diamond
5.43
  0.5  0.5  0.0
  0.0  0.5  0.5
  0.5  0.0  0.5
2
Direct
  0.00  0.00  0.00
  0.25  0.25  0.25
"""
        mktempdir() do tmpdir
            poscar_file = joinpath(tmpdir, "POSCAR")
            write(poscar_file, poscar_v4)

            data = read_poscar(poscar_file)

            @test isnothing(data.species)
            @test data.counts == [2]
        end

        # Test Cartesian coordinates
        poscar_cart = """
Si2 diamond
5.43
  0.5  0.5  0.0
  0.0  0.5  0.5
  0.5  0.0  0.5
Si
2
Cartesian
  0.00  0.00  0.00
  1.3575  1.3575  1.3575
"""
        mktempdir() do tmpdir
            poscar_file = joinpath(tmpdir, "POSCAR")
            write(poscar_file, poscar_cart)

            data = read_poscar(poscar_file)

            # Cartesian coordinates should be kept as-is
            @test data.positions[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.positions[:, 2] ≈ [1.3575, 1.3575, 1.3575]
        end

        # Test with selective dynamics
        poscar_sd = """
Si2 diamond
5.43
  0.5  0.5  0.0
  0.0  0.5  0.5
  0.5  0.0  0.5
Si
2
Selective dynamics
Direct
  0.00  0.00  0.00 T T T
  0.25  0.25  0.25 F F F
"""
        mktempdir() do tmpdir
            poscar_file = joinpath(tmpdir, "POSCAR")
            write(poscar_file, poscar_sd)

            data = read_poscar(poscar_file)

            @test data.counts == [2]
            @test size(data.positions) == (3, 2)
        end
    end

    @testset "read_vasprun" begin
        # Create a minimal vasprun.xml
        vasprun_content = """<?xml version="1.0"?>
<modeling>
 <parameters>
  <separator name="electronic" >
   <i name="NELECT">     8.00000000</i>
  </separator>
 </parameters>
 <calculation>
  <structure name="finalpos">
   <crystal>
    <varray name="basis">
     <v>   2.7150000   2.7150000   0.0000000 </v>
     <v>   0.0000000   2.7150000   2.7150000 </v>
     <v>   2.7150000   0.0000000   2.7150000 </v>
    </varray>
    <varray name="rec_basis">
     <v>   0.1841620   0.1841620  -0.1841620 </v>
     <v>  -0.1841620   0.1841620   0.1841620 </v>
     <v>   0.1841620  -0.1841620   0.1841620 </v>
    </varray>
   </crystal>
   <varray name="positions">
    <v>   0.0000000   0.0000000   0.0000000 </v>
    <v>   0.2500000   0.2500000   0.2500000 </v>
   </varray>
  </structure>
  <dos>
   <i name="efermi">     5.67890000</i>
  </dos>
  <eigenvalues>
   <array>
    <set>
     <set comment="spin 1">
      <set comment="kpoint 1">
       <r>   -5.1234    1.00000000 </r>
       <r>    6.7890    1.00000000 </r>
       <r>    6.7890    1.00000000 </r>
       <r>    6.7890    1.00000000 </r>
      </set>
      <set comment="kpoint 2">
       <r>   -3.4567    1.00000000 </r>
       <r>    4.5678    1.00000000 </r>
       <r>    5.6789    0.50000000 </r>
       <r>    7.8901    0.00000000 </r>
      </set>
     </set>
    </set>
   </array>
  </eigenvalues>
 </calculation>
 <atominfo>
  <array name="atoms">
   <set>
    <rc><c>Si </c><c>   1</c></rc>
    <rc><c>Si </c><c>   1</c></rc>
   </set>
  </array>
 </atominfo>
 <kpoints>
  <varray name="kpointlist">
   <v>   0.0000000   0.0000000   0.0000000 </v>
   <v>   0.5000000   0.5000000   0.0000000 </v>
  </varray>
  <varray name="weights">
   <v>   0.5000000 </v>
   <v>   0.5000000 </v>
  </varray>
 </kpoints>
</modeling>
"""
        mktempdir() do tmpdir
            vasprun_file = joinpath(tmpdir, "vasprun.xml")
            write(vasprun_file, vasprun_content)

            data = read_vasprun(vasprun_file)

            # Check Fermi energy
            @test data.fermi ≈ 5.6789

            # Check number of electrons
            @test data.nelect ≈ 8.0

            # Check lattice vectors
            @test size(data.lattice) == (3, 3)
            @test data.lattice[:, 1] ≈ [2.715, 2.715, 0.0]

            # Check reciprocal lattice
            @test size(data.rec_lattice) == (3, 3)

            # Check positions
            @test size(data.positions) == (3, 2)
            @test data.positions[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.positions[:, 2] ≈ [0.25, 0.25, 0.25]

            # Check species
            @test length(data.species) == 2
            @test all(s -> strip(s) == "Si", data.species)

            # Check k-points
            @test size(data.kpoints) == (3, 2)
            @test data.kpoints[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.kpoints[:, 2] ≈ [0.5, 0.5, 0.0]

            # Check weights
            @test length(data.weights) == 2
            @test all(w -> w ≈ 0.5, data.weights)

            # Check eigenvalues
            @test size(data.ebands) == (4, 2, 1)  # nbands, nkpts, nspin
            @test data.ebands[1, 1, 1] ≈ -5.1234
            @test data.ebands[2, 1, 1] ≈ 6.789

            # Check occupations
            @test size(data.occupations) == (4, 2, 1)
            @test data.occupations[1, 1, 1] ≈ 1.0
            @test data.occupations[3, 2, 1] ≈ 0.5
        end
    end

    @testset "read_eigenval" begin
        # Create a minimal EIGENVAL file
        eigenval_content = """   2   2   1   1
  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0
  CAR
   8   2   4

   0.0000000E+00  0.0000000E+00  0.0000000E+00  0.5000000E+00
    1      -5.1234    1.00000
    2       6.7890    1.00000
    3       6.7890    1.00000
    4       6.7890    1.00000

   0.5000000E+00  0.5000000E+00  0.0000000E+00  0.5000000E+00
    1      -3.4567    1.00000
    2       4.5678    1.00000
    3       5.6789    0.50000
    4       7.8901    0.00000

"""
        mktempdir() do tmpdir
            eigenval_file = joinpath(tmpdir, "EIGENVAL")
            write(eigenval_file, eigenval_content)

            data = read_eigenval(eigenval_file)

            @test data.nelect == 8
            @test data.nkpts == 2
            @test data.nbands == 4

            # Check k-points
            @test size(data.kpoints) == (3, 2)
            @test data.kpoints[:, 1] ≈ [0.0, 0.0, 0.0]
            @test data.kpoints[:, 2] ≈ [0.5, 0.5, 0.0]

            # Check weights
            @test data.weights[1] ≈ 0.5
            @test data.weights[2] ≈ 0.5

            # Check eigenvalues
            @test size(data.ebands) == (4, 2, 1)
            @test data.ebands[1, 1, 1] ≈ -5.1234
            @test data.ebands[2, 2, 1] ≈ 4.5678

            # Check occupations
            @test data.occupations[1, 1, 1] ≈ 1.0
            @test data.occupations[3, 2, 1] ≈ 0.5
        end
    end

end
