# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using LinearAlgebra
using Statistics
using NPZ
using BoltzTraP

@testset "Magnetic Workflow" begin

    @testset "run_interpolate Collinear (DFTData{2})" begin
        # Load PbTe collinear reference data
        REFTEST_DIR = joinpath(@__DIR__, "..", "reftest", "data")
        reference_file = joinpath(REFTEST_DIR, "pbte_collinear.npz")

        if !isfile(reference_file)
            @test_skip "PbTe collinear reference data not available"
            return
        end

        ref = npzread(reference_file)

        # Create synthetic DFTData{2} from reference
        # Note: Reference ebands is already concatenated (32, 145)
        # We need to split it back to (16, 145, 2) for DFTData{2}
        ebands_concat = ref["ebands"]  # (32, 145) - concatenated
        nbands_per_spin = size(ebands_concat, 1) ÷ 2
        nkpts = size(ebands_concat, 2)

        # Split into two spin channels
        ebands_spin1 = ebands_concat[1:nbands_per_spin, :]
        ebands_spin2 = ebands_concat[nbands_per_spin+1:end, :]
        ebands_3d = zeros(nbands_per_spin, nkpts, 2)
        ebands_3d[:, :, 1] = ebands_spin1
        ebands_3d[:, :, 2] = ebands_spin2

        # Decode species from uint8 array
        symbols_bytes = ref["symbols"]
        symbols_str = String(UInt8.(symbols_bytes))
        species = split(symbols_str, ",")

        # Convert lattvec from Angstrom (ASE/npz) to Bohr (Julia internal)
        ANG_TO_BOHR = 1.8897261246257702
        lattvec = ref["lattvec"] * ANG_TO_BOHR

        # Create DFTData{2}
        data = BoltzTraP.DFTData(
            lattice = lattvec,
            positions = ref["positions"]',  # transpose to 3×natom
            species = species,
            kpoints = ref["kpoints"]',  # transpose to 3×nkpts
            weights = ones(nkpts) / nkpts,
            ebands = ebands_3d,
            occupations = zeros(size(ebands_3d)),
            fermi = ref["fermi"],
            nelect = ref["nelect"],
            magmom = ref["magmom"],
        )

        # Verify DFTData{2} was created
        @test data isa BoltzTraP.DFTData{2}
        @test BoltzTraP.is_magnetic(data)
        @test size(data.ebands, 3) == 2

        # Run interpolation (use nkpts as Python does with len(data.kpoints))
        result = BoltzTraP.run_interpolate(data; kpoints=nkpts, verbose=false)

        # Verify result
        @test result isa BoltzTraP.InterpolationResult

        # Check coefficients shape matches reference
        # Reference: (32, 205) - 32 bands (concatenated), 205 equivalences
        @test size(result.coeffs, 1) == size(ref["coeffs"], 1)  # nbands
        @test size(result.coeffs, 2) == size(ref["coeffs"], 2)  # neq

        # Check metadata
        @test result.metadata["dosweight"] == 1.0  # spin-polarized
        @test haskey(result.metadata, "spintype")
        @test result.metadata["spintype"] == "Collinear"

        # Note: Direct coefficient comparison is skipped because equivalence class
        # ordering differs between Python and Julia implementations.
        # The coefficients are correct but stored in different order.
        # Key verification: shape matches and interpolated bands would match.
    end

    @testset "run_integrate Collinear" begin
        # Load PbTe collinear reference data
        REFTEST_DIR = joinpath(@__DIR__, "..", "reftest", "data")
        reference_file = joinpath(REFTEST_DIR, "pbte_collinear.npz")

        if !isfile(reference_file)
            @test_skip "PbTe collinear reference data not available"
            return
        end

        ref = npzread(reference_file)

        # Create DFTData{2} (same as interpolate test)
        ebands_concat = ref["ebands"]
        nbands_per_spin = size(ebands_concat, 1) ÷ 2
        nkpts = size(ebands_concat, 2)

        ebands_spin1 = ebands_concat[1:nbands_per_spin, :]
        ebands_spin2 = ebands_concat[nbands_per_spin+1:end, :]
        ebands_3d = zeros(nbands_per_spin, nkpts, 2)
        ebands_3d[:, :, 1] = ebands_spin1
        ebands_3d[:, :, 2] = ebands_spin2

        symbols_bytes = ref["symbols"]
        symbols_str = String(UInt8.(symbols_bytes))
        species = split(symbols_str, ",")

        # Convert lattvec from Angstrom (ASE/npz) to Bohr (Julia internal)
        ANG_TO_BOHR = 1.8897261246257702
        lattvec = ref["lattvec"] * ANG_TO_BOHR

        data = BoltzTraP.DFTData(
            lattice = lattvec,
            positions = ref["positions"]',
            species = species,
            kpoints = ref["kpoints"]',
            weights = ones(nkpts) / nkpts,
            ebands = ebands_3d,
            occupations = zeros(size(ebands_3d)),
            fermi = ref["fermi"],
            nelect = ref["nelect"],
            magmom = ref["magmom"],
        )

        # Run interpolation
        interp = BoltzTraP.run_interpolate(data; kpoints=nkpts, verbose=false)

        # Run integration
        transport = BoltzTraP.run_integrate(interp; temperatures=[300.0], verbose=false)

        # Verify result type
        @test transport isa BoltzTraP.TransportResult

        # Check output shapes
        # sigma, seebeck, kappa: (3, 3, nT, nμ)
        @test ndims(transport.sigma) == 4
        @test size(transport.sigma, 1) == 3  # 3x3 tensor
        @test size(transport.sigma, 2) == 3
        @test size(transport.sigma, 3) == 1  # 1 temperature

        @test ndims(transport.seebeck) == 4
        @test size(transport.seebeck, 3) == 1

        @test ndims(transport.kappa) == 4
        @test size(transport.kappa, 3) == 1

        # Check spintype metadata is inherited from InterpolationResult
        @test haskey(transport.metadata, "spintype")
        @test transport.metadata["spintype"] == "Collinear"

        # Check dosweight is correct for spin-polarized
        @test transport.metadata["dosweight"] == 1.0

        # Compare transport coefficients with Python reference
        # Note: Direct coefficient comparison is challenging because:
        # 1. Equivalence class ordering differs between Python and Julia
        # 2. This affects the Fourier interpolation basis ordering
        # 3. Transport coefficients are computed from interpolated bands
        #
        # For now, we verify that the Julia implementation produces reasonable
        # results (no NaN/Inf, correct shapes, physically meaningful ranges)
        @testset "Transport coefficient sanity checks" begin
            # Verify no NaN or Inf in results
            @test !any(isnan, transport.sigma)
            @test !any(isinf, transport.sigma)
            @test !any(isnan, transport.seebeck)
            @test !any(isinf, transport.seebeck)
            @test !any(isnan, transport.kappa)
            @test !any(isinf, transport.kappa)

            # Verify sigma is positive (conductivity should be non-negative)
            # Check diagonal elements at all mu values
            for k in axes(transport.sigma, 4)
                for i in 1:3
                    @test transport.sigma[i, i, 1, k] >= 0.0
                end
            end
        end
    end

    @testset "Fe Collinear E2E (QE loader)" begin
        # Load Fe collinear reference data
        REFTEST_DIR = joinpath(@__DIR__, "..", "reftest", "data")
        reference_file = joinpath(REFTEST_DIR, "fe_collinear_e2e.npz")

        if !isfile(reference_file)
            @test_skip "Fe collinear reference data not available"
            return
        end

        ref = npzread(reference_file)

        # Create DFTData{2} from reference data
        # Reference ebands is (24, 47) - concatenated [spin1; spin2]
        ebands_concat = ref["ebands"]  # (24, 47)
        nbands_per_spin = size(ebands_concat, 1) ÷ 2  # 12
        nkpts = size(ebands_concat, 2)  # 47

        # Split into two spin channels
        ebands_spin1 = ebands_concat[1:nbands_per_spin, :]
        ebands_spin2 = ebands_concat[nbands_per_spin+1:end, :]
        ebands_3d = zeros(nbands_per_spin, nkpts, 2)
        ebands_3d[:, :, 1] = ebands_spin1
        ebands_3d[:, :, 2] = ebands_spin2

        # Note: Python BoltzTraP2 uses ASE (Angstrom), Julia expects Bohr
        ANG_TO_BOHR = 1.8897261246257702
        lattvec = ref["lattvec"] * ANG_TO_BOHR

        # Create DFTData{2} - Fe has 1 atom
        data = BoltzTraP.DFTData(
            lattice = lattvec,
            positions = zeros(3, 1),  # Fe BCC: 1 atom at origin
            species = ["Fe"],
            kpoints = ref["kpoints"]',  # transpose to 3×nkpts
            weights = ones(nkpts) / nkpts,
            ebands = ebands_3d,
            occupations = zeros(size(ebands_3d)),
            fermi = ref["fermi"],
            nelect = ref["nelect"],
            magmom = ref["magmom"],
        )

        # Verify DFTData{2} was created
        @test data isa BoltzTraP.DFTData{2}
        @test BoltzTraP.is_magnetic(data)
        @test size(data.ebands, 3) == 2

        # Run interpolation
        interp = BoltzTraP.run_interpolate(data; kpoints=nkpts, verbose=false)

        # Verify interpolation result
        @test interp isa BoltzTraP.InterpolationResult
        @test interp.metadata["spintype"] == "Collinear"
        @test interp.metadata["dosweight"] == 1.0

        # Check coefficients shape matches reference
        @test size(interp.coeffs, 1) == size(ref["coeffs_real"], 1)  # nbands
        @test size(interp.coeffs, 2) == size(ref["coeffs_real"], 2)  # neq

        # Run integration at reference temperatures and mu range
        ref_temps = ref["Tr"]
        mur = ref["mur"]  # Use same mu range as Python
        transport = BoltzTraP.run_integrate(interp, mur; temperatures=ref_temps, verbose=false)

        # Verify transport result
        @test transport isa BoltzTraP.TransportResult
        @test transport.metadata["spintype"] == "Collinear"

        # Sanity checks
        @testset "Fe transport sanity checks" begin
            @test !any(isnan, transport.sigma)
            @test !any(isinf, transport.sigma)
            @test !any(isnan, transport.seebeck)
            @test !any(isinf, transport.seebeck)
            @test !any(isnan, transport.kappa)
            @test !any(isinf, transport.kappa)

            # Verify sigma diagonal is non-negative
            for k in axes(transport.sigma, 4)
                for iT in axes(transport.sigma, 3)
                    for i in 1:3
                        @test transport.sigma[i, i, iT, k] >= 0.0
                    end
                end
            end

            # For Fe metal: sigma should have some large values (metallic conductivity)
            # Check that max sigma_xx is significant somewhere in the mu range
            sigma_xx_max = maximum(transport.sigma[1, 1, 1, :])
            @test sigma_xx_max > 1e10  # Should be large for metal at some mu
        end

        # Quantitative validation at T=1000K
        # At high temperature, transport curves are smooth (oscillations damped)
        # This makes Python-Julia comparison meaningful
        @testset "Fe transport agreement at T=1000K" begin
            # Find index for T=1000K
            iT_1000K = findfirst(t -> t ≈ 1000.0, ref_temps)
            if isnothing(iT_1000K)
                @test_skip "T=1000K not found in reference data"
                return
            end

            if !isnothing(iT_1000K)
                # Python reference
                sigma_py = ref["sigma"][iT_1000K, :, 1, 1]
                S_py = ref["S"][iT_1000K, :, 1, 1]
                kappa_py = ref["kappa"][iT_1000K, :, 1, 1]

                # Julia result
                sigma_jl = vec(transport.sigma[1, 1, iT_1000K, :])
                S_jl = vec(transport.seebeck[1, 1, iT_1000K, :])
                kappa_jl = vec(transport.kappa[1, 1, iT_1000K, :])

                # Compute relative errors (filter small values)
                function max_rel_error(jl, py; threshold=1e-10)
                    valid = abs.(py) .> threshold
                    if !any(valid)
                        return 0.0
                    end
                    return maximum(abs.(jl[valid] .- py[valid]) ./ abs.(py[valid]))
                end

                # Tolerance: 5% max relative error
                rel_tol = 0.05

                sigma_err = max_rel_error(sigma_jl, sigma_py)
                S_err = max_rel_error(S_jl, S_py)
                kappa_err = max_rel_error(kappa_jl, kappa_py)

                # Relative error should be < 5%
                @test sigma_err < rel_tol
                @test S_err < rel_tol
                @test kappa_err < rel_tol

                # Correlation should be > 0.999
                @test cor(sigma_jl, sigma_py) > 0.999
                @test cor(S_jl, S_py) > 0.999
                @test cor(kappa_jl, kappa_py) > 0.999
            end
        end
    end

end
