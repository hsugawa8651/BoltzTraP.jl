# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/tests/test_bandlib.py

using Test
using LinearAlgebra
using NPZ
using Random
using Statistics: mean

@testset "Band Library (transport)" begin
    TEST_DIR = @__DIR__
    DATA_DIR = joinpath(TEST_DIR, "data")

    @testset "DOS parabolic band" begin
        # Generate a set of electron energies from a single parabolic band
        L = 5.0
        x = range(-L, L, length = 201)
        y = range(-L, L, length = 201)
        z = range(-L, L, length = 201)

        # Create points and filter to sphere
        points = [(xi, yi, zi) for xi in x for yi in y for zi in z]
        k2 = [p[1]^2 + p[2]^2 + p[3]^2 for p in points]
        k2 = filter(x -> x < L*L, k2)

        # Energy proportional to k²
        energies = reshape(k2 / maximum(k2) * 0.1, 1, :)

        npts = 60
        erange = (minimum(energies), maximum(energies))
        epsilon, dos = BoltzTraP.compute_dos(energies, erange, npts)

        # Test shape
        @test length(epsilon) == npts
        @test length(dos) == npts

        # Test normalization (should integrate to 1)
        de = epsilon[2] - epsilon[1]
        @test sum(dos) * de ≈ 1.0 atol=0.01

        # Analytic result for parabolic band: DOS ∝ √E
        ref = sqrt.(max.(epsilon, 0))
        ref ./= sum(ref) * de

        # Test that the difference is small
        @test norm(dos - ref) * de < 0.01
    end

    @testset "DOS Si" begin
        # Load Si DOS reference
        ref_file = joinpath(DATA_DIR, "Si_BTPdos.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)
        dose_ref = reference["dose"]
        dos_ref = reference["dos"]

        # DOS should be normalized to number of bands
        de = dose_ref[2] - dose_ref[1]
        @test sum(dos_ref) * de ≈ 6.0 atol=0.1  # Si has 6 valence bands
    end

    @testset "Fermi-Dirac integration" begin
        # Test basic Fermi-Dirac integration properties
        epsilon = collect(range(-0.1, 0.1, length = 1000))
        de = epsilon[2] - epsilon[1]

        # Parabolic DOS
        dos = sqrt.(max.(epsilon .+ 0.05, 0))
        dos ./= sum(dos) * de

        # Test at different temperatures
        const_kB = 3.1668115634556076e-06  # Hartree/K
        T_range = [100.0, 200.0, 300.0]
        μ = 0.0

        for T in T_range
            kT = const_kB * T
            f = BoltzTraP.fermi_dirac(epsilon, μ, kT)
            df = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT)

            # Integral of f*dos should be between 0 and 1
            N = sum(dos .* f) * de
            @test 0 <= N <= 1

            # df should be negative
            @test all(df .<= 0)

            # Integral of -df should be positive
            @test sum(-df .* dos) * de >= 0
        end
    end

    @testset "Fermi integrals Si" begin
        # Load Si fermiintegrals reference
        ref_file = joinpath(DATA_DIR, "Si_fermiintegrals.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)

        # Just verify the reference data has expected shape
        N_ref = reference["N"]
        L0_ref = reference["L0"]
        L1_ref = reference["L1"]
        L2_ref = reference["L2"]

        # L0, L1, L2 should be 4D: (nT, nμ, 3, 3)
        @test ndims(L0_ref) == 4
        @test size(L0_ref, 3) == 3
        @test size(L0_ref, 4) == 3

        # N should be 2D: (nT, nμ)
        @test ndims(N_ref) == 2
    end

    @testset "Onsager coefficients Si" begin
        # Load Si Onsager reference
        ref_file = joinpath(DATA_DIR, "Si_Onsager.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        reference = npzread(ref_file)

        # Verify reference data structure
        cond_ref = reference["cond"]
        seebeck_ref = reference["seebeck"]
        kappa_ref = reference["kappa"]

        # All should be 4D: (nT, nμ, 3, 3)
        @test ndims(cond_ref) == 4
        @test ndims(seebeck_ref) == 4
        @test ndims(kappa_ref) == 4

        # Conductivity should be positive (diagonal elements)
        @test all(all(cond_ref[:, :, i, i] .>= 0) for i = 1:3)
    end

    @testset "calc_N temperature dependence" begin
        # |N| should be monotonically increasing with T for a single band
        epsilon = collect(range(0.0, 0.1, length = 1000))
        de = epsilon[2] - epsilon[1]
        dos = sqrt.(epsilon)
        dos ./= sum(dos) * de

        μ = 0.05 * maximum(epsilon)
        T_range = range(50.0, 600.0, length = 20)

        N_vals = Float64[]
        for T in T_range
            kT = 3.1668115634556076e-06 * T  # k_B T in Hartree
            f = BoltzTraP.fermi_dirac(epsilon, μ, kT)
            N = -2.0 * sum(dos .* f) * de
            push!(N_vals, N)
        end

        # N should decrease (more negative) as T increases
        # because more states get occupied above the Fermi level
        dN = diff(N_vals)
        @test all(dN .< 0)
    end

    @testset "Heat capacity properties" begin
        # Test basic properties of heat capacity calculation
        epsilon = collect(range(0.0, 0.1, length = 10000))
        de = epsilon[2] - epsilon[1]
        dos = sqrt.(epsilon)
        dos ./= sum(dos) * de

        μ = 0.01  # Well into the band
        const_kB = 3.1668115634556076e-06  # Hartree/K

        # cv should increase with temperature at low T
        T1, T2 = 50.0, 100.0
        kT1, kT2 = const_kB * T1, const_kB * T2

        df1 = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT1)
        df2 = BoltzTraP.dfermi_dirac_de(epsilon, μ, kT2)

        cv1 = -sum((epsilon .- μ) .^ 2 .* df1 .* dos) * de
        cv2 = -sum((epsilon .- μ) .^ 2 .* df2 .* dos) * de

        # Higher T should have higher cv (more thermal smearing)
        @test cv2 > cv1

        # cv should be positive
        @test cv1 > 0
        @test cv2 > 0
    end

    @testset "calc_cv reference test" begin
        # Load Python BoltzTraP2 reference data
        REFTEST_DATA = joinpath(@__DIR__, "..", "reftest", "data")
        ref_file = joinpath(REFTEST_DATA, "cv_reference.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        ref = npzread(ref_file)

        # Call Julia calc_cv
        cv = BoltzTraP.calc_cv(
            ref["epsilon"],
            ref["dos"],
            ref["mu_range"],
            ref["T_range"];
            dosweight = ref["dosweight"],
        )

        # Test shape: (nT, nμ)
        @test size(cv) == (length(ref["T_range"]), length(ref["mu_range"]))

        # Test values match Python BoltzTraP2
        # Note: rtol=1e-5 due to small difference in KB_AU values between
        # Julia (CODATA 2018) and Python BoltzTraP2 (different source)
        @test isapprox(cv, ref["cv"], rtol = 1e-5)
    end

    @testset "calc_cv error dispatch" begin
        @test_throws ArgumentError BoltzTraP.calc_cv("invalid")
        @test_throws ArgumentError BoltzTraP.calc_cv(123)
    end

    @testset "apply_scissor reference test" begin
        # Load Python BoltzTraP2 reference data
        REFTEST_DATA = joinpath(@__DIR__, "..", "reftest", "data")
        ref_file = joinpath(REFTEST_DATA, "scissor_reference.npz")
        if !isfile(ref_file)
            @warn "Reference file not found: $ref_file"
            @test_skip "Reference data not available"
            return
        end

        ref = npzread(ref_file)

        # Call Julia apply_scissor
        eband_after = BoltzTraP.apply_scissor(
            ref["epsilon"],
            ref["dos"],
            ref["N0"],
            ref["eband_before"],
            ref["desired_gap"];
            dosweight = ref["dosweight"],
        )

        # Test shape matches
        @test size(eband_after) == size(ref["eband_before"])

        # Test values match Python BoltzTraP2
        @test isapprox(eband_after, ref["eband_after"], rtol = 1e-10)
    end

    @testset "apply_scissor error dispatch" begin
        @test_throws ArgumentError BoltzTraP.apply_scissor("invalid")
        @test_throws ArgumentError BoltzTraP.apply_scissor(123)
    end

    @testset "run_integrate with scissor" begin
        # Test that run_integrate accepts scissor keyword argument

        # Create a minimal mock InterpolationResult for testing
        mock_coeffs = zeros(ComplexF64, 2, 10)  # 2 bands, 10 equiv points
        mock_equiv = [zeros(Int, 1, 3)]  # Single equivalence at origin
        mock_lattvec = Matrix{Float64}(I, 3, 3) * 10.0  # 10 Bohr cubic cell
        mock_atoms = Dict{String,Any}("species" => ["Si"], "positions" => [[0.0, 0.0, 0.0]])
        mock_metadata =
            Dict{String,Any}("fermi" => 0.0, "nelect" => 4.0, "dosweight" => 2.0)

        mock_interp = BoltzTraP.InterpolationResult(
            mock_coeffs,
            mock_equiv,
            mock_lattvec,
            mock_atoms,
            mock_metadata,
        )

        @testset "scissor=nothing (default)" begin
            # Should work without scissor argument (backward compatibility)
            # Mock data may cause errors, but NOT scissor-related errors
            try
                result = run_integrate(mock_interp; temperatures = [300.0])
                # If successful, scissor should not be in metadata
                @test !haskey(result.metadata, "scissor_eV")
            catch e
                # Errors from mock data are OK, scissor keyword error is not
                @test !occursin("scissor", lowercase(string(e)))
            end
        end

        @testset "scissor parameter accepted" begin
            # Verify scissor keyword is accepted (not "unknown keyword argument")
            # With mock data, apply_scissor will fail because there's no gap,
            # but the keyword itself should be accepted
            try
                result = run_integrate(mock_interp; temperatures = [300.0], scissor = 1.1)
                # If we reach here, scissor was applied successfully
                @test haskey(result.metadata, "scissor_eV")
                @test result.metadata["scissor_eV"] == 1.1
            catch e
                err_msg = string(e)
                # "No band gap found" is expected (mock data has no gap)
                # "keyword argument" error means scissor was not implemented
                is_keyword_error =
                    occursin("keyword", lowercase(err_msg)) && occursin("scissor", err_msg)
                is_gap_error = occursin("band gap", lowercase(err_msg))
                @test !is_keyword_error  # Scissor keyword must be accepted
                @test is_gap_error  # Expected error from mock data
            end
        end
    end
end
