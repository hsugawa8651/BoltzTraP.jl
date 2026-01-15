# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using NPZ
using BoltzTraP: calc_N, solve_for_mu

@testset "calc_N and solve_for_mu" begin
    # Load reference data
    ref_path = joinpath(@__DIR__, "..", "reftest", "data", "calc_N_solve_mu.npz")
    ref = npzread(ref_path)

    epsilon = ref["epsilon"]
    dos = ref["dos"]
    fermi = ref["fermi"]
    dosweight = ref["dosweight"]
    mu_test = ref["mu_test"]
    T_test = ref["T_test"]
    calc_N_ref = ref["calc_N_results"]
    N0_test = ref["N0_test"]
    solve_mu_basic_ref = ref["solve_mu_basic"]
    solve_mu_refine_ref = ref["solve_mu_refine"]
    solve_mu_center_ref = ref["solve_mu_center"]

    # Histogram bin size (for tolerance estimation)
    de = epsilon[2] - epsilon[1]

    @testset "calc_N" begin
        for (iT, T) in enumerate(T_test)
            for (iμ, μ) in enumerate(mu_test)
                N_julia = calc_N(epsilon, dos, μ, T; dosweight)
                N_python = calc_N_ref[iT, iμ]
                @test isapprox(N_julia, N_python; rtol=1e-6)
            end
        end
    end

    @testset "solve_for_mu (basic)" begin
        for (iT, T) in enumerate(T_test)
            for (iN, N0) in enumerate(N0_test)
                μ_ref = solve_mu_basic_ref[iT, iN]
                if !isnan(μ_ref)
                    μ_julia = solve_for_mu(epsilon, dos, N0, T; dosweight, refine=false, try_center=false)
                    # Allow tolerance of 1 histogram bin
                    @test abs(μ_julia - μ_ref) <= de + 1e-10
                end
            end
        end
    end

    @testset "solve_for_mu (refine)" begin
        for (iT, T) in enumerate(T_test)
            for (iN, N0) in enumerate(N0_test)
                μ_ref = solve_mu_refine_ref[iT, iN]
                if !isnan(μ_ref)
                    μ_julia = solve_for_mu(epsilon, dos, N0, T; dosweight, refine=true, try_center=false)
                    # At T=0, bisection behavior may differ slightly; allow 2 bins
                    tol = T == 0.0 ? 2 * de : de / 10
                    @test abs(μ_julia - μ_ref) <= tol
                end
            end
        end
    end

    @testset "solve_for_mu (center)" begin
        for (iT, T) in enumerate(T_test)
            for (iN, N0) in enumerate(N0_test)
                μ_ref = solve_mu_center_ref[iT, iN]
                if !isnan(μ_ref)
                    μ_julia = solve_for_mu(epsilon, dos, N0, T; dosweight, refine=true, try_center=true)
                    # At T=0, bisection behavior may differ slightly; allow 2 bins
                    tol = T == 0.0 ? 2 * de : de / 10
                    @test abs(μ_julia - μ_ref) <= tol
                end
            end
        end
    end
end
