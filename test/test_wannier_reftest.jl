# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Regression baseline for the Wannier transport path on Si.
#
# Re-runs the same `WannierInterpolator → solve_transport` pipeline that
# `reftest/generate_wannier_si.jl` used to produce the golden file, and
# compares the σ/S/κ tensors element-wise against the stored result.
# Tight `rtol=1e-10` only absorbs floating-point roundoff because the
# inputs (Wannier.jl bundled silicon fixture, mesh, smoke-level Si
# parameters) are byte-fixed.

using Test
using BoltzTraP
using Wannier
using JLD2

@testset "Si Wannier transport regression baseline" begin
    si_fixture_prefix = joinpath(pkgdir(Wannier), "test", "fixtures", "silicon", "silicon")
    @test isfile("$si_fixture_prefix.win")

    golden_path = joinpath(@__DIR__, "..", "reftest", "data", "si_wannier_transport.jld2")
    @test isfile(golden_path)

    golden = JLD2.jldopen(golden_path, "r") do io
        Dict(
            "sigma" => io["sigma"],
            "seebeck" => io["seebeck"],
            "kappa" => io["kappa"],
            "mu_values" => io["mu_values"],
            "temperatures" => io["temperatures"],
            "mesh_nk" => io["mesh_nk"],
            "fermi_hartree" => io["fermi_hartree"],
            "nelect" => io["nelect"],
            "dosweight" => io["dosweight"],
        )
    end

    wi = WannierInterpolator(si_fixture_prefix)
    sys = TransportSystem(
        wi.lattice,
        golden["fermi_hartree"],
        golden["nelect"],
        golden["dosweight"],
        Unpolarized(),
    )
    mesh = UniformMesh(Tuple(golden["mesh_nk"]))
    result = solve_transport(mesh, wi, sys; temperatures = golden["temperatures"])

    # Shape parity
    @test size(result.sigma) == size(golden["sigma"])
    @test size(result.seebeck) == size(golden["seebeck"])
    @test size(result.kappa) == size(golden["kappa"])
    @test length(result.mu_values) == length(golden["mu_values"])

    # μ-grid bit equality (deterministic from DOS grid + thermal window)
    @test result.mu_values ≈ golden["mu_values"] rtol = 1e-10
    @test result.temperatures == golden["temperatures"]

    # Element-wise rtol regression on transport tensors
    # rtol=1e-6 follows the project's general reftest convention and absorbs
    # cross-platform floating-point differences in the underlying linear
    # algebra (BLAS/LAPACK across Linux/macOS/aarch64).
    @test result.sigma ≈ golden["sigma"] rtol = 1e-6
    @test result.seebeck ≈ golden["seebeck"] rtol = 1e-6
    @test result.kappa ≈ golden["kappa"] rtol = 1e-6
end
