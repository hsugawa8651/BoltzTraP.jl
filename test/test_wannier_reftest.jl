# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Regression baseline for the Wannier transport path on Si.
#
# Re-runs the same `WannierInterpolator → solve_transport` pipeline that
# `reftest/generate_wannier_si.jl` used to produce the golden file, and
# compares the σ/S/κ tensors against the stored result. The velocity tensor
# at degenerate k-points is built gauge-invariantly (averaged over approach
# directions), so the result is reproducible across LAPACK builds and the
# comparison is a tight, full-tensor `rtol=1e-8`.

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

    # Full-tensor regression. The transport tensors are built from a
    # degeneracy-aware, gauge-invariant velocity construction (the band
    # velocities at degenerate k-points are averaged over approach
    # directions), so the result no longer depends on the arbitrary
    # eigenbasis a given LAPACK build picks inside a degenerate multiplet.
    # The recomputation therefore matches the golden across platforms to
    # near machine precision (≈ 2e-13 norm, x86_64 vs aarch64), and both
    # diagonal and off-diagonal tensor entries are checked.
    # Off-diagonal tensor entries are symmetry zeros, so compare with a
    # per-tensor `atol` floor (`rtol` × the tensor's largest entry) in
    # addition to `rtol`: the diagonal is checked at 1e-8 relative, while
    # the near-zero off-diagonal entries are checked against the same
    # absolute floor instead of an ill-defined relative tolerance.
    atol_sigma = 1e-8 * maximum(abs, golden["sigma"])
    atol_seebeck = 1e-8 * maximum(abs, golden["seebeck"])
    atol_kappa = 1e-8 * maximum(abs, golden["kappa"])
    for i = 1:3, j = 1:3
        @test isapprox(
            result.sigma[i, j, :, :],
            golden["sigma"][i, j, :, :];
            rtol = 1e-8,
            atol = atol_sigma,
        )
        @test isapprox(
            result.seebeck[i, j, :, :],
            golden["seebeck"][i, j, :, :];
            rtol = 1e-8,
            atol = atol_seebeck,
        )
        @test isapprox(
            result.kappa[i, j, :, :],
            golden["kappa"][i, j, :, :];
            rtol = 1e-8,
            atol = atol_kappa,
        )
    end
end
