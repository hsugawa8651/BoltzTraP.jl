# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Test
using BoltzTraP
using Wannier

@testset "WannierInterpolator (Wannier.jl ext)" begin
    # Si fixture shipped with Wannier.jl (num_wann=8, num_bands=12, mp_grid=4×4×4)
    sifix = joinpath(pkgdir(Wannier), "test", "fixtures", "silicon", "silicon")
    @assert isfile("$sifix.win") "Wannier.jl Si fixture missing — package layout may have changed"

    @testset "construction from Wannier90 prefix" begin
        wi = WannierInterpolator(sifix)
        @test wi isa BoltzTraP.AbstractInterpolator
        @test wi isa WannierInterpolator
        @test wi.spintype isa Unpolarized
        @test size(wi.lattice) == (3, 3)
    end

    @testset "construction from Wannier.InterpModel" begin
        imodel = Wannier.read_w90_interp(sifix)
        wi = WannierInterpolator(imodel)
        @test wi isa BoltzTraP.AbstractInterpolator
        @test wi.model === imodel
    end

    @testset "construction equivalence (prefix vs InterpModel)" begin
        wi_prefix = WannierInterpolator(sifix)
        imodel = Wannier.read_w90_interp(sifix)
        wi_imodel = WannierInterpolator(imodel)

        # Both paths must produce the same lattice (Bohr conversion identical)
        @test wi_prefix.lattice == wi_imodel.lattice

        # And the same band energies on a few k-points
        kpts = [
            0.0 0.0 0.0
            0.5 0.0 0.0
            0.25 0.25 0.25
        ]
        eband_prefix = BoltzTraP.interpolate_bands(wi_prefix, kpts)
        eband_imodel = BoltzTraP.interpolate_bands(wi_imodel, kpts)
        @test eband_prefix == eband_imodel
    end

    @testset "interpolate_bands shape and units" begin
        wi = WannierInterpolator(sifix)
        kpts = [
            0.0 0.0 0.0
            0.5 0.0 0.0
            0.25 0.25 0.25
        ]
        ebands = BoltzTraP.interpolate_bands(wi, kpts)
        @test size(ebands) == (8, 3)               # (n_wann, nk)
        @test eltype(ebands) <: AbstractFloat
        # Si valence/conduction span ~25 eV total → in Hartree, |E| < 1
        @test all(abs.(ebands) .< 1.0)
    end

    @testset "interpolate_velocities shape, units, and conversion factor" begin
        wi = WannierInterpolator(sifix)
        kpts = [
            0.0 0.0 0.0
            0.5 0.0 0.0
            0.25 0.25 0.25
        ]
        vbands = BoltzTraP.interpolate_velocities(wi, kpts)
        @test size(vbands) == (3, 8, 3)            # (3, n_wann, nk)
        @test eltype(vbands) <: AbstractFloat

        # BT.jl output must equal Wannier raw output (eV·Å), permuted to the
        # BoltzTraP shape `(3, n_wann, nk)` and rescaled to Hartree·Bohr.
        # Wannier.velocity returns `(E, vᴴ)`; only `vᴴ` is needed here.
        kpts_3xN = Matrix{Float64}(transpose(kpts))
        _, v_eVA = Wannier.velocity(wi.model.kRvectors.Rvectors, wi.model.H, kpts_3xN)
        factor = 1 / (BoltzTraP.HA_TO_EV * BoltzTraP.BOHR_TO_ANG)
        @test vbands ≈ permutedims(v_eVA, (3, 1, 2)) .* factor
    end

    @testset "solve_transport smoke (UniformMesh + WannierInterpolator)" begin
        wi = WannierInterpolator(sifix)
        # Mock TransportSystem for Si: 8 valence electrons, dosweight 2.0,
        # Fermi level ≈ 5.5 eV ≈ 0.2 Hartree (per Phase 0c smoke value)
        sys = TransportSystem(wi.lattice, 0.2, 8.0, 2.0, Unpolarized())
        result = solve_transport(UniformMesh((4, 4, 4)), wi, sys; temperatures = [300.0])
        @test result isa TransportResult
        @test result.temperatures == [300.0]
        # σ shape is (3, 3, nT, nμ): leading two axes are the Cartesian
        # tensor, axis 3 is temperature, axis 4 is the auto-generated μ grid.
        @test size(result.sigma)[1:3] == (3, 3, 1)
        @test size(result.sigma, 4) > 0
    end
end
