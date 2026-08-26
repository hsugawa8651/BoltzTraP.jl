# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Absolute-value unit regression test. Self-consistency between the Fourier
# and Wannier paths (or between BT.jl and a stale reference) cannot catch a
# unit-conversion bug that both sides share — this is exactly how the
# σ/κ_el inflation bug (Fourier ×1.8897, Wannier ×6.7483) went undetected.
#
# The check here is independent of DFT and of BoltzTraP2: a parabolic band
# E(k) = |k|²/2 (atomic units, m*=1) has an exact CRTA solution
# σ/τ = n e²/m*, so the code's output can be compared directly against a
# closed-form value derived from first principles, not against another run
# of the same code.

using Test
using BoltzTraP
using BoltzTraP: TransportSystem, UniformMesh, Unpolarized, _solve_transport_from_bands
using LinearAlgebra
using StaticArrays

# σ/τ = n e²/m* (SI), with n the carrier density derived from the code's own
# BTPDOS/calc_N. Only fundamental constants and the lattice enter — no
# BoltzTraP2 comparison, no stored reference file.
const _BOHR_M = 5.29177210903e-11   # m, CODATA (exact by SI definition since 2019)
const _QE = 1.602176634e-19         # C, elementary charge (exact, SI)
const _ME = 9.1093837015e-31        # kg, electron mass (CODATA 2018)

struct _DrudeDummyInterpolator <: BoltzTraP.AbstractInterpolator end

# Build a single parabolic band E(k)=|k|²/2 (m*=1, atomic units) on a
# UniformMesh nk×nk×nk for the reciprocal lattice implied by `lattice`
# (Bohr, columns are lattice vectors). vvband is the exact analytic v⊗v = k⊗k.
function _build_parabolic_band(lattice::AbstractMatrix{Float64}, nk::Int)
    B = 2π * inv(Matrix(lattice)')   # reciprocal lattice, 1/Bohr
    npts = nk^3
    eband = zeros(1, npts)
    vvband = zeros(1, 3, 3, npts)
    i = 1
    for a1 = 0:(nk-1), a2 = 0:(nk-1), a3 = 0:(nk-1)
        f = (((a1, a2, a3) .+ nk ÷ 2) .% nk .- nk ÷ 2) ./ nk
        k = B * collect(f)
        eband[1, i] = dot(k, k) / 2
        for x = 1:3, y = 1:3
            vvband[1, x, y, i] = k[x] * k[y]
        end
        i += 1
    end
    eband, vvband
end

# Exact σ/τ [S/(m·s)] for the parabolic band, from the code's own carrier
# density (BTPDOS + calc_N) so only the Onsager-coefficient step is checked.
function _drude_ratio(lattice::AbstractMatrix{Float64}, nk::Int, mu::Float64)
    eband, vvband = _build_parabolic_band(lattice, nk)
    sys = TransportSystem(SMatrix{3,3,Float64}(lattice), 0.0, 1.0e-4, 2.0, Unpolarized())
    r = _solve_transport_from_bands(
        eband,
        vvband,
        sys,
        UniformMesh((nk, nk, nk)),
        _DrudeDummyInterpolator();
        temperatures = [300.0],
        mur = [mu],
        bins = 20000,
        scissor = nothing,
        verbose = false,
    )
    epsilon = r.dos["epsilon_Ha"]
    dos = r.dos["dos"]
    Ncell = -BoltzTraP.calc_N(epsilon, dos, mu, 300.0; dosweight = 2.0)  # electrons/cell
    Vm3 = abs(det(Matrix(lattice))) * _BOHR_M^3
    n_m3 = Ncell / Vm3
    exact = n_m3 * _QE^2 / _ME
    r.sigma[1, 1, 1, 1] / exact
end

@testset "Absolute units (Drude exact solution)" begin
    @testset "Parabolic band, cubic lattice, multiple μ" begin
        lat = 10.0 * Matrix(1.0I, 3, 3)  # Bohr
        for mu in (-0.010, -0.008, -0.006)  # Ha, below the band bottom (non-degenerate)
            ratio = _drude_ratio(lat, 40, mu)
            @test isapprox(ratio, 1.0; rtol = 1e-2)
        end
    end

    @testset "Lattice-shape independence" begin
        # The bug's signature was a lattice-dependent multiplier (1/a₀ or
        # 1/a₀³) that a single-lattice check cannot distinguish from a
        # correct, lattice-independent result. Cover cubic, tetragonal, and
        # a non-orthogonal (fcc) cell.
        mu = -0.008
        for lat in (
            8.0 * Matrix(1.0I, 3, 3),
            15.0 * Matrix(1.0I, 3, 3),
            diagm([8.0, 10.0, 12.0]),
            10.26 / 2 * [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0],
        )
            ratio = _drude_ratio(lat, 45, mu)
            @test isapprox(ratio, 1.0; rtol = 1e-2)
        end
    end
end
