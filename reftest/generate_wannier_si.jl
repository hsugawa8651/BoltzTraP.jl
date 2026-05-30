# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Generate the golden Si Wannier transport regression baseline used by
# `test/test_wannier_reftest.jl`.
#
# This script is fully self-contained: it reads the bundled silicon
# fixture shipped with Wannier.jl (num_wann=8, num_bands=12, mp_grid=4×4×4)
# and writes a JLD2 file containing the σ/S/κ tensors plus the auto-
# generated μ grid for the chosen UniformMesh and temperature set.
#
# This script assumes the current Julia environment provides
# `BoltzTraP`, `Wannier`, and `JLD2`. The usual idiom is to use a
# throw-away environment that develops the local BoltzTraP source and
# adds `Wannier` / `JLD2`:
#
#     julia -e '
#         using Pkg
#         Pkg.activate(; temp=true)
#         Pkg.develop(path = "<path to BoltzTraP.jl>")
#         Pkg.add(["Wannier", "JLD2"])
#         include("<path to BoltzTraP.jl>/reftest/generate_wannier_si.jl")
#     '
#
# The output `reftest/data/si_wannier_transport.jld2` is then loaded by
# the reftest in `test/test_wannier_reftest.jl` for an element-wise
# tight-rtol comparison against the same recomputation.

using BoltzTraP
using Wannier
using JLD2
using Dates

# -----------------------------------------------------------------------------
# Inputs (kept in sync with `test/test_wannier_reftest.jl`)
# -----------------------------------------------------------------------------

const SI_FIXTURE_PREFIX =
    joinpath(pkgdir(Wannier), "test", "fixtures", "silicon", "silicon")
const MESH_NK = (16, 16, 16)
const TEMPERATURES = [300.0, 500.0, 700.0]

# Smoke-level Si parameters, matching `test/test_wannier_ext.jl`.
# These freeze the inputs so the regression is reproducible; they are
# not intended as physically accurate values for Si.
const FERMI_HARTREE = 0.2
const NELECT = 8.0
const DOSWEIGHT = 2.0

# -----------------------------------------------------------------------------

@assert isfile("$SI_FIXTURE_PREFIX.win") "Wannier.jl Si fixture missing at $SI_FIXTURE_PREFIX.win — package layout may have changed"

@info "Building WannierInterpolator from bundled Si fixture" SI_FIXTURE_PREFIX
const wi = WannierInterpolator(SI_FIXTURE_PREFIX)
@info "WannierInterpolator built" n_wann = size(wi.model.H, 1)

const sys = TransportSystem(wi.lattice, FERMI_HARTREE, NELECT, DOSWEIGHT, Unpolarized())
const mesh = UniformMesh(MESH_NK)

@info "Solving transport" mesh = MESH_NK temperatures = TEMPERATURES
const result = solve_transport(mesh, wi, sys; temperatures = TEMPERATURES)
@info "Result shape" sigma_size = size(result.sigma) n_mu = length(result.mu_values)

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

const OUT_DIR = joinpath(@__DIR__, "data")
mkpath(OUT_DIR)
const OUT_FILE = joinpath(OUT_DIR, "si_wannier_transport.jld2")

JLD2.jldopen(OUT_FILE, "w") do io
    # Transport tensors
    io["sigma"] = result.sigma
    io["seebeck"] = result.seebeck
    io["kappa"] = result.kappa
    io["mu_values"] = result.mu_values
    io["temperatures"] = result.temperatures
    # Provenance / reproducibility
    io["mesh_nk"] = collect(MESH_NK)
    io["fermi_hartree"] = FERMI_HARTREE
    io["nelect"] = NELECT
    io["dosweight"] = DOSWEIGHT
    io["fixture_prefix"] = "Wannier.jl bundled silicon fixture (pkgdir(Wannier)/test/fixtures/silicon/silicon)"
    io["wannier_pkg_version"] = string(pkgversion(Wannier))
    io["boltztrap_pkg_version"] = string(pkgversion(BoltzTraP))
    io["generated_at"] = string(Dates.now())
end

const out_size = stat(OUT_FILE).size
@info "Saved" OUT_FILE size_kb = out_size / 1024
