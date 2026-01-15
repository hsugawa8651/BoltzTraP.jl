# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/fite.py

# Note: BOHR_TO_ANG is defined in units.jl
# Python BoltzTraP2 uses Ångström for lattvec in velocity calculation
# We must match this for σ and κ to be consistent

#=
    determine_fft_dims(equivalences)

Determine FFT grid dimensions from equivalence classes.
=#
function determine_fft_dims(equivalences)
    allvec = vcat([Matrix(equiv) for equiv in equivalences]...)
    maxabs = maximum(abs.(allvec), dims=1)
    return Tuple(2 .* vec(maxabs) .+ 1)
end

# ============================================================================
# Optimized getBTPbands implementation using function barrier pattern
# Achieves ~4x speedup by:
# 1. Pre-computing Cartesian coordinates (lattvec * R)
# 2. Filling all 4 grids in a single loop
# 3. Using function barrier for type stability
# ============================================================================

#=
    _flatten_coeffs(bandcoeff, equivalences) -> Vector{ComplexF64}

Flatten band coefficients with normalization for FFT grid filling.
Each coefficient is normalized by the number of stars in its equivalence class.
=#
function _flatten_coeffs(bandcoeff::AbstractVector, equivalences::Vector)
    total_pts = sum(size(e, 1) for e in equivalences)
    flat = Vector{ComplexF64}(undef, total_pts)
    k = 1
    @inbounds for (c, equiv) in zip(bandcoeff, equivalences)
        nstar = size(equiv, 1)
        c_norm = c / nstar
        for _ in 1:nstar
            flat[k] = c_norm
            k += 1
        end
    end
    return flat
end

#=
    _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

Fill energy and velocity FFT grids in a single loop.
Uses pre-computed Cartesian coordinates (allvec_t) and indices for efficiency.
=#
function _fill_grids!(
    egrid::Array{ComplexF64,3},
    vgrid1::Array{ComplexF64,3},
    vgrid2::Array{ComplexF64,3},
    vgrid3::Array{ComplexF64,3},
    flat_coeffs::Vector{ComplexF64},
    allvec_t::Matrix{Float64},
    all_indices::Vector{CartesianIndex{3}},
)
    fill!(egrid, 0)
    fill!(vgrid1, 0)
    fill!(vgrid2, 0)
    fill!(vgrid3, 0)

    @inbounds for k in eachindex(all_indices)
        idx = all_indices[k]
        c = flat_coeffs[k]
        egrid[idx] = c
        vgrid1[idx] = im * allvec_t[1, k] * c
        vgrid2[idx] = im * allvec_t[2, k] * c
        vgrid3[idx] = im * allvec_t[3, k] * c
    end
end

#=
    _getBTPbands_kernel!(eband, vvband, coeffs, equivalences, allvec_t, all_indices, dims, plan)

Inner kernel for getBTPbands with concrete types (function barrier pattern).
This function is type-stable and optimized for performance.
=#
function _getBTPbands_kernel!(
    eband::Matrix{Float64},
    vvband::Union{Array{Float64,4},Nothing},
    coeffs::AbstractMatrix,
    equivalences::Vector{<:AbstractMatrix{<:Integer}},
    allvec_t::Matrix{Float64},
    all_indices::Vector{CartesianIndex{3}},
    dims::NTuple{3,Int},
    plan,
)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    scale = prod(dims)

    # Work arrays (reused across bands)
    egrid = zeros(ComplexF64, dims)
    vgrid1 = zeros(ComplexF64, dims)
    vgrid2 = zeros(ComplexF64, dims)
    vgrid3 = zeros(ComplexF64, dims)

    for iband in 1:nbands
        flat_coeffs = _flatten_coeffs(coeffs[iband, :], equivalences)
        _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

        eband[iband, :] = scale .* real.(vec(plan * egrid))

        if !isnothing(vvband)
            vb1 = scale .* real.(vec(plan * vgrid1))
            vb2 = scale .* real.(vec(plan * vgrid2))
            vb3 = scale .* real.(vec(plan * vgrid3))

            @inbounds for k in 1:npts
                vvband[iband, 1, 1, k] = vb1[k] * vb1[k]
                vvband[iband, 1, 2, k] = vb1[k] * vb2[k]
                vvband[iband, 1, 3, k] = vb1[k] * vb3[k]
                vvband[iband, 2, 1, k] = vb2[k] * vb1[k]
                vvband[iband, 2, 2, k] = vb2[k] * vb2[k]
                vvband[iband, 2, 3, k] = vb2[k] * vb3[k]
                vvband[iband, 3, 1, k] = vb3[k] * vb1[k]
                vvband[iband, 3, 2, k] = vb3[k] * vb2[k]
                vvband[iband, 3, 3, k] = vb3[k] * vb3[k]
            end
        end
    end
end

#=
    getBTPbands(coeffs, equivalences, lattvec; compute_velocity=true)

Reconstruct all bands using FFT.

Uses function barrier pattern for performance optimization:
1. Pre-computes Cartesian coordinates (lattvec * R) outside hot loop
2. Fills all 4 grids (energy + 3 velocity) in a single loop
3. Calls type-stable inner kernel

# Returns
- `eband`: Band energies (nbands × npts)
- `vvband`: Velocity outer products (nbands × 3 × 3 × npts) if compute_velocity=true

Note: lattvec is expected in Bohr (atomic units) and is converted to Ångström
internally to match Python BoltzTraP2's velocity calculation conventions.
=#
function getBTPbands(coeffs, equivalences, lattvec; compute_velocity=true)
    dims = determine_fft_dims(equivalences)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    d1, d2, d3 = dims

    # Pre-computation (amortized cost)
    plan = plan_ifft(zeros(ComplexF64, dims))

    # Convert lattvec from Bohr to Ångström to match Python BoltzTraP2
    lattvec_ang = lattvec * BOHR_TO_ANG

    # Pre-compute Cartesian coordinates: allvec_t[i, k] = (equiv * lattvec)[k, i]
    all_equiv_mat = vcat([Matrix(equiv) for equiv in equivalences]...)
    allvec_t = Matrix{Float64}(permutedims(all_equiv_mat * lattvec_ang))

    # Pre-compute CartesianIndex array for grid filling
    total_pts = sum(size(e, 1) for e in equivalences)
    all_indices = Vector{CartesianIndex{3}}(undef, total_pts)
    k = 1
    for equiv in equivalences
        for ir in 1:size(equiv, 1)
            r1 = equiv[ir, 1]
            r2 = equiv[ir, 2]
            r3 = equiv[ir, 3]
            all_indices[k] = CartesianIndex(mod(r1, d1) + 1, mod(r2, d2) + 1, mod(r3, d3) + 1)
            k += 1
        end
    end

    # Output arrays
    eband = zeros(nbands, npts)
    vvband = compute_velocity ? zeros(nbands, 3, 3, npts) : nothing

    # Call type-stable inner kernel
    _getBTPbands_kernel!(eband, vvband, coeffs, equivalences, allvec_t, all_indices, dims, plan)

    return eband, vvband
end

#=
    getBTPbands_parallel(coeffs, equivalences, lattvec; compute_velocity=true)

Parallel version of getBTPbands using Julia threads.

Uses the same optimization strategy as getBTPbands:
1. Pre-computes Cartesian coordinates (lattvec * R)
2. Each thread gets its own work arrays
3. Fills all 4 grids in a single loop per band
=#
function getBTPbands_parallel(coeffs, equivalences, lattvec; compute_velocity=true)
    dims = determine_fft_dims(equivalences)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    d1, d2, d3 = dims
    scale = prod(dims)

    # Pre-computation (amortized cost, shared across threads)
    plan = plan_ifft(zeros(ComplexF64, dims))
    lattvec_ang = lattvec * BOHR_TO_ANG
    all_equiv_mat = vcat([Matrix(equiv) for equiv in equivalences]...)
    allvec_t = Matrix{Float64}(permutedims(all_equiv_mat * lattvec_ang))

    total_pts = sum(size(e, 1) for e in equivalences)
    all_indices = Vector{CartesianIndex{3}}(undef, total_pts)
    k = 1
    for equiv in equivalences
        for ir in 1:size(equiv, 1)
            r1 = equiv[ir, 1]
            r2 = equiv[ir, 2]
            r3 = equiv[ir, 3]
            all_indices[k] = CartesianIndex(mod(r1, d1) + 1, mod(r2, d2) + 1, mod(r3, d3) + 1)
            k += 1
        end
    end

    # Output arrays
    eband = zeros(nbands, npts)
    vvband = compute_velocity ? zeros(nbands, 3, 3, npts) : nothing

    # Parallel processing with thread-local work arrays
    Threads.@threads for iband in 1:nbands
        # Thread-local work arrays
        egrid = zeros(ComplexF64, dims)
        vgrid1 = zeros(ComplexF64, dims)
        vgrid2 = zeros(ComplexF64, dims)
        vgrid3 = zeros(ComplexF64, dims)

        flat_coeffs = _flatten_coeffs(coeffs[iband, :], equivalences)
        _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

        eband[iband, :] = scale .* real.(vec(plan * egrid))

        if compute_velocity
            vb1 = scale .* real.(vec(plan * vgrid1))
            vb2 = scale .* real.(vec(plan * vgrid2))
            vb3 = scale .* real.(vec(plan * vgrid3))

            @inbounds for k in 1:npts
                vvband[iband, 1, 1, k] = vb1[k] * vb1[k]
                vvband[iband, 1, 2, k] = vb1[k] * vb2[k]
                vvband[iband, 1, 3, k] = vb1[k] * vb3[k]
                vvband[iband, 2, 1, k] = vb2[k] * vb1[k]
                vvband[iband, 2, 2, k] = vb2[k] * vb2[k]
                vvband[iband, 2, 3, k] = vb2[k] * vb3[k]
                vvband[iband, 3, 1, k] = vb3[k] * vb1[k]
                vvband[iband, 3, 2, k] = vb3[k] * vb2[k]
                vvband[iband, 3, 3, k] = vb3[k] * vb3[k]
            end
        end
    end

    return eband, vvband
end
