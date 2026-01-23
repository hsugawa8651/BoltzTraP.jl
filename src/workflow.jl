# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Logging

#= 
    get_equivalences(lattvec, positions, types, magmom, nkpt_target; symprec=1e-5)

Compute equivalence classes for approximately `nkpt_target` k-points.

Convenience function combining:
- `calc_nrotations` - get number of symmetry operations
- `compute_radius` - estimate sphere radius
- `compute_bounds` - compute search bounds
- `calc_sphere_quotient_set` - build equivalence classes

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `positions`: Atomic positions (3×natoms or natoms×3)
- `types`: Vector of atom type indices
- `magmom`: Magnetic moments (`nothing` for unpolarized)
- `nkpt_target`: Target number of equivalence classes

# Returns
- `equivalences`: Vector of equivalence class matrices
- `radius`: Actual radius used
- `nrotations`: Number of symmetry operations
=#
function get_equivalences(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom,
    nkpt_target::Integer;
    symprec::Real = 1e-5,
)
    # Ensure positions are natoms×3
    if size(positions, 1) == 3 && size(positions, 2) != 3
        positions = positions'
    end

    # Convert to Vector of SVectors for internal functions
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]

    # Get number of rotations
    nrot = calc_nrotations(lattvec, pos_vec, types, magmom; symprec)

    # Compute radius for target number of equivalences
    radius = compute_radius(lattvec, nrot, nkpt_target)

    # Compute bounds
    bounds = compute_bounds(lattvec, radius)

    # Compute equivalence classes
    equivalences =
        calc_sphere_quotient_set(lattvec, pos_vec, types, magmom, radius, bounds; symprec)

    return equivalences, radius, nrot
end

"""
    run_interpolate(data::DFTData{1}; kwargs...) -> InterpolationResult

Run interpolation workflow on [`DFTData`](@ref) from loaders.

# Arguments

  - `data`: [`DFTData`](@ref) from loaders ([`load_vasp`](@ref), [`load_qe`](@ref), [`load_wien2k`](@ref), [`load_gene`](@ref), [`load_abinit`](@ref), [`load_dftk`](@ref))

See `run_interpolate(::NamedTuple; kwargs...)` for keyword arguments.

# Returns

[`InterpolationResult`](@ref) containing coefficients, equivalences, and metadata.
"""
function run_interpolate(
    data::DFTData{1};
    source::String = "unknown",
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
)
    # Convert DFTData to NamedTuple for existing implementation
    data_nt = (
        lattice = data.lattice,
        positions = data.positions,
        species = data.species,
        kpoints = data.kpoints,
        weights = data.weights,
        ebands = data.ebands,
        occupations = data.occupations,
        fermi = data.fermi,
        nelect = data.nelect,
    )
    return run_interpolate(
        data_nt;
        source,
        output,
        kpoints,
        multiplier,
        emin,
        emax,
        absolute,
        verbose,
        symprec,
    )
end

"""
    run_interpolate(data::DFTData{2}; kwargs...) -> InterpolationResult

Run interpolation workflow on spin-polarized (collinear) [`DFTData`](@ref).

Spin channels are concatenated following Python BoltzTraP2 convention:
`ebands[:,:,1]` and `ebands[:,:,2]` are vertically concatenated to form
`(2*nbands, nkpts)` before interpolation.

# Arguments
- `data`: [`DFTData{2}`](@ref) from loaders with spin-polarized data

See `run_interpolate(::NamedTuple; kwargs...)` for keyword arguments.

# Returns
[`InterpolationResult`](@ref) with `metadata["spintype"] = "Collinear"` and
`metadata["dosweight"] = 1.0`.
"""
function run_interpolate(
    data::DFTData{2};
    source::String = "unknown",
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
)
    # Concatenate spin channels (Python BoltzTraP2 convention)
    # ebands: (nbands, nkpts, 2) → (2*nbands, nkpts, 1)
    ebands_spin1 = data.ebands[:, :, 1]
    ebands_spin2 = data.ebands[:, :, 2]
    ebands_concat = vcat(ebands_spin1, ebands_spin2)

    # Reshape to 3D for compatibility with NamedTuple interface
    ebands_3d = reshape(ebands_concat, size(ebands_concat, 1), size(ebands_concat, 2), 1)

    # Convert DFTData to NamedTuple with concatenated bands and magmom
    data_nt = (
        lattice = data.lattice,
        positions = data.positions,
        species = data.species,
        kpoints = data.kpoints,
        weights = data.weights,
        ebands = ebands_3d,
        occupations = data.occupations,  # Note: occupations not concatenated (unused in interpolation)
        fermi = data.fermi,
        nelect = data.nelect,
        magmom = data.magmom,  # Pass magmom for symmetry
    )

    return run_interpolate(
        data_nt;
        source,
        output,
        kpoints,
        multiplier,
        emin,
        emax,
        absolute,
        verbose,
        symprec,
    )
end

#=
    run_interpolate(data::DFTData{N}; kwargs...) where N

Catch-all method for unsupported spin configurations (N > 2).
=#
function run_interpolate(data::DFTData{N}; kwargs...) where {N}
    if N > 2
        error(
            "Non-collinear calculations (nspin=$N) are not supported.\n" *
            "Only DFTData{1} (unpolarized) and DFTData{2} (collinear) are supported."
        )
    end
    # This should not be reached for N=1 or N=2
    error("Unexpected nspin=$N in run_interpolate")
end

"""
    run_interpolate(data::NamedTuple; kwargs...) -> InterpolationResult

Run interpolation workflow on pre-loaded band structure data (legacy interface).

# Arguments

  - `data`: NamedTuple with band structure data (for backward compatibility)

# Keyword Arguments

  - `source="unknown"`: Data source description for metadata
  - `output=nothing`: Output file path (no file saved if `nothing`)
  - `kpoints=5000`: Target number of k-points/equivalences
  - `multiplier=nothing`: Enhancement factor for k-points (alternative to `kpoints`)
  - `emin=-Inf`: Minimum energy relative to Fermi in Ha (no filter if `-Inf`)
  - `emax=+Inf`: Maximum energy relative to Fermi in Ha (no filter if `+Inf`)
  - `absolute=false`: Interpret emin/emax as absolute energies
  - `verbose=false`: Print progress information
  - `symprec=1e-5`: Symmetry precision

# Returns

[`InterpolationResult`](@ref) containing coefficients, equivalences, and metadata.

# Example

```julia
# From VASP
data = load_vasp("./Si.vasp")
result = run_interpolate(data; source = "VASP", kpoints = 5000)

# From QE
data = load_qe("./Si.qe")
result = run_interpolate(data; source = "QE")

# From DFTK (requires DFTK.jl loaded)
data = load_dftk(scfres)
result = run_interpolate(data; source = "DFTK")
```

See also: [`load_vasp`](@ref), [`load_qe`](@ref), [`load_dftk`](@ref)
"""
function run_interpolate(
    data::NamedTuple;
    source::String = "unknown",
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
)
    # Log received arguments for debugging
    @debug "run_interpolate called" source output kpoints multiplier emin emax absolute verbose symprec

    # Validate arguments
    if isnothing(kpoints) && isnothing(multiplier)
        kpoints = 5000  # Default
    elseif !isnothing(kpoints) && !isnothing(multiplier)
        error("Cannot specify both kpoints and multiplier")
    end

    @debug "run_interpolate after defaults" kpoints multiplier

    verbose && println("Processing data from $source...")

    # 1. Prepare structure data
    lattvec = data.lattice
    positions = data.positions  # 3×natoms, need to transpose
    species = data.species
    types = [findfirst(==(s), unique(species)) for s in species]
    magmom = get(data, :magmom, nothing)  # For collinear magnetic calculations

    # 2. Determine target k-points
    if !isnothing(multiplier)
        nkpt_target = size(data.kpoints, 2) * multiplier
        verbose && println("Using multiplier=$multiplier → nkpt_target=$nkpt_target")
    else
        nkpt_target = kpoints
    end

    # 3. Get space group info
    sginfo = get_spacegroup_info(lattvec, positions', types; symprec)
    verbose && println(
        "  Space group: $(sginfo.spacegroup_number) ($(sginfo.international_symbol))",
    )

    # 4. Compute equivalences
    verbose && println("Computing equivalences for ~$nkpt_target k-points...")
    equivalences, radius, nrot = get_equivalences(
        lattvec, positions', types, magmom, nkpt_target; symprec
    )
    verbose && println("  Rotations: $nrot")
    verbose && println("  Radius: $(round(radius, digits=2))")
    verbose && println("  Equivalences: $(length(equivalences))")

    # 4. Prepare band data (already in atomic units)
    kpts = data.kpoints'  # nk×3
    ebands_raw = data.ebands[:, :, 1]  # nbands×nk (spin 1), in Ha
    fermi = data.fermi  # in Ha
    nbands_total = size(ebands_raw, 1)

    # Determine dosweight from spin polarization
    # Use magmom to detect collinear (ebands may be concatenated for collinear)
    dosweight = isnothing(magmom) ? 2.0 : 1.0

    # 5. Filter bands by energy
    # emin/emax are in Ha (relative to Fermi unless absolute=true)
    if absolute
        emin_abs = emin
        emax_abs = emax
    else
        emin_abs = fermi + emin
        emax_abs = fermi + emax
    end

    # Find bands within energy range
    band_min = minimum(ebands_raw, dims = 2)[:]
    band_max = maximum(ebands_raw, dims = 2)[:]
    band_mask = (band_max .>= emin_abs) .& (band_min .<= emax_abs)
    selected_bands = findall(band_mask)

    if isempty(selected_bands)
        error("No bands found in energy range [$emin_abs, $emax_abs] Ha")
    end

    ebands = ebands_raw[selected_bands, :]
    verbose && println("  Bands: $(length(selected_bands))/$nbands_total selected")
    verbose && println(
        "  Energy range: [$(round(emin_abs, digits=6)), $(round(emax_abs, digits=6))] Ha",
    )

    # 6. Convert equivalences to matrix format for FourierInterpolator
    equiv_matrices = [hcat(eq...)' for eq in equivalences]

    # 7. Create interpolator
    verbose && println("Fitting Fourier coefficients...")
    interp = FourierInterpolator(kpts, ebands, equiv_matrices, lattvec)

    # 8. Create result with metadata
    atoms = Dict{String,Any}(
        "species" => species,
        "positions" => collect(positions'),
        "lattice" => collect(lattvec),
    )
    # Determine spintype from magmom
    spintype = isnothing(magmom) ? "Unpolarized" : "Collinear"

    metadata = Dict{String,Any}(
        "fermi" => fermi,
        "nelect" => data.nelect,
        "dosweight" => dosweight,
        "spintype" => spintype,
        "selected_bands" => selected_bands,
        "nkpt_original" => size(data.kpoints, 2),
        "nkpt_target" => nkpt_target,
        "radius" => radius,
        "nrotations" => nrot,
        "emin" => emin,
        "emax" => emax,
        "absolute" => absolute,
        "source" => source,
        "spacegroup_number" => sginfo.spacegroup_number,
        "spacegroup_symbol" => sginfo.international_symbol,
    )

    result = InterpolationResult(interp; atoms = atoms, metadata = metadata)

    # 9. Save if output specified
    if !isnothing(output)
        verbose && println("Saving to $output...")
        save_interpolation(output, result)
    end

    verbose && println("Done.")
    return result
end

"""
    run_interpolate(directory::String; kwargs...) -> InterpolationResult

Run the complete interpolation workflow on VASP data.

# Arguments

  - `directory`: Path to directory containing vasprun.xml

# Keyword Arguments

  - `output=nothing`: Output file path (no file saved if `nothing`)
  - `kpoints=5000`: Target number of k-points/equivalences
  - `multiplier=nothing`: Enhancement factor for k-points (alternative to `kpoints`)
  - `emin=-Inf`: Minimum energy relative to Fermi in Ha (no filter if `-Inf`)
  - `emax=+Inf`: Maximum energy relative to Fermi in Ha (no filter if `+Inf`)
  - `absolute=false`: Interpret emin/emax as absolute energies
  - `verbose=false`: Print progress information
  - `symprec=1e-5`: Symmetry precision

# Returns

[`InterpolationResult`](@ref) containing coefficients, equivalences, and metadata.

# Example

```julia
result = run_interpolate("./Si.vasp"; kpoints = 5000, verbose = true)
```

See also: [`run_integrate`](@ref), [`save_interpolation`](@ref)
"""
function run_interpolate(
    directory::String;
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
)
    # Log CLI arguments
    @debug "run_interpolate(directory) called" directory output kpoints multiplier emin emax absolute verbose symprec

    # Check directory exists
    if !isdir(directory)
        @debug "Directory not found, trying as file path" directory
    end

    # Load VASP data
    verbose && println("Loading VASP data from $directory...")
    data = load_vasp(directory)

    @debug "VASP data loaded" nkpts=size(data.kpoints, 2) nbands=size(data.ebands, 1) fermi=data.fermi

    return run_interpolate(
        data;
        source = directory,
        output,
        kpoints,
        multiplier,
        emin,
        emax,
        absolute,
        verbose,
        symprec,
    )
end

"""
    run_integrate(interp::InterpolationResult; kwargs...) -> TransportResult
    run_integrate(interp::InterpolationResult, mur; kwargs...) -> TransportResult
    run_integrate(input::String; kwargs...) -> TransportResult
    run_integrate(input::String, mur; kwargs...) -> TransportResult

Compute transport coefficients from interpolation result.

# Arguments

  - `interp`: [`InterpolationResult`](@ref) from [`run_interpolate`](@ref)
  - `input`: Path to interpolation result file (.jld2) - alternative to `interp`
  - `mur`: Chemical potential values in Hartree (positional argument, optional).
    Same name as Python BoltzTraP2 for compatibility.

# Keyword Arguments

  - `temperatures=[300.0]`: Vector of temperatures in K
  - `output=nothing`: Output file path (no file saved if `nothing`)
  - `bins=0`: Number of DOS histogram bins (auto if `0`)
  - `scissor=nothing`: Target band gap in eV for scissor correction.
    If specified, conduction bands are shifted to achieve this gap.
    Only applicable to semiconductors with a band gap.
  - `verbose=false`: Print progress information

# Returns

[`TransportResult`](@ref) containing σ, S, κ tensors, DOS, and metadata.

# Examples

```julia
# Direct from InterpolationResult (μ auto-generated)
interp = run_interpolate("./Si.vasp")
transport = run_integrate(interp; temperatures = [300.0])

# From file
transport = run_integrate("si_interp.jld2"; temperatures = [200.0, 300.0, 400.0])

# With explicit μ grid (positional argument, same as Python BoltzTraP2)
mur = range(-0.5, 0.5, length = 100) .* EV_TO_HA  # eV to Ha
transport = run_integrate(interp, mur; temperatures = [300.0])

# With scissor correction (shift conduction bands to 1.1 eV gap)
transport = run_integrate(interp; temperatures = [300.0], scissor = 1.1)
```

See also: [`run_interpolate`](@ref), [`save_integrate`](@ref), [`apply_scissor`](@ref)
"""
function run_integrate(
    interp::InterpolationResult;
    temperatures::AbstractVector{<:Real} = [300.0],
    output::Union{String,Nothing} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)
    # Log received arguments
    @debug "run_integrate called" temperatures output bins scissor verbose

    # Extract metadata
    # Note: fermi is stored in Ha (atomic units)
    nelect = get(interp.metadata, "nelect", 0.0)
    dosweight = get(interp.metadata, "dosweight", 2.0)
    fermi_dft = get(interp.metadata, "fermi", 0.0)  # in Ha
    source = get(interp.metadata, "source_file", "unknown")
    spintype = get(interp.metadata, "spintype", "Unpolarized")  # v0.1 compat

    if nelect == 0.0
        error("nelect not found in interpolation metadata")
    end

    verbose && println("  nelect: $nelect")
    verbose && println("  dosweight: $dosweight")
    verbose && println("  Fermi (from DFT): $(round(fermi_dft, digits=6)) Ha")

    # 1. Reconstruct bands on FFT grid
    verbose && println("Reconstructing bands via FFT...")
    eband, vvband = getBTPbands(interp.coeffs, interp.equivalences, interp.lattvec)
    nbands, npts = size(eband)
    verbose && println("  Bands: $nbands, FFT points: $npts")
    verbose && println(
        "  Energy range: [$(round(minimum(eband), digits=4)), $(round(maximum(eband), digits=4))] Ha",
    )

    # 1.5. Apply scissor correction if requested
    scissor_Ha = nothing
    if !isnothing(scissor)
        scissor_Ha = scissor * EV_TO_HA  # Convert eV to Ha
        verbose && println("Applying scissor correction (target gap: $scissor eV)...")

        # Compute temporary DOS to find current gap
        npts_dos_temp = bins > 0 ? bins : 500
        epsilon_temp, dos_temp, _ = BTPDOS(eband, vvband; npts = npts_dos_temp)

        # Apply scissor to shift conduction bands
        eband = apply_scissor(epsilon_temp, dos_temp, nelect, eband, scissor_Ha; dosweight)
        verbose && println("  Scissor applied: target gap = $scissor eV")
    end

    # 2. Compute DOS and transport DOS
    npts_dos = bins > 0 ? bins : 500
    verbose && println("Computing DOS (bins=$npts_dos)...")
    epsilon, dos, vvdos = BTPDOS(eband, vvband; npts = npts_dos)
    verbose && println(
        "  DOS energy range: [$(round(epsilon[1], digits=4)), $(round(epsilon[end], digits=4))] Ha",
    )

    # 3. Determine μ range (auto-generate from DOS grid)
    # margin = 9 * kB * T_max (ensures ~1e-4 accuracy at edges)
    kB_Ha_K = KB_AU
    margin = 9.0 * kB_Ha_K * maximum(temperatures)
    μ_min = epsilon[1] + margin
    μ_max = epsilon[end] - margin

    if μ_min >= μ_max
        error("Energy window too narrow for requested temperatures")
    end

    μ_indices = findall(e -> e > μ_min && e < μ_max, epsilon)
    μ_range = epsilon[μ_indices]
    verbose && println(
        "  μ range (auto): $(length(μ_range)) points in [$(round(μ_min, digits=4)), $(round(μ_max, digits=4))] Ha",
    )

    # 4. Solve for intrinsic chemical potential at each temperature
    verbose && println("Computing intrinsic μ for each temperature...")
    Tr = collect(Float64, temperatures)
    nT = length(Tr)
    μ0 = zeros(nT)
    for (iT, T) in enumerate(Tr)
        μ0[iT] = solve_for_mu(
            epsilon,
            dos,
            nelect,
            T;
            dosweight,
            refine = true,
            try_center = true,
        )
        verbose && println("  T=$(T)K: μ0=$(round(μ0[iT], digits=6)) Ha")
    end

    # Refined Fermi level (T=0)
    fermi_Ha =
        solve_for_mu(epsilon, dos, nelect, 0.0; dosweight, refine = true, try_center = true)
    verbose && println("  Refined Fermi: $(round(fermi_Ha, digits=6)) Ha")

    # 5. Compute Fermi integrals
    verbose && println("Computing Fermi integrals...")
    N, L0, L1, L2 = fermi_integrals(epsilon, dos, vvdos, μ_range, Tr; dosweight)
    verbose && println("  L0 shape: $(size(L0))")

    # 6. Calculate Onsager coefficients
    verbose && println("Computing Onsager coefficients...")
    # Convert lattvec from Bohr to Ångström for volume calculation
    # This matches Python BoltzTraP2's unit conventions
    lattvec_ang = interp.lattvec * BOHR_TO_ANG  # BOHR_TO_ANG from units.jl
    vuc = abs(det(lattvec_ang))  # Unit cell volume in Ų
    σ, S, κ = calc_onsager_coefficients(L0, L1, L2, Tr, vuc)
    verbose && println("  σ shape: $(size(σ))")

    # 7. Convert μ_range from Ha to eV for output (HA_TO_EV from units.jl)
    μ_range_eV = μ_range .* HA_TO_EV

    # 8. Create result
    result_metadata = Dict{String,Any}(
        "source" => source,
        "nelect" => nelect,
        "dosweight" => dosweight,
        "spintype" => spintype,
        "fermi_Ha" => fermi_Ha,
        "fermi_eV" => fermi_Ha * HA_TO_EV,
        "fermi_dft_Ha" => fermi_dft,
        "fermi_dft_eV" => fermi_dft * HA_TO_EV,
        "mu0_Ha" => μ0,
        "mu0_eV" => μ0 .* HA_TO_EV,
        "vuc_ang3" => vuc,
        "nbands" => nbands,
        "npts_fft" => npts,
        "npts_dos" => npts_dos,
    )
    # Add scissor info if applied
    if !isnothing(scissor)
        result_metadata["scissor_eV"] = scissor
        result_metadata["scissor_Ha"] = scissor_Ha
    end

    # Reshape tensors to (3, 3, nT, nμ) for TransportResult
    nμ = length(μ_range)
    σ_out = permutedims(σ, (3, 4, 1, 2))  # (nT, nμ, 3, 3) -> (3, 3, nT, nμ)
    S_out = permutedims(S, (3, 4, 1, 2))
    κ_out = permutedims(κ, (3, 4, 1, 2))

    dos_info = Dict{String,Any}(
        "epsilon_Ha" => epsilon,
        "epsilon_eV" => epsilon .* HA_TO_EV,
        "dos" => dos,
    )

    result = TransportResult(Tr, μ_range_eV, σ_out, S_out, κ_out, dos_info, result_metadata)

    # 9. Save if output specified
    if !isnothing(output)
        verbose && println("Saving to $output...")
        save_integrate(output, result)
    end

    verbose && println("Done.")
    return result
end

# Method with explicit μ grid (positional argument)
function run_integrate(
    interp::InterpolationResult,
    mur::AbstractVector{<:Real};
    temperatures::AbstractVector{<:Real} = [300.0],
    output::Union{String,Nothing} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)
    @debug "run_integrate(interp, mur) called" temperatures output bins scissor verbose

    # Extract metadata
    nelect = get(interp.metadata, "nelect", 0.0)
    dosweight = get(interp.metadata, "dosweight", 2.0)
    fermi_dft = get(interp.metadata, "fermi", 0.0)
    source = get(interp.metadata, "source_file", "unknown")
    spintype = get(interp.metadata, "spintype", "Unpolarized")  # v0.1 compat

    if nelect == 0.0
        error("nelect not found in interpolation metadata")
    end

    verbose && println("  nelect: $nelect")
    verbose && println("  dosweight: $dosweight")
    verbose && println("  Fermi (from DFT): $(round(fermi_dft, digits=6)) Ha")

    # 1. Reconstruct bands on FFT grid
    verbose && println("Reconstructing bands via FFT...")
    eband, vvband = getBTPbands(interp.coeffs, interp.equivalences, interp.lattvec)
    nbands, npts = size(eband)
    verbose && println("  Bands: $nbands, FFT points: $npts")

    # 1.5. Apply scissor correction if requested
    scissor_Ha = nothing
    if !isnothing(scissor)
        scissor_Ha = scissor * EV_TO_HA  # Convert eV to Ha
        verbose && println("Applying scissor correction (target gap: $scissor eV)...")

        # Compute temporary DOS to find current gap
        npts_dos_temp = bins > 0 ? bins : 500
        epsilon_temp, dos_temp, _ = BTPDOS(eband, vvband; npts = npts_dos_temp)

        # Apply scissor to shift conduction bands
        eband = apply_scissor(epsilon_temp, dos_temp, nelect, eband, scissor_Ha; dosweight)
        verbose && println("  Scissor applied: target gap = $scissor eV")
    end

    # 2. Compute DOS and transport DOS
    npts_dos = bins > 0 ? bins : 500
    verbose && println("Computing DOS (bins=$npts_dos)...")
    epsilon, dos, vvdos = BTPDOS(eband, vvband; npts = npts_dos)

    # 3. Use provided μ range (same as Python BoltzTraP2 `mur` argument)
    μ_range = collect(Float64, mur)
    verbose && println(
        "  μ range (provided): $(length(μ_range)) points in [$(round(minimum(μ_range), digits=4)), $(round(maximum(μ_range), digits=4))] Ha",
    )

    # 4. Solve for intrinsic chemical potential at each temperature
    verbose && println("Computing intrinsic μ for each temperature...")
    Tr = collect(Float64, temperatures)
    nT = length(Tr)
    μ0 = zeros(nT)
    for (iT, T) in enumerate(Tr)
        μ0[iT] = solve_for_mu(
            epsilon,
            dos,
            nelect,
            T;
            dosweight,
            refine = true,
            try_center = true,
        )
        verbose && println("  T=$(T)K: μ0=$(round(μ0[iT], digits=6)) Ha")
    end

    # Refined Fermi level (T=0)
    fermi_Ha =
        solve_for_mu(epsilon, dos, nelect, 0.0; dosweight, refine = true, try_center = true)
    verbose && println("  Refined Fermi: $(round(fermi_Ha, digits=6)) Ha")

    # 5. Compute Fermi integrals
    verbose && println("Computing Fermi integrals...")
    N, L0, L1, L2 = fermi_integrals(epsilon, dos, vvdos, μ_range, Tr; dosweight)

    # 6. Calculate Onsager coefficients
    verbose && println("Computing Onsager coefficients...")
    # Convert lattvec from Bohr to Ångström for volume calculation
    # This matches Python BoltzTraP2's unit conventions
    lattvec_ang = interp.lattvec * BOHR_TO_ANG  # BOHR_TO_ANG from units.jl
    vuc = abs(det(lattvec_ang))  # Unit cell volume in Ų
    σ, S, κ = calc_onsager_coefficients(L0, L1, L2, Tr, vuc)

    # 7. Convert μ_range from Ha to eV for output (HA_TO_EV from units.jl)
    μ_range_eV = μ_range .* HA_TO_EV

    # 8. Create result
    result_metadata = Dict{String,Any}(
        "source" => source,
        "nelect" => nelect,
        "dosweight" => dosweight,
        "spintype" => spintype,
        "fermi_Ha" => fermi_Ha,
        "fermi_eV" => fermi_Ha * HA_TO_EV,
        "fermi_dft_Ha" => fermi_dft,
        "fermi_dft_eV" => fermi_dft * HA_TO_EV,
        "mu0_Ha" => μ0,
        "mu0_eV" => μ0 .* HA_TO_EV,
        "vuc_ang3" => vuc,
        "nbands" => nbands,
        "npts_fft" => npts,
        "npts_dos" => npts_dos,
    )
    # Add scissor info if applied
    if !isnothing(scissor)
        result_metadata["scissor_eV"] = scissor
        result_metadata["scissor_Ha"] = scissor_Ha
    end

    # Reshape tensors to (3, 3, nT, nμ)
    nμ = length(μ_range)
    σ_out = permutedims(σ, (3, 4, 1, 2))
    S_out = permutedims(S, (3, 4, 1, 2))
    κ_out = permutedims(κ, (3, 4, 1, 2))

    dos_info = Dict{String,Any}(
        "epsilon_Ha" => epsilon,
        "epsilon_eV" => epsilon .* HA_TO_EV,
        "dos" => dos,
    )

    result = TransportResult(Tr, μ_range_eV, σ_out, S_out, κ_out, dos_info, result_metadata)

    # 9. Save if output specified
    if !isnothing(output)
        verbose && println("Saving to $output...")
        save_integrate(output, result)
    end

    verbose && println("Done.")
    return result
end

# File-based convenience method (μ auto-generated)
function run_integrate(
    input::String;
    temperatures::AbstractVector{<:Real} = [300.0],
    output::Union{String,Nothing} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)
    @debug "run_integrate(input) called" input temperatures output bins scissor verbose

    # Check file exists
    if !isfile(input)
        @debug "Input file not found" input
        error("Input file not found: $input")
    end

    verbose && println("Loading interpolation data from $input...")
    interp = load_interpolation(input)

    @debug "Interpolation data loaded" nbands=size(interp.coeffs, 1) nequiv=length(
        interp.equivalences,
    )

    return run_integrate(interp; temperatures, output, bins, scissor, verbose)
end

# File-based convenience method with explicit μ grid
function run_integrate(
    input::String,
    mur::AbstractVector{<:Real};
    temperatures::AbstractVector{<:Real} = [300.0],
    output::Union{String,Nothing} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)
    @debug "run_integrate(input, mur) called" input temperatures output bins scissor verbose

    if !isfile(input)
        error("Input file not found: $input")
    end

    verbose && println("Loading interpolation data from $input...")
    interp = load_interpolation(input)

    return run_integrate(interp, mur; temperatures, output, bins, scissor, verbose)
end

# ============================================================================ 
# Electron count and chemical potential functions
# ============================================================================ 

#= 
    calc_N(epsilon, dos, μ, T; dosweight=2.0)

Compute the electron count by integrating over the DOS.

# Arguments
- `epsilon`: Array of energies at which the DOS is available [Ha]
- `dos`: Density of states
- `μ`: Chemical potential [Ha]
- `T`: Temperature [K]
- `dosweight`: Maximum occupancy of an electron mode (2.0 for non-spin-polarized)

# Returns
Electron count (negative value for comparison with N0).
=#
function calc_N(
    epsilon::AbstractVector,
    dos::AbstractVector,
    μ::Real,
    T::Real;
    dosweight::Real = 2.0,
)
    if T == 0.0
        # Zero temperature: step function
        occ = ifelse.(epsilon .< μ, 1.0, 0.0)
        # Handle exact equality at Fermi level
        for i in eachindex(epsilon)
            if epsilon[i] ≈ μ
                occ[i] = 0.5
            end
        end
    else
        # Finite temperature: Fermi-Dirac distribution
        kBT = T * KB_AU
        occ = fermi_dirac(epsilon, μ, kBT)
    end
    de = epsilon[2] - epsilon[1]
    return -dosweight * sum(dos .* occ) * de
end

#= 
    solve_for_mu(epsilon, dos, N0, T; dosweight=2.0, refine=false, try_center=false)

Estimate the chemical potential required to have N0 electrons.

Uses binary search since N(μ) is monotonically increasing with μ.
If μ falls in a wide gap (relative to kB * T), then μ is moved to
the center of the gap.

# Arguments
- `epsilon`: Array of energies at which the DOS is available [Ha]
- `dos`: Density of states
- `N0`: Number of valence electrons in the compound
- `T`: Temperature [K]
- `dosweight`: Maximum occupancy of an electron mode
- `refine`: If `true`, use root finding for exact μ
- `try_center`: If `true`, snap μ to center of large gaps

# Returns
An estimate of the intrinsic chemical potential at a given temperature [Ha].
=#
function solve_for_mu(
    epsilon::AbstractVector,
    dos::AbstractVector,
    N0::Real,
    T::Real;
    dosweight::Real = 2.0,
    refine::Bool = false,
    try_center::Bool = false,
)
    n = length(epsilon)

    # 1. Binary search to find index where calc_N(ε) + N0 ≈ 0
    # N(μ) is monotonically increasing, so calc_N(μ) + N0 goes from + to -
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        residual = calc_N(epsilon, dos, epsilon[mid], T; dosweight) + N0
        if residual > 0
            lo = mid  # Need higher μ
        else
            hi = mid  # Need lower μ
        end
    end

    # Choose the index with smaller |residual|
    res_lo = abs(calc_N(epsilon, dos, epsilon[lo], T; dosweight) + N0)
    res_hi = abs(calc_N(epsilon, dos, epsilon[hi], T; dosweight) + N0)
    pos = res_lo < res_hi ? lo : hi
    μ = epsilon[pos]
    center = false
    lepsilon, hepsilon = 0.0, 0.0  # Will be set if in gap

    # 2. Check if μ falls in a gap
    if dos[pos] == 0.0
        # Find gap edges
        lpos = findlast(i -> dos[i] != 0.0, 1:pos)
        hpos_offset = findfirst(i -> dos[i] != 0.0, pos:n)
        hpos = isnothing(hpos_offset) ? nothing : pos + hpos_offset - 1

        if isnothing(lpos) || isnothing(hpos)
            error("μ lies outside the range of band energies")
        end

        lepsilon = epsilon[lpos]
        hepsilon = epsilon[hpos]
        kBT = T * KB_AU

        # If μ is in a gap and far enough from edges, move to center
        if try_center &&
           kBT > 0 &&
           min(hepsilon - μ, μ - lepsilon) >= FD_XMAX_GAP * kBT / 2.0
            pos = round(Int, 0.5 * (lpos + hpos))
            μ = epsilon[pos]
            center = true
        end
    end

    # 3. Refinement using root finding
    if refine
        if center
            # Exact center of gap
            μ = 0.5 * (lepsilon + hepsilon)
        else
            # Find root of calc_N(μ) + N0 = 0 using bisection
            lmu, hmu = epsilon[lo], epsilon[hi]
            μ = _bisection_root(
                μ_arg -> calc_N(epsilon, dos, μ_arg, T; dosweight) + N0,
                lmu,
                hmu,
            )
        end
    end

    return μ
end

#= 
    _bisection_root(f, a, b; tol=1e-12, maxiter=100)

Find root of f in interval [a, b] using bisection method.
Return midpoint if root is not bracketed (fa * fb > 0).
=#
function _bisection_root(f, a::Real, b::Real; tol::Real = 1e-12, maxiter::Int = 100)
    fa, fb = f(a), f(b)

    # Check if root is bracketed
    if fa * fb > 0
        # Not bracketed, return midpoint
        return 0.5 * (a + b)
    end

    for _ = 1:maxiter
        c = 0.5 * (a + b)
        if (b - a) < tol
            return c
        end
        fc = f(c)
        if fc == 0
            return c
        elseif fa * fc < 0
            b, fb = c, fc
        else
            a, fa = c, fc
        end
    end

    return 0.5 * (a + b)
end

# ============================================================================
# Scissor Correction
# ============================================================================

"""
    apply_scissor(epsilon, dos, N0, eband, desired_gap; dosweight=2.0) -> Array

Shift all conduction bands to achieve the desired value of the band gap.

This function identifies bands that lie entirely above the Fermi level
(conduction bands) and shifts them by a constant amount to match the
specified band gap.

# Arguments

  - `epsilon`: energy grid from BTPDOS (npts,), in Hartree
  - `dos`: density of states (npts,)
  - `N0`: number of valence electrons
  - `eband`: band energies (nbands, nkpoints), in Hartree
  - `desired_gap`: target band gap in Hartree
  - `dosweight`: maximum occupancy of an electron mode (2.0 for non-spin-polarized)

# Returns

Array of shifted band energies with the same shape as `eband`.

# Errors

  - `ArgumentError`: if there is no gap at the Fermi level (material is metallic)
  - `ArgumentError`: if the Fermi energy lies outside the band energy range

# Example

```julia
# Shift Si bands to achieve 1.1 eV gap
desired_gap = 1.1 * EV_TO_HA
eband_shifted = apply_scissor(epsilon, dos, nelect, eband, desired_gap)
```

See also: [`solve_for_mu`](@ref)
"""
function apply_scissor(epsilon, dos, N0, eband, desired_gap; dosweight = 2.0)
    # 1. Find current Fermi energy
    fermi =
        solve_for_mu(epsilon, dos, N0, 0.0; dosweight, refine = false, try_center = false)

    # 2. Check if Fermi lies in a gap
    pos = argmin(abs.(epsilon .- fermi))
    if dos[pos] != 0.0
        throw(ArgumentError("No band gap found at Fermi level (material is metallic)"))
    end

    # 3. Find gap edges (VBM and CBM)
    lpos = -1
    hpos = -1

    # Search downward for VBM
    for i = pos:-1:1
        if dos[i] != 0.0
            lpos = i
            break
        end
    end

    # Search upward for CBM
    for i = pos:length(dos)
        if dos[i] != 0.0
            hpos = i
            break
        end
    end

    if lpos == -1 || hpos == -1
        throw(ArgumentError("Fermi energy lies outside the band energy range"))
    end

    # 4. Calculate shift amount
    current_gap = epsilon[hpos] - epsilon[lpos]
    delta = desired_gap - current_gap

    # 5. Identify conduction bands (all states above Fermi)
    # A band is a conduction band if its minimum energy across all k-points is above Fermi
    conduction_mask = vec(minimum(eband, dims = 2)) .> fermi

    # 6. Apply shift to conduction bands
    eband_shifted = copy(eband)
    eband_shifted[conduction_mask, :] .+= delta

    return eband_shifted
end

# Error dispatch for invalid argument types
function apply_scissor(x::Any)
    throw(
        ArgumentError(
            "apply_scissor requires (epsilon, dos, N0, eband, desired_gap), got $(typeof(x))",
        ),
    )
end

# ============================================================================
# Convenience methods for TransportResult
# ============================================================================

"""
    calc_cv(transport::TransportResult) -> Matrix{Float64}

Compute the electronic contribution to the heat capacity from transport results.

This is a convenience method that extracts all necessary data from the
[`TransportResult`](@ref) struct. The `epsilon`, `dos`, `temperatures`,
`mu_values`, and `dosweight` are automatically retrieved.

# Returns

Matrix of shape (nT, nμ) with heat capacity in SI units (J/K).

# Example

```julia
transport = run_integrate("si_interp.jld2"; temperatures = [300.0])
cv = calc_cv(transport)
```

See also: [`calc_cv(epsilon, dos, μ_range, T_range; dosweight)`](@ref calc_cv)
"""
function calc_cv(transport::TransportResult)
    if isnothing(transport.dos)
        error("TransportResult does not contain DOS data")
    end

    epsilon = transport.dos["epsilon_Ha"]
    dos = transport.dos["dos"]
    dosweight = get(transport.metadata, "dosweight", 2.0)

    return calc_cv(epsilon, dos, transport.mu_values, transport.temperatures; dosweight)
end

# Error dispatch for invalid single-argument calls
function calc_cv(x::Any)
    throw(ArgumentError("calc_cv requires a TransportResult, got $(typeof(x))"))
end
