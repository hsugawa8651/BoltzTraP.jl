# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Brillouin
using Plots

#=
    _generate_manual_kpath(kpath, lattvec, npoints) -> (kvecs, kdist, labels, positions)

Generate k-path from manual specification. Returns interpolated k-points,
cumulative distances, labels and label positions.
=#
function _generate_manual_kpath(kpath::NamedTuple, lattvec::Matrix{Float64}, npoints::Int)
    # Validate kpath structure
    if !haskey(kpath, :points) || !haskey(kpath, :paths)
        error("kpath must have :points and :paths fields")
    end

    points_dict = kpath.points
    path_segments = kpath.paths

    # Get reciprocal lattice for distance calculation
    reclattvec = 2π * inv(lattvec)'

    # Interpolate along the path
    kvecs = Vector{Vector{Float64}}()
    kdist = Float64[]
    path_labels = String[]
    label_positions = Float64[]

    current_dist = 0.0

    for segment in path_segments
        for i in 1:(length(segment)-1)
            label_start = segment[i]
            label_end = segment[i+1]

            k_start = Float64.(points_dict[label_start])
            k_end = Float64.(points_dict[label_end])

            # Segment length in reciprocal space (Cartesian)
            dk_frac = k_end - k_start
            dk_cart = reclattvec * dk_frac
            segment_length = norm(dk_cart)

            # Interpolate npoints along this segment
            for j in 0:(npoints-1)
                t = j / npoints
                k_frac = k_start + t * dk_frac
                push!(kvecs, k_frac)
                push!(kdist, current_dist + t * segment_length)
            end

            # Record start label
            if i == 1 || isempty(path_labels) || label_positions[end] < current_dist - 1e-10
                push!(path_labels, label_start)
                push!(label_positions, current_dist)
            end

            current_dist += segment_length
        end

        # Add the final point of the segment
        label_final = segment[end]
        k_final = Float64.(points_dict[label_final])
        push!(kvecs, k_final)
        push!(kdist, current_dist)
        push!(path_labels, label_final)
        push!(label_positions, current_dist)
    end

    return kvecs, kdist, path_labels, label_positions
end

"""
    plot_bands(result::InterpolationResult; kwargs...) -> Plot

Plot interpolated band structure along high-symmetry k-path.

# Keyword Arguments
- `npoints::Int = 100`: Number of k-points per path segment
- `emin::Float64 = -1.0`: Minimum energy relative to Fermi [eV]
- `emax::Float64 = 1.0`: Maximum energy relative to Fermi [eV]
- `fermi_line::Bool = true`: Show Fermi level as dashed line
- `output::Union{String,Nothing} = nothing`: Save figure to file if specified
- `kpath::Union{NamedTuple,Nothing} = nothing`: Manual k-path specification (see below)

When `kpath` is provided, the spacegroup metadata is not required. Format:
```julia
kpath = (
    points = Dict("Γ" => [0,0,0], "X" => [0.5,0,0.5], "L" => [0.5,0.5,0.5]),
    paths = [["Γ", "X", "L", "Γ"]]
)
```

# Returns
- Plots.Plot object

# Example
```julia
using BoltzTraP

# Auto k-path from spacegroup (requires spacegroup in metadata)
result = load_interpolation("si_interp.jld2")
plot_bands(result; emin=-5.0, emax=5.0, output="bands.png")

# Manual k-path (for .bt2 files without spacegroup)
result = load_interpolation("si_interp.bt2")
kpath = (
    points = Dict("Γ" => [0,0,0], "X" => [0.5,0,0.5], "L" => [0.5,0.5,0.5]),
    paths = [["Γ", "X", "L", "Γ"]]
)
plot_bands(result; kpath=kpath, emin=-5.0, emax=5.0)
```
"""
function plot_bands(
    result::InterpolationResult;
    npoints::Int = 100,
    emin::Float64 = -1.0,
    emax::Float64 = 1.0,
    fermi_line::Bool = true,
    output::Union{String,Nothing} = nothing,
    kpath::Union{NamedTuple,Nothing} = nothing,
)
    # Extract data from result
    coeffs = result.coeffs
    equivalences = result.equivalences
    lattvec = result.lattvec
    metadata = result.metadata
    atoms = result.atoms

    fermi = get(metadata, "fermi", 0.0)  # in Ha

    # Generate k-path: either from manual specification or auto from spacegroup
    if !isnothing(kpath)
        # Manual k-path specification
        kvecs, kdist, path_labels, label_positions = _generate_manual_kpath(
            kpath, lattvec, npoints
        )
    else
        # Auto k-path from spacegroup
        sgnum = get(metadata, "spacegroup_number", nothing)
        if isnothing(sgnum)
            error(
                "Space group not found in metadata. " *
                "Either provide kpath argument manually or use a .jld2 file with spacegroup metadata."
            )
        end
        sgnum = Int(sgnum)  # Ensure Int64 for Brillouin.jl

        # Get conventional lattice vectors for Brillouin.jl
        # atoms["lattice"] is stored as 3×3 matrix in Bohr
        lattice = atoms["lattice"]
        # Convert Bohr to Angstrom for Brillouin.jl (BOHR_TO_ANG from units.jl)
        lattice_ang = lattice .* BOHR_TO_ANG

        # Convert to tuple format for Brillouin.jl
        Rs = (lattice_ang[:, 1], lattice_ang[:, 2], lattice_ang[:, 3])

        # Generate k-path using Brillouin.jl
        brillouin_kpath = irrfbz_path(sgnum, Rs)
        kpath_interp = Brillouin.interpolate(brillouin_kpath, npoints)

        # Collect k-points (in fractional coordinates)
        kvecs = collect(kpath_interp)

        # Compute cumulative distance along k-path
        kdist = cumdists(kpath_interp)

        # Get path labels and positions by finding high-symmetry points in interpolated path
        all_paths = paths(brillouin_kpath)
        pts = points(brillouin_kpath)

        # Find positions of high-symmetry points by matching k-vectors
        # Track search start index to handle repeated symbols across paths
        path_labels = String[]
        label_positions = Float64[]
        search_start = 1

        for path in all_paths
            for sym in path
                pt = pts[sym]
                # Find the index where this point appears, starting from search_start
                for i in search_start:length(kvecs)
                    kv = kvecs[i]
                    if isapprox(collect(kv), collect(pt), atol=1e-6)
                        push!(path_labels, String(sym))
                        push!(label_positions, kdist[i])
                        search_start = i  # Next search starts from here
                        break
                    end
                end
            end
            # After each path, move search start past the last found point
            search_start += 1
        end
    end

    # Convert kvecs to nk × 3 matrix for getBands
    kpoints = reduce(hcat, kvecs)'

    # Convert equivalences to required format
    equiv_matrices = [Matrix(eq) for eq in equivalences]

    # Interpolate bands at k-path points
    ebands, _ = getBands(kpoints, equiv_matrices, lattvec, coeffs)

    # Convert energies from Ha to eV relative to Fermi (HA_TO_EV from units.jl)
    ebands_eV = (ebands .- fermi) .* HA_TO_EV

    # Remove duplicate labels at same position (merge as "X|K")
    unique_labels = String[]
    unique_positions = Float64[]
    for (lbl, pos) in zip(path_labels, label_positions)
        if isempty(unique_positions) || abs(pos - unique_positions[end]) > 1e-6
            push!(unique_labels, lbl)
            push!(unique_positions, pos)
        else
            # Merge labels at same position (e.g., "X|K")
            unique_labels[end] = unique_labels[end] * "|" * lbl
        end
    end

    # Create plot
    nbands = size(ebands_eV, 1)

    p = Plots.plot(
        size = (600, 400),
        xlabel = "",
        ylabel = "Energy (eV)",
        legend = false,
        grid = true,
        gridalpha = 0.3,
    )

    # Plot each band
    for iband in 1:nbands
        band_energies = ebands_eV[iband, :]
        Plots.plot!(p, kdist, band_energies, color=:blue, linewidth=1.5)
    end

    # Add Fermi level
    if fermi_line
        Plots.hline!(p, [0.0], color=:red, linestyle=:dash, linewidth=1, label="E_F")
    end

    # Set y-axis limits
    Plots.ylims!(p, (emin, emax))

    # Add high-symmetry point labels
    Plots.xticks!(p, (unique_positions, unique_labels))

    # Add vertical lines at high-symmetry points
    for pos in unique_positions
        Plots.vline!(p, [pos], color=:gray, linestyle=:dot, linewidth=0.5, label="")
    end

    # Save if output specified
    if !isnothing(output)
        Plots.savefig(p, output)
        println("Band structure saved to: $output")
    end

    return p
end

"""
    plot_bands(file::String; kwargs...) -> Plot

Plot interpolated band structure from file.

Load interpolation result from file (.jld2 or .bt2) and plot band structure.
For .bt2 files without spacegroup metadata, the `kpath` argument is required.

# Arguments
- `file`: Path to interpolation result file (.jld2 or .bt2)

# Keyword Arguments
Same as `plot_bands(result::InterpolationResult; ...)`.

# Example
```julia
# From .jld2 (auto k-path from spacegroup)
plot_bands("si_interp.jld2"; emin=-5.0, emax=5.0)

# From .bt2 (manual k-path required)
kpath = (
    points = Dict("Γ" => [0,0,0], "X" => [0.5,0,0.5], "L" => [0.5,0.5,0.5]),
    paths = [["Γ", "X", "L", "Γ"]]
)
plot_bands("si_interp.bt2"; kpath=kpath, emin=-5.0, emax=5.0)
```
"""
function plot_bands(
    file::String;
    npoints::Int = 100,
    emin::Float64 = -1.0,
    emax::Float64 = 1.0,
    fermi_line::Bool = true,
    output::Union{String,Nothing} = nothing,
    kpath::Union{NamedTuple,Nothing} = nothing,
)
    result = load_interpolation(file)
    return plot_bands(
        result;
        npoints = npoints,
        emin = emin,
        emax = emax,
        fermi_line = fermi_line,
        output = output,
        kpath = kpath,
    )
end

"""
    plot_transport(result::TransportResult; kwargs...) -> Plot

Plot transport coefficients vs chemical potential or temperature.

# Keyword Arguments
- `quantity::String = "seebeck"`: Property to plot (seebeck, sigma, kappa)
- `component::String = "xx"`: Tensor component (xx, yy, zz, xy, ...)
- `abscissa::String = "mu"`: X-axis variable (mu or T)
- `temperature::Float64 = 300.0`: Temperature for μ plot [K]
- `mu_index::Int = 0`: μ index for T plot (0 = auto-select near Fermi)
- `output::Union{String,Nothing} = nothing`: Save figure to file

# Returns
- Plots.Plot object
"""
function plot_transport(
    result::TransportResult;
    quantity::String = "seebeck",
    component::String = "xx",
    abscissa::String = "mu",
    temperature::Float64 = 300.0,
    mu_index::Int = 0,
    output::Union{String,Nothing} = nothing,
)
    # Parse component
    comp_map = Dict(
        "xx" => (1, 1), "yy" => (2, 2), "zz" => (3, 3),
        "xy" => (1, 2), "xz" => (1, 3), "yz" => (2, 3),
        "yx" => (2, 1), "zx" => (3, 1), "zy" => (3, 2),
    )
    if !haskey(comp_map, lowercase(component))
        error("Unknown component: $component. Use xx, yy, zz, xy, xz, yz.")
    end
    i, j = comp_map[lowercase(component)]

    # Get data array based on quantity
    data = if quantity == "seebeck" || quantity == "S"
        result.seebeck
    elseif quantity == "sigma" || quantity == "conductivity"
        result.sigma
    elseif quantity == "kappa" || quantity == "thermal"
        result.kappa
    else
        error("Unknown quantity: $quantity. Use seebeck, sigma, or kappa.")
    end

    # Get units and labels
    units = Dict(
        "seebeck" => "μV/K",
        "S" => "μV/K",
        "sigma" => "S/m/s",
        "conductivity" => "S/m/s",
        "kappa" => "W/m/K/s",
        "thermal" => "W/m/K/s",
    )

    # Get Fermi energy
    md = result.metadata
    if haskey(md, "fermi_dft_Ha")
        fermi = md["fermi_dft_Ha"]
    elseif haskey(md, "fermi_Ha")
        fermi = md["fermi_Ha"]
    else
        fermi = 0.0
    end
    fermi_eV = fermi * HA_TO_EV  # HA_TO_EV from units.jl

    if abscissa == "mu"
        # Plot vs chemical potential at fixed temperature
        iT = findfirst(t -> isapprox(t, temperature, atol=1.0), result.temperatures)
        if isnothing(iT)
            error("Temperature $temperature K not found. Available: $(result.temperatures)")
        end

        # data shape: (3, 3, nT, nμ)
        y_data = data[i, j, iT, :]

        # Convert units
        if quantity in ["seebeck", "S"]
            y_data = y_data .* 1e6  # V/K → μV/K
        end

        # μ is already in eV, relative to refined Fermi
        # We want relative to DFT Fermi
        mu_eV = result.mu_values .- fermi_eV

        p = Plots.plot(
            mu_eV, y_data,
            xlabel = "μ - E_F (eV)",
            ylabel = "$(uppercase(quantity))_$(component) ($(units[quantity]))",
            title = "T = $(Int(temperature)) K",
            legend = false,
            linewidth = 2,
        )
    else
        # Plot vs temperature at fixed μ
        if mu_index == 0
            # Find μ closest to Fermi
            mu_index = argmin(abs.(result.mu_values .- fermi_eV))
        end

        y_data = data[i, j, :, mu_index]

        if quantity in ["seebeck", "S"]
            y_data = y_data .* 1e6
        end

        p = Plots.plot(
            result.temperatures, y_data,
            xlabel = "Temperature (K)",
            ylabel = "$(uppercase(quantity))_$(component) ($(units[quantity]))",
            title = "μ = $(round(result.mu_values[mu_index] - fermi_eV, digits=3)) eV",
            legend = false,
            linewidth = 2,
        )
    end

    if !isnothing(output)
        Plots.savefig(p, output)
        println("Plot saved to: $output")
    end

    return p
end

"""
    plot_transport(file::String; kwargs...) -> Plot

Plot transport coefficients from file.

Load integration result from file (.jld2) and plot transport coefficients.

# Arguments
- `file`: Path to integration result file (.jld2)

# Keyword Arguments
Same as `plot_transport(result::TransportResult; ...)`.

# Example
```julia
plot_transport("si_transport.jld2"; quantity="seebeck", component="xx")
plot_transport("si_transport.jld2"; quantity="sigma", abscissa="T")
```
"""
function plot_transport(
    file::String;
    quantity::String = "seebeck",
    component::String = "xx",
    abscissa::String = "mu",
    temperature::Float64 = 300.0,
    mu_index::Int = 0,
    output::Union{String,Nothing} = nothing,
)
    result = load_integrate(file)
    return plot_transport(
        result;
        quantity = quantity,
        component = component,
        abscissa = abscissa,
        temperature = temperature,
        mu_index = mu_index,
        output = output,
    )
end
