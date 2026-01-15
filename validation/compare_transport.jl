# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Compare transport coefficients between Python BoltzTraP2 and Julia BoltzTraP.jl
#
# Usage:
#   julia --project validation/compare_transport.jl <python_npz> <julia_jld2> [options]
#
# Examples:
#   julia --project validation/compare_transport.jl reftest/data/si_end2end.npz si_transport.jld2 --title Si
#   julia --project validation/compare_transport.jl reftest/data/si_end2end.npz si_transport.jld2 --title Si --alltemp
#   julia --project validation/compare_transport.jl reftest/data/pbte_end2end.npz pbte_transport.jld2 --title PbTe
#
# Output:
#   validation/transport_<title>_<T>K.png (default)
#   or specified by -o option

using NPZ
using JLD2
using Plots
using Printf
using LinearAlgebra
using LaTeXStrings
using Interpolations

# Use GR backend for high-quality output
gr()

# Constants
const Ha_to_eV = 27.211386245988
const V_to_uV = 1e6

# Property configurations
const PROPERTIES = Dict(
    "S" => (
        name = "Seebeck coefficient",
        symbol = L"S_{xx}",
        unit = L"$\mu$V/K",
        extract_py = (data, iT) -> data.S[iT, :, 1, 1] .* V_to_uV,
        extract_jl = (result, iT) -> vec(result["seebeck"][1, 1, iT, :]) .* V_to_uV,
    ),
    "sigma" => (
        name = "Electrical conductivity",
        symbol = L"\sigma_{xx}/\tau",
        unit = L"$(\Omega \cdot m \cdot s)^{-1}$",
        extract_py = (data, iT) -> data.sigma[iT, :, 1, 1],
        extract_jl = (result, iT) -> vec(result["sigma"][1, 1, iT, :]),
    ),
    "kappa" => (
        name = "Electronic thermal conductivity",
        symbol = L"\kappa_{xx}/\tau",
        unit = L"W/(m $\cdot$ K $\cdot$ s)",
        extract_py = (data, iT) -> data.kappa[iT, :, 1, 1],
        extract_jl = (result, iT) -> vec(result["kappa"][1, 1, iT, :]),
    ),
)

function load_python_reference(datafile)
    """Load Python BoltzTraP2 end-to-end transport results (npz format)."""
    data = npzread(datafile)
    return (
        Tr = data["Tr"],           # Temperatures (K)
        mur = data["mur"],         # Chemical potential (Hartree)
        fermi = data["fermi"],     # Fermi level (Hartree)
        S = data["S"],             # Seebeck (nT, nmu, 3, 3) [V/K]
        sigma = data["sigma"],     # Conductivity (nT, nmu, 3, 3)
        kappa = data["kappa"],     # Thermal conductivity (nT, nmu, 3, 3)
    )
end

function load_julia_result(datafile)
    """Load Julia BoltzTraP.jl transport results (jld2 format)."""
    return jldopen(datafile, "r") do file
        Dict(
            "temperatures" => file["temperatures"],
            "mu_values" => file["mu_values"],
            "metadata" => file["metadata"],
            "seebeck" => file["seebeck"],
            "sigma" => file["sigma"],
            "kappa" => file["kappa"],
        )
    end
end

function compare_property(py_data, jl_result, prop_name; temperature=300.0)
    """Compare a specific transport property between Python and Julia."""
    prop = PROPERTIES[prop_name]

    # Get mu values relative to Fermi level (in eV)
    jl_mu_eV = jl_result["mu_values"]
    jl_fermi_eV = jl_result["metadata"]["fermi_dft_eV"]
    jl_mu_rel = jl_mu_eV .- jl_fermi_eV

    py_mu_rel = (py_data.mur .- py_data.fermi) .* Ha_to_eV

    # Find temperature index
    jl_iT = findfirst(t -> isapprox(t, temperature, atol=1.0), jl_result["temperatures"])
    py_iT = findfirst(t -> isapprox(t, temperature, atol=1.0), py_data.Tr)

    if isnothing(jl_iT) || isnothing(py_iT)
        error("Temperature $(temperature) K not found in data")
    end

    # Extract property values
    jl_vals = prop.extract_jl(jl_result, jl_iT)
    py_vals = prop.extract_py(py_data, py_iT)

    # Check if using same mu grid
    same_mu_grid = length(jl_mu_rel) == length(py_mu_rel) &&
                   all(isapprox.(jl_mu_rel, py_mu_rel, atol=1e-6))

    return (
        jl_mu = jl_mu_rel,
        py_mu = py_mu_rel,
        jl_vals = jl_vals,
        py_vals = py_vals,
        T = jl_result["temperatures"][jl_iT],
        same_mu_grid = same_mu_grid,
        prop_name = prop_name,
    )
end

function create_joss_figure(py_data, jl_result, output_path; title="", xlims=(-0.5, 0.5), temperature=300.0)
    """Create JOSS-style 3-panel vertical layout (S, sigma, kappa)."""

    # Compare all properties
    cmp_S = compare_property(py_data, jl_result, "S"; temperature=temperature)
    cmp_sigma = compare_property(py_data, jl_result, "sigma"; temperature=temperature)
    cmp_kappa = compare_property(py_data, jl_result, "kappa"; temperature=temperature)

    # Create 3-panel vertical figure
    println("Generating 3-panel figure...")

    # Panel (a): Seebeck
    p1 = scatter(
        cmp_S.py_mu, cmp_S.py_vals,
        label="Python BoltzTraP2",
        markersize=5, markershape=:circle, markercolor=:red,
        markerstrokewidth=0, alpha=0.8,
        ylabel=L"$S_{xx}$ ($\mu$V/K)",
        xlims=xlims, xformatter=_->"",
        legend=:bottomleft, legendfontsize=8,
        title="$title @ $(@sprintf("%.0f", cmp_S.T)) K",
        grid=true, gridalpha=0.3,
        left_margin=5Plots.mm, right_margin=5Plots.mm,
        tickfontsize=9, guidefontsize=10, titlefontsize=11
    )
    scatter!(p1, cmp_S.jl_mu, cmp_S.jl_vals,
        label="Julia BoltzTraP.jl",
        markersize=6, markershape=:xcross, markercolor=:blue,
        markerstrokewidth=2, alpha=0.9
    )

    # Panel (b): Conductivity (log scale)
    py_sigma_pos = max.(abs.(cmp_sigma.py_vals), 1e10)
    jl_sigma_pos = max.(abs.(cmp_sigma.jl_vals), 1e10)
    p2 = scatter(
        cmp_sigma.py_mu, py_sigma_pos,
        label="", markersize=5, markershape=:circle, markercolor=:red,
        markerstrokewidth=0, alpha=0.8,
        ylabel=L"$\sigma_{xx}/\tau$ (1/$\Omega$ms)",
        xlims=xlims, xformatter=_->"",
        yscale=:log10,
        legend=false,
        grid=true, gridalpha=0.3,
        left_margin=5Plots.mm, right_margin=5Plots.mm,
        tickfontsize=9, guidefontsize=10
    )
    scatter!(p2, cmp_sigma.jl_mu, jl_sigma_pos,
        label="", markersize=6, markershape=:xcross, markercolor=:blue,
        markerstrokewidth=2, alpha=0.9
    )

    # Panel (c): Thermal conductivity (log scale)
    py_kappa_pos = max.(abs.(cmp_kappa.py_vals), 1e5)
    jl_kappa_pos = max.(abs.(cmp_kappa.jl_vals), 1e5)
    p3 = scatter(
        cmp_kappa.py_mu, py_kappa_pos,
        label="", markersize=5, markershape=:circle, markercolor=:red,
        markerstrokewidth=0, alpha=0.8,
        xlabel=L"$\mu - E_F$ (eV)",
        ylabel=L"$\kappa_{xx}/\tau$ (W/mKs)",
        xlims=xlims,
        yscale=:log10,
        legend=false,
        grid=true, gridalpha=0.3,
        left_margin=5Plots.mm, right_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        tickfontsize=9, guidefontsize=10
    )
    scatter!(p3, cmp_kappa.jl_mu, jl_kappa_pos,
        label="", markersize=6, markershape=:xcross, markercolor=:blue,
        markerstrokewidth=2, alpha=0.9
    )

    # Combine panels
    fig = plot(p1, p2, p3,
        layout=(3, 1),
        size=(500, 600),
        dpi=300,
        fontfamily="Computer Modern"
    )

    savefig(fig, output_path)
    println("Figure saved: $output_path")

    # Print statistics (only if same grid)
    println("\nStatistics:")
    for (name, cmp) in [("S", cmp_S), ("sigma", cmp_sigma), ("kappa", cmp_kappa)]
        if cmp.same_mu_grid
            diff = abs.(cmp.jl_vals .- cmp.py_vals)
            println("  $name: max diff = $(@sprintf("%.4g", maximum(diff))), mean = $(@sprintf("%.4g", sum(diff)/length(diff)))")
        else
            println("  $name: different mu grids (Julia: $(length(cmp.jl_vals)), Python: $(length(cmp.py_vals)))")
        end
    end

    return fig
end

function print_usage()
    println("""
    Compare transport coefficients between Python BoltzTraP2 and Julia BoltzTraP.jl

    Usage:
      julia --project validation/compare_transport.jl <python_npz> <julia_jld2> [options]

    Arguments:
      <python_npz>   Path to Python BoltzTraP2 results (npz format from reftest)
      <julia_jld2>   Path to Julia BoltzTraP.jl transport results (jld2 format)

    Options:
      -o, --output <file>   Output figure path (default: validation/transport_<title>_<T>K.png)
      --title <name>        Material name for figure title (default: untitled)
      --xlims <min,max>     X-axis limits in eV (default: -0.5,0.5)
      -t, --temperature <T> Temperature in K (default: 300)
      --alltemp             Generate figures for all common temperatures
      -h, --help            Show this help

    Examples:
      julia --project validation/compare_transport.jl reftest/data/si_end2end.npz si_transport.jld2 --title Si
      julia --project validation/compare_transport.jl reftest/data/si_end2end.npz si_transport.jld2 --title Si --alltemp
      julia --project validation/compare_transport.jl reftest/data/pbte_end2end.npz pbte_transport.jld2 --title PbTe --xlims -0.3,0.3
    """)
end

function parse_args(args)
    """Parse command line arguments."""
    if isempty(args) || "-h" in args || "--help" in args
        return nothing
    end

    # Required arguments
    if length(args) < 2
        println("Error: Missing required arguments")
        return nothing
    end

    python_npz = args[1]
    julia_jld2 = args[2]

    # Default options
    output = ""  # Will be set later based on title
    title = "untitled"
    xlims = (-0.5, 0.5)
    temperature = 300.0
    all_temperatures = false

    # Parse optional arguments
    i = 3
    while i <= length(args)
        if args[i] == "-o" || args[i] == "--output"
            i += 1
            output = args[i]
        elseif args[i] == "--title"
            i += 1
            title = args[i]
        elseif args[i] == "--xlims"
            i += 1
            parts = split(args[i], ",")
            xlims = (parse(Float64, parts[1]), parse(Float64, parts[2]))
        elseif args[i] == "-t" || args[i] == "--temperature"
            i += 1
            temperature = parse(Float64, args[i])
        elseif args[i] == "--alltemp"
            all_temperatures = true
        end
        i += 1
    end

    return (
        python_npz = python_npz,
        julia_jld2 = julia_jld2,
        output = output,
        title = title,
        xlims = xlims,
        temperature = temperature,
        all_temperatures = all_temperatures,
    )
end

function get_output_path(title, temperature, base_output="")
    """Generate output path for a given temperature."""
    safe_title = replace(title, " " => "_")
    temp_int = Int(round(temperature))
    if isempty(base_output)
        return joinpath(@__DIR__, "transport_$(safe_title)_$(temp_int)K.png")
    else
        # If output specified, append temperature before extension
        base, ext = splitext(base_output)
        return "$(base)_$(temp_int)K$(ext)"
    end
end

function main()
    args = parse_args(ARGS)

    if isnothing(args)
        print_usage()
        return
    end

    println("=" ^ 60)
    println("Comparing transport coefficients")
    println("=" ^ 60)
    println()
    println("Python data: $(args.python_npz)")
    println("Julia data:  $(args.julia_jld2)")
    println("Title:       $(args.title)")
    println("X limits:    $(args.xlims)")
    if args.all_temperatures
        println("Temperature: all (--alltemp)")
    else
        println("Temperature: $(args.temperature) K")
    end
    println()

    # Check files exist
    if !isfile(args.python_npz)
        error("Python npz file not found: $(args.python_npz)")
    end
    if !isfile(args.julia_jld2)
        error("Julia jld2 file not found: $(args.julia_jld2)")
    end

    # Load data
    println("Loading Python BoltzTraP2 data...")
    py_data = load_python_reference(args.python_npz)
    println("  Temperatures: $(py_data.Tr) K")
    println("  mu points: $(length(py_data.mur))")

    println("\nLoading Julia BoltzTraP.jl data...")
    jl_result = load_julia_result(args.julia_jld2)
    println("  Temperatures: $(jl_result["temperatures"]) K")
    println("  mu points: $(length(jl_result["mu_values"]))")

    # Determine temperatures to plot
    if args.all_temperatures
        # Find common temperatures between Python and Julia data
        py_temps = Set(py_data.Tr)
        jl_temps = Set(jl_result["temperatures"])
        common_temps = sort(collect(intersect(py_temps, jl_temps)))
        if isempty(common_temps)
            error("No common temperatures found between Python and Julia data")
        end
        println("\nCommon temperatures: $(common_temps) K")
    else
        common_temps = [args.temperature]
    end

    # Create figures for each temperature
    for temp in common_temps
        output_path = get_output_path(args.title, temp, args.output)
        println("\nGenerating figure for T=$(Int(round(temp)))K...")
        println("  Output: $output_path")
        create_joss_figure(py_data, jl_result, output_path;
            title=args.title,
            xlims=args.xlims,
            temperature=temp
        )
    end

    println("\n" * "=" ^ 60)
    println("Done ($(length(common_temps)) figure(s) generated)")
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
