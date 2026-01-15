# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# Generate Julia transport results using Python's mu range for fair comparison.
#
# Usage:
#   julia --project validation/generate_transport.jl <material> [options]
#
# Examples:
#   julia --project validation/generate_transport.jl si
#   julia --project validation/generate_transport.jl pbte
#   julia --project validation/generate_transport.jl all
#
# Output:
#   validation/<material>_transport.jld2

using NPZ
using BoltzTraP

# Material configurations
const MATERIALS = Dict(
    "si" => (
        input = "reftest/data/vasp_si.npz",
        reference = "reftest/data/si_end2end.npz",
    ),
    "pbte" => (
        input = "reftest/data/vasp_pbte.npz",
        reference = "reftest/data/pbte_end2end.npz",
    ),
)

# Constants
const HA_TO_EV = 27.211386245988

function decode_species(symbols_raw)
    """Decode species from npz format (ASCII byte array)."""
    # symbols are stored as ASCII codes in npz, comma-separated
    symbols_str = String(UInt8.(symbols_raw))
    return split(symbols_str, ",")
end

function load_dft_data(input_file)
    """Load DFT data from VASP npz file."""
    data = npzread(input_file)
    species = decode_species(data["symbols"])
    return (
        lattice = data["lattvec"],
        positions = transpose(data["positions"]),
        species = species,
        kpoints = transpose(data["kpoints"]),
        ebands = data["ebands"],
        fermi = data["fermi"],
        nelect = data["nelect"],
        dosweight = data["dosweight"],
    )
end

function load_python_params(reference_file)
    """Load Python's mur and Tr from reference npz file."""
    ref = npzread(reference_file)
    return (
        mur = ref["mur"],      # Chemical potential in Hartree
        Tr = ref["Tr"],        # Temperatures in K
        fermi = ref["fermi"],  # Fermi level in Hartree
    )
end

function generate_transport(material; kpoints=5000, verbose=true)
    """Generate Julia transport results using Python's mu range."""
    if !haskey(MATERIALS, material)
        error("Unknown material: $material. Available: $(keys(MATERIALS))")
    end

    config = MATERIALS[material]

    println("=" ^ 60)
    println("Generating transport for: $material")
    println("=" ^ 60)

    # Check files exist
    if !isfile(config.input)
        error("Input file not found: $(config.input)")
    end
    if !isfile(config.reference)
        error("Reference file not found: $(config.reference)")
    end

    # Load data
    println("\nLoading DFT data from $(config.input)...")
    dft_data = load_dft_data(config.input)
    println("  species: $(dft_data.species)")
    println("  k-points: $(size(dft_data.kpoints, 2))")
    println("  bands: $(size(dft_data.ebands, 1))")

    println("\nLoading Python reference parameters from $(config.reference)...")
    py_params = load_python_params(config.reference)
    println("  mur points: $(length(py_params.mur))")
    println("  mur range (Ha): $(minimum(py_params.mur)) to $(maximum(py_params.mur))")
    println("  mur - fermi range (Ha): $(minimum(py_params.mur) - py_params.fermi) to $(maximum(py_params.mur) - py_params.fermi)")
    println("  mur range (eV): $(minimum(py_params.mur) * HA_TO_EV) to $(maximum(py_params.mur) * HA_TO_EV)")
    println("  temperatures: $(py_params.Tr) K")

    # mur is already in Hartree (run_integrate expects Ha)
    mur_Ha = py_params.mur

    # Run interpolation
    println("\nRunning interpolation (kpoints=$kpoints)...")
    interp = run_interpolate(dft_data; source="VASP", kpoints=kpoints, verbose=verbose)
    println("  Done. Equivalences: $(length(interp.equivalences))")

    # Run integration with Python's mur (in Ha)
    println("\nRunning integration with Python's mur...")
    transport = run_integrate(interp, mur_Ha;
        temperatures=collect(Float64, py_params.Tr),
        verbose=verbose
    )
    println("  Done.")

    # Save result
    output_file = joinpath(@__DIR__, "$(material)_transport.jld2")
    println("\nSaving to $output_file...")
    save_integrate(output_file, transport)
    println("  Done.")

    # Print summary
    println("\n" * "=" ^ 60)
    println("Summary for $material:")
    println("  mu points: $(length(transport.mu_values))")
    println("  mu range (eV): $(minimum(transport.mu_values)) to $(maximum(transport.mu_values))")
    println("  temperatures: $(transport.temperatures) K")
    println("  output: $output_file")
    println("=" ^ 60)

    return transport
end

function print_usage()
    println("""
    Generate Julia transport results using Python's mu range for fair comparison.

    Usage:
      julia --project validation/generate_transport.jl <material> [options]

    Arguments:
      <material>   Material name: si, pbte, or "all"

    Options:
      -k, --kpoints <N>   Number of k-points for interpolation (default: 5000)
      -q, --quiet         Suppress verbose output
      -h, --help          Show this help

    Examples:
      julia --project validation/generate_transport.jl si
      julia --project validation/generate_transport.jl pbte -k 10000
      julia --project validation/generate_transport.jl all

    Output:
      validation/<material>_transport.jld2

    Note:
      This script reads Python's mur (chemical potential range) from the
      reference npz files to ensure fair comparison with Python BoltzTraP2.
    """)
end

function parse_args(args)
    if isempty(args) || "-h" in args || "--help" in args
        return nothing
    end

    material = args[1]
    kpoints = 5000
    verbose = true

    i = 2
    while i <= length(args)
        if args[i] == "-k" || args[i] == "--kpoints"
            i += 1
            kpoints = parse(Int, args[i])
        elseif args[i] == "-q" || args[i] == "--quiet"
            verbose = false
        end
        i += 1
    end

    return (material = lowercase(material), kpoints = kpoints, verbose = verbose)
end

function main()
    args = parse_args(ARGS)

    if isnothing(args)
        print_usage()
        return
    end

    if args.material == "all"
        for mat in sort(collect(keys(MATERIALS)))
            generate_transport(mat; kpoints=args.kpoints, verbose=args.verbose)
            println()
        end
    else
        generate_transport(args.material; kpoints=args.kpoints, verbose=args.verbose)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
