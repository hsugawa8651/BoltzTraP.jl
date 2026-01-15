# Developer Guide

This page documents internal data structures, types, and workflows for developers extending BoltzTraP.jl.

For validating results against Python BoltzTraP2, see the [Validation](validation.md) page.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                      Public API                              │
│  run_interpolate, run_integrate, load_vasp, load_qe, ...    │
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                    Workflow Layer                            │
│  workflow.jl: orchestrates interpolation and integration    │
└─────────────────────────────────────────────────────────────┘
                              │
┌──────────────┬──────────────┬──────────────┬────────────────┐
│   Symmetry   │ Interpolation│ Reconstruction│   Transport   │
│  sphere.jl   │interpolation │reconstruction │  bandlib.jl   │
│ symmetry.jl  │    .jl       │    .jl        │               │
│equivalences.jl│             │               │               │
└──────────────┴──────────────┴──────────────┴────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                      I/O Layer                               │
│  io/vasp.jl, io/qe.jl, io/serialization.jl, io/loader.jl   │
└─────────────────────────────────────────────────────────────┘
```

## Core Data Types

### InterpolationResult

Stores the result of band structure interpolation.

```julia
struct InterpolationResult
    coeffs::Matrix{ComplexF64}       # Fourier coefficients (nbands × neq)
    equivalences::Vector{Matrix{Int}} # Equivalence class representatives
    lattvec::Matrix{Float64}         # Lattice vectors (3×3, Bohr)
    atoms::Union{Nothing,Dict}       # Atomic structure (optional)
    metadata::Dict{String,Any}       # Additional information
end
```

**Metadata keys:**
- `fermi`: Fermi energy (Ha)
- `nelect`: Number of electrons
- `dosweight`: DOS weight (2.0 for non-spin-polarized)
- `selected_bands`: Band indices used
- `spacegroup_number`, `spacegroup_symbol`: Space group info

### IntegrateResult

Stores transport coefficients from BZ integration.

```julia
struct IntegrateResult
    temperatures::Vector{Float64}    # Temperatures (K)
    mu_values::Vector{Float64}       # Chemical potentials (eV)
    sigma::Array{Float64,4}          # Conductivity (3×3×nT×nμ), S/m
    seebeck::Array{Float64,4}        # Seebeck (3×3×nT×nμ), V/K
    kappa::Array{Float64,4}          # Thermal conductivity (3×3×nT×nμ), W/(m·K)
    dos::Union{Nothing,Dict}         # DOS data (optional)
    metadata::Dict{String,Any}       # Additional information
end
```

### FourierInterpolator

Internal type for band interpolation (not exported).

```julia
struct FourierInterpolator
    coeffs::Matrix{ComplexF64}       # Fourier coefficients
    equivalences::Vector{Vector{SVector{3,Int}}}  # Equivalence classes
    lattvec::SMatrix{3,3,Float64}    # Lattice vectors
end
```

**Key methods:**
- `interpolate_bands(interp, kpoints)` → band energies
- `interpolate_velocities(interp, kpoints)` → band velocities

## I/O Data Format

All loaders (`load_vasp`, `load_qe`, `load_dftk`) return a NamedTuple:

```julia
(
    lattice::Matrix{Float64},     # 3×3, Bohr
    positions::Matrix{Float64},   # 3×natoms, fractional
    species::Vector{String},      # Element symbols
    kpoints::Matrix{Float64},     # 3×nkpts, fractional
    weights::Vector{Float64},     # nkpts
    ebands::Array{Float64,3},     # nbands×nkpts×nspin, Hartree
    occupations::Array{Float64,3},# nbands×nkpts×nspin
    fermi::Float64,               # Hartree
    nelect::Float64,              # Number of electrons
)
```

## Internal Functions

### Sphere Module (`sphere.jl`)

| Function | Description |
|----------|-------------|
| `compute_bounds(lattvec, radius)` | Compute search bounds for lattice points |
| `compute_radius(lattvec, nrot, nkpt_target)` | Estimate radius for target k-points |
| `calc_nrotations(lattvec, positions, types, magmom)` | Count symmetry operations |
| `get_unique_rotations(...)` | Get unique rotation matrices |
| `calc_tensor_basis(...)` | Compute symmetry-adapted tensor basis |

### Equivalences Module (`equivalences.jl`)

| Function | Description |
|----------|-------------|
| `calc_sphere_quotient_set(...)` | Compute equivalence classes |
| `lattice_points_in_sphere(...)` | Find lattice points within radius |

### Interpolation Module (`interpolation.jl`)

| Function | Description |
|----------|-------------|
| `compute_phase_factors(kpoints, equivalences)` | Compute exp(2πi k·R) |
| `FourierInterpolator(kpts, ebands, equivs, lattvec)` | Fit Fourier coefficients |
| `interpolate_bands(interp, kpoints)` | Evaluate bands at k-points |
| `interpolate_velocities(interp, kpoints)` | Evaluate velocities |

### Transport Module (`bandlib.jl`)

| Function | Description |
|----------|-------------|
| `fermi_dirac(ε, μ, kT)` | Fermi-Dirac distribution |
| `dfermi_dirac_de(ε, μ, kT)` | Derivative of f(E) |
| `calc_Onsager_coefficients(...)` | Compute L0, L1, L2 tensors |
| `calc_transport_coefficients(...)` | Compute σ, S, κ from Onsager |
| `calc_N(...)` | Carrier concentration |
| `solve_for_mu(...)` | Find μ for target carrier density |

## Unit Constants (`units.jl`)

```julia
const BOHR_TO_ANG = 0.529177210903     # Length conversion
const ANG_TO_BOHR = 1 / BOHR_TO_ANG
const HA_TO_EV = 27.211386245988       # Energy conversion
const EV_TO_HA = 1 / HA_TO_EV
const KB_AU = 3.1668115634438576e-6    # Boltzmann constant (Ha/K)
const SIGMA_CONV = 5.2586177987809155e-24   # σ/τ: S/(m·s)
const SEEBECK_CONV = 1.9325063968640632e-25 # S intermediate
const KAPPA_CONV = 7.101830018500863e-27    # κ/τ: W/(m·K·s)
```

## Adding New Features

### Adding a New I/O Format

1. Create `src/io/newformat.jl`:

```julia
function load_newformat(directory::String)
    # Parse files and return NamedTuple
    return (
        lattice = ...,   # 3×3, Bohr
        kpoints = ...,   # 3×nkpts, fractional
        ebands = ...,    # nbands×nkpts×nspin, Hartree
        fermi = ...,     # Hartree
        # ... other fields
    )
end
```

2. Register in `src/io/loader.jl`:

```julia
function _detect_newformat(directory::String)
    isfile(joinpath(directory, "newformat.out"))
end

# In _init_loaders():
register_loader("NewFormat", _detect_newformat, load_newformat)
```

3. Export in `src/BoltzTraP.jl`:

```julia
export load_newformat
```

### Adding a New Transport Property

1. Add calculation in `src/bandlib.jl`
2. Add to `IntegrateResult` struct in `src/io/serialization.jl`
3. Update `run_integrate` in `src/workflow.jl`
4. Add CLI option in `src/cli.jl` (if needed)

## Debugging

### CLI Debug Mode

Use the `--debug` flag to enable detailed logging output:

```bash
# Interpolation with debug output
boltztrap interpolate ./Si.vasp --debug

# Integration with debug output
boltztrap integrate si_interp.jld2 -t 300 --debug
```

Debug output includes:
- CLI argument values received
- DFT format detection process
- File loading details
- Workflow function parameters

### Environment Variable

Alternatively, set the `JULIA_DEBUG` environment variable:

```bash
JULIA_DEBUG=BoltzTraP boltztrap interpolate ./Si.vasp
```

Or in Julia:

```julia
ENV["JULIA_DEBUG"] = "BoltzTraP"
using BoltzTraP
# Now all @debug statements will be printed
```

### Logging Levels

BoltzTraP.jl uses Julia's standard `Logging` module:

| Level | Purpose | Enabled by |
|-------|---------|------------|
| `@debug` | Internal details, argument values | `--debug` flag |
| `@info` | Success messages | Always |
| `@warn` | Non-fatal issues | Always |
| `@error` | Fatal errors | Always |

### Verbose Mode vs Debug Mode

| Flag | Output | Audience |
|------|--------|----------|
| `-v, --verbose` | Progress messages via `println` | End users |
| `--debug` | Internal state via `@debug` | Developers |

Example with both:
```bash
boltztrap interpolate ./Si.vasp -v --debug
```

## Reference Tests (reftest)

BoltzTraP.jl includes 75 reference tests that validate numerical equivalence with Python BoltzTraP2. These tests compare Julia output against pre-computed Python reference data.

For complete documentation on running and extending reference tests, see the dedicated **[Reference Tests](reftest.md)** page.

Quick start:

```bash
# Generate reference data (requires Python BoltzTraP2)
cd reftest
pip install boltztrap2
python generate_1_sphere.py
python generate_2_interpolation_si.py
# ... (see reftest.md for full list)

# Run tests
cd ..
julia --project -e 'using Pkg; Pkg.test()'
```

## Functional Tests (ftest)

Functional tests are located in the `ftest/` directory. Unlike unit tests in `test/`, these are standalone scripts for manual verification and debugging.

### Available Tests

| File | Purpose |
|------|---------|
| `test_cli_all_options.jl` | Comprehensive test of all CLI options |
| `test_cli_args.jl` | Basic argument propagation test |
| `test_format_detection.jl` | DFT format auto-detection test |

### Running Functional Tests

```bash
# Run from project root
cd BoltzTraP.jl

# Test all CLI options (requires VASP test data)
julia --project=. ftest/test_cli_all_options.jl

# With custom data directory
julia --project=. ftest/test_cli_all_options.jl /path/to/vasp/dir

# Test format detection
julia --project=. ftest/test_format_detection.jl

# With debug output
JULIA_DEBUG=BoltzTraP julia --project=. ftest/test_cli_args.jl
```

### Test Data

Functional tests expect VASP data in `../BoltzTraP2-public/data/Si.vasp/` by default. You can specify a different directory as a command-line argument.

Required files:
```
vasp_directory/
└── vasprun.xml
```

### Writing New Functional Tests

Template for a new functional test:

```julia
# ftest/test_new_feature.jl
using Logging

# Enable debug logging
ENV["JULIA_DEBUG"] = "BoltzTraP"
global_logger(ConsoleLogger(stderr, Logging.Debug))

using BoltzTraP

function main()
    println("=" ^ 60)
    println("Functional Test: New Feature")
    println("=" ^ 60)

    # Test 1
    println("-" ^ 60)
    println("Test 1: Description")
    println("-" ^ 60)
    @info "Running test 1"

    # ... test code ...

    @assert condition "Error message"
    println("  ✓ PASS")

    println("=" ^ 60)
    println("All tests passed!")
    println("=" ^ 60)
end

main()
```

## DFT Format Auto-Detection

### How It Works

The CLI uses `load_dft()` to automatically detect the DFT format:

1. Try each registered loader in reverse order (VASP first, then QE)
2. Each loader has a detector function that checks for format-specific files
3. First successful detection wins

### Detector Functions

| Format | Detection |
|--------|-----------|
| VASP | `vasprun.xml` exists in directory |
| QE | `*.save/data-file-schema.xml` or `*.save/data-file.xml` exists |

### Forcing a Format

Use `--format` to skip auto-detection:

```bash
boltztrap interpolate ./data -f vasp  # Force VASP
boltztrap interpolate ./data -f qe    # Force QE
boltztrap interpolate ./data -f auto  # Auto-detect (default)
```

### Adding a New Format

1. Create loader in `src/io/newformat.jl`
2. Register in `src/io/loader.jl`:

```julia
function _detect_newformat(directory::String)
    # Return true if format detected
    isfile(joinpath(directory, "newformat.out"))
end

# In _init_loaders():
register_loader("NewFormat", _detect_newformat, load_newformat)
```

## Code Structure

```
src/
├── cli.jl              # CLI commands (Comonicon)
├── workflow.jl         # High-level workflows with @debug logging
├── io/
│   ├── loader.jl       # Auto-detection logic
│   ├── vasp.jl         # VASP loader
│   └── qe.jl           # QE loader
└── ...

ftest/                  # Functional tests (manual)
test/                   # Unit tests (automated)
```

## Troubleshooting

### "No compatible DFT format found"

1. Check directory path is correct
2. Run with `--debug` to see detection attempts
3. Try explicit format: `--format vasp`
4. Verify required files exist (e.g., `vasprun.xml`)

### Debug Output Not Showing

1. Ensure `--debug` is passed to the command
2. Or set `ENV["JULIA_DEBUG"] = "BoltzTraP"` before `using BoltzTraP`
3. Check stderr (debug goes to stderr, not stdout)

### Functional Test Failures

1. Check test data path exists
2. Run with debug logging enabled
3. Check assertion error messages for details

## Migration Guide

### v0.2 (upcoming): Magnetic Material Support

The `spintype` metadata field has been added to both `InterpolationResult` and `IntegrateResult` to prepare for magnetic material support in v0.2.

**Current behavior (v0.1):**
- All calculations assume non-magnetic (unpolarized) materials
- `metadata["spintype"]` = `"Unpolarized"` (always)

**Planned v0.2 changes:**
- `SpinType` enum: `Unpolarized`, `Collinear`, `NonCollinear`
- Parametric types: `InterpolationResult{ST}` where `ST<:SpinType`
- Spin-polarized VASP/QE support

**Forward compatibility:**
Files saved with v0.1 will be loadable in v0.2. The `spintype` metadata field ensures smooth migration.
