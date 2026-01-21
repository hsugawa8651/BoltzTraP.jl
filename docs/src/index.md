# BoltzTraP.jl

Julia implementation of [BoltzTraP2](https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap2/), a band structure interpolation and semi-classical transport coefficient calculator.

## Features

!!! info "BoltzTraP.jl provides two workflows"
    - **[CLI workflow](@ref "CLI Workflow")** (Command Line): Run from terminal with `boltztrap` command. No programming required.
    - **[API workflow](@ref "API Workflow")** (Julia Session): Call functions from Julia REPL or scripts. Suitable for customization and automation.

    Both workflows produce identical results.

### Available Commands

| Feature | CLI | API | Input | Output |
|---------|-----|-----|-------|--------|
| Band interpolation | `boltztrap interpolate` | [`run_interpolate`](@ref)`()` | VASP, QE, DFTK (`scfres`) | [`InterpolationResult`](@ref) (`.jld2`, `.bt2`) |
| Transport calculation | `boltztrap integrate` | [`run_integrate`](@ref)`()` | [`InterpolationResult`](@ref) (`.jld2`, `.bt2`) | [`TransportResult`](@ref) (`.jld2`) |
| Plot bands | `boltztrap plotbands` | `plot_bands()` | [`InterpolationResult`](@ref) (`.jld2`, `.bt2`) | PNG, PDF |
| Plot transport | `boltztrap plot` | `plot_transport()` | [`TransportResult`](@ref) (`.jld2`) | PNG, PDF |
| Describe results | `boltztrap describe` | [`describe`](@ref)`()` | [`InterpolationResult`](@ref) or [`TransportResult`](@ref) | - |

See [CLI Workflow](@ref) and [API Workflow](@ref) for detailed usage. For file format details, see [Input Formats](@ref input-formats) and [Output Formats](@ref output-formats).

## Quick Start

```bash
# 1. Interpolate band structure from VASP data
boltztrap interpolate ./Si.vasp -v
# Output: Si.vasp_interp.jld2

# 2. Compute transport coefficients at 300K
boltztrap integrate Si.vasp_interp.jld2 -t 300 -v
# Output: Si.vasp_interp_transport.jld2

# 3. Plot Seebeck coefficient
boltztrap plot Si.vasp_interp_transport.jld2 -q seebeck -o seebeck.png
```

For detailed CLI usage, see [CLI Workflow](@ref). For Julia API, see [API Workflow](@ref).

## Scope and Limitations

This package currently supports **non-magnetic materials** (spin-polarized support planned).

### Features beyond Python BoltzTraP2

- **Julia API**: Programmatic workflow with `run_interpolate()` and `run_integrate()` functions
- **[DFTK.jl](https://dftk.org/) integration**: Direct coupling with Julia-native DFT package via `load_dftk()`
- **JLD2 format**: Fast HDF5-based native Julia format for result files

## Installation

**Requirements:** Julia 1.10 or later

```julia
using Pkg
Pkg.add("BoltzTraP")
```

### [CLI Setup (Optional)](@id cli-setup)

For command-line usage (`boltztrap` command), additional setup is required:

**Option 1: Shell Alias (Recommended)**

Add to `~/.bashrc` or `~/.zshrc`:

```bash
alias boltztrap='julia -e "using BoltzTraP; BoltzTraP.command_main(ARGS)" --'
```

**Option 2: Install Binary**

```julia
using BoltzTraP
BoltzTraP.comonicon_install()  # Installs to ~/.julia/bin/
```

Then add to PATH: `export PATH="$HOME/.julia/bin:$PATH"`

!!! note
    For API usage ([`run_interpolate`](@ref), [`run_integrate`](@ref)), no additional setup is needed.

## Compatibility with Python BoltzTraP2

BoltzTraP.jl aims for numerical compatibility with Python BoltzTraP2. To ensure bit-for-bit compatibility with the reference implementation, BoltzTraP.jl uses explicit numerical constants for all unit conversions, avoiding potential discrepancies from unit libraries.

Key differences:

- Julia uses atomic units (Hartree) internally, consistent with DFTK.jl
- Default output format is JLD2 (Julia native); `.bt2` format also supported for Python compatibility

## Citation

If you use BoltzTraP.jl in your research, please cite both:

**BoltzTraP.jl:**
> Sugawara, H. (2026). BoltzTraP.jl: Julia implementation of BoltzTraP2. Zenodo. [doi:10.5281/zenodo.18253186](https://doi.org/10.5281/zenodo.18253186)

**Original BoltzTraP2:**
> Madsen, G. K., Carrete, J., & Verstraete, M. J. (2018). BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients. *Computer Physics Communications*, 231, 140-145. [doi:10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

## License

GPL-3.0-or-later (same as original BoltzTraP2)

## Credits

Julia port by Hiroharu Sugawara.

Based on [BoltzTraP2 v25.11.1](https://pypi.org/project/BoltzTraP2/25.11.1/) ([source](https://gitlab.com/sousaw/BoltzTraP2)) by:
- Georg K. H. Madsen
- Jesús Carrete
- Matthieu J. Verstraete
