# CLI Workflow

This page describes the **CLI workflow** - using BoltzTraP.jl from the command line.

!!! info "CLI vs API Workflow"
    - **CLI workflow**: Run calculations via `boltztrap` shell command
    - **API workflow**: Call Julia functions ([`run_interpolate`](@ref), [`run_integrate`](@ref)) directly

    For API workflow, see [API Workflow](@ref).

## Prerequisites

- BoltzTraP.jl installed (see [Installation](@ref))
- CLI set up (see [CLI Setup](@ref cli-setup)) — requires additional step after `using BoltzTraP`
- DFT calculation results ([VASP, Quantum ESPRESSO, or other supported formats](@ref "Input Formats"))

---

## Workflow Overview

```
DFT Output (vasprun.xml)
         │
         ▼
┌─────────────────────────────────────┐
│  boltztrap interpolate ./Si.vasp   │
└─────────────────────────────────────┘
         │
         ▼
   interpolation.jld2
         │
         ▼
┌─────────────────────────────────────┐
│  boltztrap integrate interp.jld2   │
│               -t 300               │
└─────────────────────────────────────┘
         │
         ▼
   transport.jld2 (σ, S, κ)
```

---

## Commands

### `interpolate`

```@docs
BoltzTraP.interpolate
```

!!! tip "API Workflow"
    For programmatic access, see [`run_interpolate`](@ref) in the [API Workflow](@ref) page.

---

### `integrate`

```@docs
BoltzTraP.integrate
```

!!! tip "API Workflow"
    For programmatic access, see [`run_integrate`](@ref) in the [API Workflow](@ref) page.

---

### `describe`

Inspect interpolation or transport result files. See [Output Formats](@ref) for details.

```bash
boltztrap describe interpolation.jld2
boltztrap describe transport.jld2
```

---

### `plot`

```@docs
BoltzTraP.plot
```

---

### `plotbands`

```@docs
BoltzTraP.plotbands
```

---

## Complete Example

```bash
# 1. Prepare VASP calculation directory
ls ./Si.vasp/
# vasprun.xml  POSCAR

# 2. Run interpolation
boltztrap interpolate ./Si.vasp -k 5000 -v
# Output: Si.vasp_interp.jld2

# 3. Compute transport at multiple temperatures
boltztrap integrate Si.vasp_interp.jld2 -t "100:800:100" -v
# Output: Si.vasp_interp_transport.jld2

# 4. (Optional) Load results in Julia for analysis
julia -e '
using JLD2, BoltzTraP
transport = load_integrate("Si.vasp_interp_transport.jld2")
println("Temperatures: ", transport.temperatures)
println("Seebeck at 300K: ", transport.seebeck[1,1,3,:])
'
```

---

## Help

```bash
# General help
boltztrap --help

# Command-specific help
boltztrap interpolate --help
boltztrap integrate --help
boltztrap describe --help
boltztrap plot --help
boltztrap plotbands --help
```

---

## Next Steps

- [Input Formats](@ref) - DFT input formats and loaders
- [Output Formats](@ref) - Result file formats
- [API Workflow](@ref) - Using Julia API for custom analysis
