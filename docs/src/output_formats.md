# [Output Formats](@id output-formats)

BoltzTraP.jl produces result files in JLD2 format (Julia native) and supports Python BoltzTraP2's `.bt2` format for compatibility.

## Inspecting Result Files

Use the `describe` command to inspect result files:

```@docs
describe
```

**API usage:**
```julia
using BoltzTraP

# Describe interpolation result
interp = run_interpolate("./Si.vasp")
describe(interp)

# Describe transport result
transport = run_integrate(interp; temperatures=[300.0])
describe(transport)

# Describe from file
describe("si_interp.jld2")
describe("si_transport.jld2")
```

**CLI usage:**
```bash
# Describe interpolation result
boltztrap describe interpolation.jld2

# Describe transport result
boltztrap describe transport.jld2

# Describe Python .bt2 file
boltztrap describe python_result.bt2
```

---

## Interpolation Result (`.jld2`)

```@docs
InterpolationResult
```

**Loading in Julia:**
```julia
using BoltzTraP

# Using convenience function
result = load_interpolation("interpolation.jld2")

# Or using JLD2 directly
using JLD2
data = load("interpolation.jld2")
coeffs = data["coeffs"]
equivalences = data["equivalences"]
fermi = data["metadata"]["fermi_energy"]
```

**Saving:**
```julia
using BoltzTraP

interp = run_interpolate("./Si.vasp"; kpoints=5000)
save_interpolation("si_interp.jld2", interp)
```

---

## Transport Result (`.jld2`)

```@docs
TransportResult
```

**Loading in Julia:**
```julia
using BoltzTraP

# Using convenience function
transport = load_integrate("transport.jld2")

# Access results
println("Temperatures: ", transport.temperatures)
println("Seebeck at 300K: ", transport.seebeck[1, 1, 1, :])

# Or using JLD2 directly
using JLD2
data = load("transport.jld2")
sigma = data["sigma"]      # 3x3 x nT x nmu
seebeck = data["seebeck"]  # 3x3 x nT x nmu
```

**Saving:**
```julia
using BoltzTraP

transport = run_integrate(interp; temperatures=[300.0])
save_integrate("si_transport.jld2", transport)
```

---

## CSV Export

Transport results can be exported to CSV for external analysis:

```julia
using BoltzTraP

transport = load_integrate("transport.jld2")
save_integrate_csv("transport.csv", transport)
```

**CSV columns:**

| Column | Unit | Description |
|--------|------|-------------|
| `T` | K | Temperature |
| `mu` | Ha | Chemical potential |
| `N` | - | Carrier concentration |
| `sigma_xx`, `sigma_yy`, `sigma_zz` | S/(m·s) | σ/τ (electrical conductivity / relaxation time) |
| `S_xx`, `S_yy`, `S_zz` | V/K | S (Seebeck coefficient) |
| `kappa_xx`, `kappa_yy`, `kappa_zz` | W/(m·K·s) | κ₀/τ (electronic thermal conductivity / relaxation time) |

---

## Compatibility with Python BoltzTraP2

BoltzTraP.jl supports reading and writing Python BoltzTraP2's `.bt2` interpolation files, enabling seamless migration from existing Python workflows.

| Operation | Python `.bt2` | Julia `.jld2` |
|-----------|---------------|---------------|
| Read interpolation | Supported | Native |
| Write interpolation | Supported | Native (default) |
| Read transport | Not supported | Native |
| Write transport | Not supported | Native |

**Usage examples:**

```bash
# Read Python .bt2 and compute transport
boltztrap integrate python_result.bt2 -t 300

# Describe .bt2 file contents
boltztrap describe python_result.bt2

# Plot bands from .bt2 (requires --kpath for files without spacegroup)
boltztrap plotbands python_result.bt2 --kpath "G:0,0,0;X:0.5,0,0.5|G-X-G"

# Write interpolation in .bt2 format
boltztrap interpolate ./Si.vasp -o si_interp.bt2
```

**Format details:**

| Format | Structure | Use case |
|--------|-----------|----------|
| `.bt2` | LZMA-compressed JSON | Python compatibility |
| `.jld2` | HDF5-based binary | Julia native (recommended) |

!!! note "Transport results"
    Transport calculation results (`.btj` in Python) are only available in
    `.jld2` format. This is intentional as Python BoltzTraP2 also uses
    plain text formats (`.trace`, `.condtens`) for transport output.

---

## File Size Guidelines

| Material | k-points | Interpolation | Transport |
|----------|----------|---------------|-----------|
| Si (2 atoms) | 5000 | ~1 MB | ~5 MB |
| Complex (50 atoms) | 10000 | ~50 MB | ~20 MB |

---

## Next Steps

- [Input Formats](@ref input-formats) - DFT input formats and loaders
- [Interpolation](@ref) - [`InterpolationResult`](@ref) usage
- [Integration](@ref) - [`TransportResult`](@ref) usage
