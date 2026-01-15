# API Workflow

This page describes the **API workflow** - using BoltzTraP.jl as a Julia library.

!!! info "CLI vs API Workflow"
    - **CLI workflow**: Run calculations via `boltztrap` shell command
    - **API workflow**: Call Julia functions ([`run_interpolate`](@ref), [`run_integrate`](@ref)) directly

    For CLI workflow, see [CLI Workflow](@ref).

BoltzTraP.jl supports multiple input methods for transport calculations:

| Method | Input | Loader |
|--------|-------|--------|
| **VASP** | `vasprun.xml` | `load_vasp` |
| **Quantum ESPRESSO** | `data-file-schema.xml` | `load_qe` |
| **Wien2k** | `case.struct`, `case.energy` | `load_wien2k` |
| **GENE/Generic** | `case.structure`, `case.energy` | `load_gene` |
| **ABINIT** | `*_GSR.nc` (NetCDF) | `load_abinit` |
| **DFTK.jl** | `scfres` from DFTK.jl | `load_dftk` |

See [Input Formats](@ref input-formats) for detailed format specifications.

---

## Overview

All workflows follow the same pattern:

```
┌───────────────────────────────┐
│         Data Source           │
│ (VASP/QE/Wien2k/GENE/ABINIT/  │
│  DFTK)                        │
└───────────────┬───────────────┘
                │ load_*() → DFTData
                ▼
┌───────────────────────────────┐
│       run_interpolate         │
│        (Fourier fit)          │
└───────────────┬───────────────┘
                │
                ▼
        InterpolationResult
                │
                ▼
┌───────────────────────────────┐
│        run_integrate          │
│    (transport coefficients)   │
└───────────────┬───────────────┘
                │
                ▼
         TransportResult
```

[`run_integrate`](@ref) accepts both [`InterpolationResult`](@ref) directly or a file path:

```julia
# Direct (no intermediate file)
# Note: "./Si.vasp" is a directory containing vasprun.xml and POSCAR
interp = run_interpolate("./Si.vasp")
transport = run_integrate(interp; temperatures=[300.0])

# Via file (for later reuse)
interp = run_interpolate("./Si.vasp")
save_interpolation("si_interp.jld2", interp)
transport = run_integrate("si_interp.jld2"; temperatures=[300.0])
```

---

## VASP Workflow

### Step 1: Interpolate

[`run_interpolate`](@ref) reads DFT data and fits band energies to a Fourier series.

```julia
interp = run_interpolate("./Si.vasp"; kpoints=5000, verbose=true)
```

### Step 2: Integrate

[`run_integrate`](@ref) computes transport coefficients by BZ integration.

```julia
transport = run_integrate(interp; temperatures=[300.0], verbose=true)
```

### Complete Example

```julia
using BoltzTraP

# Direct workflow (no intermediate files)
interp = run_interpolate("./Si.vasp"; kpoints=5000, verbose=true)
transport = run_integrate(interp;
    temperatures = collect(100:100:800),
    verbose = true
)

# Optional: save results for later
save_interpolation("si_interp.jld2", interp)
save_integrate("si_transport.jld2", transport)

# Analyze results
using Plots
T_idx = findfirst(==(300.0), transport.temperatures)
S_xx = transport.seebeck[1, 1, T_idx, :] .* 1e6  # V/K -> uV/K
plot(transport.mu_values, S_xx, xlabel="mu (eV)", ylabel="S_xx (uV/K)")
```

---

## QE Workflow

Similar to VASP, using [`load_qe`](@ref):

```julia
using BoltzTraP

# Load QE data
data = load_qe("./Si.qe")

# Run interpolation
interp = run_interpolate(data; source="./Si.qe", kpoints=5000, verbose=true)

# Compute transport
transport = run_integrate(interp; temperatures=[300.0])
```

---

## DFTK Workflow

See [DFTK.jl Integration](@ref) for detailed documentation.

---

## Generic Data Workflow

For other data sources, prepare a NamedTuple with the required fields:

```julia
using BoltzTraP

# Prepare data (all in atomic units: Hartree, Bohr)
data = (
    lattice = ...,       # 3x3 matrix, columns are lattice vectors (Bohr)
    positions = ...,     # 3xN matrix, fractional coordinates
    species = [...],     # Vector of element symbols
    kpoints = ...,       # 3xM matrix, fractional coordinates
    weights = ...,       # M-vector, k-point weights
    ebands = ...,        # (nbands, M, nspin) array (Hartree)
    occupations = ...,   # (nbands, M, nspin) array
    fermi = ...,         # Fermi energy (Hartree)
    nelect = ...,        # Number of electrons
)

interp = run_interpolate(data; source="custom", kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])
```

---

## Utility Functions

### Format Detection

[`detected_format`](@ref) auto-detects the DFT format in a directory:

```julia
using BoltzTraP

format = detected_format("./calculation")
println("Detected: $format")  # "vasp" or "qe" or nothing
```

This is used internally by [`run_interpolate`](@ref) and [`load_dft`](@ref) for auto-detection.

---

## Next Steps

- [Interpolation](@ref) - `run_interpolate` details
- [Integration](@ref) - `run_integrate` details
- [Plotting](@ref) - Visualization functions
- [Input Formats](@ref input-formats) - DFT input formats and loaders
- [Output Formats](@ref output-formats) - Result file formats
- [Conventions](@ref) - Units and array conventions
