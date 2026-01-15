# File Formats

BoltzTraP.jl uses specific file formats for input and output.

## Input Formats

### Format Detection

See [`detected_format`](@ref) for automatic format detection.

**Example:**
```julia
format = detected_format("./calculation")  # Returns "VASP", "QE", or nothing
```

### VASP

BoltzTraP.jl reads VASP output files from a directory:

| File | Required | Description |
|------|----------|-------------|
| `vasprun.xml` | Yes | Band energies, k-points, lattice, Fermi energy |
| `POSCAR` | No | Crystal structure (also in vasprun.xml) |
| `EIGENVAL` | No | Alternative band energy source |

**Directory structure:**
```
Si.vasp/
├── vasprun.xml    # Main data source
├── POSCAR         # Optional
└── EIGENVAL       # Optional
```

### Quantum ESPRESSO (experimental)

Support for Quantum ESPRESSO output is experimental:

| File | Description |
|------|-------------|
| `*.xml` | Data file from `pw.x` |
| `bands.dat.gnu` | Band data from `bands.x` |

---

## Output Formats

### Interpolation Result (`.jld2`)

The `interpolate` command produces a JLD2 file containing:

| Key | Type | Description |
|-----|------|-------------|
| `coeffs` | `Matrix{ComplexF64}` | Fourier coefficients (nbands x nequiv) |
| `equivalences` | `Vector{Matrix{Int}}` | Equivalence class lattice points |
| `lattvec` | `Matrix{Float64}` | Lattice vectors (3x3, Bohr) |
| `metadata` | `Dict{String,Any}` | See below |

**Metadata contents:**

| Key | Type | Description |
|-----|------|-------------|
| `fermi_energy` | `Float64` | Fermi energy (Ha) |
| `nelect` | `Float64` | Number of electrons |
| `dosweight` | `Float64` | DOS weight (2.0 for non-spin-polarized) |
| `selected_bands` | `UnitRange{Int}` | Band indices used |
| `source_file` | `String` | Original DFT file |
| `creation_date` | `String` | ISO 8601 timestamp |
| `generator` | `String` | "BoltzTraP.jl" |

**Loading in Julia:**
```julia
using JLD2

data = load("interpolation.jld2")
coeffs = data["coeffs"]
equivalences = data["equivalences"]
fermi = data["metadata"]["fermi_energy"]
```

### Transport Result (`.jld2`)

The `integrate` command produces a JLD2 file containing:

| Key | Type | Description |
|-----|------|-------------|
| `temperatures` | `Vector{Float64}` | Temperatures (K) |
| `mu` | `Vector{Float64}` | Chemical potentials (Ha) |
| `N` | `Matrix{Float64}` | Carrier concentration (nT x nmu) |
| `sigma` | `Array{Float64,4}` | Conductivity/tau (3x3 x nT x nmu) |
| `seebeck` | `Array{Float64,4}` | Seebeck coefficient (3x3 x nT x nmu) |
| `kappa` | `Array{Float64,4}` | Thermal conductivity/tau (3x3 x nT x nmu) |
| `L0`, `L1`, `L2` | `Array{Float64,4}` | Fermi integrals |
| `mu0` | `Vector{Float64}` | Intrinsic chemical potential at each T |
| `epsilon` | `Vector{Float64}` | DOS energy grid (Ha) |
| `dos` | `Vector{Float64}` | Density of states |
| `vvdos` | `Array{Float64,3}` | Velocity-velocity DOS |
| `metadata` | `Dict{String,Any}` | Source info, creation date |

**Loading in Julia:**
```julia
using JLD2

transport = load("transport.jld2")
sigma = transport["sigma"]      # 3x3 x nT x nmu
seebeck = transport["seebeck"]  # 3x3 x nT x nmu
```

### CSV Export

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
| `sigma_xx`, `sigma_yy`, `sigma_zz` | S/(m*s) | Conductivity/tau diagonal |
| `S_xx`, `S_yy`, `S_zz` | V/K | Seebeck coefficient diagonal |
| `kappa_xx`, `kappa_yy`, `kappa_zz` | W/(m*K*s) | Thermal conductivity/tau diagonal |

---

## Compatibility with Python BoltzTraP2

BoltzTraP.jl supports reading Python BoltzTraP2's `.bt2` interpolation files,
enabling seamless migration from existing Python workflows.

| Operation | Python `.bt2` | Julia `.jld2` |
|-----------|---------------|---------------|
| Read interpolation | ✅ Supported | ✅ Native |
| Write interpolation | ✅ Supported | ✅ Native (default) |
| Read transport | ❌ Not supported | ✅ Native |
| Write transport | ❌ Not supported | ✅ Native |

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

- [Interpolation](@ref) - `InterpolationResult` and save/load functions
- [Integration](@ref) - `IntegrateResult` and save/load functions
