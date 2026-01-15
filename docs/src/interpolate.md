# Interpolation

Band structure interpolation using Fourier series fitting.

## run_interpolate

```@docs
run_interpolate
```

### Example

```julia
using BoltzTraP

# Basic usage
result = run_interpolate("./Si.vasp")

# With options
result = run_interpolate("./Si.vasp";
    kpoints = 10000,     # Target equivalence classes
    emin = -0.5,         # Energy filter (Ha, relative to Fermi)
    emax = 0.5,
    verbose = true
)

# Save to file
save_interpolation("si_interp.jld2", result)
```

---

## Input Data

[`run_interpolate`](@ref) accepts:
- Directory path (auto-detects format)
- [`DFTData`](@ref) from loaders ([`load_vasp`](@ref), [`load_qe`](@ref), [`load_wien2k`](@ref), [`load_gene`](@ref), [`load_abinit`](@ref), [`load_dftk`](@ref))

See [Input Formats](@ref) for loader documentation and supported file formats.

### Example

```julia
using BoltzTraP

# Option 1: Directory path (auto-detect format)
result = run_interpolate("./Si.vasp")

# Option 2: Load data explicitly
data = load_vasp("./Si.vasp")
result = run_interpolate(data; source="./Si.vasp")

# Access DFTData fields (all in atomic units: Bohr, Hartree)
data.lattice       # 3×3 matrix (Bohr)
data.kpoints       # 3×nk matrix (fractional)
data.ebands        # nbands×nk×nspin array (Hartree)
data.fermi         # Fermi energy (Hartree)
data.positions     # 3×natoms matrix (fractional)
data.species       # Vector of element symbols
```

---

## Output: save/load

```@docs
save_interpolation
load_interpolation
```

### Example

```julia
# Save interpolation result
save_interpolation("si_interp.jld2", result)

# Load for later use
loaded = load_interpolation("si_interp.jld2")
```

### Supported Formats

| Extension | Format | Description |
|-----------|--------|-------------|
| `.jld2` | JLD2 | Julia native (HDF5-based, fast) |
| `.bt2` | BT2 | Python BoltzTraP2 compatible (JSON + LZMA) |

---

## Working with InterpolationResult

See [`InterpolationResult`](@ref) in [Output Formats](@ref) for field details.

### Accessing Fields

```julia
result = run_interpolate("./Si.vasp")

# Access fields
println("Number of bands: ", size(result.coeffs, 1))
println("Number of equivalences: ", length(result.equivalences))
println("Fermi energy: ", result.metadata["fermi_energy"], " Ha")

# Coefficients for band 1
band1_coeffs = result.coeffs[1, :]
```

---

## Next Steps

- [Integration](@ref) - Compute transport coefficients
- [Output Formats](@ref) - Result file formats
