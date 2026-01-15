# Conventions

BoltzTraP.jl follows consistent conventions for arrays and physical units.

---

## Array Conventions

### K-points

K-points are stored as **3×nk matrices** (column vectors):

```julia
kpoints = [
    kx₁ kx₂ kx₃ ...
    ky₁ ky₂ ky₃ ...
    kz₁ kz₂ kz₃ ...
]  # 3×nk
```

### Band Energies

Band energies are stored as **nbands×nk×nspin arrays**:

```julia
ebands = zeros(nbands, nk, nspin)
ebands[iband, ik, ispin] = energy
```

### Transport Tensors

Transport tensors (σ, S, κ) are stored as **3×3×nT×nμ arrays**:

```julia
sigma[i, j, iT, imu]  # σ_ij at temperature iT and chemical potential imu
```

To get the trace:
```julia
trace_sigma = sigma[1,1,:,:] + sigma[2,2,:,:] + sigma[3,3,:,:]
```

To get the average:
```julia
avg_sigma = (sigma[1,1,:,:] + sigma[2,2,:,:] + sigma[3,3,:,:]) / 3
```

---

## Unit System

BoltzTraP.jl follows Python BoltzTraP2 conventions for compatibility.

### Summary by Module

| Module | Energy | Length | Notes |
|--------|--------|--------|-------|
| **I/O (VASP)** | eV | Å | Native VASP format |
| **Interpolation** | - | - | Unit-agnostic |
| **Transport** | Hartree | Bohr | Internal atomic units |
| **CLI (emin/emax)** | Hartree | - | BoltzTraP2 compatible |
| **Output** | Preserves input | Preserves input | |

### Interpolation: Unit-Agnostic

The interpolation module performs pure mathematical operations (Fourier fitting, phase factors, band reconstruction) independent of physical units.

```julia
# If input is in eV (from VASP), output is in eV
result = run_interpolate("./Si.vasp")
```

### Transport: Atomic Units

The transport module uses atomic units internally (ℏ = e = mₑ = kB = 1).

Output is converted to SI units:
- σ: S/m (electrical conductivity)
- S: V/K (Seebeck coefficient)
- κ: W/(m·K) (thermal conductivity)

### CLI Parameters

For Python BoltzTraP2 compatibility, `--emin` and `--emax` are in **Hartree** relative to Fermi:

```bash
# Energy window: Fermi ± 0.5 Ha (≈ ±13.6 eV)
boltztrap interpolate ./Si.vasp --emin -0.5 --emax 0.5
```

### Conversion Reference

| Quantity | Atomic → SI |
|----------|-------------|
| Energy | 1 Hartree (Ha) = 27.211 eV |
| Length | 1 Bohr = 0.529 Å |
| Temperature | 1 Ha/kB = 315,775 Kelvin (K) |

---

## Type Summary

| Type | Description |
|------|-------------|
| [`InterpolationResult`](@ref) | Fourier coefficients and metadata |
| [`TransportResult`](@ref) | Transport tensors and DOS |
| [`DFTData`](@ref) | Abstract type for DFT input data |
| [`NonMagneticData`](@ref) | Non-magnetic (unpolarized) DFT data |
| [`SpinPolarizedData`](@ref) | Spin-polarized DFT data (defined, full support in v0.2) |

---

## DFT Data Types

```@docs
DFTData
NonMagneticData
SpinPolarizedData
nspin
is_magnetic
```

### Usage

```julia
using BoltzTraP

# Load DFT data (returns NonMagneticData for unpolarized calculations)
data = load_vasp("./Si.vasp")

# Check spin type
nspin(data)        # Returns 1 for non-magnetic
is_magnetic(data)  # Returns false for non-magnetic

# Type dispatch
function process(data::NonMagneticData)
    # Handle non-magnetic case
end
```

!!! note "v0.2: Magnetic Material Support"
    `SpinPolarizedData` is defined but not yet fully implemented.
    Full support for spin-polarized calculations is planned for v0.2.
