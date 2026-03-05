# Magnetic Materials

BoltzTraP.jl supports spin-polarized (collinear magnetic) calculations.

## Supported Spin Types

| Spin Type | Support | Notes |
|-----------|---------|-------|
| Non-magnetic | Full | Default (dosweight=2) |
| Collinear | Full (v0.3+) | Spin-polarized (dosweight=1) |

## DFT Calculation Requirements

To use collinear magnetic calculations, prepare spin-polarized DFT data:

| DFT Code | Setting |
|----------|---------|
| VASP | `ISPIN=2` |
| Quantum ESPRESSO | `nspin=2` |
| DFTK.jl | `spin_polarization=:collinear` |

## Workflow

**No changes required** - BoltzTraP.jl automatically detects spin polarization from the band structure dimensions.

The workflow is identical to non-magnetic calculations:

```julia
using BoltzTraP

# Load spin-polarized data (automatic detection)
data = load_vasp("./Fe.vasp")

# Run interpolation and transport (same API)
interp = run_interpolate(data; kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])
```

### Verification

You can verify that magnetic data was detected correctly:

```julia
data = load_vasp("./Fe.vasp")

# Check spin type
is_magnetic(data)  # true
nspin(data)        # 2
```

### Result Metadata

[`TransportResult`](@ref) includes spin type information in its `metadata` field:

```julia
transport.metadata["spintype"]   # "Collinear"
transport.metadata["dosweight"]  # 1.0
```

See [Output Formats](@ref) for full `TransportResult` documentation.

## Internal Handling

When spin-polarized data is loaded, BoltzTraP.jl automatically:

1. **Detects spin polarization** from band structure dimensions (`nspin=2`)
2. **Concatenates spin channels** for Fourier interpolation (following Python BoltzTraP2)
3. **Applies time-reversal symmetry** considerations for equivalence classes
4. **Sets `dosweight=1.0`** (vs 2.0 for non-magnetic)

## Validation

Fe collinear magnetic transport has been validated at T=1000K against Python BoltzTraP2.

See [Validation](@ref) for details.

!!! note "Why T=1000K?"
    At lower temperatures (e.g., 300K), the Seebeck coefficient shows oscillations
    (~±250 μV/K) due to the BoltzTraP Fourier interpolation method.
    At T=1000K, thermal broadening naturally smooths these oscillations,
    enabling clean validation.
