# DFTK.jl Integration

This is a variant of the [API Workflow](@ref) that enables a complete Julia-only workflow from first-principles electronic structure calculations to Boltzmann transport coefficients.

By combining [DFTK.jl](https://dftk.org/) (Julia-native DFT) with BoltzTraP.jl, the entire workflow runs in a single Julia session without external dependencies.

```@docs
load_dftk
```

---

## Workflow

The input is the [SCF result](https://docs.dftk.org/stable/guide/self_consistent_field/) (`scfres`) from DFTK's `self_consistent_field` function.

```
DFTK.jl                    BoltzTraP.jl
┌──────────┐              ┌──────────────┐
│ scfres   │──────────────│ load_dftk    │
│ (bands,  │   Extension  │ (extract)    │
│  kpts)   │              └──────┬───────┘
└──────────┘                     │
                                 ▼
                    ┌────────────────────┐
                    │ run_interpolate    │
                    └────────┬───────────┘
                             ▼
                    ┌────────────────────┐
                    │ run_integrate      │
                    └────────┬───────────┘
                             ▼
                    ┌────────────────────┐
                    │ TransportResult    │
                    └────────────────────┘
```

### Example

```julia
using DFTK
using BoltzTraP

# Setup DFTK calculation
a = 10.26  # Silicon lattice constant (Bohr)
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
Si = ElementPsp(:Si; psp=load_psp("hgh/lda/si-q4"))
atoms = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

model = model_LDA(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut=15, kgrid=[10, 10, 10])
scfres = self_consistent_field(basis)

# Extract data for BoltzTraP
data = load_dftk(scfres)

# Run interpolation and transport
interp = run_interpolate(data; source="DFTK", kpoints=5000, verbose=true)
transport = run_integrate(interp; temperatures=[300.0])

# Access results
println("Seebeck at 300K: ", transport.seebeck[1, 1, 1, :])
```

!!! warning "Load order"
    `using DFTK` must be called **before** `using BoltzTraP` to activate the extension.

### Notes

- DFTK uses atomic units (Hartree, Bohr) internally - no conversion needed
- Dense k-grid (10x10x10 or more) recommended for accurate interpolation
- Only non-spin-polarized calculations are supported

---

## Testing

BoltzTraP.jl uses Julia package extensions for optional dependencies. DFTK extension requires special handling when running tests.

### Unit Tests

Unit tests are in `test/test_dftk_extension.jl`:

```bash
# DFTK tests are slow (require SCF), so they're gated by environment variable
TEST_DFTK=true julia --project -e 'using Pkg; Pkg.test()'
```

### Manual Testing

Test in REPL:

```julia
# Must load DFTK BEFORE BoltzTraP to trigger extension
using DFTK
using BoltzTraP

# Now load_dftk is available
@assert hasmethod(load_dftk, Tuple{Any})
```

### End-to-End Test

```bash
julia --project ftest/test_dftk_e2e.jl
```

This runs a complete workflow: DFTK SCF → [`load_dftk`](@ref) → [`run_interpolate`](@ref) → [`run_integrate`](@ref).
