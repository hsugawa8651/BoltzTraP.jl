# [Utility Functions](@id utilities)

Additional functions for band structure analysis and transport calculations.

---

## Electronic Heat Capacity

### calc_cv

Compute the electronic contribution to the heat capacity.

```@docs
calc_cv
```

### Example

```julia
using BoltzTraP

# Load interpolation result
interp = load_interpolation("si_interp.jld2")

# Get DOS from interpolation
epsilon = interp.epsilon
dos = interp.dos

# Define chemical potential and temperature ranges
mu_range = range(-0.1, 0.1, length=101)  # Ha
T_range = [100.0, 200.0, 300.0, 400.0, 500.0]  # K

# Calculate electronic heat capacity
cv = calc_cv(epsilon, dos, mu_range, T_range; dosweight=2.0)
# Returns: (nT, nμ) matrix in J/K
```

### Using with TransportResult (Recommended)

```julia
using BoltzTraP

# Run transport calculation
transport = run_integrate("si_interp.jld2"; temperatures=[300.0])

# Convenience method - extracts all parameters automatically
cv = calc_cv(transport)
```

All required parameters (`epsilon`, `dos`, `mu_values`, `temperatures`, `dosweight`) are automatically extracted from the [`TransportResult`](@ref).

### Output

| Field | Type | Description |
|-------|------|-------------|
| `cv` | `Matrix{Float64}` | Heat capacity (nT × nμ) in J/K |

---

## Scissor Correction

### apply_scissor

Apply a scissor correction to shift conduction bands relative to valence bands.

```@docs
apply_scissor
```

### Use Case

DFT calculations often underestimate band gaps. The scissor correction shifts
conduction band energies by a constant value to match experimental band gaps.

### High-Level API (Recommended)

Use the `scissor` parameter in [`run_integrate`](@ref):

```julia
using BoltzTraP

# Apply scissor correction to achieve 1.1 eV band gap
transport = run_integrate("si_interp.jld2";
    temperatures = [300.0],
    scissor = 1.1  # Target gap in eV
)

# Check if scissor was applied
if haskey(transport.metadata, "scissor_eV")
    println("Scissor applied: $(transport.metadata["scissor_eV"]) eV")
end
```

### CLI Usage

```bash
boltztrap integrate si_interp.jld2 -t 300 --scissor 1.1
```

### Low-Level API

For direct manipulation of band energies:

```julia
using BoltzTraP

# Apply scissor to raw band data
eband_shifted = apply_scissor(epsilon, dos, nelect, eband, desired_gap;
    dosweight = 2.0
)
```

Where `desired_gap` is in Hartree. Use `EV_TO_HA` for conversion from eV.

---

## DOS Smoothing (v0.3)

### smoothen_DOS

!!! note "Coming in v0.3"
    This function will be available in v0.3.

Apply Gaussian smoothing to the density of states.

---

## Validation

To validate `calc_cv` against Python BoltzTraP2:

```bash
julia --project validation/compare_cv.jl
```

See [Validation](@ref) for details.
