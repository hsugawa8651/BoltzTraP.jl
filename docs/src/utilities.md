# Utility Functions

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

### Using with TransportResult

```julia
using BoltzTraP

# Run transport calculation
transport = run_integrate("si_interp.jld2"; temperatures=[300.0])

# Use DOS from transport result
cv = calc_cv(
    transport.epsilon,
    transport.dos,
    transport.mu_values,
    transport.temperatures;
    dosweight=2.0
)
```

### Output

| Field | Type | Description |
|-------|------|-------------|
| `cv` | `Matrix{Float64}` | Heat capacity (nT × nμ) in J/K |

---

## Scissor Correction (v0.2)

### apply_scissor

!!! note "Coming in v0.2"
    This function will be available in v0.2.

Apply a scissor correction to shift conduction bands relative to valence bands.

```julia
# Planned API
interp_corrected = apply_scissor(interp, delta_gap)
```

### Use Case

DFT calculations often underestimate band gaps. The scissor correction shifts
conduction band energies by a constant value to match experimental band gaps.

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
