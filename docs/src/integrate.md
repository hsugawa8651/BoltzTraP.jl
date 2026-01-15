# Integration

Compute transport coefficients by Brillouin zone integration.

## run_integrate

```@docs
run_integrate
```

### Example

```julia
using BoltzTraP

# Basic usage
transport = run_integrate("si_interp.jld2";
    temperatures = [300.0]
)

# Temperature range
transport = run_integrate("si_interp.jld2";
    temperatures = collect(100:50:500),  # 100, 150, ..., 500 K
    bins = 2000,
    verbose = true
)

# Save results
save_integrate("si_transport.jld2", transport)
```

---

## Output: save/load

```@docs
save_integrate
load_integrate
```

### Example

```julia
# Save integration results
save_integrate("si_transport.jld2", transport)

# Load for analysis
loaded = load_integrate("si_transport.jld2")

# Access transport tensors
sigma = loaded.sigma      # 3×3×nT×nμ conductivity
seebeck = loaded.seebeck  # 3×3×nT×nμ Seebeck coefficient
kappa = loaded.kappa      # 3×3×nT×nμ thermal conductivity
```

### Supported Formats

| Extension | Format | Description |
|-----------|--------|-------------|
| `.jld2` | JLD2 | Julia native (full tensor data) |
| `.csv` | CSV | Text format (trace averages only) |

---

## Analyzing Results

### Temperature Dependence

```julia
# Get Seebeck at chemical potential μ=0 for all temperatures
idx_mu0 = argmin(abs.(transport.mu_values))
S_vs_T = transport.seebeck[1, 1, :, idx_mu0] .* 1e6  # V/K → μV/K

using Plots
plot(transport.temperatures, S_vs_T,
    xlabel = "Temperature (K)",
    ylabel = "S_xx (μV/K)",
    title = "Seebeck coefficient"
)
```

### Chemical Potential Dependence

```julia
# Get conductivity at 300K vs chemical potential
T_idx = findfirst(==(300.0), transport.temperatures)
sigma_vs_mu = transport.sigma[1, 1, T_idx, :]

plot(transport.mu_values, sigma_vs_mu,
    xlabel = "μ (eV)",
    ylabel = "σ_xx (S/m)",
    title = "Conductivity at 300K"
)
```

---

## Working with TransportResult

See [`TransportResult`](@ref) in [Output Formats](@ref) for field details.

### Accessing Fields

```julia
transport = run_integrate("si_interp.jld2"; temperatures=[300.0, 400.0, 500.0])

# Access temperatures
println("Temperatures: ", transport.temperatures)  # [300.0, 400.0, 500.0]

# Seebeck coefficient at 300K (first temperature)
S_300K = transport.seebeck[:, :, 1, :]  # 3×3×nμ tensor

# Diagonal components
S_xx = transport.seebeck[1, 1, 1, :]  # S_xx at 300K vs μ
S_yy = transport.seebeck[2, 2, 1, :]  # S_yy at 300K vs μ

# Find Seebeck at Fermi level (μ = 0)
idx = argmin(abs.(transport.mu_values))
S_at_fermi = transport.seebeck[1, 1, 1, idx]
```

---

## Next Steps

- [Output Formats](@ref) - Result file formats and `describe` command
- [Conventions](@ref) - Units and array conventions
