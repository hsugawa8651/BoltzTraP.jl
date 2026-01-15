# Plotting

BoltzTraP.jl provides plotting functions for visualizing band structures and transport coefficients.

!!! note "No additional setup required"
    [`plot_bands`](@ref) and [`plot_transport`](@ref) work directly after `using BoltzTraP`.
    For custom plots, add `using Plots` to access the full Plots.jl API.

---

## Band Structure Plotting

[`plot_bands`](@ref) plots interpolated band structure along high-symmetry k-paths.

### From InterpolationResult

```julia
using BoltzTraP

interp = run_interpolate("./Si.vasp"; kpoints=5000)
plot_bands(interp; emin=-5.0, emax=5.0, output="bands.png")
```

### From File

```julia
# From .jld2 (auto k-path from spacegroup metadata)
plot_bands("si_interp.jld2"; emin=-5.0, emax=5.0, output="bands.png")

# From .bt2 (manual k-path required if no spacegroup)
kpath = (
    points = Dict("G" => [0,0,0], "X" => [0.5,0,0.5], "L" => [0.5,0.5,0.5]),
    paths = [["G", "X", "L", "G"]]
)
plot_bands("si_interp.bt2"; kpath=kpath, emin=-5.0, emax=5.0)
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `npoints` | Int | 100 | K-points per path segment |
| `emin` | Float64 | -1.0 | Minimum energy relative to Fermi [eV] |
| `emax` | Float64 | 1.0 | Maximum energy relative to Fermi [eV] |
| `fermi_line` | Bool | true | Show Fermi level as dashed line |
| `kpath` | NamedTuple | nothing | Manual k-path specification |
| `output` | String | nothing | Save to file if specified |

### Manual K-path Format

For files without spacegroup metadata (e.g., Python .bt2 files), specify k-path manually:

```julia
kpath = (
    points = Dict(
        "G" => [0.0, 0.0, 0.0],      # Gamma
        "X" => [0.5, 0.0, 0.5],      # X
        "L" => [0.5, 0.5, 0.5],      # L
        "W" => [0.5, 0.25, 0.75],    # W
    ),
    paths = [["G", "X", "W", "L", "G"]]  # Path segments
)
```

---

## Transport Coefficient Plotting

[`plot_transport`](@ref) plots transport coefficients vs chemical potential or temperature.

### From TransportResult

```julia
using BoltzTraP

interp = run_interpolate("./Si.vasp"; kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])

# Seebeck vs chemical potential
plot_transport(transport; quantity="seebeck", component="xx", output="seebeck.png")

# Conductivity vs temperature
plot_transport(transport; quantity="sigma", abscissa="T", output="sigma_T.png")
```

### From File

```julia
plot_transport("si_transport.jld2"; quantity="seebeck", component="xx")
plot_transport("si_transport.jld2"; quantity="sigma", abscissa="T")
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `quantity` | String | "seebeck" | Property: seebeck, sigma, kappa |
| `component` | String | "xx" | Tensor component: xx, yy, zz, xy, ... |
| `abscissa` | String | "mu" | X-axis: mu (chemical potential) or T (temperature) |
| `temperature` | Float64 | 300.0 | Temperature for mu plot [K] |
| `mu_index` | Int | 0 | mu index for T plot (0 = auto-select near Fermi) |
| `output` | String | nothing | Save to file if specified |

### Available Quantities

| Quantity | Description | Units |
|----------|-------------|-------|
| `seebeck` or `S` | Seebeck coefficient | uV/K |
| `sigma` or `conductivity` | Electrical conductivity / tau | S/m/s |
| `kappa` or `thermal` | Thermal conductivity / tau | W/m/K/s |

---

## CLI Commands

For command-line usage, see [CLI Workflow](@ref):

```bash
# Band structure
boltztrap plotbands si_interp.jld2 --emin -5 --emax 5 -o bands.png
boltztrap plotbands si_interp.bt2 --kpath "G:0,0,0;X:0.5,0,0.5|G-X-G" -o bands.png

# Transport coefficients
boltztrap plot si_transport.jld2 -q seebeck -c xx -o seebeck.png
boltztrap plot si_transport.jld2 -q sigma --abscissa T -o sigma_T.png
```

---

## API Reference

```@docs
plot_bands
plot_transport
```
