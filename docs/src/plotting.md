# Plotting

BoltzTraP.jl provides plotting via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(data)` recipes for `BandPlotData` / `TransportPlotData` |
| PlotsExt | `using Plots` | `plot_bands`, `plot_transport` (full-featured with auto k-path) |
| PythonPlotExt | `using PythonPlot` | `savefig_publication` for publication-quality PDF/PNG/SVG figures |

---

## RecipesBase Recipes

These recipes are activated when any RecipesBase-compatible backend is loaded.

```julia
using BoltzTraP, Plots   # or CairoMakie, etc.
```

### Band Structure — `BandPlotData`

`BandPlotData` is a backend-independent container created internally by `plot_bands`.

```julia
interp = run_interpolate("./Si.vasp"; kpoints=5000)
bpd = BandPlotData(interp)  # created internally by plot_bands
plot(bpd)
```

The recipe draws band lines, a Fermi level (red dashed), and vertical lines at high-symmetry points.

| Attribute | Default |
|-----------|---------|
| ylabel | `"Energy (eV)"` |
| ylims | `(emin, emax)` |
| xticks | high-symmetry labels |
| framestyle | `:box` |
| linecolor | `:black` |
| linewidth | 1.5 |

### Transport Coefficients — `TransportPlotData`

```julia
transport = run_integrate(interp; temperatures=[300.0])
tpd = TransportPlotData(transport)  # created internally by plot_transport
plot(tpd)
```

| Attribute | Default |
|-----------|---------|
| xlabel | from data |
| ylabel | from data |
| title | from data |
| linewidth | 2 |

### Customization

Override any attribute with standard keyword arguments:

```julia
plot(bpd, ylabel="E − E_F (eV)", linecolor=:blue, size=(800, 600))
```

---

## Plots.jl Convenience Functions

!!! note "Requires Plots.jl"
    [`plot_bands`](@ref) and [`plot_transport`](@ref) require `using Plots`.
    They internally create `BandPlotData` / `TransportPlotData` and provide
    additional features such as automatic k-path generation via Brillouin.jl.

### Band Structure Plotting

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

## Publication-quality Figures with PythonPlot

[`savefig_publication`](@ref) writes publication-quality PDF, PNG, or SVG
figures via [PythonPlot.jl](https://github.com/JuliaPy/PythonPlot.jl).
Loading `PythonPlot` activates `BoltzTraPPythonPlotExt`, which handles
matplotlib through `CondaPkg`/`PythonCall` (no manual matplotlib install
required).

### Builder: `build_transport_plot_data`

[`build_transport_plot_data`](@ref) extracts a `TransportPlotData` from a
`TransportResult`, performing label formatting, unit conversion, and input
validation in one place. The same `TransportPlotData` is consumed by the
RecipesBase recipe, `plot_transport` (Plots), and `savefig_publication`
(PythonPlot), so all backends render identical data and labels.

```julia
using BoltzTraP

transport = load_integrate("si_transport.jld2")

tpd = build_transport_plot_data(transport;
    quantity = "seebeck",       # "seebeck"/"S", "sigma"/"conductivity", "kappa"/"thermal"
    component = "xx",           # xx, yy, zz, xy, xz, yz, yx, zx, zy
    abscissa = "mu",            # "mu" or "T"
    temperature = 300.0,        # K — used when abscissa = "mu"
)
```

`build_transport_plot_data` raises `ArgumentError` for unknown quantity,
component, or abscissa values, and when the requested temperature is not
present in `transport.temperatures`.

### `savefig_publication`

```julia
using BoltzTraP, PythonPlot

# Single panel (PDF inferred from the file extension)
savefig_publication(tpd, "seebeck.pdf";
    axis_width_cm = 8.0, axis_height_cm = 5.0)

# 1×3 subplot grid for S, σ, κ vs μ at T = 300 K
tpds = [
    build_transport_plot_data(transport;
        quantity = q, component = "xx", abscissa = "mu", temperature = 300.0)
    for q in ("seebeck", "sigma", "kappa")
]
savefig_publication(tpds, "transport_combined.pdf";
    axis_width_cm = 6.0, axis_height_cm = 4.5,
    layout = (1, 3))
```

Layout, margins, and gaps are specified in centimeters; the figure size in
inches is derived internally. Output format is inferred from the file
extension (`.pdf`, `.png`, `.svg`, ...).

### Demo: Si transport at T = 300 K

The figures below were produced from `examples/110_si_vasp_pbe_paw.jl`
followed by `examples/119_si_vasp_pbe_paw_publication.jl`:

- [Combined S/σ/κ vs μ (1×3 layout)](assets/publication_Si_combined_300K.pdf)
- [Seebeck vs μ (single panel)](assets/publication_Si_seebeck_300K.pdf)

### Band figures

`savefig_publication(::BandPlotData, path; ...)` is also provided. Construct
a `BandPlotData` from interpolated band energies and a k-path before
calling it; an end-to-end builder for band data is planned for a future
release.

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
build_transport_plot_data
savefig_publication
```
