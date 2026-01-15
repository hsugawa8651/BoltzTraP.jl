# Validation

Visual comparison tools for validating BoltzTraP.jl against Python BoltzTraP2.

## Prerequisites

1. Generate Python reference data:
```bash
cd reftest
python generate_5_e2e.py
```

2. Generate Julia transport results:
```bash
cd BoltzTraP.jl
boltztrap interpolate reftest/data/vasp_si.npz -k 5000 -o si_interp.jld2
boltztrap integrate si_interp.jld2 -t 300 -o si_transport.jld2
```

## Usage

```bash
julia --project validation/compare_transport.jl <python_npz> <julia_jld2> [options]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `<python_npz>` | Path to Python BoltzTraP2 results (npz format from reftest) |
| `<julia_jld2>` | Path to Julia BoltzTraP.jl transport results (jld2 format) |

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output <file>` | Output figure path | `validation/transport_<title>_<T>K.png` |
| `--title <name>` | Material name for figure title | `untitled` |
| `--xlims <min,max>` | X-axis limits in eV | `-0.5,0.5` |
| `-t, --temperature <T>` | Temperature in K | `300` |
| `--alltemp` | Generate figures for all common temperatures | - |
| `-h, --help` | Show help | - |

## Examples

### Generate validation figures

```bash
cd BoltzTraP.jl

# Silicon (all temperatures)
julia --project validation/compare_transport.jl \
    reftest/data/si_end2end.npz \
    si_transport.jld2 \
    --title Si --alltemp

# PbTe (all temperatures)
julia --project validation/compare_transport.jl \
    reftest/data/pbte_end2end.npz \
    pbte_transport.jld2 \
    --title PbTe --alltemp

```

### JOSS Paper Figures

For JOSS paper submission, generate figures with specific output names:

```bash
# Figure 1: Silicon transport coefficients
julia --project validation/compare_transport.jl \
    reftest/data/si_end2end.npz \
    si_transport.jld2 \
    --title Si \
    -o paper/figure1.png

# PbTe transport coefficients
julia --project validation/compare_transport.jl \
    reftest/data/pbte_end2end.npz \
    pbte_transport.jld2 \
    --title PbTe \
    -o paper/figure_pbte.png
```

## Output

Figures are saved as `validation/transport_<title>_<T>K.png`:

| Material | 200K | 300K | 400K |
|----------|------|------|------|
| Si | `transport_Si_200K.png` | `transport_Si_300K.png` | `transport_Si_400K.png` |
| PbTe | `transport_PbTe_200K.png` | `transport_PbTe_300K.png` | `transport_PbTe_400K.png` |

Custom path can be specified with `-o` option.

## Figure Contents

3-panel vertical layout showing:
1. **Seebeck coefficient** $S_{xx}$ (μV/K)
2. **Electrical conductivity** $\sigma_{xx}/\tau$ (1/Ωms, log scale)
3. **Thermal conductivity** $\kappa_{xx}/\tau$ (W/mKs, log scale)

Legend shows:
- Red circles: Python BoltzTraP2
- Blue crosses: Julia BoltzTraP.jl

## Test Data Sources

Validation uses data from [BoltzTraP2](https://gitlab.com/souza-group/BoltzTraP2) test materials.

Additional loader tests use synthetic structures generated with [pymatgen](https://pymatgen.org/) (MIT), [ASE](https://wiki.fysik.dtu.dk/ase/) (LGPL-2.1), and [NumPy](https://numpy.org/) (BSD-3). These generated files are included in the distribution. See `reftest/README.md` for details.
