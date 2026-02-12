# BoltzTraP.jl Examples

Last-Modified: 2026-02-12

## Numbering System

3-digit `xyz` format:
- `x`: Category (physics / method)
- `y`: Material (y=0 reserved)
- `z`: Sequential number within material

| x | Category | Version |
|---|----------|---------|
| 0 | Documentation (this index) | v0.3 |
| 1 | Non-magnetic, DFT direct | v0.3 |
| 3 | Magnetic, DFT direct | v0.3 |
| 5 | Wannier.jl integration | v0.4 (planned) |
| 7 | Visualization | v0.3 |
| 9 | Validation / Benchmarks | v0.3 |

**Material Codes** (2nd digit, y=0 reserved):

| xy | Material | Notes |
|----|----------|-------|
| 11z | Si | Diamond, Fd-3m — all DFT formats |
| 12z | PbTe | Rock salt, Fm-3m — thermoelectric |
| 13z | Li | BCC, Im-3m |
| 15z | Bi2Te3 | Topological insulator, SOC |
| 16z | CoSb3 | Skutterudite |
| 31z | Fe | BCC, Im-3m — collinear magnetic |
| 32z | CrI3 | Antiferromagnetic |

**DFT Format Codes** (3rd digit, for loader examples):

| z | Format |
|---|--------|
| 0 | VASP |
| 1 | Quantum ESPRESSO |
| 2 | Wien2k |
| 3 | GENE |
| 4 | ABINIT |
| 5 | DFTK.jl |
| 7+ | Applied features (scissor, cv, etc.) |

## File List

### 1yz: Non-magnetic

| File | Description |
|------|-------------|
| `110_si_vasp.jl` | Si VASP — basic workflow (interpolation → transport) |

### 3yz: Magnetic

(planned)

### 7xx: Visualization

(planned)

### 9xx: Validation

(planned)

## How to Run

```bash
cd BoltzTraP.jl

# Run a single example
julia --project=. examples/110_si_vasp.jl
```

## Dependencies

All examples use:
- `BoltzTraP` (this package)

Some examples additionally require:
- `DFTK` (for DFTK.jl integration examples, 115/315)
- `Plots` (for visualization examples, 7xx)

## See Also

- [Workflow Guide](https://hsugawa8651.github.io/BoltzTraP.jl/workflow/)
- [Magnetic Materials](https://hsugawa8651.github.io/BoltzTraP.jl/magnetic/)
- [Input Formats](https://hsugawa8651.github.io/BoltzTraP.jl/input_formats/)
