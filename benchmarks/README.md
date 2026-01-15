# Benchmarks

End-to-end performance benchmarks comparing BoltzTraP.jl (Julia) and BoltzTraP2 (Python).

## Overview

The benchmark measures the complete workflow:
1. **Interpolation** (`fitde3D` / `run_interpolate`) - FFT-based band structure fitting
2. **Integration** (`getBTPbands` + transport coefficients) - Transport property calculation

I/O time is excluded from measurements to focus on computational performance.

## Input Data

- **Material**: Silicon (Si)
- **Source**: `benchmarks/data/Si.vasp/` (from [BoltzTraP2 data](https://gitlab.com/souza-group/boltztrap2/-/tree/master/data/Si.vasp), GPL-3.0)
- **Files**: `vasprun.xml`, `POSCAR`

## Scripts

| Script | Description |
|--------|-------------|
| `scaling_benchmark.py` | Python BoltzTraP2 benchmark |
| `scaling_benchmark.jl` | Julia BoltzTraP.jl benchmark |
| `generate_figure.jl` | Generate comparison figure |

## Quick Start

### Local Execution

```bash
# Python benchmark
python benchmarks/scaling_benchmark.py

# Julia benchmark
julia --project benchmarks/scaling_benchmark.jl

# Generate comparison figure (requires both JSON results)
julia --project benchmarks/generate_figure.jl
```

### GitHub Actions

The benchmark runs automatically on:
- Release publication
- Manual trigger via `workflow_dispatch`

```bash
# Manual trigger (CLI)
gh workflow run benchmarks.yml

# With PR creation
gh workflow run benchmarks.yml -f create_pr=true
```

## Requirements

### Python
- BoltzTraP2: `pip install boltztrap2` or `conda install boltztrap2`
- numpy
- pyFFTW (optional, for FFTW backend): `pip install pyfftw`

### Julia
- BoltzTraP.jl (this package)
- JSON.jl, PythonPlot.jl: `] add JSON PythonPlot`

## FFT Backend

| Language | Default | With pyFFTW |
|----------|---------|-------------|
| Python | numpy.fft | FFTW (via pyFFTW) |
| Julia | FFTW.jl | FFTW.jl |

For fair comparison in CI, both use FFTW (pyFFTW installed in GitHub Actions).
Local benchmarks without pyFFTW will use numpy.fft.

## Parameters

| Parameter | Value |
|-----------|-------|
| Temperature | 300 K |
| K-points | 1000, 2000, 4000 |
| DOS bins | 500 |
| Runs | 5 (median reported) |

## Output Files

| File | Content |
|------|---------|
| `scaling_results_python.json` | Python benchmark results |
| `scaling_results_julia.json` | Julia benchmark results |
| `benchmark.png` | End-to-end performance scaling figure |
| `benchmark_breakdown.png` | Performance breakdown (interpolate/integrate) |

## JSON Format

```json
{
  "language": "Julia",
  "version": "1.11.x",
  "timestamp": "2026-01-14T12:00:00Z",
  "parameters": {
    "temperature_K": 300,
    "dos_bins": 500,
    "n_runs": 5
  },
  "results": [
    {
      "kpoints": 1000,
      "total_ms": 1700,
      "interpolate_ms": 500,
      "integrate_ms": 1200,
      "std_ms": 50,
      "peak_memory_MiB": 256
    }
  ]
}
```
