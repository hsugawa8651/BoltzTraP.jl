# BoltzTraP.jl

[![CI](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/docs.yml/badge.svg)](https://hsugawa8651.github.io/BoltzTraP.jl)
[![codecov](https://codecov.io/gh/hsugawa8651/BoltzTraP.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hsugawa8651/BoltzTraP.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18253186.svg)](https://doi.org/10.5281/zenodo.18253186)

Julia port of [BoltzTraP2 v25.11.1](https://pypi.org/project/BoltzTraP2/25.11.1/) ([source](https://gitlab.com/sousaw/BoltzTraP2)), a band structure interpolation and transport coefficient calculator.

## Features (v0.1)

* Nonmagnetic materials (spin-polarized support planned)

### Available Commands

| Feature | CLI | API | Input | Output |
|---------|-----|-----|-------|--------|
| Band interpolation | `boltztrap interpolate` | `run_interpolate()` | VASP, QE, DFTK (`scfres`) | `InterpolationResult` (`.jld2`, `.bt2`) |
| Transport calculation | `boltztrap integrate` | `run_integrate()` | `InterpolationResult` (`.jld2`, `.bt2`) | `TransportResult` (`.jld2`) |
| Plot bands | `boltztrap plotbands` | `plot_bands()` | `InterpolationResult` (`.jld2`, `.bt2`) | PNG, PDF |
| Plot transport | `boltztrap plot` | `plot_transport()` | `TransportResult` (`.jld2`) | PNG, PDF |
| Describe results | `boltztrap describe` | `describe()` | `InterpolationResult` or `TransportResult` | - |

## Installation

```julia
using Pkg
Pkg.add("BoltzTraP")
```

## Documentation

See [hsugawa8651.github.io/BoltzTraP.jl](https://hsugawa8651.github.io/BoltzTraP.jl) for full documentation.

## Validation

BoltzTraP.jl is validated against Python BoltzTraP2 through 75 reference tests. To run validation tests locally:

```bash
# Generate reference data (requires Python BoltzTraP2)
cd reftest
pip install boltztrap2
python generate_1_sphere.py
python generate_2_interpolation_si.py
# ... (see reftest/README.md for all scripts)

# Run Julia tests
cd ..
julia --project -e 'using Pkg; Pkg.test()'

# Run with DFTK extension tests (slow, requires DFTK.jl)
TEST_DFTK=true julia --project -e 'using Pkg; Pkg.test()'
```

See [Developer Guide](https://hsugawa8651.github.io/BoltzTraP.jl/developer/) for details.

## Test Data

This package includes test data from two sources:

| Source | Content | License |
|--------|---------|---------|
| [BoltzTraP2](https://gitlab.com/souza-group/BoltzTraP2) | Real DFT calculations (Si, Li, PbTe, etc.) | GPL-3.0 |
| Synthetic (via [pymatgen](https://pymatgen.org/)) | Low-symmetry structures (monoclinic, triclinic) | Data only |

The synthetic test data was generated using pymatgen (MIT), ASE (LGPL-2.1), and NumPy (BSD-3). See `reftest/README.md` for details.

## Contributing

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/hsugawa8651/BoltzTraP.jl/issues).

## Citation

If you use BoltzTraP.jl in your research, please cite both:

**BoltzTraP.jl:**
> Sugawara, H. (2026). BoltzTraP.jl: Julia port of BoltzTraP2 for band structure interpolation and transport coefficient calculation (Version 0.1.1) [Computer software]. https://doi.org/10.5281/zenodo.18328895

**Original BoltzTraP2:**
> Madsen, G. K., Carrete, J., & Verstraete, M. J. (2018). BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients. *Computer Physics Communications*, 231, 140-145. [doi:10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

## License

GPL-3.0-or-later (same as original BoltzTraP2)
