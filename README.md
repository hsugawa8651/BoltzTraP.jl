# BoltzTraP.jl

[![CI](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/hsugawa8651/BoltzTraP.jl/actions/workflows/docs.yml/badge.svg)](https://hsugawa8651.github.io/BoltzTraP.jl)
[![codecov](https://codecov.io/gh/hsugawa8651/BoltzTraP.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hsugawa8651/BoltzTraP.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18253185.svg)](https://doi.org/10.5281/zenodo.18253185)

Julia port of [BoltzTraP2 v25.11.1](https://pypi.org/project/BoltzTraP2/25.11.1/) ([source](https://gitlab.com/sousaw/BoltzTraP2)), a band structure interpolation and transport coefficient calculator.

## Features

| Material Type | Support |
|--------------|---------|
| Non-magnetic | ✅ Full |
| Collinear magnetic | ✅ Full (v0.3) |
| Non-collinear magnetic | ⬜ Planned |

* Electronic heat capacity (`calc_cv`) and scissor correction (`apply_scissor`)

### Available Commands

| Feature | CLI | API | Input | Output |
|---------|-----|-----|-------|--------|
| Band interpolation | `boltztrap interpolate` | `run_interpolate()` | VASP, QE, Wien2k, GENE, ABINIT, DFTK (`scfres`) | `InterpolationResult` (`.jld2`, `.bt2`) |
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

## Quick Start

```julia
using BoltzTraP

# Non-magnetic: Si (VASP)
interp = run_interpolate("./Si.vasp"; kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0, 500.0])

# Collinear magnetic: Fe (QE)
# No code changes needed — spin polarization is detected automatically
interp = run_interpolate("./Fe.qe"; kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0, 1000.0])
```

## Validation

BoltzTraP.jl is validated against Python BoltzTraP2 through 75 reference tests covering symmetry, interpolation, transport, I/O, and end-to-end workflows (< 10⁻⁶ relative error). To run validation tests locally:

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

If you use BoltzTraP.jl in your research, please cite the following:

**BoltzTraP.jl:**
> Sugawara, H. (2026). BoltzTraP.jl: Julia port of BoltzTraP2 for band structure interpolation and transport coefficient calculation (Version 0.4.1) [Computer software]. https://doi.org/10.5281/zenodo.20469759

**Original BoltzTraP2:**
> Madsen, G. K., Carrete, J., & Verstraete, M. J. (2018). BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients. *Computer Physics Communications*, 231, 140-145. [doi:10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

**BoltzWann (Wannier-basis transport methodology, reimplemented in BoltzTraP.jl v0.4+):**
> Pizzi, G., Volja, D., Kozinsky, B., Fornari, M., & Marzari, N. (2014). BoltzWann: A code for the evaluation of thermoelectric and electronic transport properties with a maximally-localized Wannier functions basis. *Computer Physics Communications*, 185, 422-429. [doi:10.1016/j.cpc.2013.09.015](https://doi.org/10.1016/j.cpc.2013.09.015)

**Wannier90 v3 (MLWF construction, used via Wannier.jl):**
> Pizzi, G., et al. (2020). Wannier90 as a community code: new features and applications. *Journal of Physics: Condensed Matter*, 32, 165902. [doi:10.1088/1361-648X/ab51ff](https://doi.org/10.1088/1361-648X/ab51ff)

**Wannier.jl (Julia MLWF construction):**
> Qiao, J. Wannier.jl: Julia package for Wannier functions [Computer software]. https://github.com/qiaojunfeng/Wannier.jl

## License

GPL-3.0-or-later (same as original BoltzTraP2)
