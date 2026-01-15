# Reference Tests

BoltzTraP.jl provides reference tests (reftest) to validate numerical equivalence with Python BoltzTraP2.

These reference tests were used by the developers during implementation to ensure correctness. The same tests and data generation scripts are provided for users who wish to independently verify the implementation.

## Overview

Reference testing ensures BoltzTraP.jl produces results identical to the original Python implementation:

| Aspect | Description |
|--------|-------------|
| Location | [`reftest/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/reftest) directory |
| Tolerance | < 10⁻⁶ relative error |
| Test count | 75+ reference tests |
| Data source | [BoltzTraP2 v25.11.1](https://gitlab.com/souza-group/BoltzTraP2) |

## What is Tested

Reference tests validate numerical equivalence for:

| Category | What is Compared | Tolerance |
|----------|------------------|-----------|
| Symmetry | Equivalence classes, star functions, rotation matrices | exact |
| Interpolation | Phase factors, Fourier coefficients, reconstructed bands | < 10⁻¹⁰ |
| Transport | DOS, transport DOS, Fermi integrals (L₀, L₁, L₂), Onsager coefficients (σ, S, κ) | < 10⁻⁶ |
| I/O | Lattice vectors, atomic positions, k-points, band energies, Fermi level | exact |
| End-to-end | Full workflow from DFT data to transport coefficients | < 10⁻⁶ |

### Test Materials

| Material | Structure | Type | Purpose |
|----------|-----------|------|---------|
| Si | Diamond (Fd-3m) | Semiconductor | High-symmetry cubic |
| PbTe | Rock salt (Fm-3m) | Thermoelectric | Heavy elements |
| Synthetic | Monoclinic, Triclinic | - | Low-symmetry validation |

### Test Categories

#### Symmetry and Equivalences

Tests for k-point symmetry and equivalence class generation.

| Test | Material | Description |
|------|----------|-------------|
| `unit_cube` | Generic | Unit cube lattice vectors |
| `simple_cubic` | Generic | Simple cubic with atoms |
| `si_diamond` | Si | Diamond structure symmetry |
| `monoclinic` | Generic | Monoclinic space group |
| `triclinic_p1` | Generic | P1 (no symmetry) |

#### Band Interpolation

Tests for Fourier coefficient fitting and band reconstruction.

| Test | Material | k-points | Bands | Equivalences |
|------|----------|----------|-------|--------------|
| `si_interpolation` | Si | 165 | 6 | 1,102 |

#### Transport Coefficients

Tests for Fermi integrals and Onsager coefficients.

| Test | Material | Temperatures | Chemical potentials |
|------|----------|--------------|---------------------|
| `simple_transport` | Si | 200, 300, 400 K | 50 points |

#### File I/O

Tests for reading DFT output files.

| Format | Materials | Test Count |
|--------|-----------|------------|
| VASP | Si, PbTe | 18 |
| QE | Si, Fe, CrI3 | 27 |
| Wien2k | Si, CoSb3, Bi2Te3 | 27 |
| GENE | Si, Synthetic | 18 |
| ABINIT | Si | 9 |

#### End-to-End Tests

Complete workflow tests from DFT data to transport coefficients.

| Test | Material | Description |
|------|----------|-------------|
| `si_end2end` | Si | Full Si workflow |
| `pbte_end2end` | PbTe | Full PbTe workflow |

## Running Reference Tests

For instructions on running reference tests (prerequisites, data generation, execution), see:

**[reftest/README.md](https://github.com/hsugawa8651/BoltzTraP.jl/blob/main/reftest/README.md)**

## See Also

- [Validation](validation.md) - Visual comparison with Python results
- [Benchmarks](benchmarks.md) - Performance comparison
