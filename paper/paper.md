---
title: 'BoltzTraP.jl: A Julia package for band structure interpolation and transport coefficient calculations'
tags:
  - Julia
  - condensed matter physics
  - transport coefficients
  - Boltzmann transport
  - thermoelectrics
  - band structure
  - Seebeck coefficient
  - electrical conductivity
authors:
  - name: Hiroharu Sugawara
    orcid: 0000-0002-0071-2396
    corresponding: true
    email: hsugawa@tmu.ac.jp
    affiliation: 1
affiliations:
  - name: Graduate School of Systems Design, Tokyo Metropolitan University, Japan
    index: 1
date: 21 January 2026
bibliography: paper.bib
repository: https://github.com/hsugawa8651/BoltzTraP.jl
archive_doi: 10.5281/zenodo.18328895
---

# Summary

Thermoelectric materials convert heat directly into electricity and vice versa, enabling applications in waste heat recovery, solid-state cooling, and remote power generation. The efficiency of thermoelectric devices depends on the figure of merit ZT, which is determined by the Seebeck coefficient, electrical conductivity, and thermal conductivity. Accurate prediction of these transport coefficients from first-principles calculations is essential for discovering and optimizing new thermoelectric materials.

`BoltzTraP.jl` is a Julia implementation of the BoltzTraP2 algorithm [@madsen2018boltztrap2] for calculating semi-classical transport coefficients from density functional theory (DFT) band structures. The package interpolates band energies using a Fourier expansion in reciprocal space, reconstructs the bands on a dense k-point grid via FFT, and computes transport tensors using the Boltzmann transport equation in the constant relaxation time approximation.

# Statement of Need

BoltzTraP2 [@madsen2018boltztrap2], implemented with C++ and Python, is the de facto standard for computing transport coefficients from DFT calculations, with over 2,000 citations. The algorithm—Fourier interpolation of band structures followed by Boltzmann transport integration—is well-established and widely validated.

`BoltzTraP.jl` is a faithful Julia port that implements the exact same algorithm as BoltzTraP2, producing numerically equivalent results. The motivation for this port is integration with the growing Julia ecosystem for materials science:

- **Julia ecosystem integration**: `BoltzTraP.jl` integrates natively with Julia packages like DFTK.jl [@herbst2021dftk] and Wannier.jl [@qiao2026wannierjl], enabling direct transport calculations from DFTK self-consistent field results without intermediate files or language bridges.
- **Pure Julia, no compilation**: BoltzTraP2 requires compiling a C++ extension, which can fail on some systems—particularly on HPC clusters with non-standard compiler configurations. `BoltzTraP.jl` has no external compiled dependencies.
- **HPC-friendly**: Julia's package manager handles dependencies without conda/pip conflicts common on shared HPC systems. No C++ compilation means no compiler version mismatches. As pure Julia code, `BoltzTraP.jl` is designed to be highly portable and should run on any system where Julia is available, from personal laptops to high-performance computing (HPC) clusters.
- **Performance**: BLAS-optimized algorithms achieve 1.8-3.4x end-to-end speedup over BoltzTraP2 (depending on problem size), enabling rapid screening of thermoelectric materials.


The package is not intended to replace BoltzTraP2 for existing Python/Wien2k workflows, but to serve users who prefer or require a Julia-based toolchain, particularly those working within the Julia DFT ecosystem or on HPC systems where Python environment management is problematic.

# Comparison with BoltzTraP2

`BoltzTraP.jl` implements the same algorithm as BoltzTraP2, with differences in implementation language and ecosystem integration:

| Aspect | BoltzTraP2 (Python) | BoltzTraP.jl |
|--------|:-------------------:|:------------:|
| Algorithm | Fourier interpolation | Identical |
| Numerical results | Reference | Equivalent (< 10⁻⁶ error) |
| Language | Python + C++ | Pure Julia |
| Compilation required | Yes (C++ extension) | No |
| Julia DFT integration | No | Yes (DFTK.jl) |
| Performance (end-to-end) | 1x (baseline) | 1.8-3.4x faster |
| Spin-polarized calculations | Yes | Not yet (planned) |
| Unit handling | Manual conversion | Atomic units (Hartree) |
| Input formats | VASP, QE, Wien2k, GENE, ABINIT | VASP, QE, Wien2k, GENE, ABINIT, DFTK |
| Output format | `.bt2`/`.btj` (LZMA+JSON) | `.jld2` (HDF5, default)[^bt2] |

[^bt2]: BoltzTraP.jl can also read and write `.bt2` files (BoltzTraP2's native format) for compatibility with existing workflows.

# Example Usage

## From DFT Output Files

```julia
using BoltzTraP

# Interpolate band structure from VASP calculation
interp = run_interpolate("./Si.vasp"; kpoints=5000, verbose=true)

# Compute transport coefficients
transport = run_integrate(interp;
    temperatures = [200.0, 300.0, 400.0, 500.0],
    verbose = true
)

# Access results (in constant relaxation time approximation)
# transport.sigma: σ/τ, electrical conductivity / relaxation time (S/m/s)
# transport.seebeck: S, Seebeck coefficient (V/K)
# transport.kappa: κ₀/τ, electronic thermal conductivity / relaxation time (W/m/K/s)
```

## DFTK.jl Integration

```julia
using DFTK
using BoltzTraP

# Run DFT calculation with DFTK
a = 10.26  # Silicon lattice constant (Bohr)
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
Si = ElementPsp(:Si; psp=load_psp("hgh/lda/si-q4"))
atoms = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

model = model_LDA(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut=15, kgrid=[10, 10, 10])
scfres = self_consistent_field(basis)

# Direct transport calculation
data = load_dftk(scfres)
interp = run_interpolate(data; source="DFTK", kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])
```

## Command Line Interface

```bash
# Interpolate band structure
boltztrap interpolate -k 5000 -o si_interp.jld2 ./Si.vasp

# Compute transport at multiple temperatures
boltztrap integrate si_interp.jld2 -t 100:100:500 -o si_transport.jld2
```

# Performance

Benchmarks on a MacBook Pro (Apple M2) comparing end-to-end performance (interpolation + integration) for silicon (1102 equivalence classes, 6 bands):

| k-points | BoltzTraP2 (Python) | BoltzTraP.jl | Speedup |
|----------|:-------------------:|:------------:|:-------:|
| 1000 | 1.80 s | 0.53 s | 3.4x |
| 2000 | 2.48 s | 1.14 s | 2.2x |
| 4000 | 3.68 s | 2.08 s | 1.8x |

The performance improvement comes from reformulating the algorithm as BLAS matrix operations while preserving numerical equivalence with the original Python implementation. Both implementations use BLAS with default multithreading settings for fair comparison. \autoref{fig:benchmark} shows the performance breakdown by operation type.

![End-to-end performance comparison between BoltzTraP2 and BoltzTraP.jl. Stacked bars show time spent in interpolation and integration phases. BoltzTraP.jl achieves 1.8-3.4x speedup with efficient memory allocation scaling linearly with the number of k-points.](benchmark.png){#fig:benchmark}

# Validation

Since `BoltzTraP.jl` is a port of BoltzTraP2, validation focuses on numerical equivalence with the original implementation.

## Reference Testing

The package was developed using reference testing to guarantee identical results with BoltzTraP2:

1. Reference data was generated by running BoltzTraP2 on test systems and saving intermediate results (k-points, eigenvalues, interpolation coefficients, transport tensors).
2. Each Julia function was validated to reproduce the reference output within numerical precision.
3. All reference tests must pass before changes are merged.

This approach ensures that `BoltzTraP.jl` produces the same results as BoltzTraP2—the correctness of the algorithm is inherited from the extensively-validated original.

## Numerical Comparison

| System | Type | Bands | K-points | Max error (hartree) |
|--------|------|-------|----------|----------------|
| Si | Semiconductor | 6 | 165 | < 10⁻¹² |
| PbTe | Thermoelectric | 14 | 145 | < 10⁻⁶ |

To demonstrate that `BoltzTraP.jl` faithfully reproduces the original Python implementation, \autoref{fig:validation} compares transport coefficients computed by both codes for silicon at 300 K. The results are visually indistinguishable, confirming numerical equivalence across the entire chemical potential range.

![End-to-end validation of transport coefficients for silicon at 300 K, using VASP [@kresse1996efficient] band structure data distributed with BoltzTraP2 [@madsen2018boltztrap2]. Three panels show Seebeck coefficient (S), electrical conductivity (σ/τ), and thermal conductivity (κ/τ) as functions of chemical potential. BoltzTraP2 (red circles) and BoltzTraP.jl (blue markers) produce numerically equivalent results.](transport_Si_300K.png){#fig:validation}

The figure exhibits the expected semiconductor physics: the Seebeck coefficient is positive for p-type carriers (μ in the valence band region) and negative for n-type carriers (μ in the conduction band region), with sign reversal occurring within the band gap. The electrical and thermal conductivities vanish within the gap and increase as the chemical potential enters the bands.

The test suite includes:
- Reference tests verifying numerical equivalence with BoltzTraP2
- Over 1000 unit tests covering all package functionality
- Continuous integration on Linux, macOS, and Windows

# Documentation and Installation

Full documentation is available at: https://hsugawa8651.github.io/BoltzTraP.jl/

Installation via Julia's package manager:

```julia
using Pkg
Pkg.add("BoltzTraP")
```

# Community Guidelines

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/hsugawa8651/BoltzTraP.jl/issues).

# Acknowledgements

This package is a Julia port of BoltzTraP2, and we gratefully acknowledge the original authors: Georg K. H. Madsen, Jesus Carrete, and Matthieu J. Verstraete. We also thank the Julia community and the developers of DFTK.jl, Wannier.jl, Spglib.jl, FFTW.jl. `BoltzTraP.jl` is released under the GPL-3.0-or-later license, the same as BoltzTraP2.

# References
