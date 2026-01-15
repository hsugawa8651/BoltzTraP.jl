# Reference Test (reftest)

Reference testing environment for validating BoltzTraP.jl against Python BoltzTraP2.

Reference tests ensure that BoltzTraP.jl produces numerically identical results to the original Python implementation. Each test compares Julia output against pre-computed Python reference data.

**Note:** These reference tests and data generation scripts were used by the developers during the implementation of BoltzTraP.jl. The same tools are provided for users who wish to independently verify the numerical equivalence.

## What is Tested

Reference tests validate numerical equivalence for:

| Category | What is Compared | Tolerance |
|----------|------------------|-----------|
| Symmetry | Equivalence classes, star functions, rotation matrices | exact |
| Interpolation | Phase factors, Fourier coefficients, reconstructed bands | < 10⁻¹⁰ |
| Transport | DOS, transport DOS, Fermi integrals (L₀, L₁, L₂), Onsager coefficients (σ, S, κ) | < 10⁻⁶ |
| I/O | Lattice vectors, atomic positions, k-points, band energies, Fermi level | exact |
| End-to-end | Full workflow from DFT data to transport coefficients | < 10⁻⁶ |

**Test Materials:**
- Si (diamond, semiconductor)
- PbTe (rock salt, thermoelectric)
- Synthetic monoclinic/triclinic (low-symmetry validation)

## Quick Start

```bash
cd reftest

# Step 1: Install Python BoltzTraP2 (if not already installed)
pip install boltztrap2

# Step 2: Generate Si.GENE (GENE loader test data)
python convert_dft_to_gene.py /path/to/BoltzTraP2/data/Si.vasp

# Step 3: Generate reference data (~6GB)
python generate_1_sphere.py
python generate_2_interpolation_si.py
python generate_3_transport.py
python generate_5_e2e.py
python generate_6_loaders.py gene  # Requires Si.GENE from Step 2

# Step 4: Run comparison tests
julia compare.jl
```

**Note:** Step 2 creates `Si.GENE` in the BoltzTraP2 data directory, which is required for GENE loader tests.

## Prerequisites

### Python

Install Python BoltzTraP2:

```bash
pip install boltztrap2
```

Or from source:

```bash
git clone https://gitlab.com/souza-group/BoltzTraP2.git BoltzTraP2-public
cd BoltzTraP2-public
pip install -e .
```

### Julia

Ensure NPZ.jl is available:

```julia
using Pkg
Pkg.add("NPZ")
```

## Step 1: Generate Reference Data

Reference data files (`data/*.npz`, ~6GB) are **not included** in the repository. Generate them locally using Python scripts.

### Data Directory Option

By default, generator scripts use the BoltzTraP2 data directory from the pip-installed package. You can specify a custom data directory:

```bash
# Default: use pip-installed BoltzTraP2 data directory
python generate_5_e2e.py

# Custom: use a local clone of BoltzTraP2
python generate_5_e2e.py --data-dir /path/to/BoltzTraP2/data
```

### Generator Scripts

```bash
cd reftest

# Category 1: Sphere, equivalences, symmetry
python generate_1_sphere.py

# Category 2: Interpolation
python generate_2_interpolation_si.py

# Category 3: Transport
python generate_3_transport.py

# Category 4: I/O
python generate_4_io_poscar.py
python generate_4_io_bt2.py

# Category 5: End-to-end (Si, PbTe)
python generate_5_e2e.py

# Category 6: DFT loader validation
python generate_6_loaders.py vasp      # VASP loader
python generate_6_loaders.py qe        # Quantum ESPRESSO loader
python generate_6_loaders.py wien2k    # Wien2k loader
python generate_6_loaders.py gene      # GENE loader (requires Si.GENE, see below)
python generate_6_loaders.py abinit    # ABINIT loader
python generate_6_loaders.py all       # All loaders (requires Si.GENE for GENE)

# Category 7: Synthetic low-symmetry test data
python generate_synthetic_lowsym.py    # All formats (GENE, VASP, QE, Wien2k, ABINIT)
```

### Utility: DFT to GENE Converter

Convert any DFT data to GENE format:

```bash
# Convert VASP to GENE
python convert_dft_to_gene.py /path/to/Si.vasp

# Convert QE to GENE
python convert_dft_to_gene.py /path/to/Si.ESPRESSO/out --name Si

# Convert Wien2k to GENE
python convert_dft_to_gene.py /path/to/Si --output ./Si.GENE
```

This utility was used to create `Si.GENE` from `Si.vasp` for GENE loader testing.

## Step 2: Run Comparison Tests

### Basic Usage

```bash
cd reftest
julia compare.jl
```

### Specific Tests

```bash
julia compare.jl si_diamond
julia compare.jl transport
julia compare.jl integrate_e2e
```

### Visual Comparison (optional)

Generate transport comparison figures. First, generate Julia transport results, then compare with Python results:

```bash
cd ..  # BoltzTraP.jl directory

# Generate Julia transport (see docs/src/validation.md for details)
# Then compare:
julia --project validation/compare_transport.jl \
    reftest/data/si_end2end.npz si_transport.jld2 --title Si

julia --project validation/compare_transport.jl \
    reftest/data/pbte_end2end.npz pbte_transport.jld2 --title PbTe
```

Output figures are saved in `validation/` as `transport_<title>_<T>K.png`.

## Directory Structure

```
reftest/
├── README.md
├── common.py                         # Shared utilities for generator scripts
├── compare.jl                        # Julia comparison script
├── convert_dft_to_gene.py            # DFT → GENE format converter
├── generate_0_unit_test_data.py      # [0] Unit test fixtures (test/data/)
├── generate_1_sphere.py              # [1] Sphere, equivalences, symmetry
├── generate_2_interpolation_si.py    # [2] Interpolation (Si, cubic)
├── generate_3_transport.py           # [3] Transport coefficients
├── generate_4_io_poscar.py           # [4] POSCAR I/O
├── generate_4_io_bt2.py              # [4] BT2 file format
├── generate_5_e2e.py                 # [5] End-to-end (Si, PbTe)
├── generate_6_loaders.py             # [6] DFT loader validation
├── generate_synthetic_lowsym.py      # [7] Synthetic low-symmetry test data
├── generate_synthetic_p1.py          # [7] Synthetic P1 triclinic data
└── data/                             # Generated .npz files (gitignored, ~6GB)
    ├── unit_cube.npz
    ├── simple_cubic.npz
    ├── si_diamond.npz
    ├── si_interpolation.npz
    ├── simple_transport.npz
    ├── si_end2end.npz
    ├── pbte_end2end.npz
    └── ...
```

## Data Source

### BoltzTraP2 Test Materials

Reference data originates from [BoltzTraP2](https://gitlab.com/souza-group/BoltzTraP2) test materials.

**Package Provenance:**

| Item | Value |
|------|-------|
| Package | `boltztrap2-25.11.1.tar.gz` |
| PyPI | https://pypi.org/project/BoltzTraP2/25.11.1/ |
| SHA256 | `1a5493cdc9a9d834f1f33f7efad89896e3f8a3f8e36eee0c1699f41d5666234a` |

**Test Data Checksums (VASP vasprun.xml SHA256):**

| Material | File Date | SHA256 |
|----------|-----------|--------|
| Si.vasp | 2018-01-08 | `09492bee9109f3e82deb1592458a9b2229eac2b1b0d4ad3751360fb544efedb7` |
| PbTe.vasp.unpolarized | 2018-08-08 | `d2cee34d033706f9c988442944d2c0c9bb339981c365b5c3d3f2446a990619a5` |
| Li.vasp | 2018-01-08 | `e5c4aa251aa307c79e1c82c92ee8946ffd24c0d37b91db7e5ecc708e011cfdcf` |

**Test Data Checksums (Wien2k case.struct SHA256):**

| Material | Energy File | NSpin | SHA256 |
|----------|-------------|-------|--------|
| Si | .energy | 1 | `e3108b1ff3844bc8214a88757656a08c36e5f60ec21aa75aef5b3250eb5d1e2d` |
| Li.W2K | .energy | 1 | `0b408714755dad8bb2a81302ce2b0384d8859df62c4f5dcc76ebf6e5c93021f7` |
| CoSb3 | .energy | 1 | `b52f0bf5f8b33ee8d16d0e8cd34a5bfd242b84c03bb6ee940725f741761831ea` |
| Bi2Te3 | .energyso | 2* | `615ec78b87c52c560a712758934432c9540b1f4f7c9aee1dd0d43e2546232f6e` |

*Bi2Te3 uses spin-orbit coupling (.energyso), which sets dosweight=1.0 (nspin=2 in DFTData).

**Verify package checksum:**

```bash
pip download boltztrap2==25.11.1 --no-deps
shasum -a 256 boltztrap2-25.11.1.tar.gz
```

**Verify test data checksum:**

```bash
python -c "import BoltzTraP2; print(BoltzTraP2.__path__[0])"
# Then verify:
shasum -a 256 <path>/data/Si.vasp/vasprun.xml
```

**Citation:**

> Madsen, G. K., Carrete, J., & Verstraete, M. J. (2018). BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients. *Computer Physics Communications*, 231, 140-145. [doi:10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

### Synthetic Low-Symmetry Test Data

The synthetic monoclinic and triclinic test structures were generated using:

| Tool | Version | License | Purpose |
|------|---------|---------|---------|
| [pymatgen](https://pymatgen.org/) | 2024.x | MIT | Structure creation, format conversion |
| [ASE](https://wiki.fysik.dtu.dk/ase/) | 3.x | LGPL-2.1 | Structure manipulation |
| [NumPy](https://numpy.org/) | 1.x | BSD-3 | Numerical operations |
| [netCDF4](https://unidata.github.io/netcdf4-python/) | 1.x | MIT | ABINIT NetCDF output |

These files contain **synthetic band structure data** (free-electron model) for testing DFT loader functionality with low-symmetry crystal structures.

**Generator script:** `reftest/generate_synthetic_lowsym.py`

## Notes

- Tolerance: rtol=1e-10, atol=1e-12
- Reference data must be regenerated if Python BoltzTraP2 version changes
- Reference data (~6GB) is gitignored; generate locally

---

## Appendix A: Synthetic Low-Symmetry Test Data Details

The `generate_synthetic_lowsym.py` script creates synthetic DFT output files for testing low-symmetry crystal structure handling across all loaders.

### Purpose

Original BoltzTraP2 test data contains only high-symmetry structures (space groups 166-229). This script creates:
- **Monoclinic** structures (β = 100°, α = γ = 90°)
- **Triclinic** structures (α = 85°, β = 80°, γ = 75°)

### Generated Formats

All 5 DFT formats are generated:

| Directory | Loader | Files Created |
|-----------|--------|---------------|
| `Synthetic.GENE.{monoclinic,triclinic}/` | GENE | `.structure`, `.energy` |
| `Synthetic.vasp.{monoclinic,triclinic}/` | VASP | `vasprun.xml` |
| `Synthetic.ESPRESSO.{monoclinic,triclinic}/out/` | QE | `synthetic.xml` |
| `Synthetic.W2K.{monoclinic,triclinic}/` | Wien2k | `.struct`, `.energy`, `.scf` |
| `Synthetic.abinit.{monoclinic,triclinic}/` | ABINIT | `synthetic_GSR.nc` |

### Physical Model

- **Structure**: 2 atoms (C, N) at general positions to break all symmetry
- **Band Model**: Free-electron approximation E(k) = |k + G|² / 2
- **K-points**: 4×4×4 Monkhorst-Pack grid (64 points)
- **Bands**: 8 bands per k-point

### Usage

```bash
cd reftest
python generate_synthetic_lowsym.py

# Or with custom BoltzTraP2 data directory (for writing):
python generate_synthetic_lowsym.py --data-dir /path/to/BoltzTraP2/data
```

**Note:** This script writes synthetic DFT files to the BoltzTraP2 data directory. When using the pip-installed package, ensure you have write permissions or use `--data-dir` to specify a writable location.

Generated data is used by `test/test_loaders_reftest.jl` in the "Synthetic Low-Symmetry Tests" testset.

### Adding New Test Structures

To add a new cell type:

1. Define cell parameters in `generate_synthetic_lowsym.py`:
   ```python
   ("new_type", 5.1, 4.8, 6.2, 75.0, 80.0, 85.0),  # a, b, c, α, β, γ
   ```

2. Regenerate all formats:
   ```bash
   python generate_synthetic_lowsym.py
   ```

3. Add test case in `test/test_loaders_reftest.jl`:
   ```julia
   ("new_type", (α = 75.0, β = 80.0, γ = 85.0)),
   ```

## Appendix B: Unit Test Fixture Data

The `test/data/` directory contains pre-computed reference data (~7MB) used by `Pkg.test()`. These fixtures are **included** in the repository for convenience.

### Files

| File | Description | Shape |
|------|-------------|-------|
| `Si_equivalences.npz` | Equivalence classes for Si | 5048 classes |
| `Si_fitde3D.npz` | Fourier coefficients | (6, 5048) |
| `Si_BTPdos.npz` | Density of states | 93 energy points |
| `Si_BTPdos_lambda.npz` | DOS (with λ=0) | 93 energy points |
| `Si_fermiintegrals.npz` | Fermi integrals L0, L1, L2 | (40, 100, 3, 3) |
| `Si_Onsager.npz` | Onsager coefficients | (40, 100, 3, 3) |
| `Si_cv.npz` | Heat capacity | (40, 100) |
| `kpoints.npz` | K-points from Si | (165, 3) |
| `Si_old_mommat_ref.npz` | Momentum matrix (old format) | - |
| `Si_new_mommat_ref.npz` | Momentum matrix (new format) | - |

### Generation Script

To regenerate these fixtures (e.g., for a different BoltzTraP2 version):

```bash
cd reftest

# Generate to /tmp for verification (default)
python generate_0_unit_test_data.py

# Or write directly to test/data/ (overwrites existing)
python generate_0_unit_test_data.py --write
```

**Parameters:**
- 40 temperatures (100-800 K)
- 100 chemical potentials (±0.15 Ha around Fermi level)
- 93 DOS points
- Si diamond structure from BoltzTraP2 data

**Note:** The exact output depends on the BoltzTraP2 version and DFT data. Minor numerical differences may occur.

## Appendix C: Pure-Julia Synthetic Data Generator

For users who want to verify loaders without Python, a pure-Julia generator is available.

### Location

`test/generate_synthetic_data.jl`

### Usage

```bash
# Generate synthetic test data
julia --project=. test/generate_synthetic_data.jl

# Verify all loaders
julia --project=. test/generate_synthetic_data.jl --verify
```

### Features

- Generates monoclinic and triclinic structures
- Supports all 5 formats: GENE, VASP, QE, Wien2k, ABINIT
- ABINIT requires NCDatasets.jl (optional)
- Uses free-electron band model for synthetic band energies

### Output

Files are generated in `test/synthetic_data/`:

```
test/synthetic_data/
├── Synthetic.GENE.monoclinic/
├── Synthetic.GENE.triclinic/
├── Synthetic.vasp.monoclinic/
├── Synthetic.vasp.triclinic/
├── Synthetic.ESPRESSO.monoclinic/out/
├── Synthetic.ESPRESSO.triclinic/out/
├── Synthetic.W2K.monoclinic/
├── Synthetic.W2K.triclinic/
├── Synthetic.abinit.monoclinic/
└── Synthetic.abinit.triclinic/
```

### Example Output

```
$ julia --project=. test/generate_synthetic_data.jl --verify
Verifying loaders...
============================================================

[monoclinic] Expected: α=90.0°, β=100.0°, γ=90.0°
------------------------------------------------------------
  GENE    : PASS (α=90.0°, β=100.0°, γ=90.0°)
  VASP    : PASS (α=90.0°, β=100.0°, γ=90.0°)
  QE      : PASS (α=90.0°, β=100.0°, γ=90.0°)
  Wien2k  : PASS (α=90.0°, β=100.0°, γ=90.0°)
  ABINIT  : PASS (α=90.0°, β=100.0°, γ=90.0°)

[triclinic] Expected: α=85.0°, β=80.0°, γ=75.0°
------------------------------------------------------------
  GENE    : PASS (α=85.0°, β=80.0°, γ=75.0°)
  ...
============================================================
All tests PASSED
```
