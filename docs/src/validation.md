# Validation

This page explains how BoltzTraP.jl was validated against Python BoltzTraP2, and how you can verify the results yourself.

## Overview

BoltzTraP.jl is a faithful Julia port of Python BoltzTraP2. All transport coefficients (electrical conductivity σ, Seebeck coefficient S, thermal conductivity κ) match within numerical precision (< 10⁻⁶ relative error).

## How BoltzTraP.jl Was Validated

### Reference Testing Methodology

BoltzTraP.jl was developed using **reference testing** - a methodology where:

1. Python BoltzTraP2 computes reference outputs for various inputs
2. Julia BoltzTraP.jl computes the same quantities from the same inputs
3. Results are compared to verify numerical equivalence

This approach guarantees that BoltzTraP.jl produces results identical to the original Python implementation.

### Validation Summary

| Aspect | Details |
|--------|---------|
| Test cases | 727 (5004 assertions) |
| Tolerance | < 10⁻⁶ relative error |
| Data source | BoltzTraP2 v25.11.1 (PyPI) |
| Materials tested | Si, PbTe, Fe (collinear) |

### Reference Data Source

All reference data originates from the test materials included in the [BoltzTraP2 package](https://gitlab.com/souza-group/BoltzTraP2) (v25.11.1). These DFT calculations were performed by the BoltzTraP2 developers and are distributed with the package for testing purposes.

**Repository and Citation:**

- **GitLab**: [https://gitlab.com/souza-group/BoltzTraP2](https://gitlab.com/souza-group/BoltzTraP2)
- **PyPI**: [https://pypi.org/project/BoltzTraP2/](https://pypi.org/project/BoltzTraP2/)
- **Paper**: Madsen, G. K., Carrete, J., & Verstraete, M. J. (2018). BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients. *Computer Physics Communications*, 231, 140-145. [doi:10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

**DFT Calculation Details:**

| Material | Structure | DFT Code | XC Functional | k-grid | Pseudopotential |
|----------|-----------|----------|---------------|--------|-----------------|
| Si | Diamond (Fd-3m) | VASP | PBE | 17×17×17 | PAW_PBE Si |
| PbTe | Rock salt (Fm-3m) | VASP | PBE | 16×16×16 | PAW_PBE Pb_d, Te |
| Fe | BCC (Im-3m) | QE | PBE | - | ONCV |

**Data Provenance and Verification:**

For package checksums, test data checksums, and verification instructions, see [`reftest/README.md`](https://github.com/hsugawa8651/BoltzTraP.jl/blob/main/reftest/README.md#data-source).

### Components Tested

| Component | Test Cases | Description |
|-----------|------------|-------------|
| Symmetry | 37 | Space group detection, rotations, equivalences |
| Types | 46 | DFTData, SpinTypes, BandData |
| Interpolation | 57 | Fourier coefficient fitting, reconstruction |
| Transport | 31 | Fermi integrals, Onsager coefficients |
| I/O | 130 | VASP/QE file parsing |
| Loaders | 92 | Multi-format loader comparison |
| Workflow | 63 | End-to-end magnetic validation |

### Reference Data Scripts

All reference data is generated from Python BoltzTraP2 using scripts in the [`reftest/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/reftest) directory.
See [`reftest/README.md`](https://github.com/hsugawa8651/BoltzTraP.jl/blob/main/reftest/README.md) for complete details.

**Directory structure:**

```
reftest/
├── README.md                         # Complete usage instructions
├── common.py                         # Shared utilities for generators
├── compare.jl                        # Julia comparison script
├── convert_dft_to_gene.py            # DFT → GENE format converter
├── generate_0_unit_test_data.py      # [0] Unit test fixtures (test/data/)
├── generate_1_sphere.py              # [1] Sphere, equivalences, symmetry
├── generate_2_interpolation_si.py    # [2] Interpolation (Si, cubic)
├── generate_3_transport.py           # [3] Transport coefficients
├── generate_4_io_poscar.py           # [4] POSCAR I/O
├── generate_4_io_bt2.py              # [4] BT2 file format
├── generate_5_e2e.py                 # [5] End-to-end (Si, Li, PbTe)
├── generate_6_loaders.py             # [6] DFT loader validation
├── generate_synthetic_lowsym.py      # [7] Synthetic low-symmetry test data
├── generate_synthetic_p1.py          # [7] Synthetic P1 triclinic data
└── data/                             # Generated .npz files (gitignored)
    ├── si_interpolation.npz
    ├── si_end2end.npz
    ├── pbte_end2end.npz
    └── ...
```

**Key scripts:**

| Category | Script | Description |
|----------|--------|-------------|
| 0 | `generate_0_unit_test_data.py` | Unit test fixtures |
| 1 | `generate_1_sphere.py` | Symmetry and equivalences |
| 2 | `generate_2_interpolation_*.py` | Interpolation coefficients |
| 3 | `generate_3_transport.py` | Transport coefficients |
| 4 | `generate_4_io_*.py` | File I/O (POSCAR, BT2) |
| 5 | `generate_5_e2e.py` | End-to-end workflows |
| 6 | `generate_6_loaders.py` | DFT loader validation (VASP, QE, Wien2k, GENE, ABINIT) |
| 7 | `generate_synthetic_*.py` | Synthetic low-symmetry test data |
| - | `compare.jl` | Julia comparison against reference data |

### CoSb₃ Fourier documented baseline (manual)

The Fourier path was exercised on CoSb₃ (conventional cubic skutterudite,
space group 204, a = 17.080296 Bohr ≈ 9.04 Å; lattice data shipped at
`reftest/data/geometries/cosb3.jl`) using a DFTK SCF result followed by a
non-self-consistent band evaluation (`compute_bands`) on a 16×16×16 mesh,
then BoltzTraP.jl's Fourier (star-function) interpolation and transport.
The baseline records σ_xx × τ at T = 300 K (τ = 10 fs) for 17
chemical-potential offsets μ_eff in eV relative to the intrinsic Fermi
level; the full data is shipped at `reftest/data/cosb3_fourier_baseline.toml`.

This is a manual documented baseline — a frozen BoltzTraP.jl
self-consistency record, not an automated CI regression test (unlike the
Si Fourier and Si Wannier reftests, which run in CI) and not a
Pizzi-conformance test. It exists so future lattice or interpolation
regressions in the Fourier path are detectable by hand. The two anchor
points are:

| μ_eff | σ_xx × τ | reference |
|------:|---------:|:---|
| (eV) | (×10⁶ S/m) | (×10⁶ S/m) |
| **−1.25** | **0.9779** | ≈ 1.0 (Pizzi 2014 valence anchor, ~3% agreement) |
| **0.00** | **0.0394** | ≈ 0 (mid-gap, expected) |

The DFTK SCF checkpoint used to produce these numbers (≈ 295 MB) is not
shipped with the repository; the values in
`reftest/data/cosb3_fourier_baseline.toml` stand on their own as the
documented baseline.

## Verify It Yourself

You can reproduce the validation using the tools provided in the [`reftest/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/reftest) and [`validation/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/validation) directories.

## Prerequisites

### 1. Install Python BoltzTraP2

```bash
pip install boltztrap2
```

Or from source:

```bash
git clone https://gitlab.com/souza-group/BoltzTraP2.git
cd BoltzTraP2
pip install -e .
```

### 2. Prepare Test Data

BoltzTraP2 includes test materials. After installation, locate the data directory:

```bash
# Find BoltzTraP2 installation
python -c "import BoltzTraP2; print(BoltzTraP2.__path__[0])"

# Test data is in the 'data' subdirectory:
#   BoltzTraP2/data/Si.vasp/
#   BoltzTraP2/data/PbTe.vasp.unpolarized/
```

## Quick Validation

### Step 1: Generate Python Reference Data

```bash
cd BoltzTraP.jl/reftest

# Generate all reference data
python generate_1_sphere.py
python generate_2_interpolation_si.py
python generate_3_transport.py
python generate_4_io_poscar.py
python generate_4_io_bt2.py
python generate_5_e2e.py
```

This creates `.npz` files in `reftest/data/` containing Python-computed results:
- `si_end2end.npz` - Si transport coefficients
- `pbte_end2end.npz` - PbTe transport coefficients
- And others for intermediate validation

### Step 2: Run Automated Comparison

```bash
cd BoltzTraP.jl/reftest

# Run all comparisons
julia compare.jl

# Or run specific tests
julia compare.jl si_diamond      # Symmetry test
julia compare.jl transport       # Transport test
julia compare.jl integrate_e2e   # End-to-end test
```

### Step 3: Generate Julia Transport Results

Generate Julia transport results from the same input data:

```bash
cd BoltzTraP.jl

# Generate Si transport (from VASP data in reftest)
julia --project -e '
using BoltzTraP, NPZ
data = npzread("reftest/data/vasp_si.npz")
dft_data = (
    lattice = data["lattvec"],
    positions = transpose(data["positions"]),
    species = ["Si", "Si"],
    kpoints = transpose(data["kpoints"]),
    ebands = data["ebands"],
    fermi = data["fermi"],
    nelect = data["nelect"],
    dosweight = data["dosweight"],
)
interp = run_interpolate(dft_data; source="VASP", kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])
save_integrate("si_transport.jld2", transport)
'

# Similarly for PbTe
julia --project -e '
using BoltzTraP, NPZ
data = npzread("reftest/data/vasp_pbte.npz")
dft_data = (
    lattice = data["lattvec"],
    positions = transpose(data["positions"]),
    species = ["Pb", "Te"],
    kpoints = transpose(data["kpoints"]),
    ebands = data["ebands"],
    fermi = data["fermi"],
    nelect = data["nelect"],
    dosweight = data["dosweight"],
)
interp = run_interpolate(dft_data; source="VASP", kpoints=5000)
transport = run_integrate(interp; temperatures=[300.0])
save_integrate("pbte_transport.jld2", transport)
'
```

### Step 4: Visual Comparison

Compare Python and Julia results visually:

```bash
cd BoltzTraP.jl

# Compare Si transport coefficients
julia --project validation/compare_transport.jl \
    reftest/data/si_end2end.npz \
    si_transport.jld2 \
    --title Si

# Compare PbTe transport coefficients
julia --project validation/compare_transport.jl \
    reftest/data/pbte_end2end.npz \
    pbte_transport.jld2 \
    --title PbTe
```

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--title <name>` | Material name for figure title | `untitled` |
| `--xlims <min,max>` | X-axis limits in eV | `-0.5,0.5` |
| `-t, --temperature <T>` | Temperature in K | `300` |
| `--alltemp` | Generate figures for all common temperatures | - |
| `-o, --output <file>` | Output figure path | `validation/transport_<title>_<T>K.png` |

### Step 5: View Results

Comparison figures are saved in `validation/`:

| File | Description |
|------|-------------|
| `transport_Si_300K.png` | 3-panel comparison for Si at 300K (S, σ, κ) |
| `transport_PbTe_300K.png` | 3-panel comparison for PbTe at 300K (S, σ, κ) |

## Available Materials

| Material | Type | Description |
|----------|------|-------------|
| Si | Semiconductor | Silicon (diamond structure) |
| PbTe | Thermoelectric | Lead telluride |
| Fe | Ferromagnetic metal | BCC Iron (collinear, T=1000K) |

## Available Properties

| Property | Symbol | Units |
|----------|--------|-------|
| Seebeck coefficient | S | V/K |
| Electrical conductivity | σ | S/m (divided by τ) |
| Thermal conductivity | κ | W/(m·K) (divided by τ) |

## Detailed Validation

For more comprehensive validation including intermediate quantities (Fourier coefficients, DOS, Fermi integrals), see the [Reference Tests](reftest.md) page.

## Expected Results

When validation is successful, you should see:

1. **Visual overlap**: Python (dashed) and Julia (solid) curves should be indistinguishable
2. **Numerical agreement**: Relative differences < 10⁻⁶

Example output:
```
Comparing S for Si at T=300K
  Python: S_xx = 6.575e-06 V/K
  Julia:  S_xx = 6.575e-06 V/K
  Relative difference: 2.3e-10
  ✓ PASS
```

## Troubleshooting

### "Reference data not found"

Generate the reference data first:

```bash
cd reftest
python generate_5_e2e.py
```

### "BoltzTraP2 module not found"

Install Python BoltzTraP2:

```bash
pip install boltztrap2
```

### Large numerical differences

If differences exceed 10⁻⁶:

1. Check BoltzTraP2 version (v25.11.1 recommended)
2. Ensure reference data was regenerated after any BoltzTraP2 update
3. Report an issue on GitHub with details

## See Also

- [Reference Tests](reftest.md) - Technical details of reference tests (for developers extending tests)
- [Benchmarks](benchmarks.md) - Performance comparison with Python
