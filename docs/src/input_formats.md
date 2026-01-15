# [Input Formats](@id input-formats)

BoltzTraP.jl supports multiple DFT output formats. Each format has a dedicated loader function that returns a [`DFTData`](@ref) struct.

| Format | Input Files | Loader |
|--------|-------------|--------|
| [VASP](@ref input-vasp) | `vasprun.xml` | [`load_vasp`](@ref) |
| [Quantum ESPRESSO](@ref input-qe) | `data-file-schema.xml` | [`load_qe`](@ref) |
| [Wien2k](@ref input-wien2k) | `case.struct`, `case.energy` | [`load_wien2k`](@ref) |
| [GENE/Generic](@ref input-gene) | `case.structure`, `case.energy` | [`load_gene`](@ref) |
| [ABINIT](@ref input-abinit) | `*_GSR.nc` (NetCDF) | [`load_abinit`](@ref) |
| [DFTK.jl](@ref input-dftk) | `scfres` from DFTK.jl | [`load_dftk`](@ref) |

---

## Format Detection

```@docs
detected_format
```

```@docs
load_dft
```

---

## [VASP](@id input-vasp)

```@docs
load_vasp
```

**Required files:**

| File | Required | Description |
|------|----------|-------------|
| `vasprun.xml` | Yes | Band energies, k-points, lattice, Fermi energy |
| `POSCAR` | No | Crystal structure (also in vasprun.xml) |
| `EIGENVAL` | No | Alternative band energy source |

**Directory structure:**
```
Si.vasp/
├── vasprun.xml    # Main data source
├── POSCAR         # Optional
└── EIGENVAL       # Optional
```

---

## [Quantum ESPRESSO](@id input-qe)

```@docs
load_qe
```

**Required files:**

| File | Required | Description |
|------|----------|-------------|
| `*.xml` or `data-file-schema.xml` | Yes | Main XML from `pw.x` |
| `bands.dat.gnu` | Optional | Band data from `bands.x` |

**Supported calculations**:
- Non-magnetic (nspin=1)
- Collinear spin-polarized (nspin=2, `lsda=true`)

**Not supported**: Non-collinear / spin-orbit calculations

**Directory structure:**
```
Si.ESPRESSO/out/
└── silicon.xml    # or data-file-schema.xml
```

---

## [Wien2k](@id input-wien2k)

```@docs
load_wien2k
```

**Required files:**

| File | Required | Description |
|------|----------|-------------|
| `case.struct` | Yes | Crystal structure |
| `case.energy` or `case.energyso` | Yes | Band energies |
| `case.scf` | Yes | Fermi level |
| `case.mommat2` | No | Momentum matrix elements |

**Energy file variants:**

| Files | Spin | Notes |
|-------|------|-------|
| `case.energy` | Non-spin-polarized | dosweight = 2.0 |
| `case.energyso` | Spin-orbit coupling | dosweight = 1.0, NSpin = 1 |
| `case.energyup` + `case.energydn` | Spin-polarized | dosweight = 1.0, NSpin = 2 |
| `case.energysoup` + `case.energysodn` | SOC + spin | dosweight = 1.0, NSpin = 2 |

**Supported lattice types:**

| Code | Lattice Type |
|------|--------------|
| `P` | Primitive |
| `H` | Hexagonal |
| `R` | Rhombohedral |
| `F` | Face-centered |
| `B` | Body-centered (I) |
| `CXY`, `CXZ`, `CYZ` | Base-centered (C, B, A) |

**Directory structure:**
```
Si/
├── Si.struct      # Crystal structure
├── Si.scf         # SCF output (Fermi level)
├── Si.energy      # Band energies
└── Si.mommat2     # Optional momentum matrix
```

---

## [GENE/Generic](@id input-gene)

```@docs
load_gene
```

**Required files:**

| File | Required | Description |
|------|----------|-------------|
| `case.structure` | Yes | Crystal structure (lattice + atomic positions) |
| `case.energy` | Yes | Band energies and Fermi level |

**File formats:**

- `.structure` format:
  - Line 1: Title (ignored)
  - Lines 2-4: Lattice vectors (3x3, Bohr, row-major)
  - Line 5: Number of atoms
  - Lines 6+: Element symbol + Cartesian coordinates (Bohr)

- `.energy` format:
  - Line 1: Title (ignored)
  - Line 2: `nk nspin efermi(Ry)`
  - For each spin channel, for each k-point:
    - Header: `kx ky kz nband`
    - Band lines: `energy [vx vy vz]` (Rydberg units)

**Spin handling:**
- For spin-polarized calculations (`nspin=2`), spin channels are concatenated
- Stored as `DFTData{1}` with doubled band count and `dosweight=1.0`

**Directory structure:**
```
Si.GENE/
├── Si.structure   # Crystal structure
└── Si.energy      # Band energies
```

---

## [ABINIT](@id input-abinit)

```@docs
load_abinit
```

**Required files:**

| File | Required | Description |
|------|----------|-------------|
| `*_GSR.nc` | Yes | Ground State Results (NetCDF format) |

**Installation:**
```julia
using Pkg
Pkg.add("NCDatasets")
```

**NetCDF variables used:**

| Variable | Shape | Description |
|----------|-------|-------------|
| `primitive_vectors` | (3, 3) | Lattice vectors (Bohr) |
| `reduced_atom_positions` | (3, natom) | Fractional coordinates |
| `atom_species` | (natom,) | Species index (1-based) |
| `atom_species_names` | (80, ntypat) | Element names |
| `reduced_coordinates_of_kpoints` | (3, nkpt) | K-points (fractional) |
| `eigenvalues` | (nband, nkpt, nspin) | Band energies (Hartree) |
| `fermie` | scalar | Fermi energy (Hartree) |
| `nelect` | scalar | Number of electrons |

**Directory structure:**
```
Si.abinit/
└── outsi_DS1_GSR.nc   # NetCDF Ground State Results
```

---

## [DFTK.jl](@id input-dftk)

See [DFTK.jl Integration](@ref) for detailed documentation.

---

## Test Data

BoltzTraP.jl includes synthetic test data for all supported formats in [`test/synthetic_data/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/test/synthetic_data):

| Structure | Format |
|-----------|--------|
| Monoclinic | ABINIT, GENE, QE, VASP, Wien2k |
| Triclinic | ABINIT, GENE, QE, VASP, Wien2k |

These files are generated by [`test/generate_synthetic_data.jl`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/test/generate_synthetic_data.jl) (pure Julia, no external dependencies). For the Python-based generator, see [Reference Tests](@ref) page.

!!! note "Validation Reference Data"
    For real material data (Si, PbTe, Li, etc.) used in validation against Python BoltzTraP2, see [Validation](@ref) page.

---

## Next Steps

- [Output Formats](@ref output-formats) - Result file formats and `describe` command
- [API Workflow](@ref) - Using loaders in Julia code
- [CLI Workflow](@ref) - Command-line interface
