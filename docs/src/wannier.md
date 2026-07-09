# Wannier path

BoltzTraP.jl supports band and velocity interpolation through a Wannier
basis via the [Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl)
package. The Wannier path is loaded as a weak-dependency extension
(`BoltzTraPWannierExt`) and becomes available once `using Wannier` is
in scope alongside `using BoltzTraP`.

The extension exposes a [`WannierInterpolator`](@ref) type that
implements the same [`AbstractInterpolator`](@ref) interface as the
default Fourier (Shankland) path. Once constructed, it goes through the
same entry point, [`run_integrate`](@ref)`(interp, sys; sampling)`.

```julia
using BoltzTraP
using Wannier

wi = WannierInterpolator("path/to/silicon")          # prefix entry point
sys = TransportSystem(wi; fermi = 0.2, nelect = 8.0)
result = run_integrate(wi, sys; sampling = UniformMesh((8, 8, 8)), temperatures = [300.0])
```

Only the interpolator changes between the two paths. What the Wannier
path additionally needs — the Fermi energy, the electron count, and an
explicit mesh — is visible in the call: the Wannier90 files carry the
Hamiltonian, not the state of the electron system.

## Construction

`WannierInterpolator` provides two construction entry points. Both
produce a value of the same type and feed the same downstream
interpolation methods.

### From a Wannier90 file prefix

```julia
wi = WannierInterpolator(prefix::String; spintype = Unpolarized())
```

The prefix points to a Wannier90 run on disk: `<prefix>.win` plus the
unitary matrices stored in `<prefix>.chk` (or `<prefix>.chk.fmt`). The
constructor calls `Wannier.read_w90_interp` internally and then
delegates to the in-memory entry point below, centralizing the
Ångström → Bohr lattice conversion.

This is the usual entry point when the Wannierization step has already
been performed by an external Wannier90 binary (typical for VASP, QE,
Wien2k, or ABINIT pipelines).

### From an in-memory `Wannier.InterpModel`

```julia
wi = WannierInterpolator(imodel::Wannier.InterpModel; spintype = Unpolarized())
```

This entry point accepts an `InterpModel` produced inside the running
Julia session, without first writing a `.chk` binary to disk. Use it
when the Wannierization is performed in-process — for example a
`Wannier.read_w90` followed by `disentangle` and `Wannier.InterpModel`,
or a DFTK → Wannier.jl pipeline that builds an `InterpModel`
directly from an SCF result.

The prefix-string constructor is itself a thin wrapper around this
method.

## Running a transport calculation

### Building the system

`TransportSystem` collects the material constants that the transport
integrals need. The lattice and the spin type come from the
interpolator; the Fermi energy and the electron count do not exist in
the Wannier90 files and have to be supplied:

```julia
sys = TransportSystem(wi; fermi = 0.2, nelect = 8.0)   # fermi in Hartree
```

`dosweight` defaults to the spin degeneracy implied by the spin type,
`2.0` when unpolarized and `1.0` otherwise.

### Basic usage

The Wannier path goes through the same entry point as the Fourier path:

```julia
result = run_integrate(wi, sys; sampling = UniformMesh((8, 8, 8)),
                       temperatures = [300.0])
```

[`solve_transport`](@ref) remains available for callers that want to
address the three dispatch axes directly.

### An explicit mesh size is required

The Fourier path derives its integration grid from the star functions
that accompany the interpolation coefficients. The Wannier path has no
such structure — the `Wannier.InterpModel` carries only the real-space
Hamiltonian and the associated R-vectors — so the caller must supply an
explicit k-mesh.

In practice this means `UniformMesh((nk1, nk2, nk3))` with concrete
integers; a `UniformMesh` constructed with `nk === nothing` is rejected
when paired with a `WannierInterpolator`.

```julia
result = run_integrate(wi, sys; sampling = UniformMesh((8, 8, 8)))  # OK
run_integrate(wi, sys)                                              # errors
```

The grid that was used is recorded in `result.metadata["sampling_nk"]`.

Passing an explicit size on the Fourier path has the opposite outcome:
the size is ignored, with a warning, because the interpolation basis
already fixes the grid.

## Unit conventions

Internally BoltzTraP.jl works in atomic units (Hartree for energy,
Bohr for length), to match the rest of the code base. Wannier.jl, in
turn, follows the Wannier90 convention of eV for energies and Ångström
for lengths.

The extension performs the conversions inside its
[`interpolate_bands`](@ref) and [`interpolate_velocities`](@ref)
methods, so user-facing arrays returned by these functions are already
in BoltzTraP's Hartree/Bohr units:

| Quantity              | Wannier.jl native | After extension       |
|-----------------------|-------------------|------------------------|
| Band energies         | eV                | Hartree                |
| Group velocities      | eV·Å              | Hartree·Bohr           |
| Lattice (constructor) | Å                 | Bohr (stored in `wi.lattice`) |

The numerical conversion factors (`HA_TO_EV ≈ 27.2114`,
`BOHR_TO_ANG ≈ 0.5292`) are the same constants used by the Fourier
path, so the Wannier and Fourier interpolators round-trip with the
same definition of "atomic units" as the rest of the package.

The `k`-point convention also differs by transposition: BoltzTraP.jl
takes `kpoints` as an `(nk, 3)` matrix, whereas Wannier.jl expects
`(3, nk)`. The extension transposes internally so callers continue to
use the BoltzTraP shape.

## Available methods

```@docs
WannierInterpolator
```

The methods implemented on top of [`WannierInterpolator`](@ref) — band
interpolation, velocity interpolation, and the
[`solve_transport`](@ref) entry point — share their docstrings with
the generic [`AbstractInterpolator`](@ref) interface, so the same API
documentation applies regardless of which interpolator is in use.

## Wannier path validation

The Wannier path is exercised through two regression baselines and one
external secondary reference.

### Silicon regression baseline (CI)

`test/test_wannier_reftest.jl` re-runs a pipeline from
`WannierInterpolator` to `solve_transport` on the silicon fixture
shipped with [Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl)
(`pkgdir(Wannier)/test/fixtures/silicon/silicon`, 8 Wannier functions,
12 bands, 4×4×4 Monkhorst-Pack grid) and compares the transport
tensors to a pre-recorded JLD2 file
(`reftest/data/si_wannier_transport.jld2`). The full σ, S and κ tensors
(all entries `i, j = 1, 2, 3`) are compared at `rtol = 1e-8` across
`T = 300, 500, 700` K. The band velocities at degenerate k-points are
constructed gauge-invariantly (averaged over approach directions), so the
transport tensors no longer depend on the arbitrary eigenbasis a given
LAPACK build selects within a degenerate multiplet and the recomputation
reproduces the baseline across platforms to near machine precision.

This test runs in CI when `TEST_WANNIER=true` is set for the extension
test group.

### CoSb₃ documented baseline (manual)

The Wannier path was also exercised on CoSb₃ (conventional cubic
skutterudite, space group 204, a = 17.080296 Bohr ≈ 9.04 Å; lattice
data shipped at `reftest/data/geometries/cosb3.jl`) using a DFTK SCF
result followed by external Wannierization with `wannier90.x`
([Pizzi et al. 2020](https://doi.org/10.1088/1361-648X/ab51ff)) and the
default random-projection setup. The baseline records σ_xx × τ at
T = 300 K (τ = 10 fs) for 17 chemical-potential offsets μ_eff in eV
relative to the intrinsic Fermi level; the full data is shipped at
`reftest/data/cosb3_wannier_baseline.toml`. The two anchor points
below are compared against the
[Pizzi et al. (2014)](https://doi.org/10.1016/j.cpc.2013.09.015)
BoltzWann reference.

The two anchor points are:

| μ_eff | σ_xx × τ | reference |
|------:|---------:|:---|
| (eV) | (×10⁶ S/m) | (×10⁶ S/m) |
| **−1.25** | **1.0481** | ≈ 1.0 (Pizzi 2014 valence anchor, ~5% agreement) |
| **0.00** | **0.0208** | ≈ 0 (mid-gap, expected) |

Conduction-side σ_xx deviates from the Pizzi reference by roughly a
factor of two due to the default `HydrogenicWannierProjection` setup
(random initial centers, generic s-like Gaussians), which under-
represents the Co d / Sb p character of the conduction manifold. The
remaining 15 baseline values document the current state of the Wannier
path output, not a regression target.

The frozen Wannier model files used to produce these numbers
(`cosb3.win`, `cosb3.chk`, `cosb3.eig`, `cosb3.mmn`, total ≈ 209 MB)
come from a fixed random-projection realization that is not reproduced
as part of the test suite. They are therefore not shipped with the
repository; the values in `reftest/data/cosb3_wannier_baseline.toml`
stand on their own as the documented baseline.

### External QE + BoltzWann cross-check

An independent reference was produced through Quantum ESPRESSO +
`pw2wannier90` + `postw90.x boltzwann`
([Pizzi et al. 2020](https://doi.org/10.1088/1361-648X/ab51ff)) on the
same CoSb₃ structure, with Co d and Sb p initial projections chosen
to match the relevant orbital character. That pipeline yields an
isotropic σ(μ) curve for the Wannier path, computed using the
BoltzWann methodology
([Pizzi et al. 2014](https://doi.org/10.1016/j.cpc.2013.09.015)), and
is recorded as an additional cross-check, not as a primary regression
target. Quantitative agreement with Pizzi et al. (2014) is limited to
the valence plateau under the projection choices available today; the
external reference is included for the reader who wants an alternative
Wannierization perspective on the same material.
