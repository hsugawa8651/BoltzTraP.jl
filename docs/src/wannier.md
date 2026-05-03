# Wannier.jl Integration

BoltzTraP.jl supports band and velocity interpolation through a Wannier
basis via the [Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl)
package. The Wannier path is loaded as a weak-dependency extension
(`BoltzTraPWannierExt`) and becomes available once `using Wannier` is
in scope alongside `using BoltzTraP`.

The extension exposes a [`WannierInterpolator`](@ref) type that
implements the same [`AbstractInterpolator`](@ref) interface as the
default Fourier (Shankland) path. Once constructed, it slots into the
existing transport pipeline through
[`solve_transport`](@ref)`(::UniformMesh, ::WannierInterpolator, sys)`.

```julia
using BoltzTraP
using Wannier

wi = WannierInterpolator("path/to/silicon")          # prefix entry point
result = solve_transport(UniformMesh((8, 8, 8)), wi, sys; temperatures = [300.0])
```

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

## Use with `solve_transport`

The Wannier path participates in the standard transport pipeline:

```julia
result = solve_transport(mesh::UniformMesh, wi::WannierInterpolator, sys::TransportSystem;
                         temperatures = ..., ...)
```

### `UniformMesh.nk` is required

The Fourier path can derive a natural integration mesh from the
symmetry-equivalent points that accompany the interpolation
coefficients. The Wannier path has no such structure — the
`Wannier.InterpModel` carries only the real-space Hamiltonian and the
associated R-vectors — so the caller must supply an explicit k-mesh
through `UniformMesh.nk`.

In practice this means `UniformMesh((nk1, nk2, nk3))` with concrete
integers; a `UniformMesh` constructed with `nk === nothing` is rejected
when paired with a `WannierInterpolator`.

```julia
mesh = UniformMesh((8, 8, 8))                            # OK
result = solve_transport(mesh, wi, sys; temperatures = [300.0])

mesh_nothing = UniformMesh(nothing)                      # Fourier-only
solve_transport(mesh_nothing, wi, sys; ...)              # errors
```

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
