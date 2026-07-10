# Architecture

This page describes BoltzTraP.jl's transport API design, in particular
the dispatch axes that allow swapping the BZ sampling and band
interpolation strategies independently.

## Method-dispatched transport API

The primary entry point for transport coefficient computation is
`solve_transport`, which dispatches on three independent axes:
`sampling::AbstractBZSampling`, `interp::AbstractInterpolator`, and
`sys::TransportSystem`.

```@docs
solve_transport
```

| Axis | Abstract type | Concrete subtype(s) | Role |
|---|---|---|---|
| `sampling` | `AbstractBZSampling` | `UniformMesh` | How the BZ is sampled |
| `interp` | `AbstractInterpolator` | `FourierInterpolator`, `WannierInterpolator` | How band energies and velocities are evaluated |
| `sys` | `TransportSystem` | (concrete struct) | Material data carrier |

A new sampling or interpolator strategy can be added by introducing a
new concrete subtype and a new `solve_transport` method dispatched on
that combination — no changes to existing types or to the entry point.

## BZ sampling

```@docs
AbstractBZSampling
UniformMesh
```

## Band interpolation

`WannierInterpolator` is documented on the [Wannier path](@ref) page, since it
only becomes available once `Wannier.jl` is loaded.

```@docs
AbstractInterpolator
FourierInterpolator
interpolate_velocities
```

## Material data carrier

`TransportSystem` collects material-dependent data that was previously
scattered across `InterpolationResult.metadata` and `BandData{ST}`:

```@docs
TransportSystem
```

Two convenience constructors assemble it from what a given interpolation
path already has. `TransportSystem(::InterpolationResult)` reads the
fields out of the interpolation metadata. `TransportSystem(::WannierInterpolator; fermi, nelect)`
takes the lattice and the spin type from the interpolator and requires
the rest, because the Wannier90 files do not carry it.

## Entry API

`run_integrate` is the entry point for both interpolation paths:

```julia
run_integrate(interp::AbstractInterpolator, sys::TransportSystem; sampling = UniformMesh(), kwargs...)
run_integrate(interp::AbstractInterpolator, sys::TransportSystem, mur::AbstractVector; kwargs...)
```

The call site is the same whichever interpolator is passed; only
`typeof(interp)` selects the backend, and `result.metadata["interpolator"]`
records which one ran. `solve_transport` stays public for callers that
want to address the three dispatch axes directly.

The Fourier path keeps its convenience forms

```julia
run_integrate(interp::InterpolationResult; kwargs...)
run_integrate(interp::InterpolationResult, mur::AbstractVector; kwargs...)
run_integrate(input::String; kwargs...)
```

which build the `TransportSystem` and the `FourierInterpolator` from the
`InterpolationResult` themselves. They restore the original
`metadata["source"]` value before the result is written out, which is why
they do not route through the two-argument form. Output is bit-equal to
the pre-refactor implementation under default arguments.

## Provenance metadata

`TransportResult.metadata` records which sampling and interpolator types
produced the output:

| Key | Example value | Source |
|---|---|---|
| `"sampling"` | `"UniformMesh"` | type name of the sampling argument |
| `"sampling_nk"` | `(23, 23, 23)` | the k-grid the integration actually used |
| `"interpolator"` | `"Fourier"` | type name of the interpolator argument with the `Interpolator` suffix stripped |

These keys are stamped by `solve_transport` and survive JLD2 round-trip
through `save_integrate` / `load_integrate`.

`"sampling_nk"` is resolved per interpolator: the Fourier path reports the
grid its interpolation basis fixes, the Wannier path reports the grid it
was given.
