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
| `interp` | `AbstractInterpolator` | `FourierInterpolator` | How band energies and velocities are evaluated |
| `sys` | `TransportSystem` | (concrete struct) | Material data carrier |

A new sampling or interpolator strategy can be added by introducing a
new concrete subtype and a new `solve_transport` method dispatched on
that combination — no changes to existing types or to the entry point.

## BZ sampling

```@docs
AbstractBZSampling
UniformMesh
```

## Material data carrier

`TransportSystem` collects material-dependent data that was previously
scattered across `InterpolationResult.metadata` and `BandData{ST}`:

```@docs
TransportSystem
```

The `TransportSystem(::InterpolationResult)` convenience constructor
extracts these fields from existing interpolation output, so callers
that already have an `InterpolationResult` need not assemble the struct
manually.

## Backward compatibility

The legacy entry points

```julia
run_integrate(interp::InterpolationResult; kwargs...)
run_integrate(interp::InterpolationResult, mur::AbstractVector; kwargs...)
```

are preserved as thin wrappers. They build a `TransportSystem` and
`FourierInterpolator` from the `InterpolationResult` and delegate to
`solve_transport(UniformMesh(), ...)`. The wrapper restores the original
`metadata["source"]` value after the delegate returns. Output is
bit-equal to the pre-refactor implementation under default arguments.

## Provenance metadata

`TransportResult.metadata` records which sampling and interpolator types
produced the output:

| Key | Example value | Source |
|---|---|---|
| `"sampling"` | `"UniformMesh"` | type name of the sampling argument |
| `"interpolator"` | `"Fourier"` | type name of the interpolator argument with the `Interpolator` suffix stripped |

These keys are stamped by `solve_transport` and survive JLD2 round-trip
through `save_integrate` / `load_integrate`.
