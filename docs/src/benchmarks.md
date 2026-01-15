# Performance Benchmarks

BoltzTraP.jl includes optimized implementations for computationally intensive functions.
This page documents the performance improvements and shows how to reproduce the benchmarks.

## End-to-End Performance (vs Python BoltzTraP2)

BoltzTraP.jl achieves **1.8-3.4x speedup** over Python BoltzTraP2 for end-to-end transport calculations.

**Hardware:** MacBook Pro (Apple M2)

| k-points | BoltzTraP2 (Python) | BoltzTraP.jl | Speedup |
|----------|---------------------|--------------|---------|
| 1000 | 1.80 s | 0.53 s | **3.4x** |
| 2000 | 2.48 s | 1.14 s | **2.2x** |
| 4000 | 3.68 s | 2.08 s | **1.8x** |

The speedup is larger for smaller problems due to Julia's efficient compilation.
For larger problems, both implementations become memory-bound with similar scaling.

## Running Benchmarks

Benchmark scripts are available in the [`benchmarks/`](https://github.com/hsugawa8651/BoltzTraP.jl/tree/main/benchmarks) directory:

```bash
# End-to-end scaling benchmarks (Python vs Julia)
python benchmarks/scaling_benchmark.py   # Run Python BoltzTraP2
julia --project benchmarks/scaling_benchmark.jl  # Run Julia BoltzTraP.jl

# Generate comparison figure
julia --project benchmarks/generate_figure.jl
```

See [`benchmarks/README.md`](https://github.com/hsugawa8651/BoltzTraP.jl/blob/main/benchmarks/README.md) for details.
