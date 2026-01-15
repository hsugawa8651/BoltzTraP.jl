#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#
# End-to-end benchmark: measure interpolate + integrate performance.
# Compares Python BoltzTraP2 with Julia BoltzTraP.jl.
#
# Usage:
#   python benchmarks/scaling_benchmark.py
#
# Output:
#   benchmarks/scaling_results_python.json

import json
import time
import tracemalloc
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    from BoltzTraP2 import dft as BTP
    from BoltzTraP2 import sphere, fite
    from BoltzTraP2 import bandlib as BL
    from BoltzTraP2 import units
    HAVE_BOLTZTRAP2 = True
except ImportError:
    print("ERROR: BoltzTraP2 not installed. Install with: pip install boltztrap2")
    HAVE_BOLTZTRAP2 = False

# Check if pyFFTW is available
try:
    import pyfftw
    HAVE_PYFFTW = True
except ImportError:
    HAVE_PYFFTW = False


def run_single_benchmark(data, nkpt, temperature=300.0, npts_dos=500):
    """Run a single interpolate + integrate cycle and return elapsed time.

    Args:
        data: BoltzTraP2 DFTData object (already loaded)
        nkpt: Target number of k-points for equivalences
        temperature: Temperature in K
        npts_dos: Number of DOS bins

    Returns:
        Tuple of (total_time_ms, interpolate_time_ms, integrate_time_ms)
    """
    # Pre-compute lattice info (not part of timing)
    lattvec = data.get_lattvec()
    vuc = data.get_volume()

    # Interpolate (includes get_equivalences + fitde3D)
    t_interp_start = time.perf_counter()
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt)
    coeffs = fite.fitde3D(data, equivalences)
    t_interp_end = time.perf_counter()
    interpolate_ms = (t_interp_end - t_interp_start) * 1000

    # Integrate
    t_integ_start = time.perf_counter()
    eband, vvband, _ = fite.getBTPbands(equivalences, coeffs, lattvec)
    epsilon, dos, vvdos, _ = BL.BTPDOS(eband, vvband, npts=npts_dos)

    # Set up chemical potential range
    Tr = np.array([temperature])
    margin = 9.0 * units.BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(
        epsilon > epsilon.min() + margin,
        epsilon < epsilon.max() - margin
    )
    mur = epsilon[mur_indices]

    # Fermi integrals and Onsager coefficients
    N, L0, L1, L2, _ = BL.fermiintegrals(
        epsilon, dos, vvdos, mur=mur, Tr=Tr, dosweight=data.dosweight
    )
    sigma, seebeck, kappa, _ = BL.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, vuc
    )
    t_integ_end = time.perf_counter()
    integrate_ms = (t_integ_end - t_integ_start) * 1000

    total_ms = interpolate_ms + integrate_ms
    return total_ms, interpolate_ms, integrate_ms


def measure_memory(data, nkpt, temperature=300.0, npts_dos=500):
    """Measure peak memory usage for a single run.

    Args:
        data: BoltzTraP2 DFTData object
        nkpt: Target number of k-points
        temperature: Temperature in K
        npts_dos: Number of DOS bins

    Returns:
        Peak memory in MiB
    """
    # Pre-compute lattice info (not part of measurement)
    lattvec = data.get_lattvec()
    vuc = data.get_volume()

    tracemalloc.start()

    # Run the full calculation (includes get_equivalences)
    equivalences = sphere.get_equivalences(data.atoms, data.magmom, nkpt)
    coeffs = fite.fitde3D(data, equivalences)
    eband, vvband, _ = fite.getBTPbands(equivalences, coeffs, lattvec)
    epsilon, dos, vvdos, _ = BL.BTPDOS(eband, vvband, npts=npts_dos)

    Tr = np.array([temperature])
    margin = 9.0 * units.BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(
        epsilon > epsilon.min() + margin,
        epsilon < epsilon.max() - margin
    )
    mur = epsilon[mur_indices]

    N, L0, L1, L2, _ = BL.fermiintegrals(
        epsilon, dos, vvdos, mur=mur, Tr=Tr, dosweight=data.dosweight
    )
    sigma, seebeck, kappa, _ = BL.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, vuc
    )

    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return peak / (1024 * 1024)  # Convert to MiB


def run_scaling_benchmark(
    datadir="benchmarks/data/Si.vasp",
    nkpt_values=[1000, 2000, 4000],
    n_runs=5,
    temperature=300.0,
    npts_dos=500,
):
    """Run the full scaling benchmark.

    Args:
        datadir: Path to VASP data directory
        nkpt_values: List of target k-point counts
        n_runs: Number of runs for timing (median reported)
        temperature: Temperature in K
        npts_dos: Number of DOS bins

    Returns:
        Results dictionary
    """
    print("=" * 60)
    print("Python BoltzTraP2 End-to-End Benchmark")
    print("=" * 60)
    print()

    if not HAVE_BOLTZTRAP2:
        print("ERROR: BoltzTraP2 not available.")
        return None

    print(f"pyFFTW available: {HAVE_PYFFTW}")
    print()

    # Load data (I/O excluded from measurement)
    print(f"Loading VASP data: {datadir}")
    data = BTP.DFTData(datadir)
    print(f"  Bands: {data.ebands.shape[0]}")
    print(f"  K-points (DFT): {len(data.kpoints)}")
    print(f"  Fermi: {data.fermi:.6f} Ha")
    print()

    # Warm up
    print("Warming up...")
    run_single_benchmark(data, 500, temperature, npts_dos)
    print()

    results = {
        "language": "Python",
        "version": platform.python_version(),
        "package": "BoltzTraP2",
        "pyfftw": HAVE_PYFFTW,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "temperature_K": temperature,
            "dos_bins": npts_dos,
            "n_runs": n_runs,
        },
        "results": []
    }

    print("-" * 60)
    print(f"Benchmarking (T={temperature}K, DOS bins={npts_dos})")
    print("-" * 60)
    print()

    for nkpt in nkpt_values:
        print(f"nkpt={nkpt}:")

        # Time measurements
        times_total = []
        times_interp = []
        times_integ = []

        for run in range(n_runs):
            total, interp, integ = run_single_benchmark(
                data, nkpt, temperature, npts_dos
            )
            times_total.append(total)
            times_interp.append(interp)
            times_integ.append(integ)
            print(f"  Run {run+1}: {total:.1f} ms (interp: {interp:.1f}, integ: {integ:.1f})")

        # Memory measurement (separate run)
        peak_mem = measure_memory(data, nkpt, temperature, npts_dos)

        median_total = np.median(times_total)
        median_interp = np.median(times_interp)
        median_integ = np.median(times_integ)
        std_total = np.std(times_total)

        print(f"  Median: {median_total:.1f} ms (interp: {median_interp:.1f}, integ: {median_integ:.1f})")
        print(f"  Peak memory: {peak_mem:.1f} MiB")
        print()

        results["results"].append({
            "kpoints": nkpt,
            "total_ms": median_total,
            "interpolate_ms": median_interp,
            "integrate_ms": median_integ,
            "std_ms": std_total,
            "peak_memory_MiB": peak_mem,
        })

    # Save results
    script_dir = Path(__file__).parent
    output_file = script_dir / "scaling_results_python.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved: {output_file}")
    print()

    # Summary table
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print()
    print(f"{'kpoints':>8} | {'Total (ms)':>12} | {'Interp (ms)':>12} | {'Integ (ms)':>12} | {'Memory (MiB)':>12}")
    print("-" * 70)
    for r in results["results"]:
        print(f"{r['kpoints']:>8} | {r['total_ms']:>12.1f} | {r['interpolate_ms']:>12.1f} | {r['integrate_ms']:>12.1f} | {r['peak_memory_MiB']:>12.1f}")
    print()

    return results


if __name__ == "__main__":
    run_scaling_benchmark()
