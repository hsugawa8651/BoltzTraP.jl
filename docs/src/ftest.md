# Functional Tests

BoltzTraP.jl provides functional tests (ftest) for verifying CLI commands and end-to-end workflows. These are standalone scripts for manual execution, useful for:

- **Installation verification**: Confirm BoltzTraP.jl works correctly after installation
- **Troubleshooting**: Diagnose issues with CLI commands or format detection
- **Integration testing**: Verify DFTK integration works in your environment

## Overview

| Aspect | Description |
|--------|-------------|
| Location | `ftest/` directory |
| Execution | Manual (`julia --project ftest/script.jl`) |
| Included in `Pkg.test()` | No |
| Requirements | BoltzTraP2-public test data (for most scripts) |

!!! note "ftest vs test"
    - `test/`: Automated unit tests, run via `Pkg.test()`
    - `ftest/`: Manual functional tests, run individually

## Available Scripts

| Script | Purpose | Requirements |
|--------|---------|--------------|
| `test_cli_all_options.jl` | Test all CLI options for interpolate/integrate | VASP data |
| `test_cli_args.jl` | Test CLI argument propagation | VASP data |
| `test_cli_bt2.jl` | Test CLI with Python-generated `.bt2` files | `.bt2` file |
| `test_format_detection.jl` | Test DFT format auto-detection | Multiple formats |
| `test_dftk_e2e.jl` | End-to-end DFTK integration | DFTK.jl |

## Quick Start

### Verify Installation

After installing BoltzTraP.jl, run the format detection test:

```bash
cd BoltzTraP.jl
julia --project ftest/test_format_detection.jl
```

Expected output:
```
============================================================
Functional Test: DFT Format Auto-Detection
============================================================

------------------------------------------------------------
Test 1: Detect VASP format
------------------------------------------------------------
  Directory: ../BoltzTraP2-public/data/Si.vasp
  Detected: VASP
  ✓ PASS
...
```

### Test CLI Commands

```bash
# Test all CLI options (comprehensive)
julia --project ftest/test_cli_all_options.jl

# Test with custom data directory
julia --project ftest/test_cli_all_options.jl /path/to/vasp/data

# Basic argument test
julia --project ftest/test_cli_args.jl
```

### Test DFTK Integration

```bash
# Requires DFTK.jl (takes ~1 minute for SCF)
julia --project ftest/test_dftk_e2e.jl
```

This runs a complete workflow: DFTK SCF → [`load_dftk`](@ref) → [`run_interpolate`](@ref) → [`run_integrate`](@ref).

### Test .bt2 File Support

```bash
# Test Python BoltzTraP2 .bt2 file compatibility
julia --project ftest/test_cli_bt2.jl
```

## Debug Output

Enable debug logging for detailed output:

```bash
JULIA_DEBUG=BoltzTraP julia --project ftest/test_cli_args.jl
```

This shows:
- CLI argument values received
- DFT format detection process
- File loading details
- Workflow function parameters

## Test Data

Most functional tests expect BoltzTraP2 test data at `../BoltzTraP2-public/data/`:

```
BoltzTraP2-public/
└── data/
    ├── Si.vasp/
    │   └── vasprun.xml
    ├── Si.ESPRESSO/
    │   └── out/
    └── ...
```

You can specify a different directory as a command-line argument:

```bash
julia --project ftest/test_cli_all_options.jl /my/custom/vasp/dir
```

## Troubleshooting

### "No compatible DFT format found"

1. Run format detection test with debug output:
   ```bash
   JULIA_DEBUG=BoltzTraP julia --project ftest/test_format_detection.jl
   ```

2. Check if required files exist:
   - VASP: `vasprun.xml`
   - QE: `*.save/data-file-schema.xml`
   - Wien2k: `*.struct`, `*.energy`, `*.scf`

3. Try explicit format:
   ```julia
   data = load_vasp("/path/to/dir")  # Instead of load_dft
   ```

### DFTK Test Fails

1. Ensure DFTK is installed:
   ```julia
   using Pkg; Pkg.add("DFTK")
   ```

2. Load DFTK before BoltzTraP:
   ```julia
   using DFTK
   using BoltzTraP
   ```

3. Check DFTK version compatibility

### Test Script Not Found

Ensure you're in the BoltzTraP.jl directory:

```bash
cd /path/to/BoltzTraP.jl
julia --project ftest/test_format_detection.jl
```

## Writing Custom Functional Tests

Template for a new functional test:

```julia
# ftest/test_my_feature.jl
using Logging

# Enable debug logging (optional)
ENV["JULIA_DEBUG"] = "BoltzTraP"
global_logger(ConsoleLogger(stderr, Logging.Debug))

using BoltzTraP

function main()
    println("=" ^ 60)
    println("Functional Test: My Feature")
    println("=" ^ 60)
    println()

    # Test 1
    println("-" ^ 60)
    println("Test 1: Description")
    println("-" ^ 60)

    result = some_function()
    @assert condition "Error message"
    println("  ✓ PASS")
    println()

    # Summary
    println("=" ^ 60)
    println("All tests passed!")
    println("=" ^ 60)
end

main()
```

## See Also

- [Validation](@ref) - Visual comparison with Python BoltzTraP2
- [Reference Tests](@ref) - Numerical equivalence tests
