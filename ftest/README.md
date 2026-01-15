# Functional Tests (ftest)

Manual functional tests for CLI commands and end-to-end workflows.

These tests are **not** included in `Pkg.test()` and must be run manually. They test CLI functionality, argument propagation, and integration with external tools like DFTK.

## Prerequisites

```bash
# BoltzTraP2 test data
ls ../BoltzTraP2-public/data/Si.vasp/vasprun.xml

# For DFTK tests
julia -e 'using Pkg; Pkg.add("DFTK")'
```

## Scripts

| Script | Description | Requirements |
|--------|-------------|--------------|
| `test_cli_all_options.jl` | Test all CLI options for interpolate/integrate | VASP data |
| `test_cli_args.jl` | Test CLI argument propagation to workflow functions | VASP data |
| `test_cli_bt2.jl` | Test CLI commands with Python-generated .bt2 files | .bt2 file |
| `test_format_detection.jl` | Test DFT format auto-detection | Multiple formats |
| `test_dftk_e2e.jl` | End-to-end DFTK integration test | DFTK.jl |

## Usage

```bash
cd BoltzTraP.jl

# Run all CLI option tests
julia --project ftest/test_cli_all_options.jl

# Run with custom VASP directory
julia --project ftest/test_cli_all_options.jl /path/to/vasp/data

# Test CLI argument propagation
julia --project ftest/test_cli_args.jl

# Test .bt2 file handling
julia --project ftest/test_cli_bt2.jl

# Test format auto-detection
julia --project ftest/test_format_detection.jl

# DFTK end-to-end test (requires DFTK.jl, slow)
julia --project ftest/test_dftk_e2e.jl
```

### Debug Output

Enable debug logging to see detailed output:

```bash
JULIA_DEBUG=BoltzTraP julia --project ftest/test_cli_args.jl
```

## Difference from test/

| Directory | Purpose | Execution | Included in CI |
|-----------|---------|-----------|----------------|
| `test/` | Unit tests, integration tests | `Pkg.test()` | Yes |
| `ftest/` | Manual functional tests, CLI tests | Manual | No |

Functional tests are kept separate because:
- They require manual verification of output
- Some tests are slow (e.g., DFTK SCF calculation)
- They test CLI behavior that's hard to automate

## Notes

- All scripts print `✓ PASS` or `✗ FAIL` for each test case
- Exit code is non-zero on failure
- DFTK test takes ~1 minute due to SCF calculation
