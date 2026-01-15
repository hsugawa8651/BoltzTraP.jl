# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""Common utilities for reftest scripts."""

import argparse
import os
import sys
from pathlib import Path

# Output directory for generated .npz files
OUTPUT_DIR = Path(__file__).parent / "data"
OUTPUT_DIR.mkdir(exist_ok=True)


def get_boltztrap2_data_dir():
    """Get BoltzTraP2 data directory.

    Priority:
    1. --data-dir command line argument
    2. pip installed BoltzTraP2 package location

    Returns:
        Path to BoltzTraP2 data directory

    Usage:
        python generate_*.py                          # Use pip installed location
        python generate_*.py --data-dir /path/to/data # Use specified path
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--data-dir', type=str, default=None,
                        help='Path to BoltzTraP2 data directory')
    args, _ = parser.parse_known_args()

    if args.data_dir:
        data_dir = Path(args.data_dir)
        if not data_dir.exists():
            print(f"Error: Specified data directory does not exist: {data_dir}")
            sys.exit(1)
        return data_dir

    # Default: use pip installed BoltzTraP2
    try:
        import BoltzTraP2
        data_dir = Path(BoltzTraP2.__path__[0]) / "data"
        if not data_dir.exists():
            print(f"Error: BoltzTraP2 data directory not found: {data_dir}")
            print("Try: pip install boltztrap2")
            sys.exit(1)
        return data_dir
    except ImportError:
        print("Error: BoltzTraP2 not installed.")
        print("Install with: pip install boltztrap2")
        print("Or specify data directory: --data-dir /path/to/BoltzTraP2/data")
        sys.exit(1)


def get_material_path(material_name):
    """Get path to a specific material directory.

    Args:
        material_name: e.g., "Si.vasp", "PbTe.vasp.unpolarized", "Li.W2K"

    Returns:
        Path to material directory
    """
    data_dir = get_boltztrap2_data_dir()
    return data_dir / material_name
