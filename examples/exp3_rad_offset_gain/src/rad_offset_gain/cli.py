#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os

from rad_offset_gain.lib import add_radiance_offset


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply a radiance offset to L1B NetCDF data.")

    parser.add_argument("--input", required=True, help="Input NetCDF file path")
    parser.add_argument("--output", required=True, help="Output NetCDF file path")
    parser.add_argument(
        "--offset", type=float, default=0.01, help="Offset fraction (e.g., 0.01 for 1%%)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    try:
        add_radiance_offset(args.input, args.output, args.offset)
    except Exception as e:
        raise RuntimeError(f"Error processing file {args.input}: {e}")


if __name__ == "__main__":
    main()
