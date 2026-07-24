#!/usr/bin/env python3
"""Compatibility wrapper for ``busco-alg-painter paint``."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from busco_alg_painter.cli import paint_main

if __name__ == "__main__":
    paint_main()
