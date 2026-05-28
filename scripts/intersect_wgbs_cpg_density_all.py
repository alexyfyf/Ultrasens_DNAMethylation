#!/usr/bin/env python3
"""Compatibility wrapper for intersect_known_wgbs_cpg_density.py."""
from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("intersect_known_wgbs_cpg_density.py")), run_name="__main__")
