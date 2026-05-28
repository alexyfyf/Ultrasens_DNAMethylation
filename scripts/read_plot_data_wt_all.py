#!/usr/bin/env python3
"""Compatibility wrapper for analyze_cpg_density_switch_batch.py."""
from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("analyze_cpg_density_switch_batch.py")), run_name="__main__")
