#!/usr/bin/env python3
"""Compatibility wrapper for plot_cpg_switch_parameters.py."""
from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("plot_cpg_switch_parameters.py")), run_name="__main__")
