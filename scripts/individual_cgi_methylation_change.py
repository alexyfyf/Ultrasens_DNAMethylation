#!/usr/bin/env python3
"""Compatibility wrapper for compare_cgi_methylation_change.py."""
from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("compare_cgi_methylation_change.py")), run_name="__main__")
