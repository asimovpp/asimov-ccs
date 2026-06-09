#!/usr/bin/env python3
"""Generate reference values for distorted finite-volume test cases.

Geometry is built from a clean two-cell stencil and then transformed to
isolate skewness, non-orthogonality, and uneven spacing.  Field values are
manufactured point solutions, so the Fortran tests can check the discrete
kernel and equation algebra without relying on hand-written constants.
SymPy is a generator-only dependency used for exact algebra and Fortran
printing; the generated Fortran module has no Python runtime dependency.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path

from refgen.fields import build_fields
from refgen.fortran import render_module
from refgen.geometry import build_cases
from refgen.refs import field_refs, geometry_refs


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "distorted_refs.f90"


def render_refs() -> str:
    cases = build_cases()
    fields = build_fields()
    case_data = [geometry_refs(case) for case in cases]
    field_data = [[field_refs(case_ref, field) for field in fields] for case_ref in case_data]
    return render_module(cases, fields, case_data, field_data)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate distorted-stencil Fortran references.")
    parser.add_argument("--write", action="store_true", help="write the generated Fortran module")
    parser.add_argument("--check", action="store_true", help="check that the generated module is current")
    args = parser.parse_args()

    if not args.write and not args.check:
        parser.error("choose --write or --check")

    generated = render_refs()
    if args.write:
        OUTPUT.write_text(generated, encoding="utf-8")

    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != generated:
            diff = difflib.unified_diff(
                current.splitlines(),
                generated.splitlines(),
                fromfile=str(OUTPUT),
                tofile="generated",
                lineterm="",
            )
            print("\n".join(diff))
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
