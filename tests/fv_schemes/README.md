# Distorted FV References

This directory owns the generated references for the distorted
finite-volume tests.

`generate_distorted_refs.py` is the source for `distorted_refs.f90`. The
Fortran module is generated output, checked in with the tests, and consumed
directly by the Fortran test programs. Do not edit generated constants by
hand. Change the Python generator or its inputs, regenerate the module, and
review the generator diff and generated Fortran diff together.

## Files

- `generate_distorted_refs.py`: command-line entry point.
- `refgen/geometry.py`: two-cell stencil geometry and distortion cases.
- `refgen/fields.py`: analytic scalar and pressure fields with exact
  gradients.
- `refgen/refs.py`: finite-volume kernel and equation reference formulas.
- `refgen/fortran.py`: Fortran module rendering.
- `distorted_refs.f90`: generated Fortran constants used by the tests.
- `pyproject.toml`: Python dependency record for this harness.

## Dependencies

Use Python 3 with NumPy and SymPy available in the active environment.
`pyproject.toml` records those dependencies for the generator. In this
repository, use plain `python3` commands for the documented workflow.

## Cases and Fields

The generator writes five cases:

- `base`: undistorted two-cell stencil.
- `skew`: face centre shifted tangentially from the cell-centre line.
- `nonorth`: face normal rotated away from the cell-centre line.
- `uneven`: unequal owner and neighbour spacing and volumes.
- `all`: skewness, non-orthogonality, and uneven spacing together.

The generator writes three field families:

- `lin`: linear scalar and pressure variation.
- `quad`: quadratic variation, including x/y/z cross terms.
- `trig`: smooth trigonometric variation in all three coordinates.

Together these references exercise skewness, non-orthogonality, uneven
spacing, and 3-D field and gradient terms.

## Update Workflow

Run from the repository root:

```sh
python3 tests/fv_schemes/generate_distorted_refs.py --check
python3 tests/fv_schemes/generate_distorted_refs.py --write
python3 tests/fv_schemes/generate_distorted_refs.py --check
python3 -m compileall -q tests/fv_schemes/generate_distorted_refs.py tests/fv_schemes/refgen
```

Use `--check` before editing to confirm the checked-in Fortran is current.
Use `--write` only after an intentional generator change. The second
`--check` confirms that `distorted_refs.f90` matches the generator output.

## Test Consumers

- `tests/kernels/test_distorted_kernels.f90` checks kernel-level distorted
  stencil terms against the generated references.
- `tests/fv_schemes/test_distorted_schemes.f90` checks equation assembly,
  mass-flux behavior, and residual reductions on the same generated stencils.
