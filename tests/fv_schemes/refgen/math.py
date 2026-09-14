"""Small symbolic helpers shared by the reference generator."""

from __future__ import annotations

import numpy as np
import sympy as sp


def vec(*values: sp.Expr | int | float) -> sp.Matrix:
    """Build a SymPy column vector from scalar-like values."""

    return sp.Matrix([sp.sympify(value) for value in values])


def unit(value: sp.Matrix) -> sp.Matrix:
    """Return a unit vector with exact symbolic components."""

    return sp.simplify(value / value.norm())


def pick_min(*values: sp.Expr) -> sp.Expr:
    """Choose the smaller symbolic value using high precision evaluation."""

    return min(values, key=lambda value: float(sp.N(value, 30)))


def as_array(value: sp.Expr | sp.Matrix) -> np.ndarray:
    """Convert scalar or matrix reference values to an object array."""

    if isinstance(value, sp.MatrixBase):
        return np.array(value.tolist(), dtype=object)
    return np.asarray(value, dtype=object)
