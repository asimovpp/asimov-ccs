"""Manufactured point fields sampled by the distorted stencil tests.

The fields are analytic so SymPy can derive exact gradients before the
generator writes decimal Fortran constants.  Keeping the fields analytic makes
the references easier to audit while still exercising linear, curved, and
trigonometric variation across the same distorted control volumes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import sympy as sp

from .math import vec

X, Y, Z = sp.symbols("x y z")
ScalarSampler = Callable[[sp.Matrix], tuple[sp.Expr, sp.Matrix]]


@dataclass(frozen=True)
class FieldSpec:
    name: str
    phi: ScalarSampler
    pressure: ScalarSampler


def symbolic_sampler(expr: sp.Expr) -> ScalarSampler:
    """Return a sampler for a symbolic scalar field and its gradient."""

    grad_expr = vec(*(sp.diff(expr, var) for var in (X, Y, Z)))

    def sample(point: sp.Matrix) -> tuple[sp.Expr, sp.Matrix]:
        subs = {X: point[0], Y: point[1], Z: point[2]}
        value = sp.simplify(expr.subs(subs))
        gradient = vec(*(sp.simplify(component.subs(subs)) for component in grad_expr))
        return value, gradient

    return sample


def field_spec(name: str, phi: sp.Expr, pressure: sp.Expr) -> FieldSpec:
    """Pair scalar and pressure fields with their exact symbolic gradients."""

    return FieldSpec(name, symbolic_sampler(phi), symbolic_sampler(pressure))


def build_fields() -> list[FieldSpec]:
    """Return field families used to probe the same distorted stencils."""

    phi_lin = (
        sp.Integer(1)
        + sp.Rational(7, 10) * X
        - sp.Rational(1, 5) * Y
        + sp.Rational(13, 100) * Z
    )
    p_lin = sp.Rational(1, 4) + sp.Rational(4, 5) * X - sp.Rational(2, 5) * Y + sp.Rational(3, 20) * Z

    phi_quad = (
        phi_lin
        + sp.Rational(3, 20) * X**2
        + sp.Rational(1, 2) * X * Y
        - sp.Rational(3, 10) * X * Z
        + sp.Rational(1, 5) * Y**2
        + sp.Rational(7, 25) * Y * Z
        - sp.Rational(11, 50) * Z**2
    )
    p_quad = (
        sp.Rational(1, 4)
        + sp.Rational(4, 5) * X
        - sp.Rational(2, 5) * Y
        - sp.Rational(1, 20) * X**2
        + sp.Rational(1, 5) * X * Y
        - sp.Rational(7, 50) * X * Z
        + sp.Rational(1, 10) * Y**2
        + sp.Rational(3, 25) * Y * Z
        + sp.Rational(9, 100) * Z**2
    )

    phi_trig = (
        sp.Rational(9, 10)
        + sp.sin(sp.pi * X / 3)
        - sp.Rational(1, 3) * sp.cos(sp.pi * Y / 4)
        + sp.Rational(1, 5) * sp.sin(sp.pi * Z / 6)
    )
    p_trig = (
        sp.Rational(2, 5)
        + sp.Rational(3, 5) * sp.sin(sp.pi * X / 4)
        + sp.Rational(1, 4) * sp.cos(sp.pi * Y / 5)
        - sp.Rational(1, 6) * sp.cos(sp.pi * Z / 7)
    )

    return [
        field_spec("lin", phi_lin, p_lin),
        field_spec("quad", phi_quad, p_quad),
        field_spec("trig", phi_trig, p_trig),
    ]


def sample_field(field: FieldSpec, xp: sp.Matrix, xn: sp.Matrix) -> dict[str, sp.Expr | sp.Matrix]:
    """Sample scalar and pressure fields at the owner and neighbour cells."""

    phi_p, gphi_p = field.phi(xp)
    phi_n, gphi_n = field.phi(xn)
    p_p, gp_p = field.pressure(xp)
    p_n, gp_n = field.pressure(xn)

    return {
        "phi_val": vec(phi_p, phi_n),
        "grad_phi": sp.Matrix.hstack(gphi_p, gphi_n),
        "pressure": vec(p_p, p_n),
        "grad_p": sp.Matrix.hstack(gp_p, gp_n),
    }
