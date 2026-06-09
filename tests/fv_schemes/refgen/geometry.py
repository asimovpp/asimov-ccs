"""Geometry cases for the two-cell distorted stencil.

Ferziger, Peric, and Street discuss finite-volume face fluxes in Chapter 4
and the deferred correction terms needed on non-orthogonal grids in Chapter 9.
These cases start from one clean stencil and apply one transform at a time so
the Fortran tests can isolate skewness, non-orthogonality, and uneven spacing.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .math import vec, unit

ScalarPair = tuple[sp.Expr, sp.Expr]


@dataclass(frozen=True)
class CaseSpec:
    name: str
    xp: sp.Matrix
    xn: sp.Matrix
    xf: sp.Matrix
    normal: sp.Matrix
    volume: ScalarPair


def base_case() -> CaseSpec:
    """Base two-cell stencil with an orthogonal centred face."""

    return CaseSpec(
        "base",
        vec(0, 0, 0),
        vec(1, 0, 0),
        vec(sp.Rational(1, 2), 0, 0),
        vec(1, 0, 0),
        (sp.Rational(2, 5), sp.Rational(2, 5)),
    )


def skew(case: CaseSpec, offset: sp.Matrix) -> CaseSpec:
    """Move the face centre tangentially without changing cell centres."""

    normal_component = offset.dot(case.normal) / case.normal.dot(case.normal) * case.normal
    tangent_offset = sp.simplify(offset - normal_component)
    return CaseSpec(case.name, case.xp, case.xn, case.xf + tangent_offset, case.normal, case.volume)


def nonorth(case: CaseSpec, normal: sp.Matrix) -> CaseSpec:
    """Replace the face normal without moving cell or face centres."""

    return CaseSpec(case.name, case.xp, case.xn, case.xf, unit(normal), case.volume)


def uneven(case: CaseSpec) -> CaseSpec:
    """Stretch neighbour spacing and use unequal control-volume sizes."""

    return CaseSpec(
        case.name,
        case.xp,
        vec(sp.Rational(6, 5), 0, 0),
        vec(sp.Rational(2, 5), 0, 0),
        case.normal,
        (sp.Rational(31, 100), sp.Rational(47, 100)),
    )


def named(name: str, case: CaseSpec) -> CaseSpec:
    return CaseSpec(name, case.xp, case.xn, case.xf, case.normal, case.volume)


def build_cases() -> list[CaseSpec]:
    """Return base, isolated distortion, and combined cases."""

    base = base_case()
    skew_offset = vec(0, sp.Rational(1, 5), sp.Rational(1, 7))
    nonorth_normal = vec(sp.Rational(4, 5), sp.Rational(3, 5), sp.Rational(1, 5))

    return [
        named("base", base),
        named("skew", skew(base, skew_offset)),
        named("nonorth", nonorth(base, nonorth_normal)),
        named("uneven", uneven(base)),
        named("all", nonorth(skew(uneven(base), skew_offset), nonorth_normal)),
    ]
