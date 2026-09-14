"""Render generated reference data as a Fortran module."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import sympy as sp
from sympy.codegen.ast import Declaration, FloatBaseType, IntBaseType, String, Variable, value_const
from sympy.codegen.fnodes import ArrayConstructor, array, reshape
from sympy.codegen.futils import FCodePrinter

from .geometry import CaseSpec
from .fields import FieldSpec
from .math import as_array

REAL_CCS = FloatBaseType(String("real(ccs_real)"))
INT_DEFAULT = IntBaseType(String("integer"))
N_CASES = sp.Symbol("n_cases")
N_FIELDS = sp.Symbol("n_fields")
ReferenceValue = sp.Expr | sp.Matrix | tuple[sp.Expr, ...]
CASE_SCALARS = (
    "area",
    "schmidt",
    "mflux",
    "lf",
    "dx_mag",
    "dx_orth",
    "dx_face",
    "gamma_f",
    "diff_plain",
    "diff_corr",
    "adv_flux",
    "vol_f",
    "inv_a_f",
    "pc_coeff",
)
CASE_ARRAYS = (
    ("xp", (3, N_CASES)),
    ("xn", (3, N_CASES)),
    ("xf", (3, N_CASES)),
    ("normal", (3, N_CASES)),
    ("volume", (2, N_CASES)),
    ("mu", (2, N_CASES)),
    ("rho", (2, N_CASES)),
    ("inv_a", (2, N_CASES)),
    ("xp_orth", (3, N_CASES)),
    ("xn_orth", (3, N_CASES)),
    ("rvec", (3, 2, N_CASES)),
)
FIELD_SCALARS = (
    "diff_exp",
    "cd_face",
    "corr_face",
    "cd_exp",
    "upwind_exp",
    "luds_exp",
    "face_corr",
    "mom_adv_res",
    "mom_diff_res",
    "mom_res",
    "pc_res",
    "rc_flux",
)
FIELD_ARRAYS = (
    ("phi_val", (2, N_FIELDS, N_CASES)),
    ("pressure", (2, N_FIELDS, N_CASES)),
    ("grad_phi", (3, 2, N_FIELDS, N_CASES)),
    ("grad_p", (3, 2, N_FIELDS, N_CASES)),
)


class CCSPrinter(FCodePrinter):
    """Fortran printer that emits CCS kinded real literals."""

    def _print_Float(self, expr: sp.Float) -> str:
        return f"{float(expr):.17e}_ccs_real"


def real_lit(value: sp.Expr) -> sp.Float:
    return sp.Float(sp.N(value, 30), 17)


def int_param(name: str, value: int) -> Declaration:
    var = Variable(name, type=INT_DEFAULT, attrs=[value_const], value=sp.Integer(value))
    return var.as_Declaration()


def real_array(name: str, values: np.ndarray, decl_shape: tuple[sp.Expr | int, ...]) -> Declaration:
    data = np.asarray(values, dtype=object)
    flat = [real_lit(value) for value in data.reshape(-1, order="F")]
    constructor = ArrayConstructor(flat)
    if data.ndim == 1:
        value = constructor
    else:
        value = reshape(constructor, ArrayConstructor([sp.sympify(dim) for dim in decl_shape]))

    var = array(name, decl_shape, type=REAL_CCS, attrs=[value_const], value=value)
    return var.as_Declaration()


def case_array(value: ReferenceValue) -> np.ndarray:
    if isinstance(value, tuple):
        return np.asarray(value, dtype=object)
    return as_array(value)


def collect_case(refs: Sequence[dict[str, ReferenceValue]], name: str) -> np.ndarray:
    values = [case_array(case[name]) for case in refs]
    if values[0].shape == ():
        return np.asarray([value.item() for value in values], dtype=object)
    return np.stack(values, axis=-1)


def collect_field(refs: Sequence[Sequence[dict[str, sp.Expr | sp.Matrix]]], name: str) -> np.ndarray:
    values = [[as_array(field[name]) for field in case] for case in refs]
    if values[0][0].shape == ():
        return np.asarray([[value.item() for value in case] for case in values], dtype=object).T
    return np.stack([np.stack(case, axis=-1) for case in values], axis=-1)


def declarations(
    cases: Sequence[CaseSpec],
    fields: Sequence[FieldSpec],
    case_refs: Sequence[dict[str, ReferenceValue]],
    field_refs: Sequence[Sequence[dict[str, sp.Expr | sp.Matrix]]],
) -> list[str | Declaration]:
    decls: list[str | Declaration] = [int_param(f"case_{case.name}", i + 1) for i, case in enumerate(cases)]
    decls.extend(int_param(f"field_{field.name}", i + 1) for i, field in enumerate(fields))
    decls.append(int_param("n_cases", len(cases)))
    decls.append(int_param("n_fields", len(fields)))

    decls.append("")
    decls.append("! Non-unit material and flow values make interpolation mistakes visible.")
    for name in CASE_SCALARS:
        decls.append(real_array(name, collect_case(case_refs, name), (N_CASES,)))

    decls.append("")
    decls.append("! Geometry for the owner cell, neighbour cell, and shared face.")
    for name, shape in CASE_ARRAYS:
        decls.append(real_array(name, collect_case(case_refs, name), shape))

    decls.append("")
    decls.append("! Kernel explicit terms and assembled residual references.")
    for name in FIELD_SCALARS:
        decls.append(real_array(name, collect_field(field_refs, name), (N_FIELDS, N_CASES)))

    decls.append("")
    decls.append("! Sampled values and exact gradients for the analytic fields.")
    for name, shape in FIELD_ARRAYS:
        decls.append(real_array(name, collect_field(field_refs, name), shape))

    return decls


def render_module(
    cases: Sequence[CaseSpec],
    fields: Sequence[FieldSpec],
    case_refs: Sequence[dict[str, ReferenceValue]],
    field_refs: Sequence[Sequence[dict[str, sp.Expr | sp.Matrix]]],
) -> str:
    printer = CCSPrinter({"standard": 2003, "source_format": "free"})
    lines = [
        "!> Reference values generated by generate_distorted_refs.py.",
        "!> Do not edit this file by hand; update the generator instead.",
        "!> SymPy is used only by the generator for exact algebra and printing.",
        "module distorted_refs",
        "  use kinds, only: ccs_real",
        "  implicit none",
        "  public",
        "",
    ]

    for decl in declarations(cases, fields, case_refs, field_refs):
        if isinstance(decl, str):
            lines.append(f"  {decl}" if decl else "")
        else:
            printed = printer.doprint(decl).replace("integer*4", "integer")
            lines.extend(f"  {line}" if line else "" for line in printed.splitlines())

    lines.extend(["", "end module distorted_refs", ""])
    return "\n".join(lines)
