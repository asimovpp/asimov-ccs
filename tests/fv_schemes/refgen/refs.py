"""Finite-volume reference formulas for distorted two-cell stencils."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .fields import FieldSpec, sample_field
from .geometry import CaseSpec
from .math import pick_min

ScalarPair = tuple[sp.Expr, sp.Expr]
ReferenceValue = sp.Expr | sp.Matrix | ScalarPair


@dataclass(frozen=True)
class MaterialSpec:
    """Material and solver values shared by every distorted geometry case."""

    area: sp.Expr
    mu: ScalarPair
    rho: ScalarPair
    schmidt: sp.Expr
    mflux: sp.Expr
    inv_a: ScalarPair


def default_materials() -> MaterialSpec:
    """Return deliberately asymmetric values for interpolation checks."""

    # These non-unit, non-integer values make coefficient and interpolation
    # mistakes visible instead of hiding them behind symmetric test data.
    return MaterialSpec(
        area=sp.Rational(73, 100),
        mu=(sp.Rational(13, 10), sp.Rational(19, 20)),
        rho=(sp.Rational(9, 10), sp.Rational(23, 20)),
        schmidt=sp.Rational(37, 50),
        mflux=sp.Rational(37, 100),
        inv_a=(sp.Rational(6, 5), sp.Rational(17, 20)),
    )


def geometry_refs(spec: CaseSpec) -> dict[str, ReferenceValue]:
    """Compute geometry and material terms that do not depend on a field."""

    materials = default_materials()
    dx_pn = spec.xn - spec.xp
    dx_mag = sp.simplify(dx_pn.norm())
    dx_mag_sq = sp.simplify(dx_mag**2)
    lf = sp.simplify(1 - dx_pn.dot(spec.xf - spec.xp) / dx_mag_sq)
    dx_orth = pick_min((spec.xf - spec.xp).dot(spec.normal), (spec.xn - spec.xf).dot(spec.normal))
    dx_face = sp.simplify(2 * dx_orth)
    xp_orth = spec.xf - dx_orth * spec.normal
    xn_orth = spec.xf + dx_orth * spec.normal
    rvec = sp.Matrix.hstack(xp_orth - spec.xp, xn_orth - spec.xn)

    mu_f = lf * materials.mu[0] + (1 - lf) * materials.mu[1]
    rho_f = lf * materials.rho[0] + (1 - lf) * materials.rho[1]
    gamma_f = sp.simplify(mu_f / (rho_f * materials.schmidt))
    diff_plain = sp.simplify(materials.area * gamma_f / dx_mag)
    diff_corr = sp.simplify(materials.area * gamma_f / dx_face)
    adv_flux = sp.simplify(materials.mflux * materials.area)
    vol_f = sp.simplify(lf * spec.volume[0] + (1 - lf) * spec.volume[1])
    inv_a_f = sp.simplify(lf * materials.inv_a[0] + (1 - lf) * materials.inv_a[1])
    pc_coeff = sp.simplify(materials.area * vol_f * inv_a_f / dx_face)

    return {
        "xp": spec.xp,
        "xn": spec.xn,
        "xf": spec.xf,
        "normal": spec.normal,
        "area": materials.area,
        "volume": spec.volume,
        "mu": materials.mu,
        "rho": materials.rho,
        "schmidt": materials.schmidt,
        "inv_a": materials.inv_a,
        "mflux": materials.mflux,
        "lf": lf,
        "dx_mag": dx_mag,
        "dx_orth": dx_orth,
        "dx_face": dx_face,
        "xp_orth": xp_orth,
        "xn_orth": xn_orth,
        "rvec": rvec,
        "gamma_f": gamma_f,
        "diff_plain": diff_plain,
        "diff_corr": diff_corr,
        "adv_flux": adv_flux,
        "vol_f": vol_f,
        "inv_a_f": inv_a_f,
        "pc_coeff": pc_coeff,
    }


def field_refs(case_ref: dict[str, ReferenceValue], field: FieldSpec) -> dict[str, sp.Expr | sp.Matrix]:
    """Compute field-dependent kernel and equation references."""

    sampled = sample_field(field, case_ref["xp"], case_ref["xn"])
    phi_val = sampled["phi_val"]
    grad_phi = sampled["grad_phi"]
    pressure = sampled["pressure"]
    grad_p = sampled["grad_p"]
    rvec = case_ref["rvec"]

    lf = case_ref["lf"]
    diff_corr = case_ref["diff_corr"]
    adv_flux = case_ref["adv_flux"]
    pc_coeff = case_ref["pc_coeff"]
    vol_f = case_ref["vol_f"]
    inv_a_f = case_ref["inv_a_f"]
    dx_face = case_ref["dx_face"]
    xp = case_ref["xp"]
    xn = case_ref["xn"]
    dx_pn = xn - xp

    diff_exp = sp.simplify(-diff_corr * (grad_phi[:, 1].dot(rvec[:, 1]) - grad_phi[:, 0].dot(rvec[:, 0])))
    cd_face = sp.simplify(lf * phi_val[0] + (1 - lf) * phi_val[1])
    corr_face = sp.simplify(
        sp.Rational(1, 2) * (phi_val[0] + phi_val[1])
        + sp.Rational(1, 2) * (grad_phi[:, 0].dot(rvec[:, 0]) + grad_phi[:, 1].dot(rvec[:, 1]))
    )
    cd_exp = sp.simplify(adv_flux * cd_face - adv_flux * phi_val[0])
    luds_exp = sp.simplify(adv_flux * grad_phi[:, 0].dot(rvec[:, 0]))
    face_corr = sp.simplify(corr_face - sp.Rational(1, 2) * (phi_val[0] + phi_val[1]))
    pc_explicit = sp.simplify(
        -pc_coeff * (grad_phi[:, 1].dot(rvec[:, 1]) - grad_phi[:, 0].dot(rvec[:, 0]))
    )

    # Rhie-Chow mass flux correction: pressure jump, non-orthogonal
    # gradient correction, then volume/inverse-diagonal interpolation.
    pressure_flux = sp.simplify(
        (pressure[1] - pressure[0])
        + (grad_p[:, 1].dot(rvec[:, 1]) - grad_p[:, 0].dot(rvec[:, 0]))
        - sp.Rational(1, 2) * (grad_p[:, 0] + grad_p[:, 1]).dot(dx_pn)
    )
    pressure_flux = sp.simplify(-pressure_flux / dx_face)
    mom_adv_rhs = sp.simplify(-cd_exp - face_corr * adv_flux)
    mom_diff_rhs = sp.simplify(-diff_exp)
    mom_adv_res = sp.simplify(mom_adv_rhs - adv_flux * phi_val[0])
    mom_diff_res = sp.simplify(mom_diff_rhs - diff_corr * phi_val[0] + diff_corr * phi_val[1])

    return {
        **sampled,
        "diff_exp": diff_exp,
        "cd_face": cd_face,
        "corr_face": corr_face,
        "cd_exp": cd_exp,
        "upwind_exp": sp.Integer(0),
        "luds_exp": luds_exp,
        "face_corr": face_corr,
        "mom_adv_res": mom_adv_res,
        "mom_diff_res": mom_diff_res,
        "mom_res": sp.simplify(mom_adv_res + mom_diff_res),
        "pc_res": sp.simplify(pc_explicit - pc_coeff * phi_val[0] + pc_coeff * phi_val[1]),
        "rc_flux": sp.simplify(vol_f * inv_a_f * pressure_flux),
    }
