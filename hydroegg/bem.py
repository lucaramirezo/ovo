"""Blade Element Momentum (BEM) solver and constraint-aware blade design for the
HydroEgg rotor.

Run ``python -m hydroegg.bem`` to regenerate ``data/blade_design.json`` and the
``docs/figures/bem_*.png`` figures. The solver half of this module (sections,
polars, ``Rotor``, ``solve``) is reusable on its own; the design study is below
the ``Design study`` banner.

Conventions
-----------
Frame: X = body axis, flow travels +X, nose at x = 0, rotor downstream at the tail,
Z up. ``r`` is always the radius measured perpendicular to the body axis.
Units: mm and deg for geometry, SI for physics (function arguments say which).

Cone ``beta``: aft (+X, downstream) lean of the blade pitch axis out of the plane
of rotation. Twist ``theta``: angle between the chord line and the tangential
(azimuthal) direction, measured in the section plane normal to the pitch axis,
positive when the leading edge sits toward upstream (-X). The inflow angle
``phi`` is measured the same way, so the angle of attack is ``alpha = phi - theta``.

Momentum model for a coned element: axial momentum through the projected annulus
2 pi r dr, the spanwise (radial) component of the inflow along the coned blade is
neglected (independence principle):

    tan(phi)     = U (1 - a) cos(beta) / (Omega r (1 + a'))
    a / (1 - a)  = sigma Cn cos^2(beta) / (4 F sin^2 phi)    Buhl branch above a = 0.4
    a' / (1 + a')= sigma Ct / (4 F sin(phi) cos(phi))
    sigma = B c / (2 pi r),  Cn = Cl cos phi + Cd sin phi,  Ct = Cl sin phi - Cd cos phi

F is the product of the Prandtl tip and hub loss factors. With beta = 0 this is the
textbook BEM. The high-induction branch is Buhl (2005, NREL/TP-500-36834) written
in the closed form of Ning (2014, Wind Energy 17:1327). Stations that the relaxed
fixed-point iteration cannot converge fall back to Ning's bracketed residual in
phi; anything still unconverged is flagged in the result.
"""

from __future__ import annotations

import functools
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from . import params as P

ROOT = Path(__file__).resolve().parents[1]
DATA_JSON = ROOT / "data" / "blade_design.json"
FIG_DIR = ROOT / "docs" / "figures"

RHO = P.WATER_RHO            # kg/m^3
NU = P.WATER_NU              # m^2/s
G = 9.81                     # m/s^2

# Section family: NACA 44xx (4 % camber at 40 % chord), thickness blended root -> tip.
CAMBER_M = 0.04
CAMBER_P = 0.40
TC_ROOT = 0.18
TC_TIP = 0.12

# Polar model assumptions (see docs/notes/blade-design.md).
N_CRIT = 5.0                 # disturbed inflow (river turbulence, printed surface)
MODEL_SIZE = "xlarge"

# Geometry inputs owned by this study (proposed for params.py, see the JSON).
R_AERO_ROOT = 30.0           # mm, start of the aerodynamic span, perpendicular radius
MAX_CHORD = 28.0             # mm, three stowed blades side by side around the tail
LUG_CHORD = 22.0             # mm, rounded root lug that sits in the hub clevis
LUG_THICKNESS = 8.0          # mm, 3 mm pin bore + 2.5 mm walls
CONE_NOMINAL = 12.0          # deg, straight pitch axis, design value (8 to 15 expected)
PITCH_AXIS = 0.30            # chord fraction of the stacking (pitch) axis


# ==========================================================================
# Sections
# ==========================================================================

def naca44_section(t_c: float, te_c: float = 0.0, n_per_side: int = 81,
                   m: float = CAMBER_M, p: float = CAMBER_P) -> np.ndarray:
    """Unit-chord NACA 4-digit section with a blunt trailing edge.

    Thickness uses the closed-TE coefficient (-0.1036 x^4) applied normal to the
    mean line (standard construction). The blunt TE is then added linearly: each
    surface is offset by 0.5 * te_c * (x/c) normal to the chord, which is exactly
    the trailing-edge term of a CST/Kulfan section, so NeuralFoil sees the same
    shape the CAD builds. Returns (2n-1, 2) points in Selig order (upper TE -> LE
    -> lower TE). The upper surface is the convex suction side.
    """
    b = np.linspace(0.0, math.pi, n_per_side)
    x = 0.5 * (1.0 - np.cos(b))
    yt = 5 * t_c * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                    + 0.2843 * x**3 - 0.1036 * x**4)
    yc = np.where(x <= p, m / p**2 * (2 * p * x - x**2),
                  m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))
    dyc = np.where(x <= p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th) + 0.5 * te_c * x
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th) - 0.5 * te_c * x
    upper = np.column_stack([xu[::-1], yu[::-1]])
    lower = np.column_stack([xl[1:], yl[1:]])
    return np.vstack([upper, lower])


def section_area(t_c: float, te_c: float = 0.0) -> float:
    """Area of the unit-chord section (multiply by c^2)."""
    xy = naca44_section(t_c, te_c)
    x, y = xy[:, 0], xy[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def section_max_thickness(t_c: float, te_c: float = 0.0) -> float:
    """Max thickness / chord measured normal to the chord (includes the TE add-on)."""
    xy = naca44_section(t_c, te_c, n_per_side=201)
    n = (len(xy) + 1) // 2
    up, lo = xy[:n][::-1], xy[n - 1:]
    xs = np.linspace(0, 1, 401)
    return float(np.max(np.interp(xs, up[:, 0], up[:, 1]) - np.interp(xs, lo[:, 0], lo[:, 1])))


def naca_code(t_c: float) -> str:
    return f"44{int(round(100 * t_c)):02d}"


@functools.lru_cache(maxsize=8192)
def _kulfan(t_c: float, te_c: float) -> tuple[dict, float, float]:
    """Kulfan parameters of the section, prepared exactly as
    neuralfoil.get_aero_from_airfoil does (normalise, then 8-weight CST fit).
    Returns (params, rotation_deg, scale): alpha_nf = alpha + rotation, Re_nf = Re / scale."""
    import aerosandbox as asb

    af = asb.Airfoil(name="hydroegg-naca44", coordinates=naca44_section(t_c, te_c))
    norm = af.normalize(return_dict=True)
    kaf = norm["airfoil"].to_kulfan_airfoil(n_weights_per_side=8, normalize_coordinates=False)
    kp = {k: np.asarray(v, dtype=float) for k, v in kaf.kulfan_parameters.items()}
    return kp, float(norm["rotation_angle"]), float(norm["scale_factor"])


# ==========================================================================
# Polars
# ==========================================================================

class NeuralFoilPolar:
    """Direct NeuralFoil polars for a set of stations (one t/c per station).

    Call signature shared by all polar classes:
        cl, cd = polar(idx, alpha_deg, re, te_c)
    idx: station index per point, te_c: trailing-edge thickness / chord per point.
    """

    def __init__(self, t_c, n_crit: float = N_CRIT, xtr: float = 1.0,
                 model_size: str = MODEL_SIZE):
        self.t_c = np.atleast_1d(np.asarray(t_c, dtype=float))
        self.n_crit, self.xtr, self.model_size = n_crit, xtr, model_size

    def aero(self, idx, alpha_deg, re, te_c, chunk: int = 1000) -> dict[str, np.ndarray]:
        import neuralfoil as nf

        shape = np.shape(alpha_deg)
        idx, al, re, te = (np.broadcast_to(v, shape).ravel()
                           for v in (idx, alpha_deg, re, te_c))
        keys = np.column_stack([np.round(self.t_c[idx], 4), np.round(np.clip(te, 0, None), 4)])
        uniq, inv = np.unique(keys, axis=0, return_inverse=True)
        inv = inv.ravel()
        fits = [_kulfan(float(a), float(b)) for a, b in uniq]
        up = np.stack([f[0]["upper_weights"] for f in fits])
        lo = np.stack([f[0]["lower_weights"] for f in fits])
        le = np.array([float(f[0]["leading_edge_weight"]) for f in fits])
        tet = np.array([float(f[0]["TE_thickness"]) for f in fits])
        rot = np.array([f[1] for f in fits])
        scl = np.array([f[2] for f in fits])
        out = {k: np.empty(al.size) for k in ("cl", "cd", "confidence", "cp_min")}
        for s in range(0, al.size, chunk):
            sl = slice(s, s + chunk)
            j = inv[sl]
            kp = {"upper_weights": up[j].T, "lower_weights": lo[j].T,
                  "leading_edge_weight": le[j], "TE_thickness": tet[j]}
            r = nf.get_aero_from_kulfan_parameters(
                kp, al[sl] + rot[j], np.clip(re[sl], 1e3, None) / scl[j],
                n_crit=self.n_crit, xtr_upper=self.xtr, xtr_lower=self.xtr,
                model_size=self.model_size)
            ue = np.column_stack([np.abs(r[k]) for k in r if "ue/vinf" in k])
            out["cl"][sl], out["cd"][sl] = r["CL"], r["CD"]
            out["confidence"][sl] = r["analysis_confidence"]
            out["cp_min"][sl] = 1.0 - np.max(ue, axis=1) ** 2
        return {k: v.reshape(shape) for k, v in out.items()}

    def __call__(self, idx, alpha_deg, re, te_c):
        r = self.aero(idx, alpha_deg, re, te_c)
        return r["cl"], r["cd"]


ALPHA_GRID = np.r_[np.arange(-14.0, -2.0, 1.0), np.arange(-2.0, 16.0, 0.5),
                   np.arange(16.0, 32.1, 2.0)]
RE_GRID = 10.0 ** (np.arange(23, 38) / 6.0)        # 6.8e3 .. 1.5e6, 6 points per decade
TE_GRID = np.array([0.0, 0.02, 0.03, 0.045, 0.07, 0.12])


class TablePolar:
    """NeuralFoil polars tabulated per station on (te_c, alpha, log Re).

    Interpolation: linear in te_c and alpha, cubic (Catmull-Rom) in log Re, and Cd is
    interpolated as log(Cd). Measured against direct NeuralFoil on random points in
    the design envelope: rms error about 0.005 in Cl and 1 % in Cd (see the notes).
    Used inside sweeps and optimisation; reported numbers are re-solved with
    ``NeuralFoilPolar`` directly (``solve(..., init=...)``). Outside the grid the
    values are held constant (alpha beyond +32 deg is deep stall anyway).
    """

    def __init__(self, base: NeuralFoilPolar, te_grid=TE_GRID, alpha_grid=ALPHA_GRID,
                 re_grid=RE_GRID):
        self.base = base
        self.te, self.al, self.lre = (np.asarray(te_grid, float), np.asarray(alpha_grid, float),
                                      np.log(np.asarray(re_grid, float)))
        S, T, A, R = np.meshgrid(np.arange(base.t_c.size), np.arange(self.te.size),
                                 np.arange(self.al.size), np.arange(self.lre.size), indexing="ij")
        res = base.aero(S.ravel(), self.al[A.ravel()], np.exp(self.lre[R.ravel()]),
                        self.te[T.ravel()])
        self.cl = res["cl"].reshape(S.shape)
        self.lcd = np.log(res["cd"]).reshape(S.shape)
        self.confidence = res["confidence"].reshape(S.shape)
        self.t_c = base.t_c

    @staticmethod
    def _bracket(grid, v, lo=0, hi_off=2):
        i = np.clip(np.searchsorted(grid, v) - 1, lo, grid.size - hi_off)
        w = np.clip((v - grid[i]) / (grid[i + 1] - grid[i]), 0.0, 1.0)
        return i, w

    def __call__(self, idx, alpha_deg, re, te_c):
        shape = np.broadcast_shapes(np.shape(idx), np.shape(alpha_deg), np.shape(re), np.shape(te_c))
        idx, al, re, te = (np.broadcast_to(v, shape) for v in (idx, alpha_deg, re, te_c))
        it, wt = self._bracket(self.te, te)
        ia, wa = self._bracket(self.al, al)
        ir, t = self._bracket(self.lre, np.log(np.clip(re, 1.0, None)), lo=1, hi_off=3)
        wr = ((-t**3 + 2 * t**2 - t) / 2, (3 * t**3 - 5 * t**2 + 2) / 2,
              (-3 * t**3 + 4 * t**2 + t) / 2, (t**3 - t**2) / 2)
        cl = np.zeros(shape)
        lcd = np.zeros(shape)
        for dt, ft in ((0, 1 - wt), (1, wt)):
            for da, fa in ((0, 1 - wa), (1, wa)):
                for off, fr in zip((-1, 0, 1, 2), wr):
                    f = ft * fa * fr
                    cl += f * self.cl[idx, it + dt, ia + da, ir + off]
                    lcd += f * self.lcd[idx, it + dt, ia + da, ir + off]
        return cl, np.exp(lcd)


class IdealPolar:
    """Thin-airfoil polar for solver verification: Cl = 2 pi (alpha - alpha0), Cd const."""

    def __init__(self, alpha0_deg: float = 0.0, cd: float = 0.0):
        self.alpha0, self.cd0 = alpha0_deg, cd

    def __call__(self, idx, alpha_deg, re, te_c):
        cl = 2 * math.pi * np.radians(np.asarray(alpha_deg) - self.alpha0)
        return cl, np.full(np.shape(cl), self.cd0)


# ==========================================================================
# Rotor and solver
# ==========================================================================

@dataclass
class Rotor:
    """Blade discretised into annular elements. All lengths in mm, angles in deg."""

    r_mm: np.ndarray              # element centres, radius perpendicular to the axis
    dr_mm: np.ndarray             # element widths, perpendicular to the axis
    chord_mm: np.ndarray
    twist_deg: np.ndarray
    t_c: np.ndarray
    te_mm: float = P.MIN_TE       # absolute trailing-edge thickness (0 = sharp)
    cone_deg: float | np.ndarray = 0.0   # local aft cone of the pitch axis
    x_aft_mm: np.ndarray | None = None   # pitch-axis aft offset from the hinge; None = straight
    r_tip_mm: float = P.ROTOR_R_TIP_TARGET
    r_hub_mm: float = R_AERO_ROOT        # hub-loss radius = start of the aerodynamic span
    r_hinge_mm: float = P.HINGE_R
    n_blades: int = P.N_BLADES
    tip_loss: bool = True
    hub_loss: bool = True
    station: np.ndarray | None = None    # polar station index per element; None = 0..n-1

    @property
    def cone(self) -> np.ndarray:
        return np.broadcast_to(np.asarray(self.cone_deg, float), self.r_mm.shape)

    @property
    def x_aft(self) -> np.ndarray:
        if self.x_aft_mm is not None:
            return np.asarray(self.x_aft_mm, float)
        return (self.r_mm - self.r_hinge_mm) * np.tan(np.radians(self.cone))

    @property
    def te_c(self) -> np.ndarray:
        return np.broadcast_to(self.te_mm, self.r_mm.shape) / self.chord_mm


def elements(n: int = 24, r_root: float = R_AERO_ROOT,
             r_tip: float = P.ROTOR_R_TIP_TARGET) -> tuple[np.ndarray, np.ndarray]:
    """Element centres and widths (mm), cosine-clustered toward root and tip."""
    s = 0.5 * (1 - np.cos(np.linspace(0, math.pi, n + 1)))
    edges = r_root + (r_tip - r_root) * s
    return 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)


A_FLOOR = -5.0     # axial induction below this is treated as a non-physical solution


def _axial_induction(k: np.ndarray, F: np.ndarray) -> np.ndarray:
    """a from k = sigma Cn cos^2(beta) / (4 F sin^2 phi): momentum for a <= 0.4,
    Buhl's empirical turbulent-wake branch above (Ning 2014 closed form)."""
    k = np.maximum(k, -0.9)
    a_mom = k / (1 + k)
    g1 = 2 * F * k - (10 / 9 - F)
    g2 = np.maximum(2 * F * k - F * (4 / 3 - F), 1e-12)
    g3 = 2 * F * k - (25 / 9 - 2 * F)
    small = np.abs(g3) < 1e-6
    a_buhl = np.where(small, 1 - 1 / (2 * np.sqrt(g2)),
                      (g1 - np.sqrt(g2)) / np.where(small, 1.0, g3))
    return np.where(k <= 2 / 3, a_mom, a_buhl)


def solve(rotor: Rotor, polar, u_m_s, tsr, *, relax: float = 0.35, tol: float = 1e-7,
          max_iter: int = 400, method: str = "auto", init: dict | None = None
          ) -> dict[str, np.ndarray]:
    """Solve the BEM equations for one or more operating points.

    u_m_s, tsr: scalars or 1-D arrays, broadcast together to n_op operating points.
    TSR is based on the perpendicular tip radius. method: "auto" (relaxed fixed point,
    Ning bracketing fallback), or "ning" (bracketing for every element).

    Returns element arrays of shape (n_op, n_el) and integrated loads of shape (n_op,):
    cp, ct, power_W, torque_Nm, thrust_N, rpm, m_flap_Nm and m_edge_Nm (per blade,
    about the hinge pin; flap positive = pushes the blade aft onto its deployed stop),
    plus converged flags and the method used per element (0 fixed point, 1 fallback,
    -1 not converged). ``init`` (a previous result with the same shape) seeds a and
    a', so a cheap table-polar solve can be refined with direct NeuralFoil polars.
    """
    u, tsr = np.broadcast_arrays(np.atleast_1d(np.asarray(u_m_s, float)),
                                 np.atleast_1d(np.asarray(tsr, float)))
    n_op, n_el = u.size, rotor.r_mm.size
    B = rotor.n_blades
    R, Rh = rotor.r_tip_mm * 1e-3, rotor.r_hub_mm * 1e-3
    r = np.broadcast_to(rotor.r_mm * 1e-3, (n_op, n_el))
    c = np.broadcast_to(rotor.chord_mm * 1e-3, (n_op, n_el))
    th = np.broadcast_to(np.radians(rotor.twist_deg), (n_op, n_el))
    beta = np.broadcast_to(np.radians(rotor.cone), (n_op, n_el))
    cb = np.cos(beta)
    U = np.broadcast_to(u[:, None], (n_op, n_el))
    omega = np.broadcast_to((tsr * u / R)[:, None], (n_op, n_el))
    sigma = B * c / (2 * math.pi * r)
    st = np.arange(n_el) if rotor.station is None else np.asarray(rotor.station)
    idx = np.broadcast_to(st, (n_op, n_el))
    te_c = np.broadcast_to(rotor.te_c, (n_op, n_el))
    lam = omega * r / (U * cb)                       # effective local speed ratio

    def loss(phi, rs):
        """Prandtl tip x hub loss at radius rs (m)."""
        sphi = np.maximum(np.abs(np.sin(phi)), 1e-6)
        F = np.ones_like(phi)
        if rotor.tip_loss:
            F = F * 2 / math.pi * np.arccos(np.clip(np.exp(-B / 2 * (R - rs) / (rs * sphi)), 0, 1))
        if rotor.hub_loss:
            F = F * 2 / math.pi * np.arccos(np.clip(np.exp(-B / 2 * (rs - Rh) / (Rh * sphi)), 0, 1))
        return np.clip(F, 1e-4, 1.0)

    def update(phi, W, sel=Ellipsis):
        """New (a, kp, F) given phi and W (for Re); sel picks a subset of elements."""
        alpha = np.degrees(phi - th[sel])
        re = W * c[sel] / NU
        cl, cd = polar(idx[sel], alpha, re, te_c[sel])
        sp, cp_ = np.sin(phi), np.cos(phi)
        cn, ct = cl * cp_ + cd * sp, cl * sp - cd * cp_
        F = loss(phi, r[sel])
        k = sigma[sel] * cn * cb[sel] ** 2 / (4 * F * np.maximum(sp**2, 1e-12))
        spc = np.where(np.abs(sp * cp_) < 1e-12, 1e-12, sp * cp_)
        kp = np.clip(sigma[sel] * ct / (4 * F * spc), -0.9, 0.95)
        return _axial_induction(k, F), kp, k

    # ---- relaxed fixed point on (a, a') -------------------------------------
    a = np.full((n_op, n_el), 0.3) if init is None else np.array(init["a"], float)
    ap = np.zeros((n_op, n_el)) if init is None else np.array(init["ap"], float)
    done = np.zeros((n_op, n_el), bool)
    if method == "auto":
        # Annuli are independent in BEM, so converged elements are frozen and only
        # the active set is iterated (exact, and much cheaper with NeuralFoil).
        act = np.ones((n_op, n_el), bool)
        for _ in range(max_iter):
            sel = np.nonzero(act)
            vn, vt = U[sel] * cb[sel] * (1 - a[sel]), omega[sel] * r[sel] * (1 + ap[sel])
            a_new, kp, _k = update(np.arctan2(vn, vt), np.hypot(vn, vt), sel)
            ap_new = kp / (1 - kp)
            da, dap = a_new - a[sel], ap_new - ap[sel]
            conv = (np.abs(da) < tol) & (np.abs(dap) < tol)
            a[sel] = a[sel] + relax * da
            ap[sel] = ap[sel] + relax * dap
            done[sel] = conv
            act[sel] = ~conv
            if not act.any():
                break
        done &= a > A_FLOOR                     # a fixed point on the k clamp is not physical
    used = np.where(done, 0, -1)

    # ---- Ning bracketing fallback in phi --------------------------------------
    bad = ~done
    if bad.any():
        sel = np.nonzero(bad)
        lam_s, U_s, cb_s, om_s, r_s = lam[sel], U[sel], cb[sel], omega[sel], r[sel]
        W_s = np.hypot(U_s * cb_s * (1 - a[sel]), om_s * r_s * (1 + ap[sel]))

        def resid(phi):
            """Ning (2014) residual. In the momentum region sin(phi)/(1-a) is written
            as sin(phi)(1+k), which is smooth through k = -1; dividing by 1-a with a
            clamped k created false roots at a = -9 (found at the hub-loss element)."""
            nonlocal W_s
            for _ in range(2):                  # settle Re (W) at this phi
                a_s, kp_s, k_s = update(phi, W_s, sel)
                W_s = np.hypot(U_s * cb_s * (1 - a_s), om_s * r_s / (1 - kp_s))
            sp = np.sin(phi)
            lhs = np.where(k_s <= 2 / 3, sp * (1 + k_s), sp / np.maximum(1 - a_s, 1e-9))
            f = lhs - np.cos(phi) * (1 - kp_s) / lam_s
            return f, a_s, kp_s

        eps = 1e-6
        lo = np.full(lam_s.shape, eps)
        hi = np.full(lam_s.shape, math.pi / 2)
        f_lo, *_ = resid(lo)
        f_hi, *_ = resid(hi)
        alt = f_lo * f_hi > 0                    # try the other branch (phi > 90 deg)
        lo = np.where(alt, math.pi / 2, lo)
        hi = np.where(alt, math.pi - eps, hi)
        f_lo2, *_ = resid(lo)
        f_hi2, *_ = resid(hi)
        ok = f_lo2 * f_hi2 <= 0
        f_lo = f_lo2
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f_mid, a_m, kp_m = resid(mid)
            left = f_lo * f_mid <= 0
            hi = np.where(left, mid, hi)
            lo = np.where(left, lo, mid)
            f_lo = np.where(left, f_lo, f_mid)
        phi_b = 0.5 * (lo + hi)
        f_fin, a_m, kp_m = resid(phi_b)
        conv = ok & (np.abs(f_fin) < 1e-5) & (a_m > A_FLOOR)   # a root, not a jump; physical a
        # bracketed elements keep the bisection state even when flagged (best estimate)
        a[sel] = np.where(ok, a_m, a[sel])
        ap[sel] = np.where(ok, kp_m / (1 - kp_m), ap[sel])
        used[sel] = np.where(conv, 1, -1)
        phi_fix = np.full((n_op, n_el), np.nan)
        phi_fix[sel] = np.where(ok, phi_b, np.nan)

    # ---- consolidated state and loads -----------------------------------------
    vn, vt = U * cb * (1 - a), omega * r * (1 + ap)
    phi, W = np.arctan2(vn, vt), np.hypot(vn, vt)
    if bad.any():                               # identical at a true root; consistent if flagged
        phi = np.where(np.isnan(phi_fix), phi, phi_fix)
    alpha = np.degrees(phi - th)
    re = W * c / NU
    cl, cd = polar(idx, alpha, re, te_c)
    F = loss(phi, r)
    cn = cl * np.cos(phi) + cd * np.sin(phi)
    ct_ = cl * np.sin(phi) - cd * np.cos(phi)
    q = 0.5 * RHO * W**2 * c
    ds = np.broadcast_to(rotor.dr_mm * 1e-3, (n_op, n_el)) / cb    # span length of element
    n_el_force, t_el_force = q * cn * ds, q * ct_ * ds              # per blade, N
    thrust = B * np.sum(n_el_force * cb, axis=1)
    torque = B * np.sum(t_el_force * r, axis=1)
    om = tsr * u / R
    power = om * torque
    area = math.pi * R**2
    dx = np.broadcast_to(rotor.x_aft * 1e-3, (n_op, n_el))
    drh = r - rotor.r_hinge_mm * 1e-3
    m_flap = np.sum(n_el_force * (dx * np.sin(beta) + drh * cb), axis=1)
    m_edge = np.sum(t_el_force * np.hypot(dx, drh), axis=1)
    return {
        "u_m_s": u, "tsr": tsr, "rpm": om * 60 / (2 * math.pi), "omega_rad_s": om,
        "cp": power / (0.5 * RHO * u**3 * area), "ct": thrust / (0.5 * RHO * u**2 * area),
        "power_W": power, "torque_Nm": torque, "thrust_N": thrust,
        "m_flap_Nm": m_flap, "m_edge_Nm": m_edge,
        "a": a, "ap": ap, "phi_deg": np.degrees(phi), "alpha_deg": alpha, "re": re,
        "cl": cl, "cd": cd, "F": F, "W_m_s": W,
        "fn_N_m": q * cn, "ft_N_m": q * ct_,
        "method": used, "converged": used >= 0,
        "n_unconverged": np.sum(used < 0, axis=1), "n_fallback": np.sum(used == 1, axis=1),
    }


# ==========================================================================
# Design study
# ==========================================================================

DESIGN_U = P.DESIGN_FLOW_V                      # m/s
U_EVAL = (1.0, 1.5, 2.0, 2.5)                   # m/s
TSR_CANDIDATES = (3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
TSR_TIE = 0.005             # relative Cp band treated as a tie when choosing the TSR
TSR_SWEEP = np.round(np.arange(1.5, 11.01, 0.25), 2)
N_EL = 24
N_STATIONS = 21
R_TIP = P.ROTOR_R_TIP_TARGET
TE_GRID_NARROW = np.array([0.015, 0.025, 0.035, 0.05])   # covers 0.6 mm on 14-28 mm chords

# Estimates used only for the structural and cavitation side checks.
PLA_RHO = 1240.0            # kg/m^3, solid PLA, assumed
P_ATM = 101325.0            # Pa
P_VAP = 2339.0              # Pa, water at 20 C
DEPTH_M = 0.10              # m, assumed immersion of the blade tip at the top of its circle


def thickness_ratio(r_mm, r_root: float = R_AERO_ROOT, r_tip: float = R_TIP):
    """Base NACA t/c, linear from TC_ROOT at the aero root to TC_TIP at the tip."""
    s = (np.asarray(r_mm, float) - r_root) / (r_tip - r_root)
    return TC_ROOT + (TC_TIP - TC_ROOT) * np.clip(s, 0, 1)


def linear(r_mm, y_root: float, y_tip: float, r_root: float = R_AERO_ROOT, r_tip: float = R_TIP):
    return y_root + (y_tip - y_root) * (np.asarray(r_mm, float) - r_root) / (r_tip - r_root)


def hull_following_axis(r_mm) -> tuple[np.ndarray, np.ndarray]:
    """Estimate of a pitch axis that copies the hull when folded 90 deg forward about
    the hinge: the blade point at deployed radius r sits at x = HINGE_X - (r - HINGE_R)
    when stowed, at the hull radius there. Returns (x_aft_mm, local cone_deg)."""
    from . import profile

    xs = P.HINGE_X - (np.asarray(r_mm, float) - P.HINGE_R)
    x_aft = np.asarray(profile.radius(xs)) - profile.radius(P.HINGE_X)
    cone = np.degrees(np.arctan([-profile.slope(float(x)) for x in xs]))
    return x_aft, np.asarray(cone)


def verify_glauert(tsr: float = 4.0, n: int = 40) -> tuple[float, float]:
    """Solver check: Glauert's optimum rotor (wake rotation, no drag, no tip loss)
    over the aero span, flown with an ideal polar. Returns (Cp_bem, Cp_theory)."""
    from scipy.integrate import quad
    from scipy.optimize import brentq

    def a_of(x):
        return brentq(lambda a: 16 * a**3 - 24 * a**2 + a * (9 - 3 * x**2) - 1 + x**2,
                      0.25 + 1e-12, 1 / 3)

    def integrand(x):
        a = a_of(x)
        return (1 - 3 * a) / (4 * a - 1) * (1 - a) * x**3

    cp_th = 8 / tsr**2 * quad(integrand, tsr * R_AERO_ROOT / R_TIP, tsr, limit=200)[0]
    r, dr = elements(n)
    x = tsr * r / R_TIP
    a = np.array([a_of(v) for v in x])
    ap = (1 - 3 * a) / (4 * a - 1)
    phi = np.arctan((1 - a) / (x * (1 + ap)))
    alpha_d = 5.0
    cl = 2 * math.pi * math.radians(alpha_d)
    sigma = 4 * a * np.sin(phi) ** 2 / ((1 - a) * cl * np.cos(phi))
    rot = Rotor(r, dr, sigma * 2 * math.pi * r / P.N_BLADES, np.degrees(phi) - alpha_d,
                np.full(n, 0.15), te_mm=0.0, tip_loss=False, hub_loss=False)
    return float(solve(rot, IdealPolar(), DESIGN_U, tsr)["cp"][0]), float(cp_th)


class Study:
    """Element discretisation + polars shared by every design evaluated here."""

    def __init__(self, n_el: int = N_EL, n_crit: float = N_CRIT, xtr: float = 1.0,
                 te_grid=TE_GRID):
        self.r, self.dr = elements(n_el)
        self.t_c = thickness_ratio(self.r)
        self.direct = NeuralFoilPolar(self.t_c, n_crit=n_crit, xtr=xtr)
        self.table = TablePolar(self.direct, te_grid=te_grid)

    def rotor(self, chord, twist, te_mm: float = P.MIN_TE, cone=0.0, **kw) -> Rotor:
        return Rotor(self.r, self.dr, np.asarray(chord, float), np.asarray(twist, float),
                     self.t_c, te_mm=te_mm, cone_deg=cone, **kw)

    def linear_rotor(self, x, cone=CONE_NOMINAL, **kw) -> Rotor:
        """x = (chord_root_mm, chord_tip_mm, twist_root_deg, twist_tip_deg)."""
        return self.rotor(linear(self.r, x[0], x[1]), linear(self.r, x[2], x[3]), cone=cone, **kw)

    def perf(self, rot: Rotor, u, tsr, direct: bool = True, **kw) -> dict:
        res = solve(rot, self.table, u, tsr, **kw)
        return solve(rot, self.direct, u, tsr, init=res, **kw) if direct else res

    def cp_max(self, rot: Rotor, u: float = DESIGN_U, direct: bool = True,
               tsr_grid=TSR_SWEEP) -> tuple[float, dict]:
        """(tsr_opt, result at tsr_opt): table sweep, parabolic refinement, then a
        direct-polar solve at the refined TSR."""
        sw = solve(rot, self.table, u, tsr_grid)
        j = int(np.clip(np.argmax(sw["cp"]), 1, len(tsr_grid) - 2))
        x, y = tsr_grid[j - 1:j + 2], sw["cp"][j - 1:j + 2]
        den = y[0] - 2 * y[1] + y[2]
        t_opt = x[1] + 0.5 * (x[1] - x[0]) * (y[0] - y[2]) / den if den < 0 else x[1]
        t_opt = float(np.clip(t_opt, x[0], x[2]))
        return t_opt, self.perf(rot, u, t_opt, direct=direct)

    def runaway(self, rot: Rotor, u: float, sweep: dict) -> dict:
        """Zero-torque (free-wheeling) point above the Cp peak, by secant on torque."""
        q, t = sweep["torque_Nm"], sweep["tsr"]
        j0 = int(np.argmax(sweep["cp"]))
        if not np.any(q[j0:] < 0):
            raise ValueError("no zero-torque crossing in the TSR sweep; extend TSR_SWEEP")
        k = j0 + int(np.argmax(q[j0:] < 0))
        t0, t1 = t[k - 1], t[k]
        for _ in range(3):
            r0, r1 = self.perf(rot, u, t0), self.perf(rot, u, t1)
            q0, q1 = r0["torque_Nm"][0], r1["torque_Nm"][0]
            t0, t1 = t1, t1 - q1 * (t1 - t0) / (q1 - q0)
        return self.perf(rot, u, t1)

    def schmitz(self, tsr: float, te_mm: float = 0.0, u: float = DESIGN_U) -> dict:
        """Schmitz optimum, phi = 2/3 atan(1/x), c = 16 pi r sin^2(phi/2) / (B Cl).

        The angle of attack maximises L/D at the Re the chord itself produces: since
        c = K / Cl, a lower Cl means a longer chord and a higher Re. At these low Re
        that coupling moves the best angle well below the fixed-Re best-L/D angle."""
        n = self.r.size
        x = tsr * self.r / R_TIP
        phi = (2 / 3) * np.arctan(1 / x)
        K = 16 * math.pi * self.r / P.N_BLADES * np.sin(phi / 2) ** 2   # chord_mm * Cl
        W = u * np.sqrt((2 / 3) ** 2 + x**2)                              # a ~ 1/3 estimate
        al = np.arange(-2.0, 14.01, 0.25)
        na = al.size
        I, AL = np.repeat(np.arange(n), na), np.tile(al, n)
        Kp, Wp = np.repeat(K, na), np.repeat(W, na)
        cl = np.ones(n * na)
        for _ in range(4):
            c = Kp / np.maximum(cl, 0.2)
            cl, cd = self.direct(I, AL, Wp * c * 1e-3 / NU, te_mm / c)
        ld = np.where(cl > 0.2, cl / cd, -np.inf).reshape(n, na)
        j = np.argmax(ld, axis=1)
        cl_o = cl.reshape(n, na)[np.arange(n), j]
        return {"chord": K / cl_o, "twist": np.degrees(phi) - al[j], "alpha": al[j],
                "cl": cl_o, "ld": ld[np.arange(n), j]}

    def _annulus_search(self, C, T, te_mm, cone, u, tsr):
        n, k = C.shape
        cone_el = np.broadcast_to(np.asarray(cone, float), (n,))
        rot = Rotor(np.repeat(self.r, k), np.repeat(self.dr, k), C.reshape(-1), T.reshape(-1),
                    np.repeat(self.t_c, k), te_mm=te_mm, cone_deg=np.repeat(cone_el, k),
                    station=np.repeat(np.arange(n), k))
        res = solve(rot, self.table, u, tsr, max_iter=300)
        ok = ((res["method"][0] == 0) & (res["a"][0] >= 0) & (res["a"][0] <= 0.7)
              & (res["alpha_deg"][0] > -4) & (res["alpha_deg"][0] < 14))
        ft = np.where(ok, res["ft_N_m"][0], -np.inf).reshape(n, k)
        j = np.argmax(ft, axis=1)
        return C[np.arange(n), j], T[np.arange(n), j]

    def bem_optimal(self, tsr: float, te_mm: float = 0.0, cone=0.0, u: float = DESIGN_U) -> dict:
        """BEM-optimal chord and twist. Annuli are independent in BEM, so the rotor
        optimum is every element maximising its own tangential force at (u, tsr).
        Search: coarse grid in (log chord, twist), then a fine grid around the best
        cell. Only fixed-point solutions with 0 <= a <= 0.7 and -4 < alpha < 14 deg
        are admitted (no turbulent-wake or post-stall branches)."""
        n = self.r.size
        cg, tg = np.meshgrid(np.geomspace(4.0, 60.0, 24), np.arange(-6.0, 36.1, 1.0), indexing="ij")
        c0, t0 = self._annulus_search(np.tile(cg.ravel(), (n, 1)), np.tile(tg.ravel(), (n, 1)),
                                      te_mm, cone, u, tsr)
        cf = c0[:, None, None] * np.geomspace(0.89, 1.12, 9)[None, :, None]
        tf = t0[:, None, None] + np.linspace(-1.0, 1.0, 21)[None, None, :]
        C = np.broadcast_to(cf, (n, 9, 21)).reshape(n, -1)
        T = np.broadcast_to(tf, (n, 9, 21)).reshape(n, -1)
        c1, t1 = self._annulus_search(C, T, te_mm, cone, u, tsr)
        return {"chord": c1, "twist": t1}

    def linearise(self, y) -> tuple[float, float]:
        """Root and tip values of the line through y(0.7 R) and y(0.9 R)."""
        y7, y9 = np.interp([0.7 * R_TIP, 0.9 * R_TIP], self.r, y)
        slope = (y9 - y7) / (0.2 * R_TIP)
        return (float(y7 + slope * (R_AERO_ROOT - 0.7 * R_TIP)),
                float(y7 + slope * (R_TIP - 0.7 * R_TIP)))

    def optimise_linear(self, x0, tsr: float, cone=CONE_NOMINAL, u: float = DESIGN_U,
                        maxfev: int = 300) -> tuple[np.ndarray, float]:
        """Max Cp(tsr) over linear chord and twist, x = (c_root, c_tip, tw_root, tw_tip),
        inside the printability and folding bounds. Nelder-Mead on table polars."""
        from scipy.optimize import minimize

        lo = np.array([P.MIN_CHORD, P.MIN_CHORD, -5.0, -10.0])
        hi = np.array([MAX_CHORD, MAX_CHORD, P.MAX_ROOT_TWIST, 15.0])
        x0 = np.clip(np.asarray(x0, float), lo, hi)
        steps = np.array([2.0, 1.0, 2.0, 1.0])
        simplex = [x0]
        for i in range(4):
            v = x0.copy()
            v[i] = x0[i] - steps[i] if x0[i] - steps[i] >= lo[i] else x0[i] + steps[i]
            simplex.append(v)

        def f(x):
            return -float(solve(self.linear_rotor(x, cone=cone), self.table, u, tsr)["cp"][0])

        res = minimize(f, x0, method="Nelder-Mead", bounds=list(zip(lo, hi)),
                       options={"initial_simplex": np.array(simplex), "xatol": 0.02,
                                "fatol": 1e-6, "maxfev": maxfev})
        return res.x, -res.fun


def _ops(res: dict, j: int = 0) -> dict:
    """Scalar operating-point summary with unit-suffixed keys."""
    return {
        "flow_m_s": round(float(res["u_m_s"][j]), 3), "tsr": round(float(res["tsr"][j]), 3),
        "rpm": round(float(res["rpm"][j]), 1), "cp": round(float(res["cp"][j]), 4),
        "ct": round(float(res["ct"][j]), 4), "power_W": round(float(res["power_W"][j]), 3),
        "torque_Nm": round(float(res["torque_Nm"][j]), 4),
        "thrust_N": round(float(res["thrust_N"][j]), 3),
        "flap_moment_hinge_per_blade_Nm": round(float(res["m_flap_Nm"][j]), 4),
        "edge_moment_hinge_per_blade_Nm": round(float(res["m_edge_Nm"][j]), 4),
        "unconverged_elements": int(res["n_unconverged"][j]),
    }


def run_study(log=print) -> dict:
    """Full design study. Returns every number written to the JSON and figures."""
    out: dict = {}
    st = Study()
    out["study"] = st
    cp_g, cp_t = verify_glauert(4.0)
    out["glauert"] = {"cp_bem": cp_g, "cp_theory": cp_t}
    log(f"Glauert check TSR 4: BEM {cp_g:.5f} vs theory {cp_t:.5f}")

    # ---- 1. design TSR -------------------------------------------------------
    rows, designs = [], {}
    for tsr in TSR_CANDIDATES:
        opt = st.bem_optimal(tsr)
        r_opt = st.rotor(opt["chord"], opt["twist"], te_mm=0.0)
        sch = st.schmitz(tsr)
        r_sch = st.rotor(sch["chord"], sch["twist"], te_mm=0.0)
        x0 = [*st.linearise(opt["chord"]), *st.linearise(opt["twist"])]
        x, _ = st.optimise_linear(x0, tsr)
        r_fin = st.linear_rotor(x)
        p_opt, p_sch, p_fin = (st.perf(r, DESIGN_U, tsr) for r in (r_opt, r_sch, r_fin))
        t_pk, p_pk = st.cp_max(r_fin)
        rows.append({
            "design_tsr": tsr, "cp_bem_optimum": float(p_opt["cp"][0]),
            "cp_schmitz": float(p_sch["cp"][0]), "cp_constrained": float(p_fin["cp"][0]),
            "cp_constrained_peak": float(p_pk["cp"][0]), "tsr_constrained_peak": t_pk,
            "schmitz_tip_chord_mm": float(sch["chord"][-1]), "schmitz_root_chord_mm": float(sch["chord"][0]),
            "schmitz_root_twist_deg": float(sch["twist"][0]),
            "optimum_chord_0p9R_mm": float(np.interp(0.9 * R_TIP, st.r, opt["chord"])),
            "constrained_x": [float(v) for v in x],
            "rpm_at_1p5_m_s": float(p_fin["rpm"][0]), "torque_at_1p5_m_s_Nm": float(p_fin["torque_Nm"][0]),
        })
        designs[tsr] = {"opt": opt, "sch": sch, "x": x}
        log(f"TSR {tsr}: optimum {rows[-1]['cp_bem_optimum']:.4f} schmitz {rows[-1]['cp_schmitz']:.4f} "
            f"constrained {rows[-1]['cp_constrained']:.4f} (peak {rows[-1]['cp_constrained_peak']:.4f} "
            f"at {t_pk:.2f}) x={np.round(x, 2)}")
    best = max(r["cp_constrained"] for r in rows)
    tied = [r for r in rows if r["cp_constrained"] >= (1 - TSR_TIE) * best]
    tsr_d = max(r["design_tsr"] for r in tied)
    out["tsr_rows"], out["tsr_d"] = rows, tsr_d
    log(f"design TSR -> {tsr_d}")

    # ---- 2. constraint chain at the design TSR --------------------------------
    opt = designs[tsr_d]["opt"]
    opt_te = st.bem_optimal(tsr_d, te_mm=P.MIN_TE)
    cl_r, cl_t = st.linearise(opt_te["chord"])
    tw_r, tw_t = st.linearise(opt_te["twist"])
    c_lin, t_lin = linear(st.r, cl_r, cl_t), linear(st.r, tw_r, tw_t)
    c_clip = np.clip(c_lin, P.MIN_CHORD, MAX_CHORD)
    t_cap = np.minimum(t_lin, P.MAX_ROOT_TWIST)
    x_fin, _ = st.optimise_linear([cl_r, cl_t, tw_r, tw_t], tsr_d)
    x_fin = np.array([round(x_fin[0], 1), round(x_fin[1], 1), round(x_fin[2], 1), round(x_fin[3], 1)])
    chain_defs = [
        ("optimum", "BEM optimum (sharp TE, no cone)", opt["chord"], opt["twist"], 0.0, 0.0),
        ("te", "+ blunt TE 0.6 mm (re-optimised)", opt_te["chord"], opt_te["twist"], P.MIN_TE, 0.0),
        ("linear", "+ linear chord and twist (0.7R-0.9R fit)", c_lin, t_lin, P.MIN_TE, 0.0),
        ("chord", f"+ chord in [{P.MIN_CHORD:g}, {MAX_CHORD:g}] mm (clipped)", c_clip, t_lin, P.MIN_TE, 0.0),
        ("twist", f"+ twist <= {P.MAX_ROOT_TWIST:g} deg (clipped)", c_clip, t_cap, P.MIN_TE, 0.0),
        ("cone", f"+ {CONE_NOMINAL:g} deg aft cone", c_clip, t_cap, P.MIN_TE, CONE_NOMINAL),
        ("final", "final: linear re-optimised in all bounds", linear(st.r, x_fin[0], x_fin[1]),
         linear(st.r, x_fin[2], x_fin[3]), P.MIN_TE, CONE_NOMINAL),
    ]
    iso_defs = [
        ("iso_chord", "chord clip alone", np.clip(opt_te["chord"], P.MIN_CHORD, MAX_CHORD), opt_te["twist"], P.MIN_TE, 0.0),
        ("iso_twist", "twist cap alone", opt_te["chord"], np.minimum(opt_te["twist"], P.MAX_ROOT_TWIST), P.MIN_TE, 0.0),
        ("iso_cone", "cone alone", opt_te["chord"], opt_te["twist"], P.MIN_TE, CONE_NOMINAL),
    ]
    sch = designs[tsr_d]["sch"]
    ref_defs = [("schmitz", "Schmitz (sharp TE, no cone)", sch["chord"], sch["twist"], 0.0, 0.0)]
    chain = []
    for key, label, c, t, te, cone in chain_defs + iso_defs + ref_defs:
        rot = st.rotor(c, t, te_mm=te, cone=cone)
        p_d = st.perf(rot, DESIGN_U, tsr_d)
        t_pk, p_pk = st.cp_max(rot)
        chain.append({"key": key, "label": label, "chord": np.asarray(c), "twist": np.asarray(t),
                      "te_mm": te, "cone_deg": cone, "cp_design_tsr": float(p_d["cp"][0]),
                      "cp_peak": float(p_pk["cp"][0]), "tsr_peak": t_pk,
                      "chord_root_mm": float(np.interp(R_AERO_ROOT, st.r, c)),
                      "chord_tip_mm": float(np.interp(R_TIP, st.r, c)),
                      "twist_root_deg": float(np.interp(R_AERO_ROOT, st.r, t)),
                      "twist_tip_deg": float(np.interp(R_TIP, st.r, t))})
        log(f"  {label:48s} Cp(TSR {tsr_d}) {chain[-1]['cp_design_tsr']:.4f}  peak {chain[-1]['cp_peak']:.4f} at {t_pk:.2f}")
    out["chain"], out["x_final"] = chain, x_fin
    out["schmitz"], out["optimum"], out["optimum_te"] = sch, opt, opt_te

    # ---- 3. final design performance -------------------------------------------
    rot = st.linear_rotor(x_fin)
    out["rotor"] = rot
    perf = {"sweep": {}, "mppt": [], "design_tsr": [], "runaway": [], "envelope": []}
    n_unconv = n_fallback = 0
    for u in U_EVAL:
        sw = st.perf(rot, u, TSR_SWEEP)
        t_pk, p_pk = st.cp_max(rot, u)
        p_d = st.perf(rot, u, tsr_d)
        p_run = st.runaway(rot, u, sw)
        perf["sweep"][u] = sw
        perf["mppt"].append(p_pk)
        perf["design_tsr"].append(p_d)
        perf["runaway"].append(p_run)
        valid = sw["tsr"] <= p_run["tsr"][0] + 1e-9
        perf["envelope"].append({
            "flow_m_s": u,
            "max_thrust_N": float(max(np.max(sw["thrust_N"][valid]), p_run["thrust_N"][0])),
            "max_flap_moment_hinge_per_blade_Nm": float(max(np.max(sw["m_flap_Nm"][valid]), p_run["m_flap_Nm"][0])),
            "max_torque_Nm": float(np.max(sw["torque_Nm"][valid])),
        })
        for r_ in (sw, p_pk, p_d, p_run):
            n_unconv += int(np.sum(r_["n_unconverged"]))
            n_fallback += int(np.sum(r_["n_fallback"]))
        log(f"  U {u}: MPPT TSR {t_pk:.2f} Cp {p_pk['cp'][0]:.4f} P {p_pk['power_W'][0]:.2f} W "
            f"rpm {p_pk['rpm'][0]:.0f} | runaway TSR {p_run['tsr'][0]:.2f} rpm {p_run['rpm'][0]:.0f} "
            f"T {p_run['thrust_N'][0]:.1f} N")
    perf["n_unconverged"], perf["n_fallback"] = n_unconv, n_fallback
    out["perf"] = perf
    i15 = U_EVAL.index(DESIGN_U)
    dp = perf["mppt"][i15]

    # ---- 4. checks -------------------------------------------------------------
    chk: dict = {}
    tab_only = solve(rot, st.table, DESIGN_U, dp["tsr"][0])
    ning = solve(rot, st.direct, DESIGN_U, dp["tsr"][0], method="ning", init=tab_only)
    chk["cp_table_vs_direct"] = (float(tab_only["cp"][0]), float(dp["cp"][0]))
    chk["cp_fixed_point_vs_ning"] = (float(dp["cp"][0]), float(ning["cp"][0]))
    chk["max_abs_a_diff_fixed_point_vs_ning"] = float(np.max(np.abs(ning["a"] - dp["a"])))
    all_cp = np.concatenate([perf["sweep"][u]["cp"] for u in U_EVAL])
    chk["max_cp_any_sweep"] = float(np.max(all_cp))
    mid = dp["F"][0] > 0.9
    chk["a_midspan_range"] = (float(dp["a"][0][mid].min()), float(dp["a"][0][mid].max()))
    chk["a_max_design"] = float(dp["a"][0].max())
    chk["a_max_design_r_mm"] = float(st.r[np.argmax(dp["a"][0])])
    chk["n_elements_buhl_design"] = int(np.sum(dp["a"][0] > 0.4))
    chk["ap_max_design"] = float(dp["ap"][0].max())
    aero = st.direct.aero(np.arange(st.r.size), dp["alpha_deg"][0], dp["re"][0], rot.te_c)
    chk["nf_confidence_min_design"] = float(aero["confidence"].min())
    al = np.arange(0.0, 20.01, 0.25)
    n = st.r.size
    cl_s, _ = st.direct(np.repeat(np.arange(n), al.size), np.tile(al, n),
                        np.repeat(dp["re"][0], al.size), np.repeat(rot.te_c, al.size))
    a_stall = al[np.argmax(cl_s.reshape(n, al.size), axis=1)]
    margin = a_stall - dp["alpha_deg"][0]
    chk["stall_margin_min_deg"] = float(margin[mid].min())
    chk["stall_margin_r_mm"] = float(st.r[mid][np.argmin(margin[mid])])
    chk["alpha_design_range_deg"] = (float(dp["alpha_deg"][0][mid].min()), float(dp["alpha_deg"][0][mid].max()))
    chk["re_design_range"] = (float(dp["re"][0].min()), float(dp["re"][0].max()))
    chk["re_range_all_mppt"] = (float(min(p["re"][0].min() for p in perf["mppt"])),
                                float(max(p["re"][0].max() for p in perf["mppt"])))
    cav = {}
    for tag, res in (("mppt_2p5_m_s", perf["mppt"][-1]), ("runaway_2p5_m_s", perf["runaway"][-1])):
        ae = st.direct.aero(np.arange(n), res["alpha_deg"][0], res["re"][0], rot.te_c)
        sig = (P_ATM + RHO * G * DEPTH_M - P_VAP) / (0.5 * RHO * res["W_m_s"][0] ** 2)
        m = sig + ae["cp_min"]
        k = int(np.argmin(m))
        cav[tag] = {"sigma_min": float(sig.min()), "cp_min_at_worst": float(ae["cp_min"][k]),
                    "margin_min": float(m[k]), "r_mm": float(st.r[k])}
    chk["cavitation"] = cav
    out["checks"] = chk

    # ---- 5. sensitivities (Cp peak at 1.5 m/s) ----------------------------------
    sens = []
    base_t, base_p = dp["tsr"][0], dp["cp"][0]

    def add(label, s_, r_):
        t_, p_ = s_.cp_max(r_)
        sens.append({"case": label, "cp_peak": float(p_["cp"][0]), "tsr_peak": t_,
                     "delta_cp_pct": 100 * (float(p_["cp"][0]) / base_p - 1)})
        log(f"  sens {label:44s} Cp {p_['cp'][0]:.4f} at {t_:.2f} ({sens[-1]['delta_cp_pct']:+.1f} %)")

    sens.append({"case": f"nominal ({CONE_NOMINAL:g} deg cone, n_crit {N_CRIT:g})", "cp_peak": float(base_p),
                 "tsr_peak": float(base_t), "delta_cp_pct": 0.0})
    for cone in (0.0, 8.0, 15.0):
        add(f"cone {cone:g} deg (straight axis)", st, replace(rot, cone_deg=cone))
    x_aft, cone_h = hull_following_axis(st.r)
    add("hull-following curved axis (estimate)", st, replace(rot, cone_deg=cone_h, x_aft_mm=x_aft))
    out["hull_axis"] = {"r_mm": st.r, "x_aft_mm": x_aft, "cone_deg": cone_h}
    from . import profile
    add(f"hub-loss radius {profile.radius(P.HINGE_X):.1f} mm (hub surface)", st,
        replace(rot, r_hub_mm=float(profile.radius(P.HINGE_X))))
    add("no hub loss", st, replace(rot, hub_loss=False))
    for label, kw in (("n_crit 9 (clean inflow, smooth blade)", {"n_crit": 9.0}),
                      ("forced transition x/c 0.05 (rough print)", {"xtr": 0.05}),
                      ("48 elements (discretisation)", {"n_el": 48})):
        s2 = Study(te_grid=TE_GRID_NARROW, **kw)
        add(label, s2, s2.linear_rotor(x_fin))
    out["sensitivity"] = sens

    return out


# ==========================================================================
# Post-processing: curves, structure estimates, JSON, figures
# ==========================================================================

WORK_DIR = ROOT / "out" / "work" / "bem"
CACHE = WORK_DIR / "study.pkl"
POST_CACHE = WORK_DIR / "post.pkl"
U_FINE = np.round(np.arange(1.0, 2.501, 0.1), 2)


def stadium(width: float, thick: float) -> dict:
    """Area (mm^2) and section moduli (mm^3) of a stadium (rounded slot) section:
    a width x thick rectangle with semicircular ends, width along its long axis."""
    a, r = width - thick, thick / 2
    area = a * thick + math.pi * r**2
    i_thin = a * thick**3 / 12 + math.pi * thick**4 / 64        # depth = thick
    d = 4 * r / (3 * math.pi)                                      # half-disc centroid
    i_half = math.pi * r**4 / 8 - (math.pi * r**2 / 2) * d**2
    i_wide = thick * a**3 / 12 + 2 * (i_half + (math.pi * r**2 / 2) * (a / 2 + d) ** 2)
    return {"area_mm2": area, "z_thin_mm3": i_thin / r, "z_wide_mm3": i_wide / (width / 2)}


def structure_estimates(out: dict) -> dict:
    """Blade mass, centrifugal relief of the flap moment and shank bending stress.
    Estimates: solid PLA, straight pitch axis at the nominal cone, the hinge flap
    and edge moments applied to the lug section (an upper bound for the shank)."""
    x = out["x_final"]
    perf = out["perf"]
    beta = math.radians(CONE_NOMINAL)
    rs = np.linspace(R_AERO_ROOT, R_TIP, 121)
    cs = linear(rs, x[0], x[1])
    area = np.array([section_area(float(t), P.MIN_TE / c) * c**2 for t, c in zip(thickness_ratio(rs), cs)])
    lin_mass = PLA_RHO * area * 1e-9 / math.cos(beta)               # kg per mm of radius
    lug = stadium(LUG_CHORD, LUG_THICKNESS)
    shank_len = (R_AERO_ROOT - P.HINGE_R) / math.cos(beta)
    m_shank = PLA_RHO * 0.5 * (lug["area_mm2"] + area[0]) * shank_len * 1e-9
    x_aft = (rs - P.HINGE_R) * math.tan(beta)
    i_cf = float(np.trapezoid(lin_mass * rs * x_aft, rs)) * 1e-6     # kg m^2, sum dm r x_aft
    res_cases = [("mppt_1p5_m_s", perf["mppt"][U_EVAL.index(DESIGN_U)]),
                 ("mppt_2p5_m_s", perf["mppt"][-1]),
                 ("runaway_2p5_m_s", perf["runaway"][-1])]
    cases = []
    for tag, res in res_cases:
        om = float(res["omega_rad_s"][0])
        mf, me = float(res["m_flap_Nm"][0]), float(res["m_edge_Nm"][0])
        m_cf = i_cf * om**2
        cases.append({
            "case": tag, "rpm": float(res["rpm"][0]),
            "flap_moment_hinge_Nm": mf, "edge_moment_hinge_Nm": me,
            "centrifugal_relief_Nm": m_cf, "centrifugal_relief_pct": 100 * m_cf / mf,
            "flap_stress_lug_thin_depth_MPa": mf * 1e3 / lug["z_thin_mm3"],
            "flap_stress_lug_wide_depth_MPa": mf * 1e3 / lug["z_wide_mm3"],
            "edge_stress_lug_thin_depth_MPa": me * 1e3 / lug["z_thin_mm3"],
            "edge_stress_lug_wide_depth_MPa": me * 1e3 / lug["z_wide_mm3"],
        })
    env = perf["envelope"][-1]
    return {
        "assumptions": "solid PLA 1240 kg/m^3, straight pitch axis at the nominal cone, "
                       "hinge moments applied to the 22 x 8 mm stadium lug (upper bound)",
        "blade_mass_g": 1e3 * (float(np.trapezoid(lin_mass, rs)) + m_shank),
        "blade_mass_aero_g": 1e3 * float(np.trapezoid(lin_mass, rs)),
        "lug_section": lug,
        "cases": cases,
        "max_flap_moment_2p5_m_s_Nm": env["max_flap_moment_hinge_per_blade_Nm"],
        "max_thrust_2p5_m_s_N": env["max_thrust_N"],
    }


def post(out: dict, st: Study, log=print) -> dict:
    """Curves and tables for the JSON and figures (cheap, run after run_study)."""
    tsr_d = out["tsr_d"]
    chain = {c["key"]: c for c in out["chain"]}
    rot = out["rotor"]
    curves = {"tsr": TSR_SWEEP}
    for key in ("optimum", "schmitz", "linear"):
        c = chain[key]
        r_ = st.rotor(c["chord"], c["twist"], te_mm=c["te_mm"], cone=c["cone_deg"])
        curves[key] = st.perf(r_, DESIGN_U, TSR_SWEEP)["cp"]
    for u in U_EVAL:
        curves[f"final_{u:g}"] = out["perf"]["sweep"][u]["cp"]
    out["curves"] = curves
    log("  curves done")
    fine = {"u": U_FINE, "mppt": [], "fixed": None, "runaway": []}
    for u in U_FINE:
        _, p = st.cp_max(rot, float(u))
        fine["mppt"].append(p)
        sw = solve(rot, st.table, float(u), TSR_SWEEP)
        fine["runaway"].append(st.runaway(rot, float(u), sw))
    fine["fixed"] = st.perf(rot, U_FINE, tsr_d)
    out["fine"] = fine
    log("  fine speed grid done")
    out["structure"] = structure_estimates(out)
    return out


# ---- JSON ---------------------------------------------------------------------

def _f(v, nd: int = 4) -> float:
    return round(float(v), nd)


def _l(v, nd: int = 4) -> list:
    return [round(float(x), nd) for x in np.ravel(v)]


def station_table(x_fin, n: int = N_STATIONS, cone: float = CONE_NOMINAL) -> dict:
    r = np.linspace(R_AERO_ROOT, R_TIP, n)
    c = linear(r, x_fin[0], x_fin[1])
    tw = linear(r, x_fin[2], x_fin[3])
    tc = thickness_ratio(r)
    te_c = P.MIN_TE / c
    tmax = np.array([section_max_thickness(float(a), float(b)) for a, b in zip(tc, te_c)]) * c
    tb = math.tan(math.radians(cone))
    return {
        "radius_mm": _l(r, 3), "chord_mm": _l(c, 3), "twist_deg": _l(tw, 3),
        "thickness_ratio": _l(tc, 4), "naca": [f"NACA {naca_code(t)}" for t in tc],
        "te_thickness_mm": [P.MIN_TE] * n, "te_thickness_ratio": _l(te_c, 4),
        "max_thickness_mm": _l(tmax, 3),
        "pitch_axis_x_mm": _l(P.HINGE_X + (r - P.HINGE_R) * tb, 3),
        "pitch_axis_aft_of_hinge_mm": _l((r - P.HINGE_R) * tb, 3),
        "span_from_hinge_mm": _l((r - P.HINGE_R) / math.cos(math.radians(cone)), 3),
    }


def build_json(out: dict) -> dict:
    x = out["x_final"]
    tsr_d = out["tsr_d"]
    perf, fine, chk = out["perf"], out["fine"], out["checks"]
    i15 = U_EVAL.index(DESIGN_U)
    dp = perf["mppt"][i15]
    chain = {c["key"]: c for c in out["chain"]}
    st_r = out["rotor"].r_mm

    def table(results):
        return [_ops(res) for res in results]

    def fine_table(results):
        return [_ops(res) for res in results]

    fixed = fine["fixed"]
    q_d = next(r_["torque_at_1p5_m_s_Nm"] for r_ in out["tsr_rows"] if r_["design_tsr"] == tsr_d)
    tsr_rows = []
    for row in out["tsr_rows"]:
        rr = {k: (_l(v, 2) if isinstance(v, list) else _f(v, 4)) for k, v in row.items()}
        rr["copper_loss_index"] = _f((row["torque_at_1p5_m_s_Nm"] / q_d) ** 2, 3)
        tsr_rows.append(rr)
    sw = perf["sweep"][DESIGN_U]
    q_pk = float(dp["torque_Nm"][0])
    gen_rows = [dict(_ops(dp), copper_loss_index=1.0, label="Cp peak (MPPT)")]
    for t_ in (4.0, 4.5, 5.0, 5.5):
        j = int(np.argmin(np.abs(sw["tsr"] - t_)))
        gen_rows.append(dict(_ops(sw, j), copper_loss_index=_f((float(sw["torque_Nm"][j]) / q_pk) ** 2, 3),
                             label=f"held at TSR {t_:g}"))
    ha = out["hull_axis"]
    return {
        "schema": "hydroegg.blade_design/1",
        "generated_by": "python -m hydroegg.bem",
        "conventions": {
            "frame": "X = body axis, flow +X (nose at x = 0), Z up, right-handed",
            "radius": "radius_mm is perpendicular to the body axis (not along the coned span)",
            "rotation": "right-handed about +X (Omega vector +X): CLOCKWISE seen from upstream "
                        "(observer at the nose looking downstream). A blade pointing +Z moves toward -Y.",
            "twist": "angle between chord line and the plane of rotation, measured in the section "
                     "plane normal to the pitch axis; POSITIVE = leading edge rotated upstream (-X). "
                     "Angle of attack = inflow angle phi - twist.",
            "section_orientation": "suction (convex, upper) side faces downstream (+X); leading "
                                   "edge points in the direction of rotation",
            "section_placement": "for the blade whose axis is in the +Z half-plane: pitch-axis point "
                                 "A = (pitch_axis_x_mm, 0, radius_mm); e_n = (cos b, 0, -sin b) with "
                                 "b = cone; d = cos(tw) * (+Y) + sin(tw) * e_n (LE -> TE); "
                                 "n = cos(tw) * e_n - sin(tw) * (+Y) (toward suction side). "
                                 "Section point = A + chord * ((xs - pitch_axis) * d + ys * n) for "
                                 "unit-chord coordinates (xs, ys) of hydroegg.bem.naca44_section.",
            "cone": "aft (+X) lean of the pitch axis out of the plane of rotation, about the "
                    "tangential hinge pin",
        },
        "rotor": {
            "n_blades": P.N_BLADES, "r_tip_mm": R_TIP, "r_aero_root_mm": R_AERO_ROOT,
            "hinge_r_mm": P.HINGE_R, "hinge_x_mm": P.HINGE_X, "cone_deg": CONE_NOMINAL,
            "cone_note": "straight pitch axis at the nominal cone; the hull-following curved axis "
                         "is evaluated as a sensitivity (see sensitivity)",
            "design_tsr": tsr_d, "design_flow_m_s": DESIGN_U,
            "tsr_at_cp_peak_1p5_m_s": _f(dp["tsr"][0], 3),
            "swept_diameter_mm": 2 * R_TIP,
            "rotation_seen_from_upstream": "clockwise",
            "omega_vector": "+X",
        },
        "section": {
            "family": "NACA 44xx (4 % camber at 40 % chord), t/c linear root to tip",
            "camber_m": CAMBER_M, "camber_p": CAMBER_P,
            "thickness_ratio_root": TC_ROOT, "thickness_ratio_tip": TC_TIP,
            "te_thickness_mm": P.MIN_TE,
            "te_construction": "closed-TE NACA thickness (-0.1036 x^4) plus a linear add-on "
                               "0.5 * te * (x/c) on each surface normal to the chord (CST TE term)",
            "pitch_axis_chord_fraction": PITCH_AXIS,
            "naca_code_note": "codes are rounded labels; use thickness_ratio for geometry",
        },
        "laws": {
            "chord_mm": {"type": "linear in radius", "root": _f(x[0], 2), "tip": _f(x[1], 2)},
            "twist_deg": {"type": "linear in radius", "root": _f(x[2], 2), "tip": _f(x[3], 2)},
            "thickness_ratio": {"type": "linear in radius", "root": TC_ROOT, "tip": TC_TIP},
        },
        "stations": station_table(x),
        "pitch_axis_hull_following_estimate": {
            "note": "estimate only, for the layout study: pitch axis that would copy the hull when "
                    "folded 90 deg forward about the hinge; Cp change vs the straight 12 deg axis is "
                    "in sensitivity",
            "element_radius_mm": _l(ha["r_mm"], 2), "x_aft_of_hinge_mm": _l(ha["x_aft_mm"], 2),
            "local_cone_deg": _l(ha["cone_deg"], 2),
        },
        "shank": {
            "r_start_mm": P.HINGE_R, "r_end_mm": R_AERO_ROOT,
            "length_perpendicular_mm": R_AERO_ROOT - P.HINGE_R,
            "length_along_axis_mm": _f((R_AERO_ROOT - P.HINGE_R) / math.cos(math.radians(CONE_NOMINAL)), 2),
            "lug_section": "stadium (rounded slot)", "lug_chord_mm": LUG_CHORD,
            "lug_thickness_mm": LUG_THICKNESS, "pin_d_mm": P.PIN_D,
            "blend": "loft from the lug at the hinge to the NACA root section at r_end_mm; "
                     "aerodynamically inactive (BEM span starts at r_aero_root_mm, hub loss there)",
            "max_chord_mm": LUG_CHORD,
        },
        "performance": {
            "note": "MPPT = generator load held at the Cp-peak TSR; fixed = load held at the design "
                    "TSR; runaway = zero torque (free-wheeling); flap/edge moments per blade about "
                    "the hinge pin, flap positive pushes the blade aft onto its deployed stop",
            "mppt": table(perf["mppt"]),
            "design_tsr": table(perf["design_tsr"]),
            "runaway": table(perf["runaway"]),
            "envelope_up_to_runaway": [{k: (_f(v) if isinstance(v, float) else v) for k, v in e.items()}
                                       for e in perf["envelope"]],
            "fine_mppt": fine_table(fine["mppt"]),
            "fine_fixed_design_tsr": [_ops(fixed, j) for j in range(len(fine["u"]))],
            "fine_runaway": fine_table(fine["runaway"]),
        },
        "cp_tsr": {
            "tsr": _l(out["curves"]["tsr"], 2),
            "flow_m_s": DESIGN_U,
            "cp_final": _l(out["curves"]["final_1.5"], 4),
            "cp_bem_optimum": _l(out["curves"]["optimum"], 4),
            "cp_schmitz": _l(out["curves"]["schmitz"], 4),
            "cp_final_by_flow": {f"{u:g}_m_s": _l(out["curves"][f"final_{u:g}"], 4) for u in U_EVAL},
        },
        "design_point_1p5_m_s": {
            "tsr": _f(dp["tsr"][0], 3), "element_radius_mm": _l(st_r, 2),
            "a": _l(dp["a"][0]), "a_prime": _l(dp["ap"][0]), "alpha_deg": _l(dp["alpha_deg"][0], 2),
            "phi_deg": _l(dp["phi_deg"][0], 2), "re": _l(dp["re"][0], 0), "cl": _l(dp["cl"][0]),
            "cd": _l(dp["cd"][0], 5), "loss_factor_F": _l(dp["F"][0]),
        },
        "tsr_study": tsr_rows,
        "tsr_study_note": "cp_constrained = Cp of the printable linear blade optimised for that TSR, "
                          "flown at that TSR; tsr_constrained_peak = where its Cp actually peaks. "
                          "copper_loss_index = (torque / torque of the chosen design)^2 at 1.5 m/s, "
                          "same generator (copper loss scales with current squared, current with torque)",
        "generator_tradeoff_1p5_m_s": gen_rows,
        "generator_tradeoff_note": "final blade at 1.5 m/s held at TSR at or above its Cp peak: "
                                   "copper_loss_index = (torque / torque at the Cp peak)^2",
        "tsr_rule": f"highest candidate TSR whose constrained Cp is within {100 * TSR_TIE:g} % of the best",
        "constraint_chain": [{"key": c["key"], "label": c["label"], "te_mm": c["te_mm"],
                              "cone_deg": c["cone_deg"], "cp_design_tsr": _f(c["cp_design_tsr"]),
                              "cp_peak": _f(c["cp_peak"]), "tsr_peak": _f(c["tsr_peak"], 3),
                              "chord_root_mm": _f(c["chord_root_mm"], 2), "chord_tip_mm": _f(c["chord_tip_mm"], 2),
                              "twist_root_deg": _f(c["twist_root_deg"], 2), "twist_tip_deg": _f(c["twist_tip_deg"], 2)}
                             for c in out["chain"]],
        "optimum_distribution": {
            "element_radius_mm": _l(st_r, 2),
            "bem_optimum_chord_mm": _l(chain["optimum"]["chord"], 2),
            "bem_optimum_twist_deg": _l(chain["optimum"]["twist"], 2),
            "schmitz_chord_mm": _l(chain["schmitz"]["chord"], 2),
            "schmitz_twist_deg": _l(chain["schmitz"]["twist"], 2),
        },
        "sensitivity": [{k: (_f(v) if isinstance(v, float) else v) for k, v in s.items()}
                        for s in out["sensitivity"]],
        "checks": json.loads(json.dumps(chk, default=float)),
        "glauert_verification": out["glauert"],
        "structure_estimates": json.loads(json.dumps(out["structure"], default=float)),
        "model": {
            "polars": f"NeuralFoil {MODEL_SIZE}, n_crit {N_CRIT:g}, free transition",
            "fluid": {"rho_kg_m3": RHO, "nu_m2_s": NU},
            "elements": N_EL, "tip_loss": "Prandtl", "hub_loss": f"Prandtl at r = {R_AERO_ROOT:g} mm",
            "high_induction": "Buhl (2005), closed form of Ning (2014)",
            "unconverged_elements_study_total": int(perf["n_unconverged"]),
            "fallback_elements_study_total": int(perf["n_fallback"]),
        },
    }


def write_json(out: dict, path: Path = DATA_JSON) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(build_json(out), indent=1) + "\n")


# ---- figures --------------------------------------------------------------------

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e2e1dc", "#fcfcfb"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"                    # categorical slots 1-3
SEQ = ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")               # blue ramp 250/400/550/700


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelsize": 9, "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.8,
        "lines.solid_capstyle": "round", "legend.frameon": False, "legend.fontsize": 8.5,
        "font.size": 9, "text.color": INK,
    })
    return plt


def _ref(ax, y=None, x=None, text="", where="right", va="bottom", xpos=None):
    """Thin grey reference line with a small label (at xpos if given, else at the edge)."""
    kw = dict(color=INK2, lw=0.8, alpha=0.7, zorder=1)
    if y is not None:
        ax.axhline(y, **kw)
        xl = ax.get_xlim()
        if xpos is not None:
            ax.text(xpos, y, text, ha="left", va=va, fontsize=7.5, color=INK2)
        else:
            ax.text(xl[1] if where == "right" else xl[0], y, text, ha=where, va=va, fontsize=7.5, color=INK2)
    else:
        ax.axvline(x, **kw)
        yl = ax.get_ylim()
        ax.text(x, yl[1], text, ha="left", va="top", fontsize=7.5, color=INK2, rotation=90)


def make_figures(out: dict, fig_dir: Path = FIG_DIR) -> list[Path]:
    plt = _style()
    fig_dir.mkdir(parents=True, exist_ok=True)
    chain = {c["key"]: c for c in out["chain"]}
    r = out["rotor"].r_mm
    x = out["x_final"]
    tsr_d = out["tsr_d"]
    paths = []
    rr = np.linspace(P.HINGE_R, R_TIP, 200)

    # 1. chord and twist ------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
    for ax, key, unit in ((axs[0], "chord", "chord (mm)"), (axs[1], "twist", "twist (deg)")):
        ax.axvspan(P.HINGE_R, R_AERO_ROOT, color=GRID, alpha=0.6, lw=0, zorder=0)
        ax.plot(r, chain["optimum"][key], color=S2, label="BEM optimum (sharp TE, no cone)")
        ax.plot(r, chain["schmitz"][key], color=S3, label="Schmitz")
        ax.plot(r, chain["linear"][key], color=INK2, lw=1.0, ls=(0, (4, 3)),
                label="linear fit (0.7R-0.9R) of the blunt-TE optimum")
        lo, hi = (x[0], x[1]) if key == "chord" else (x[2], x[3])
        ax.plot([R_AERO_ROOT, R_TIP], [lo, hi], color=S1, lw=2.2, label="final (printable)")
        ax.set_xlabel("radius r (mm)")
        ax.set_ylabel(unit)
        ax.set_xlim(P.HINGE_R, R_TIP + 2)
    axs[0].set_ylim(0, 45)
    axs[1].set_ylim(0, 29)
    axs[0].set_title(f"Chord, design TSR {tsr_d:g}")
    axs[1].set_title("Twist (positive = LE upstream)")
    _ref(axs[0], y=P.MIN_CHORD, text=f"min chord {P.MIN_CHORD:g}", va="top", xpos=34.0)
    _ref(axs[0], y=MAX_CHORD, text=f"max chord {MAX_CHORD:g}", xpos=62.0)
    axs[0].plot([P.HINGE_R, R_AERO_ROOT], [LUG_CHORD, LUG_CHORD], color=INK2, lw=1.4)
    axs[0].text(0.5 * (P.HINGE_R + R_AERO_ROOT), LUG_CHORD + 0.8, f"shank\n<= {LUG_CHORD:g}",
                ha="center", va="bottom", fontsize=7.5, color=INK2)
    _ref(axs[1], y=P.MAX_ROOT_TWIST, text=f"root twist cap {P.MAX_ROOT_TWIST:g}", xpos=62.0)
    axs[1].legend(loc="lower left", bbox_to_anchor=(0.17, 0.0))
    p = fig_dir / "bem_chord_twist.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths.append(p)

    # 2. Cp - TSR ---------------------------------------------------------------
    cv = out["curves"]
    t = cv["tsr"]
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True, sharey=True)
    ax = axs[0]
    ax.plot(t, cv["optimum"], color=S2, label="BEM optimum")
    ax.plot(t, cv["schmitz"], color=S3, label="Schmitz")
    ax.plot(t, cv["final_1.5"], color=S1, lw=2.2, label="final (printable)")
    ax.set_title(f"Cp vs TSR at {DESIGN_U:g} m/s")
    ax.legend(loc="lower center")
    ax2 = axs[1]
    for u, col in zip(U_EVAL, SEQ):
        ax2.plot(t, cv[f"final_{u:g}"], color=col, label=f"{u:g} m/s")
    ax2.set_title("Final blade, flow speed (Re) effect")
    ax2.legend(loc="lower center", ncol=2)
    for a_ in axs:
        a_.set_xlim(1.5, 10)
        a_.set_ylim(0, 0.62)
        a_.set_xlabel("tip speed ratio")
        _ref(a_, y=16 / 27, text="Betz 0.593", where="right", va="bottom")
        _ref(a_, x=tsr_d, text=f" design TSR {tsr_d:g}")
    axs[0].set_ylabel("power coefficient Cp")
    p = fig_dir / "bem_cp_tsr.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths.append(p)

    # 3. power, rpm, torque, thrust vs flow speed -------------------------------
    fine = out["fine"]
    u = fine["u"]
    mp = {k: np.array([float(res[k][0]) for res in fine["mppt"]]) for k in ("power_W", "rpm", "torque_Nm", "thrust_N")}
    rw = {k: np.array([float(res[k][0]) for res in fine["runaway"]]) for k in ("rpm", "thrust_N")}
    fx = fine["fixed"]
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 6.2), constrained_layout=True, sharex=True)
    ax = axs[0, 0]
    ax.plot(u, mp["power_W"], color=S1, lw=2.2, label="MPPT (Cp peak)")
    ax.plot(u, fx["power_W"], color=S2, label=f"fixed TSR {tsr_d:g}")
    ax.set_title("Shaft power (W)")
    ax.legend(loc="upper left")
    ax = axs[0, 1]
    ax.plot(u, mp["rpm"], color=S1, lw=2.2, label="MPPT")
    ax.plot(u, fx["rpm"], color=S2, label=f"fixed TSR {tsr_d:g}")
    ax.plot(u, rw["rpm"], color=S3, label="runaway (no load)")
    ax.set_title("Rotor speed (rpm)")
    ax.legend(loc="upper left")
    ax = axs[1, 0]
    ax.plot(u, mp["torque_Nm"], color=S1, lw=2.2, label="MPPT")
    ax.plot(u, fx["torque_Nm"], color=S2, label=f"fixed TSR {tsr_d:g}")
    ax.set_title("Shaft torque (N m)")
    ax.legend(loc="upper left")
    ax = axs[1, 1]
    ax.plot(u, mp["thrust_N"], color=S1, lw=2.2, label="MPPT")
    ax.plot(u, rw["thrust_N"], color=S3, label="runaway (no load)")
    ax.set_title("Rotor thrust (N)")
    ax.legend(loc="upper left")
    for a_ in axs[1]:
        a_.set_xlabel("flow speed (m/s)")
    for a_ in axs.ravel():
        a_.set_ylim(bottom=0)
        a_.axvline(DESIGN_U, color=INK2, lw=0.8, alpha=0.5, zorder=1)
    for a_, key in ((axs[0, 0], "power_W"), (axs[0, 1], "rpm"), (axs[1, 0], "torque_Nm"), (axs[1, 1], "thrust_N")):
        j = int(np.argmin(np.abs(u - DESIGN_U)))
        a_.annotate(f"{mp[key][j]:.3g}", (u[j], mp[key][j]), xytext=(4, -10), textcoords="offset points",
                    fontsize=7.5, color=INK)
        a_.annotate(f"{mp[key][-1]:.3g}", (u[-1], mp[key][-1]), xytext=(-4, 4), textcoords="offset points",
                    fontsize=7.5, color=INK, ha="right")
    p = fig_dir / "bem_power.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths.append(p)

    # 4. constraint cost ------------------------------------------------------------
    keys = ["optimum", "te", "linear", "chord", "twist", "cone", "final"]
    labels = [chain[k]["label"] for k in keys]
    cps = np.array([chain[k]["cp_design_tsr"] for k in keys])
    fig, ax = plt.subplots(figsize=(9.2, 3.4), constrained_layout=True)
    yy = np.arange(len(keys))[::-1]
    ax.barh(yy, cps, height=0.5, color=[S1 if k == "final" else "#86b6ef" for k in keys], zorder=2)
    for yv, cpv, k in zip(yy, cps, keys):
        d = 100 * (cpv / cps[0] - 1)
        ax.text(cpv + 0.004, yv, f"{cpv:.3f}" + ("" if k == "optimum" else f"  ({d:+.1f} %)"),
                va="center", fontsize=8, color=INK)
    ax.set_yticks(yy, labels)
    ax.set_xlim(0, max(cps) * 1.25)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(f"Cp at design TSR {tsr_d:g}, {DESIGN_U:g} m/s")
    ax.set_title("Cp cost of each constraint (cumulative, top to bottom)")
    p = fig_dir / "bem_constraint_cost.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths.append(p)
    return paths


def main(argv=None) -> None:
    import argparse
    import pickle

    ap = argparse.ArgumentParser(description="HydroEgg BEM blade design study")
    ap.add_argument("--cached", action="store_true",
                    help=f"reuse the cached study ({CACHE.relative_to(ROOT)}) instead of re-running it")
    ap.add_argument("--cached-post", action="store_true", help="also reuse the cached post-processing")
    args = ap.parse_args(argv)
    import time
    t0 = time.time()
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    if args.cached_post and POST_CACHE.exists():
        out = pickle.loads(POST_CACHE.read_bytes())
    else:
        if args.cached and CACHE.exists():
            out = pickle.loads(CACHE.read_bytes())
            st = Study()
        else:
            out = run_study()
            st = out.pop("study")
            CACHE.write_bytes(pickle.dumps(out))
        post(out, st)
        POST_CACHE.write_bytes(pickle.dumps(out))
    write_json(out)
    print(f"wrote {DATA_JSON.relative_to(ROOT)}")
    for p in make_figures(out):
        print(f"wrote {p.relative_to(ROOT)}")
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    # Run through the importable module so cached pickles reference hydroegg.bem classes.
    from hydroegg import bem as _bem
    _bem.main()
