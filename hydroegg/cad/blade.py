"""Folding blade: 3 identical blades at psi = 0 / 120 / 240 (builder: blade).

Geometry source: the rotor layout (hydroegg/layout.py, data/layout_rotor.json). The aero
span is a smooth C2 B-spline surface interpolated through the layout's own sections
(``layout.AeroStations`` on the solved curved pitch axis of data/layout_rotor.json), so the
CAD is the surface the static-tail pockets and the hub slots were cut for. Everything is
built in the BLADE FRAME of hydroegg/layout.py and carried to the world with
``layout.blade_to_world``:

    local X = e_s  span (deployed: radially out)
    local Y = e_c  chord, along the hinge pin (toward the leading edge)
    local Z = e_n  thickness (deployed: aft, the suction side faces +n)
    origin  = hinge point H, on the pin axis

Construction: every section is interpolated with the same parameters (same knots), then
the pole columns are interpolated along the span (``_bspline_surface``), so each layout
section is an exact, planar iso-curve of one bicubic skin; a ruled strip closes the blunt
trailing edge and two planar caps close the ends. The skin is fused with the knuckle and
heel (analytic solids), then the end chamfers, the pin bore and the magnet seats are cut.

Parts of the blade (blade frame):

* aero span r 30..91: the layout's NACA 44xx sections (data/blade_design.json) at the
  layout's own 41-station grid, blunt TE of MIN_TE.
* tip r 91..97: sections shrink elliptically in chord toward TIP_ANCHOR (the LE sweeps
  back, the TE stays nearly straight so it prints on the same support line) with a
  thickness floor, ending in a 2.8 x 0.9 mm blunt cap. Every tip section lies inside the
  nominal section at its radius (checked), so the pockets cut for the layout still fit.
* shank: layout.shank_sections (knuckle rounded-rectangle morphing into the root airfoil,
  55 deg cone at the knuckle) with CAD-only changes that stay inside the layout envelope
  or inside the knuckle: flanks SHANK_W wide (clear of the knuckle end chamfers), rounder
  rectangle corners (RECT_EXP 3 instead of 4), the rectangle's end opened to the blunt TE
  (x^8, local to its end), the first sections shrunk so the skin starts inside the
  knuckle, and the shank sections within SHANK_JOIN_GAP of the aero root plane dropped:
  layout.shank_sections ends 0.63 mm off the aero root section, the skin bridges it.
* knuckle: D KNUCKLE_D x KNUCKLE_L about the pin, 0.3 x 45 deg end chamfers, pin bore
  ROTOR_HUB knuckle_bore_d with 0.3 mm lead-in chamfers.
* heel: HEEL_REACH x HEEL_T x HEEL_W under the knuckle; its -n face is the flat deployed
  stop (x 113.5 deployed, ROTOR_HUB stop_face_x). Tip and planform corners rounded 1 mm,
  +n root filleted 0.6 mm into the knuckle; the stop face, its perimeter and its
  junction with the knuckle stay sharp (the hub stop block sits there). OCC propagates
  fillets along tangent chains, so the fillet order is fixed and asserted.
* magnet seats (4 x 2 N52, D 4.1 glue / light press, depth MAGNET H + 0.1 so a seated
  magnet sits 0.1 below its face): deployed seat on the DETENTS axis in the heel stop
  face; stowed seat where the DETENTS axis meets the shank's inner face but along the
  face NORMAL (that face is inclined ~51 deg to the radial DETENTS axis, so a radial seat
  cannot be flush), with a 0.2 mm mouth chamfer.

parts(): blade_0..2 in the WORLD STOWED pose (fold FOLD_DEG) with a pivot such that
rotating the stowed shape by deployed_deg about the axis through the hinge point gives
the deployed pose (proved numerically by ``checks()``).

Run ``scripts/cadrun.sh .venv/bin/python -m hydroegg.cad.blade`` to build, check and
export (out/cad/step/blade.step, out/cad/stl/blade.stl, out/cad/checks/blade.json and the
preview out/work/cad/blade/final.png).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from build123d import (
    Axis,
    Edge,
    Face,
    GeomType,
    Location,
    Plane,
    Shape,
    Shell,
    Solid,
    Vector,
    Wire,
    export_step,
    export_stl,
)
from OCP.Approx import Approx_ParametrizationType
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCP.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf
from OCP.gp import gp_Pnt, gp_Trsf
from OCP.collections import Array2_gp_Pnt as TColgp_Array2OfPnt

from .. import layout as LY
from .. import params as P
from . import common as C

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "out" / "work" / "cad" / "blade"
STEP_OUT = ROOT / "out" / "cad" / "step" / "blade.step"
STL_OUT = ROOT / "out" / "cad" / "stl" / "blade.stl"
CHECK_OUT = ROOT / "out" / "cad" / "checks" / "blade.json"

# --------------------------------------------------------------------------
# Blade-builder choices (CAD detail below the layout's resolution; documented)
# --------------------------------------------------------------------------
TIP_ZONE = 6.0            # mm of span (r 91..97) over which the tip rounds
TIP_G_MIN = 0.2           # chord of the last tip section / nominal tip chord
TIP_ANCHOR = 0.85        # chord fraction the tip sections shrink toward (TE nearly straight)
TIP_T_MIN = 0.9           # mm, max thickness of the last tip section (>= 1.5 x MIN_TE)
TIP_N = 10                # tip sections (sin-spaced, dense at the very tip)
SHANK_N = 16              # shank sections, as layout.blade_surface
SHANK_JOIN_GAP = 1.8      # last kept shank section this far before the aero root plane
KNUCKLE_CHAMFER = 0.3     # knuckle end edges (also elephant-foot relief on the bed)
SHANK_W = P.KNUCKLE_L - 2 * KNUCKLE_CHAMFER - 0.2  # 13.2: shank flanks 0.1 inside the chamfer start
RECT_TE_EXP = 8           # knuckle-rectangle TE opening ~ x^8 (local to the TE end)
RECT_EXP = 3.0            # knuckle-rectangle super-ellipse exponent (layout 4: ~0.4 mm corners;
                          # 3 rounds the shank flank corners, removing material only)
ROOT_SHRINK = 0.15        # first shank section shrunk 15 % in thickness, so the loft
ROOT_SHRINK_LEN = 2.5     # starts strictly inside the knuckle; back to 1 at s = 2.5
BORE_D = float(P.ROTOR_HUB["knuckle_bore_d"])  # 3.3 free running on the 3 mm pin
BORE_CHAMFER = 0.3
MAG_D, MAG_H = float(P.MAGNET["D"]), float(P.MAGNET["H"])
MAG_SEAT_D = MAG_D + 0.1  # glue / light press (printed holes close up ~0.1)
MAG_SEAT_DEPTH = MAG_H + 0.1
MAG_MOUTH_CHAMFER = 0.2   # stowed seat only (see HEEL_SEAT_CHAMFER)
HEEL_SEAT_CHAMFER = 0.0   # heel seat: no chamfer, its rim is only 0.39 mm from the heel tip
HEEL_TIP_FILLET = 1.0     # heel tip, +n edge and the two planform corners
HEEL_SIDE_FILLET = 0.0    # heel long edges stay sharp: a fillet there chains into the root fillet
HEEL_ROOT_FILLET = 0.6    # heel +n face to knuckle (never the stop-face side)
PLA_G_MM3 = 1.24e-3
STL_TOL = 0.05


# ==========================================================================
# Layout context
# ==========================================================================

def axis_from_json(blade: LY.BladeDef) -> LY.BladeAxis:
    """The solved curved pitch axis exactly as published in data/layout_rotor.json."""
    d = json.loads(LY.DATA_JSON.read_text())["axis"]
    return LY.BladeAxis(r=np.asarray(d["r_mm"], float), s=np.asarray(d["s_mm"], float),
                        lead=np.asarray(d["lead_mm"], float), aft=np.asarray(d["aft_mm"], float),
                        slope=np.asarray(d["slope_daft_ds"], float),
                        aft_coef=np.asarray(d["aft_poly_coef_mm"], float),
                        lead_coef=np.asarray(d["lead_poly_coef_mm"], float),
                        r_root=float(blade.r_root), r_tip=float(blade.r_tip))


@lru_cache(maxsize=1)
def context() -> tuple[LY.BladeDef, LY.Layout, LY.BladeAxis]:
    blade = LY.load_blade()
    lay = LY.chosen_layout()
    return blade, lay, axis_from_json(blade)


# ==========================================================================
# Sections (blade frame, (Np, 3) point loops: upper TE -> LE -> lower TE)
# ==========================================================================

def _pitch_point(r: float):
    blade, lay, ax = context()
    ld = float(ax.lead_at(r))
    return (float(LY.s_of_r(r, ld, lay.hinge_r)), ld, float(ax.aft_at(r)), float(ax.slope_at(r, ld)))


def section_unit(tc: float, te_mm: float, c: float, g: float = 1.0, xa: float = TIP_ANCHOR,
                 t_mm: float | None = None, n: int = LY.N_PER_SIDE) -> np.ndarray:
    """(2n-1, 2) section in NOMINAL chord units (x from the nominal LE, y toward the
    suction side): a NACA 44xx of chord g * c shrunk toward chord fraction xa and laid on
    the nominal camber line, max thickness t_mm (default tc g c), blunt TE te_mm.
    g = 1 reproduces layout.naca44(tc, te_mm / c) exactly."""
    blade, _, _ = context()
    m, p = blade.camber_m, blade.camber_p
    b = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(b))
    t_small = tc if t_mm is None else t_mm / (g * c)
    yt = 5 * t_small * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    xp = xa + g * (x - xa)
    yc = np.where(xp <= p, m / p**2 * (2 * p * xp - xp**2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * xp - xp**2))
    dyc = np.where(xp <= p, 2 * m / p**2 * (p - xp), 2 * m / (1 - p) ** 2 * (p - xp))
    th = np.arctan(dyc)
    h = g * yt
    ramp = 0.5 * (te_mm / c) * x
    xu, yu = xp - h * np.sin(th), yc + h * np.cos(th) + ramp
    xl, yl = xp + h * np.sin(th), yc - h * np.cos(th) - ramp
    return np.vstack([np.column_stack([xu[::-1], yu[::-1]]), np.column_stack([xl[1:], yl[1:]])])


def place_section(xy: np.ndarray, r: float) -> np.ndarray:
    """Blade-frame points of a nominal-chord-unit section at deployed radius r."""
    blade, _, _ = context()
    s, ld, aft, sl = _pitch_point(r)
    c, tw = float(blade.chord(r)), float(blade.twist(r))
    L, U = LY.section_frame(sl, tw)
    Pp = np.array([s, ld, aft])
    return Pp + np.outer((blade.pitch_axis - xy[:, 0]) * c, L) + np.outer(xy[:, 1] * c, U)


def tip_law(u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Chord factor g and thickness add-on weight w over the tip zone (u 0 -> 1)."""
    u = np.asarray(u, float)
    g = TIP_G_MIN + (1 - TIP_G_MIN) * np.sqrt(np.clip(1 - u**2, 0, None))
    return g, u**3


def tip_sections() -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Tip-zone sections (u = 0 is a nominal section, u = 1 the end cap)."""
    blade, _, _ = context()
    k = np.arange(TIP_N + 1)
    u = np.sin(0.5 * math.pi * k / TIP_N)
    r = blade.r_tip - TIP_ZONE * (1 - u)
    g, w = tip_law(u)
    c_tip = float(blade.chord(blade.r_tip))
    t_add = max(0.0, TIP_T_MIN - float(blade.t_c(blade.r_tip)) * c_tip * TIP_G_MIN)
    out, meta = [], []
    for rk, gk, wk in zip(r, g, w):
        c, tc = float(blade.chord(rk)), float(blade.t_c(rk))
        t_mm = tc * c * gk + t_add * wk
        xy = section_unit(tc, blade.te_mm, c, gk, TIP_ANCHOR, t_mm)
        out.append(place_section(xy, float(rk)))
        meta.append({"r": float(rk), "g": float(gk), "chord_mm": c * gk, "t_max_mm": t_mm, "xy": xy, "c": c})
    return np.array(out), r, meta


def shank_sections(width: float = SHANK_W, te_mm: float | None = None, root_shrink: float = ROOT_SHRINK,
                   align_root: bool = False, n: int = SHANK_N, rect_exp: float | None = None
                   ) -> tuple[np.ndarray, np.ndarray]:
    """layout.shank_sections with the CAD changes listed in the module docstring.
    shank_sections(P.KNUCKLE_L, 0.0, 0.0, False) reproduces layout.shank_sections.
    align_root=True shifts the spine along s so the last section is the aero root section
    (not used: it moves the shank up to 0.6 mm off the layout; section_grid drops the
    shank sections that reach the aero root plane instead)."""
    blade, lay, ax = context()
    te_mm = blade.te_mm if te_mm is None else te_mm
    rect_exp = RECT_EXP if rect_exp is None else rect_exp
    nps = LY.N_PER_SIDE
    r0 = blade.r_root
    s0, ld0, aft0, sl0 = float(ax.s[0]), float(ax.lead[0]), float(ax.aft[0]), float(ax.slope[0])
    c0, tw0 = float(blade.chord(r0)), float(blade.twist(r0))
    L0, _ = LY.section_frame(sl0, tw0)
    M0 = np.array([s0, ld0, aft0]) + (blade.pitch_axis - 0.5) * c0 * L0
    m_k = math.tan(math.radians(lay.shank_cone_deg))
    s = np.linspace(0.0, s0, n)
    n_mid = LY._hermite(s, 0.0, 0.0, m_k, s0, M0[2], sl0)
    sl = LY._hermite_d(s, 0.0, 0.0, m_k, s0, M0[2], sl0)
    sb = min(lay.shank_blend_start, 0.8 * s0)
    c_mid = np.where(s <= sb, 0.0, LY._hermite(np.maximum(s, sb), sb, 0.0, 0.0, s0, M0[1], 0.0))
    tt = np.clip((s - sb) / (s0 - sb), 0, 1)
    u = tt * tt * (3 - 2 * tt)
    s_mid = s + (u * (M0[0] - s0) if align_root else 0.0)
    tw = u * tw0
    foil = LY.naca44(float(blade.t_c(r0)), blade.te_mm / c0, blade.camber_m, blade.camber_p, nps)
    rect = LY.rounded_rect(2 * lay.knuckle_r / lay.knuckle_l, rect_exp, nps)
    sgn = np.r_[np.ones(nps), -np.ones(nps - 1)]
    # blunt end of the knuckle rectangle: te_mm opening concentrated at its TE end
    # (x^RECT_TE_EXP), so the shank does not grow outside the layout envelope elsewhere
    rect_ramp = sgn * 0.5 * te_mm * rect[:, 0] ** RECT_TE_EXP
    q = np.clip(s / ROOT_SHRINK_LEN, 0, 1)
    ys = 1.0 - root_shrink * (1 - q * q * (3 - 2 * q))          # thickness only, not the TE
    X = (1 - u)[:, None] * (0.5 - rect[:, 0])[None, :] * width + u[:, None] * (0.5 - foil[:, 0])[None, :] * c0
    Y = ((1 - u)[:, None] * (rect[:, 1][None, :] * lay.knuckle_l * ys[:, None] + rect_ramp[None, :])
         + u[:, None] * foil[:, 1][None, :] * c0)
    L, U = LY.section_frame(sl, tw)
    M = np.stack([s_mid, c_mid, n_mid], axis=-1)
    pts = M[:, None, :] + X[..., None] * L[:, None, :] + Y[..., None] * U[:, None, :]
    return pts, s


def aero_main_radii() -> np.ndarray:
    """Layout station radii (layout.blade_surface grid) up to the tip zone, plus its start."""
    blade, _, _ = context()
    r = np.linspace(blade.r_root, blade.r_tip, 41)
    r_tz = blade.r_tip - TIP_ZONE
    return r[r < r_tz - 0.3]


@lru_cache(maxsize=1)
def section_grid() -> dict:
    """All loft sections root -> tip: (Nv, Np, 3) blade-frame grid plus bookkeeping."""
    blade, lay, ax = context()
    sh, s_sh = shank_sections()
    r_main = aero_main_radii()
    st = LY.AeroStations(blade, r_main)
    ld = ax.lead_at(r_main)
    aero = st.points(LY.s_of_r(r_main, ld, lay.hinge_r), ld, ax.aft_at(r_main), ax.slope_at(r_main, ld))
    tip, r_tip, meta = tip_sections()
    join_err = float(np.abs(sh[-1] - aero[0]).max())
    # layout.shank_sections ends at s = s0 but the aero root section is centred at
    # s0 + (pitch_axis - 0.5) c0 L0_s (0.63 mm nearer the hinge): keep the shank sections
    # that lie at least SHANK_JOIN_GAP before the aero root plane, the skin bridges the rest
    nrm = np.cross(aero[0][len(aero[0]) // 2] - aero[0][0], aero[0][-1] - aero[0][0])
    nrm /= np.linalg.norm(nrm)
    o = aero[0].mean(axis=0)
    if np.dot(aero[1].mean(axis=0) - o, nrm) < 0:
        nrm = -nrm
    dist = np.array([float(np.dot(sec.mean(axis=0) - o, nrm)) for sec in sh])
    keep = np.nonzero(dist < -SHANK_JOIN_GAP)[0]
    gap = -dist[keep[-1]]
    G = np.concatenate([sh[keep], aero, tip], axis=0)
    root_plane = (o, nrm)
    zone = np.array(["shank"] * len(keep) + ["aero"] * len(aero) + ["tip"] * len(tip))
    return {"G": G, "zone": zone, "r_main": r_main, "r_tip_zone": r_tip, "tip_meta": meta,
            "shank_last_vs_aero_root_mm": join_err, "shank_sections_used": int(len(keep)),
            "shank_gap_to_aero_root_mm": float(gap), "root_plane": root_plane}


# ==========================================================================
# Solids
# ==========================================================================

def _chord_params(G: np.ndarray, axis: int) -> np.ndarray:
    """Normalised chord-length parameters along ``axis`` of a (nv, nu, 3) grid, averaged
    over the other direction (the same for every row, so every row gets the same knots)."""
    d = np.linalg.norm(np.diff(G, axis=axis), axis=-1)
    cum = np.concatenate([np.zeros_like(np.take(d, [0], axis=axis)), np.cumsum(d, axis=axis)], axis=axis)
    cum = cum / np.take(cum, [-1], axis=axis)
    t = cum.mean(axis=1 - axis)
    t[0], t[-1] = 0.0, 1.0
    return t


def _knots_mults(t: np.ndarray):
    from OCP.collections import Array1_double, Array1_int
    vals, mults = np.unique(np.round(t, 14), return_counts=True)
    K, M = Array1_double(1, len(vals)), Array1_int(1, len(vals))
    for k, (v, m) in enumerate(zip(vals, mults)):
        K.SetValue(k + 1, float(v))
        M.SetValue(k + 1, int(m))
    return K, M


def _bspline_surface(G: np.ndarray):
    """Skinned C2 bicubic B-spline surface through every grid point (U around the section,
    V along the span). Every section is first interpolated with the SAME parameters (so
    the same knots), then the pole columns are interpolated along the span: each section
    is then an exact iso-curve of the surface and stays in its own plane (OCC's
    PointsToBSplineSurface end conditions bend the end sections out of plane)."""
    from OCP.Geom import Geom_BSplineSurface
    from scipy.interpolate import make_interp_spline
    nv, nu = G.shape[:2]
    tu, tv = _chord_params(G, 1), _chord_params(G, 0)
    su = make_interp_spline(tu, np.transpose(G, (1, 0, 2)), k=3)          # coeffs (nu, nv, 3)
    sv = make_interp_spline(tv, np.transpose(su.c, (1, 0, 2)), k=3)       # coeffs (nv, nu, 3)
    poles = sv.c                                                          # [j_v, i_u]
    arr = TColgp_Array2OfPnt(1, poles.shape[1], 1, poles.shape[0])
    for j in range(poles.shape[0]):
        for i in range(poles.shape[1]):
            arr.SetValue(i + 1, j + 1, gp_Pnt(*map(float, poles[j, i])))
    uk, um = _knots_mults(su.t)
    vk, vm = _knots_mults(sv.t)
    return Geom_BSplineSurface(arr, uk, vk, um, vm, 3, 3, False, False)


def _edge(curve) -> Edge:
    return Edge(BRepBuilderAPI_MakeEdge(curve).Edge())


@lru_cache(maxsize=1)
def loft_surface():
    return _bspline_surface(section_grid()["G"])


@lru_cache(maxsize=1)
def span_solid() -> Solid:
    """Shank + aero + tip as one closed solid: the interpolated skin, the ruled blunt-TE
    strip and two planar end caps (root cap buried in the knuckle)."""
    S = loft_surface()
    u0, u1, v0, v1 = S.Bounds()
    skin = Face(BRepBuilderAPI_MakeFace(S, 1e-7).Face())
    e_up, e_lo = _edge(S.UIso(u0)), _edge(S.UIso(u1))
    te = Face.make_surface_from_curves(e_up, e_lo)
    caps = []
    for v in (v0, v1):
        e = _edge(S.VIso(v))
        caps.append(Face(Wire([e, Edge.make_line(e.position_at(1), e.position_at(0))])))
    shell = Shell([skin, te, *caps])
    sol = Solid(shell)
    if sol.volume < 0:
        sol = Solid(sol.wrapped.Reversed())
    return sol


def volume(shape: Shape) -> float:
    """Volume (mm3) with adaptive integration: OCC's default Gauss rule reads the B-spline
    skin ~0.9 % low (checked against a 0.005 mm mesh)."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape.wrapped, pr, 1e-7)
    return float(pr.Mass())


def _cyl(p0, p1, r: float) -> Solid:
    return C.cyl(p0, p1, r)


def _cone(p0, p1, r0: float, r1: float) -> Solid:
    """Cone frustum from p0 (radius r0) to p1 (radius r1)."""
    a, b = Vector(*p0), Vector(*p1)
    return Solid.make_cone(r0, r1, (b - a).length, Plane(origin=a, z_dir=(b - a).normalized()))


def knuckle_heel() -> Solid:
    """Knuckle cylinder about the pin (local Y) united with the heel tab."""
    _, lay, _ = context()
    R, Lk = P.KNUCKLE_D / 2, P.KNUCKLE_L
    kn = _cyl((0, -Lk / 2, 0), (0, Lk / 2, 0), R)    # end chamfers cut later (_end_chamfer_cutter)
    reach, t, w = P.HEEL_REACH, P.HEEL_T, P.HEEL_W
    heel = Solid.make_box(reach, w, t, Plane(origin=(-reach, -w / 2, -t / 2)))
    e = heel.edges()
    tip_top = e.filter_by(Axis.Y).filter_by(lambda ed: ed.center().X < -reach + 1e-6 and ed.center().Z > 0)
    tip_corners = e.filter_by(Axis.Z).filter_by(lambda ed: ed.center().X < -reach + 1e-6)
    heel = heel.fillet(HEEL_TIP_FILLET, list(tip_top) + list(tip_corners))
    kh = kn.fuse(heel).clean()
    # OCC propagates a fillet along tangent edge chains, so the order matters: the +n root
    # fillet goes first (its ends meet the heel flanks at sharp corners and stop there),
    # then the +n long edges; the -n stop-face perimeter and its knuckle junction stay sharp
    # (flat stop on the hub block, magnet seat rim 0.39 mm from the heel tip)
    s_j = -math.sqrt(R**2 - (t / 2) ** 2)
    root = [ed for ed in kh.edges() if ed.geom_type == GeomType.LINE
            and abs(ed.center().X - s_j) < 1e-3 and ed.center().Z > 0]
    if root and HEEL_ROOT_FILLET > 0:
        kh = kh.fillet(HEEL_ROOT_FILLET, root)
    sides = [ed for ed in kh.edges() if ed.geom_type == GeomType.LINE and abs(ed.tangent_at(0.5).X) > 0.999
             and abs(abs(ed.center().Y) - w / 2) < 1e-6 and abs(ed.center().Z - t / 2) < 1e-6]
    if sides and HEEL_SIDE_FILLET > 0:
        kh = kh.fillet(HEEL_SIDE_FILLET, sides)
    stray = [f for f in kh.faces() if f.geom_type not in (GeomType.PLANE, GeomType.CYLINDER)
             and f.bounding_box().min.Z < -t / 2 + 0.3 and f.bounding_box().max.X > -reach + 1.5]
    stop = [f for f in kh.faces() if f.geom_type == GeomType.PLANE and abs(f.center().Z + t / 2) < 1e-6
            and f.bounding_box().max.X < 0]
    assert kh.is_valid and not stray, "a heel fillet propagated onto the stop face side"
    assert stop and stop[0].bounding_box().size.Y <= w + 1e-6, "heel stop face no longer sharp-edged"
    return kh


def _end_chamfer_cutter() -> Solid:
    """45 deg rings that chamfer the knuckle's outer end edges. Cut after the union: OCC's
    boolean of the shank skin with a pre-chamfered knuckle fails silently (no overlap)."""
    from build123d import BuildLine, BuildPart, BuildSketch, Polyline, make_face, revolve
    R, Lk, ch = P.KNUCKLE_D / 2, P.KNUCKLE_L, KNUCKLE_CHAMFER
    y0, y1 = Lk / 2 - ch, Lk / 2 + 0.6
    pts = [(y0, R), (y1, R - (y1 - y0)), (y1, R + 1.0), (y0, R + 1.0)]
    with BuildPart() as bp:
        with BuildSketch(Plane.XY):
            with BuildLine():
                Polyline(*[(y, r) for y, r in pts], (y0, R), close=False)
            make_face()
        revolve(axis=Axis.X)
    ring = bp.part.rotate(Axis.Z, 90)          # local X (revolve axis) -> local Y (pin axis)
    return ring.fuse(ring.mirror(Plane.XZ))


def _bore_cutter() -> Solid:
    Lk, r = P.KNUCKLE_L, BORE_D / 2
    cut = _cyl((0, -Lk / 2 - 1, 0), (0, Lk / 2 + 1, 0), r)
    for sg in (-1, 1):
        y_face = sg * Lk / 2
        cut = cut.fuse(_cone((0, y_face - sg * BORE_CHAMFER, 0), (0, y_face + sg * 0.5, 0), r,
                             r + BORE_CHAMFER + 0.5))
    return cut


def _seat_cutter(center, axis, t_face: float, chamfer: float) -> Solid:
    """Magnet seat along unit ``axis`` (pointing OUT of the blade) whose face sits at
    parameter t_face (mm along axis from ``center``): bore MAG_SEAT_D from the seat bottom
    (t_face - MAG_SEAT_DEPTH) to 2 mm outside, plus an optional 45 deg mouth chamfer."""
    c, a = np.asarray(center, float), np.asarray(axis, float)
    r = MAG_SEAT_D / 2
    bot = c + (t_face - MAG_H) * a           # magnet face at t_face, its back on the seat floor
    out = c + (t_face + 2.0) * a
    cut = _cyl(bot, out, r)
    if chamfer > 0:
        p0 = c + (t_face - chamfer) * a
        p1 = c + (t_face + 1.5) * a
        cut = cut.fuse(_cone(p0, p1, r, r + chamfer + 1.5))
    return cut


def _crossing(sol: Solid, p_in, p_out, n_iter: int = 40) -> np.ndarray:
    """Boundary point of ``sol`` on the segment p_in (inside) -> p_out (outside)."""
    a, b = np.asarray(p_in, float), np.asarray(p_out, float)
    for _ in range(n_iter):
        m = 0.5 * (a + b)
        if sol.is_inside(Vector(*m), tolerance=1e-7):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def _skin_normal(p: np.ndarray) -> np.ndarray:
    """Unit normal of the loft skin at its point nearest p (sign arbitrary)."""
    from OCP.gp import gp_Vec
    S = loft_surface()
    pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, p)), S)
    uu, vv = pr.LowerDistanceParameters()
    pt, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
    S.D1(uu, vv, pt, du, dv)
    n = du.Crossed(dv)
    v = np.array([n.X(), n.Y(), n.Z()])
    return v / np.linalg.norm(v)


def stowed_seat_geometry(sol: Solid) -> dict:
    """Stowed magnet seat. The DETENTS centre lies on a RADIAL axis, but there the shank's
    inner face is inclined ~50 deg to the radial direction, so a radial 4 mm seat cannot be
    flush. The seat is centred where the DETENTS axis meets the blade's inner face and runs
    along the face normal; its face sits 0.1 below the lowest rim point."""
    c = np.asarray(P.DETENTS["stowed_blade_magnet_blade_frame"], float)
    a = np.asarray(stowed_axis_detents(), float)
    p = _crossing(sol, c - 2 * a, c + 4 * a)
    t_axis = float(np.dot(p - c, a))
    n = _skin_normal(p)
    if sol.is_inside(Vector(*(p + 0.05 * n)), tolerance=1e-7):
        n = -n                                               # outward
    e1 = np.cross(n, [1.0, 0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    t_rim = []
    for ang in np.linspace(0, 2 * math.pi, 24, endpoint=False):
        q = (MAG_SEAT_D / 2) * (math.cos(ang) * e1 + math.sin(ang) * e2)
        t_rim.append(float(np.dot(_crossing(sol, p + q - 1.5 * n, p + q + 1.5 * n) - p, n)))
    t_rim = np.array(t_rim)
    t_face = float(t_rim.min()) - (MAG_SEAT_DEPTH - MAG_H)
    seated = p + (t_face - MAG_H / 2) * n
    return {"mouth_blade_frame_mm": p.tolist(), "axis_out_blade_frame": n.tolist(),
            "detents_axis_t_at_face_mm": t_axis,
            "angle_face_normal_vs_detents_axis_deg": float(math.degrees(math.acos(min(1.0, abs(float(n @ a)))))),
            "rim_height_range_mm": [float(t_rim.min()), float(t_rim.max())], "t_face_mm": t_face,
            "seated_magnet_centre_blade_frame_mm": seated.tolist(),
            "seated_vs_detents_centre_mm": float(np.linalg.norm(seated - c))}


def stowed_axis_detents() -> np.ndarray:
    """Stowed magnet axis in the blade frame (pointing out of the inner face), from
    data/layout_rotor.json detents (params.py DETENTS carries only the centre)."""
    d = json.loads(LY.DATA_JSON.read_text())["detents"]["stowed"]
    a = np.asarray(d["blade_magnet_axis_blade_frame"], float)
    return a / np.linalg.norm(a)


@lru_cache(maxsize=1)
def _build() -> tuple[Solid, dict]:
    t0 = time.time()
    span = span_solid()
    kh = knuckle_heel()
    body = span.fuse(kh).clean()
    info = {"t_span_union_s": round(time.time() - t0, 1),
            "overlap_span_knuckle_mm3": round(volume(span) + volume(kh) - volume(body), 1)}
    body = body.cut(_end_chamfer_cutter())
    body = body.cut(_bore_cutter())
    # deployed seat: heel stop face is the plane n = -HEEL_T/2, axis -n
    cd = np.asarray(P.DETENTS["deployed_blade_magnet_blade_frame"], float)
    ad = np.array([0.0, 0.0, -1.0])
    t_face_d = (P.HEEL_T / 2 + cd[2]) - (MAG_SEAT_DEPTH - MAG_H)   # face plane at t = 2.5 - 1.3 = 1.2
    body = body.cut(_seat_cutter(cd, ad, t_face_d, HEEL_SEAT_CHAMFER))
    sg = stowed_seat_geometry(body)
    body = body.cut(_seat_cutter(sg["mouth_blade_frame_mm"], sg["axis_out_blade_frame"], sg["t_face_mm"],
                                 MAG_MOUTH_CHAMFER))
    cs = np.asarray(P.DETENTS["stowed_blade_magnet_blade_frame"], float)
    body = body.clean()
    sols = body.solids()
    if len(sols) == 1:
        body = sols[0]
    info.update(deployed_seat={"centre_detents": cd.tolist(), "axis": ad.tolist(), "t_face_mm": t_face_d},
                stowed_seat={"centre_detents": cs.tolist(), "axis_detents": stowed_axis_detents().tolist(), **sg},
                t_build_s=round(time.time() - t0, 1))
    return body, info


def blade_frame_solid() -> Solid:
    """One blade in the blade frame (X span, Y chord along the pin, Z thickness)."""
    return _build()[0]


def _loc(T: np.ndarray) -> Location:
    tr = gp_Trsf()
    tr.SetValues(*(float(T[i, j]) for i in range(3) for j in range(4)))
    return Location(tr)


def blade_transform(i: int, fold_deg: float) -> np.ndarray:
    _, lay, _ = context()
    return LY.blade_to_world(fold_deg, lay.psi(i), lay)


def blade_world(i: int, fold_deg: float = P.FOLD_DEG) -> Solid:
    """Blade i (psi = 120 i) in the world frame at a fold angle (0 deployed, FOLD_DEG stowed)."""
    return blade_frame_solid().moved(_loc(blade_transform(i, fold_deg)))


def pivot(i: int) -> dict:
    """Hinge pivot for blade i: rotating the STOWED blade by deployed_deg (right-hand rule)
    about ``axis`` through ``origin`` gives the deployed blade. Folding is +fold about e_t,
    so deploying is +FOLD_DEG about -e_t."""
    _, lay, _ = context()
    _, _, e_t = LY.rotor_axes(lay.psi(i))
    return {"origin": [float(v) for v in LY.hinge_point(lay, lay.psi(i))],
            "axis": [float(v) + 0.0 for v in -e_t], "deployed_deg": float(P.FOLD_DEG)}


# blade frame (s, c, n) -> print (x, y, z) = (s, -n, c): chord vertical, TE side down,
# knuckle end face on the bed, pin bore vertical
PRINT_R = np.array([[1.0, 0, 0, 0], [0, 0, -1.0, 0], [0, 1.0, 0, 0]])


def print_solid() -> Solid:
    return C.on_bed(blade_frame_solid().moved(_loc(PRINT_R)))


PRINT_NOTE = ("Print chord-vertical: trailing-edge side down, the knuckle's TE-side end face on the bed, "
              "pin bore vertical (round, no teardrop), layers in the span/thickness plane so flapwise "
              "bending, the knuckle lug and the heel stop load all run inside the layers. Supports: tree, "
              "build plate only, under the trailing edge (it floats {te_lo:.1f} to {te_hi:.1f} mm above the "
              "bed), the heel underside and the shank's TE ramp; nothing else overhangs > 45 deg. Magnet seats "
              "print as horizontal 4.1 mm holes (bridged). PLA or PETG, 0.2 mm layers, 4 walls, 40 % gyroid, "
              "brim on the knuckle.")


def parts() -> list[dict]:
    """Three identical blades, WORLD STOWED pose, with hinge pivots (common.py contract)."""
    _, lay, _ = context()
    pr = print_solid()
    out = []
    for i in range(lay.n_blades):
        _, e_r, _ = LY.rotor_axes(lay.psi(i))
        out.append({
            "name": f"blade_{i}", "shape": blade_world(i, P.FOLD_DEG), "print": pr,
            "material": "blade", "role": "body", "rotor": True, "qty": 1, "bought": False,
            "explode": [float(v) for v in (18.0 * e_r + np.array([12.0, 0, 0]))],
            "pivot": pivot(i),
            "notes": "Blade %d of 3 (identical prints). " % i + PRINT_NOTE.format(**_print_te_range()),
        })
    return out


def preview_parts() -> list[dict]:
    """Blades plus 3 mm pins (steel) for the visual check."""
    _, lay, _ = context()
    ps = parts()
    out = [{k: v for k, v in p.items() if k not in ("print",)} for p in ps]
    for i in range(lay.n_blades):
        H = LY.hinge_point(lay, lay.psi(i))
        _, _, e_t = LY.rotor_axes(lay.psi(i))
        pin = _cyl(H - P.PIN_L / 2 * e_t, H + P.PIN_L / 2 * e_t, P.PIN_D / 2)
        out.append({"name": f"pin_{i}", "shape": pin, "material": "steel", "role": "body", "rotor": True})
    return out


def preview_root() -> list[dict]:
    """Close-up target: blade root (s < 22) in the blade frame, with its pin."""
    root = blade_frame_solid() & Solid.make_box(30, 40, 40, Plane(origin=(-9, -15, -10)))
    pin = _cyl((0, -P.PIN_L / 2, 0), (0, P.PIN_L / 2, 0), P.PIN_D / 2)
    return [{"name": "root", "shape": root, "material": "blade"},
            {"name": "pin", "shape": pin, "material": "steel"}]


def preview_heel() -> list[dict]:
    """Close-up target: knuckle and heel (s < 6) in the blade frame."""
    part = blade_frame_solid() & Solid.make_box(14, 20, 20, Plane(origin=(-9, -10, -10)))
    return [{"name": "heel", "shape": part, "material": "blade"}]


def preview_tip() -> list[dict]:
    """Close-up target: last 20 mm of span in the blade frame."""
    tip = blade_frame_solid() & Solid.make_box(25, 40, 40, Plane(origin=(62, -15, 0)))
    return [{"name": "tip", "shape": tip, "material": "blade"}]


def preview_print() -> list[dict]:
    """The print orientation on a 120 mm plate patch (z = 0)."""
    plate = Solid.make_box(120, 50, 0.6, Plane(origin=(-60, -25, -0.6)))
    return [{"name": "blade_print", "shape": print_solid(), "material": "blade"},
            {"name": "plate", "shape": plate, "material": "stand", "role": "body"}]


@lru_cache(maxsize=1)
def _print_te_range() -> dict:
    """Height of the trailing edge above the bed in the print orientation."""
    G = section_grid()["G"]
    te = np.vstack([G[:, 0], G[:, -1]])
    zmin = -P.KNUCKLE_L / 2
    z = te[:, 1] - zmin
    return {"te_lo": float(z.min()), "te_hi": float(z.max())}


# ==========================================================================
# Checks
# ==========================================================================

def _project_surface(pts: np.ndarray) -> np.ndarray:
    """Distance of points to the loft skin (exact projection on the B-spline surface)."""
    S = loft_surface()
    out = np.full(len(pts), np.nan)
    for k, p in enumerate(pts):
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, p)), S)
        if pr.NbPoints():
            out[k] = pr.LowerDistance()
    return out


def _dist_to_faces(sol: Shape, pts: np.ndarray) -> np.ndarray:
    """Exact distance from points to the boundary of a solid (BRepExtrema on its shell)."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    shell = sol.shells()[0].wrapped if sol.shells() else sol.wrapped
    out = []
    for p in pts:
        v = BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, p))).Vertex()
        d = BRepExtrema_DistShapeShape(v, shell)
        out.append(d.Value() if d.IsDone() else np.nan)
    return np.array(out)


def _mesh(shape: Shape, tol: float = STL_TOL):
    import trimesh
    verts, faces = shape.tessellate(tol, 0.1)
    return trimesh.Trimesh(np.array([(v.X, v.Y, v.Z) for v in verts]), np.array(faces), process=True)


def _cyl_faces(sol: Shape, radius: float, tol: float = 1e-3) -> list:
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Cylinder
    found = []
    for f in sol.faces():
        ad = BRepAdaptor_Surface(f.wrapped)
        if ad.GetType() == GeomAbs_Cylinder:
            cy = ad.Cylinder()
            if abs(cy.Radius() - radius) < tol:
                ax = cy.Axis()
                loc, d = ax.Location(), ax.Direction()
                found.append({"face": f, "r": cy.Radius(), "loc": np.array([loc.X(), loc.Y(), loc.Z()]),
                              "dir": np.array([d.X(), d.Y(), d.Z()])})
    return found


def _seat_measure(sol: Solid, centre, axis) -> dict:
    """Measured magnet seat: its cylinder face, the seat-bottom position along the axis
    and the seated magnet centre (pressed to the bottom), vs the DETENTS centre."""
    c, a = np.asarray(centre, float), np.asarray(axis, float)
    best = None
    for cf in _cyl_faces(sol, MAG_SEAT_D / 2):
        off = cf["loc"] - c
        lateral = float(np.linalg.norm(off - np.dot(off, cf["dir"]) * cf["dir"]))
        par = abs(abs(float(np.dot(cf["dir"], a))) - 1.0)
        if lateral < 0.05 and par < 1e-6:
            best = cf
    if best is None:
        return {"found": False}
    ts = [float(np.dot(np.array([v.X, v.Y, v.Z]) - c, a)) for v in best["face"].vertices()]
    bb = best["face"].bounding_box()
    corners = np.array([[x, y, z] for x in (bb.min.X, bb.max.X) for y in (bb.min.Y, bb.max.Y)
                        for z in (bb.min.Z, bb.max.Z)])
    t_lo = float(min(np.dot(corners - c, a))) if not ts else min(ts)
    # seat bottom: probe along the axis from inside the seat toward the floor
    p_air = c + (t_lo + 0.3) * a
    t_bot = float(np.dot(_crossing(sol, c + (t_lo - 1.0) * a, p_air) - c, a)) if not sol.is_inside(Vector(*p_air)) else t_lo
    seated = c + (t_bot + MAG_H / 2) * a
    off = best["loc"] - c
    lateral = float(np.linalg.norm(off - np.dot(off, best["dir"]) * best["dir"]))
    return {"found": True, "seat_d_mm": 2 * best["r"], "axis_lateral_offset_mm": lateral,
            "seat_bottom_t_mm": t_bot, "seated_magnet_centre_blade_frame_mm": seated.tolist(),
            "seated_vs_detents_mm": float(np.linalg.norm(seated - c)),
            "seated_vs_detents_along_axis_mm": float(t_bot + MAG_H / 2)}


def checks(write: bool = True) -> dict:
    """Numerical checks listed in the task, returned (and written) as a dict."""
    import trimesh
    t0 = time.time()
    blade, lay, ax = context()
    sol, info = _build()
    res: dict = {"module": "hydroegg.cad.blade", "generated": time.strftime("%Y-%m-%d %H:%M"),
                 "build": info}
    # ---- consistency of the inputs
    ax_mem = LY.solve_axis(blade, lay, n=41, iters=3)
    rr = np.linspace(blade.r_root, blade.r_tip, 300)
    res["inputs"] = {
        "blade_source": blade.source,
        "axis_json_vs_resolved_max_mm": float(max(np.abs(ax.aft_at(rr) - ax_mem.aft_at(rr)).max(),
                                                  np.abs(ax.lead_at(rr) - ax_mem.lead_at(rr)).max())),
        "params_vs_layout_chosen": {
            "HINGE_X": [P.HINGE_X, lay.hinge_x], "HINGE_R": [P.HINGE_R, lay.hinge_r],
            "FOLD_DEG": [P.FOLD_DEG, lay.fold_deg], "KNUCKLE_D": [P.KNUCKLE_D, 2 * lay.knuckle_r],
            "KNUCKLE_L": [P.KNUCKLE_L, lay.knuckle_l], "HEEL_REACH": [P.HEEL_REACH, lay.heel_reach],
            "HEEL_T": [P.HEEL_T, lay.heel_t], "HEEL_W": [P.HEEL_W, lay.heel_w],
            "SHANK_CONE_DEG": [P.SHANK_CONE_DEG, lay.shank_cone_deg]},
        "params_match_layout": bool(all(abs(a - b) < 1e-9 for a, b in [
            (P.HINGE_X, lay.hinge_x), (P.HINGE_R, lay.hinge_r), (P.FOLD_DEG, lay.fold_deg),
            (P.KNUCKLE_D, 2 * lay.knuckle_r), (P.KNUCKLE_L, lay.knuckle_l), (P.HEEL_REACH, lay.heel_reach),
            (P.HEEL_T, lay.heel_t), (P.HEEL_W, lay.heel_w), (P.SHANK_CONE_DEG, lay.shank_cone_deg)])),
    }
    # reproduction of the layout generators by this module's variants
    sh_l, _ = LY.shank_sections(blade, lay, ax, SHANK_N)
    sh_m, _ = shank_sections(P.KNUCKLE_L, 0.0, 0.0, False, rect_exp=4.0)
    naca_err = max(float(np.abs(section_unit(float(blade.t_c(r)), blade.te_mm, float(blade.chord(r)), 1.0)
                                - LY.naca44(float(blade.t_c(r)), blade.te_mm / float(blade.chord(r)))).max())
                   for r in (30.0, 63.5, 97.0))
    sg = section_grid()
    res["generators"] = {"shank_sections_reproduce_layout_max_mm": float(np.abs(sh_l - sh_m).max()),
                         "section_unit_reproduces_naca44_max": naca_err,
                         "layout_shank_last_section_vs_aero_root_mm": sg["shank_last_vs_aero_root_mm"],
                         "shank_sections_used": sg["shank_sections_used"],
                         "shank_gap_to_aero_root_plane_mm": sg["shank_gap_to_aero_root_mm"],
                         "n_sections": int(sg["G"].shape[0]), "points_per_section": int(sg["G"].shape[1])}
    # ---- solid
    chk = C.check(sol, "blade (blade frame)")
    chk["volume_cm3_adaptive"] = round(volume(sol) / 1000, 4)
    res["solid"] = chk
    # ---- CAD vs layout surface
    S = LY.blade_surface(blade, lay, ax)
    r_tz = blade.r_tip - TIP_ZONE
    rng = np.random.default_rng(7)
    on = S.aero[S.aero_r <= r_tz + 1e-9]
    on_pts = on.reshape(-1, 3)[rng.choice(on.shape[0] * on.shape[1], 240, replace=False)]
    r_off = rng.uniform(blade.r_root, r_tz, 40)
    st = LY.AeroStations(blade, r_off)
    ld = ax.lead_at(r_off)
    off = st.points(LY.s_of_r(r_off, ld, lay.hinge_r), ld, ax.aft_at(r_off), ax.slope_at(r_off, ld))
    off_pts = off[:, ::3].reshape(-1, 3)
    d_on = _project_surface(on_pts)
    d_off = _project_surface(off_pts)
    d_faces = _dist_to_faces(sol, np.vstack([on_pts[:40], off_pts[:40]]))
    tipz = S.aero[S.aero_r > r_tz + 1e-9].reshape(-1, 3)
    d_tip = _project_surface(tipz)
    sh_pts = S.shank.reshape(-1, 3)
    out_kn = np.hypot(sh_pts[:, 0], sh_pts[:, 2]) > P.KNUCKLE_D / 2 + 0.3
    d_sh = _project_surface(sh_pts[out_kn])
    res["surface_vs_layout"] = {
        "aero_station_samples": int(len(on_pts)), "aero_station_max_dev_mm": float(np.nanmax(d_on)),
        "aero_between_station_samples": int(len(off_pts)),
        "aero_between_station_max_dev_mm": float(np.nanmax(d_off)),
        "aero_final_solid_faces_max_dev_mm (80 samples, BRepExtrema)": float(np.nanmax(d_faces)),
        "aero_zone_r_mm": [float(blade.r_root), r_tz],
        "tip_zone_layout_points": int(len(tipz)),
        "tip_zone_max_dev_mm (intentional rounding, material only removed)": float(np.nanmax(d_tip)),
        "shank_points_outside_knuckle": int(out_kn.sum()),
        "shank_max_dev_mm (intentional: 13.2 mm flanks, rounder corners, blunt TE, bridge to the aero root)":
            float(np.nanmax(d_sh)),
        "shank_mean_dev_mm": float(np.nanmean(d_sh)),
    }
    # tip sections inside the nominal section outline (2D, section plane)
    from matplotlib.path import Path as MPath
    worst = -np.inf
    for m in sg["tip_meta"]:
        nom = section_unit(float(blade.t_c(m["r"])), blade.te_mm, m["c"], 1.0)
        poly = MPath(nom)
        xy = m["xy"]
        inside = poly.contains_points(xy, radius=1e-9)
        dmin = np.min(np.hypot(*(xy[:, None, :] - nom[None, :, :]).transpose(2, 0, 1)), axis=1) * m["c"]
        outside = ~inside & (dmin > 1e-4)          # on-boundary points (u = 0 section) are inside
        if outside.any():
            worst = max(worst, float(dmin[outside].max()))
    res["tip"] = {"zone_r_mm": [r_tz, float(blade.r_tip)], "chord_end_mm": sg["tip_meta"][-1]["chord_mm"],
                  "t_max_end_mm": sg["tip_meta"][-1]["t_max_mm"],
                  "sections_inside_nominal": bool(worst == -np.inf),
                  "max_outside_nominal_mm": None if worst == -np.inf else worst}
    # ---- trailing edge thickness (ruled TE strip width) and tip cap
    Sx = loft_surface()
    u0, u1, v0, v1 = Sx.Bounds()
    vv = np.linspace(v0, v1, 400)
    te_w = [Sx.Value(u0, v).Distance(Sx.Value(u1, v)) for v in vv]
    capxy = sg["tip_meta"][-1]["xy"] * sg["tip_meta"][-1]["c"]
    nps = LY.N_PER_SIDE
    up, lo = capxy[:nps][::-1], capxy[nps - 1:]
    res["trailing_edge"] = {"min_mm": float(min(te_w)), "max_mm": float(max(te_w)), "MIN_TE": P.MIN_TE,
                            "pass": bool(min(te_w) >= P.MIN_TE - 1e-3),
                            "note": "exactly MIN_TE at every section; the cubic skin dips < 0.001 between them",
                            "tip_cap_chord_mm": float(np.ptp(capxy[:, 0])),
                            "tip_cap_max_thickness_mm": float(np.max(up[:, 1] - np.interp(up[:, 0], lo[:, 0], lo[:, 1]))),
                            "tip_cap_te_mm": float(capxy[0, 1] - capxy[-1, 1])}
    # ---- knuckle bore
    bores = _cyl_faces(sol, BORE_D / 2)
    res["knuckle_bore"] = {
        "found_faces": len(bores), "d_mm": [round(2 * b["r"], 4) for b in bores],
        "axis_vs_pin_axis_deg": [float(math.degrees(math.acos(min(1.0, abs(float(b["dir"][1])))))) for b in bores],
        "axis_offset_from_pin_mm": [float(np.hypot(b["loc"][0], b["loc"][2])) for b in bores],
        "pin_d_mm": P.PIN_D, "diametral_clearance_mm": BORE_D - P.PIN_D,
        "knuckle_od_mm": P.KNUCKLE_D, "knuckle_len_mm": P.KNUCKLE_L,
        "wall_around_bore_mm": P.KNUCKLE_D / 2 - BORE_D / 2}
    # ---- magnet seats
    res["magnet_seats"] = {
        "seat_d_mm": MAG_SEAT_D, "seat_depth_mm": MAG_SEAT_DEPTH,
        "deployed": _seat_measure(sol, P.DETENTS["deployed_blade_magnet_blade_frame"], [0, 0, -1]),
        "stowed": {**_seat_measure(sol, info["stowed_seat"]["seated_magnet_centre_blade_frame_mm"],
                                   info["stowed_seat"]["axis_out_blade_frame"]),
                   **{k: v for k, v in info["stowed_seat"].items() if k not in ("centre_detents",)}},
        "heel_tip_wall_mm": float(P.HEEL_REACH - (-P.DETENTS["deployed_blade_magnet_blade_frame"][0])
                                  - MAG_SEAT_D / 2),
        "note": ("seat depth MAGNET H + 0.1 from the contact face (the lowest rim point on the curved "
                 "stowed face), so a seated magnet sits 0.1 below its face; layout magnet_recess is 0.2, "
                 "so the seated centre is 0.1 nearer the face than the DETENTS centre")}
    # ---- mass (adaptive volume integration) and where it differs from the layout estimate
    from build123d import Keep
    vol = volume(sol)
    span, kh = span_solid(), knuckle_heel()
    o, nrm = sg["root_plane"]
    aero_part = span.split(Plane(origin=tuple(o), z_dir=tuple(nrm)), keep=Keep.TOP)
    v_aero = volume(aero_part)
    v_overlap = volume(span & kh)
    v_shank = volume(span) - v_aero - v_overlap
    lm = LY.blade_mass_props(lay, S)
    res["mass"] = {
        "volume_mm3": vol, "mass_g_PLA_solid": vol * PLA_G_MM3, "density_g_cm3": PLA_G_MM3 * 1e3,
        "volume_occ_default_rule_mm3": float(sol.volume),
        "layout_blade_mass_g": lm["mass_g"], "layout_volume_mm3": lm["volume_mm3"],
        "delta_pct": 100 * (vol / lm["volume_mm3"] - 1),
        "cad_parts_mm3": {"aero_r30_to_tip": v_aero, "shank_outside_knuckle": v_shank,
                          "knuckle_heel_net_of_bore_seats_chamfers": vol - v_aero - v_shank},
        "layout_parts_mm3": {"aero": lm["volume_aero_mm3"], "shank_loft": lm["volume_shank_mm3"],
                             "knuckle_minus_3mm_pin": lm["volume_knuckle_mm3"], "heel": lm["volume_heel_mm3"]},
        "note": ("aero matches; the layout's shank estimate integrates section area along the centroid path "
                 "(sections tilted up to 55 deg to it) and includes the part inside the knuckle "
                 "(%.0f mm3 here), which its knuckle term counts again; the CAD shank also has 13.2 mm "
                 "flanks; the CAD bore is 3.3 (layout 3.0) and adds magnet seats and chamfers" % v_overlap)}
    # ---- stowed blade vs pocket floor (blade 0, world)
    pm = LY.pocket_map(blade, lay, S)
    tb = json.loads(LY.DATA_JSON.read_text())["pocket_depth"]
    mw = _mesh(blade_world(0, P.FOLD_DEG), 0.03)
    W = np.asarray(mw.vertices)
    x, rho, az = LY.cyl(W, lay.psi(0))
    stat = x < lay.static_tail_x - 1e-6
    fl = pm.lookup(x[stat], az[stat])
    hr = LY.hull_r(x[stat])
    below = rho[stat] < hr - 1e-6
    allowed = np.where(np.isnan(fl), hr, fl)
    clear2d = rho[stat] - allowed
    xt = np.asarray(tb["x_mm"])
    fr_tab = np.interp(x[stat], xt, np.asarray(tb["floor_r_mm"]), left=np.nan, right=np.nan)
    clear_tab = rho[stat] - np.where(np.isnan(fr_tab), hr, fr_tab)
    az_lo = np.interp(x[stat], xt, np.asarray(tb["az_lo_deg"]), left=np.nan, right=np.nan)
    az_hi = np.interp(x[stat], xt, np.asarray(tb["az_hi_deg"]), left=np.nan, right=np.nan)
    inb = below & ~np.isnan(az_lo)
    az_margin = np.minimum(az[stat][inb] - az_lo[inb], az_hi[inb] - az[stat][inb])
    h = LY.height_above_hull(W)
    k2 = int(np.argmin(clear2d))
    res["stowed_vs_pockets"] = {
        "samples": int(stat.sum()),
        "min_clear_to_pocket_floor_2d_raster_mm": float(clear2d.min()),
        "at_x_rho_az": [float(x[stat][k2]), float(rho[stat][k2]), float(az[stat][k2])],
        "min_clear_to_floor_r_table_mm": float(np.nanmin(clear_tab)),
        "POCKET_CLEAR_design_mm": P.POCKET_CLEAR,
        "pass_never_below_floor": bool(clear2d.min() >= -0.02 and np.nanmin(clear_tab) >= -0.02),
        "min_az_margin_inside_pocket_deg": float(az_margin.min()) if az_margin.size else None,
        "max_height_over_hull_mm": float(h.max()),
        "max_height_over_hull_aero_x_lt_100_mm": float(h[x < 100].max()),
        "note": ("2D raster = layout.pocket_map (x, azimuth) floor incl. clearances, the floor the static "
                 "tail is cut to; table = pocket_depth.floor_r (deepest point per x). The stowed outer face "
                 "is flush by design (layout residual +0.065).")}
    # ---- stowed detent interface: layout boss (radial D4 at r 19.02..21.196) and static magnet
    det = json.loads(LY.DATA_JSON.read_text())["detents"]["stowed"]
    az_b = math.radians(lay.psi(0) + float(P.DETENTS["stowed_az_deg"]))
    e_rb = np.array([0.0, -math.sin(az_b), math.cos(az_b)])
    xb = float(P.DETENTS["stowed_x"])
    r_top = float(P.DETENTS["stowed_boss_top_r"])
    r_floor = r_top - float(det["boss_height_over_pocket_floor_mm"])
    boss = _cyl(np.array([xb, 0, 0]) + r_floor * e_rb, np.array([xb, 0, 0]) + r_top * e_rb, MAG_D / 2)
    st0 = blade_world(0, P.FOLD_DEG)
    T0 = blade_transform(0, P.FOLD_DEG)
    ss = info["stowed_seat"]
    mag_w = LY.apply(T0, np.asarray(ss["seated_magnet_centre_blade_frame_mm"]))
    ax_w = T0[:3, :3] @ np.asarray(ss["axis_out_blade_frame"])
    stat_w = np.array([xb, 0, 0]) + float(P.DETENTS["stowed_static_magnet_r"]) * e_rb
    gap_pair = (MAG_SEAT_DEPTH - MAG_H) + float(json.loads(LY.DATA_JSON.read_text())["chosen"]["magnet_recess_mm"])
    res["stowed_detent_interface"] = {
        "proposed_static_magnet_centre_world_psi0_mm": (mag_w + (MAG_H + gap_pair) * ax_w).tolist(),
        "proposed_static_magnet_axis_world_psi0": ax_w.tolist(),
        "proposed_pair_gap_mm": gap_pair,
        "layout_boss_vs_blade_overlap_mm3": C.interference(st0, boss),
        "layout_boss_to_blade_min_distance_mm": C.min_distance(st0, boss),
        "blade_magnet_centre_world_psi0_mm": mag_w.tolist(), "blade_magnet_axis_world_psi0": ax_w.tolist(),
        "static_magnet_centre_world_psi0_mm (DETENTS)": stat_w.tolist(),
        "blade_to_static_magnet_centre_distance_mm": float(np.linalg.norm(mag_w - stat_w)),
        "angle_blade_axis_vs_radial_deg": float(math.degrees(math.acos(min(1.0, abs(float(ax_w @ e_rb)))))),
        "note": "DETENTS puts the stowed pair on a radial axis with a flat boss top at r 21.196, but the "
                "blade's inner face there is inclined ~51 deg to radial; the blade seat follows the face "
                "normal. The static side must put its magnet on this axis, behind a boss top that follows "
                "the blade face (see open issues)."}
    try:   # pairing with the static seat when hydroegg.cad.shells is present (other builder)
        from . import shells as SH
        g = SH.stow_detent_geometry()
        sc_, sn_ = np.asarray(g["magnet_centre"], float), np.asarray(g["n"], float)
        off_ = sc_ - mag_w
        res["stowed_detent_interface"]["pairing_with_shells_static_seat"] = {
            "static_magnet_centre_world_psi0_mm": sc_.tolist(), "static_axis_world_psi0": sn_.tolist(),
            "axis_angle_deg": float(math.degrees(math.acos(min(1.0, abs(float(ax_w @ -sn_)))))),
            "lateral_offset_mm": float(np.linalg.norm(off_ - np.dot(off_, ax_w) * ax_w)),
            "face_centre_gap_mm": float(np.linalg.norm((mag_w + MAG_H / 2 * ax_w) - (sc_ + MAG_H / 2 * sn_)))}
    except Exception as e:  # noqa: BLE001
        res["stowed_detent_interface"]["pairing_with_shells_static_seat"] = {"error": f"{type(e).__name__}: {e}"}
    # ---- deployed stop face position (heel -n face) vs ROTOR_HUB stop_face_x
    stop = [f for f in sol.faces() if f.geom_type == GeomType.PLANE and abs(f.center().Z + P.HEEL_T / 2) < 1e-6
            and f.center().X < 0]
    sc = np.array([stop[0].center().X, stop[0].center().Y, stop[0].center().Z])
    bbs = stop[0].bounding_box()
    res["deployed_stop"] = {
        "stop_face_x_world_mm": float(LY.apply(blade_transform(0, 0.0), sc)[0]),
        "ROTOR_HUB_stop_face_x": P.ROTOR_HUB["stop_face_x"],
        "flat_face_from_pin_mm": [float(-bbs.max.X), float(-bbs.min.X)],
        "flat_face_width_mm": float(bbs.size.Y), "stop_face_area_mm2": float(stop[0].area)}
    # ---- root inside the hub slot envelope (swept fold 0..FOLD_DEG)
    slot = json.loads(LY.DATA_JSON.read_text())["hub_slot"]
    root_mesh = _mesh(sol, 0.05)
    R = np.asarray(root_mesh.vertices)
    R = R[(R[:, 0] < 20)]
    xs, fl_s, hw_s = (np.asarray(slot[k]) for k in ("x_mm", "slot_floor_re_mm", "slot_half_width_mm"))
    worst_floor, worst_w = np.inf, -np.inf
    for f in np.arange(0.0, P.FOLD_DEG + 1e-9, 2.0):
        Wr = LY.apply(blade_transform(0, f), R)
        xh, re, tau = LY.hinge_coords(Wr, lay)
        rh = np.hypot(Wr[:, 1], Wr[:, 2])
        m = (xh >= lay.static_tail_x + lay.hub_gap) & (rh < LY.hull_r(xh))
        if not m.any():
            continue
        i = np.clip(np.round((xh[m] - xs[0]) / (xs[1] - xs[0])).astype(int), 0, xs.size - 1)
        worst_floor = min(worst_floor, float((re[m] - (fl_s[i] + lay.run_clear)).min()))
        worst_w = max(worst_w, float((np.abs(tau[m]) - (hw_s[i] - lay.run_clear)).max()))
    res["root_vs_hub_slot"] = {"min_re_minus_swept_floor_mm": worst_floor,
                               "max_tau_minus_swept_halfwidth_mm": worst_w,
                               "note": "swept CAD root vs the layout's swept envelope (slot table minus its 0.3 "
                                       "run clearance, 0.5 mm x bins): >= -0.3 keeps the running clearance >= 0; "
                                       "the tessellation adds ~0.05"}
    # ---- pivot proof
    pv = pivot(0)
    st_ = blade_world(0, P.FOLD_DEG)
    rot = st_.rotate(Axis(pv["origin"], pv["axis"]), pv["deployed_deg"])
    dep = blade_world(0, 0.0)
    b1, b2 = rot.bounding_box(), dep.bounding_box()
    bb_err = max(abs(a - b) for a, b in zip((b1.min.X, b1.min.Y, b1.min.Z, b1.max.X, b1.max.Y, b1.max.Z),
                                            (b2.min.X, b2.min.Y, b2.min.Z, b2.max.X, b2.max.Y, b2.max.Z)))
    samp = np.asarray(_mesh(rot, 0.05).vertices)
    samp = samp[rng.choice(len(samp), 250, replace=False)]
    dev = _dist_to_faces(dep, samp)
    res["pivot_proof"] = {"pivot_blade_0": pv, "bbox_max_diff_mm": bb_err, "surface_samples": int(len(samp)),
                          "surface_max_dev_mm": float(np.nanmax(dev)),
                          "pass": bool(bb_err < 0.05 and np.nanmax(dev) < 0.05)}
    for i in (1, 2):
        pvi = pivot(i)
        roti = blade_world(i, P.FOLD_DEG).rotate(Axis(pvi["origin"], pvi["axis"]), pvi["deployed_deg"])
        bi, bd = roti.bounding_box(), blade_world(i, 0.0).bounding_box()
        res["pivot_proof"][f"bbox_max_diff_blade_{i}_mm"] = max(
            abs(a - b) for a, b in zip((bi.min.X, bi.min.Y, bi.min.Z, bi.max.X, bi.max.Y, bi.max.Z),
                                       (bd.min.X, bd.min.Y, bd.min.Z, bd.max.X, bd.max.Y, bd.max.Z)))
    # ---- deployed tip radius and world poses
    dep_bb = []
    for i in range(3):
        mv = np.asarray(_mesh(blade_world(i, 0.0), 0.05).vertices)
        dep_bb.append(float(np.hypot(mv[:, 1], mv[:, 2]).max()))
    res["deployed"] = {"tip_radius_max_mm": dep_bb, "ROTOR_R_TIP_TARGET": P.ROTOR_R_TIP_TARGET}
    # ---- print orientation
    pr = print_solid()
    pm_ = _mesh(pr, STL_TOL)
    nz = pm_.face_normals[:, 2]
    zc = pm_.triangles_center[:, 2]
    over = (nz < -math.cos(math.radians(P.MAX_OVERHANG))) & (zc > 0.3)
    res["print"] = {"orientation": "blade frame (s, c, n) -> print (x, y, z) = (s, -n, c): chord vertical, "
                                   "TE side down, knuckle end face on the bed, pin bore vertical",
                    "bbox_mm": C.check(pr)["size_mm"], "fits_bed": bool(C.fits_bed(pr)),
                    "te_height_above_bed_mm": [_print_te_range()["te_lo"], _print_te_range()["te_hi"]],
                    "overhang_area_needing_support_mm2": float(pm_.area_faces[over].sum()),
                    "bed_contact_area_mm2": float(pm_.area_faces[(nz < -0.999) & (zc < 0.05)].sum()),
                    "bore_axis_print": "vertical (z)", "notes": PRINT_NOTE.format(**_print_te_range())}
    res["t_checks_s"] = round(time.time() - t0, 1)
    if write:
        CHECK_OUT.parent.mkdir(parents=True, exist_ok=True)
        CHECK_OUT.write_text(json.dumps(LY._clean(res), indent=1) + "\n")
    return res


# ==========================================================================
# Export and preview
# ==========================================================================

def close_mesh_slivers(m, max_perimeter: float = 0.2) -> int:
    """Close tiny boundary loops (OCC mesher artifacts a few microns wide where faces meet)
    with fans oriented against their neighbours. Returns the number of loops closed;
    larger holes are left alone so a real defect still fails the watertight check."""
    faces = np.asarray(m.faces)
    directed = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    have = {tuple(e) for e in directed}
    nxt = {}
    for a, b in directed:
        if (b, a) not in have:
            nxt[b] = a                       # the hole boundary runs opposite to the face edge
    new, done, closed = [], set(), 0
    for start in list(nxt):
        if start in done:
            continue
        loop, v = [start], nxt[start]
        while v != start and v in nxt and len(loop) < 50:
            loop.append(v)
            v = nxt[v]
        done.update(loop)
        if v != start or len(loop) < 3:
            continue
        P_ = m.vertices[loop]
        per = float(np.linalg.norm(P_ - np.roll(P_, -1, axis=0), axis=1).sum())
        if per <= max_perimeter:
            new += [[loop[0], loop[k], loop[k + 1]] for k in range(1, len(loop) - 1)]
            closed += 1
    if new:
        m.faces = np.vstack([faces, np.asarray(new)])
    return closed


def export() -> dict:
    import trimesh
    from OCP.BRepTools import BRepTools
    pr = print_solid()
    BRepTools.Clean_s(pr.wrapped)          # drop cached finer triangulations from the checks
    STEP_OUT.parent.mkdir(parents=True, exist_ok=True)
    STL_OUT.parent.mkdir(parents=True, exist_ok=True)
    export_step(pr, STEP_OUT)
    export_stl(pr, STL_OUT, tolerance=STL_TOL, angular_tolerance=0.3)
    m = trimesh.load(STL_OUT)
    raw_wt, cleaned, filled = bool(m.is_watertight), 0, 0
    if not raw_wt:
        # OCC's mesher leaves zero-length segments where three fillets meet at the heel
        # tip corners; drop the degenerate slivers and re-weld
        n0 = len(m.faces)
        m.update_faces(m.nondegenerate_faces())
        m.remove_unreferenced_vertices()
        m.merge_vertices()
        cleaned = n0 - len(m.faces)
        filled = close_mesh_slivers(m)
        m.export(STL_OUT)
        m = trimesh.load(STL_OUT)
    return {"step": str(STEP_OUT.relative_to(ROOT)), "stl": str(STL_OUT.relative_to(ROOT)),
            "stl_tolerance_mm": STL_TOL, "stl_raw_watertight": raw_wt, "stl_degenerate_removed": int(cleaned),
            "stl_sliver_holes_closed": int(filled), "stl_winding_consistent": bool(m.is_winding_consistent),
            "stl_watertight": bool(m.is_watertight), "stl_volume_mm3": float(m.volume),
            "stl_volume_vs_solid_pct": float(100 * (m.volume / volume(pr) - 1)),
            "stl_triangles": int(len(m.faces)), "stl_bytes": STL_OUT.stat().st_size,
            "step_bytes": STEP_OUT.stat().st_size}


def render_final(out: Path = WORK / "final.png") -> Path:
    """Iso of the deployed rotor + a stowed view, side by side (scripts/preview.py renders)."""
    WORK.mkdir(parents=True, exist_ok=True)
    jobs = [(WORK / "final_iso.png", "view=55,18&fold=1"), (WORK / "final_stowed.png", "view=140,24&fold=0")]
    for path, q in jobs:
        subprocess.run([sys.executable, str(ROOT / "scripts/preview.py"), "hydroegg.cad.blade:preview_parts",
                        "--out", str(path), "--q", q, "-W", "1000", "-H", "800"], check=True, cwd=ROOT)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    ims = [mpimg.imread(p) for p, _ in jobs]
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.9), dpi=100)
    for a, im, t in zip(axs, ims, ("deployed (fold 0), seen from downstream",
                                   "stowed (fold 90), aft quarter (no hull; pins in steel)")):
        a.imshow(im)
        a.set_title(t, fontsize=12, loc="left", pad=6)
        a.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01, wspace=0.03)
    fig.savefig(out)
    plt.close(fig)
    for p, _ in jobs:
        p.unlink(missing_ok=True)
        for ext in (".glb", ".json"):                     # throwaway preview scenes (~20 MB)
            (ROOT / "out/work/preview" / (p.stem + ext)).unlink(missing_ok=True)
    return out


def main() -> None:
    t0 = time.time()
    sol = blade_frame_solid()
    print("built", C.check(sol, "blade"), f"{time.time() - t0:.0f} s", flush=True)
    ex = export()
    print("exported", ex, flush=True)
    res = checks(write=False)
    res["export"] = ex
    CHECK_OUT.parent.mkdir(parents=True, exist_ok=True)
    CHECK_OUT.write_text(json.dumps(LY._clean(res), indent=1) + "\n")
    print(f"wrote {CHECK_OUT.relative_to(ROOT)} ({time.time() - t0:.0f} s)", flush=True)
    if "--no-render" not in sys.argv:
        print("render", render_final(), flush=True)


if __name__ == "__main__":
    main()
