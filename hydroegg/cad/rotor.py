"""Rotor: hub, spinner, hinge pins, hub detent magnets, hub grub screw and spinner screw.

    scripts/cadrun.sh .venv/bin/python -m hydroegg.cad.rotor      # build, export, check, preview

Sources: params.py (ROTOR_HUB, SPINNER_ACCESS, DETENTS, SHAFT, HUB_GAP, clearances),
data/layout_rotor.json (solved blade axis, hub_slot, detents) and hydroegg/layout.py (blade
surface model and fold transforms). Frame as params.py: X = body axis (flow +X), Z up. Blade
azimuth psi as layout.py: e_r = (0, -sin psi, cos psi), e_t = (0, -cos psi, -sin psi), blade i at
psi = 120 i; a +Z feature is carried to psi by ``shape.rotate(Axis.X, psi)``. Hinge-plane
coordinates of blade 0: x, r_e (along e_r = +Z) and tau (along the pin, e_t = -Y, toward the LE).

Hub (printed, rotating, x 108.6 to 128)
    Body of revolution: the hull inset by spinner_wall + spinner_gap (normal to the surface) from
    the front face to a 45 deg aft cone that ends on a flat aft face (r <= HUB_AFT_FLAT_R) at the
    shaft end. The cone is what lets the spinner print base down with a short bridge instead of a
    full-diameter ceiling; HUB_AFT_FLAT_R (11) keeps the cone clear of the ears (it would touch
    the 2 mm ear cover around the pin holes at 10.3).
    Bore D8.15 (close running) with a D-flat (x 109.12 to 128, flat 0.075 off the shaft flat)
    facing the M3 grub at psi 60; the shaft's flat must run out through its aft end to x <= 109.12
    so the hub slides on (SHAFT.hub_flat_x gives only the grub station).
    Blade slots: the swept volume of the blade root (knuckle, heel, shank and first aero rows of
    the layout's blade surface, fold 0 to FOLD_DEG) cut as four tau bands, each the epigraph of the
    swept floor r_e(x) dilated by CLEAR_ROTATE:
        b0 |tau| <= HEEL_W/2 + run     knuckle + heel + all shank points (the heel's own slot)
        b1 |tau| <= KNUCKLE_L/2 + run  knuckle + shank points with |tau| > HEEL_W/2 (ear inner faces)
        b2  tau in [KNUCKLE_L/2, PIN_L/2]          shank points with tau > KNUCKLE_L/2 (over the ear)
        b3  tau in [PIN_L/2 - run, max tau + run]  shank points with tau > PIN_L/2 - run
    The swept root leans toward the LE (+tau, up to 14.7) and never passes tau -7, so the slot is
    asymmetric (the layout's hub_slot table is the symmetric |tau| envelope binned in 0.5 mm, too
    coarse for the stop: it would put the stop face at x 113.1). The heel band ends forward in an
    exact, zero-clearance stop face at x = HINGE_X - HEEL_T/2 = 113.5 (the deployed stop). Floors
    rise at most 45 deg going aft (the hub prints front face down, the spinner base down, both +X
    up); steeper rises are cut as 45 deg ramps, a true vertical jump of the swept floor is left as
    a flat bridged step.
    Ears: what is left of the hub at |tau| 7.3 to 11 around the pin; ear holes D2.90 (PRESS_FIT)
    for the ISO 8734 3 x 22 dowel, insertion channels D3.5 along the pin from tau 11 to the hub
    surface on both sides. Pin retention: press fit in both ears, and a loose ISO 8734 3 x 4 filler
    dowel in each channel, trapped by the spinner shell, so a loosened pin can move 0.4 mm at most
    and stays in both ears (without the fillers it could travel 4.4 mm, more than the 3.7 mm ear).
    Deployed-detent magnets: D3.90 x 2.2 pockets in the stop faces (DETENTS).
    M3 x 4 cup-point grub, radial at x 111.6, psi 60: tapped M3 in the PLA (2.5 drill, r 4.08 to
    9.5) with a D3.4 clearance bore out to the surface. A heat-set insert does not fit: its D4.0
    hole needs 1.6 mm of wall and leaves 1.0 to the front face and about 1.3 to the swept roots,
    and it would sit 15 mm deep. PLA tapping is enough here: the D-flat carries the torque, the
    rotor thrust goes through the spinner screw, the grub only sets the hub's axial position.
    Front face: flat slinger face 0.6 mm (HUB_GAP) from the static tail end face, with two 90 deg
    V-grooves (labyrinth chambers) around the lip-seal counterbore (r 6.3 and 8.7, 0.5 deep, kept
    inside the grub's tapped zone so 1.25 mm of wall stays to the grub bore); nothing stands proud,
    so the hub prints on it. Three D3 x 1.5 pegs on the aft face key the spinner.

Spinner (printed, rotating, x 108.6 to the tip)
    1.2 mm shell (normal) over the hub with a 0.2 mm gap; per blade the same tau bands, but cut
    through the full depth as windows (no feathered skin where a band floor would graze the inner
    cone). Inner aft wall (3 mm): a flat bridge D22.6 around the screw hole, then a 45 deg cone that
    matches the hub's aft cone; clamped against the hub aft face and the shaft end by one ISO 4762
    M3 x 8 into the tapped shaft end (rotor thrust: hub -> spinner wall -> screw -> shaft). Aft of
    the wall the spinner is solid (perimeters + infill) around the D6 access bore, which opens at
    the tip. The access bore leaves no room for a pointed tip: the shell ends at x 151.7, where the
    R = SPINNER_TIP_R round leaves MIN_WALL + chamfer around the bore. Three holes in the aft wall
    take the hub's pegs (keying). Slide the spinner on with the blades folded >= 15 deg (stowed is
    simplest): fully deployed blades block it 10.5 mm before home.

Printing: hub on its front face, spinner base (front edge) down; see SUPPORT_NOTES.
"""

from __future__ import annotations

import functools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from build123d import Axis, Edge, Face, Location, Shape, Solid, Vector, Wire, export_step

from .. import layout as L
from .. import params as P
from .. import profile
from .common import axis_cyl, check, cyl, fits_bed, on_bed

REPO = Path(__file__).resolve().parents[2]
OUT_STEP = REPO / "out/cad/step"
OUT_STL = REPO / "out/cad/stl"
OUT_CHK = REPO / "out/cad/checks"
WORK = REPO / "out/work/cad/rotor"

# --------------------------------------------------------------------------
# Derived dimensions (all from params; nothing here is a free design value unless marked)
# --------------------------------------------------------------------------
H = P.ROTOR_HUB
X_FRONT = float(H["x_front"])                  # 108.6 = STATIC_TAIL_X + HUB_GAP
X_AFT = float(H["x_aft_face"])                 # 128 = shaft end
BORE_R = H["bore_d"] / 2                       # 4.075
FLAT_R = BORE_R - H["d_flat_depth"]            # 3.575 (shaft flat at 3.5: 0.075 clearance)
FLAT_X = (109.12, X_AFT)                       # DESIGN: D-flat from just aft of the bore chamfer (1.23 mm
                                               # ahead of the grub bore) out through the aft face
GRUB_X, GRUB_PSI = float(H["grub_x"]), float(H["grub_psi_deg"])
STOP_X = float(H["stop_face_x"])               # 113.5
RUN = P.CLEAR_ROTATE                           # 0.3 blade <-> hub/spinner running clearance
SAMPLING = 0.05                                # extra clearance on cloud-sampled floors (point spacing)
SP_WALL, SP_GAP = float(H["spinner_wall"]), float(H["spinner_gap"])
HUB_INSET = SP_WALL + SP_GAP                   # hub surface = hull inset 1.4 (normal)
EAR_HOLE_D = P.PIN_D + 2 * P.PRESS_FIT         # 2.90 printed, pin pressed
CHANNEL_D = P.PIN_D + 2 * P.CLEAR_SLIDE        # 3.5 insertion channel
PIN_HALF = P.PIN_L / 2                         # 11 = ear outer face
EAR_IN = P.KNUCKLE_L / 2 + RUN                 # 7.3 = ear inner face
HEEL_HW = P.HEEL_W / 2 + RUN                   # 5.3 heel band
HUB_AFT_FLAT_R = 11.0                          # DESIGN: flat aft face radius, then a 45 deg cone
MAG_D_HOLE = P.MAGNET["D"] + 2 * P.PRESS_FIT   # 3.90
MAG_H = P.MAGNET["H"]
MAG_RECESS = 0.2                               # layout magnet_recess (DETENTS hub magnet x 112.3)
GRUB_TAP_D, GRUB_CLEAR_D, GRUB_TAP_R1 = 2.5, 3.4, 9.5
GRUB_L = 4.0                                   # M3 x 4 (ROTOR_HUB.grub)
SCREW_L, SCREW_HEAD_D, SCREW_HEAD_H = 8.0, 5.5, 3.0   # ISO 4762 M3 x 8
SCREW_CLEAR_D = 3.4
ACCESS_R = P.SPINNER_ACCESS["d"] / 2           # 3.0
AFT_WALL = float(P.SPINNER_ACCESS["aft_wall"]) # 3.0
X_HEAD = X_AFT + AFT_WALL                      # 131: screw head seat
TIP_CHAMFER = 0.3                              # access bore edge at the tip
PEG_D, PEG_H, PEG_R = 3.0, 1.5, 8.0            # DESIGN: anti-rotation pegs on the hub aft face
PEG_PSI = (60.0, 180.0, 300.0)                 # between the blades
PEG_HOLE_D, PEG_HOLE_DEPTH = PEG_D + 2 * P.CLEAR_SLIDE, PEG_H + 0.3
FILLER_L = 4.0                                 # DESIGN: ISO 8734 3 x 4 filler dowel in each pin channel
GROOVES_R, GROOVE_DEPTH = (6.3, 8.7), 0.5      # DESIGN: front-face labyrinth V-grooves (90 deg)
FRONT_CHAMFER = 0.4                            # elephant-foot relief on the bed edges
BORE_CHAMFER = 0.5
DX = 0.02                                      # floor sampling step along x
R_TOP = 40.0                                   # slot cut tops (outside everything)

MATERIAL_PRINTED = "PETG (v2, wet) / PLA (display)"

SUPPORT_NOTES = {
    "hub": ("front face on the bed (+X up). No supports: outer flank leans <= 24 deg; slot floors "
            "rise <= 45 deg; slot ends are flat bridges across the slot (<= 14.6 mm); the ear holes, "
            "pin channels and grub bore are horizontal holes D <= 3.5 (arched tops); magnet pockets "
            "and the bore open upward; the aft cone is a 45 deg roof; the V-grooves have 45 deg flanks."),
    "spinner": ("front edge (base) on the bed, tip up. No supports: shell leans <= 26 deg; the inner "
                "aft wall is a 45 deg cone plus a flat bridge of D 22.6 around the screw hole; the slot "
                "windows have vertical side walls and flat bridged aft ends (<= 14.6 mm across); the "
                "access bore is vertical; brim recommended (bed contact is the 1.2 mm front edge)."),
}


# --------------------------------------------------------------------------
# Layout model (blade surface reconstructed from the published axis table)
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def layout_model():
    """(data, lay, blade, axis, surf): layout_rotor.json, the chosen Layout, the BEM blade, the
    solved curved pitch axis from the JSON table (no re-solve) and the blade-frame surface."""
    d = json.loads(L.DATA_JSON.read_text())
    lay = L.chosen_layout()
    ch = d["chosen"]
    for k_lay, k_json, k_par in (("hinge_x", "hinge_x_mm", P.HINGE_X), ("hinge_r", "hinge_r_mm", P.HINGE_R),
                                 ("static_tail_x", "static_tail_x_mm", P.STATIC_TAIL_X),
                                 ("fold_deg", "fold_deg", P.FOLD_DEG), ("pin_l", "pin_l_mm", P.PIN_L),
                                 ("knuckle_l", "knuckle_l_mm", P.KNUCKLE_L), ("heel_t", "heel_t_mm", P.HEEL_T),
                                 ("heel_reach", "heel_reach_mm", P.HEEL_REACH), ("heel_w", "heel_w_mm", P.HEEL_W),
                                 ("run_clear", "run_clear_mm", P.CLEAR_ROTATE)):
        v = getattr(lay, k_lay)
        assert abs(v - ch[k_json]) < 1e-9 and abs(v - k_par) < 1e-9, (k_lay, v, ch[k_json], k_par)
    assert abs(lay.knuckle_r * 2 - P.KNUCKLE_D) < 1e-9
    blade = L.load_blade()
    a = d["axis"]
    ax = L.BladeAxis(np.array(a["r_mm"]), np.array(a["s_mm"]), np.array(a["lead_mm"]), np.array(a["aft_mm"]),
                     np.array(a["slope_daft_ds"]), np.array(a["aft_poly_coef_mm"]),
                     np.array(a["lead_poly_coef_mm"]), blade.r_root, blade.r_tip)
    surf = L.blade_surface(blade, lay, ax)
    return d, lay, blade, ax, surf


def hull_r(x):
    return np.asarray(profile.radius(np.asarray(x, float)), float)


def hull_in(x, wall):
    """Hull inset by ``wall`` along the normal (profile.inner_radius is 0.5 um short of the exact
    normal offset on this tail, so 1 um is added to keep every shell >= its nominal wall)."""
    return np.asarray(profile.inner_radius(np.asarray(x, float), wall=wall + 0.001), float)


def shell_thickness() -> dict:
    """Exact minimum normal distance between the spinner's outer and inner profile curves over the
    shell (front face to the start of the inner cone), and between hub and spinner (gap)."""
    xx = np.linspace(100.0, 140.0, 400001)
    ro = hull_r(xx)
    xs = np.linspace(X_FRONT, x_cone_spinner(), 60)
    ri = hull_in(xs, SP_WALL)
    rh = hull_in(np.linspace(X_FRONT, x_cone_hub(), 60), HUB_INSET)
    wall = min(float(np.min(np.hypot(xx - x, ro - r))) for x, r in zip(xs, ri))
    xi = np.linspace(X_FRONT - 2, X_AFT, 200001)
    rin = hull_in(xi, SP_WALL)
    gap = min(float(np.min(np.hypot(xi - x, rin - r))) for x, r in zip(np.linspace(X_FRONT, x_cone_hub(), 60), rh))
    return {"spinner_shell_min_mm": wall, "hub_spinner_gap_min_mm": gap}


# --------------------------------------------------------------------------
# Slot bands (numpy): swept floors in the hinge plane of blade 0
# --------------------------------------------------------------------------

def _dilated_floor(xs: np.ndarray, x: np.ndarray, re: np.ndarray, c: float) -> np.ndarray:
    """Floor r_e(x) of the epigraph of a point set, dilated by a disc of radius c (+ half a bin)."""
    f = np.full(xs.size, np.inf)
    i = np.floor((x - xs[0]) / DX + 0.5).astype(int)
    ok = (i >= 0) & (i < xs.size)
    np.minimum.at(f, i[ok], re[ok])
    ce = c + DX / 2
    k = int(math.ceil(ce / DX))
    out = np.full(xs.size, np.inf)
    for j in range(-k, k + 1):
        dd = j * DX
        if abs(dd) > ce:
            continue
        sh = np.roll(f, j)
        if j > 0:
            sh[:j] = np.inf
        elif j < 0:
            sh[j:] = np.inf
        out = np.minimum(out, sh - math.sqrt(ce * ce - dd * dd))
    return out


def _printable(f: np.ndarray, jump: float = 0.3) -> np.ndarray:
    """Scan aft (+x = up when printed): the floor may rise at most 45 deg (DX per sample); a raw
    jump larger than ``jump`` per sample is a vertical step of the swept floor and stays a flat
    bridged step. A floor can only get lower (more clearance)."""
    g = f.copy()
    k = DX * math.tan(math.radians(P.MAX_OVERHANG))
    for i in range(1, g.size):
        if not (np.isfinite(g[i]) and np.isfinite(g[i - 1])):
            continue
        if f[i] - f[i - 1] > jump:
            continue
        g[i] = min(g[i], g[i - 1] + k)
    return g


def _rdp(pts: np.ndarray, tol: float) -> np.ndarray:
    """Ramer-Douglas-Peucker polyline simplification."""
    if len(pts) < 3:
        return pts
    keep = np.zeros(len(pts), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        p, q = pts[a], pts[b]
        d = q - p
        n = math.hypot(*d)
        seg = pts[a + 1:b] - p
        dist = np.abs(seg[:, 0] * d[1] - seg[:, 1] * d[0]) / n if n > 0 else np.hypot(*seg.T)
        j = int(np.argmax(dist))
        if dist[j] > tol:
            keep[a + 1 + j] = True
            stack += [(a, a + 1 + j), (a + 1 + j, b)]
    return pts[keep]


@functools.lru_cache(maxsize=1)
def slot_bands() -> dict:
    """Band floors and outline polygons for blade 0 (hinge-plane x, r_e; tau extents)."""
    _, lay, _, _, surf = layout_model()
    hx, hr = lay.hinge_x, lay.hinge_r
    xs = 107.6 + DX * np.arange(int(round((134.0 - 107.6) / DX)) + 1)
    folds = np.arange(0.0, lay.fold_deg + 1e-9, 0.5)
    # shank + first aero rows: dense surface cloud, swept
    cloud = np.vstack([L.densify(surf.shank, 0.15).reshape(-1, 3),
                       L.densify(surf.aero[:12], 0.2).reshape(-1, 3)])
    X, R, T = [], [], []
    for f in folds:
        W = L.apply(L.blade_to_world(f, 0.0, lay), cloud)
        x, re, tau = L.hinge_coords(W, lay)
        rho = np.hypot(W[:, 1], W[:, 2])
        m = (x >= xs[0]) & (rho < hull_r(x) + 0.6)
        X.append(x[m]); R.append(re[m]); T.append(tau[m])
    X, R, T = np.concatenate(X), np.concatenate(R), np.concatenate(T)
    eps = 1e-6
    # knuckle: exact circle about the pin, dilated
    fk = np.full(xs.size, np.inf)
    rk = lay.knuckle_r + RUN
    m = np.abs(xs - hx) <= rk
    fk[m] = hr - np.sqrt(rk**2 - (xs[m] - hx) ** 2)
    # heel: outline (s in [-reach, 0], n in [-t/2, t/2]) swept in 0.1 deg steps, dilated
    s1 = np.arange(-lay.heel_reach, 1e-9, 0.01)
    n1 = np.arange(-lay.heel_t / 2, lay.heel_t / 2 + 1e-9, 0.01)
    ss = np.concatenate([s1, s1, np.full(n1.size, -lay.heel_reach)])
    nn = np.concatenate([np.full(s1.size, -lay.heel_t / 2), np.full(s1.size, lay.heel_t / 2), n1])
    fr = np.radians(np.arange(0.0, lay.fold_deg + 1e-9, 0.1))[:, None]
    hxs = (hx - ss[None] * np.sin(fr) + nn[None] * np.cos(fr)).ravel()
    hrs = (hr + ss[None] * np.cos(fr) + nn[None] * np.sin(fr)).ravel()
    fh = _dilated_floor(xs, hxs, hrs, RUN)
    fh[xs < STOP_X - 1e-9] = np.inf              # the stop face: no clearance ahead of the heel
    assert hxs.min() >= STOP_X - 1e-9, "heel crosses the stop plane"
    c = RUN + SAMPLING
    t_max, t_min = float(T.max()), float(T.min())
    bands = []
    raw0 = np.minimum.reduce([fk, fh, _dilated_floor(xs, X, R, c)])
    bands.append(dict(name="b0_heel", tau=(-HEEL_HW, HEEL_HW), raw=raw0))
    m1 = np.abs(T) > P.HEEL_W / 2 + eps
    bands.append(dict(name="b1_knuckle", tau=(-EAR_IN, EAR_IN),
                      raw=np.minimum(fk, _dilated_floor(xs, X[m1], R[m1], c))))
    for sgn in (1, -1):
        m2 = sgn * T > P.KNUCKLE_L / 2 + eps
        if m2.any():
            hi = min(PIN_HALF, float(np.abs(T[m2]).max()) + c)
            lo = P.KNUCKLE_L / 2
            bands.append(dict(name=f"b2_ear{'+' if sgn > 0 else '-'}",
                              tau=(lo, hi) if sgn > 0 else (-hi, -lo),
                              raw=_dilated_floor(xs, X[m2], R[m2], c)))
        m3 = sgn * T > PIN_HALF - RUN + eps
        if m3.any():
            hi = float(np.abs(T[m3]).max()) + c
            lo = PIN_HALF - RUN
            bands.append(dict(name=f"b3_outer{'+' if sgn > 0 else '-'}",
                              tau=(lo, hi) if sgn > 0 else (-hi, -lo),
                              raw=_dilated_floor(xs, X[m3], R[m3], c)))
    for b in bands:
        g = _printable(b["raw"])
        b["floor"] = g
        fin = np.isfinite(g)
        idx = np.nonzero(fin)[0]
        runs = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)
        polys = []
        for r_ in runs:
            pts = np.column_stack([xs[r_], g[r_]])
            if b["name"] == "b0_heel":                     # exact vertical stop step at STOP_X
                left = pts[pts[:, 0] < STOP_X - 1e-9]
                right = pts[pts[:, 0] >= STOP_X - 1e-9]
                fk_stop = hr - math.sqrt(rk**2 - (STOP_X - hx) ** 2)
                r_right = float(right[0, 1])
                right = right[right[:, 0] > STOP_X + 1e-9]
                pts = np.vstack([left, [[STOP_X, fk_stop], [STOP_X, r_right]], right])
                b["stop_face_r_e"] = (r_right, fk_stop)
            pts = _rdp(pts, 0.004)
            poly = np.vstack([[pts[0, 0], R_TOP], pts, [pts[-1, 0], R_TOP]])
            keep = np.ones(len(poly), bool)
            keep[1:] = np.hypot(*np.diff(poly, axis=0).T) > 1e-6
            polys.append(poly[keep])
        b["polys"] = polys
        b["x_range"] = (float(xs[idx].min()), float(xs[idx].max()))
        b["floor_min"] = float(np.min(g[fin]))
    return {"xs": xs, "bands": bands, "tau_range": (t_min, t_max), "n_cloud_pts": int(X.size),
            "fk": fk, "fh": fh}


# --------------------------------------------------------------------------
# Solid helpers
# --------------------------------------------------------------------------

def _V(x, r):
    return Vector(float(x), 0.0, float(r))


def _revolve(segments: list) -> Solid:
    """Revolve a closed (x, r) outline made of ('line', [p0, p1, ...]), ('spline', [...]) and
    ('arc', [p0, pm, p1]) segments, each starting where the previous one ended."""
    edges = []
    for kind, pts in segments:
        vs = [_V(*p) for p in pts]
        if kind == "line":
            edges += [Edge.make_line(a, b) for a, b in zip(vs[:-1], vs[1:]) if (a - b).length > 1e-7]
        elif kind == "spline":
            edges.append(Edge.make_spline(vs))
        elif kind == "arc":
            edges.append(Edge.make_three_point_arc(*vs))
    face = Face(Wire(edges))
    return Solid.revolve(face, 360, Axis.X)


def _spline_pts(fn, x0, x1, n=40):
    s = (1 - np.cos(np.linspace(0, math.pi, n))) / 2
    xs = x0 + (x1 - x0) * s
    return [(float(x), float(fn(x))) for x in xs]


def _solve(fn, lo, hi, it=80):
    flo = fn(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = fn(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _fuse(shapes):
    out = shapes[0]
    for s in shapes[1:]:
        out = out.fuse(s)
    return out.clean()


def _single(shape) -> Solid:
    sols = shape.solids()
    if len(sols) == 1:
        return sols[0]
    return max(sols, key=lambda s: s.volume)


def _isect(a: Shape, b: Shape) -> Solid:
    r = a.intersect(b)
    sols = [x for x in (r if isinstance(r, (list, tuple)) else [r]) for x in x.solids()]
    return _fuse(sols) if len(sols) > 1 else sols[0]


def rot_psi(shape: Shape, psi: float) -> Shape:
    """Carry a shape built for blade 0 (psi = 0) to azimuth psi (layout convention)."""
    return shape.rotate(Axis.X, psi) if psi else shape


def psi_list():
    return [360.0 * i / P.N_BLADES for i in range(P.N_BLADES)]


def hinge_frame(psi: float):
    _, lay, _, _, _ = layout_model()
    e_x, e_r, e_t = L.rotor_axes(psi)
    return L.hinge_point(lay, psi), e_r, e_t


# --------------------------------------------------------------------------
# Slot cut (blade 0), hub, spinner
# --------------------------------------------------------------------------

def _band_prism(poly, lo: float, hi: float) -> Solid:
    """Polygon in the blade-0 hinge plane (x, r_e) extruded along tau from lo to hi."""
    vs = [Vector(float(x), -lo, float(r)) for x, r in poly]          # tau = -y at psi 0
    face = Face(Wire([Edge.make_line(a, c) for a, c in zip(vs, vs[1:] + vs[:1])]))
    return Solid.extrude(face, Vector(0, -1, 0) * (hi - lo))


@functools.lru_cache(maxsize=1)
def slot_cut0() -> Solid:
    """Hub slot for blade 0: union of the band prisms (floor polygons extruded along tau)."""
    return _fuse([_band_prism(poly, *b["tau"]) for b in slot_bands()["bands"] for poly in b["polys"]])


@functools.lru_cache(maxsize=1)
def spinner_window0() -> Solid:
    """Spinner slot for blade 0: each band cut through the full depth (r_e >= 0) over its x range.
    The spinner never needs material under a slot, and a full-depth window leaves no feathered
    skin where a band floor crosses the spinner's inner surface at a grazing angle."""
    out = []
    for b in slot_bands()["bands"]:
        for poly in b["polys"]:
            x0, x1 = float(poly[:, 0].min()), float(poly[:, 0].max())
            rect = np.array([[x0, -1.0], [x1, -1.0], [x1, R_TOP], [x0, R_TOP]])
            out.append(_band_prism(rect, *b["tau"]))
    return _fuse(out)


def x_cone_hub() -> float:
    return _solve(lambda x: hull_in(x, HUB_INSET) - (HUB_AFT_FLAT_R + X_AFT - x), 112.0, X_AFT)


@functools.lru_cache(maxsize=1)
def hub_body() -> Solid:
    """Hub body of revolution (no bore, no slots)."""
    rc = lambda x: float(hull_in(x, HUB_INSET))  # noqa: E731
    xc = x_cone_hub()
    x1 = X_FRONT + FRONT_CHAMFER
    segs = [("line", [(X_FRONT, 0.0), (X_FRONT, rc(X_FRONT) - FRONT_CHAMFER), (x1, rc(x1))]),
            ("spline", _spline_pts(rc, x1, xc, 40)),
            ("line", [(xc, rc(xc)), (X_AFT, HUB_AFT_FLAT_R), (X_AFT, 0.0), (X_FRONT, 0.0)])]
    return _revolve(segs)


def _d_bore() -> Solid:
    x0, x1 = FLAT_X
    full = axis_cyl(X_FRONT - 1.0, x0 + 0.01, BORE_R)
    d = axis_cyl(x0, x1 + 1.0, BORE_R)
    cutter = Solid.make_box(x1 - x0 + 4, 20, 10).moved(Location((x0 - 2, -10, FLAT_R)))
    d = d.cut(rot_psi(cutter, GRUB_PSI))
    cham = _revolve([("line", [(X_FRONT - 0.5, BORE_R + BORE_CHAMFER + 0.5), (X_FRONT + BORE_CHAMFER, BORE_R - 0.01),
                               (X_FRONT - 0.5, BORE_R - 0.01), (X_FRONT - 0.5, BORE_R + BORE_CHAMFER + 0.5)])])
    return _fuse([full, d, cham])


def _grooves() -> Solid:
    out = []
    for r in GROOVES_R:
        h = GROOVE_DEPTH
        out.append(_revolve([("line", [(X_FRONT - 0.3, r - h - 0.3), (X_FRONT + h, r), (X_FRONT - 0.3, r + h + 0.3),
                                       (X_FRONT - 0.3, r - h - 0.3)])]))
    return _fuse(out)


def _pin_holes0() -> Solid:
    """Ear hole (press) and both insertion channels for blade 0, along e_t."""
    Hp, e_r, e_t = hinge_frame(0.0)
    Hp = np.asarray(Hp)
    ear = cyl(Hp - (PIN_HALF + 0.05) * e_t, Hp + (PIN_HALF + 0.05) * e_t, EAR_HOLE_D / 2)
    ch1 = cyl(Hp + (PIN_HALF - 0.05) * e_t, Hp + 24.0 * e_t, CHANNEL_D / 2)
    ch2 = cyl(Hp - (PIN_HALF - 0.05) * e_t, Hp - 24.0 * e_t, CHANNEL_D / 2)
    return _fuse([ear, ch1, ch2])


def magnet_center0() -> np.ndarray:
    xm, re, tau = P.DETENTS["deployed_hub_magnet_x_re_tau"]
    _, e_r, e_t = hinge_frame(0.0)
    return np.array([xm, 0, 0]) + re * np.asarray(e_r) + tau * np.asarray(e_t)


def _magnet_pocket0() -> Solid:
    c = magnet_center0()
    x0 = STOP_X - MAG_RECESS - MAG_H
    return cyl((x0, c[1], c[2]), (STOP_X + 0.3, c[1], c[2]), MAG_D_HOLE / 2)


def _grub_hole() -> Solid:
    tap = cyl((GRUB_X, 0, 3.0), (GRUB_X, 0, GRUB_TAP_R1), GRUB_TAP_D / 2)
    clr = cyl((GRUB_X, 0, GRUB_TAP_R1 - 0.01), (GRUB_X, 0, 30.0), GRUB_CLEAR_D / 2)
    return rot_psi(_fuse([tap, clr]), GRUB_PSI)


def _pegs() -> list:
    out = []
    for psi in PEG_PSI:
        p = _revolve([("line", [(X_AFT - 0.2, 0.0), (X_AFT - 0.2, PEG_D / 2), (X_AFT + PEG_H - 0.3, PEG_D / 2),
                                (X_AFT + PEG_H, PEG_D / 2 - 0.3), (X_AFT + PEG_H, 0.0), (X_AFT - 0.2, 0.0)])])
        out.append(rot_psi(p.moved(Location((0, 0, PEG_R))), psi))
    return out


@functools.lru_cache(maxsize=1)
def hub() -> Solid:
    body = hub_body()
    tools = [_d_bore(), _grooves(), _grub_hole()]
    cut0 = _fuse([slot_cut0(), _pin_holes0(), _magnet_pocket0()])
    tools += [rot_psi(cut0, psi) for psi in psi_list()]
    body = body.cut(_fuse(tools))
    body = _fuse([body] + _pegs())
    return _single(body)


def tip_geometry() -> dict:
    """Spinner tip: end plane x_e where the R = SPINNER_TIP_R round, tangent to the hull and to the
    end face, leaves MIN_WALL + chamfer of flat face outside the access bore."""
    R = P.SPINNER_TIP_R
    target = ACCESS_R + TIP_CHAMFER + P.MIN_WALL      # radius where the round meets the end face

    def fillet(xe):
        xc = xe - R
        xx = np.linspace(xc - 2 * R, xe + 0.5, 4000)
        rr = hull_r(xx)

        def gap(rc):
            return float(np.min(np.hypot(xx - xc, rr - rc))) - R
        rc = _solve(lambda r: gap(r), 0.0, float(hull_r(xc)) - 1e-6)
        k = int(np.argmin(np.hypot(xx - xc, rr - rc)))
        return xc, rc, float(xx[k]), float(rr[k])

    xe = _solve(lambda x: fillet(x)[1] - target, 140.0, 160.0, it=50)
    xc, rc, xt, rt = fillet(xe)
    a0 = math.atan2(rt - rc, xt - xc)
    am = 0.5 * a0
    mid = (xc + R * math.cos(am), rc + R * math.sin(am))
    return {"x_end": xe, "r_face_top": rc, "tangent": (xt, rt), "mid": mid, "center": (xc, rc),
            "hull_r_at_end": float(hull_r(xe))}


def x_cone_spinner() -> float:
    off = HUB_AFT_FLAT_R + SP_GAP * math.sqrt(2)
    return _solve(lambda x: hull_in(x, SP_WALL) - (off + X_AFT - x), 112.0, X_AFT)


@functools.lru_cache(maxsize=1)
def spinner_body() -> Solid:
    ro = lambda x: float(hull_r(x))              # noqa: E731
    ri = lambda x: float(hull_in(x, SP_WALL))    # noqa: E731
    tg = tip_geometry()
    xe, rft = tg["x_end"], tg["r_face_top"]
    xt, rt = tg["tangent"]
    xci = x_cone_spinner()
    r_ceil = HUB_AFT_FLAT_R + SP_GAP * math.sqrt(2)
    x1 = X_FRONT + 0.3
    segs = [("line", [(X_FRONT, ri(X_FRONT)), (X_FRONT, ro(X_FRONT) - 0.3), (x1, ro(x1))]),
            ("spline", _spline_pts(ro, x1, xt, 60)),
            ("arc", [(xt, rt), tg["mid"], (xe, rft)]),
            ("line", [(xe, rft), (xe, ACCESS_R + TIP_CHAMFER), (xe - TIP_CHAMFER, ACCESS_R),
                      (X_HEAD, ACCESS_R), (X_HEAD, SCREW_CLEAR_D / 2), (X_AFT, SCREW_CLEAR_D / 2),
                      (X_AFT, r_ceil), (xci, ri(xci))]),
            ("spline", list(reversed(_spline_pts(ri, X_FRONT, xci, 40))))]
    return _revolve(segs)


@functools.lru_cache(maxsize=1)
def spinner() -> Solid:
    body = spinner_body()
    tools = [rot_psi(cyl((X_AFT - 0.5, 0, PEG_R), (X_AFT + PEG_HOLE_DEPTH, 0, PEG_R), PEG_HOLE_D / 2), psi)
             for psi in PEG_PSI]
    tools += [rot_psi(spinner_window0(), psi) for psi in psi_list()]
    return _single(body.cut(_fuse(tools)))


# --------------------------------------------------------------------------
# Bought parts (simple, dimensionally accurate models)
# --------------------------------------------------------------------------

def _round_bar(L_, r, ch0=0.0, ch1=0.0) -> Solid:
    """Cylinder along +X from 0 to L_ with end chamfers (axis profile)."""
    pts = [(0.0, 0.0), (0.0, r - ch0), (ch0, r), (L_ - ch1, r), (L_, r - ch1), (L_, 0.0), (0.0, 0.0)]
    return _revolve([("line", pts)])


def _hex_prism(s: float, depth: float) -> Solid:
    """Hex socket cutter (across flats s) along +X from 0 to depth."""
    R = s / math.sqrt(3)
    vs = [Vector(0, R * math.cos(math.radians(60 * k + 30)), R * math.sin(math.radians(60 * k + 30)))
          for k in range(6)]
    face = Face(Wire([Edge.make_line(a, b) for a, b in zip(vs, vs[1:] + vs[:1])]))
    return Solid.extrude(face, Vector(depth, 0, 0))


def _place_along(shape: Shape, p0, direction) -> Shape:
    """Move a part built along +X from the origin so it starts at p0 and points along direction."""
    d = Vector(*direction).normalized()
    ax = Vector(1, 0, 0)
    ang = math.degrees(math.acos(max(-1.0, min(1.0, ax.dot(d)))))
    s = shape
    if ang > 1e-9:
        axis = ax.cross(d)
        if axis.length < 1e-9:
            axis = Vector(0, 0, 1)
        s = s.rotate(Axis((0, 0, 0), axis.normalized()), ang)
    return s.moved(Location(tuple(p0)))


def pin_model() -> Solid:
    return _round_bar(P.PIN_L, P.PIN_D / 2, 0.3, 0.3)


def pin(i: int) -> Solid:
    Hp, _, e_t = hinge_frame(psi_list()[i])
    Hp = np.asarray(Hp)
    return _place_along(pin_model(), Hp - PIN_HALF * np.asarray(e_t), e_t)


def filler(i: int, side: int) -> Solid:
    """ISO 8734 3 x 4 filler dowel in the insertion channel of pin i, side +1 / -1, against the pin end."""
    Hp, _, e_t = hinge_frame(psi_list()[i])
    Hp, e_t = np.asarray(Hp), np.asarray(e_t)
    p0 = Hp + side * PIN_HALF * e_t if side > 0 else Hp - (PIN_HALF + FILLER_L) * e_t
    return _place_along(_round_bar(FILLER_L, P.PIN_D / 2, 0.3, 0.3), p0, e_t)


def magnet_model() -> Solid:
    return _round_bar(MAG_H, P.MAGNET["D"] / 2, 0.1, 0.1)


def hub_magnet(i: int) -> Solid:
    c = magnet_center0()
    m = magnet_model().moved(Location((c[0] - MAG_H / 2, c[1], c[2])))
    return rot_psi(m, psi_list()[i])


def grub_model() -> Solid:
    """ISO 4029 M3 x 4 cup point, built along +X from the cup (x 0) to the socket end (x 4)."""
    body = _round_bar(GRUB_L, 1.5, 0.35, 0.2)
    cup = Solid.make_cylinder(0.7, 0.4, plane=_plane_x(-0.001))
    sock = _hex_prism(1.5, 1.6).moved(Location((GRUB_L - 1.6 + 0.001, 0, 0)))
    return _single(body.cut(cup).cut(sock))


def _plane_x(x0: float):
    from build123d import Plane
    return Plane(origin=(x0, 0, 0), x_dir=(0, 1, 0), z_dir=(1, 0, 0))


def grub() -> Solid:
    e = np.array([0.0, -math.sin(math.radians(GRUB_PSI)), math.cos(math.radians(GRUB_PSI))])
    r_flat = P.SHAFT_D / 2 - P.SHAFT["hub_flat_depth"]           # 3.5: the shaft's flat
    return _place_along(grub_model(), np.array([GRUB_X, 0, 0]) + r_flat * e, e)


def screw_model() -> Solid:
    """ISO 4762 M3 x 8 along +X: head bearing face at x 0, head to x -3, shank to x +8 (thread)."""
    shank = _round_bar(SCREW_L, 1.5, 0.0, 0.3)
    head = _round_bar(SCREW_HEAD_H, SCREW_HEAD_D / 2, 0.3, 0.0).moved(Location((-SCREW_HEAD_H, 0, 0)))
    sock = _hex_prism(2.5, 1.5).moved(Location((-SCREW_HEAD_H - 0.001, 0, 0)))
    return _single(_fuse([shank, head]).cut(sock))


def spinner_screw() -> Solid:
    """Head on the spinner's counterbore floor (x 131), thread into the shaft end (x 128 -> 123)."""
    return screw_model().rotate(Axis.Y, 180).moved(Location((X_HEAD, 0, 0)))


# --------------------------------------------------------------------------
# Part dicts
# --------------------------------------------------------------------------

def _print_pose(s: Shape) -> Shape:
    return on_bed(s, (0, -90, 0))                # +X up: front face / front edge on the bed


@functools.lru_cache(maxsize=1)
def _parts_cached() -> tuple:
    out = []
    h = hub()
    out.append(dict(name="hub", shape=h, print=_print_pose(h), material="hub", role="body", rotor=True,
                    qty=1, bought=False, explode=[0, 0, 0],
                    notes="Print front face down, " + SUPPORT_NOTES["hub"] + " " + MATERIAL_PRINTED
                          + ", 4 perimeters, 40 % infill. Tap the grub bore M3 after printing; press "
                          "the pins and magnets."))
    sp = spinner()
    out.append(dict(name="spinner", shape=sp, print=_print_pose(sp), material="hub", role="body", rotor=True,
                    qty=1, bought=False, explode=[40, 0, 0],
                    notes="Print base (front edge) down, " + SUPPORT_NOTES["spinner"] + " "
                          + MATERIAL_PRINTED + ", 3 perimeters, 30 % infill in the solid tail."))
    for i, psi in enumerate(psi_list()):
        _, e_r, e_t = hinge_frame(psi)
        out.append(dict(name=f"hinge_pin_{i}", shape=pin(i), print=None, material="steel", role="internal",
                        rotor=True, qty=1, bought=True,
                        explode=[round(float(v), 3) for v in 16.0 * np.asarray(e_t) + 10.0 * np.asarray(e_r)],
                        notes="ISO 8734 parallel pin 3 x 22 m6, stainless (A1); pressed into the ears "
                              "(D2.90 printed), blade knuckle runs free on it. VERIFY stock length."))
    for i, psi in enumerate(psi_list()):
        _, e_r, e_t = hinge_frame(psi)
        for k, side in enumerate((1, -1)):
            out.append(dict(name=f"pin_filler_{i}{'ab'[k]}", shape=filler(i, side), print=None, material="steel",
                            role="internal", rotor=True, qty=1, bought=True,
                            explode=[round(float(v), 3) for v in side * 24.0 * np.asarray(e_t) + 10.0 * np.asarray(e_r)],
                            notes="ISO 8734 parallel pin 3 x 4 m6 (A1) used as a loose filler in the pin "
                                  "insertion channel: pin + filler can move at most 0.4 mm before the "
                                  "spinner stops them, so a loosened hinge pin stays in both ears."))
    for i, psi in enumerate(psi_list()):
        out.append(dict(name=f"hub_magnet_{i}", shape=hub_magnet(i), print=None, material="magnet",
                        role="internal", rotor=True, qty=1, bought=True,
                        explode=[8.0, 0, 0],
                        notes="N52 disc 4 x 2, axial; pressed (or CA-glued) into the stop face pocket, "
                              "0.2 mm below the face; attracting the blade's heel magnet."))
    eg = [0.0, -math.sin(math.radians(GRUB_PSI)), math.cos(math.radians(GRUB_PSI))]
    out.append(dict(name="hub_grub", shape=grub(), print=None, material="steel", role="internal", rotor=True,
                    qty=1, bought=True, explode=[0.0, round(22 * eg[1], 3), round(22 * eg[2], 3)],
                    notes="ISO 4029 M3 x 4 cup point, A2; on the shaft D-flat, 1.5 mm hex key."))
    out.append(dict(name="spinner_screw", shape=spinner_screw(), print=None, material="steel", role="internal",
                    rotor=True, qty=1, bought=True, explode=[70, 0, 0],
                    notes="ISO 4762 M3 x 8 socket head, A2 (A4 for v2), medium threadlocker; clamps spinner "
                          "and hub to the shaft end, carries the rotor thrust."))
    return tuple(out)


def parts() -> list[dict]:
    return [dict(p) for p in _parts_cached()]


def _strip(p: dict) -> dict:
    return {k: v for k, v in p.items() if k != "print"}


def preview_parts() -> list[dict]:
    return [_strip(p) for p in parts()]


# ---- preview helpers (blade-root proxy meshes from the layout surface, sections) ----

def _grid_mesh(G: np.ndarray):
    """Closed-ish quad-strip mesh of a (Ns, Np, 3) loft grid with fan end caps."""
    import trimesh
    ns, npt = G.shape[:2]
    V = G.reshape(-1, 3)
    F = []
    for i in range(ns - 1):
        for j in range(npt):
            a, b = i * npt + j, i * npt + (j + 1) % npt
            c_, d_ = (i + 1) * npt + (j + 1) % npt, (i + 1) * npt + j
            F += [[a, b, c_], [a, c_, d_]]
    cen0, cen1 = len(V), len(V) + 1
    V = np.vstack([V, G[0].mean(0), G[-1].mean(0)])
    for j in range(npt):
        F.append([cen0, (j + 1) % npt, j])
        F.append([cen1, (ns - 1) * npt + j, (ns - 1) * npt + (j + 1) % npt])
    m = trimesh.Trimesh(V, np.array(F), process=True)
    trimesh.repair.fix_normals(m)
    return m


@functools.lru_cache(maxsize=4)
def blade_root_mesh_blade_frame(n_aero: int | None = None):
    """Blade root (knuckle, heel, shank + aero span) as trimeshes in the blade frame (s, c, n).
    Preview proxy built from the layout surface; not the blade part."""
    import trimesh
    _, lay, _, _, surf = layout_model()
    loft = np.concatenate([surf.shank[:-1], surf.aero[:n_aero]], axis=0)
    parts_ = [_grid_mesh(loft)]
    kn = trimesh.creation.cylinder(radius=lay.knuckle_r, height=lay.knuckle_l, sections=96)
    kn.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [1, 0, 0]))   # axis -> y (c)
    parts_.append(kn)
    hb = trimesh.creation.box(extents=[lay.heel_reach, lay.heel_w, lay.heel_t])
    hb.apply_translation([-lay.heel_reach / 2, 0, 0])
    parts_.append(hb)
    return parts_


def blade_proxy(i: int, fold: float, n_aero: int | None = None):
    import trimesh
    _, lay, _, _, _ = layout_model()
    T = L.blade_to_world(fold, psi_list()[i], lay)
    ms = []
    for m in blade_root_mesh_blade_frame(n_aero):
        mm = m.copy()
        mm.apply_transform(T)
        ms.append(mm)
    return trimesh.util.concatenate(ms)


def _half(shape: Shape, keep_neg_y: bool = True) -> Shape:
    box = Solid.make_box(200, 100, 200).moved(Location((40, 0 if keep_neg_y else -100, -100)))
    return shape.cut(box)


def preview_assembly_deployed() -> list[dict]:
    ps = preview_parts()
    for i in range(P.N_BLADES):
        ps.append({"name": f"blade_proxy_{i}", "mesh": blade_proxy(i, 0.0), "material": "blade"})
    return ps


def preview_assembly_stowed() -> list[dict]:
    ps = preview_parts()
    for i in range(P.N_BLADES):
        ps.append({"name": f"blade_proxy_{i}", "mesh": blade_proxy(i, P.FOLD_DEG), "material": "blade"})
    return ps


def preview_hub_only() -> list[dict]:
    keep = {"hub", "hinge_pin_0", "hinge_pin_1", "hinge_pin_2", "hub_magnet_0", "hub_magnet_1",
            "hub_magnet_2", "hub_grub"}
    out = []
    for p in preview_parts():
        if p["name"] in keep:
            q = dict(p)
            q["role"] = "body"
            out.append(q)
    return out


def preview_hub_bare() -> list[dict]:
    return [p for p in preview_parts() if p["name"] == "hub"]


def preview_spinner_only() -> list[dict]:
    return [p for p in preview_parts() if p["name"] == "spinner"]


def preview_section() -> list[dict]:
    """Hub + spinner + pin + screw cut on the blade-0 hinge plane (keep y <= 0, the +tau side)."""
    out = []
    for p in preview_parts():
        if p["name"] in ("hub", "spinner", "spinner_screw", "hinge_pin_0", "hub_magnet_0"):
            q = dict(p)
            q["shape"] = _half(p["shape"], keep_neg_y=False)
            q["role"] = "body"
            q["name"] = p["name"] + "_half"
            out.append(q)
    for f in (0.0, P.FOLD_DEG):
        out.append({"name": f"blade_proxy_f{int(f)}", "mesh": blade_proxy(0, f, 4), "material": "blade"})
    return out


def _light(ps: list[dict]) -> list[dict]:
    """Inspection colours: printed parts in light materials so edges read in the render."""
    out = []
    for p in ps:
        q = dict(p)
        if q.get("material") == "hub":
            q["material"] = "generator" if q["name"].startswith("hub") else "default"
        out.append(q)
    return out


def inspect_hub() -> list[dict]:
    return _light([p for p in preview_parts() if p["name"] == "hub"])


def inspect_spinner() -> list[dict]:
    return _light([p for p in preview_parts() if p["name"] == "spinner"])


def inspect_section() -> list[dict]:
    return _light(preview_section())


def inspect_roots(fold: float = 0.0) -> list[dict]:
    ps = _light([p for p in preview_parts() if p["name"] in ("hub", "spinner")])
    for i in range(P.N_BLADES):
        ps.append({"name": f"blade_proxy_{i}", "mesh": blade_proxy(i, fold, 4), "material": "blade"})
    return ps


def preview_print_poses() -> list[dict]:
    out = []
    for k, p in enumerate([q for q in parts() if q["print"] is not None]):
        s = p["print"].moved(Location((k * 70.0, 0, 0)))
        out.append({"name": p["name"] + "_print", "shape": s, "material": p["material"]})
    return out


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

def to_mesh(shape: Shape, tol: float = 0.02, angular: float = 0.1):
    import trimesh
    v, f = shape.tessellate(tol, angular)
    return trimesh.Trimesh(np.array([(p.X, p.Y, p.Z) for p in v]), np.array(f), process=True)


def to_manifold(mesh):
    import manifold3d as mf
    m = mf.Manifold(mf.Mesh(np.asarray(mesh.vertices, np.float32), np.asarray(mesh.faces, np.uint32)))
    if m.is_empty():
        raise RuntimeError(f"mesh is not a closed manifold: {m.status()}")
    return m


def _mf_transform(m, T: np.ndarray):
    return m.transform(np.asarray(T[:3, :4], float))


@functools.lru_cache(maxsize=1)
def root_manifolds(n_aero: int = 14):
    """Blade root proxy as manifolds in the blade frame (s, c, n): (loft, knuckle, heel). The loft
    (shank + first aero rows of layout.blade_surface) is the union of the convex hulls of
    consecutive sections, so slight section concavities are filled (conservative)."""
    import manifold3d as mf
    _, lay, _, _, surf = layout_model()
    loft = np.concatenate([surf.shank[:-1], surf.aero[:n_aero]], axis=0)
    hulls = [mf.Manifold.hull_points(np.vstack([loft[k], loft[k + 1]])) for k in range(len(loft) - 1)]
    loft_m = mf.Manifold.batch_boolean(hulls, mf.OpType.Add)
    kn = mf.Manifold.cylinder(lay.knuckle_l, lay.knuckle_r, -1.0, 180, True).rotate((90.0, 0.0, 0.0))
    heel = mf.Manifold.cube((lay.heel_reach, lay.heel_w, lay.heel_t), True).translate((-lay.heel_reach / 2, 0, 0))
    return loft_m, kn, heel


def swept_clearance(hub_m, spin_m, d_fold: float = 1.0) -> dict:
    """Min gap of the blade-0 root proxy (knuckle, heel, shank, aero span) to the hub and spinner
    meshes over fold 0..FOLD_DEG. The heel is checked separately: at fold 0 its leading face lands
    on the stop face by design (contact)."""
    import manifold3d as mf
    _, lay, _, _, _ = layout_model()
    loft_m, kn_m, heel_m = root_manifolds()
    root_noheel = loft_m + kn_m
    res = {"folds_deg": [0.0, lay.fold_deg, d_fold]}
    worst = {"hub": (9e9, None), "spinner": (9e9, None)}
    worst_heel = {"hub": (9e9, None), "spinner": (9e9, None)}
    for f in np.arange(0.0, lay.fold_deg + 1e-9, d_fold):
        T = L.blade_to_world(f, 0.0, lay)
        r = _mf_transform(root_noheel, T)
        h = _mf_transform(heel_m, T)
        for name, tgt in (("hub", hub_m), ("spinner", spin_m)):
            g = tgt.min_gap(r, 2.0)
            if g < worst[name][0]:
                worst[name] = (g, float(f))
            if f >= 6.0 - 1e-9:
                gh = tgt.min_gap(h, 2.0)
                if gh < worst_heel[name][0]:
                    worst_heel[name] = (gh, float(f))
    T0 = L.blade_to_world(0.0, 0.0, lay)
    h0 = _mf_transform(heel_m, T0)
    inter = (hub_m ^ h0).volume()
    res.update(root_min_gap_hub_mm=worst["hub"][0], root_min_gap_hub_fold=worst["hub"][1],
               root_min_gap_spinner_mm=worst["spinner"][0], root_min_gap_spinner_fold=worst["spinner"][1],
               heel_min_gap_hub_fold_ge6_mm=worst_heel["hub"][0], heel_min_gap_hub_fold=worst_heel["hub"][1],
               heel_min_gap_spinner_fold_ge6_mm=worst_heel["spinner"][0],
               heel_stop_contact_gap_fold0_mm=hub_m.min_gap(h0, 2.0),
               heel_stop_overlap_fold0_mm3=float(inter),
               note="blade root proxy from layout.blade_surface (shank + aero loft, knuckle cylinder R4 x 14, "
                    "heel box 8 x 10 x 5) vs tessellated hub/spinner (0.02 mm); stop contact at fold 0 is the "
                    "design (heel face on the stop face); heel run clearance reported from fold 6 deg, where "
                    "the opening stop exceeds the 0.3 mm run clearance")
    return res


def spinner_slide_on(spin_m, folds=(0.0, 15.0, 30.0, 45.0, 60.0, 90.0), stroke: float = 50.0) -> dict:
    """Spinner slid on from aft over the hub with the blades fitted: min gap blade roots (all three)
    vs the spinner translated aft by d in [0, stroke]."""
    _, lay, _, _, _ = layout_model()
    loft_m, kn_m, heel_m = root_manifolds()
    root = loft_m + kn_m + heel_m
    out = {}
    for f in folds:
        roots = None
        for i in range(P.N_BLADES):
            r = _mf_transform(root, L.blade_to_world(f, psi_list()[i], lay))
            roots = r if roots is None else roots + r
        worst = (9e9, None)
        for d in np.arange(0.0, stroke + 1e-9, 0.5):
            g = spin_m.translate((float(d), 0.0, 0.0)).min_gap(roots, 2.0)
            if g < worst[0]:
                worst = (g, float(d))
        out[f"fold_{int(f)}"] = {"min_gap_mm": worst[0], "at_offset_mm": worst[1]}
    return out


def _ray_hit(m, o, e):
    """(distance, hit normal) of the first hit along the segment o -> e, or None."""
    hits = m.ray_cast(tuple(map(float, o)), tuple(map(float, e)))
    best = None
    for h in hits:
        d = float(np.linalg.norm(np.asarray(h.position) - o))
        if d > 1e-3 and (best is None or d < best[0]):
            best = (d, np.asarray(h.normal, float))
    return best


def _ray_len(m, o, e) -> float | None:
    h = _ray_hit(m, o, e)
    return None if h is None else h[0]


def wall_thickness(mesh, n_rays: int = 30000, seed: int = 1) -> dict:
    """Local wall thickness by ray casting along the inward normal from random surface points
    (manifold ray_cast). A hit counts as a wall when the opposite face is within 25 deg of
    parallel (|n_hit . n| > 0.9); hits on faces meeting at a corner (45 deg chamfers, 90 deg
    edges) measure the distance to an edge, not a wall, and are counted separately."""
    import trimesh
    m = to_manifold(mesh)
    pts, fi = trimesh.sample.sample_surface(mesh, n_rays, seed=seed)
    nrm = mesh.face_normals[fi]
    best, corner = [], 0
    for p, n in zip(pts, nrm):
        o = p - 2e-3 * n
        h = _ray_hit(m, o, p - 10.0 * n)
        if h is None:
            continue
        if abs(float(np.dot(h[1], n))) < 0.9:
            corner += 1
            continue
        best.append((h[0] + 2e-3, p))
    best.sort(key=lambda t: t[0])
    return {"min_mm": float(best[0][0]), "at": [round(float(v), 2) for v in best[0][1]],
            "p01_mm": float(np.percentile([b[0] for b in best], 1)), "n_wall_rays": len(best),
            "n_corner_rays": corner,
            "thinnest": [[round(b[0], 3)] + [round(float(v), 2) for v in b[1]] for b in best[:8]]}


def overhang_report(mesh_print, tol_deg: float = 1.0) -> dict:
    """Downward-facing area in the print pose steeper than MAX_OVERHANG (+ tol_deg for tessellation
    noise on exact 45 deg faces), ignoring the bed face (z ~ 0), split into horizontal ceilings
    (bridges) and inclined overhangs; inclined faces are clustered (1 mm grid) and listed."""
    fn = mesh_print.face_normals
    a = mesh_print.area_faces
    c = mesh_print.triangles_center
    ang = np.degrees(np.arcsin(np.clip(-fn[:, 2], -1, 1)))          # overhang from vertical
    down = (ang > P.MAX_OVERHANG + tol_deg) & (c[:, 2] > 0.05)
    flat = down & (ang > 89.0)
    incl = down & ~flat
    out = {"bed_contact_mm2": float(a[(ang > 89.0) & (c[:, 2] <= 0.05)].sum()),
           "bridge_area_mm2": float(a[flat].sum()), "inclined_overhang_area_mm2": float(a[incl].sum()),
           "inclined_overhang_max_deg": float(ang[incl].max()) if incl.any() else None}
    if incl.any():
        keys = np.round(c[incl] / 2.0).astype(int)
        uniq, inv = np.unique(keys, axis=0, return_inverse=True)
        areas = np.bincount(inv.ravel(), weights=a[incl])
        order = np.argsort(-areas)[:8]
        out["inclined_clusters_print_xyz_area"] = [[int(v) * 2 for v in uniq[k]] + [round(float(areas[k]), 2)]
                                                   for k in order]
    if flat.any():
        out["bridge_z_mm"] = sorted({round(float(z), 1) for z in c[flat, 2]})[:12]
    return out


def pin_retention(spin_m) -> dict:
    """Axial travel of each hinge pin along its axis (both directions) until its end face meets the
    spinner (ray cast from 37 points of the pin end disc against the spinner mesh), and what is
    then left of the pin in the far ear."""
    _, lay, _, _, _ = layout_model()
    res = []
    for psi in psi_list():
        Hp, e_r, e_t = (np.asarray(v) for v in hinge_frame(psi))
        e_x = np.array([1.0, 0, 0])
        for sgn in (1, -1):
            best = 9e9
            for rr in (0.0, 0.75, 1.5):
                for a in np.radians(np.arange(0, 360, 30 if rr else 360)):
                    o = Hp + sgn * PIN_HALF * e_t + rr * (math.cos(a) * e_x + math.sin(a) * e_r)
                    d = _ray_len(spin_m, o, o + sgn * 30.0 * e_t)
                    if d is not None:
                        best = min(best, d)
            res.append(best)
    t_shell = float(min(res))
    t = t_shell - FILLER_L                     # the 3 x 4 filler rides between the pin end and the spinner
    ear_t = PIN_HALF - EAR_IN
    core = float(hull_in(np.array([lay.hinge_x]), HUB_INSET)[0])
    return {"channel_tau_mm": [PIN_HALF, math.sqrt(core**2 - lay.hinge_r**2)],
            "pin_end_to_spinner_mm": t_shell, "filler_l_mm": FILLER_L,
            "free_travel_to_spinner_mm": t, "pin_end_to_spinner_all_mm": [round(v, 3) for v in res],
            "far_ear_thickness_mm": ear_t, "far_ear_engagement_left_mm": ear_t - t,
            "pin_stays_in_knuckle_and_near_ear": bool(t < PIN_HALF + P.KNUCKLE_L / 2),
            "note": "primary retention: press fit in both ears; secondary: a 3 x 4 filler dowel in each "
                    "insertion channel, trapped by the spinner shell, limits a loosened pin to this travel, "
                    "with both ears still engaged"}


def _rot_about(origin, axis, deg: float) -> np.ndarray:
    """4x4 rotation by deg (right-hand) about axis through origin."""
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    t = math.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    Rm = np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * K @ K
    T = np.eye(4)
    T[:3, :3] = Rm
    T[:3, 3] = np.asarray(origin, float) - Rm @ np.asarray(origin, float)
    return T


def cross_checks(hub_m, spin_m, d_deg: float = 2.0) -> dict:
    """Checks against the other builders' parts when their modules import: the real blade swept
    through its pivot (stowed -> deployed) vs hub and spinner, the shaft vs the hub bore, the
    spinner screw in the shaft end, and the aft shell vs the hub front (HUB_GAP)."""
    out = {}
    try:
        from . import blade as B
        bp = [p for p in B.parts() if p["name"] == "blade_0"][0]
        bm = to_manifold(to_mesh(bp["shape"], 0.02, 0.1))
        pv = bp["pivot"]
        worst = {"hub": (9e9, None), "spinner": (9e9, None)}
        for a in np.arange(0.0, pv["deployed_deg"] - 6.0 + 1e-9, d_deg):      # stowed .. 6 deg before stop
            m = _mf_transform(bm, _rot_about(pv["origin"], pv["axis"], a))
            for k, tgt in (("hub", hub_m), ("spinner", spin_m)):
                g = tgt.min_gap(m, 2.0)
                if g < worst[k][0]:
                    worst[k] = (g, float(pv["deployed_deg"] - a))
        dep = _mf_transform(bm, _rot_about(pv["origin"], pv["axis"], pv["deployed_deg"]))
        out["blade"] = {"source": "hydroegg.cad.blade blade_0 (stowed pose + pivot)",
                        "min_gap_hub_mm": worst["hub"][0], "at_fold_deg": worst["hub"][1],
                        "min_gap_spinner_mm": worst["spinner"][0], "spinner_at_fold_deg": worst["spinner"][1],
                        "deployed_overlap_hub_mm3": float((hub_m ^ dep).volume()),
                        "deployed_gap_hub_mm": float(hub_m.min_gap(dep, 2.0)),
                        "deployed_gap_spinner_mm": float(spin_m.min_gap(dep, 2.0)),
                        "folds_checked_deg": [6.0, float(pv["deployed_deg"]), d_deg],
                        "note": "fold < 6 deg excluded from the run clearance: the heel closes onto the stop "
                                "face there (deployed: contact, overlap must be 0)"}
    except Exception as e:  # noqa: BLE001
        out["blade"] = {"error": f"{type(e).__name__}: {e}"}
    try:
        from . import drivetrain as D
        dps = {p["name"]: p for p in D.parts()}
        sh = dps.get("shaft")
        if sh is not None:
            sm = to_manifold(to_mesh(sh["shape"], 0.005, 0.05))
            bb = sh["shape"].bounding_box()
            out["shaft"] = {"overlap_hub_mm3": float((hub_m ^ sm).volume()),
                            "min_gap_hub_mm": float(hub_m.min_gap(sm, 1.0)),
                            "shaft_aft_end_x_mm": float(bb.max.X), "hub_aft_face_x_mm": X_AFT}
            scr = to_manifold(to_mesh(spinner_screw(), 0.01))
            out["shaft"]["spinner_screw_in_shaft_overlap_mm3"] = float((scr ^ sm).volume())
            out["shaft"]["note"] = ("screw overlap = M3 thread modelled as D3 in the D2.5 tap drill (thread "
                                    "engagement), expected > 0")
    except Exception as e:  # noqa: BLE001
        out["shaft"] = {"error": f"{type(e).__name__}: {e}"}
    try:
        from . import shells as SH
        aft = [p for p in SH.parts() if p["name"] == "aft_shell"]
        if aft:
            am = to_manifold(to_mesh(aft[0]["shape"], 0.05, 0.2))
            out["aft_shell"] = {"min_gap_hub_mm": float(hub_m.min_gap(am, 2.0)),
                                "min_gap_spinner_mm": float(spin_m.min_gap(am, 2.0)),
                                "overlap_hub_mm3": float((hub_m ^ am).volume()),
                                "expected_gap_mm": P.HUB_GAP}
    except Exception as e:  # noqa: BLE001
        out["aft_shell"] = {"error": f"{type(e).__name__}: {e}"}
    return out


def run_checks(part_list: list[dict], meshes: dict) -> dict:
    t0 = time.time()
    _, lay, _, _, _ = layout_model()
    sb = slot_bands()
    out = {"schema": "hydroegg.cad.rotor/1", "generated_by": "python -m hydroegg.cad.rotor"}
    # ---- slot bands
    out["slot_bands"] = [{"name": b["name"], "tau_mm": [round(b["tau"][0], 3), round(b["tau"][1], 3)],
                          "x_range_mm": [round(v, 3) for v in b["x_range"]],
                          "floor_min_r_e_mm": round(b["floor_min"], 4),
                          "n_vertices": int(sum(len(p) for p in b["polys"]))} for b in sb["bands"]]
    b0 = sb["bands"][0]
    out["swept_root_tau_range_mm"] = [round(v, 3) for v in sb["tau_range"]]
    # ---- interfaces
    ic = []

    def add(name, value, limit, ok, note=""):
        ic.append({"check": name, "value": value, "limit": limit, "ok": bool(ok), "note": note})

    add("hub front x = STATIC_TAIL_X + HUB_GAP", X_FRONT, P.STATIC_TAIL_X + P.HUB_GAP,
        abs(X_FRONT - (P.STATIC_TAIL_X + P.HUB_GAP)) < 1e-9)
    add("hub aft face x = shaft end (SHAFT.x1)", X_AFT, P.SHAFT["x1"], abs(X_AFT - P.SHAFT["x1"]) < 1e-9)
    add("hub bore vs SHAFT_D: radial clearance (mm)", round(BORE_R - P.SHAFT_D / 2, 4), "0.05..0.1",
        0.05 <= BORE_R - P.SHAFT_D / 2 <= 0.1, "bore D%.2f close running on the D8 g6 shaft" % (2 * BORE_R))
    r_flat_shaft = P.SHAFT_D / 2 - P.SHAFT["hub_flat_depth"]
    add("D-flat clearance to the shaft flat (mm)", round(FLAT_R - r_flat_shaft, 4), "0.05..0.1",
        0.05 <= FLAT_R - r_flat_shaft <= 0.1, "D-flat x %.1f to %.1f faces the grub (psi %g)" % (*FLAT_X, GRUB_PSI))
    add("grub on the D-flat (x inside the flat, both)", GRUB_X, list(FLAT_X),
        FLAT_X[0] + 1.5 <= GRUB_X <= FLAT_X[1] - 1.5,
        "shaft flat must run from x <= %.1f out through the aft end (x %.0f) so the D-bore slides on" % FLAT_X)
    add("grub thread engagement in tapped PLA (mm)", round(r_flat_shaft + GRUB_L - BORE_R, 3), ">= 3",
        r_flat_shaft + GRUB_L - BORE_R >= 3.0, "M3 x 4 cup point from the flat (r 3.5) to r 7.5; tapped to r %.1f"
        % GRUB_TAP_R1)
    add("stop face x = HINGE_X - HEEL_T/2 = ROTOR_HUB.stop_face_x", STOP_X, lay.hinge_x - lay.heel_t / 2,
        abs(STOP_X - (lay.hinge_x - lay.heel_t / 2)) < 1e-9)
    add("stop face span r_e (mm)", [round(v, 3) for v in b0["stop_face_r_e"]], "heel face r_e 6.0..10.88",
        b0["stop_face_r_e"][0] <= lay.hinge_r - lay.heel_reach and b0["stop_face_r_e"][1] > 10.0,
        "flat stop face; its top is limited by the knuckle running clearance")
    add("ear inner face tau = KNUCKLE_L/2 + CLEAR_ROTATE", EAR_IN, 7.3, abs(EAR_IN - 7.3) < 1e-9)
    add("ear thickness = PIN_L/2 - ear inner face = HUB_EAR_T", round(PIN_HALF - EAR_IN, 4), P.HUB_EAR_T,
        abs(PIN_HALF - EAR_IN - P.HUB_EAR_T) < 1e-9)
    add("ear hole D (press, PRESS_FIT per side)", round(EAR_HOLE_D, 3), P.PIN_D, True,
        "ROTOR_HUB.ear_hole_d 3.0 is the nominal pin; printed bore %.2f" % EAR_HOLE_D)
    add("heel slot floor r_e min / wall over the shaft bore (mm)",
        [round(b0["floor_min"], 4), round(b0["floor_min"] - BORE_R, 4)], P.MIN_WALL,
        b0["floor_min"] - BORE_R >= P.MIN_WALL - 1e-6)
    add("magnet pocket: D / depth / face recess (mm)", [MAG_D_HOLE, MAG_H + MAG_RECESS, MAG_RECESS],
        "DETENTS hub magnet x 112.3", abs((STOP_X - MAG_RECESS - MAG_H / 2) - P.DETENTS["deployed_hub_magnet_x_re_tau"][0]) < 1e-6)
    eng = SCREW_L - AFT_WALL
    add("spinner screw M3 x 8: thread engagement in the shaft (mm)", eng, ">= 1.5 d = 4.5 and <= 8 (tap depth)",
        4.5 <= eng <= 8.0)
    add("access bore D vs M3 head D5.5 (radial clearance mm)", round(ACCESS_R - SCREW_HEAD_D / 2, 3), ">= 0.2",
        ACCESS_R - SCREW_HEAD_D / 2 >= 0.2)
    tg = tip_geometry()
    add("spinner tip: end x / face radius / tip round R", [round(tg["x_end"], 3), round(tg["r_face_top"], 3),
                                                          P.SPINNER_TIP_R], "r_face >= access r + chamfer + MIN_WALL",
        tg["r_face_top"] >= ACCESS_R + TIP_CHAMFER + P.MIN_WALL - 1e-6)
    hub_m = to_manifold(meshes["hub"])
    spin_m = to_manifold(meshes["spinner"])
    add("D-flat start to the grub bore (wall, mm)", round(GRUB_X - GRUB_TAP_D / 2 - FLAT_X[0], 3), ">= MIN_WALL",
        GRUB_X - GRUB_TAP_D / 2 - FLAT_X[0] >= P.MIN_WALL - 1e-9,
        "and %.2f mm aft of the 0.5 mm bore chamfer" % (FLAT_X[0] - (X_FRONT + BORE_CHAMFER)))
    add("front V-grooves inside the grub tapped zone: wall to the grub bore (mm)",
        round(GRUB_X - GRUB_TAP_D / 2 - (X_FRONT + GROOVE_DEPTH), 3), ">= MIN_WALL",
        max(GROOVES_R) + GROOVE_DEPTH <= GRUB_TAP_R1 and GRUB_X - GRUB_TAP_D / 2 - (X_FRONT + GROOVE_DEPTH) >= P.MIN_WALL)
    pr = pin_retention(spin_m)
    out["pin_retention"] = pr
    add("pins retained: free axial travel (pin + 3 x 4 filler to the spinner) / far-ear engagement left (mm)",
        [round(pr["free_travel_to_spinner_mm"], 3), round(pr["far_ear_engagement_left_mm"], 3)],
        "engagement left >= 1.0", pr["far_ear_engagement_left_mm"] >= 1.0, pr["note"])
    out["interface_checks"] = ic
    # ---- geometry checks on the meshes
    out["hub_spinner_overlap_mm3"] = float((hub_m ^ spin_m).volume())
    out["hub_spinner_min_gap_mm"] = float(hub_m.min_gap(spin_m, 1.0))
    print("  sweep check...", flush=True)
    out["swept_clearance"] = swept_clearance(hub_m, spin_m)
    print("  slide-on check...", flush=True)
    out["spinner_slide_on"] = so = spinner_slide_on(spin_m)
    ok_folds = [k for k, v in so.items() if v["min_gap_mm"] >= RUN - 0.01]
    add("spinner slides on from aft over the fitted blades (folds with >= run clearance all the way)",
        ok_folds, "assembly fold", "fold_90" in ok_folds,
        "deployed blades (fold 0) block it %.1f mm before home: fold the blades (stowed) first"
        % so["fold_0"]["at_offset_mm"] if so["fold_0"]["min_gap_mm"] < RUN - 0.01 else "")
    print("  wall thickness...", flush=True)
    out["wall_thickness"] = {k: wall_thickness(to_mesh(p["shape"], 0.005, 0.05))
                             for k, p in ((q["name"], q) for q in part_list) if k in ("hub", "spinner")}
    out["shell_profile"] = shell_thickness()
    if "--no-cross" not in sys.argv:
        print("  cross checks (blade, shaft, aft shell)...", flush=True)
        out["cross_checks"] = cross_checks(hub_m, spin_m)
    # bought parts vs printed parts (pins in ears, magnets in pockets, grub, screw)
    bought = {}
    for p in part_list:
        if p["bought"]:
            bm = to_manifold(to_mesh(p["shape"], 0.01))
            bought[p["name"]] = {"overlap_hub_mm3": round(float((bm ^ hub_m).volume()), 4),
                                 "overlap_spinner_mm3": round(float((bm ^ spin_m).volume()), 4)}
    out["bought_vs_printed_overlap"] = bought
    out["elapsed_s"] = round(time.time() - t0, 1)
    return out


# --------------------------------------------------------------------------
# Main: export STEP / STL, checks, previews
# --------------------------------------------------------------------------

def _round(o, nd=4):
    if isinstance(o, float):
        return round(o, nd)
    if isinstance(o, dict):
        return {k: _round(v, nd) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_round(v, nd) for v in o]
    if isinstance(o, (np.floating,)):
        return round(float(o), nd)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main(argv: list[str]) -> None:
    import trimesh
    t0 = time.time()
    for d in (OUT_STEP, OUT_STL, OUT_CHK, WORK):
        d.mkdir(parents=True, exist_ok=True)
    print("slot bands...", flush=True)
    sb = slot_bands()
    print("  %d cloud points, tau %.2f..%.2f (%.1f s)" % (sb["n_cloud_pts"], *sb["tau_range"], time.time() - t0))
    print("building parts...", flush=True)
    ps = parts()
    print("  built (%.1f s)" % (time.time() - t0))
    report = {"parts": []}
    meshes = {}
    for p in ps:
        s = p["shape"]
        rec = check(s, p["name"])
        export_step(s, OUT_STEP / f"{p['name']}.step")
        if p["print"] is not None:
            m = to_mesh(p["print"], 0.02, 0.1)
            m.export(OUT_STL / f"{p['name']}.stl")
            mr = trimesh.load(OUT_STL / f"{p['name']}.stl")
            rec.update(watertight=bool(mr.is_watertight), stl_volume_cm3=round(float(mr.volume) / 1000, 3),
                       stl_faces=int(len(mr.faces)), print_size_mm=check(p["print"])["size_mm"],
                       fits_bed=fits_bed(p["print"]), overhang=overhang_report(mr))
            meshes[p["name"]] = to_mesh(s, 0.02, 0.1)
        else:
            mr = to_mesh(s, 0.01)
            rec.update(watertight=bool(mr.is_watertight))
        rec.update(bought=p["bought"], material=p["material"], qty=p["qty"])
        report["parts"].append(rec)
        print("  %-16s valid=%s solids=%d vol=%.3f cm3 watertight=%s" % (p["name"], rec["valid"], rec["n_solids"],
                                                                         rec["volume_cm3"], rec["watertight"]))
    print("checks...", flush=True)
    report.update(run_checks(ps, meshes))
    (OUT_CHK / "rotor.json").write_text(json.dumps(_round(report), indent=1))
    print("wrote", OUT_CHK / "rotor.json", "(%.1f s)" % (time.time() - t0))
    if "--no-preview" not in argv:
        render_previews(PREVIEWS)


PREVIEWS = [("preview_assembly_deployed", "rotor_deployed.png", "view=-35,22&fold=0&zoom=1.5"),
            ("preview_assembly_stowed", "rotor_stowed.png", "view=-150,20&fold=0&zoom=1.4"),
            ("inspect_hub", "hub_front.png", "view=-140,25&zoom=1.5"),
            ("inspect_hub", "hub_aft.png", "view=40,35&zoom=1.5"),
            ("inspect_spinner", "spinner.png", "view=-150,30&zoom=1.4"),
            ("inspect_section", "section.png", "view=side&zoom=1.3"),
            ("preview_print_poses", "print_poses.png", "view=-30,30&zoom=1.2")]


def render_previews(jobs, out_dir: Path = WORK) -> list[Path]:
    """Render preview functions of this module in-process (reuses the build cache): writes a GLB +
    manifest under out/work/preview and screenshots it through scripts/render.py."""
    import subprocess
    from ..export_web import write_scene
    outs = []
    for fn, png, q in jobs:
        ps = globals()[fn]()
        stem = Path(png).stem
        write_scene(ps, REPO / "out/work/preview", glb_name=f"rotor_{stem}.glb")
        qq = f"m=../out/work/preview/rotor_{stem}.json&{q}"
        out = out_dir / png
        subprocess.run([sys.executable, str(REPO / "scripts/render.py"), str(out), qq, "-W", "1400", "-H", "900"],
                       check=True, cwd=REPO)
        outs.append(out)
    return outs

if __name__ == "__main__":
    main(sys.argv[1:])
