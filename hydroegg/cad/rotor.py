"""Rotor: hub, spinner, spinner tip cap, hinge pins, hub detent magnets, hub grub screw and spinner screw.

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
    The pin goes in LAST: drop the blade into its slot, then press the pin through the channel.
    Deployed-detent magnets (finding T03): a D4 x 2 disc did not fit the 4.5 mm stop zone (0.15 mm
    to the knuckle arc). Glued 6.1 x 2.1 x 2.2 seats for N52 blocks 6 x 2 x 2 (6 along the pin, 2
    radial, magnetised along x) centred at r_e 8.2 (common.DEP_MAGNET), 1.27 mm to the knuckle arc.
    M3 x 4 cup-point grub, radial at x 111.6, psi 60: tapped M3 in the PLA (2.5 drill, r 4.08 to
    9.5) with a D3.4 clearance bore out to the surface. A heat-set insert does not fit: its D4.0
    hole needs 1.6 mm of wall and leaves 1.0 to the front face and about 1.3 to the swept roots,
    and it would sit 15 mm deep. PLA tapping is enough here: the D-flat carries the torque, the
    rotor thrust goes through the spinner screw, the grub only sets the hub's axial position.
    Front face: flat slinger face 0.6 mm (HUB_GAP) from the static tail end face; the v1 labyrinth
    V-grooves were removed (finding T11: first-layer squish closes them and the r 8.7 groove stood
    1.0 mm from the grub thread). Nothing stands proud, so the hub prints on it. Three D3 x 1.5
    pegs on the aft face key the spinner. (A 2 mm trim of the hub's outer skin to the window
    outline was tried and dropped, see hub_trim0.)

Spinner (printed, rotating, x 108.6 to the joint plane X_JOINT 151.725)
    1.2 mm shell (normal) over the hub with a 0.2 mm gap. Inner aft wall (3 mm): a flat D22.6 seat
    ceiling around the screw hole (the clamp face on the hub and the 3 peg holes), then a 45 deg
    cone that matches the hub's aft cone; clamped against the hub aft face and the shaft end by
    one ISO 4762 M3 x 8 into the tapped shaft end (rotor thrust: hub -> spinner wall -> screw ->
    shaft). Windows (finding T01, see spinner_window0): one smooth outline per blade (convex hull
    of the slot bands + 1.5 offset, arc joins) with a 45 deg gable top; full depth only where the
    heel and knuckle sweep (it closes with 45 deg flanks and a <= 15 mm bridge at x ~125.9, leaving
    a complete ring of the cone under the seat ceiling), elsewhere a recess above a 1.2 mm skin
    whose floor rises at 45 deg out through the surface. Aft of the wall the spinner is solid
    around the D6 access bore, which opens into a 45 deg countersink at the joint plane.
    Tip cap (finding T02): a printed push-in cap continues profile.radius from X_JOINT to an
    R SPINNER_TIP_R round on the axis (built tip x ~161, L/D 2.15): D5.7 x 12 spigot with three
    crush ribs to D6.1 in the access bore, 45 deg shoulder in the countersink, and a fingernail
    V-groove at the joint (0.8 x 45 deg chamfers on both sides) to pull it for the M3 screw.
    Slide the spinner on with the blades folded >= 15 deg (stowed is simplest).

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
from .common import DEP_MAGNET, axis_cyl, check, cyl, fits_bed, on_bed

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
MAG_H = P.MAGNET["H"]
MAG_RECESS = 0.2                               # layout magnet_recess (DETENTS hub magnet x 112.3)
# Deployed-detent pair (finding T03): a D4 x 2 disc does not fit the 4.5 mm stop zone (0.15 mm to
# the knuckle arc on the hub, 0.39 to the heel tip on the blade). Both sides take an N52 BLOCK
# 6 x 2 x 2 (6 along the pin, 2 radial, magnetised through 2 along x) centred at r_e DEP_MAG_RE,
# CA-glued in DEP_SEAT seats; each face sits MAG_RECESS below its stop face, so the deployed pair
# gap is 2 x MAG_RECESS = 0.4 (DETENTS design gap). PROPOSED for params DETENTS (not edited).
DEP_MAG = tuple(DEP_MAGNET["size_tau_re_x"])   # tau, r_e, x (magnetisation) mm
DEP_SEAT = tuple(DEP_MAGNET["seat_tau_re"])     # glued seat (tau, r_e)
DEP_MAG_RE = DEP_MAGNET["r_e"]                  # r_e of the block centre (DETENTS had 8.4388 for the disc)
assert abs(MAG_RECESS - DEP_MAGNET["recess"]) < 1e-12
GRUB_TAP_D, GRUB_CLEAR_D, GRUB_TAP_R1 = 2.5, 3.4, 9.5
GRUB_L = 4.0                                   # M3 x 4 (ROTOR_HUB.grub)
SCREW_L, SCREW_HEAD_D, SCREW_HEAD_H = 8.0, 5.5, 3.0   # ISO 4762 M3 x 8
SCREW_CLEAR_D = 3.4
ACCESS_R = P.SPINNER_ACCESS["d"] / 2           # 3.0
AFT_WALL = float(P.SPINNER_ACCESS["aft_wall"]) # 3.0
X_HEAD = X_AFT + AFT_WALL                      # 131: screw head seat
# Tip (finding T02): the spinner ends at the joint plane X_JOINT in a 45 deg countersink around
# the D6 access bore; a printed push-in tip cap continues profile.radius to a SPINNER_TIP_R round.
X_JOINT = 151.725                              # DESIGN: old blunt end plane (L/D unchanged up to here)
CAP_SPIGOT_D, CAP_SPIGOT_L = 5.7, 12.0         # 0.15 radial in the D6 bore, 3 crush ribs to D6.1
CAP_RIB_R, CAP_RIB_W = 3.05, 0.6               # rib tips (0.05 radial interference) and width
CAP_GROOVE = 0.8                               # fingernail V-groove: 45 deg chamfers on both sides of the joint
PEG_D, PEG_H, PEG_R = 3.0, 1.5, 8.0            # DESIGN: anti-rotation pegs on the hub aft face
PEG_PSI = (60.0, 180.0, 300.0)                 # between the blades
PEG_HOLE_D, PEG_HOLE_DEPTH = PEG_D + 2 * P.CLEAR_SLIDE, PEG_H + 0.3
FILLER_L = 4.0                                 # DESIGN: ISO 8734 3 x 4 filler dowel in each pin channel
# (front-face labyrinth V-grooves removed for v1, finding T11: first-layer squish closes them and
#  the r 8.7 groove stood 1.0 mm from the grub thread)
FRONT_CHAMFER = 0.4                            # elephant-foot relief on the bed edges
BORE_CHAMFER = 0.5
DX = 0.02                                      # floor sampling step along x
R_TOP = 40.0                                   # slot cut tops (outside everything)

MATERIAL_PRINTED = "PETG (v2, wet) / PLA (display)"

SUPPORT_NOTES = {
    "hub": ("front face on the bed (+X up). No supports: outer flank leans <= 24 deg; slot floors "
            "rise <= 45 deg; slot ends are flat bridges across the slot (<= 14.6 mm); the ear holes, "
            "pin channels and grub bore are horizontal holes D <= 3.5 (arched tops); magnet pockets "
            "and the bore open upward; the aft cone is a 45 deg roof; the front face is flat (no grooves)."),
    "spinner": ("front edge (base) on the bed, tip up. No supports: shell leans <= 26 deg; the inner "
                "aft wall is a 45 deg cone plus a flat D22.6 bridge (the hub clamp face) that sits on a "
                "complete ring of that cone; each window's full-depth part closes with 45 deg flanks "
                "and a flat bridge <= 15 mm anchored on both flanks, its shallow part ends in a 45 deg "
                "floor that runs out through the surface; the access bore and its 45 deg countersink "
                "are vertical; brim recommended (bed contact is the 1.2 mm front edge)."),
    "spinner_cap": ("No supports: D5.7 spigot (3 crush ribs with 45 deg lead-ins) on the bed, 45 deg "
                    "shoulder, outer surface leans <= 31 deg; 3 mm brim (25 mm2 footprint, 21 mm tall)."),
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


# ---- spinner windows ----------------------------------------------------------------------
# One outline per blade in the (x, tau) plane, extruded along e_r. The spinner prints base down
# (+X up), so a window's aft end is the part that must be support-free.
WIN_OFFSET = 1.5          # DESIGN: outline = convex hull of the band rectangles, offset outward (arc joins)
WIN_DEEP_AFT = 0.5        # DESIGN: aft clearance of the full-depth part (its gable) beyond the deep band ends
WIN_BRIDGE = 15.0         # DESIGN: longest flat bridge across a window top (print)
WIN_SKIN_N = P.MIN_WALL   # skin left (normal) under the shallow part of a window
WIN_X_THIN = 119.5        # forward of this the spinner is too thin for a skin + recess: full depth
WIN_CHAMFER = 0.4         # outer window edge
HUB_TRIM = 2.0            # hub outer skin trimmed with the full-depth outline (no slot teeth show)


def r_inner_spinner(x):
    """Spinner inner surface radius: hull inset SP_WALL, then the 45 deg cone to the aft wall."""
    x = np.asarray(x, float)
    xci = x_cone_spinner()
    r_ceil = HUB_AFT_FLAT_R + SP_GAP * math.sqrt(2)
    return np.where(x < xci, hull_in(np.minimum(x, xci), SP_WALL), r_ceil + (X_AFT - x))


def _xt_face(pts) -> Face:
    """Polygon of (x, tau) points in the plane r_e = -1 of blade 0 (tau = -y at psi 0)."""
    vs = [Vector(float(x), -float(t), -1.0) for x, t in pts]
    return Face(Wire.make_polygon(vs, close=True))


def _faces_of(s) -> list:
    """Faces of a boolean result (Shape, Sketch, Compound or ShapeList)."""
    if isinstance(s, (list, tuple)):
        return [f for x in s for f in _faces_of(x)]
    return list(s.faces()) if hasattr(s, "faces") else [s]


def _gable(pts_rad, x_min: float, bridge: float) -> dict:
    """Smallest 45 deg gable roof with a flat top of at most ``bridge`` (mm) over points (x, tau, r)
    grown by r: left flank x = cL + tau, right flank x = cR - tau, flat top at xa."""
    cL = max(x - t + r * math.sqrt(2) for x, t, r in pts_rad)
    cR = max(x + t + r * math.sqrt(2) for x, t, r in pts_rad)
    xa = max((cL + cR - bridge) / 2, x_min)
    return {"cL": cL, "cR": cR, "xa": xa, "tau_a": (xa - cL, cR - xa)}


def _gable_face(g: dict, x_low: float = 90.0, half: float = 40.0) -> Face:
    """Region under a gable (big polygon) in the window plane."""
    cL, cR, xa = g["cL"], g["cR"], g["xa"]
    return _xt_face([(x_low, -half), (cL - half, -half), (xa, xa - cL), (xa, cR - xa), (cR - half, half),
                     (x_low, half)])


@functools.lru_cache(maxsize=1)
def window_layout() -> dict:
    """Window geometry of blade 0 (see spinner_window0): outlines, gables, skin floor profile."""
    from scipy.spatial import ConvexHull
    from build123d import Kind, offset
    sb = slot_bands()
    xs = sb["xs"]
    rects, deep = [], []
    rS0 = r_inner_spinner(xs) + WIN_SKIN_N * np.sqrt(1 + np.gradient(r_inner_spinner(xs), xs) ** 2)
    T = np.full(xs.size, np.inf)                       # shallow (non-deep) band floors, as radius
    for b in sb["bands"]:
        t0, t1 = b["tau"]
        tmin = 0.0 if t0 <= 0.0 <= t1 else min(abs(t0), abs(t1))
        g = b["floor"]
        fin = np.isfinite(g)
        fr = np.where(fin, np.sqrt(np.maximum(g, 0.0) ** 2 + tmin * tmin), np.inf)
        is_deep = fin & (fr < rS0 + 0.05)
        for poly in b["polys"]:
            rects.append((b["name"], float(poly[:, 0].min()), float(poly[:, 0].max()), t0, t1))
        idx = np.nonzero(is_deep)[0]
        shallow = fin & ~is_deep
        # a band whose shallow tail stands above the hull surface needs no cut there (b3)
        shallow &= fr < hull_r(xs) - 0.05
        if idx.size:
            x1d = float(xs[idx.max()]) if shallow.any() else float(xs[np.nonzero(fin)[0].max()])
            deep.append((b["name"], float(xs[idx.min()]), x1d, t0, t1))
        T = np.where(shallow, np.minimum(T, fr), T)
    # skin floor of the shallow part: >= the skin, <= every shallow band floor, rising at most 45 deg
    # going aft (up in print); after the last band it keeps rising and leaves through the surface
    has = np.isfinite(T)
    first = int(np.argmax(has))
    Tt = np.where(np.arange(xs.size) < first, rS0, T)
    F = Tt.copy()
    for i in range(1, xs.size):
        F[i] = min(F[i], F[i - 1] + DX)
    rS = np.maximum(rS0, F)
    # outer outline P (hull + offset) and its roof
    corners = np.array([(x, t) for _, x0, x1, t0, t1 in rects for x in (x0, x1) for t in (t0, t1)])
    hv = corners[ConvexHull(corners).vertices]
    Pf = offset(_xt_face(hv), WIN_OFFSET, kind=Kind.ARC)
    gS = _gable([(x, t, WIN_OFFSET) for x, t in hv], float(hv[:, 0].max()) + WIN_OFFSET, WIN_BRIDGE)
    # full-depth part D: hull of the deep band rectangles (and of the whole outline forward of
    # WIN_X_THIN, where the spinner is too thin for a skin) + offset, aft end cut by its own gable
    dc = [(x, t) for _, x0, x1, t0, t1 in deep for x in (x0, x1) for t in (t0, t1)]
    for a, b in zip(hv, np.roll(hv, -1, axis=0)):                  # hull edges crossing x = WIN_X_THIN
        if (a[0] - WIN_X_THIN) * (b[0] - WIN_X_THIN) < 0:
            u = (WIN_X_THIN - a[0]) / (b[0] - a[0])
            dc.append((WIN_X_THIN, float(a[1] + u * (b[1] - a[1]))))
    dc = np.array(dc)
    dv = dc[ConvexHull(dc).vertices]
    gD = _gable([(x, t, WIN_DEEP_AFT) for x, t in dv], float(dv[:, 0].max()) + WIN_DEEP_AFT, WIN_BRIDGE)
    return dict(xs=xs, rects=rects, deep=deep, hull=hv, deep_hull=dv, P=Pf, gable_S=gS, gable_D=gD,
                rS0=rS0, rS=rS, T=T)


def _slab(face_or_faces) -> Solid:
    out = [Solid.extrude(f, Vector(0, 0, R_TOP + 1.0)) for f in _faces_of(face_or_faces)]
    return _fuse(out) if len(out) > 1 else out[0]


@functools.lru_cache(maxsize=1)
def window_parts0() -> dict:
    """Solids of the blade-0 window: 'deep' (full-depth slab D) and 'shallow' (the outer outline
    slab above the skin floor), plus their 2D faces."""
    W = window_layout()
    gS, gD = W["gable_S"], W["gable_D"]
    hv = W["hull"]
    # outer outline: P plus the gable's flat top. The flanks touch P where its offset arcs turn
    # 45 deg (support points of the hull grown along (1, -+1)/sqrt 2).
    vL = hv[int(np.argmax(hv[:, 0] - hv[:, 1]))]
    vR = hv[int(np.argmax(hv[:, 0] + hv[:, 1]))]
    s2 = WIN_OFFSET / math.sqrt(2)
    tL, tR = (vL[0] + s2, vL[1] - s2), (vR[0] + s2, vR[1] + s2)
    cen = (float(hv[:, 0].mean()), float(hv[:, 1].mean()))
    roof = _xt_face([tL, (gS["xa"], gS["tau_a"][0]), (gS["xa"], gS["tau_a"][1]), tR, cen])
    outer = _faces_of(W["P"].fuse(roof).clean())
    deep_f = _faces_of(offset_face(_xt_face(W["deep_hull"]), WIN_OFFSET).intersect(_gable_face(gD)))
    # skin floor (revolved keep region) under the shallow part
    xs, rS = W["xs"], W["rS"]
    k = np.nonzero(rS < hull_r(xs) + 2.0)[0]
    x_end = float(xs[k.max()])
    sel = (xs >= WIN_X_THIN - 1.0) & (xs <= x_end)
    pts = _rdp(np.column_stack([xs[sel], rS[sel]]), 0.003)
    prof = [(float(pts[0, 0]), 0.0)] + [(float(x), float(r)) for x, r in pts] + \
           [(x_end, 45.0), (145.0, 45.0), (145.0, 0.0)]
    keep = _revolve([("line", prof + [prof[0]])])
    shallow = _slab(outer).cut(keep)
    deep = _slab(deep_f)
    return dict(outer_face=outer, deep_face=deep_f, deep=deep, shallow=shallow, keep=keep,
                tangent_points=(tL, tR), x_skin_exit=x_end)


def offset_face(f: Face, d: float):
    from build123d import Kind, offset
    return offset(f, d, kind=Kind.ARC)


@functools.lru_cache(maxsize=1)
def spinner_window0() -> Solid:
    """Spinner window for blade 0 (finding T01: one smooth outline, support-free aft end).

    Outline at the outer surface: the convex hull of the hub-slot band rectangles (x, tau), offset
    WIN_OFFSET outward with arc joins (no corner gives material back to the swept root), with the
    aft side closed by a 45 deg gable whose flat top is at most WIN_BRIDGE wide. Depth:
      * full depth (r_e >= -1) only where a band floor lies below the spinner's inner surface + the
        skin (the heel and knuckle sweep, x <= ~124.7) and everywhere forward of WIN_X_THIN; that
        region D is the hull of the deep band parts + WIN_OFFSET, its aft end cut by a 45 deg
        gable (flat top <= WIN_BRIDGE) over the deep corners + WIN_DEEP_AFT, so it closes at
        x ~125.9 and a complete ring of the 45 deg cone stands under the D22.6 seat ceiling;
      * elsewhere the window cuts only above a revolved skin floor: >= WIN_SKIN_N (normal) over the
        inner surface, <= every shallow band floor (run clearance kept), rising at most 45 deg going
        aft, and after the last band it rises on through the outer surface, which closes the
        window (no ceiling at all)."""
    w = window_parts0()
    return _fuse([w["deep"], w["shallow"]])


def hub_trim0() -> Solid:
    """Hub outer HUB_TRIM mm inside the full-depth window outline (blade 0), proposed by finding T01
    so no slot steps show through the window. Built and checked, then NOT applied to the hub: the
    band slots run the full depth, so their steps still showed below the 2 mm trim, while the trim's
    constant-tau side walls met the hub cone at a grazing angle (1.08 mm walls) and the guard that
    kept it off the pin channel left a 12 mm2 cantilever ledge. The smooth spinner outline, 1.5 mm
    outside the slots, is the visible fix."""
    w = window_parts0()
    xc = x_cone_hub()
    rc = lambda x: float(hull_in(x, HUB_INSET + HUB_TRIM))  # noqa: E731
    xs = np.linspace(X_FRONT - 1.0, xc, 30)
    r_cc = float(hull_in(xc, HUB_INSET)) - HUB_TRIM * math.sqrt(2)
    prof = [(X_FRONT - 1.0, 0.0)] + [(float(x), rc(max(x, X_FRONT))) for x in xs[:-1]]
    prof += [(xc, min(r_cc, rc(xc))), (X_AFT + 2.0, min(r_cc, rc(xc)) - (X_AFT + 2.0 - xc)), (X_AFT + 2.0, 0.0)]
    keep = _revolve([("line", prof + [prof[0]])])
    # keep the hub whole around the +tau pin insertion channel (D3.5 along the pin, tau >= PIN_HALF):
    # a trim surface crossing the channel at a grazing angle would leave a membrane over it
    Hp, _, _ = hinge_frame(0.0)
    m = CHANNEL_D / 2 + P.MIN_WALL
    guard = _box_hinge(Hp[0] - m, Hp[0] + m, -5.0, R_TOP, PIN_HALF - m, R_TOP)
    return w["deep"].cut(keep).cut(guard)


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


def _pin_holes0() -> Solid:
    """Ear hole (press) and both insertion channels for blade 0, along e_t."""
    Hp, e_r, e_t = hinge_frame(0.0)
    Hp = np.asarray(Hp)
    ear = cyl(Hp - (PIN_HALF + 0.05) * e_t, Hp + (PIN_HALF + 0.05) * e_t, EAR_HOLE_D / 2)
    ch1 = cyl(Hp + (PIN_HALF - 0.05) * e_t, Hp + 24.0 * e_t, CHANNEL_D / 2)
    ch2 = cyl(Hp - (PIN_HALF - 0.05) * e_t, Hp - 24.0 * e_t, CHANNEL_D / 2)
    return _fuse([ear, ch1, ch2])


def magnet_center0() -> np.ndarray:
    """Seated hub magnet centre (blade 0): x from STOP_X - recess - H/2, r_e DEP_MAG_RE, tau 0."""
    _, e_r, e_t = hinge_frame(0.0)
    xm = STOP_X - MAG_RECESS - DEP_MAG[2] / 2
    return np.array([xm, 0, 0]) + DEP_MAG_RE * np.asarray(e_r)


def _box_hinge(x0: float, x1: float, re0: float, re1: float, t0: float, t1: float) -> Solid:
    """Box in the blade-0 hinge frame (x, r_e = +Z, tau = -Y)."""
    return Solid.make_box(x1 - x0, t1 - t0, re1 - re0).moved(Location((x0, -t1, re0)))


def _magnet_pocket0() -> Solid:
    """Glued seat for the 6 x 2 x 2 block (6 along the pin, 2 radial, magnetised along x): DEP_SEAT
    (+0.1 on each size), depth H + MAG_RECESS from the stop face, open 0.3 past it."""
    c = magnet_center0()
    sx, sr = DEP_SEAT[0] / 2, DEP_SEAT[1] / 2
    return _box_hinge(STOP_X - MAG_RECESS - DEP_MAG[2], STOP_X + 0.3, c[2] - sr, c[2] + sr, -sx, sx)


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
    tools = [_d_bore(), _grub_hole()]                      # front-face V-grooves dropped for v1 (T11)
    cut0 = _fuse([slot_cut0(), _pin_holes0(), _magnet_pocket0()])    # (hub_trim0 not applied: see its docstring)
    tools += [rot_psi(cut0, psi) for psi in psi_list()]
    body = body.cut(_fuse(tools))
    body = _fuse([body] + _pegs())
    return _single(body)


@functools.lru_cache(maxsize=1)
def tip_geometry() -> dict:
    """Spinner / tip-cap joint (finding T02) and the cap's tip round.

    Joint plane X_JOINT. The spinner's access bore (D6) opens into a 45 deg countersink that
    reaches the hull radius at X_JOINT; the cap's 45 deg shoulder seats in it (support-free on
    both prints). A fingernail V-groove: the spinner's outer edge is chamfered 45 deg from
    X_JOINT - CAP_GROOVE down to the countersink, the cap's edge 45 deg back out to the hull.
    The cap's outer surface continues profile.radius to a SPINNER_TIP_R round on the axis."""
    xj = X_JOINT
    rj = float(hull_r(xj))
    x_cs0 = xj - (rj - ACCESS_R)                          # countersink meets the bore
    xa = xj - CAP_GROOVE
    ra = float(hull_r(xa))
    xb = (ra + xa - ACCESS_R + x_cs0) / 2                 # spinner chamfer meets the countersink
    rb = ACCESS_R + xb - x_cs0
    xc_ = _solve(lambda x: float(hull_r(x)) - (rb + x - xb), xb, xb + 3.0)   # cap chamfer meets the hull
    rc_ = float(hull_r(xc_))
    Rt = P.SPINNER_TIP_R
    xx = np.linspace(150.0, P.HULL_L, 20001)
    rr = hull_r(xx)
    xcen = _solve(lambda xc: float(np.min(np.hypot(xx - xc, rr))) - Rt, 150.0, P.HULL_L - Rt * 0.5)
    k = int(np.argmin(np.hypot(xx - xcen, rr)))
    xt, rt = float(xx[k]), float(rr[k])
    a0 = math.atan2(rt, xt - xcen)
    mid = (xcen + Rt * math.cos(a0 / 2), Rt * math.sin(a0 / 2))
    return {"x_joint": xj, "r_joint": rj, "x_countersink": x_cs0, "spinner_chamfer": [(xa, ra), (xb, rb)],
            "cap_chamfer": [(xb, rb), (xc_, rc_)], "tip_center_x": xcen, "tangent": (xt, rt), "mid": mid,
            "x_tip": xcen + Rt, "x_end": xj, "groove_depth": float(hull_r(xb)) - rb}


def x_cone_spinner() -> float:
    off = HUB_AFT_FLAT_R + SP_GAP * math.sqrt(2)
    return _solve(lambda x: hull_in(x, SP_WALL) - (off + X_AFT - x), 112.0, X_AFT)


@functools.lru_cache(maxsize=1)
def spinner_body() -> Solid:
    ro = lambda x: float(hull_r(x))              # noqa: E731
    ri = lambda x: float(hull_in(x, SP_WALL))    # noqa: E731
    tg = tip_geometry()
    (xa, ra), (xb, rb) = tg["spinner_chamfer"]
    xci = x_cone_spinner()
    r_ceil = HUB_AFT_FLAT_R + SP_GAP * math.sqrt(2)
    x1 = X_FRONT + 0.3
    segs = [("line", [(X_FRONT, ri(X_FRONT)), (X_FRONT, ro(X_FRONT) - 0.3), (x1, ro(x1))]),
            ("spline", _spline_pts(ro, x1, xa, 60)),
            ("line", [(xa, ra), (xb, rb), (tg["x_countersink"], ACCESS_R),
                      (X_HEAD, ACCESS_R), (X_HEAD, SCREW_CLEAR_D / 2), (X_AFT, SCREW_CLEAR_D / 2),
                      (X_AFT, r_ceil), (xci, ri(xci))]),
            ("spline", list(reversed(_spline_pts(ri, X_FRONT, xci, 40))))]
    return _revolve(segs)


@functools.lru_cache(maxsize=1)
def spinner_cap() -> Solid:
    """Printed push-in tip cap (T02): D5.7 spigot with 3 crush ribs (D6.1 over the tips) in the
    access bore, 45 deg shoulder in the spinner's countersink, outer surface = profile.radius to an
    R SPINNER_TIP_R round. Pull it with a fingernail in the joint groove to reach the M3 screw."""
    tg = tip_geometry()
    r_sp = CAP_SPIGOT_D / 2
    x_sh = tg["x_countersink"] - (ACCESS_R - r_sp)        # spigot meets the 45 deg shoulder
    x_s0 = x_sh - CAP_SPIGOT_L
    (xb, rb), (xc_, rc_) = tg["cap_chamfer"]
    xt, rt = tg["tangent"]
    ro = lambda x: float(hull_r(x))              # noqa: E731
    # R SPINNER_TIP_R round; it ends in a D0.5 flat on the axis (a surface of revolution with a pole on
    # the axis tessellates with a crack there)
    Rt, xcen = P.SPINNER_TIP_R, tg["tip_center_x"]
    rt = math.sqrt(Rt ** 2 - (xt - xcen) ** 2)                 # exactly on the round
    a0, a1 = math.atan2(rt, xt - xcen), math.asin(0.25 / Rt)
    pend = (xcen + Rt * math.cos(a1), 0.25)
    pmid = (xcen + Rt * math.cos(0.5 * (a0 + a1)), Rt * math.sin(0.5 * (a0 + a1)))
    hull_pts = _spline_pts(ro, xc_, xt, 24)
    hull_pts[-1] = (xt, rt)
    segs = [("line", [(x_s0, 0.0), (x_s0, r_sp - 0.3), (x_s0 + 0.3, r_sp), (x_sh, r_sp), (xc_, rc_)]),   # (xb, rb) is collinear
            ("spline", hull_pts),
            ("arc", [(xt, rt), pmid, pend]),
            ("line", [pend, (pend[0], 0.0), (x_s0, 0.0)])]
    body = _revolve(segs)
    xr0, xr1, ri0 = x_s0 + 1.2, x_sh - 0.2, r_sp - 0.3
    pts = [(xr0, ri0), (xr1, ri0), (xr1, CAP_RIB_R), (xr0 + 0.55, CAP_RIB_R), (xr0, CAP_RIB_R - 0.55)]
    vs = [Vector(float(x), -CAP_RIB_W / 2, float(r)) for x, r in pts]
    rib = Solid.extrude(Face(Wire.make_polygon(vs, close=True)), Vector(0, CAP_RIB_W, 0))   # 45 deg lead-in
    ribs = [rot_psi(rib, a) for a in (60.0, 180.0, 300.0)]      # clear of the revolve seam (+Z)
    return _single(_fuse([body] + ribs))




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
    """N52 block 6 x 2 x 2 (DEP_MAG): x along the magnetisation, centred at the origin, 0.1 edge chamfer."""
    t, r, h = DEP_MAG
    b = Solid.make_box(h, t, r).moved(Location((-h / 2, -t / 2, -r / 2)))
    try:
        b = b.chamfer(0.1, None, b.edges())
    except Exception:  # noqa: BLE001
        pass
    return b


def hub_magnet(i: int) -> Solid:
    """Block seated on its seat floor (face MAG_RECESS below the stop face)."""
    c = magnet_center0()
    return rot_psi(magnet_model().moved(Location(tuple(float(v) for v in c))), psi_list()[i])


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
                          + ", 4 perimeters, 40 % infill. Tap the grub bore M3 after printing; glue the "
                          "three 6 x 2 x 2 magnets into the stop-face seats (CA). The hinge pins go in "
                          "only after the blades: drop each blade into its slot, then press its pin "
                          "through the D3.5 channel into both ears."))
    sp = spinner()
    out.append(dict(name="spinner", shape=sp, print=_print_pose(sp), material="hub", role="body", rotor=True,
                    qty=1, bought=False, explode=[40, 0, 0],
                    notes="Print base (front edge) down, " + SUPPORT_NOTES["spinner"] + " "
                          + MATERIAL_PRINTED + ", 3 perimeters, 30 % infill in the solid tail."))
    cp = spinner_cap()
    out.append(dict(name="spinner_cap", shape=cp, print=_print_pose(cp), material="hub", role="body", rotor=True,
                    qty=1, bought=False, explode=[55, 0, 0],
                    notes="Print spigot down, " + SUPPORT_NOTES["spinner_cap"] + " " + MATERIAL_PRINTED
                          + ", 3 perimeters, 100 % infill (small part)."))
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
                        notes="N52 block 6 x 2 x 2, magnetised through the 2 mm along x; CA-glued in the "
                              "6.1 x 2.1 stop-face seat, face 0.2 mm below the stop face; attracting the "
                              "blade's heel magnet (pair gap 0.4)."))
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


def to_mesh_closed(shape: Shape, tol: float = 0.005, angular: float = 0.05):
    """Finest closed tessellation: OCC leaves cracks on some B-spline seams at very fine tolerances,
    so the tolerance is relaxed (x2 steps, up to 0.04) until the mesh is a closed manifold."""
    from OCP.BRepTools import BRepTools
    t = tol
    while True:
        BRepTools.Clean_s(shape.wrapped)
        m = to_mesh(shape, t, angular)
        if m.is_watertight or t >= 0.04:
            m.metadata["tol"] = t
            return m
        t *= 2


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
    return {"mesh_tol_mm": mesh.metadata.get("tol"), "min_mm": float(best[0][0]), "at": [round(float(v), 2) for v in best[0][1]],
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


def _occ_hits(shape, o, d) -> list:
    """Sorted line parameters of every crossing of the line o + t d with the shape's faces."""
    from OCP.gp import gp_Dir, gp_Lin, gp_Pnt
    from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
    it = IntCurvesFace_ShapeIntersector()
    it.Load(shape.wrapped, 1e-7)
    it.Perform(gp_Lin(gp_Pnt(*map(float, o)), gp_Dir(*map(float, d))), -1e3, 1e3)
    return sorted(it.WParameter(i) for i in range(1, it.NbPnt() + 1))


def seat_walls() -> dict:
    """Material around the deployed-magnet seat of blade 0 on the real hub BREP (T03): lines along
    e_r (+Z at psi 0) and along tau from inside the seat, and from the mouth rim (x STOP_X - 0.02)
    outward along e_r and diagonally; the wall is the next material run after the seat face."""
    h = hub()
    c = magnet_center0()
    x0 = STOP_X - MAG_RECESS - DEP_MAG[2]
    sr, st = DEP_SEAT[1] / 2, DEP_SEAT[0] / 2
    above, below, side, rim = [], [], [], []
    for x in np.linspace(x0 + 0.05, STOP_X - 0.02, 12):
        for t in (-st + 0.05, 0.0, st - 0.05):
            o = np.array([x, -t, c[2]])
            hp = [v for v in _occ_hits(h, o, (0, 0, 1)) if v > 0]
            if len(hp) >= 2:
                above.append(hp[1] - hp[0])
            hn = sorted(-v for v in _occ_hits(h, o, (0, 0, 1)) if v < 0)
            if len(hn) >= 2:
                below.append(hn[1] - hn[0])
        for sg in (1, -1):
            o = np.array([x, 0.0, c[2]])
            hs = [v for v in _occ_hits(h, o, (0, -sg, 0)) if v > 0]
            if len(hs) >= 2:
                side.append(hs[1] - hs[0])
    for t in np.linspace(-st, st, 9):
        for d in ((0, 0, 1), (0.0, 0.0, -1.0), (-0.7071, 0, 0.7071)):
            for zr in (c[2] + sr, c[2] - sr):
                if (d[2] > 0) != (zr > c[2]):
                    continue
                o = np.array([STOP_X - 0.02, -t, zr]) + 0.01 * np.array(d)
                hh = [v for v in _occ_hits(h, o, d) if v > 0]
                if hh:
                    rim.append(hh[0])
    return {"above_min_mm": float(min(above)), "below_min_mm": float(min(below)), "side_min_mm": float(min(side)),
            "mouth_rim_min_mm": float(min(rim)), "n_lines": len(above) + len(below) + len(side) + len(rim),
            "note": "rays from inside the seat along e_r (up: the knuckle arc of the slot, down: the heel slot) "
                    "and along the pin (sides), and from the seat mouth rim on the stop face outward"}


def layer_scan(stl: Path, layer: float = 0.2, min_area: float = 1.0) -> dict:
    """Layer-by-layer unsupported regions of a print STL (45 deg rule, 0.2 mm layers): area, inscribed
    width and the fraction of each region's boundary within 0.25 mm of the layer below (anchored).
    A bridge anchored on both sides reads >= ~0.4; a cantilever tread ~0.3 or less."""
    import manifold3d as mf
    import trimesh
    from matplotlib.path import Path as MP
    tm = trimesh.load(stl, force="mesh")
    M = mf.Manifold(mf.Mesh(vert_properties=np.asarray(tm.vertices, np.float32),
                            tri_verts=np.asarray(tm.faces, np.uint32)))
    zs = np.arange(0.1 + layer, float(tm.bounds[1][2]), layer)
    prev = M.slice(0.1)
    regs = []
    for z in zs:
        cur = M.slice(float(z))
        un = cur - prev.offset(layer + 1e-3, mf.JoinType.Round)
        for cp in un.decompose():
            a = cp.area()
            if a < min_area:
                continue
            lo, hi = 0.0, 30.0
            for _ in range(14):
                m_ = 0.5 * (lo + hi)
                if cp.offset(-m_, mf.JoinType.Round).is_empty():
                    hi = m_
                else:
                    lo = m_
            pts = []
            for poly in cp.to_polygons():
                q = np.vstack([poly, poly[:1]])
                for a_, b_ in zip(q[:-1], q[1:]):
                    n = max(1, int(np.linalg.norm(b_ - a_) / 0.1))
                    pts.append(a_ + np.outer((np.arange(n) + 0.5) / n, b_ - a_))
            pts = np.vstack(pts)
            w = np.zeros(len(pts))
            for poly in prev.offset(0.25, mf.JoinType.Round).to_polygons():
                poly = np.asarray(poly, float)
                if len(poly) < 3:
                    continue
                sgn = 1 if 0.5 * (np.dot(poly[:, 0], np.roll(poly[:, 1], -1)) - np.dot(poly[:, 1], np.roll(poly[:, 0], -1))) > 0 else -1
                w += sgn * MP(poly).contains_points(pts)
            regs.append({"z_mm": round(float(z), 2), "area_mm2": round(float(a), 2), "width_mm": round(2 * lo, 2),
                         "anchored": round(float((w > 0.5).mean()), 3)})
        prev = cur
    return {"n_regions": len(regs), "total_mm2": round(sum(r["area_mm2"] for r in regs), 1), "regions": regs[:20],
            "method": "manifold slices every %.1f mm; region = slice minus the previous slice grown by one layer "
                      "(45 deg); regions < %.0f mm2 ignored" % (layer, min_area)}


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
    c_m = magnet_center0()
    add("deployed magnet seat (block 6 x 2 x 2, glued): seat tau x r_e / depth / face recess / centre r_e (mm)",
        [list(DEP_SEAT), round(DEP_MAG[2] + MAG_RECESS, 3), MAG_RECESS, DEP_MAG_RE], "inside the stop face",
        b0["stop_face_r_e"][0] + 1.0 <= DEP_MAG_RE - DEP_SEAT[1] / 2 and DEP_MAG_RE + DEP_SEAT[1] / 2
        <= b0["stop_face_r_e"][1] - P.MIN_WALL,
        "PROPOSED DETENTS.deployed_hub_magnet_x_re_tau = [%.2f, %.2f, 0] (was the 4 x 2 disc at r_e 8.4388)"
        % (c_m[0], DEP_MAG_RE))
    eng = SCREW_L - AFT_WALL
    add("spinner screw M3 x 8: thread engagement in the shaft (mm)", eng, ">= 1.5 d = 4.5 and <= 8 (tap depth)",
        4.5 <= eng <= 8.0)
    add("access bore D vs M3 head D5.5 (radial clearance mm)", round(ACCESS_R - SCREW_HEAD_D / 2, 3), ">= 0.2",
        ACCESS_R - SCREW_HEAD_D / 2 >= 0.2)
    tg = tip_geometry()
    r_sp = CAP_SPIGOT_D / 2
    x_s0 = tg["x_countersink"] - (ACCESS_R - r_sp) - CAP_SPIGOT_L
    head_top = X_HEAD + SCREW_HEAD_H
    add("tip cap: spigot end to the M3 head (axial clearance, mm)", round(x_s0 - head_top, 3), ">= 1.0",
        x_s0 - head_top >= 1.0, "spigot D%.1f x %.1f from x %.2f; head top x %.1f" % (CAP_SPIGOT_D, CAP_SPIGOT_L,
                                                                                     x_s0, head_top))
    add("tip cap: joint x / tip x / built length / L over D", [round(tg["x_joint"], 3), round(tg["x_tip"], 3),
                                                              round(tg["x_tip"], 3), round(tg["x_tip"] / P.HULL_D, 3)],
        "HULL_L %g (Myring point), L/D %.2f" % (P.HULL_L, P.HULL_L / P.HULL_D), tg["x_tip"] > X_JOINT + 5.0,
        "the R%g round of the Myring point ends the egg at x %.2f" % (P.SPINNER_TIP_R, tg["x_tip"]))
    add("tip cap: fingernail groove depth below the hull (mm)", round(tg["groove_depth"], 3), ">= 0.2",
        tg["groove_depth"] >= 0.2)
    xx = np.linspace(tg["x_countersink"], X_JOINT - CAP_GROOVE, 200)
    wall_cs = float(np.min((hull_r(xx) - (ACCESS_R + xx - tg["x_countersink"])) / math.sqrt(2)))
    add("spinner wall between the 45 deg countersink and the hull (normal, mm, before the groove chamfer)",
        round(wall_cs, 3), ">= 0.3 at the groove edge", wall_cs >= 0.3)
    hub_m = to_manifold(meshes["hub"])
    spin_m = to_manifold(meshes["spinner"])
    add("D-flat start to the grub bore (wall, mm)", round(GRUB_X - GRUB_TAP_D / 2 - FLAT_X[0], 3), ">= MIN_WALL",
        GRUB_X - GRUB_TAP_D / 2 - FLAT_X[0] >= P.MIN_WALL - 1e-9,
        "and %.2f mm aft of the 0.5 mm bore chamfer" % (FLAT_X[0] - (X_FRONT + BORE_CHAMFER)))
    add("grub M3 thread (major D3) to the hub front face (wall, mm); front-face grooves removed (T11)",
        round(GRUB_X - 1.5 - X_FRONT, 3), ">= MIN_WALL", GRUB_X - 1.5 - X_FRONT >= P.MIN_WALL - 1e-9)
    sw = seat_walls()
    out["deployed_seat_walls"] = sw
    add("deployed magnet seat walls (OCC rays on the hub): to the knuckle arc above / below / to the slot sides, "
        "and from the mouth rim (mm)", [round(sw["above_min_mm"], 3), round(sw["below_min_mm"], 3),
                                       round(sw["side_min_mm"], 3), round(sw["mouth_rim_min_mm"], 3)],
        ">= MIN_WALL", min(sw["above_min_mm"], sw["below_min_mm"], sw["side_min_mm"],
                           sw["mouth_rim_min_mm"]) >= P.MIN_WALL - 0.02, sw["note"])
    W = window_layout()
    wp = window_parts0()
    out["windows"] = {"outline_offset_mm": WIN_OFFSET, "bridge_max_mm": WIN_BRIDGE,
                      "outer_gable": {k: (list(v) if isinstance(v, tuple) else v) for k, v in W["gable_S"].items()},
                      "deep_gable": {k: (list(v) if isinstance(v, tuple) else v) for k, v in W["gable_D"].items()},
                      "deep_band_parts": [list(d) for d in W["deep"]],
                      "skin_floor_exit_x_mm": wp["x_skin_exit"],
                      "skin_n_mm": WIN_SKIN_N, "full_depth_forward_of_x_mm": WIN_X_THIN}
    add("window full-depth part: flat top x / bridge width / ring of cone left under the seat ceiling (mm)",
        [round(W["gable_D"]["xa"], 3), round(W["gable_D"]["tau_a"][1] - W["gable_D"]["tau_a"][0], 3),
         round(X_AFT - W["gable_D"]["xa"], 3)], "bridge <= %g, ring >= 1.0" % WIN_BRIDGE,
        W["gable_D"]["tau_a"][1] - W["gable_D"]["tau_a"][0] <= WIN_BRIDGE + 1e-6 and X_AFT - W["gable_D"]["xa"] >= 1.0)
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
    out["wall_thickness"] = {k: wall_thickness(to_mesh_closed(p["shape"], 0.005, 0.05))
                             for k, p in ((q["name"], q) for q in part_list) if k in ("hub", "spinner", "spinner_cap")}
    print("  layer scan...", flush=True)
    out["layer_scan"] = {k: layer_scan(OUT_STL / f"{k}.stl") for k in ("hub", "spinner", "spinner_cap")}
    cap_m = to_manifold(meshes["spinner_cap"])
    out["cap_vs_spinner"] = {"overlap_mm3": float((cap_m ^ spin_m).volume()), "min_gap_mm": float(cap_m.min_gap(spin_m, 1.0)),
                             "note": "3 crush ribs (0.05 mm radial) pressed in the D6 bore; the 45 deg shoulder "
                                     "sits on the countersink (contact)"}
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
