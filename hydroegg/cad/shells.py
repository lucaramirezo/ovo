"""Static hull of HydroEgg: the front shell (nose to the equator, socket side) and the aft
shell (equator to the static tail end, spigot side), plus their cutaway variants.

    scripts/cadrun.sh .venv/bin/python -m hydroegg.cad.shells      # build, export, check, render

Frame and conventions as hydroegg/params.py: mm, X = body axis (flow +X, nose at x = 0), Z up.
Azimuth ``theta`` is measured about +X from +Z toward +Y (common.at_azimuth). The rotor study's
pocket table (data/layout_rotor.json pocket_depth) uses ``psi`` = -theta with its azimuths
relative to the blade (+ toward e_t); a pocket of the blade at theta_b therefore spans
theta_b + [-az_hi, -az_lo]. pocket_profile() converts once and asserts a known sample point.

Both shells print with +X pointing DOWN (front shell on its equator face, nose up; aft shell
on its tail end face, spigot tip up), so "down" for the 45 deg overhang rule is +X on both.

Design decisions taken here (the interface values themselves come from params.py and the
layout JSON; nothing below changes them):

* Pockets: ONE cutter per blade lofted (ruled) through the rows of the pocket_depth table:
  floor_r(x) over theta_b + [-az_hi, -az_lo]. The table floor already rises <= 45 deg (rotor
  study printable_floor). The side walls follow the table and are widened only where they
  would widen aft faster than 45 deg (an overhang when printed +X down) and stand >= 0.5 mm. A matching "liner" (the cutter grown by
  POCKET_FLOOR_WALL normal to the floor and the side walls) is added inside each shell first,
  so >= 2.0 mm of shell stays under every floor and beside every side wall.
* Socket / spigot: JOINT and BAYONET as given. The lip's forward face is a 45 deg cone that
  matches the lug's 45 deg aft face (0.2 mm axial play), the lip's aft face is a 45 deg chamfer
  (it is an overhang when the front shell prints), the twist channel's forward wall is a 45 deg
  ceiling. The channel ends exactly at the locked lug, so the hard stop aligns the pockets of
  both shells. The seal bore gets a 20 deg x 0.4 mm O-ring lead-in at the equator face.
* Carrier seat: the 1.6 mm ring (CARRIER.ring) follows the shell wall; the twist lock is
  defined in CARRIER_LOCK (shell side, the carrier plate builder matches it): 3 tabs enter
  through 3 notches in the ring's aft lip and turn 30 deg toward -theta (the generator's
  reaction torque drives them into the hard stop) inside a groove whose forward wall is a 45
  deg seat cone.
* Fins: pads and half-dovetail grooves per FIN (yawed with the fin). Retention: the reviewer is
  right that an M3 heat-set insert does not fit the 10 mm pad beside a 4 mm deep groove, so the
  fin is held by a 2 mm cross pin through the pad and the tongue (FIN_RAIL), hole printed as a
  teardrop; the groove is closed at the end stop, so drag pushes the fin home and the pin only
  keeps it from sliding forward in handling.
* Stowed detent: DETENTS gives a radial magnet axis at (x 104.25, az 18.75), but the stowed
  blade's inner face there leans ~49-51 deg from radial, so a radial 4 mm pair cannot meet
  face to face (the blade module puts its magnet on the face normal). The static magnet shares
  the blade magnet's mouth point on that face and leans 45 deg from radial toward the normal
  (the steepest support-free direction; boss face, hole and boss side all print at 45 deg),
  5.7 deg off the blade magnet's axis. The boss face clears the stowed blade by 0.1 mm over its
  whole 6.5 mm diameter (layout.py surface), which leaves ~1.0 mm between the magnet faces
  (~1.3 N per the layout force table, ~7x the blade's gravity moment). See STOW_DETENT.
* Bearing tube: crush ribs (3 per seat) reach D21.9 / D15.9 (PRESS_FIT per side) in relieved
  bores D22.4 / D16.3; the seal shoulder's underside is a 0.6 mm flat ring (the seal's stop)
  plus a 45 deg cone, the only intentional sub-45 deg-rule face of the tube (reported).
* Nose: plug forward of NOSE_PLUG.x_45 with a flat bridge (38.6 mm, reported as a bridge),
  tow eye ring with a teardrop hole, raised gland seat block, 10.5 mm tap-drill bore through
  block, wall and plug (tapped M12x1.5 after printing, no modelled thread).
"""

from __future__ import annotations

import functools
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from build123d import (
    Axis,
    Edge,
    Face,
    Location,
    Plane,
    Solid,
    Vector,
    Wire,
    export_brep,
    export_step,
    export_stl,
    import_brep,
)
from scipy.ndimage import maximum_filter1d

from .. import params as P
from .. import profile
from . import common as C

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
WORK = REPO / "out" / "work" / "cad" / "shells"
OUT_STEP = REPO / "out" / "cad" / "step"
OUT_STL = REPO / "out" / "cad" / "stl"
OUT_CHECKS = REPO / "out" / "cad" / "checks"

TAN45 = 1.0
DX_TAB = 0.5                               # pocket table row spacing (mm)
FLOOR_PRINT_DEG = 44.5                     # CAD pocket floor rise limit (table: 45, see pocket_profile)

# --------------------------------------------------------------------------
# Interface values (params.py; nothing here re-derives them)
# --------------------------------------------------------------------------
EQ = P.JOINT["x_eq"]                       # 57.5 equator plane
J = P.JOINT
BY = P.BAYONET
R_BORE = J["seal_bore_d"] / 2              # 32.045 socket seal bore
R_SPIGOT = J["spigot_seal_od"] / 2         # 31.945 spigot over the seal zone
R_GROOVE = J["groove_bottom_d"] / 2        # 29.94 O-ring groove bottom
R_LUG = J["lug_od"] / 2                    # 31.795 lug tip
R_SPIG_LUG = J["spigot_lug_zone_od"] / 2   # 30.295 spigot OD under the lugs
R_LIP = J["lip_d"] / 2                     # 30.545 socket lip
R_CH = J["channel_d"] / 2                  # 32.045 twist channel
R_SPIG_IN = J["spigot_bore_d"] / 2         # 27.94 spigot bore
X_TIP = J["x_spigot_tip"]                  # 43.595
X_LIP0, X_SEAL0 = J["x_lip"]               # 46.595, 49.595
X_CH0 = J["x_channel"][0]                  # 43.195
X_RING0 = J["x_ring_ceiling_start"]        # 41.695
GR = J["groove_x"]                         # x_aft 56, x_fwd_bottom 53.1, x_fwd_top 51.095

LUG_W, LUG_H, LUG_T = BY["lug_w_deg"], BY["lug_h"], BY["lug_t"]
LOCK = list(BY["lock_theta_deg"])          # [72, 192, 312]
ENTRY = list(BY["entry_theta_deg"])        # [37, 157, 277]
GAP_W = BY["gap_w_deg"]                    # 20
BAY_AX_CLEAR = 0.2                         # lug aft cone <-> lip cone, axial play
LUG_CHAMFER = 0.4                          # lug forward-outer edge (O-ring and entry lead)
SEAL_LEAD = (0.4, 20.0)                    # seal bore lead-in at the equator: radial mm, deg
EQ_CHAMFER = 0.5                           # outer edge of the equator faces (parting V, elephant foot)
TAIL_CHAMFER = 0.4                         # outer edge of the static tail end face

BLADE_THETA = [(-120.0 * i) % 360.0 for i in range(P.N_BLADES)]   # theta = -psi: 0, 240, 120

# Carrier twist lock, shell side. The carrier plate (drivetrain builder) carries 3 tabs that
# match this: x TAB_X, tip radius TAB_R, 45 deg chamfer TAB_CH on the front-outer corner,
# TAB_W wide, locked at lock_theta. Entry = lock + twist (tabs turn toward -theta to lock).
CARRIER_LOCK = dict(
    n_tabs=3, tab_w_deg=16.0, tab_x=[20.5, 22.5], tab_r_out=27.2, tab_front_chamfer_45=1.0,
    plate_r=P.CARRIER["r"], plate_x=[P.CARRIER["x0"], P.CARRIER["x0"] + P.CARRIER["t"]],
    lock_theta_deg=[90.0, 210.0, 330.0], twist_deg=P.CARRIER["twist_deg"],
    clear_radial=P.CLEAR_SLIDE, clear_axial=0.2, clear_cone_normal=0.2, notch_side_clear_deg=1.0,
    ring_x=list(P.CARRIER["ring"]["x"]), ring_t=P.CARRIER["ring"]["t_radial"],
    note="tabs pass the ring's aft lip through 3 notches at entry_theta, then turn twist_deg toward "
         "-theta (the generator reaction torque on the stator acts toward -theta and drives the tabs "
         "into the groove's end wall); forward seat = the tab's 45 deg chamfer on the groove's 45 deg "
         "ceiling, aft retention = the ring's lip (tab overlap ~1.1 mm). Plate front face nominal x "
         "18.5 with +-0.2 play.",
)
CARRIER_LOCK["entry_theta_deg"] = [(t + CARRIER_LOCK["twist_deg"]) % 360.0 for t in CARRIER_LOCK["lock_theta_deg"]]

# Fin rail (shell side): groove section in the yawed pad, cross pin retention.
FIN = P.FIN
FIN_RAIL = dict(
    pad_w=FIN["rail_w"], groove_depth=FIN["groove_depth"],
    groove_w_top=4.5, groove_neck_h=2.5, groove_w_bottom=6.0,
    groove_w_offset=-2.25,              # w of the flat (vertical) wall, on the fin's flat side (-w)
    end_stop_u=FIN["chord"],            # groove closed at u = chord (fin TE at x_te 55.5)
    pin_d=2.0, pin_hole_d=2.0, pin_l=10.0, pin_u=32.0, pin_depth=2.0,
    tongue_pin_hole_d=2.3,
    note="u: along the yawed root chord from the LE (x_le 16.37 at w = 0), w: across, + toward +theta; "
         "depth: radial at the fin angle below the root line (pad top). Groove: flat wall at w = -2.25 "
         "full depth (the fin's flat Clark Y side faces -theta), +w wall vertical to depth 2.5 then a 45 "
         "deg undercut to w = +3.75 at depth 4 (mean width 4.78). Cross pin: ISO 8734 2 x 10 (or a 1.75 "
         "mm filament stub) along w at u 32, depth 2.0, snug in the pad (D2.0 teardrop), D2.3 in the "
         "tongue.",
)

STOW_DETENT = dict(axis_lean_deg=45.0, boss_r=3.25, magnet_hole_d=P.MAGNET["D"] + 0.1,
                   magnet_hole_depth=P.MAGNET["H"] + 0.2, face_clear=0.1)

CUTAWAY_THETA = (-73.0, 17.0)              # removed sector (theta), see cutaway_sector()
# Seam azimuths of the full solids of revolution (windows free of every cutter that touches them):
SEAM_OUTER = 50.0      # hulls: between pocket 0 (theta <= ~26) and the theta 93 fin pad (>= ~79)
SEAM_CAVITY = 137.0    # front cavity: between carrier groove 1 (<= 129) and bayonet channel 2 (>= 146.5)
SEAM_TUBE = 30.0       # tube and bore: between the theta 0 rib and the 60 deg crush rib
SEAM_SPIGOT = 137.0    # spigot: between the lugs (64..80, 184..200, 304..320)


# --------------------------------------------------------------------------
# Small geometry helpers
# --------------------------------------------------------------------------

def r_out(x) -> float:
    return float(profile.radius(x))


def r_in(x) -> float:
    return float(profile.inner_radius(x))


def pt(x: float, r: float, th: float) -> Vector:
    """World point at axial x, radius r, azimuth theta (deg, from +Z toward +Y)."""
    t = math.radians(th)
    return Vector(x, r * math.sin(t), r * math.cos(t))


def _solve(f, a: float, b: float, n: int = 60) -> float:
    """Root of f on [a, b] by bisection (f(a), f(b) of opposite sign)."""
    fa = f(a)
    for _ in range(n):
        m = 0.5 * (a + b)
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def _profile_edges(nodes: list) -> list:
    """nodes: list of ('L', (x, r)) and ('S', [(x, r), ...]) items forming a closed chain in the
    meridian half plane (y = 0, z = r). 'S' items are splines whose first point must equal the
    previous node's last point."""
    edges, last = [], None
    first = None
    for kind, val in nodes:
        if kind == "L":
            p = Vector(val[0], 0.0, val[1])
            if last is not None and (p - last).length > 1e-9:
                edges.append(Edge.make_line(last, p))
            last = p
            first = first or p
        else:
            pts = [Vector(x, 0.0, r) for x, r in val]
            if last is not None and (pts[0] - last).length > 1e-9:
                edges.append(Edge.make_line(last, pts[0]))
            edges.append(Edge.make_spline(pts))
            last = pts[-1]
            first = first or pts[0]
    if (first - last).length > 1e-9:
        edges.append(Edge.make_line(last, first))
    return edges


def revolve_nodes(nodes: list, th0: float | None = None, th1: float | None = None,
                  seam: float = 50.0) -> Solid:
    """Revolve a closed meridian profile about X: full turn (seam at azimuth `seam`), or the
    sector [th0, th1] (deg). Seams must sit where no later boolean cuts across them: a seam
    inside the theta 0 rib made OCC return a null fuse, and a channel cut across a seam left a
    solid that passed BRepCheck but failed every later boolean through that sector."""
    face = Face(Wire(_profile_edges(nodes)))
    if th0 is None:
        return C.at_azimuth(Solid.revolve(face, 360.0, Axis.X), seam)
    s = Solid.revolve(face, th1 - th0, Axis.X)       # covers [-(th1 - th0), 0]
    return C.at_azimuth(s, th1)


def poly_sector(pts_xr, th0: float, th1: float) -> Solid:
    return revolve_nodes([("L", p) for p in pts_xr], th0, th1)


def curve(fn, x0: float, x1: float, n: int = 60, cos: bool = False) -> list:
    s = (1 - np.cos(np.linspace(0, math.pi, n))) / 2 if cos else np.linspace(0, 1, n)
    xs = x0 + (x1 - x0) * s
    return [(float(x), float(fn(x))) for x in xs]


def halfspace(point, normal, size: float = 600.0) -> Solid:
    """Big box whose face lies on the plane (point, normal), on the -normal side."""
    pl = Plane(origin=Vector(*point), z_dir=Vector(*normal))
    return pl * Solid.make_box(size, size, size).moved(Location((-size / 2, -size / 2, -size)))


def fuse(a, *others):
    out = a
    for o in others:
        out = out.fuse(o)
    return out


def inter(a, b) -> Solid:
    """Boolean common of two solids as one Solid (largest piece)."""
    r = a.intersect(b)
    sols = [s for item in (r or []) for s in item.solids()]
    if not sols:
        raise RuntimeError("empty intersection")
    return max(sols, key=lambda s: s.volume) if len(sols) == 1 else sorted(sols, key=lambda s: -s.volume)[0]


def teardrop(center, axis, up, r: float, L: float) -> Solid:
    """Horizontal-hole cutter: cylinder of radius r and length L centred on `center` along
    `axis`, plus the 45 deg point toward `up` (both unit vectors, up normal to axis)."""
    c, ax, up = Vector(*center), Vector(*axis).normalized(), Vector(*up).normalized()
    up = (up - ax * up.dot(ax)).normalized()
    right = up.cross(ax)
    o = c - ax * (L / 2)
    cyl = Solid.make_cylinder(r, L, Plane(origin=o, z_dir=ax))
    sq = Solid.make_box(r, r, L, Plane(origin=o, x_dir=(up + right).normalized(), z_dir=ax))
    return cyl.fuse(sq)


def single(shape, name: str) -> Solid:
    """Clean a boolean result and return its largest solid (asserting nothing else of size)."""
    shape = shape.clean()
    sols = sorted(shape.solids(), key=lambda s: -s.volume)
    if len(sols) > 1:
        extra = sum(s.volume for s in sols[1:])
        if extra > 1.0:
            raise RuntimeError(f"{name}: {len(sols)} solids, {extra:.2f} mm3 outside the main one")
    return sols[0]


# --------------------------------------------------------------------------
# Pockets (shared by both shells)
# --------------------------------------------------------------------------

def _printable_bounds(lo, hi, r_arc, dx, active=None):
    """Widen the angular bounds so no side wall widens aft (down in print) faster than 45 deg:
    going aft, a bound may move outward by at most dx (arc at radius r_arc) per row. Rows with
    active False are left as they are (the pocket tip, where the side walls are < 0.5 mm tall)."""
    lo, hi = lo.copy(), hi.copy()
    for i in range(len(lo) - 2, -1, -1):
        if active is not None and not active[i]:
            continue
        d = math.degrees(dx * TAN45 / r_arc[i])
        hi[i] = max(hi[i], hi[i + 1] - d)
        lo[i] = min(lo[i], lo[i + 1] + d)
    return lo, hi


@functools.lru_cache(maxsize=1)
def pocket_profile() -> dict:
    """Pocket and liner sections for the blade at theta 0, from layout_rotor.json pocket_depth.

    Returns arrays over x: floor (radius), th_lo / th_hi (theta, deg), and the liner's floor and
    bounds. Asserts the table's azimuth convention on the DETENTS stowed sample point."""
    tab = json.loads((DATA / "layout_rotor.json").read_text())
    pd = tab["pocket_depth"]
    x = np.array(pd["x_mm"], float)
    fl = np.array(pd["floor_r_mm"], float)
    az_lo = np.array(pd["az_lo_deg"], float)
    az_hi = np.array(pd["az_hi_deg"], float)
    assert np.allclose(np.diff(x), DX_TAB), "pocket table rows are not on the 0.5 mm grid"
    assert "theta = -psi" in P.POCKETS["tables"]
    assert "+ toward e_t" in tab["conventions"]["pocket_azimuth"]

    # --- convention check on a known sample: the stowed-detent point (x, az, r) of DETENTS
    from .. import layout as RL
    det = tab["detents"]["stowed"]
    xs_, az_, rs_ = det["x_mm"], det["az_deg"], det["blade_inner_r_mm"]
    _, e_r, e_t = RL.rotor_axes(0.0 + az_)                      # rotor frame at psi = az
    p_rot = np.array([xs_, 0.0, 0.0]) + rs_ * e_r
    p_cad = pt(xs_, rs_, -az_)                                 # theta = -psi
    assert np.allclose(p_rot, [p_cad.X, p_cad.Y, p_cad.Z], atol=1e-9), (p_rot, p_cad)
    assert RL.rotor_axes(0.0)[2][1] < 0                        # +az (e_t at psi 0) is -Y: theta < 0
    i = int(np.argmin(np.abs(x - xs_)))
    assert az_lo[i] <= az_ <= az_hi[i] and fl[i] < rs_         # the point lies inside the pocket

    # --- theta bounds of the pocket at theta_b = 0: [-az_hi, -az_lo], exactly the table (its 0.5
    # deg raster steps stay: smoothing them into an envelope cost the fin pads 0.4 mm of land),
    # widened only where a side wall would widen aft faster than 45 deg and is >= 0.5 mm tall
    th_lo, th_hi = -az_hi, -az_lo
    r_arc = np.array([r_out(v) for v in x]) + 1.0
    deep = (np.array([r_out(v) for v in x]) - fl) >= 0.5
    th_lo, th_hi = _printable_bounds(th_lo, th_hi, r_arc, DX_TAB, active=deep)
    # square tip: the rows ahead of the first 0.5 mm deep row (the table's narrow 0.01 mm deep tip
    # row) take that row's bounds, so the floor's 45 deg ramp closes the pocket with parallel side
    # walls instead of sub-millimetre ledges
    i0 = int(np.argmax(deep))
    th_lo[:i0], th_hi[:i0] = th_lo[i0], th_hi[i0]

    # --- floor rule check (the table must already rise <= 45 deg going forward), then the same
    # rule at FLOOR_PRINT_DEG: the ruled loft twists the table's exact 45 deg ramps to 45.34 deg,
    # so the CAD floor is cut a hair deeper (<= ~0.16 mm at the top of the root-notch ramp)
    rise = fl[:-1] - fl[1:]
    assert rise.max() <= DX_TAB * TAN45 + 1e-6, f"table floor rises {rise.max():.3f} per row"
    fl_table = fl.copy()
    kf = DX_TAB * math.tan(math.radians(FLOOR_PRINT_DEG))
    fl = fl.copy()
    for k in range(len(fl) - 2, -1, -1):
        fl[k] = min(fl[k], fl[k + 1] + kf)

    # --- extend half a row past both ends so the cutter leaves the hull cleanly
    x_e = np.concatenate([[x[0] - DX_TAB], x, [x[-1] + DX_TAB]])
    fl_e = np.concatenate([[fl[0] + kf], fl, [fl[-1]]])
    lo_e = np.concatenate([[th_lo[0] + 0.0], th_lo, [th_lo[-1]]])
    hi_e = np.concatenate([[th_hi[0]], th_hi, [th_hi[-1]]])
    assert fl_e[0] > r_out(x_e[0]), "forward pocket row must start outside the hull"

    # --- liner: floor - W along the floor normal, bounds widened by W at the floor radius
    W = P.POCKET_FLOOR_WALL + 0.15
    s_floor = np.abs(np.gradient(fl_e, x_e))
    s_hull = np.abs(np.array([profile.slope(v) for v in x_e]))
    s_max = maximum_filter1d(np.maximum(s_floor, s_hull), 3, mode="nearest")
    lf = fl_e - W * np.sqrt(1 + s_max ** 2)
    # liner inner surface: material may not step inward going up (-x) faster than 45 deg
    for k in range(1, len(lf)):
        lf[k] = min(lf[k], lf[k - 1] + DX_TAB * TAN45)
    dth = np.degrees(W / fl_e)
    l_lo, l_hi = lo_e - dth, hi_e + dth
    l_lo, l_hi = _printable_bounds(l_lo, l_hi, np.array([r_in(v) for v in x_e]), DX_TAB)
    return dict(x=x_e, floor=fl_e, th_lo=lo_e, th_hi=hi_e, liner_floor=lf, liner_lo=l_lo, liner_hi=l_hi,
                table=dict(x=x, floor=fl_table, az_lo=az_lo, az_hi=az_hi), cad_floor=fl, tip_rows=i0,
                sample=dict(x=xs_, az=az_, r=rs_))


def _sector_wire(x: float, r0: float, r1: float, th0: float, th1: float) -> Wire:
    thm = 0.5 * (th0 + th1)
    return Wire([
        Edge.make_line(pt(x, r0, th0), pt(x, r1, th0)),
        Edge.make_three_point_arc(pt(x, r1, th0), pt(x, r1, thm), pt(x, r1, th1)),
        Edge.make_line(pt(x, r1, th1), pt(x, r0, th1)),
        Edge.make_three_point_arc(pt(x, r0, th1), pt(x, r0, thm), pt(x, r0, th0)),
    ])


def _decimate(x, cols, tols) -> np.ndarray:
    """Indices of the rows to keep so that linear interpolation between kept rows reproduces
    every dropped row within tols (per column)."""
    keep = [0]
    i = 0
    n = len(x)
    while i < n - 1:
        j = i + 1
        while j + 1 < n:
            ok = True
            for c, tol in zip(cols, tols):
                t = (x[i + 1:j + 1] - x[i]) / (x[j + 1] - x[i])
                lin = c[i] + t * (c[j + 1] - c[i])
                if np.abs(lin - c[i + 1:j + 1]).max() > tol:
                    ok = False
                    break
            if not ok:
                break
            j += 1
        keep.append(j)
        i = j
    return np.array(keep)


def pocket_cutter(theta_b: float = 0.0, liner: bool = False, r_top: float = 41.0) -> Solid:
    """The reusable pocket cutter (or its liner) of the blade at theta_b: a ruled loft through
    the table rows (rows that linear interpolation reproduces within 0.01 mm / 0.02 deg are
    skipped)."""
    pp = pocket_profile()
    if liner:
        f, lo, hi = pp["liner_floor"], pp["liner_lo"], pp["liner_hi"]
    else:
        f, lo, hi = pp["floor"], pp["th_lo"], pp["th_hi"]
    k = _decimate(pp["x"], (f, lo, hi), (0.01, 0.02, 0.02))
    ws = [_sector_wire(float(pp["x"][i]), float(f[i]), r_top, float(lo[i]), float(hi[i])) for i in k]
    s = Solid.make_loft(ws, ruled=True)
    assert s.is_valid
    return C.at_azimuth(s, theta_b) if theta_b else s


def pocket_cutters(liner: bool = False) -> list:
    base = pocket_cutter(0.0, liner)
    return [base if t == 0 else C.at_azimuth(base, t) for t in BLADE_THETA]


# --------------------------------------------------------------------------
# Front shell
# --------------------------------------------------------------------------

def _front_hull() -> Solid:
    """Hull of revolution x 0 .. equator, with the 45 deg parting chamfer at the equator."""
    xe = EQ - EQ_CHAMFER
    return revolve_nodes([("S", curve(r_out, 0.0, xe, 140, cos=True)),
                          ("L", (EQ, P.HULL_D / 2 - EQ_CHAMFER)), ("L", (EQ, 0.0)), ("L", (0.0, 0.0))])


def _front_cavity() -> Solid:
    """Inside of the front shell: nose-plug bridge face, wall, carrier ring, socket ring, seal
    bore with its O-ring lead-in."""
    x45 = P.NOSE_PLUG["x_45"]
    rx0, rx1 = P.CARRIER["ring"]["x"]
    t_ring = P.CARRIER["ring"]["t_radial"]
    r_ri1 = r_in(rx1) - t_ring
    x_ch1 = _solve(lambda x: r_ri1 + (x - rx1) - r_in(x), rx1, rx1 + 10)        # ring aft 45 deg chamfer
    x_sr = _solve(lambda x: R_LIP + (X_RING0 - x) - r_in(x), X_RING0 - 8, X_RING0)  # socket ring ceiling
    lead_r, lead_deg = SEAL_LEAD
    x_lead = EQ - lead_r / math.tan(math.radians(lead_deg))
    nodes = [
        ("L", (x45, 0.0)), ("L", (x45, r_in(x45))),
        ("S", curve(r_in, x45, rx0, 30)),
        ("L", (rx0, r_in(rx0) - t_ring)),
        ("S", curve(lambda x: r_in(x) - t_ring, rx0, rx1, 12)),
        ("L", (x_ch1, r_in(x_ch1))),
        ("S", curve(r_in, x_ch1, x_sr, 24)),
        ("L", (X_RING0, R_LIP)),
        ("L", (X_SEAL0 - (R_BORE - R_LIP), R_LIP)),        # lip aft face: 45 deg chamfer
        ("L", (X_SEAL0, R_BORE)),
        ("L", (x_lead, R_BORE)), ("L", (EQ, R_BORE + lead_r)), ("L", (EQ + 1.0, R_BORE + lead_r)),
        ("L", (EQ + 1.0, 0.0)),
    ]
    return revolve_nodes(nodes, seam=SEAM_CAVITY)


def _socket_cutters() -> list:
    """Twist channels (entry gap edge to the locked lug, 45 deg ceiling, 45 deg lip cone) and the
    entry gaps through the lip, one per lug."""
    out = []
    # lip forward cone: parallel to the lug's 45 deg aft face, BAY_AX_CLEAR behind it
    x_cone = lambda r: X_TIP + LUG_T + BAY_AX_CLEAR - (r - R_SPIG_LUG)       # noqa: E731
    ch = [(X_RING0, R_LIP - 0.6), (X_RING0, R_LIP), (X_CH0, R_CH), (x_cone(R_CH), R_CH),
          (x_cone(R_LIP), R_LIP), (x_cone(R_LIP), R_LIP - 0.6)]
    # entry gap through the lip: its outer face (r = R_CH = the seal bore radius) stops exactly at
    # the bore edge and the cutter then drops 45 deg into the open bore, so it never shares a face
    # with the bore (a 0.4 mm overlap made OCC return an invalid solid for the 277 deg lug)
    gap = [(x_cone(R_CH) - 0.3, R_LIP - 0.6), (x_cone(R_CH) - 0.3, R_CH), (X_SEAL0, R_CH),
           (X_SEAL0 + 1.5, R_CH - 1.5), (X_SEAL0 + 1.5, R_LIP - 0.6)]
    for e, l in zip(ENTRY, LOCK):
        # one cutter per lug: channel (dead end 0.5 deg past the gap edge, so no face of the two
        # pieces is coplanar) fused with the entry gap
        out.append(poly_sector(ch, e - GAP_W / 2 - 0.5, l + LUG_W / 2).fuse(poly_sector(gap, e - GAP_W / 2, e + GAP_W / 2)))
    return out


def _carrier_cutters() -> list:
    """Twist-lock groove (with its 45 deg seat ceiling) and the entry notches in the ring."""
    L = CARRIER_LOCK
    t0, t1 = L["tab_x"]
    Rt, h = L["tab_r_out"], L["tab_front_chamfer_45"]
    Rg = Rt + L["clear_radial"]
    xa = t1 + L["clear_axial"]                         # groove aft wall
    b_tab = Rt - t0 - h                                # tab chamfer line r = x + b_tab
    b_g = b_tab + L["clear_cone_normal"] * math.sqrt(2.0)
    x_c = Rg - b_g                                     # ceiling meets the groove top
    r_low = 22.0
    groove = [(r_low - b_g, r_low), (x_c, Rg), (xa, Rg), (xa, r_low)]
    notch = [(t1 - 0.1, r_low), (t1 - 0.1, Rg), (t1 + 6.0, Rg), (t1 + 6.0, r_low)]
    hw = L["tab_w_deg"] / 2
    out = []
    for lk, en in zip(L["lock_theta_deg"], L["entry_theta_deg"]):
        lo = lk - hw                                   # hard stop wall (locked tab's -theta side)
        hi = en + hw + L["notch_side_clear_deg"]
        if hi < lo:
            hi += 360.0
        out.append(poly_sector(groove, lo, hi))
        out.append(poly_sector(notch, en - hw - L["notch_side_clear_deg"], en + hw + L["notch_side_clear_deg"]))
    return out


def _eye() -> Solid:
    """Tow eye: ring in the XZ plane (hole axis Y) ahead of the nose, teardrop hole (tip toward
    -X, up when the front shell prints nose-up), outer edges rounded."""
    E = P.TOW_EYE
    xc, rh = E["hole_centre_x"], E["hole_d"] / 2
    ro, t = rh + E["ring_w"], E["ring_t"]
    ring = Solid.make_cylinder(ro, t, Plane(origin=(xc, -t / 2, 0), z_dir=(0, 1, 0)))
    try:
        ring = ring.fillet(1.2, ring.edges())
    except Exception:  # noqa: BLE001
        pass
    # teardrop hole, point toward -X (up when the front shell prints nose-up)
    hole = teardrop((xc, 0.0, 0.0), (0, 1, 0), (-1, 0, 0), rh, t + 2.0)
    assert abs(hole.bounding_box().min.X - (xc - rh * math.sqrt(2.0))) < 1e-6
    return ring.cut(hole)


def gland_frame() -> dict:
    """GLAND pose: spot-face centre S, outward axis A (unit), from params."""
    G = P.GLAND
    a = math.radians(G["axis_angle_from_x_deg"])
    th = math.radians(G["theta_deg"])
    e_r = np.array([0.0, math.sin(th), math.cos(th)])
    A = math.cos(a) * np.array([1.0, 0.0, 0.0]) + math.sin(a) * e_r
    S = np.array([G["spot_face_centre_x"], 0.0, 0.0]) + G["spot_face_centre_r"] * e_r
    return dict(S=S, A=A, e_r=e_r)


GLAND_DRILL_D = 10.5                       # M12x1.5 tap drill
GLAND_BORE_LEN = 9.6                       # from the spot face: fully through the plug into the cavity


def _gland_block() -> Solid:
    G, F = P.GLAND, gland_frame()
    S, A = F["S"], F["A"]
    L = 14.0
    pl = Plane(origin=Vector(*(S - L * A)), z_dir=Vector(*A))
    blk = Solid.make_cylinder(G["seat_d"] / 2, L, pl)
    top = [e for e in blk.edges() if abs((e.position_at(0) - Vector(*S)).dot(Vector(*A))) < 1e-6
           and abs((e.position_at(0) - Vector(*S)).length - G["seat_d"] / 2) < 1e-3]
    try:
        blk = blk.chamfer(0.5, None, top)
    except Exception:  # noqa: BLE001
        pass
    return blk


def _gland_bore() -> Solid:
    F = gland_frame()
    S, A = F["S"], F["A"]
    r = GLAND_DRILL_D / 2
    pl = Plane(origin=Vector(*(S - GLAND_BORE_LEN * A)), z_dir=Vector(*A))
    bore = Solid.make_cylinder(r, GLAND_BORE_LEN + 3.0, pl)
    cs = 0.75                                     # 45 deg entry chamfer for the tap
    cone = Solid.make_cone(r, r + cs + 3.0, cs + 3.0, Plane(origin=Vector(*(S - cs * A)), z_dir=Vector(*A)))
    return bore.fuse(cone)


# ---- fins ---------------------------------------------------------------

def fin_frame() -> dict:
    """Yawed fin root frame at theta = 0 (lat = +Y toward +theta, along = +Z radial)."""
    i = math.radians(FIN["incidence_deg"])
    sense = -1.0 if FIN["te_side"] == "-theta" else 1.0
    ci, si = math.cos(i), sense * math.sin(i)
    c = FIN["chord"]
    x_te = FIN["x_te"]
    x_le = x_te - c * ci
    x_p = x_le + FIN["yaw_pivot_chord_frac"] * c * ci
    y_le, y_te = FIN["root_r_le"], FIN["root_r_te"]
    k = (y_te - y_le) / (x_te - x_le)
    u_dir = np.array([ci, si, 0.0])
    w_dir = np.array([-si, ci, 0.0])
    o_le = np.array([x_p, 0.0, 0.0]) - FIN["yaw_pivot_chord_frac"] * c * u_dir   # LE on the hull line
    assert abs(o_le[0] - FIN["x_le"]) < 1e-3, (o_le, FIN["x_le"])

    def root_z(x):
        return y_le + k * (x - x_le)

    def P_(u, w=0.0, depth=0.0):
        """Point on (or depth below) the root line at chord station u and lateral w."""
        q = o_le + u * u_dir + w * w_dir
        return np.array([q[0], q[1], root_z(q[0]) - depth])

    return dict(ci=ci, si=si, u=u_dir, w=w_dir, o=o_le, k=k, x_le=x_le, x_te=x_te, root_z=root_z, P=P_,
                phi=math.atan(k))


def _pad_local() -> Solid:
    """Rail pad at theta 0 (before at_azimuth): yawed 10 mm prism under the root line, from the
    LE end to the equator plane (its aft end stands on the bed), top edges rounded."""
    F = fin_frame()
    W = FIN["rail_w"]
    u_len = 48.0
    o = F["o"] - (W / 2) * F["w"] + np.array([0.0, 0.0, 18.0])
    box = Solid.make_box(u_len, W, 40.0, Plane(origin=Vector(*o), x_dir=Vector(*F["u"]), z_dir=(0, 0, 1)))
    phi = F["phi"]
    top = halfspace((F["x_le"], 0.0, F["root_z"](F["x_le"])), (-math.sin(phi), 0.0, math.cos(phi)))
    eq = halfspace((EQ, 0.0, 0.0), (1.0, 0.0, 0.0))
    pad = inter(inter(box, top), eq)
    # round the front vertical edges first, then the whole top boundary (except the bed edge at the
    # equator) in one rolling-ball pass: filleting the top edges first left B-spline corner blends
    # that tessellated with cracks
    ud = Vector(*F["u"])
    front = [e for e in pad.edges() if abs((e.position_at(0) - Vector(*F["o"])).dot(ud)) < 1e-3
             and abs((e.position_at(1) - Vector(*F["o"])).dot(ud)) < 1e-3
             and abs((e.position_at(1) - e.position_at(0)).normalized().Z) > 0.95]
    pad = pad.fillet(1.5, front)
    nt = Vector(-math.sin(phi), 0.0, math.cos(phi))
    topf = max(pad.faces(), key=lambda f: f.normal_at(f.center()).dot(nt))
    edges = [e for e in topf.edges() if min(e.position_at(0).X, e.position_at(1).X, e.center().X) < EQ - 0.2]
    pad = pad.fillet(1.0, edges)
    return pad


def _groove_local() -> Solid:
    """Half-dovetail groove (open at the pad's front, closed at the end stop) and the cross-pin
    teardrop hole, at theta 0."""
    F, R = fin_frame(), FIN_RAIL
    w0 = R["groove_w_offset"]
    wt, wb, hn, d = R["groove_w_top"], R["groove_w_bottom"], R["groove_neck_h"], R["groove_depth"]
    sec = [(w0, 1.5), (w0 + wt, 1.5), (w0 + wt, -hn), (w0 + wb, -d), (w0, -d)]
    u0, u1 = -1.0, R["end_stop_u"]
    base = F["P"](u0)
    pts = [Vector(*(base + w * F["w"] + np.array([0.0, 0.0, z]))) for w, z in sec]
    face = Face(Wire.make_polygon(pts, close=True))
    vec = F["P"](u1) - F["P"](u0)
    groove = Solid.extrude(face, Vector(*vec))
    # cross pin: along w through the whole pad, teardrop tip toward -X (up in print)
    c = F["P"](R["pin_u"], 0.0, R["pin_depth"])
    pin = teardrop(c, F["w"], (-1.0, 0.0, 0.0), R["pin_hole_d"] / 2, R["pad_w"] + 4.0)
    return groove, pin


def pads_added() -> list:
    """The pad material the front shell really gains (pads outside the hull's 1 mm skin core)."""
    core = _front_core()
    return [p.cut(core) for p in _pads()[0]]


def _pads() -> tuple[list, list]:
    pad0 = _pad_local()
    groove0, pin0 = _groove_local()
    pads = [C.at_azimuth(pad0, t) for t in FIN["theta_deg"]]
    cuts = [C.at_azimuth(g, t) for t in FIN["theta_deg"] for g in (groove0, pin0)]
    return pads, cuts


# --------------------------------------------------------------------------
# Cutaway sector
# --------------------------------------------------------------------------

def cutaway_sector(front: bool) -> Solid:
    """Removed sector CUTAWAY_THETA. Front shell: from a cone through the wall at x 16 (the nose
    keeps its plug, eye and gland; the cut end overhangs < 45 deg) to the equator. Aft shell: the
    whole length."""
    th0, th1 = CUTAWAY_THETA
    if front:
        # boundary x_b(r) = xw - 1.25 (r - rw): the cut face leans 38.7 deg from vertical when
        # printed nose-up (|dx/dr| >= 1 is the 45 deg rule) and crosses the wall almost square
        xw = 16.0
        rw = 0.5 * (r_out(xw) + r_in(xw))
        xb = lambda r: xw - 1.25 * (r - rw)         # noqa: E731
        pts = [(xb(0.0), 0.0), (xb(60.0), 60.0), (EQ + 2.0, 60.0), (EQ + 2.0, 0.0)]
    else:
        pts = [(X_TIP - 2.0, 0.0), (X_TIP - 2.0, 60.0), (P.STATIC_TAIL_X + 2.0, 60.0), (P.STATIC_TAIL_X + 2.0, 0.0)]
    return poly_sector(pts, th0, th1)


def _front_core(skin: float = 1.0) -> Solid:
    """Hull of revolution shrunk by `skin` (radially): external features are clipped to stay
    outside it, so they merge with the outer skin only and never meet the inner features."""
    x0 = _solve(lambda x: r_out(x) - skin - 0.5, 0.0, 5.0)
    return revolve_nodes([("L", (x0, 0.0)), ("S", curve(lambda x: r_out(x) - skin, x0, EQ + 1.0, 100, cos=True)),
                          ("L", (EQ + 1.0, 0.0))])


def build_front_shell() -> Solid:
    t0 = time.time()
    hull = _front_hull()
    core = _front_core()
    pads, pad_cuts = _pads()
    ext = [p.cut(core) for p in pads] + [_gland_block().cut(core)]
    outer = hull.fuse(_eye(), *ext)
    body = outer.cut(_front_cavity())
    liners = [inter(l, hull) for l in pocket_cutters(liner=True)]
    body = body.fuse(*liners)
    cutters = _socket_cutters() + _carrier_cutters() + [_gland_bore()] + pad_cuts + pocket_cutters()
    body = body.cut(*cutters)
    body = single(body, "front_shell")
    print(f"front shell {time.time() - t0:.1f} s", flush=True)
    return body


# --------------------------------------------------------------------------
# Aft shell
# --------------------------------------------------------------------------

def _aft_hull() -> Solid:
    """Hull of revolution from the equator to the static tail end (parting and tail chamfers)."""
    x0, x1 = EQ + EQ_CHAMFER, P.STATIC_TAIL_X - TAIL_CHAMFER
    return revolve_nodes([("L", (EQ, 0.0)), ("L", (EQ, P.HULL_D / 2 - EQ_CHAMFER)),
                          ("S", curve(r_out, x0, x1, 120)),
                          ("L", (P.STATIC_TAIL_X, r_out(x1) - TAIL_CHAMFER)), ("L", (P.STATIC_TAIL_X, 0.0))])


def _spigot() -> Solid:
    """Spigot ring: lug zone, step, seal zone with the O-ring groove (aft flank radial, tip-side
    flank 45 deg), tip chamfer."""
    ch = 0.5
    r0 = R_SPIG_IN - 0.5                       # the cavity cut trims it back to the spigot bore
    nodes = [("L", (X_TIP, r0)), ("L", (X_TIP, R_SPIG_LUG - ch)), ("L", (X_TIP + ch, R_SPIG_LUG)),
             ("L", (X_SEAL0, R_SPIG_LUG)), ("L", (X_SEAL0, R_SPIGOT)),
             ("L", (GR["x_fwd_top"], R_SPIGOT)), ("L", (GR["x_fwd_bottom"], R_GROOVE)),
             ("L", (GR["x_aft"], R_GROOVE)), ("L", (GR["x_aft"], R_SPIGOT)),
             ("L", (EQ + 0.01, R_SPIGOT)), ("L", (EQ + 0.01, r0))]
    return revolve_nodes(nodes, seam=SEAM_SPIGOT)


def _lugs() -> list:
    pts = [(X_TIP, R_SPIG_LUG - 0.3), (X_TIP, R_LUG - LUG_CHAMFER), (X_TIP + LUG_CHAMFER, R_LUG),
           (X_TIP + LUG_T - LUG_H, R_LUG), (X_TIP + LUG_T, R_SPIG_LUG), (X_TIP + LUG_T, R_SPIG_LUG - 0.3)]
    return [poly_sector(pts, l - LUG_W / 2, l + LUG_W / 2) for l in LOCK]


def _aft_cavity() -> Solid:
    """Spigot bore, 45 deg root ring, hull inner wall down to the tail end block (x 101)."""
    xb = P.TAIL_END_BLOCK["x0"]
    x_m = _solve(lambda x: R_SPIG_IN + (x - EQ) - r_in(x), EQ, EQ + 12)
    nodes = [("L", (X_TIP - 1.0, 0.0)), ("L", (X_TIP - 1.0, R_SPIG_IN)), ("L", (EQ, R_SPIG_IN)),
             ("L", (x_m, r_in(x_m))), ("S", curve(r_in, x_m, xb, 60)), ("L", (xb, 0.0))]
    return revolve_nodes(nodes)


BT = P.BEARING_TUBE
R608 = 11.2          # 608 seat bore between the crush ribs (D22.4)
R608_RIB = P.BEARING["D"] / 2 + P.PRESS_FIT     # 10.95 (D21.9)
R688 = 8.15          # 688 seat bore between the crush ribs (D16.3)
R688_RIB = P.BEARING_ALT["D"] / 2 + P.PRESS_FIT  # 7.95 (D15.9)
R_SHOULDER = BT["shoulder_front_id"] / 2         # 9.0
R_SEAL_SH = BT["seal_shoulder_id"] / 2           # 6.25
R_SEAL = BT["seal_bore_d"] / 2                   # 8.0
SEAL_STOP_FLAT = 0.6


def _tube() -> Solid:
    px, pr = BT["profile_x"], BT["profile_r"]
    ch = 0.5
    nodes = [("L", (px[0], 0.0)), ("L", (px[0], pr[0] - ch)), ("L", (px[0] + ch, pr[0]))]
    nodes += [("L", (x, r)) for x, r in zip(px[1:-1], pr[1:-1])]
    xe = P.TAIL_END_BLOCK["x0"] + 0.5
    nodes += [("L", (xe, pr[-1])), ("L", (xe, 0.0))]
    return revolve_nodes(nodes, seam=SEAM_TUBE)


def _bore() -> Solid:
    """Tube bore: mouth lead-in, 608 seat, shoulder, 688 seat, seal shoulder (flat stop ring + 45
    deg), lip-seal bore with its press-in lead."""
    x608 = P.BEARING["x0"]                        # 63
    x_sh0, x_sh1 = BT["shoulder_front_x"]         # 70, 96
    xs0, xs1 = BT["seal_shoulder_x"]              # 101, 102.5
    li = BT["mouth_lead_in"]
    xm = BT["x0"]
    xe = P.STATIC_TAIL_X
    cone_dx = (R_SEAL - SEAL_STOP_FLAT) - R_SEAL_SH
    nodes = [("L", (xm - 1.0, 0.0)), ("L", (xm - 1.0, R608 + li)), ("L", (xm, R608 + li)),
             ("L", (xm + li, R608)), ("L", (x_sh0, R608)), ("L", (x_sh0, R_SHOULDER)),
             ("L", (x_sh1, R_SHOULDER)), ("L", (x_sh1, R688)), ("L", (xs0, R688)),
             ("L", (xs0, R_SEAL_SH)), ("L", (xs1 - cone_dx, R_SEAL_SH)),
             ("L", (xs1, R_SEAL - SEAL_STOP_FLAT)), ("L", (xs1, R_SEAL)),
             ("L", (xe - 0.5, R_SEAL)), ("L", (xe, R_SEAL + 0.5)), ("L", (xe + 1.0, R_SEAL + 0.5)),
             ("L", (xe + 1.0, 0.0))]
    assert x608 >= xm + li
    return revolve_nodes(nodes, seam=SEAM_TUBE)


def _crush_ribs() -> list:
    """3 crush ribs per bearing seat (between the shell ribs), lead ramp at the entry end."""
    out = []
    for r_seat, r_tip, x0, x1 in ((R608, R608_RIB, P.BEARING["x0"], BT["shoulder_front_x"][0]),
                                  (R688, R688_RIB, BT["shoulder_front_x"][1] + 0.3, BT["seal_shoulder_x"][0])):
        prof = [(x0, r_seat + 0.3), (x0, r_seat), (x0 + 0.6, r_tip), (x1, r_tip), (x1, r_seat + 0.3)]
        f = Face(Wire.make_polygon([Vector(x, -0.6, r) for x, r in prof], close=True))
        rib = Solid.extrude(f, Vector(0, 1.2, 0))
        out += [C.at_azimuth(rib, t) for t in (60.0, 180.0, 300.0)]
    return out


def _ribs(hull: Solid) -> list:
    rb = BT["ribs"]
    t = rb["t"]
    x0, x1 = rb["x"][0], P.TAIL_END_BLOCK["x0"] + 0.5
    plate = Solid.make_box(x1 - x0, t, 40.0 - 9.0).moved(Location((x0, -t / 2, 9.0)))
    return [inter(C.at_azimuth(plate, th), hull) for th in rb["theta_deg"]]


# ---- stowed detent ------------------------------------------------------

_BLADE_PTS = None


def stowed_blade_points(i: int = 0, step: float = 0.3) -> np.ndarray:
    """World points (stowed, fold FOLD_DEG) of blade i from hydroegg.layout's surface model."""
    global _BLADE_PTS
    from .. import layout as RL
    if _BLADE_PTS is None:
        blade = RL.load_blade()
        lay = RL.chosen_layout()
        ax = RL.solve_axis(blade, lay, n=41, iters=3)
        S = RL.blade_surface(blade, lay, ax)
        _BLADE_PTS = dict(lay=lay, S=S)
    lay, S = _BLADE_PTS["lay"], _BLADE_PTS["S"]
    T = RL.blade_to_world(lay.fold_deg, lay.psi(i), lay)
    clouds = [RL.densify(S.aero, step).reshape(-1, 3), RL.densify(S.shank, step).reshape(-1, 3),
              S.knuckle.reshape(-1, 3), S.heel]
    return RL.apply(T, np.vstack(clouds))


def stow_detent_geometry() -> dict:
    """Static stowed-detent seat for blade 0 (theta 0); the other blades are copies about X.

    The blade's stowed magnet sits where the DETENTS radial line (x 104.25, az 18.75) meets the
    blade's inner face and runs along that face's normal, ~51 deg from radial (the blade
    module does the same). The static magnet shares the mouth point and leans toward that normal
    as far as printing allows: 45 deg from radial (the boss face, the hole and the boss side are
    then all at the 45 deg limit when printed +X down). The boss face stands face_clear below the
    stowed blade (layout.py surface) over the whole boss."""
    det = P.DETENTS
    x_c, az = det["stowed_x"], det["stowed_az_deg"]
    th = math.radians(-az)
    e_r = np.array([0.0, math.sin(th), math.cos(th)])
    e_t = np.array([0.0, math.cos(th), -math.sin(th)])
    ex = np.array([1.0, 0.0, 0.0])
    B = stowed_blade_points(0, 0.15)
    o = np.array([x_c, 0.0, 0.0])
    rel = B - o
    al = rel @ e_r
    near = np.linalg.norm(rel - np.outer(al, e_r), axis=1) < 0.5
    p = o + float(al[near].min()) * e_r                                   # first guess
    for _ in range(3):                                                    # local plane fit, line crossing
        loc = B[np.linalg.norm(B - p, axis=1) < 1.5]
        cen = loc.mean(axis=0)
        _, _, vt = np.linalg.svd(loc - cen)
        nb = vt[2] if vt[2] @ e_r < 0 else -vt[2]                         # blade face normal, into the pocket
        loc = loc[np.abs((loc - cen) @ nb) < 0.3]                         # inner sheet only
        cen = loc.mean(axis=0)
        _, _, vt = np.linalg.svd(loc - cen)
        nb = vt[2] if vt[2] @ e_r < 0 else -vt[2]
        p = o + float(((cen - o) @ nb) / (e_r @ nb)) * e_r                # mouth on the inner face
    v = -nb                                                               # static -> blade
    a, b, c = v @ ex, v @ e_r, v @ e_t
    lean_face = math.degrees(math.atan2(a, b))
    lean = min(lean_face, STOW_DETENT["axis_lean_deg"])
    m = math.hypot(a, b)
    n = m * (math.sin(math.radians(lean)) * ex + math.cos(math.radians(lean)) * e_r) + c * e_t
    n /= np.linalg.norm(n)
    rel = B - p
    s = rel @ n
    perp = np.linalg.norm(rel - np.outer(s, n), axis=1)
    near = perp <= STOW_DETENT["boss_r"] + 0.3
    s_face = float(s[near].min()) - STOW_DETENT["face_clear"]
    face = p + s_face * n
    mag_c = face - (0.2 + P.MAGNET["H"] / 2) * n
    return dict(contact=p, mouth=p, n=n, blade_face_normal=nb, blade_face_lean_deg=lean_face, lean_deg=lean,
                axis_mismatch_deg=float(math.degrees(math.acos(min(1.0, float(n @ -nb))))),
                s_face=s_face, face=face, magnet_centre=mag_c, blade_min_s=float(s[near].min()),
                n_blade_pts=int(near.sum()))


def _detent_bosses(pockets: list) -> tuple[list, list]:
    g = stow_detent_geometry()
    n, face = g["n"], g["face"]
    D = STOW_DETENT
    pl = Plane(origin=Vector(*(face - 6.0 * n)), z_dir=Vector(*n))
    boss = Solid.make_cylinder(D["boss_r"], 6.0, pl)
    # no chamfers on the boss rim or the magnet-hole mouth: around an axis that already leans
    # 45 deg, a 45 deg chamfer cone has a strip facing straight down when printed
    rh, dh = D["magnet_hole_d"] / 2, D["magnet_hole_depth"]
    hole = Solid.make_cylinder(rh, dh + 1.0, Plane(origin=Vector(*(face - dh * n)), z_dir=Vector(*n)))
    bosses, holes = [], []
    for t, pk in zip(BLADE_THETA, pockets):
        b = boss if t == 0 else C.at_azimuth(boss, t)
        bosses.append(inter(b, pk))
        holes.append(hole if t == 0 else C.at_azimuth(hole, t))
    return bosses, holes


def build_aft_shell() -> Solid:
    t0 = time.time()
    hull = _aft_hull()
    outer = fuse(hull, _spigot())
    body = outer.cut(_aft_cavity())
    liners = [inter(l, hull) for l in pocket_cutters(liner=True)]
    body = body.fuse(_tube(), *_ribs(hull), *_lugs(), *liners)
    pockets = pocket_cutters()
    body = body.cut(_bore(), *pockets)
    bosses, holes = _detent_bosses(pockets)
    body = fuse(body, *_crush_ribs(), *bosses)
    body = body.cut(*holes)
    body = single(body, "aft_shell")
    print(f"aft shell {time.time() - t0:.1f} s", flush=True)
    return body


# --------------------------------------------------------------------------
# Parts, cache, previews
# --------------------------------------------------------------------------

PRINT_ROT = (0.0, 90.0, 0.0)               # +X -> -Z: equator / tail end face on the bed
NAMES = ("front_shell", "aft_shell", "front_shell_cutaway", "aft_shell_cutaway")


def _cache_key() -> str:
    h = hashlib.sha1()
    for p in (Path(__file__), REPO / "hydroegg/params.py", REPO / "hydroegg/profile.py",
              REPO / "hydroegg/cad/common.py", DATA / "layout_rotor.json"):
        h.update(p.read_bytes())
    return h.hexdigest()[:12]


def build_all(force: bool = False) -> dict:
    """The four shell solids (world pose), cached as BREP under out/work/cad/shells/cache."""
    cdir = WORK / "cache" / _cache_key()
    files = {n: cdir / f"{n}.brep" for n in NAMES}
    if not force and all(f.exists() for f in files.values()):
        return {n: import_brep(str(f)) for n, f in files.items()}
    import shutil
    for old in (WORK / "cache").glob("*"):                 # stale keys (older sources)
        if old.is_dir() and old != cdir:
            shutil.rmtree(old, ignore_errors=True)
    front = build_front_shell()
    aft = build_aft_shell()
    fc = single(front.cut(cutaway_sector(True)), "front_shell_cutaway")
    ac = single(aft.cut(cutaway_sector(False)), "aft_shell_cutaway")
    out = dict(front_shell=front, aft_shell=aft, front_shell_cutaway=fc, aft_shell_cutaway=ac)
    cdir.mkdir(parents=True, exist_ok=True)
    for n, s in out.items():
        export_brep(s, str(files[n]))
    return out


NOTES = {
    "front_shell": "PLA (display) / PETG (v2). Prints nose-up on its equator face, no supports: the nose "
                   "plug's 38.6 mm flat bridge is the only bridge; eye and pin holes are teardrops; socket "
                   "ceiling, lip chamfer, carrier groove ceiling and pocket ramps are 45 deg; fin pads stand "
                   "on the bed. Tap the M12x1.5 after printing (10.5 mm drill modelled).",
    "aft_shell": "PLA / PETG. Prints on its tail end face, spigot tip up, no supports: lugs, O-ring groove "
                 "flank, root ring, notch ramps and detent bosses at 45 deg; the 0.6 mm seal stop ring under "
                 "the seal shoulder is the only sub-rule face. Bearings pressed with a draw bolt; crush ribs "
                 "set the press fit.",
}


def parts() -> list[dict]:
    S = build_all()
    out = []
    for n in NAMES:
        base = n.replace("_cutaway", "")
        cut = n.endswith("_cutaway")
        out.append(dict(
            name=n, shape=S[n], print=C.on_bed(S[n], PRINT_ROT), material="shell",
            role="shell_cutaway" if cut else "shell", rotor=False, qty=1, bought=False,
            explode=[-80.0, 0.0, 0.0] if base == "front_shell" else [30.0, 0.0, 0.0],
            notes=NOTES[base] + (" Cutaway: sector theta %g..%g removed (display only)." % CUTAWAY_THETA if cut else ""),
        ))
    return out


def _view(dicts):
    """Preview dicts (no 'print' entry, it is not JSON). SHELLS_PREVIEW_MAT overrides the material
    (a darker one shows edges better while modelling)."""
    import os
    mat = os.environ.get("SHELLS_PREVIEW_MAT")
    out = []
    for d in dicts:
        e = {k: v for k, v in d.items() if k in ("name", "shape", "material", "role")}
        if mat:
            e["material"] = mat
        out.append(e)
    return out


def preview_parts() -> list[dict]:
    """Both shells, assembled pose, outside."""
    return _view([d for d in parts() if d["role"] == "shell"])


def preview_cutaway() -> list[dict]:
    return _view([dict(d, role="body") for d in parts() if d["role"] == "shell_cutaway"])


def preview_front() -> list[dict]:
    return _view([d for d in parts() if d["name"] == "front_shell"])


def preview_aft() -> list[dict]:
    return _view([d for d in parts() if d["name"] == "aft_shell"])


def preview_section() -> list[dict]:
    """Shells cut by the XZ plane, the -Y half removed (view=side looks from -Y)."""
    keep = halfspace((0.0, 0.0, 0.0), (0.0, -1.0, 0.0))
    out = []
    for d in parts():
        if d["role"] == "shell":
            out.append(dict(name=d["name"] + "_sec", shape=inter(d["shape"], keep), material="shell"))
    return _view(out)


def preview_print() -> list[dict]:
    """Both shells in print orientation, side by side."""
    out = []
    for k, d in enumerate([d for d in parts() if d["role"] == "shell"]):
        out.append(dict(name=d["name"] + "_print", shape=d["print"].moved(Location((k * 110.0 - 55.0, 0, 0))),
                        material="shell"))
    return _view(out)


# --------------------------------------------------------------------------
# Checks (measured on the solids), exports, figures
# --------------------------------------------------------------------------

def _occ_rays(shape, origins, dirs, far: float = 200.0) -> np.ndarray:
    """Distance from each origin along its unit direction to the first face of `shape` (inf if
    none). Origins inside the material give the material thickness along the ray."""
    from OCP.gp import gp_Dir, gp_Lin, gp_Pnt
    from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
    isec = IntCurvesFace_ShapeIntersector()
    isec.Load(shape.wrapped, 1e-7)
    out = np.full(len(origins), np.inf)
    for k, (o, d) in enumerate(zip(origins, dirs)):
        isec.Perform(gp_Lin(gp_Pnt(*map(float, o)), gp_Dir(*map(float, d))), 1e-6, far)
        if isec.IsDone() and isec.NbPnt() > 0:
            out[k] = min(isec.WParameter(i) for i in range(1, isec.NbPnt() + 1))
    return out


def _occ_inside(shape, pts) -> np.ndarray:
    """Point-in-solid by a majority vote of ray parities (3 fixed skew directions). OCC's
    BRepClass3d classifier (and distance_to, which uses it) reported points 2.5 mm inside an
    empty pocket of the aft shell as IN; ray parity in 7 directions said OUT."""
    from OCP.gp import gp_Dir, gp_Lin, gp_Pnt
    from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
    isec = IntCurvesFace_ShapeIntersector()
    isec.Load(shape.wrapped, 1e-7)
    dirs = [(0.5774, 0.6513, 0.4924), (-0.3122, 0.8331, -0.4566), (0.7071, -0.4083, 0.5774)]
    out = np.zeros(len(pts), bool)
    for k, p in enumerate(pts):
        votes = 0
        for d in dirs:
            isec.Perform(gp_Lin(gp_Pnt(*map(float, p)), gp_Dir(*d)), 1e-6, 1e3)
            votes += isec.NbPnt() % 2
        out[k] = votes >= 2
    return out


def _mesh(shape, tol: float = 0.03):
    import trimesh
    v, f = shape.tessellate(tol, 0.1)
    return trimesh.Trimesh(np.array([(q.X, q.Y, q.Z) for q in v]), np.array(f), process=True)


def overhang_scan(shape, x_bed: float, limit: float = P.MAX_OVERHANG, tol_deg: float = 0.5) -> dict:
    """Faces that point down (+X: both shells print +X down) more than `limit` from vertical and
    are not on the bed (x = x_bed). Angles use the exact B-rep normal at every tessellation
    triangle's centroid (a tessellated 45 deg cone would otherwise read 45.7 deg). Faces within
    tol_deg of the limit (the 45 deg ramps of ruled lofts read up to 45.34) are counted as "at the
    limit" with their largest angle. Reports the offending B-rep faces with area, angle, position."""
    from build123d import GeomType
    bad_faces, area_bad, area_45, area_bed, max_at_limit = [], 0.0, 0.0, 0.0, 0.0
    for f in shape.faces():
        v, t = f.tessellate(0.03, 0.1)
        if not t:
            continue
        V = np.array([(q.X, q.Y, q.Z) for q in v])
        T = np.array(t)
        cr = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])
        ar = 0.5 * np.linalg.norm(cr, axis=1)
        ok = ar > 1e-12
        if not ok.any():
            continue
        T, cr, ar = T[ok], cr[ok], ar[ok]
        nt = cr / (2 * ar[:, None])
        cen = V[T].mean(axis=1)
        if f.geom_type == GeomType.PLANE:
            n0 = f.normal_at(f.center())
            ne = np.tile(np.array([n0.X, n0.Y, n0.Z]), (len(T), 1))
        else:
            ne = nt.copy()
            for k_, c_ in enumerate(cen):
                try:
                    q_ = f.normal_at(Vector(*c_))
                    ne[k_] = (q_.X, q_.Y, q_.Z)
                except Exception:  # noqa: BLE001 (projection did not converge: keep the facet normal)
                    pass
        ne *= np.sign(np.sum(ne * nt, axis=1))[:, None]
        ang = np.degrees(np.arcsin(np.clip(ne[:, 0], -1.0, 1.0)))       # + = facing down (+X)
        is_bed = f.geom_type == GeomType.PLANE and abs(cen[:, 0].mean() - x_bed) < 0.02 and ne[0, 0] > 0.999
        if is_bed:
            area_bed += float(ar.sum())
            continue
        b = ang > limit + tol_deg
        lim = (ang > limit - 0.5) & ~b
        area_45 += float(ar[lim].sum())
        if lim.any():
            max_at_limit = max(max_at_limit, float(ang[lim].max()))
        if b.any():
            area_bad += float(ar[b].sum())
            cc = cen[b]
            r = np.hypot(cc[:, 1], cc[:, 2])
            th = np.degrees(np.arctan2(cc[:, 1], cc[:, 2])) % 360
            bad_faces.append(dict(area_mm2=round(float(ar[b].sum()), 3), max_deg=round(float(ang[b].max()), 2),
                                  geom=str(f.geom_type).split(".")[-1],
                                  x_mm=[round(float(cc[:, 0].min()), 2), round(float(cc[:, 0].max()), 2)],
                                  r_mm=[round(float(r.min()), 2), round(float(r.max()), 2)],
                                  theta_deg=round(float(np.median(th)), 1)))
    bad_faces.sort(key=lambda d: -d["area_mm2"])
    return dict(limit_deg=limit, tol_deg=tol_deg, area_mm2=round(area_bad, 3), area_at_45_mm2=round(area_45, 1),
                at_limit_max_deg=round(max_at_limit, 2),
                bed_area_mm2=round(area_bed, 1), n_faces=len(bad_faces), faces=bad_faces[:15])


def stl_check(path: Path) -> dict:
    import trimesh
    m = trimesh.load(path)
    zmin = m.bounds[0][2]
    n, a, c = m.face_normals, m.area_faces, m.triangles_center
    bed = (c[:, 2] - zmin < 0.02) & (n[:, 2] < -0.999)
    return dict(watertight=bool(m.is_watertight), volume_cm3=round(float(m.volume) / 1000, 3),
                euler=int(m.euler_number), bbox_mm=[round(float(v), 2) for v in m.extents],
                z_min=round(float(zmin), 4), bed_contact_mm2=round(float(a[bed].sum()), 1),
                n_faces=int(len(m.faces)))


def _volume(shape, eps: float = 1e-7) -> float:
    """Volume with a tight integration tolerance (the default reads ~2e-4 off on these solids)."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape.wrapped, g, eps, True)
    return g.Mass()


def robustness_probe(shape) -> dict:
    """Split the solid by 6 planes through the axis; the halves must add up to the whole (a
    valid-but-broken solid fails this, see revolve_nodes). Precise volumes."""
    V = _volume(shape)
    errs = []
    for th in range(0, 180, 30):
        t = math.radians(th + 7.0)
        nrm = (0.0, -math.cos(t), math.sin(t))
        parts_ = [shape.intersect(halfspace((0, 0, 0), nrm)), shape.intersect(halfspace((0, 0, 0), tuple(-v for v in nrm)))]
        if any(p_ is None for p_ in parts_):
            errs.append(1.0)
            continue
        v = sum(_volume(x) for p_ in parts_ for x in p_)
        errs.append(abs(v - V) / V)
    return dict(max_rel_err=float(max(errs)), ok=bool(max(errs) < 1e-4))


def floor_wall_check(shape, x_range) -> dict:
    """Shell left under the pocket floors: rays from just under the floor along its inward
    normal, at every table row in x_range and 5 azimuths per row inside the side walls."""
    pp = pocket_profile()
    x, fl = pp["table"]["x"], pp["cad_floor"]
    lo, hi = pp["th_lo"][1:-1], pp["th_hi"][1:-1]
    sf = np.gradient(fl, x)
    O, D, meta = [], [], []
    for k in range(len(x)):
        if not (x_range[0] + 0.3 <= x[k] <= x_range[1] - 0.3):
            continue
        if r_out(x[k]) - fl[k] < 0.3:
            continue
        m_deg = math.degrees(1.5 / fl[k])
        for tb_ in BLADE_THETA:
            for f_ in np.linspace(0.0, 1.0, 5):
                th = tb_ + lo[k] + m_deg + f_ * (hi[k] - lo[k] - 2 * m_deg)
                e = np.array([0.0, math.sin(math.radians(th)), math.cos(math.radians(th))])
                nn = np.array([sf[k], 0.0, 0.0]) - e
                nn /= np.linalg.norm(nn)
                p0 = np.array([x[k], 0.0, 0.0]) + fl[k] * e
                O.append(p0 + 0.02 * nn)
                D.append(nn)
                meta.append((x[k], th))
    d = _occ_rays(shape, O, D) + 0.02
    k = int(np.argmin(d))
    return dict(n_rays=len(d), min_mm=round(float(d.min()), 3), at_x_theta=[round(meta[k][0], 2), round(meta[k][1] % 360, 1)],
                p05_mm=round(float(np.percentile(d, 5)), 3))


def gland_engagement(shape) -> dict:
    """Thread engagement in solid material along the 72 generatrices of the M12 major diameter,
    from the spot face inward (the bore is tapped through block, wall and plug)."""
    F = gland_frame()
    S, A = F["S"], F["A"]
    u = np.cross(A, [0.0, 1.0, 0.0])
    u /= np.linalg.norm(u)
    v = np.cross(A, u)
    O, D = [], []
    for ph in np.linspace(0, 2 * math.pi, 72, endpoint=False):
        O.append(S - 0.05 * A + 6.0 * (math.cos(ph) * u + math.sin(ph) * v))
        D.append(-A)
    d = np.minimum(_occ_rays(shape, O, D) + 0.05, P.GLAND["thread_length"])
    # where the drill actually removes material (the drill tip itself ends in the open cavity)
    plain = _front_hull().fuse(_gland_block().cut(_front_core())).cut(_front_cavity())
    xs = inter(_gland_bore(), plain).bounding_box().max.X
    return dict(min_mm=round(float(d.min()), 2), mean_mm=round(float(d.mean()), 2),
                layout_min_mm=P.GLAND["thread_engagement_min"], bore_x_max_mm=round(float(xs), 2))


def blade_pocket_check(S: dict, step: float = 0.4) -> dict:
    """Stowed blades (hydroegg.layout surface model, all 3) against the CAD pockets: analytic
    clearances to the cutter's floor and side walls, then OCC point classification of every
    suspect point (outside a pocket's angular range, within 0.25 mm of a floor, or near a detent
    boss) against both shells."""
    pp = pocket_profile()
    xt, fl, lo, hi = pp["x"], pp["floor"], pp["th_lo"], pp["th_hi"]
    g = stow_detent_geometry()
    n_pts = n_sus = n_in = 0
    fa = fr = sc = np.inf
    for i, tb in enumerate(BLADE_THETA):
        B = stowed_blade_points(i, step)
        x, r = B[:, 0], np.hypot(B[:, 1], B[:, 2])
        rel = (np.degrees(np.arctan2(B[:, 1], B[:, 2])) - tb + 180.0) % 360.0 - 180.0
        m = (x <= P.STATIC_TAIL_X + 0.5) & (r < profile.radius(np.clip(x, 0, P.HULL_L)) + 0.05)
        f, a, b = np.interp(x, xt, fl), np.interp(x, xt, lo), np.interp(x, xt, hi)
        inr = (rel >= a) & (rel <= b) & (x >= xt[0])
        q = m & inr
        fc = r - f
        if (q & (x <= 100)).any():
            fa = min(fa, float(fc[q & (x <= 100)].min()))
        if (q & (x > 100)).any():
            fr = min(fr, float(fc[q & (x > 100)].min()))
        if q.any():
            sc = min(sc, float((np.minimum(rel - a, b - rel)[q] * np.radians(1.0) * r[q]).min()))
        # near a detent boss: blade-0 frame, distance along the boss axis to the boss face plane
        B0 = np.column_stack([x, r * np.sin(np.radians(rel)), r * np.cos(np.radians(rel))])   # blade 0 frame
        rel0 = B0 - g["contact"]
        s_ax = rel0 @ g["n"]
        perp = np.linalg.norm(rel0 - np.outer(s_ax, g["n"]), axis=1)
        boss = (perp < STOW_DETENT["boss_r"] + 0.6) & (s_ax - g["s_face"] < 0.3)
        sus = m & (~inr | (fc < 0.25) | boss)
        Bs = B[sus]
        n_pts += int(m.sum())
        n_sus += int(sus.sum())
        n_in += int(_occ_inside(S["aft_shell"], Bs[Bs[:, 0] > X_TIP - 0.5]).sum())
        n_in += int(_occ_inside(S["front_shell"], Bs[Bs[:, 0] < EQ + 0.5]).sum())
    return dict(n_points=n_pts, n_suspects=n_sus, n_inside=n_in, floor_clear_aero_mm=round(fa, 3),
                floor_clear_root_mm=round(fr, 3), side_clear_mm=round(sc, 3))


def pad_pocket_arc(step: float = 0.05) -> tuple:
    """Smallest arc (mm, at the hull radius) between a yawed rail pad side and a CAD pocket bound
    over the pad length: the layout's definition (pad walls are planes parallel to the radial line
    at the fin angle, so a wall at lateral offset lat meets the hull at theta + asin(lat / R))."""
    pp = pocket_profile()
    F = fin_frame()
    x_p = F["o"][0] + FIN["yaw_pivot_chord_frac"] * FIN["chord"] * F["ci"]
    best = (np.inf, None)
    for x in np.arange(P.POCKETS["x"][0], EQ + 1e-9, step):
        R = r_out(x)
        dx = x - x_p
        a = math.degrees(math.asin((-FIN["rail_w"] / 2 + dx * F["si"]) / F["ci"] / R))
        b = math.degrees(math.asin((FIN["rail_w"] / 2 + dx * F["si"]) / F["ci"] / R))
        for tf in FIN["theta_deg"]:
            for tb in BLADE_THETA:
                lo = float(np.interp(x, pp["x"], pp["th_lo"])) + tb
                hi = float(np.interp(x, pp["x"], pp["th_hi"])) + tb
                for off in (-360.0, 0.0, 360.0):
                    g = max(lo + off - (tf + b), (tf + a) - (hi + off))
                    if g > -90 and math.radians(g) * R < best[0]:
                        best = (math.radians(g) * R, (tf, tb, round(float(x), 2)))
    return best


def interface_checks(S: dict) -> list:
    """Interface checks measured on the solids (value, limit, pass)."""
    out = []

    def add(check, value, ok, limit=None):
        out.append(dict(check=check, value=value if isinstance(value, str) else f"{value}",
                        limit="" if limit is None else f"{limit}", ok=bool(ok)))

    front, aft = S["front_shell"], S["aft_shell"]
    pp = pocket_profile()
    smp = pp["sample"]
    add("pocket table convention: DETENTS stowed point (x %.2f, az %.2f psi) = theta %.2f inside the pocket"
        % (smp["x"], smp["az"], -smp["az"]), "asserted", True)
    # stowed blades vs the shells
    bp = blade_pocket_check(S)
    add("stowed blades (layout.py surface, 3 blades): points inside either shell (OCC on %d suspects of %d)"
        % (bp["n_suspects"], bp["n_points"]), bp["n_inside"], bp["n_inside"] == 0, 0)
    add("stowed aero span above the CAD pocket floor (min radial clearance, x <= 100), mm", bp["floor_clear_aero_mm"],
        bp["floor_clear_aero_mm"] >= P.POCKET_CLEAR - 0.01, P.POCKET_CLEAR)
    add("stowed root above the CAD pocket floor (min radial clearance, x > 100), mm", bp["floor_clear_root_mm"],
        bp["floor_clear_root_mm"] >= 0.0, "> 0 (floor cut deeper by the 45 deg rule)")
    add("stowed blade to the CAD pocket side walls (min arc), mm", bp["side_clear_mm"], bp["side_clear_mm"] >= 0.5, 0.5)
    g = stow_detent_geometry()
    add("stowed-detent boss face to the stowed blade (normal, min over the face), mm", round(g["blade_min_s"] - g["s_face"], 3),
        g["blade_min_s"] - g["s_face"] >= 0.05, ">= 0.05")
    add("stowed-detent static magnet axis vs the blade face normal at the shared mouth point (45 deg print limit), deg",
        round(g["axis_mismatch_deg"], 2), g["axis_mismatch_deg"] <= 8.0, "<= 8")
    # floor walls
    fw = floor_wall_check(front, (P.POCKETS["x"][0], EQ))
    aw = floor_wall_check(aft, (EQ, P.STATIC_TAIL_X))
    add("front shell wall under the pocket floors (rays along the floor normal), min mm", fw["min_mm"],
        fw["min_mm"] >= P.POCKET_FLOOR_WALL - 0.05, P.POCKET_FLOOR_WALL)
    add("aft shell wall under the pocket floors (rays along the floor normal), min mm", aw["min_mm"],
        aw["min_mm"] >= P.POCKET_FLOOR_WALL - 0.05, P.POCKET_FLOOR_WALL)
    # fin pads vs pockets
    arc, where = pad_pocket_arc()
    add("fin rail pads to the pocket edges (arc at the hull radius, layout definition; at x %.2f), mm" % where[2],
        round(arc, 3), arc >= 2.0, "2.0 (layout 2.22 sampled, 2.20 independent)")
    hull = _front_hull()
    pads = [p_.cut(hull) for p_ in _pads()[0]]
    pk = [inter(q, hull) for q in pocket_cutters()]
    dpad = min(pd.distance_to(q) for pd in pads for q in pk)
    add("fin rail pads (above the hull) to the pocket void (3D min distance, information), mm", round(dpad, 3),
        dpad >= 1.5, "1.5")
    # generator bay
    G = P.GENERATOR
    bay = C.axis_cyl(G["x_bay_fwd"], G["x_bay_fwd"] + G["bay_L"], G["bay_D"] / 2)
    db = min(front.distance_to(bay), aft.distance_to(bay))
    add("generator bay D%.1f x %g to the shells, mm" % (G["bay_D"], G["bay_L"]), round(db, 3), db >= 0.5, 0.5)
    # carrier lock
    L = CARRIER_LOCK
    ring = lambda x: r_in(x) - P.CARRIER["ring"]["t_radial"]                     # noqa: E731
    xa = L["tab_x"][1] + L["clear_axial"]
    add("carrier tab overlap behind the ring lip (tab tip - lip inner r at the groove aft wall), mm",
        round(L["tab_r_out"] - ring(xa), 3), L["tab_r_out"] - ring(xa) >= 0.8, 0.8)
    plate_cl = min(ring(x) - L["plate_r"] for x in np.linspace(L["plate_x"][0] - 0.2, L["plate_x"][1] + 0.2, 50))
    add("carrier plate rim to the ring (radial, plate x +-0.2), mm", round(plate_cl, 3), plate_cl >= 0.1, 0.1)
    car = _carrier_cutters()
    hull_skin = min(c_.distance_to(_front_outer_skin()) for c_ in car[::2])
    add("carrier groove to the outer hull surface (wall left), mm", round(hull_skin, 3), hull_skin >= P.MIN_WALL, P.MIN_WALL)
    # gland
    ge = gland_engagement(front)
    add("gland M12 thread engagement in solid material (min over 72 generatrices), mm", ge["min_mm"],
        ge["min_mm"] >= 6.0, "6.0 (layout 6.05)")
    add("gland bore (cable exit) x max, ahead of the carrier ring at x %.1f" % P.CARRIER["ring"]["x"][0], ge["bore_x_max_mm"],
        ge["bore_x_max_mm"] <= P.CARRIER["ring"]["x"][0] - 1.0, "<= %.1f" % (P.CARRIER["ring"]["x"][0] - 1.0))
    # joint
    vi = sum(x.volume for x in (front.intersect(aft) or []))
    add("front shell / aft shell interference volume (locked), mm3", round(vi, 4), vi < 0.5, 0.5)
    lugs = _lugs()
    dl = min(front.distance_to(l_) for l_ in lugs)
    add("locked lugs to the front shell (hard stop contact), mm", round(dl, 4), dl < 0.01, "contact")
    add("seal bore / spigot radial clearance, mm", round(R_BORE - R_SPIGOT, 3), abs(R_BORE - R_SPIGOT - P.ORING_EQUATOR["radial_clearance"]) < 1e-6,
        P.ORING_EQUATOR["radial_clearance"])
    # tube vs rotor keep-out, tail end
    kp = json.loads((DATA / "layout_rotor.json").read_text())["pocket_summary"]["keepout_min_mm"]
    add("bearing tube outer r under the tail notch keep-out, mm", BT["r_out_rear"], BT["r_out_rear"] < kp, "< %.2f" % kp)
    bb = aft.bounding_box()
    add("aft shell ends at STATIC_TAIL_X (hub gap %.1f to the hub front)" % P.HUB_GAP, round(bb.max.X, 4),
        abs(bb.max.X - P.STATIC_TAIL_X) < 1e-6, P.STATIC_TAIL_X)
    add("claw B step OD %.0f to the tube mouth (radial), mm" % P.COUPLER["claw_b_step_od"],
        round(R608 - P.COUPLER["claw_b_step_od"] / 2, 3), R608 - P.COUPLER["claw_b_step_od"] / 2 >= 0.3, 0.3)
    # cutaway keeps the fin pads
    sec = cutaway_sector(True)
    vp = sum(x.volume for pd in pads for x in (pd.intersect(sec) or []))
    add("cutaway sector theta %g..%g cuts a fin pad (volume), mm3" % CUTAWAY_THETA, round(vp, 3), vp < 0.01, 0)
    return out


def _front_outer_skin() -> Solid:
    """Solid outside the front hull (a thick shell around it) for wall-thickness distances."""
    big = revolve_nodes([("L", (-30.0, 0.0)), ("L", (-30.0, 80.0)), ("L", (EQ + 5, 80.0)), ("L", (EQ + 5, 0.0))])
    return big.cut(_front_hull())


def section_figure(S: dict, path: Path) -> None:
    """2D meridian sections of both shells with the internals' envelopes from params."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, PathPatch, Rectangle
    from matplotlib.path import Path as MPath

    planes = [(0.0, "theta 0 (top: blade 0 pocket, rib) / 180 (bottom: gland)"),
              (LOCK[0], "theta %g (top: locked lug, channel stop) / %g (bottom)" % (LOCK[0], LOCK[0] + 180)),
              (FIN["theta_deg"][0], "theta %g (top: fin pad, groove, cross pin) / %g (bottom: fin pad)"
               % tuple(FIN["theta_deg"]))]
    fig, axs = plt.subplots(len(planes), 1, figsize=(13, 5.0 * len(planes)))
    cols = {"front_shell": "#C9B48A", "aft_shell": "#8F7E5C"}
    for ax, (th, title) in zip(axs, planes):
        t = math.radians(th)
        e_s = Vector(0.0, math.sin(t), math.cos(t))
        pl = Plane(origin=(0, 0, 0), x_dir=(1, 0, 0), z_dir=(0, math.cos(t), -math.sin(t)))
        for n in ("front_shell", "aft_shell"):
            res = S[n].intersect(pl) or []
            for f in [ff for item in res for ff in item.faces()]:
                verts, codes = [], []
                for w in [f.outer_wire(), *f.inner_wires()]:
                    q = w.positions(list(np.linspace(0, 1, 1500)))
                    ring = [(p.X, p.dot(e_s)) for p in q]
                    verts += ring + [ring[0]]
                    codes += [MPath.MOVETO] + [MPath.LINETO] * (len(ring) - 1) + [MPath.CLOSEPOLY]
                ax.add_patch(PathPatch(MPath(verts, codes), fc=cols[n], ec="#3b3428", lw=0.5, zorder=3))
        # references from params (dashed envelopes)
        G = P.GENERATOR
        for sgn in (1, -1):
            kw = dict(fill=False, lw=0.8, ls=(0, (3, 2)), zorder=5)
            ax.add_patch(Rectangle((G["x_bay_fwd"], 0), G["bay_L"], sgn * G["bay_D"] / 2, ec="#2f6fb5", **kw))
            ax.add_patch(Rectangle((P.CARRIER["x0"], 0), P.CARRIER["t"], sgn * P.CARRIER["r"], ec="#5a3fb0", **kw))
            for b in (P.BEARING, P.BEARING_ALT, P.LIP_SEAL):
                ax.add_patch(Rectangle((b["x0"], sgn * b["d"] / 2), b["B"], sgn * (b["D"] - b["d"]) / 2, ec="#222", **kw))
            ax.add_patch(Rectangle((P.SHAFT["x0"], 0), P.SHAFT["L"], sgn * P.SHAFT["d"] / 2, ec="#777", **kw))
            cx = 0.5 * (GR["x_aft"] + GR["x_fwd_bottom"])
            ax.add_patch(Circle((cx, sgn * (R_GROOVE + P.ORING_EQUATOR["cs"] / 2)), P.ORING_EQUATOR["cs"] / 2,
                                fc="#111", ec="none", zorder=6))
        xs = np.linspace(0, P.HULL_L, 600)
        for sgn in (1, -1):
            ax.plot(xs, sgn * profile.radius(xs), color="#999", lw=0.6, ls=":", zorder=1)
        ax.axvline(EQ, color="#999", lw=0.6, ls=":")
        ax.axhline(0, color="#bbb", lw=0.5)
        ax.set_aspect("equal")
        ax.set_xlim(-16, 112)
        ax.set_ylim(-48, 48)
        ax.set_title("Shell section, plane through the axis at " + title, loc="left", fontsize=10)
        ax.set_xlabel("x (mm), flow ->")
        ax.set_ylabel("r (mm)")
    fig.text(0.01, 0.005, "front shell light, aft shell dark; dashed: generator bay, carrier plate, 608 / 688 / lip seal, "
             "shaft (params.py); black dot: O-ring cross-section; dotted: hull profile", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main(argv: list[str]) -> None:
    t0 = time.time()
    force = "--force" in argv
    S = build_all(force=force)
    print(f"built {time.time() - t0:.1f} s", flush=True)
    for d in (OUT_STEP, OUT_STL, OUT_CHECKS, WORK):
        d.mkdir(parents=True, exist_ok=True)
    report = dict(module="hydroegg.cad.shells", cache_key=_cache_key(), parts=[], conventions=dict(
        theta="about +X from +Z toward +Y; pocket table theta = -psi (asserted)",
        print="both shells print +X down (front: equator face on the bed, nose up; aft: tail end face on the bed, spigot up)",
        step="world pose (assembled)", stl="print orientation, z = 0 bed"))
    for d in parts():
        n = d["name"]
        export_step(d["shape"], str(OUT_STEP / f"{n}.step"))
        export_stl(d["print"], str(OUT_STL / f"{n}.stl"), tolerance=0.03, angular_tolerance=0.1)
        x_bed = d["shape"].bounding_box().max.X
        rec = dict(C.check(d["shape"], n), print_check=C.check(d["print"], n), fits_bed=C.fits_bed(d["print"]),
                   stl=stl_check(OUT_STL / f"{n}.stl"), overhang=overhang_scan(d["shape"], x_bed),
                   robustness=robustness_probe(d["shape"]), role=d["role"], notes=d["notes"])
        report["parts"].append(rec)
        print(n, json.dumps({k: rec[k] for k in ("valid", "n_solids", "volume_cm3")}), "stl", rec["stl"]["watertight"],
              "overhang", round(rec["overhang"]["area_mm2"], 2), "probe", rec["robustness"]["ok"],
              f"{time.time() - t0:.1f} s", flush=True)
    report["interface_checks"] = interface_checks(S)
    for c in report["interface_checks"]:
        print(("PASS " if c["ok"] else "FAIL ") + c["check"] + ": " + c["value"], flush=True)
    g = stow_detent_geometry()
    report["stow_detent"] = dict({k: (np.round(v, 4).tolist() if isinstance(v, np.ndarray) else v) for k, v in g.items()},
                                 **STOW_DETENT)
    report["carrier_lock"] = CARRIER_LOCK
    report["fin_rail"] = FIN_RAIL
    report["cutaway_theta_deg"] = list(CUTAWAY_THETA)
    pp = pocket_profile()
    tb_, i0_ = pp["table"], pp["tip_rows"]
    wid = np.maximum(pp["th_hi"][1:-1] - (-tb_["az_lo"]), (-tb_["az_hi"]) - pp["th_lo"][1:-1])
    report["pockets"] = dict(
        table_rows=len(tb_["x"]), loft_rows_kept=int(len(_decimate(pp["x"], (pp["floor"], pp["th_lo"], pp["th_hi"]),
                                                                     (0.01, 0.02, 0.02)))),
        side_wall_widening_deg_max=round(float(wid[i0_:].max()), 3),
        squared_tip_rows_x=[float(v) for v in tb_["x"][:i0_]],
        floor_deepened_mm_max=round(float((tb_["floor"] - pp["cad_floor"]).max()), 3),
        floor_print_deg=FLOOR_PRINT_DEG,
        note="bounds = table az (theta = -psi); widened only where a >= 0.5 mm side wall would widen aft faster "
             "than 45 deg; rows ahead of the first 0.5 mm deep row take its bounds (square tip); floor = table "
             "floor re-limited to %g deg rise; liner = floor - (POCKET_FLOOR_WALL + 0.15) along the floor normal"
             % FLOOR_PRINT_DEG)
    section_figure(S, WORK / "section_2d.png")
    (OUT_CHECKS / "shells.json").write_text(json.dumps(report, indent=1, default=float))
    print(f"checks written {time.time() - t0:.1f} s", flush=True)
    if "--no-render" not in argv:
        # mid-gray material: the ivory 'shell' material washes out edges in these review renders
        renders = [("preview_parts", "shells_iso", "view=iso&zoom=1.3", "default"),
                   ("preview_section", "shells_section", "view=-25,20&zoom=1.3", "default"),
                   ("preview_cutaway", "shells_cutaway", "view=hero&zoom=1.3", "default"),
                   ("preview_print", "shells_print", "view=-30,25", "default")]
        for fn, tag, q, mat in renders:
            for attempt in (1, 2):                      # headless Chrome sometimes times out on navigate
                r = subprocess.run([sys.executable, str(REPO / "scripts/preview.py"), f"hydroegg.cad.shells:{fn}",
                                    "--out", str(WORK / f"{tag}.png"), "--q", q], cwd=REPO,
                                   env=dict(__import__("os").environ, SHELLS_PREVIEW_MAT=mat))
                if r.returncode == 0:
                    break
                print(f"render {tag} failed (attempt {attempt})", flush=True)
    print(f"done {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
