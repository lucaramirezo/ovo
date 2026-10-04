"""Drivetrain of the HydroEgg: shaft, claw coupling, motor carrier, GL40 replica and the
bought parts of the power path (bearings, lip seal, shaft collar, O-ring, cable gland,
screws and heat-set inserts).

Run ``scripts/cadrun.sh .venv/bin/python -m hydroegg.cad.drivetrain`` to build every part,
write out/cad/step/<part>.step (world pose), out/cad/stl/<part>.stl (print pose, printed
parts only) and out/cad/checks/drivetrain.json.

Sources: hydroegg/params.py (SHAFT, COUPLER, CARRIER, GENERATOR, BEARING, BEARING_ALT,
LIP_SEAL, ORING_EQUATOR, JOINT, GLAND, SHAFT_STOP, INSERT_M3, clearances) and
data/hardware.json (GL40 hole patterns and solder-pad keep-out, collar, gland).

Frame (params.py): X = body axis, flow +X, nose at x = 0, Z up; theta about +X from +Z
toward +Y. Rotor azimuth psi = -theta (rotor half convention).

Design notes (choices made here, see DESIGN below and the checks JSON):

* Claw coupling. The layout gives the stack between the generator rotor face (x 46.5) and
  the 608 front face (x 63): claw A, 1.0 mm axial float, claw B with an OD 20 step over its
  aft end and a 0.5 mm OD 12 shim on the 608 inner race. Here the two claws interleave:
  claw A is a 2 mm flange on the rotor face with two dogs (theta 0 / 180) that carry the
  counterbores of its two M3 x 4 screws (rotor-face holes at 0 / 180; the holes at 90 / 270
  stay free because claw B's dogs pass over them). Claw B is a printed sleeve pressed (3
  crush ribs) over the Ruland MSC-8-SS set-screw collar, which bears on the shim and carries
  the thrust; its two dogs (theta 90 / 270) reach forward between claw A's dogs. Dog tips
  stop 1.0 mm short of the opposite claw (axial float) and the flanks run 0.3 mm apart
  (0.15 per flank), V lead-in on every dog tip.
* The collar replaces the layout's alternative (retaining-ring groove + 2 grub screws): no
  groove to grind into a 56 HRC shaft. The collar set screw lands on a 0.5 mm flat.
* Display build: the replica is two printed parts, a static stator (bolted to the carrier
  through the adapter spacer, M3 heat-set inserts on PCD 24 like the real motor's threads)
  and a bell that turns on a printed journal, so the rotor of the display model really
  spins through the real coupling. In the closed hull the bell can move aft only by the
  1 mm claw float, far less than its 7.5 mm journal, so it cannot come off.
* Lip seal (v2) goes in AFTER the shaft, from the tail end, over a thin brass installation
  sleeve (K&S 9 x 0.45 mm tube, ID 8.1) slid over the protruding shaft end; that covers the
  hub D-flat and the tapped end. A sleeve cannot be used the layout's way (seal first, shaft
  pushed through it from the front), because nothing thicker than 8 mm passes the 688 bore.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from build123d import (
    Align,
    Axis,
    Box,
    Circle,
    Ellipse,
    Location,
    Plane,
    Polygon,
    Rectangle,
    RegularPolygon,
    Shape,
    Solid,
    Compound,
    Vector,
    revolve,
    export_step,
    export_stl,
)

from .. import params as P
from .common import at_azimuth, check, fits_bed, interference, min_distance, on_bed

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out" / "cad"
WORK = ROOT / "out" / "work" / "cad" / "drivetrain"
HW_JSON = ROOT / "data" / "hardware.json"

S, CP, CR, BT, B8 = P.SHAFT, P.COUPLER, P.CARRIER, P.BEARING, P.BEARING_ALT
G, LS, OR, J = P.GENERATOR, P.LIP_SEAL, P.ORING_EQUATOR, P.JOINT


def _hw():
    try:
        return json.loads(HW_JSON.read_text())
    except Exception:  # noqa: BLE001
        return {}


_HW = _hw()


def _find_key(d, key):
    if isinstance(d, dict):
        if key in d:
            return d[key]
        for v in d.values():
            r = _find_key(v, key)
            if r is not None:
                return r
    elif isinstance(d, list):
        for v in d:
            r = _find_key(v, key)
            if r is not None:
                return r
    return None


_PAD = (_find_key(_HW, "solder_pad_keepout") or {}).get("value") or {}
PAD_WINDOW = _PAD.get("cad_window") or {"radial_from_mm": 12.0, "radial_to_mm": 20.5,
                                         "tangential_half_width_mm": 6.0}
PAD_BLOCK = {"r0": _PAD.get("radial_from_mm", 12.93), "r1": _PAD.get("r_range_mm", [12.93, 19.72])[1],
             "hw": _PAD.get("tangential_half_width_mm", 5.09)}

# --------------------------------------------------------------------------
# Design values owned by this module (choices, not catalogue data)
# --------------------------------------------------------------------------
DESIGN = dict(
    # claw coupling (interleaved dogs)
    claw_a_flange_t=2.0,          # flange on the rotor face, under the screw heads
    claw_a_od=32.4,               # flange / dog lobes: counterbore D6 at r 12 + 1.2 wall
    claw_dog_r_in=4.6,            # dog inner radius (shaft r 4 + 0.6)
    claw_dog_half_deg=45.0,       # 2 dogs per claw, 90 deg pitch between A and B dogs
    claw_flank_gap=0.30,          # total gap between facing flanks (0.15 per flank)
    claw_tip_chamfer=0.6,         # V lead-in on every dog tip
    claw_b_floor_t=1.2,           # claw B floor ahead of the collar (MIN_WALL)
    claw_b_cone_x=(60.5, 61.5),   # OD 22 -> OD 20 cone (45 deg, prints aft face down)
    collar_pocket_d=16.2,         # collar pocket, 3 crush ribs to D 15.9 (PRESS_FIT per side)
    crush_rib_r=0.5,
    setscrew_access_d=4.0,        # teardrop access hole for the collar set screw key
    setscrew_theta=0.0,           # assembled azimuth of the collar set screw / shaft flat
    collar_flat_x=(55.0, 62.0),   # shaft flat under the collar set screw (x 58.5)
    hub_flat_x0=108.4,            # hub D-flat from here to the shaft end (hub front 108.6)
    shaft_chamfer=0.5,
    # carrier
    carrier_cbore=(6.2, 1.7),     # ISO 7380 M3 button head (dk 5.7, k 1.65): D x depth
    carrier_clear_d=3.4,          # M3 clearance
    carrier_slot=(38.14, 40.0, 3.2),   # radial slots for PCD 38.14 (GM4108H) .. PCD 40 (GL40 M2.5)
    carrier_window_corner_r=2.0,
    # carrier twist-lock tabs: taken from hydroegg.cad.shells.CARRIER_LOCK (shell side owns the
    # ring groove and notches); this copy is the fallback if shells.py cannot be imported
    carrier_lock_fallback=dict(n_tabs=3, tab_w_deg=16.0, tab_x=[20.5, 22.5], tab_r_out=27.2,
                               tab_front_chamfer_45=1.0, lock_theta_deg=[90.0, 210.0, 330.0], twist_deg=30.0),
    # GL40 replica (display)
    stator_base=(44.0, 3.0),      # D x t of the stator flange (front face x 25)
    stator_core=(38.0, 18.2),     # stator stack OD / ID
    stator_core_x1=40.0,
    stator_slots=24,
    post_d=(15.0, 12.0),          # stator post lower / journal diameter
    post_shoulder_x=37.2,
    post_top_x=45.0,
    bell_wall=1.6,
    bell_face_t=2.5,
    bell_rim_gap=0.5,             # bell rim to stator flange (the visible air-gap line)
    bell_sleeve_od=17.0,
    bell_sleeve_x0=37.5,
    bell_boss_d=7.2,              # M3 insert boss (hole 4.0 + 2 x 1.6)
    bell_boss_x0=40.5,
    insert_hole_depth=6.0,
    spacer_od=46.0,
)

def carrier_lock() -> dict:
    """Twist-lock tab spec shared with the front shell (shells.CARRIER_LOCK is the source)."""
    try:
        from .shells import CARRIER_LOCK as L
        src = "hydroegg.cad.shells.CARRIER_LOCK"
    except Exception:  # noqa: BLE001
        L, src = DESIGN["carrier_lock_fallback"], "drivetrain fallback copy"
    L = dict(L)
    L.setdefault("entry_theta_deg", [(t + L["twist_deg"]) % 360.0 for t in L["lock_theta_deg"]])
    L["source"] = src
    return L


THETA_SCREWS = [45.0, 135.0, 225.0, 315.0]     # motor pattern clocked 45 deg (CARRIER)
THETA_PAD = 180.0                              # solder-pad block faces the cable notch

# --------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------


def rev(pts) -> Solid:
    """Revolve a closed (x, r) polygon (r >= 0) about the X axis."""
    face = Plane.XZ * Polygon(*[(float(x), float(r)) for x, r in pts], align=None)
    return revolve(face, Axis.X)


def tube(x0: float, x1: float, ri: float, ro: float, ch: float = 0.0) -> Solid:
    """Tube on the axis; ch chamfers the four corners."""
    c = min(ch, (ro - ri) / 2.5, (x1 - x0) / 2.5) if ch else 0.0
    if c <= 0:
        return rev([(x0, ri), (x1, ri), (x1, ro), (x0, ro)])
    return rev([(x0, ri + c), (x0 + c, ri), (x1 - c, ri), (x1, ri + c), (x1, ro - c), (x1 - c, ro),
                (x0 + c, ro), (x0, ro - c)])


def PX(x: float) -> Plane:
    """Plane normal to +X at x; local (u, v) = world (y, z)."""
    return Plane(origin=(x, 0, 0), x_dir=(0, 1, 0), z_dir=(1, 0, 0))


def uv(r: float, th_deg: float) -> tuple[float, float]:
    t = math.radians(th_deg)
    return (r * math.sin(t), r * math.cos(t))


def xcyl(x0: float, x1: float, r: float, th_deg: float = 0.0, rr: float = 0.0) -> Solid:
    """Cylinder parallel to X from x0 to x1, radius r, axis at radius rr / azimuth th."""
    u, v = uv(rr, th_deg)
    return Solid.make_cylinder(r, x1 - x0, Plane(origin=(x0, u, v), x_dir=(0, 1, 0), z_dir=(1, 0, 0)))


def blind_hole(x_open: float, x_end: float, r: float, th: float = 0.0, rr: float = 0.0) -> Solid:
    """Hole parallel to X from the open face x_open to x_end with a 90 deg (45 deg wall) point at the
    end, so a blind hole printed with its open end on the bed needs no bridge."""
    sgn = 1.0 if x_end > x_open else -1.0
    a, b = (x_open - 0.1 * sgn), x_end
    body = xcyl(min(a, b), max(a, b), r, th, rr)
    u, v = uv(rr, th)
    cone = Solid.make_cone(r, 0.02, r, Plane(origin=(x_end, u, v), x_dir=(0, 1, 0), z_dir=(sgn, 0, 0)))
    return body.fuse(cone)


def radial_cyl(x: float, th_deg: float, r: float, r0: float, r1: float) -> Solid:
    """Cylinder of radius r on a radial line at (x, theta) from radius r0 to r1."""
    t = math.radians(th_deg)
    d = Vector(0, math.sin(t), math.cos(t))
    o = Vector(x, 0, 0) + d * r0
    return Solid.make_cylinder(r, r1 - r0, Plane(origin=o, x_dir=(1, 0, 0), z_dir=d))


def radial_teardrop(x: float, th_deg: float, d: float, r0: float, r1: float, down=(-1, 0, 0)) -> Solid:
    """Radial hole of diameter d with a 45 deg teardrop roof pointing away from `down`
    (the print 'down' direction in world), so it prints without support."""
    t = math.radians(th_deg)
    n = Vector(0, math.sin(t), math.cos(t))
    up = -Vector(*down)
    up = (up - n * up.dot(n)).normalized()
    side = n.cross(up).normalized()
    rr = d / 2
    # teardrop: circle + a square rotated 45 deg whose corner points 'up'
    pl = Plane(origin=Vector(x, 0, 0) + n * r0, x_dir=side, z_dir=n)
    tip = rr * math.sqrt(2)
    poly = [(0, tip), (rr / math.sqrt(2), rr / math.sqrt(2)), (-rr / math.sqrt(2), rr / math.sqrt(2))]
    f = (pl * Circle(rr)).faces()[0]
    tri = (pl * Polygon(*poly, align=None)).faces()[0]
    s1 = Solid.extrude(f, n * (r1 - r0))
    s2 = Solid.extrude(tri, n * (r1 - r0))
    return s1.fuse(s2).clean()


def halfspace(normal, offset: float, size: float = 400.0) -> Solid:
    """Solid half-space {p : p . n >= offset} clipped to a big box."""
    n = Vector(*normal).normalized()
    pl = Plane(origin=n * offset, z_dir=n)
    return Solid.make_box(size, size, size, Plane(origin=pl.origin - pl.x_dir * (size / 2) - pl.y_dir * (size / 2),
                                                  x_dir=pl.x_dir, z_dir=n))


def flat_cut(x0: float, x1: float, n: Vector, r_plane: float, w: float = 20.0, h: float = 10.0) -> Solid:
    """Box that removes everything beyond the plane at distance r_plane along n (n normal to X),
    between x0 and x1: used for shaft flats."""
    pl = Plane(origin=Vector(x0, 0, 0) + n * r_plane, x_dir=(1, 0, 0), z_dir=n)
    return Solid.make_box(x1 - x0, w, h, Plane(origin=pl.origin - pl.y_dir * (w / 2), x_dir=(1, 0, 0), z_dir=n))


def dog(x0: float, x1: float, r_in: float, r_out: float, th_c: float, half_deg: float,
        flank_off: float, tip_x: float, tip_ch: float) -> Solid:
    """Claw dog: annular sector about theta th_c (+-half_deg) between x0 and x1, each flank
    moved inward by flank_off (parallel offset), V lead-in chamfer on the tip at tip_x."""
    ring = tube(x0, x1, r_in, r_out)
    out = ring
    for sgn in (+1, -1):
        a = math.radians(th_c + sgn * half_deg)
        # flank plane through the axis at angle a; inward normal toward the dog centre
        e = Vector(0, math.sin(a), math.cos(a))            # radial direction of the flank
        tdir = Vector(0, math.cos(a), -math.sin(a))        # +theta direction at the flank
        n_in = tdir * (-sgn)                               # toward the dog centre
        out = out & halfspace(n_in, flank_off)
    # keep only the sector around th_c (the two half-spaces also keep the opposite sector)
    c = math.radians(th_c)
    out = out & halfspace(Vector(0, math.sin(c), math.cos(c)), 0.0)
    if tip_ch > 0:
        edges = [e for e in out.edges() if str(e.geom_type).endswith("LINE")
                 and abs(e.center().X - tip_x) < 1e-6 and abs(e.start_point().X - e.end_point().X) < 1e-6]
        if edges:
            out = out.chamfer(tip_ch, tip_ch, edges)
    return out


def m3_socket_screw(x_head: float, length: float, th: float, r: float, direction: float = -1.0) -> Solid:
    """ISO 4762 M3 (dk 5.5, k 3) with its head bearing face at x_head; shank goes toward
    direction (-1 = -X)."""
    xh0, xh1 = (x_head, x_head + 3.0) if direction < 0 else (x_head - 3.0, x_head)
    head = xcyl(xh0, xh1, 2.75, th, r)
    xs0, xs1 = (x_head - length, x_head) if direction < 0 else (x_head, x_head + length)
    shank = xcyl(xs0, xs1, 1.45, th, r)
    sock_x = (xh1 - 1.6, xh1) if direction < 0 else (xh0, xh0 + 1.6)
    u, v = uv(r, th)
    hexa = Solid.extrude((Plane(origin=(sock_x[0], u, v), x_dir=(0, 1, 0), z_dir=(1, 0, 0))
                          * RegularPolygon(1.25 / math.cos(math.pi / 6), 6)).faces()[0], (1.6, 0, 0))
    return head.fuse(shank).clean() - hexa


def m3_button_screw(x_head: float, length: float, th: float, r: float) -> Solid:
    """ISO 7380 M3 button head (dk 5.7, k 1.65), bearing face at x_head, shank toward +X."""
    u, v = uv(r, th)
    head = rev([(x_head - 1.65, 0.0), (x_head - 1.65, 1.0), (x_head - 1.2, 2.0), (x_head - 0.5, 2.7),
                (x_head, 2.85), (x_head, 0.0)]).moved(Location((0, u, v)))
    shank = xcyl(x_head, x_head + length, 1.45, th, r)
    hexa = Solid.extrude((Plane(origin=(x_head - 1.65, u, v), x_dir=(0, 1, 0), z_dir=(1, 0, 0))
                          * RegularPolygon(1.0 / math.cos(math.pi / 6), 6)).faces()[0], (1.1, 0, 0))
    return head.fuse(shank).clean() - hexa


def heat_set_insert(x_face: float, th: float, r: float, direction: float = 1.0) -> Solid:
    """CNC Kitchen M3 x 5.7 insert modelled at the printed hole size (D 3.98, no overlap with
    the hole wall) with its M3 bore, top flush with x_face, going toward direction."""
    L = P.INSERT_M3["L"]
    x0, x1 = (x_face, x_face + L) if direction > 0 else (x_face - L, x_face)
    body = xcyl(x0, x1, 3.98 / 2, th, r)
    bore = xcyl(x0 - 0.1, x1 + 0.1, 1.25, th, r)
    # knurl rings for the look: two shallow grooves
    g = []
    for k in (0.3, 0.65):
        xm = x0 + k * L
        g.append(tube(xm - 0.35, xm + 0.35, 1.8, 2.5).moved(Location((0, *uv(r, th)))))
    out = body - bore
    for gg in g:
        out = out - gg
    return out


# --------------------------------------------------------------------------
# Shaft and the parts on it
# --------------------------------------------------------------------------


def shaft() -> Solid:
    d, x0, x1 = S["d"], S["x0"], S["x1"]
    r, c = d / 2, DESIGN["shaft_chamfer"]
    body = rev([(x0, 0), (x0, r - c), (x0 + c, r), (x1 - c, r), (x1, r - c), (x1, 0)])
    depth = S["hub_flat_depth"]
    # hub D-flat facing the grub (psi 60 -> theta -60)
    th_g = -P.ROTOR_HUB["grub_psi_deg"]
    n = Vector(0, math.sin(math.radians(th_g)), math.cos(math.radians(th_g)))
    fx0 = DESIGN["hub_flat_x0"]
    body = body - flat_cut(fx0, x1 + 1.0, n, r - depth)
    # collar flat at the set-screw azimuth
    th_c = DESIGN["setscrew_theta"]
    n2 = Vector(0, math.sin(math.radians(th_c)), math.cos(math.radians(th_c)))
    fa, fb = DESIGN["collar_flat_x"]
    body = body - flat_cut(fa, fb, n2, r - 0.5)
    # M3 x 8 tapped aft end: 2.5 tap drill 10 deep with a 118 deg point (thread not modelled)
    # M3 x 8 tapped end: thread modelled at the major diameter (3.02, so the spinner screw fits with no
    # overlap), 2.5 tap drill to 10 deep with a 118 deg point (0.05 flat, no apex)
    tip = x1 - 10.0 - 1.2 / math.tan(math.radians(59))
    drill = rev([(x1 + 0.1, 0), (x1 + 0.1, 1.81), (x1 - 0.3, 1.51), (x1 - 8.0, 1.51), (x1 - 8.0, 1.25),
                 (x1 - 10.0, 1.25), (tip, 0.05), (tip, 0)])
    return body - drill


def shaft_collar() -> Solid:
    """Ruland MSC-8-SS one-piece set-screw collar 8 x 16 x 8 (stainless), aft face on the shim."""
    x1 = 63.0 - CP["shim_t"]
    x0 = x1 - 8.0
    body = tube(x0, x1, 4.0, 8.0, ch=0.4)
    hole = radial_cyl((x0 + x1) / 2, DESIGN["setscrew_theta"], 1.65, 3.5, 8.5)  # M4 tap (VERIFY size)
    screw = radial_cyl((x0 + x1) / 2, DESIGN["setscrew_theta"], 1.65, 4.0, 6.5)
    sock = radial_cyl((x0 + x1) / 2, DESIGN["setscrew_theta"], 1.0, 5.6, 6.6)
    return (body - hole).fuse(screw - sock).clean()


def shim() -> Solid:
    x1 = 63.0
    return tube(x1 - CP["shim_t"], x1, 4.05, CP["shim_od"] / 2)


def claw_b() -> Solid:
    """Printed claw B: sleeve pressed over the shaft collar (thrust collar), two dogs forward."""
    D = DESIGN
    x_end = 63.0 - CP["shim_t"]                     # aft face on the shim (62.5)
    x_collar0 = x_end - 8.0                         # 54.5
    x_base = x_collar0 - D["claw_b_floor_t"]        # 53.3
    f_ = CP["float"]
    x_a_tip = x_base - f_                           # claw A dog tips (52.3)
    x_b_tip = 46.5 + D["claw_a_flange_t"] + f_      # claw B dog tips (49.5)
    ro, rs = CP["od"] / 2, CP["claw_b_step_od"] / 2
    cx0, cx1 = D["claw_b_cone_x"]
    body = rev([(x_base, 4.2), (x_base, ro - 0.5), (x_base + 0.5, ro), (cx0, ro), (cx1, rs),
                (x_end - 0.4, rs), (x_end, rs - 0.4), (x_end, 4.2)])
    pocket = rev([(x_collar0, 0), (x_collar0, D["collar_pocket_d"] / 2), (x_end - 0.5, D["collar_pocket_d"] / 2),
                  (x_end + 0.01, D["collar_pocket_d"] / 2 + 0.51), (x_end + 0.01, 0)])
    body = body - pocket
    # crush ribs: 3 half-round ribs reaching D 15.9 (PRESS_FIT -0.05 per side on the collar OD 16)
    r_rib = D["crush_rib_r"]
    rc = 8.0 + P.PRESS_FIT + r_rib                  # rib inner extent = 8 - 0.05
    for k in range(3):
        th = 60.0 + 120.0 * k
        u_, v_ = uv(rc, th)
        rib = xcyl(x_collar0, x_end - 1.1, r_rib, th, rc).fuse(
            Solid.make_cone(r_rib, 0.02, 0.5, Plane(origin=(x_end - 1.1, u_, v_), x_dir=(0, 1, 0), z_dir=(1, 0, 0))))
        body = body.fuse(rib & tube(x_collar0, x_end, 7.0, D["collar_pocket_d"] / 2 + 0.01))
    body = body.clean()
    # set-screw access (print: aft face down, so 'down' is +X)
    body = body - radial_teardrop((x_collar0 + x_end) / 2, D["setscrew_theta"], D["setscrew_access_d"],
                                  7.0, ro + 1.0, down=(1, 0, 0))
    for th in (90.0, 270.0):
        body = body.fuse(dog(x_b_tip, x_base + 0.01, D["claw_dog_r_in"], ro, th, D["claw_dog_half_deg"],
                             D["claw_flank_gap"] / 2, x_b_tip, D["claw_tip_chamfer"]))
    return body.clean()


def claw_a() -> Solid:
    """Printed claw A on the generator rotor face: flange + two dogs carrying the screw counterbores."""
    D = DESIGN
    x0 = G["x_rotor_face"]                           # 46.5
    xf = x0 + D["claw_a_flange_t"]                   # 48.5
    x_tip = 63.0 - CP["shim_t"] - 8.0 - D["claw_b_floor_t"] - CP["float"]   # 52.3
    ro = D["claw_a_od"] / 2
    body = rev([(x0, 4.2), (x0, ro - 0.4), (x0 + 0.4, ro), (xf - 0.5, ro), (xf, ro - 0.5), (xf, 4.2)])
    for th in (0.0, 180.0):
        d = dog(xf - 0.01, x_tip, D["claw_dog_r_in"], ro - 0.5, th, D["claw_dog_half_deg"],
                D["claw_flank_gap"] / 2, x_tip, D["claw_tip_chamfer"])
        body = body.fuse(d)
    body = body.clean()
    r_pcd = G["rotor_face_pcd"] / 2
    for th in (0.0, 180.0):
        body = body - xcyl(x0 - 0.1, xf, 1.7, th, r_pcd)
        body = body - xcyl(xf, x_tip + 0.1, 3.0, th, r_pcd)
    return body


def claw_a_screws() -> Compound:
    """2 x ISO 4762 M3 x 4 (clamp 2.0 mm flange, 2.0 mm into the rotor-face thread / insert)."""
    xf = G["x_rotor_face"] + DESIGN["claw_a_flange_t"]
    r = G["rotor_face_pcd"] / 2
    return Compound([m3_socket_screw(xf, 4.0, th, r) for th in (0.0, 180.0)])


# --------------------------------------------------------------------------
# Bearings, seal, O-ring
# --------------------------------------------------------------------------


def _bearing(x0: float, d: float, D: float, B: float, r_ir: float, r_or: float, rec: float, seal_t: float,
             gap: float):
    """Steel rings (one solid: the seal grooves are cut between the rings) + the two rubber
    seals (compound of 2) recessed rec below each face, a gap ring next to the inner race."""
    x1 = x0 + B
    steel = tube(x0, x1, d / 2, D / 2, ch=0.3)
    g = rec + seal_t
    for xa, xb in ((x0 - 0.1, x0 + g), (x1 - g, x1 + 0.1)):
        steel = steel - tube(xa, xb, r_ir, r_or)
    # inner-race land: a small chamfer on the inner ring shoulder (reads as two rings)
    seals = Compound([tube(x0 + rec, x0 + g, r_ir + gap, r_or), tube(x1 - g, x1 - rec, r_ir + gap, r_or)])
    return steel, seals


def bearing_608():
    return _bearing(BT["x0"], BT["d"], BT["D"], BT["B"], 6.05, 9.6, 0.25, 0.35, 0.25)


def bearing_688():
    return _bearing(B8["x0"], B8["d"], B8["D"], B8["B"], 5.0, 7.0, 0.2, 0.3, 0.2)


def lip_seal():
    """TC 8x16x5 NBR double-lip seal, main lip facing aft (water side), at x 102.5 - 107.5.
    OD modelled at 15.98 and lips at r 4.005 (no overlap with the bore and the shaft)."""
    x0 = LS["x0"]
    ro = LS["D"] / 2 - 0.01
    s = [(0.0, ro), (5.0, ro), (5.0, 7.05), (1.6, 7.05), (1.6, 6.35), (4.1, 6.05), (4.6, 5.6),
         (4.6, 5.0), (3.7, 4.005), (1.9, 4.45), (1.0, 4.6), (0.5, 4.005), (0.0, 4.7)]
    body = rev([(x0 + a, b) for a, b in s])
    groove = Solid.make_torus(6.25, 0.55, Plane(origin=(x0 + 3.45, 0, 0), x_dir=(0, 1, 0), z_dir=(1, 0, 0)))
    body = body - groove
    spring = Solid.make_torus(6.25, 0.42, Plane(origin=(x0 + 3.45, 0, 0), x_dir=(0, 1, 0), z_dir=(1, 0, 0)))
    return body, spring


def oring() -> Solid:
    """ISO 3601 / AS568 -141 as installed: stretched onto the groove bottom (D 59.88) and
    squeezed into the gland (2.105 mm): elliptical section of the same area as the
    stretched round section, against the radial aft flank."""
    cs = OR["cs"] * math.sqrt(OR["id"] / (OR["groove_d"]))     # stretched section (volume kept)
    area = math.pi / 4 * cs ** 2
    b = (OR["gland_depth"] - 0.04) / 2                         # radial semi-axis, 0.02 air each side so
    a = area / (math.pi * b)                                   # faceted meshes never clash; axial semi-axis
    r_c = OR["groove_d"] / 2 + 0.02 + b
    x_c = J["groove_x"]["x_aft"] - a - 0.02
    face = Plane.XZ * Ellipse(a, b).moved(Location((x_c, r_c)))
    return revolve(face, Axis.X), dict(cs_stretched=cs, a=a, b=b, x_c=x_c, r_c=r_c)


# --------------------------------------------------------------------------
# Carrier, spacer, GL40 replica
# --------------------------------------------------------------------------


def _rounded_window(pl: Plane, th: float, r0: float, r1: float, hw: float, rc: float, depth: float):
    """Radial window at azimuth th: from r0 to r1, half-width hw, corners rounded rc, extruded
    depth along the plane normal (+X)."""
    rect = Rectangle(2 * hw, r1 - r0, align=(Align.CENTER, Align.MIN))
    from build123d import fillet as _f
    sk = _f(rect.vertices(), rc) if rc > 0 else rect
    sk = sk.moved(Location((0, r0))).rotate(Axis.Z, -th)   # local v = +Z (theta 0); rotate toward +Y
    return Solid.extrude((pl * sk).faces()[0], pl.z_dir * depth)


def carrier() -> Solid:
    D = DESIGN
    x0, t = CR["x0"], CR["t"]
    x1 = x0 + t
    R = CR["r"]
    plate = tube(x0, x1, CR["centre_bore_d"] / 2, R, ch=0.0)
    plate = rev([(x0, CR["centre_bore_d"] / 2 + 0.4), (x0 + 0.4, CR["centre_bore_d"] / 2),
                 (x1, CR["centre_bore_d"] / 2), (x1, R - 0.5), (x1 - 0.5, R), (x0 + 0.3, R), (x0, R - 0.3)])
    L = carrier_lock()
    t0, t1 = L["tab_x"]
    Rt, h = L["tab_r_out"], L["tab_front_chamfer_45"]
    hw = L["tab_w_deg"] / 2
    # tab section in (x, r): front-outer corner chamfered 45 deg (seats on the groove's 45 deg ceiling)
    sec = [(t0, R - 1.0), (t0, Rt - h), (t0 + h, Rt), (t1, Rt), (t1, R - 1.0)]
    for th in L["lock_theta_deg"]:
        face = Plane.XZ * Polygon(*sec, align=None)
        # revolve about +X by +angle turns +Z toward -Y (theta 0 .. -2 hw); turn back by hw so the
        # sector spans theta -hw .. +hw, then carry it to th
        tab = revolve(face, Axis.X, 2 * hw).rotate(Axis.X, -hw)
        plate = plate.fuse(at_azimuth(tab, th))
    plate = plate.clean()
    # GL40 pattern: 4 x M3 on PCD 24 (counterbored front face) + radial slots PCD 38.14-40
    r24 = G["mount_pcd"] / 2
    cb_d, cb_h = D["carrier_cbore"]
    s0, s1, sw = D["carrier_slot"]
    for th in THETA_SCREWS:
        plate = plate - xcyl(x0 - 0.1, x1 + 0.1, D["carrier_clear_d"] / 2, th, r24)
        plate = plate - xcyl(x0 - 0.1, x0 + cb_h, cb_d / 2, th, r24)
        slot = Rectangle(sw, s1 / 2 - s0 / 2 + sw, align=(Align.CENTER, Align.MIN))
        from build123d import fillet as _f
        slot = _f(slot.vertices(), sw / 2 - 0.01).moved(Location((0, s0 / 2 - sw / 2))).rotate(Axis.Z, -th)
        plate = plate - Solid.extrude((PX(x0 - 0.1) * slot).faces()[0], (t + 0.2, 0, 0))
    # solder-pad keep-out window (hardware.json cad_window) merged with the cable notch at theta 180
    w = PAD_WINDOW
    plate = plate - _rounded_window(PX(x0 - 0.1), THETA_PAD, w["radial_from_mm"], w["radial_to_mm"],
                                    w["tangential_half_width_mm"], D["carrier_window_corner_r"], t + 0.2)
    cn = CR["cable_notch"]
    plate = plate - _rounded_window(PX(x0 - 0.1), cn["theta_deg"], cn["bottom_r"], R + 2.0, cn["width"] / 2,
                                    0.0, t + 0.2)
    # assembly aid on the front face, readable from the nose: "LOCK" and a chevron pointing the twist
    # direction (-theta = clockwise seen from the nose); engraved 0.4 deep (the face prints on top)
    from build123d import FontStyle, Text
    front = Plane(origin=(x0 - 0.1, 0, 0), x_dir=(0, -1, 0), z_dir=(-1, 0, 0))   # local x = -Y, y = +Z
    txt = Text("LOCK", 2.6, font="DejaVu Sans", font_style=FontStyle.BOLD).moved(Location((0, 17.2)))
    plate = plate - Compound([Solid.extrude(f, Vector(0.5, 0, 0)) for f in (front * txt).faces()])
    r_a, hw_a = 17.4, 0.45
    arc = lambda r, t: (r * math.sin(math.radians(t)), r * math.cos(math.radians(t)))   # noqa: E731
    shaft_pts = [arc(r_a + hw_a, t) for t in np.linspace(-15.0, -24.0, 10)] + \
                [arc(r_a - hw_a, t) for t in np.linspace(-24.0, -15.0, 10)]
    head_pts = [arc(r_a, -31.0), arc(r_a + 1.3, -23.5), arc(r_a - 1.3, -23.5)]
    arrow = (PX(x0 - 0.1) * Polygon(*shaft_pts, align=None)).faces()[0]
    head = (PX(x0 - 0.1) * Polygon(*head_pts, align=None)).faces()[0]
    plate = plate - Solid.extrude(head, (0.5, 0, 0))
    plate = plate - Solid.extrude(arrow, (0.5, 0, 0))
    return plate


def carrier_screws() -> Compound:
    """4 x ISO 7380 M3 x 8 button head: head in the 1.7 mm counterbore (0.05 below the face),
    through the carrier (2.3 under the head) and the spacer (2.5) into the stator: 3.2 thread."""
    x_head = CR["x0"] + DESIGN["carrier_cbore"][1]
    r24 = G["mount_pcd"] / 2
    return Compound([m3_button_screw(x_head, 8.0, th, r24) for th in THETA_SCREWS])


def adapter_spacer() -> Solid:
    D = DESIGN
    x0 = CR["x0"] + CR["t"]
    x1 = G["motor_body_x"][0]
    ro = D["spacer_od"] / 2
    sp = rev([(x0, 8.0), (x1, 8.0), (x1, ro - 0.4), (x1 - 0.4, ro), (x0 + 0.3, ro), (x0, ro - 0.3)])
    r24 = G["mount_pcd"] / 2
    for th in THETA_SCREWS:
        sp = sp - xcyl(x0 - 0.1, x1 + 0.1, D["carrier_clear_d"] / 2, th, r24)
        sp = sp - xcyl(x0 - 0.1, x1 + 0.1, 1.45, th, 20.0)          # GL40 M2.5 on PCD 40 (optional)
    w = PAD_WINDOW
    sp = sp - _rounded_window(PX(x0 - 0.1), THETA_PAD, w["radial_from_mm"], ro + 2.0,
                              w["tangential_half_width_mm"], 0.0, x1 - x0 + 0.2)
    return sp


def gen_stator() -> Solid:
    """Static half of the printed GL40 stand-in: flange (front face x 25 = the real stator face),
    stator stack with 24 slots, journal post for the bell; M3 insert holes on PCD 24 and the
    PCD 40 M2.5 holes as on the real motor; the solder-pad block is engraved (not raised) so the
    face prints flat on the bed."""
    D = DESIGN
    xf = G["motor_body_x"][0]                        # 25.0
    fd, ft = D["stator_base"]
    core_od, core_id = D["stator_core"]
    xb = xf + ft
    base = rev([(xf, 4.0), (xf, fd / 2 - 0.6), (xf + 0.6, fd / 2), (xb, fd / 2), (xb, 4.0)])
    core = tube(xb - 0.01, D["stator_core_x1"], core_id / 2, core_od / 2)
    core = rev([(xb - 0.01, core_id / 2), (xb - 0.01, core_od / 2), (D["stator_core_x1"] - 0.5, core_od / 2),
                (D["stator_core_x1"], core_od / 2 - 0.5), (D["stator_core_x1"], core_id / 2)])
    slot_w, slot_d = 1.4, 3.0
    for k in range(D["stator_slots"]):
        th = (k + 0.5) * 360.0 / D["stator_slots"]
        sl = Box(D["stator_core_x1"] - xb + 1.0, slot_w, slot_d + 1.0,
                 align=(Align.MIN, Align.CENTER, Align.MIN)).moved(Location((xb + 0.6, 0, core_od / 2 - slot_d)))
        core = core - at_azimuth(sl, th)
    pd0, pd1 = D["post_d"]
    post = rev([(xb - 0.01, 4.0), (xb - 0.01, pd0 / 2), (D["post_shoulder_x"], pd0 / 2),
                (D["post_shoulder_x"], pd1 / 2), (D["post_top_x"] - 0.4, pd1 / 2), (D["post_top_x"], pd1 / 2 - 0.4),
                (D["post_top_x"], 4.0)])
    st = base.fuse(core).fuse(post).clean()
    # inserts (PCD 24) and the M2.5 PCD 40 holes on the front face, all at theta 45 + 90 k
    r24 = G["mount_pcd"] / 2
    for th in THETA_SCREWS:
        st = st - blind_hole(xf, xf + D["insert_hole_depth"] + 0.5, P.INSERT_M3["hole_d"] / 2, th, r24)
        st = st - rev([(xf - 0.1, 0), (xf - 0.1, 2.4), (xf + 0.4, 2.0), (xf + 0.4, 0)]).moved(Location((0, *uv(r24, th))))
        st = st - blind_hole(xf, xf + 3.5, 1.25, th, 20.0)
    # engraved solder-pad block outline and its three pads (0.4 deep) at theta 180
    pb = PAD_BLOCK
    outer = _rounded_window(PX(xf - 0.1), THETA_PAD, pb["r0"], pb["r1"], pb["hw"], 1.0, 0.5)
    inner = _rounded_window(PX(xf - 0.2), THETA_PAD, pb["r0"] + 0.5, pb["r1"] - 0.5, pb["hw"] - 0.5, 0.6, 0.8)
    st = st - (outer - inner)
    for tc in (-3.0, 0.0, 3.0):
        pad = Rectangle(2.0, 2.0, align=(Align.CENTER, Align.MIN)).moved(Location((tc, 16.6)))
        pad = pad.rotate(Axis.Z, -THETA_PAD)
        st = st - Solid.extrude((PX(xf - 0.1) * pad).faces()[0], (0.5, 0, 0))
    return st


def gen_bell() -> Solid:
    """Turning half of the printed GL40 stand-in: the outer rotor bell (D 46.5, face at x 46.5)
    with its hub sleeve on the stator journal, 4 M3 insert bosses on PCD 24 (claw A uses 0 / 180)."""
    D = DESIGN
    xr = G["x_rotor_face"]                            # 46.5
    ro = G["D"] / 2
    ri = ro - D["bell_wall"]
    xface = xr - D["bell_face_t"]                     # 44.0
    x_rim = G["motor_body_x"][0] + D["stator_base"][1] + D["bell_rim_gap"]      # 28.5
    cup = rev([(x_rim, ri), (x_rim, ro - 0.3), (x_rim + 0.3, ro), (xr - 1.0, ro), (xr, ro - 1.0),
               (xr, 4.0), (xface, 4.0), (xface, ri)])
    # decorative magnet-ring line on the can
    xg = x_rim + 3.0                                  # V groove, faces 37 deg from vertical when printed
    cup = cup - rev([(xg, ro + 0.1), (xg + 0.4, ro - 0.3), (xg + 0.8, ro + 0.1)])
    jb = D["post_d"][1] / 2 + P.CLEAR_ROTATE         # journal bore radius
    sl = rev([(D["bell_sleeve_x0"], jb), (D["bell_sleeve_x0"], D["bell_sleeve_od"] / 2 - 0.4),
              (D["bell_sleeve_x0"] + 0.4, D["bell_sleeve_od"] / 2), (xface + 0.01, D["bell_sleeve_od"] / 2),
              (xface + 0.01, jb)])
    bell = cup.fuse(sl)
    r24 = G["rotor_face_pcd"] / 2
    for th in (0.0, 90.0, 180.0, 270.0):
        boss = xcyl(D["bell_boss_x0"], xface + 0.01, D["bell_boss_d"] / 2, th, r24)
        bell = bell.fuse(boss)
    bell = bell.clean()
    # journal bore up to a 1.2 mm skin under the face, D 8 cosmetic centre hole (real 8 mm bore)
    bell = bell - xcyl(D["bell_sleeve_x0"] - 0.1, xr - P.MIN_WALL, jb)
    bell = bell - xcyl(xr - P.MIN_WALL - 0.1, xr + 0.1, 4.0)
    for th in (0.0, 90.0, 180.0, 270.0):
        bell = bell - blind_hole(xr, xr - D["insert_hole_depth"], P.INSERT_M3["hole_d"] / 2, th, r24)
        bell = bell - rev([(xr + 0.1, 0), (xr + 0.1, 2.4), (xr - 0.4, 2.0), (xr - 0.4, 0)]).moved(
            Location((0, *uv(r24, th))))
    return bell


def gen_inserts() -> tuple[Compound, Compound]:
    """Heat-set inserts: 4 in the stator face (carrier screws), 2 in the bell (claw A screws)."""
    xf = G["motor_body_x"][0]
    r24 = G["mount_pcd"] / 2
    st = Compound([heat_set_insert(xf, th, r24, +1) for th in THETA_SCREWS])
    xr = G["x_rotor_face"]
    bl = Compound([heat_set_insert(xr, th, G["rotor_face_pcd"] / 2, -1) for th in (0.0, 180.0)])
    return st, bl


# --------------------------------------------------------------------------
# Cable gland and tools
# --------------------------------------------------------------------------


def gland_frame() -> tuple[Vector, Vector]:
    """Spot-face centre and the outward gland axis (unit) in world coordinates."""
    g = P.GLAND
    t = math.radians(g["theta_deg"])
    er = Vector(0, math.sin(t), math.cos(t))
    c = Vector(g["spot_face_centre_x"], 0, 0) + er * g["spot_face_centre_r"]
    a = math.radians(g["axis_angle_from_x_deg"])
    u = Vector(math.cos(a), 0, 0) + er * math.sin(a)
    return c, u.normalized()


def cable_gland() -> Solid:
    """Hummel HSK-K M12x1.5 (1.209.1200.50), simplified: thread core at the M12 minor diameter
    (10.4, no overlap with the 10.5 tap-drill bore) 8 mm into the seat, base hex SW 15, cap nut
    SW 15 with a dome, cable bore 6.5 (clamping 3 - 6.5). Height 21 above the seat face."""
    g = P.GLAND
    L = g["thread_length"]
    # local profile along +s (outward), built on the X axis then placed
    prof = [(-L, 0), (-L, 4.7), (-L + 0.5, 5.2), (0.0, 5.2), (0.0, 0)]
    core = rev(prof)
    hx = lambda s0, s1: Solid.extrude((PX(s0) * RegularPolygon(15.0 / math.sqrt(3), 6)).faces()[0], (s1 - s0, 0, 0))
    base = hx(0.0, 5.0) & rev([(0.0, 0), (0.0, 8.4), (4.4, 8.4), (5.0, 7.6), (5.0, 0)])
    neck = rev([(4.99, 0), (4.99, 6.2), (6.5, 6.2), (6.5, 0)])
    cap = hx(6.5, 12.5) & rev([(6.5, 0), (6.5, 7.6), (7.1, 8.4), (12.5, 8.4), (12.5, 0)])
    dome = rev([(12.49, 0), (12.49, 7.4), (14.5, 7.0), (17.0, 6.0), (19.2, 4.8), (21.0, 3.6), (21.0, 0)])
    body = core.fuse(base).fuse(neck).fuse(cap).fuse(dome).clean()
    body = body - rev([(-L - 0.1, 0), (-L - 0.1, 3.25), (21.1, 3.25), (21.1, 0)])
    c, u = gland_frame()
    return _place_on_axis(body, c, u)


def _place_on_axis(shape: Shape, origin: Vector, axis: Vector) -> Shape:
    """Move a shape built along +X (origin at 0) so +X maps to `axis` through `origin`."""
    ax = axis.normalized()
    ref = Vector(0, 1, 0) if abs(ax.Y) < 0.9 else Vector(0, 0, 1)
    z = ax.cross(ref).normalized()
    pl = Plane(origin=origin, x_dir=ax, z_dir=z)
    return shape.moved(pl.location)


def seal_sleeve() -> Solid:
    """Installation sleeve for the lip seal (v2): K&S brass tube 9 x 0.45 (ID 8.1) x 25,
    chamfered lead-in. A tool, shown parked aft of the spinner (x 175)."""
    x0 = 175.0
    return rev([(x0, 4.05), (x0 + 25, 4.05), (x0 + 25, 4.5), (x0 + 1.2, 4.5), (x0, 4.15)])


# --------------------------------------------------------------------------
# Part dicts
# --------------------------------------------------------------------------


def _bom(name, shape, print_shape, material, role, rotor, qty, bought, explode, notes, **extra):
    d = dict(name=name, shape=shape, print=print_shape, material=material, role=role, rotor=rotor,
             qty=qty, bought=bought, explode=explode, notes=notes)
    d.update(extra)
    return d


_CACHE: dict = {}


def parts() -> list[dict]:
    if "parts" in _CACHE:
        return _CACHE["parts"]
    out = []
    sh = shaft()
    out.append(_bom("shaft", sh, None, "steel", "internal", True, 1, True, [62, 0, 0],
                    "8 x 75.5 stainless ground rod, MISUMI PSSFJ8 class (SUS440C, 56 HRC, g6), ordered "
                    "machined: hub D-flat 0.5 deep x 108.4-128 (grub side, psi 60), collar flat 0.5 deep "
                    "x 55-62 (theta 0), M3 x 8 tapped aft end, 0.5 chamfers. No ring groove (collar is the stop)."))
    out.append(_bom("shaft_collar", shaft_collar(), None, "steel", "internal", True, 1, True, [8, 0, 0],
                    "Ruland MSC-8-SS set-screw collar 8 x 16 x 8, 303 stainless; thrust stop on the 608 inner "
                    "race through the shim; pressed into claw B. Set-screw size not confirmed (M4 assumed, VERIFY)."))
    out.append(_bom("shim_8x12x0p5", shim(), None, "steel", "internal", True, 1, True, [14, 0, 0],
                    "stainless shim 8.1 x 12 x 0.5 on the 608 inner race (clears the seals). DIN 988 lists 8 x 14; "
                    "an 8 x 12 precision shim or a DIN 988 S 8x14x0.5 both work (the 608 seals are recessed)."))
    cb = claw_b()
    out.append(_bom("claw_b", cb, on_bed(cb, (0, 90, 0)), "hub", "internal", True, 1, False, [-2, 0, 0],
                    "PETG, aft face on the bed (dogs up), no supports: 45 deg OD cone at the step, teardrop "
                    "set-screw access, the floor bridges the collar pocket (D 16.2 annulus). Press the collar in "
                    "(3 crush ribs to D 15.9) with its set screw under the access hole; a drop of CA locks it."))
    ca = claw_a()
    out.append(_bom("claw_a", ca, on_bed(ca, (0, -90, 0)), "hub", "internal", True, 1, False, [-12, 0, 0],
                    "PETG, flange on the bed (rotor-face side down), dogs up, no supports. 2 x M3 x 4 ISO 4762 "
                    "into the rotor-face holes at 0 / 180 (heads in the dog counterbores)."))
    out.append(_bom("screws_claw_a_m3x4", claw_a_screws(), None, "steel", "internal", True, 2, True, [-12, 0, 0],
                    "2 x ISO 4762 M3 x 4 A2 (2.0 mm thread engagement; the GL40 rotor holes are 2.5 deep)."))
    b6s, b6r = bearing_608()
    out.append(_bom("bearing_608", b6s, None, "steel", "internal", False, 1, True, [20, 0, 0],
                    "608-2RS stainless (S608-2RS, AISI 440C) 8 x 22 x 7, x 63-70, locating: outer race on the "
                    "18 mm tube shoulder, inner race thrust from the collar through the shim."))
    out.append(_bom("bearing_608_seals", b6r, None, "rubber", "internal", False, 1, True, [20, 0, 0],
                    "the two 2RS seals of bearing_608 (part of it, not a separate BOM line).",
                    bom_parent="bearing_608"))
    b8s, b8r = bearing_688()
    out.append(_bom("bearing_688", b8s, None, "steel", "internal", False, 1, True, [40, 0, 0],
                    "688-2RS 8 x 16 x 5, x 96-101, floating (slip fit on the shaft), outer race on the seal shoulder."))
    out.append(_bom("bearing_688_seals", b8r, None, "rubber", "internal", False, 1, True, [40, 0, 0],
                    "the two 2RS seals of bearing_688 (part of it).", bom_parent="bearing_688"))
    ls, spring = lip_seal()
    out.append(_bom("lip_seal_tc8x16x5", ls, None, "rubber", "internal", False, 1, True, [48, 0, 0],
                    "TC 8x16x5 NBR double lip (v2), x 102.5-107.5, main lip aft; installed after the shaft from "
                    "the tail end over the brass sleeve (tool_seal_sleeve)."))
    out.append(_bom("lip_seal_spring", spring, None, "steel", "internal", False, 1, True, [48, 0, 0],
                    "garter spring of the TC seal (part of it).", bom_parent="lip_seal_tc8x16x5"))
    orr, _ = oring()
    out.append(_bom("oring_141", orr, None, "rubber", "internal", False, 1, True, [30, 0, 0],
                    "ISO 3601 / AS568 -141 NBR 70 (ID 58.42 x 2.62), modelled as installed: stretched onto the "
                    "groove bottom D 59.88 and squeezed into the 2.105 gland (elliptic section, same area)."))
    car = carrier()
    out.append(_bom("motor_carrier", car, on_bed(car, (0, 90, 0)), "hub", "internal", False, 1, False,
                    [-56, 0, 0],
                    "PETG, AFT face on the bed (tabs flush with it, counterbores open upward), no supports; the "
                    "tabs' 45 deg front chamfer faces up. Twist-lock per "
                    "shells.CARRIER_LOCK: 3 tabs 16 deg wide enter at theta 120 / 240 / 0 and lock at 90 / 210 / "
                    "330 after a 30 deg twist toward -theta (the generator reaction torque drives them into the stop)."))
    out.append(_bom("screws_carrier_m3x8", carrier_screws(), None, "steel", "internal", False, 4, True,
                    [-66, 0, 0], "4 x ISO 7380 M3 x 8 A2 button head (3.2 mm into the motor / insert)."))
    sp = adapter_spacer()
    out.append(_bom("adapter_spacer", sp, on_bed(sp, (0, -90, 0)), "hub", "internal", False, 1, False,
                    [-46, 0, 0], "PETG, flat, no supports. GL40-specific: PCD 24 M3 + PCD 40 M2.5 clearance, "
                    "solder-pad window open to the rim at theta 180 (leads run out to the carrier notch)."))
    stt = gen_stator()
    out.append(_bom("generator_stator", stt, on_bed(stt, (0, -90, 0)), "generator", "internal", False, 1,
                    False, [-36, 0, 0],
                    "PLA display stand-in for the CubeMars GL40 stator (real envelope with the bell). Front face "
                    "on the bed, no supports. 4 M3 heat-set inserts on PCD 24 take the carrier screws."))
    bel = gen_bell()
    out.append(_bom("generator_bell", bel, on_bed(bel, (0, 90, 0)), "generator", "internal", True, 1, False,
                    [-22, 0, 0],
                    "PLA display stand-in for the GL40 rotor bell; turns on the stator journal (0.3 mm per side). "
                    "Rotor face on the bed, cup up, no supports. 2 M3 inserts take claw A."))
    ins_s, ins_b = gen_inserts()
    out.append(_bom("inserts_m3_stator", ins_s, None, "magnet", "internal", False, 4, True, [-36, 0, 0],
                    "4 x CNC Kitchen M3 x 5.7 heat-set inserts in the replica stator face (take the carrier "
                    "screws); display build only, the real GL40 has M3 threads. Modelled at the hole size."))
    out.append(_bom("inserts_m3_bell", ins_b, None, "magnet", "internal", True, 2, True, [-22, 0, 0],
                    "2 x CNC Kitchen M3 x 5.7 heat-set inserts in the replica bell (take claw A); display build "
                    "only. Modelled at the hole size."))
    gl = cable_gland()
    out.append(_bom("cable_gland_m12", gl, None, "default", "body", False, 1, True, [-90, 0, -14],
                    "Hummel HSK-K M12x1.5 PA grey (1.209.1200.50), thread in the tapped nose seat (v2)."))
    out.append(_bom("tool_seal_sleeve", seal_sleeve(), None, "steel", "hidden", False, 1, True, [0, 0, 0],
                    "assembly tool: K&S brass tube 9 x 0.45 mm (ID 8.1) x 25 with a filed lead-in; slid over the "
                    "protruding shaft end while the lip seal is pushed home. Parked aft of the model."))
    _CACHE["parts"] = out
    return out


def preview_parts() -> list[dict]:
    """Every drivetrain part with role 'body' so the viewer shows them without the shells."""
    res = []
    for p in parts():
        if p["role"] == "hidden":
            continue
        q = {k: v for k, v in p.items() if k != "print"}
        q["role"] = "body"
        res.append(q)
    return res


def preview_coupling() -> list[dict]:
    """Close-up: claw A, claw B, collar, shim, 608 and the shaft front, pulled apart along X."""
    keep = {"claw_a": -10.0, "screws_claw_a_m3x4": -10.0, "generator_bell": -16.0, "claw_b": 6.0,
            "shaft_collar": 18.0, "shim_8x12x0p5": 22.0, "bearing_608": 26.0, "bearing_608_seals": 26.0}
    res = []
    for p in parts():
        if p["name"] in keep:
            res.append(dict(name=p["name"], shape=p["shape"].moved(Location((keep[p["name"]], 0, 0))),
                            material=p["material"], role="body"))
    return res


def preview_generator() -> list[dict]:
    """Close-up: carrier, spacer, stator, bell pulled apart along X (shows faces and holes)."""
    keep = {"motor_carrier": -24.0, "screws_carrier_m3x8": -34.0, "adapter_spacer": -12.0, "generator_stator": 0.0,
            "inserts_m3_stator": 0.0, "inserts_m3_bell": 22.0, "generator_bell": 22.0, "claw_a": 34.0}
    res = []
    for p in parts():
        if p["name"] in keep:
            res.append(dict(name=p["name"], shape=p["shape"].moved(Location((keep[p["name"]], 0, 0))),
                            material=p["material"], role="body"))
    return res


def _latest_shells(names=("front_shell_cutaway", "aft_shell_cutaway")) -> list[dict]:
    """Shell solids from the newest shells BREP cache (display only, no rebuild)."""
    import glob
    import os
    from build123d import import_brep
    dirs = glob.glob(str(ROOT / "out/work/cad/shells/cache/*"))
    if not dirs:
        return []
    d = max(dirs, key=os.path.getmtime)
    return [dict(name=n, shape=import_brep(f"{d}/{n}.brep"), material="shell", role="body")
            for n in names if Path(f"{d}/{n}.brep").exists()]


def preview_in_shells() -> list[dict]:
    """Drivetrain (and the fins) inside the shells' cutaway variants, for a context check."""
    out = preview_parts() + _latest_shells()
    try:
        from .fins import preview_parts as fin_parts
        out += fin_parts()
    except Exception:  # noqa: BLE001
        pass
    return out


def preview_carrier() -> list[dict]:
    """The carrier alone (front face toward the camera with view=front)."""
    return [dict(name="motor_carrier", shape=carrier(), material="hub", role="body")]


def preview_print() -> list[dict]:
    """Printed parts in their print orientation, laid out side by side on a virtual bed."""
    res, x = [], 0.0
    for p in parts():
        if p["print"] is None:
            continue
        s = p["print"]
        bb = s.bounding_box()
        res.append(dict(name=p["name"] + "_print", shape=s.moved(Location((x - bb.min.X, -bb.center().Y, 0))),
                        material=p["material"], role="body"))
        x += bb.size.X + 8.0
    return res


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def _mesh_check(shape, tol=0.03):
    import trimesh
    v, f = shape.tessellate(tol, 0.2)
    m = trimesh.Trimesh(np.array([(p.X, p.Y, p.Z) for p in v]), np.array(f), process=True)
    m.merge_vertices()
    return m


def _overhang_area(m, limit_deg=46.0, z_tol=0.05, bridges: bool = False):
    """Area (mm2) of faces off the bed that overhang more than limit_deg from vertical (46 = the
    45 deg rule with 1 deg for tessellation of exact 45 deg chamfers). bridges=True returns only the
    flat downward faces (normal within 10 deg of straight down), which print as bridges."""
    n = m.face_normals
    zc = m.triangles_center[:, 2]
    lim = -math.cos(math.radians(10.0)) if bridges else -math.sin(math.radians(limit_deg))
    down = n[:, 2] < lim
    off_bed = zc > z_tol
    return float(m.area_faces[down & off_bed].sum())


def checks(parts_list=None) -> dict:
    ps = parts_list or parts()
    by = {p["name"]: p for p in ps}
    res = {"parts": [], "interfaces": [], "assembly": []}
    for p in ps:
        c = check(p["shape"], p["name"])
        c.update(bought=p["bought"], material=p["material"], role=p["role"])
        if p["print"] is not None:
            m = _mesh_check(p["print"])
            c.update(print_size_mm=[round(v, 2) for v in p["print"].bounding_box().size],
                     print_min_z=round(p["print"].bounding_box().min.Z, 4),
                     watertight=bool(m.is_watertight), fits_bed=fits_bed(p["print"]),
                     overhang_area_mm2_gt46=round(_overhang_area(m), 2),
                     of_which_bridges_mm2=round(_overhang_area(m, bridges=True), 2))
        else:
            m = _mesh_check(p["shape"])
            c.update(watertight=bool(m.is_watertight))
        res["parts"].append(c)
    I = res["interfaces"]

    def add(name, value, ok, **kw):
        I.append(dict(check=name, value=value, ok=bool(ok), **kw))

    ca, cb = by["claw_a"]["shape"], by["claw_b"]["shape"]
    add("claw A vs claw B overlap (mm3)", round(interference(ca, cb), 4), interference(ca, cb) < 1e-6)
    d_ab = min_distance(ca, cb)
    add("claw A to claw B min gap (dog flanks, mm)", round(d_ab, 3), abs(d_ab - DESIGN["claw_flank_gap"]) < 0.02)
    bbA, bbB = ca.bounding_box(), cb.bounding_box()
    x_b_base = 63.0 - CP["shim_t"] - 8.0 - DESIGN["claw_b_floor_t"]
    add("axial float: claw A dog tips to claw B base (mm)", round(x_b_base - bbA.max.X, 3),
        abs(x_b_base - bbA.max.X - CP["float"]) < 1e-6)
    add("axial float: claw B dog tips to claw A flange (mm)",
        round(bbB.min.X - (G["x_rotor_face"] + DESIGN["claw_a_flange_t"]), 3),
        abs(bbB.min.X - (G["x_rotor_face"] + DESIGN["claw_a_flange_t"]) - CP["float"]) < 1e-6)
    add("claw engagement depth (mm)", round(bbA.max.X - bbB.min.X, 3), bbA.max.X - bbB.min.X >= 2.5)
    add("claw B aft face + shim = 608 front face (x)", round(bbB.max.X + CP["shim_t"], 3),
        abs(bbB.max.X + CP["shim_t"] - BT["x0"]) < 1e-6)
    seat_r = (BT["D"] + 2 * P.PRESS_FIT) / 2
    add("claw B OD 20 step to the 608 seat bore in the tube mouth, radial (mm)",
        round(seat_r - CP["claw_b_step_od"] / 2, 3), seat_r - CP["claw_b_step_od"] / 2 >= 0.3)
    add("claw B cone ends ahead of the tube mouth (x)", DESIGN["claw_b_cone_x"][1],
        DESIGN["claw_b_cone_x"][1] <= P.BEARING_TUBE["x0"] - P.BEARING_TUBE["mouth_lead_in"])
    col = by["shaft_collar"]["shape"]
    add("collar in claw B: press overlap at the crush ribs (mm3, expected > 0)",
        round(interference(col, cb), 3), interference(col, cb) > 0)
    sh = by["shaft"]["shape"]
    bbS = sh.bounding_box()
    add("shaft x0 / x1", [round(bbS.min.X, 3), round(bbS.max.X, 3)],
        abs(bbS.min.X - S["x0"]) < 1e-6 and abs(bbS.max.X - S["x1"]) < 1e-6)
    add("shaft to claw A (mm)", round(min_distance(sh, ca), 3), min_distance(sh, ca) > 0.3)
    add("shaft vs claw B overlap (mm3)", round(interference(sh, cb), 4), interference(sh, cb) < 1e-6)
    for nm in ("bearing_608", "bearing_688", "lip_seal_tc8x16x5", "shaft_collar", "shim_8x12x0p5"):
        v = interference(sh, by[nm]["shape"])
        add(f"shaft vs {nm} overlap (mm3)", round(v, 4), v < 1e-3)
    ls = by["lip_seal_tc8x16x5"]["shape"].bounding_box()
    add("lip seal x / recess behind the static tail face (mm)", [round(ls.min.X, 3), round(ls.max.X, 3),
                                                                 round(P.STATIC_TAIL_X - ls.max.X, 3)],
        abs(P.STATIC_TAIL_X - ls.max.X - LS["recess"]) < 1e-6)
    add("lip seal outer face to hub front (mm)", round(P.ROTOR_HUB["x_front"] - ls.max.X, 3),
        P.ROTOR_HUB["x_front"] - ls.max.X >= P.HUB_GAP)
    car, spc = by["motor_carrier"]["shape"], by["adapter_spacer"]["shape"]
    stt, bel = by["generator_stator"]["shape"], by["generator_bell"]["shape"]
    bbc = car.bounding_box()
    add("carrier x / plate r", [round(bbc.min.X, 3), round(bbc.max.X, 3), CR["r"]],
        abs(bbc.min.X - CR["x0"]) < 1e-6 and abs(bbc.max.X - CR["x0"] - CR["t"]) < 1e-6)
    gen = stt.fuse(bel)
    bbg = gen.bounding_box()
    add("GL40 replica envelope D x L (x from-to)", [round(bbg.size.Y, 3), round(bbg.size.X, 3),
                                                    round(bbg.min.X, 3), round(bbg.max.X, 3)],
        abs(bbg.size.Y - G["D"]) < 0.01 and abs(bbg.min.X - G["motor_body_x"][0]) < 1e-6
        and abs(bbg.max.X - G["motor_body_x"][1]) < 1e-6)
    bay = Solid.make_cylinder(G["bay_D"] / 2, G["bay_L"], Plane(origin=(G["x_bay_fwd"], 0, 0), z_dir=(1, 0, 0)))
    for nm in ("adapter_spacer", "generator_stator", "generator_bell"):
        s = by[nm]["shape"]
        out_v = s.volume - interference(s, bay)
        add(f"{nm} inside the bay D50.5 x 24 (volume outside, mm3)", round(out_v, 3), out_v < 1e-3)
    for a_, b_ in (("motor_carrier", "adapter_spacer"), ("adapter_spacer", "generator_stator"),
                   ("generator_stator", "generator_bell"), ("generator_bell", "claw_a")):
        v = interference(by[a_]["shape"], by[b_]["shape"])
        add(f"{a_} vs {b_} overlap (mm3)", round(v, 4), v < 1e-3)
    add("bell journal: radial running clearance per side (mm)", P.CLEAR_ROTATE, True)
    add("stator to bell min gap (mm)", round(min_distance(stt, bel), 3), min_distance(stt, bel) >= 0.29)
    for nm in ("screws_carrier_m3x8", "inserts_m3_stator", "inserts_m3_bell", "screws_claw_a_m3x4"):
        for host in ("motor_carrier", "adapter_spacer", "generator_stator", "generator_bell", "claw_a"):
            v = interference(by[nm]["shape"], by[host]["shape"])
            if v > 1e-3:
                add(f"{nm} vs {host} overlap (mm3)", round(v, 4), False)
    # carrier window / notch / screw webs (2D distances at the plate)
    r24 = G["mount_pcd"] / 2
    w = PAD_WINDOW
    cbr = DESIGN["carrier_cbore"][0] / 2
    pts = []
    for th in (135.0, 225.0):
        p = np.array(uv(r24, th))
        # window: rounded rectangle at theta 180 (tangential |u| <= hw, radial r0..r1), corner r
        q = np.array([abs(p[0]), -p[1]])          # local (tangential, radial along theta 180)
        hw, r0, r1, rc = w["tangential_half_width_mm"], w["radial_from_mm"], w["radial_to_mm"], DESIGN[
            "carrier_window_corner_r"]
        cx = np.clip(q[0], -(hw - rc), hw - rc)
        cy = np.clip(q[1], r0 + rc, r1 - rc)
        d = math.hypot(q[0] - cx, q[1] - cy) - rc
        pts.append(d - cbr)
    add("carrier: counterbore (PCD 24) to the pad window web (mm)", round(min(pts), 3), min(pts) >= P.MIN_WALL)
    add("carrier: pad window to the centre bore web (mm)", round(w["radial_from_mm"] - CR["centre_bore_d"] / 2, 3),
        w["radial_from_mm"] - CR["centre_bore_d"] / 2 >= P.MIN_WALL)
    s0, s1, sw = DESIGN["carrier_slot"]
    add("carrier: slot end to the plate rim (mm)", round(CR["r"] - (s1 / 2 + sw / 2), 3),
        CR["r"] - (s1 / 2 + sw / 2) >= P.MIN_WALL)
    L = carrier_lock()
    add("carrier tabs per " + L["source"] + ": tip r / ring lip r_in / overlap (mm)",
        [L["tab_r_out"], CR["ring"]["r_in_at_plate"], round(L["tab_r_out"] - CR["ring"]["r_in_at_plate"], 3)],
        L["tab_r_out"] - CR["ring"]["r_in_at_plate"] >= 0.8)
    # gland pose
    c, u = gland_frame()
    inner = c - u * P.GLAND["thread_length"]
    add("gland spot-face centre (x, r)", [round(c.X, 4), round(math.hypot(c.Y, c.Z), 4)],
        abs(c.X - P.GLAND["spot_face_centre_x"]) < 1e-6)
    add("gland thread inner end (x, r) vs GLAND inner_end", [round(inner.X, 4), round(math.hypot(inner.Y, inner.Z), 4)],
        abs(inner.X - P.GLAND["inner_end_x"]) < 0.01 and abs(math.hypot(inner.Y, inner.Z) - P.GLAND["inner_end_r"]) < 0.01)
    gl = by["cable_gland_m12"]["shape"]
    add("gland body to the carrier plate (mm)", round(min_distance(gl, car), 3), min_distance(gl, car) >= 1.0)
    _, oi = oring()
    add("O-ring installed section a x b (mm), centre (x, r)", [round(2 * oi["a"], 3), round(2 * oi["b"], 3),
                                                              round(oi["x_c"], 3), round(oi["r_c"], 3)],
        oi["x_c"] + oi["a"] <= J["groove_x"]["x_aft"] and oi["r_c"] + oi["b"] <= OR["bore_d"] / 2)
    return res


def assembly_sequence() -> list[dict]:
    """Order of assembly with the physical check of each step (numbers from params / this module)."""
    D = DESIGN
    L = carrier_lock()
    entry = L["entry_theta_deg"]
    return [
        dict(step=1, sub="generator (bench)", action="Press 4 M3 inserts into the replica stator face and 2 into "
             "the bell bosses (soldering iron, 4.0 mm holes).", check="display build only; the real GL40 has M3 "
             "threads (stator 3.5 deep, rotor 2.5 deep).", ok=True),
        dict(step=2, sub="generator (bench)", action="Stack carrier + adapter spacer + stator and fit 4 ISO 7380 "
             "M3 x 8 from the carrier counterbores (2 mm key).", check="thread engagement 3.2 mm (insert 5.7, real "
             "motor hole 3.5); heads 0.05 below the carrier face; motor pattern clocked 45 deg so the pad window "
             "and the cable notch face theta 180.", ok=True),
        dict(step=3, sub="generator (bench)", action="Slide the bell onto the stator journal; bolt claw A to the bell "
             "face with 2 x M3 x 4 ISO 4762 at 0 / 180 (2.5 mm key down the dog counterbores).",
             check="journal 12.0 in 12.6 (CLEAR_ROTATE per side), 7.5 mm long; v2: claw A bolts to the GL40 rotor "
             "face (rotor holes 2.5 deep, 2.0 engaged).", ok=True),
        dict(step=4, sub="front shell", action="Screw the M12 gland into the tapped seat, pass the cable, leave the "
             "slack loop; solder the leads (v2) and lay the cable in the carrier notch.", check="cable OD <= 6.5 in "
             "a 7.5 notch; the pad window (r 12-20.5) opens into the notch.", ok=True),
        dict(step=5, sub="front shell", action=f"Push carrier + generator + claw A in through the socket with the tabs "
             f"at theta {entry[0]:.0f} / {entry[1]:.0f} / {entry[2]:.0f}, through the ring notches, then twist "
             f"{L['twist_deg']:.0f} deg toward -theta to lock at {L['lock_theta_deg'][0]:.0f} / "
             f"{L['lock_theta_deg'][1]:.0f} / {L['lock_theta_deg'][2]:.0f}.",
             check=f"largest part to pass the socket lip D {J['lip_d']}: carrier tabs D {2 * L['tab_r_out']:.1f}, "
             f"claw A D {D['claw_a_od']}, bell D {G['D']}. The reaction torque (rotor omega +X) pushes the "
             "stator toward -theta, i.e. into the lock stop.", ok=2 * L["tab_r_out"] < J["lip_d"]),
        dict(step=6, sub="aft shell", action="Press the 688 (through the 22 seat and the 18 bore) onto the seal "
             "shoulder, then the 608 onto its 18 shoulder (draw bolt).", check="688 OD 16 < 18: passes (layout).",
             ok=True),
        dict(step=7, sub="aft shell", action="Slide the shaft in from the front, aft end first, through 608 and 688 "
             f"until it stands {S['x1'] - P.STATIC_TAIL_X:.1f} mm proud of the static tail face (calipers).",
             check="nothing larger than 8 mm passes the bearings: the hub D-flat goes through the 608 / 688 bores "
             "(slip fit, no seal yet).", ok=True),
        dict(step=8, sub="aft shell", action="Shim on the shaft against the 608 inner race; claw B (collar pressed "
             "in) onto the shaft front until the collar sits on the shim; tighten the collar set screw on the "
             "flat through claw B's access hole (theta 0).", check=f"claw B access hole at x "
             f"{(63.0 - CP['shim_t'] - 4.0):.1f}, 1 mm aft of the equator: the open spigot (bore "
             f"D {J['spigot_bore_d']}) leaves {J['spigot_bore_d'] / 2 - CP['od'] / 2:.1f} mm of radial room for a "
             "short hex key; the shaft end stays at the gauge setting.", ok=J["spigot_bore_d"] / 2 - CP["od"] / 2 > 12),
        dict(step=9, sub="aft shell (v2)", action="Lip seal: slide the brass sleeve (ID 8.1, OD 9) over the "
             "protruding shaft end (covers the D-flat and the tapped end), push the seal over it, main lip aft, "
             "into the 16 bore onto the shoulder (a 14 mm socket as a driver), remove the sleeve.",
             check="fixes the layout flag: the seal goes on from the free shaft end, so the flat edges never "
             "touch the lip; the sleeve never has to pass a bearing.", ok=True),
        dict(step=10, sub="rotor", action="Rotor per hydroegg.cad.rotor (hub on the D-flat with a 0.6 feeler, "
             "spinner, M3 screw into the tapped shaft end).", check="thrust path: hub -> M3 screw -> shaft -> "
             "collar -> shim -> 608 inner race.", ok=True),
        dict(step=11, sub="aft shell", action="Grease the -141 O-ring and roll it into the spigot groove.",
             check="installed: stretched 2.5 % onto D 59.88, squeezed into the 2.105 gland.", ok=True),
        dict(step=12, sub="close", action="Align the bayonet entry marks and push the shells together; the claws "
             "find each other on the V lead-ins (turn the rotor a little if two dog tips meet), twist 35 deg.",
             check=f"dogs 2 x 90 deg each side, flank gap {D['claw_flank_gap']} mm, tips chamfered "
             f"{D['claw_tip_chamfer']} x 45; float {CP['float']} mm both ways.", ok=True),
    ]


def integration_check() -> dict:
    """Mesh (manifold3d) interference of the drivetrain and fins against the newest shells BREP cache
    and the rotor parts, the carrier's insertion path / twist and the generator stack going in.
    Exact OCC booleans against the full shells take minutes each; meshes at 0.03-0.05 mm give the
    same answer to the tessellation tolerance (overlaps below ~1 mm3 on curved contacts are facets)."""
    import glob
    import os
    import manifold3d as m3
    import trimesh
    from build123d import import_brep
    from scipy.spatial import cKDTree

    dirs = glob.glob(str(ROOT / "out/work/cad/shells/cache/*"))
    if not dirs:
        return {"skipped": "no shells cache"}
    d = max(dirs, key=os.path.getmtime)

    def mesh(sh, tol=0.04):
        v, f = sh.tessellate(tol, 0.25)
        return trimesh.Trimesh(np.array([(q.X, q.Y, q.Z) for q in v]), np.array(f), process=True)

    def mf(t):
        mm = m3.Mesh(vert_properties=np.asarray(t.vertices, np.float32), tri_verts=np.asarray(t.faces, np.uint32))
        mm.merge()
        return m3.Manifold(mm)

    def inter(a, b):
        return float((mf(a) ^ mf(b)).volume())

    def gap(a, b):
        pa = np.vstack([trimesh.sample.sample_surface(a, 60000, seed=1)[0], a.vertices])
        pb = np.vstack([trimesh.sample.sample_surface(b, 150000, seed=2)[0], b.vertices])
        lo, hi = np.maximum(a.bounds[0], b.bounds[0]) - 3, np.minimum(a.bounds[1], b.bounds[1]) + 3
        sel = np.all((pa >= lo) & (pa <= hi), axis=1)
        return float(cKDTree(pb).query(pa[sel])[0].min()) if sel.any() else float("inf")

    SH = {n: mesh(import_brep(f"{d}/{n}.brep"), 0.03) for n in ("front_shell", "aft_shell")}
    mine = {p["name"]: p["shape"] for p in parts()}
    try:
        from .fins import parts as fin_parts
        mine.update({p["name"]: p["shape"] for p in fin_parts()})
    except Exception:  # noqa: BLE001
        pass
    expect = {"bearing_608": "press (crush ribs)", "bearing_688": "press (crush ribs)",
              "cable_gland_m12": "seated on the spot face", "fin_0": "on the end stop", "fin_1": "on the end stop",
              "fin_cross_pins": "snug in the D2.0 pad hole", "motor_carrier": "tabs on the hard stop"}
    rows = []
    for nm, sh in (("front_shell", ("motor_carrier", "adapter_spacer", "generator_stator", "generator_bell",
                                    "claw_a", "screws_carrier_m3x8", "cable_gland_m12", "fin_0", "fin_1",
                                    "fin_cross_pins", "oring_141")),
                   ("aft_shell", ("claw_b", "shaft_collar", "shim_8x12x0p5", "shaft", "oring_141",
                                  "lip_seal_tc8x16x5", "bearing_608", "bearing_688", "claw_a", "generator_bell"))):
        for n in sh:
            if n not in mine:
                continue
            m = mesh(mine[n])
            v = inter(m, SH[nm])
            rows.append(dict(part=n, shell=nm, overlap_mm3=round(v, 3), min_gap_mm=round(gap(m, SH[nm]), 3),
                             expected=expect.get(n, "clear")))
    car = mine["motor_carrier"]
    L = carrier_lock()
    path = []
    for x in (40.0, 32.0, 26.0, 23.0, 21.0, 19.5, 18.5):
        path.append(dict(x_front=x, twist_deg=L["twist_deg"], overlap_mm3=round(inter(mesh(
            at_azimuth(car, L["twist_deg"]).moved(Location((x - 18.5, 0, 0)))), SH["front_shell"]), 3)))
    for t in (20.0, 10.0, 5.0, 0.0):
        path.append(dict(x_front=18.5, twist_deg=t, overlap_mm3=round(inter(mesh(at_azimuth(car, t)), SH["front_shell"]), 3)))
    stack = []
    for dx in (40.0, 20.0, 10.0, 4.0, 0.0):
        v = sum(inter(mesh(at_azimuth(mine[n], L["twist_deg"]).moved(Location((dx, 0, 0)))), SH["front_shell"])
                for n in ("adapter_spacer", "generator_stator", "generator_bell", "claw_a"))
        stack.append(dict(dx=dx, overlap_mm3=round(v, 3)))
    rot = []
    try:
        from .rotor import parts as rotor_parts
        for p in rotor_parts():
            if p["name"] in ("hub", "spinner", "hub_grub", "spinner_screw"):
                m = mesh(p["shape"])
                for n in ("shaft", "lip_seal_tc8x16x5"):
                    a = mesh(mine[n])
                    rot.append(dict(mine=n, rotor=p["name"], overlap_mm3=round(inter(a, m), 3),
                                    min_gap_mm=round(gap(a, m), 3)))
    except Exception as ex:  # noqa: BLE001
        rot.append(dict(error=str(ex)))
    return dict(shells_cache=str(Path(d).relative_to(ROOT)), method="manifold3d mesh booleans (0.03-0.04 mm) and "
                "surface-sample gaps (+-0.03 mm)", rows=rows, carrier_entry_and_twist=path,
                generator_stack_entry=stack, rotor=rot)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def export_all(ps=None) -> list[str]:
    ps = ps or parts()
    (OUT / "step").mkdir(parents=True, exist_ok=True)
    (OUT / "stl").mkdir(parents=True, exist_ok=True)
    written = []
    for p in ps:
        if p["role"] == "hidden" and p["bought"]:
            continue
        sp = OUT / "step" / f"{p['name']}.step"
        export_step(p["shape"], str(sp))
        written.append(str(sp.relative_to(ROOT)))
        if p["print"] is not None:
            st = OUT / "stl" / f"{p['name']}.stl"
            export_stl(p["print"], str(st), tolerance=0.02, angular_tolerance=0.15)
            written.append(str(st.relative_to(ROOT)))
    return written


def _stl_watertight(path: Path) -> bool:
    import trimesh
    m = trimesh.load(path)
    return bool(m.is_watertight)


def render_previews(module: str = "hydroegg.cad.drivetrain", views=None) -> list[str]:
    """Preview PNGs through scripts/preview.py (it takes the render lock)."""
    import subprocess
    views = views or [("preview_parts", "drivetrain_iso", "view=-35,22"),
                      ("preview_coupling", "drivetrain_coupling", "view=60,25"),
                      ("preview_generator", "drivetrain_generator", "view=-60,25"),
                      ("preview_print", "drivetrain_print", "view=-30,40"),
                      ("preview_in_shells", "drivetrain_in_shells", "view=-35,22")]
    out = []
    for fn, stem, q in views:
        png = WORK / f"{stem}.png"
        r = subprocess.run([sys.executable, str(ROOT / "scripts/preview.py"), f"{module}:{fn}", "--out", str(png),
                            "--q", q], cwd=ROOT, capture_output=True, text=True)
        if r.returncode == 0:
            out.append(str(png.relative_to(ROOT)))
    return out


def main() -> None:
    import time
    t0 = time.time()
    ps = parts()
    res = checks(ps)
    files = export_all(ps)
    for c in res["parts"]:
        stl = OUT / "stl" / f"{c['name']}.stl"
        if stl.exists():
            c["stl_watertight"] = _stl_watertight(stl)
    res["assembly"] = assembly_sequence()
    res["design"] = DESIGN
    if "--integration" in sys.argv:
        res["integration"] = integration_check()
    res["files"] = files
    res["elapsed_s"] = round(time.time() - t0, 1)
    if "--no-render" not in sys.argv:
        res["renders"] = render_previews()
    (OUT / "checks").mkdir(parents=True, exist_ok=True)
    (OUT / "checks" / "drivetrain.json").write_text(json.dumps(res, indent=1, default=str))
    bad = [c["name"] for c in res["parts"] if not c["valid"]]
    print(f"parts {len(res['parts'])}, invalid {bad}, interfaces failing "
          f"{[i['check'] for i in res['interfaces'] if not i['ok']]}, {res['elapsed_s']} s")
    for c in res["parts"]:
        print(f"  {c['name']:24s} valid={c['valid']} solids={c['n_solids']} vol={c['volume_cm3']} "
              f"wt={c.get('watertight')} ovh={c.get('overhang_area_mm2_gt46', '-')} bridges={c.get('of_which_bridges_mm2', '-')}")


if __name__ == "__main__":
    main()

