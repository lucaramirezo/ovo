"""Tolerance coupon: one small print (60 x 60 mm footprint) that calibrates every printed fit of
the Ovo before the real parts are printed.

    scripts/cadrun.sh .venv/bin/python -m ovo.cad.coupon

The coupon prints flat (as modelled, z up) together with two small keys that sit in its window
(separate bodies, 1 mm or more of air around them, they come off the bed loose):

* 3 stepped bores, each a 608 seat over a 688 seat, as in the bearing tube (the 688 passes the
  608 seat and stops on the step): 21.85 / 15.85, 21.95 / 15.95, 22.05 / 16.05 (labels).
* 6 pin holes for the ISO 8734 3 mm dowel: press 2.9 / 3.0 / 3.1, running 3.2 / 3.3 / 3.4.
* M3 heat-set insert hole 4.0 x 6.5 deep (CNC Kitchen M3 x 5.7).
* M12x1.5 tap drill 10.5 through 10 mm (tap it and try the cable gland).
* Fin dovetail: the shell's half-dovetail groove (shells.FIN_RAIL) cut vertically into the -X
  face, as it prints in the front shell (nose up); key 'tongue' = the fin tongue inset by
  CLEAR_SLIDE, printed root down like the fin.
* Bayonet: the socket section of the front shell (lip, 45 deg lip cone, twist channel with its
  45 deg ceiling, entry gap) unrolled into the +X face, printed in the shell's orientation; key
  'lug' = the spigot lug section on a short spigot plate, printed tip up like the aft shell. Push
  the lug up through the gap, slide it along the channel to the lock end.

Every value is engraved on the top face (DejaVu Sans Bold, 0.6 mm deep).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from build123d import (
    Axis,
    FontStyle,
    Location,
    Plane,
    RectangleRounded,
    Solid,
    Text,
    Vector,
    export_step,
    export_stl,
    extrude,
)

from .. import params as P
from .common import check, fits_bed, interference, min_distance

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out" / "cad"
FONT = "DejaVu Sans"

W, H = 60.0, 10.0
J, BY = P.JOINT, P.BAYONET

BORES = [  # (x, y, 608 seat D, 688 seat D)
    (-17.5, 17.5, 21.85, 15.85),
    (0.0, 0.0, 21.95, 15.95),
    (17.5, 17.5, 22.05, 16.05),
]
BORE_SPLIT_Z = 5.0           # 688 seat z 0..5 (its B), 608 seat z 5..10
PINS = [  # (x, y, D)
    (-26.0, -16.5, 2.9), (-20.5, -16.5, 3.0), (-15.0, -16.5, 3.1),
    (-26.0, -25.5, 3.2), (-20.5, -25.5, 3.3), (-15.0, -25.5, 3.4),
]
INSERT = (11.0, -23.5, P.INSERT_M3["hole_d"], 6.5)
M12 = (21.0, -21.0, 10.5)
WINDOW = (-8.0, 7.0, -28.5, -15.0)         # x0, x1, y0, y1 (through)
GROOVE_Y = -3.0                             # fin groove centre (w = 0) on the -X face
BAY_Z0 = 50.6                               # z = BAY_Z0 - x_shell maps the socket section
LABEL_H = 2.6


def _rail() -> dict:
    try:
        from .shells import FIN_RAIL as R
        return dict(R, source="ovo.cad.shells.FIN_RAIL")
    except Exception:  # noqa: BLE001
        from .fins import _RAIL_FALLBACK
        return dict(_RAIL_FALLBACK, source="fins fallback copy")


def _shell_consts() -> dict:
    """Socket / lug section numbers, the same expressions as shells._socket_cutters / _lugs."""
    X_TIP, LUG_T, LUG_H = J["x_spigot_tip"], BY["lug_t"], BY["lug_h"]
    R_SPIG_LUG, R_LUG = J["spigot_lug_zone_od"] / 2, J["lug_od"] / 2
    R_LIP, R_CH = J["lip_d"] / 2, J["channel_d"] / 2
    X_RING0, X_CH0, X_SEAL0 = J["x_ring_ceiling_start"], J["x_channel"][0], J["x_lip"][1]
    ax_clear = 0.2                                   # shells.BAY_AX_CLEAR
    x_cone = lambda r: X_TIP + LUG_T + ax_clear - (r - R_SPIG_LUG)      # noqa: E731
    return dict(X_TIP=X_TIP, LUG_T=LUG_T, LUG_H=LUG_H, R_SPIG_LUG=R_SPIG_LUG, R_LUG=R_LUG, R_LIP=R_LIP,
                R_CH=R_CH, X_RING0=X_RING0, X_CH0=X_CH0, X_SEAL0=X_SEAL0, x_cone=x_cone, lug_ch=0.4)


def _bx(r: float) -> float:
    """Coupon x of a shell radius on the +X face (face x = W/2 is the lip, r = R_LIP)."""
    c = _shell_consts()
    return W / 2 - (r - c["R_LIP"])


def _bz(x_shell: float) -> float:
    return BAY_Z0 - x_shell


def _prism_xz(pts_xz, y0: float, y1: float) -> Solid:
    """Prism of an (x, z) polygon between y0 and y1."""
    from build123d import Face, Wire
    f = Face(Wire.make_polygon([Vector(x, y0, z) for x, z in pts_xz], close=True))
    return Solid.extrude(f, Vector(0, y1 - y0, 0))


def _prism_xy(pts_xy, z0: float, z1: float) -> Solid:
    from build123d import Face, Wire
    f = Face(Wire.make_polygon([Vector(x, y, z0) for x, y in pts_xy], close=True))
    return Solid.extrude(f, Vector(0, 0, z1 - z0))


def _zcyl(x, y, r, z0, z1) -> Solid:
    return Solid.make_cylinder(r, z1 - z0, Plane(origin=(x, y, z0)))


def _cone(x, y, r0, r1, z0, z1) -> Solid:
    return Solid.make_cone(r0, r1, z1 - z0, Plane(origin=(x, y, z0)))


def _label(txt: str, x: float, y: float, h: float = LABEL_H, depth: float = 0.6) -> Solid:
    t = Text(txt, h, font=FONT, font_style=FontStyle.BOLD)
    t = t.moved(Location((x, y, H)))
    return extrude(t, -depth)


def _bay_y(theta: float, r: float) -> float:
    """Unrolled azimuth -> y on the +X face (entry gap centre theta 37 at y 0, channel toward -y)."""
    e = BY["entry_theta_deg"][0]
    return -math.radians(theta - e) * r


def bayonet_cutters() -> list[Solid]:
    c = _shell_consts()
    e, l = BY["entry_theta_deg"][0], BY["lock_theta_deg"][0]
    gw, lw = BY["gap_w_deg"], BY["lug_w_deg"]
    xo = W / 2 + 1.0
    zc0, zc1 = _bz(c["X_RING0"]), _bz(c["X_CH0"])
    zk0, zk1 = _bz(c["x_cone"](c["R_CH"])), _bz(c["x_cone"](c["R_LIP"]))
    xf = _bx(c["R_CH"])
    channel = [(W / 2, zc0), (xf, zc1), (xf, zk0), (W / 2, zk1), (xo, zk1), (xo, zc0)]
    y_ch = sorted([_bay_y(e - gw / 2 - 0.5, c["R_CH"]), _bay_y(l + lw / 2, c["R_CH"])])
    za = _bz(c["X_SEAL0"] - (c["R_CH"] - c["R_LIP"]))           # lip aft chamfer top (x 48.095)
    zb = _bz(c["X_SEAL0"])
    bore = [(W / 2, za), (xf, zb), (xf, -1.0), (xo, -1.0), (xo, za)]
    gap = [(xf, _bz(c["x_cone"](c["R_CH"]) - 0.3)), (xf, -1.0), (xo, -1.0), (xo, _bz(c["x_cone"](c["R_CH"]) - 0.3))]
    y_gap = sorted([_bay_y(e - gw / 2, c["R_LIP"]), _bay_y(e + gw / 2, c["R_LIP"])])
    return [_prism_xz(channel, *y_ch), _prism_xz(bore, y_ch[0] - 2.0, y_ch[1] + 2.0),
            _prism_xz(gap, *y_gap)], dict(y_channel=y_ch, y_gap=y_gap, z_ceiling=(zc0, zc1), z_lip_cone=(zk0, zk1),
                                          z_lip_aft=(za, zb), x_floor=xf)


def groove_section_xy(inset: float = 0.0) -> list[tuple[float, float]]:
    """The shell's fin groove (FIN_RAIL) on the -X face, as (x, y): depth into +x, w along y,
    inset > 0 gives the tongue (straight rail: plain normal offsets)."""
    R = _rail()
    w0, wt, wb, hn, d = R["groove_w_offset"], R["groove_w_top"], R["groove_w_bottom"], R["groove_neck_h"], R[
        "groove_depth"]
    c = inset
    w_lo, w_hi = w0 + c, w0 + wt - c
    d_bot = d - c
    slope = (wb - wt) / (d - hn)                    # dw/d(depth) along the flank (> 0)
    delta = -c * math.sqrt(1 + slope ** 2)          # shift toward -w
    d_n = hn + (w_hi - (w0 + wt) - delta) / slope
    w_b = (w0 + wt) + (d_bot - hn) * slope + delta
    x = lambda dep: -W / 2 + dep                    # noqa: E731
    top = -1.0 if inset == 0 else -c
    pts = [(x(top), w_lo), (x(top), w_hi), (x(d_n), w_hi), (x(d_bot), w_b), (x(d_bot), w_lo)]
    return [(px, GROOVE_Y + py) for px, py in pts]


def coupon_block(labels: bool = True) -> Solid:
    blk = extrude(RectangleRounded(W, W, 3.0), H)
    top = blk.edges().group_by(Axis.Z)[-1]
    blk = blk.chamfer(0.6, None, top)
    bot = blk.edges().group_by(Axis.Z)[0]
    blk = blk.chamfer(0.4, None, bot)
    cuts = []
    for x, y, d608, d688 in BORES:
        cuts += [_zcyl(x, y, d688 / 2, -1.0, BORE_SPLIT_Z + 0.01), _zcyl(x, y, d608 / 2, BORE_SPLIT_Z, H + 1.0),
                 _cone(x, y, d608 / 2, d608 / 2 + 1.3, H - 0.3, H + 1.0), _cone(x, y, d688 / 2 + 1.3, d688 / 2, -1.0, 0.3)]
    for x, y, d in PINS:
        cuts += [_zcyl(x, y, d / 2, -1.0, H + 1.0), _cone(x, y, d / 2 + 1.2, d / 2, -1.0, 0.2),
                 _cone(x, y, d / 2, d / 2 + 1.2, H - 0.2, H + 1.0)]
    x, y, d, dep = INSERT
    cuts += [_zcyl(x, y, d / 2, H - dep, H + 1.0), _cone(x, y, d / 2, d / 2 + 1.3, H - 0.3, H + 1.0)]
    x, y, d = M12
    cuts += [_zcyl(x, y, d / 2, -1.0, H + 1.0), _cone(x, y, d / 2, d / 2 + 1.75, H - 0.75, H + 1.0),
             _cone(x, y, d / 2 + 1.4, d / 2, -1.0, 0.4)]
    x0, x1, y0, y1 = WINDOW
    win = extrude(RectangleRounded(x1 - x0, y1 - y0, 1.0).moved(Location(((x0 + x1) / 2, (y0 + y1) / 2, -1.0))),
                  H + 2.0)
    cuts.append(win)
    cuts.append(_prism_xy(groove_section_xy(0.0), -1.0, H + 1.0))
    bc, _ = bayonet_cutters()
    cuts += bc
    # labels (top face)
    lab = [("21.85", -19.0, 3.0), ("15.85", -19.0, -0.6), ("21.95  15.95", 0.0, -12.9),
           ("22.05", 19.5, 3.0), ("16.05", 19.5, -0.6),
           ("2.9", -26.0, -12.6), ("3.0", -20.5, -12.6), ("3.1", -15.0, -12.6),
           ("3.2", -26.0, -21.6), ("3.3", -20.5, -21.6), ("3.4", -15.0, -21.6),
           ("4.0", 11.0, -18.7), ("M12 10.5", 20.0, -12.6), ("FIN", -24.8, -8.6), ("BAY", 23.2, -7.2)]
    for t, lx, ly in (lab if labels else []):
        cuts.append(_label(t, lx, ly, 2.2 if len(t) <= 3 else LABEL_H))
    out = blk
    for c in cuts:
        out = out - c
    return out


# --------------------------------------------------------------------------
# Coupon B (finding T26): the fits coupon A does not calibrate, on a second 60 x 34 block
# --------------------------------------------------------------------------
CB_W, CB_D, CB_H = 60.0, 34.0, 10.0
DBORES = [(-21.0, 8.0, 8.10), (-9.0, 8.0, 8.15), (3.0, 8.0, 8.25)]       # hub D-bore, 0.5 flat, vertical
DBORE_FLAT = P.ROTOR_HUB["d_flat_depth"]
VSEATS = [(-21.5, -4.0, (6.0, 2.0)), (-13.0, -4.0, (6.1, 2.1)), (-4.5, -4.0, (6.2, 2.2))]   # hub block seat
VDISCS = [(5.0, -4.0, 4.0), (12.0, -4.0, 4.1), (19.0, -4.0, 4.2)]      # disc seats, vertical
SEAT_DEPTH = 2.2
HSEATS = [(-22.0, (6.0, 2.0)), (-15.0, (6.1, 2.1)), (-8.0, (6.2, 2.2))]   # blade heel slot, in the -Y face
HDISCS = [(2.0, 4.0), (9.0, 4.1), (16.0, 4.2)]                             # blade shank seat, -Y face
EAR = dict(x0=12.0, x1=28.0, y=8.0, h=12.0, holes=((15.5, 2.9), (20.0, 3.0), (24.5, 3.1)))  # hub ear (3.7 thick)


def coupon_b(labels: bool = True) -> Solid:
    """Second tolerance block: hub D-bore D8.10 / 8.15 / 8.25 with the 0.5 flat (vertical, as the hub
    prints), the deployed-magnet block seats 6.0 x 2.0 / 6.1 x 2.1 / 6.2 x 2.2 and disc seats D4.0 / 4.1
    / 4.2 both vertical (hub, tail boss) and horizontal in the -Y face (blade heel and shank seats, as
    the blade prints), and a 3.7 mm upright ear with horizontal press holes D2.9 / 3.0 / 3.1 (the hub
    ears print with horizontal holes)."""
    from build123d import Box
    blk = extrude(RectangleRounded(CB_W, CB_D, 2.0), CB_H)
    blk = blk.chamfer(0.4, None, blk.edges().group_by(Axis.Z)[0])
    t = P.HUB_EAR_T
    ear = Box(EAR["x1"] - EAR["x0"], t, EAR["h"]).moved(
        Location(((EAR["x0"] + EAR["x1"]) / 2, EAR["y"], CB_H + EAR["h"] / 2 - 0.01)))
    blk = blk.fuse(ear).clean()
    cuts = []
    for x, y, d in DBORES:
        hole = _zcyl(x, y, d / 2, -1.0, CB_H + 1.0)
        flat_r = d / 2 - DBORE_FLAT
        seg = Box(d, d, CB_H + 4.0).moved(Location((x, y + flat_r + d / 2, CB_H / 2)))
        cuts += [hole - seg, _cone(x, y, d / 2 + 1.0, d / 2, -1.0, 0.0)]
    for x, y, (a, b) in VSEATS:
        cuts.append(Box(a, b, SEAT_DEPTH + 1.0).moved(Location((x, y, CB_H - SEAT_DEPTH + (SEAT_DEPTH + 1.0) / 2))))
    for x, y, d in VDISCS:
        cuts.append(_zcyl(x, y, d / 2, CB_H - SEAT_DEPTH, CB_H + 1.0))
    y_face = -CB_D / 2
    for x, (a, b) in HSEATS:                       # 6.x along z, 2.x along x, depth along +y
        cuts.append(Box(b, SEAT_DEPTH + 1.0, a).moved(Location((x, y_face + (SEAT_DEPTH - 1.0) / 2, CB_H / 2))))
    for x, d in HDISCS:
        cuts.append(Solid.make_cylinder(d / 2, SEAT_DEPTH + 1.0, Plane(origin=(x, y_face - 1.0, CB_H / 2),
                                                                     z_dir=(0, 1, 0))))
    for x, d in EAR["holes"]:
        cuts.append(Solid.make_cylinder(d / 2, t + 2.0, Plane(origin=(x, EAR["y"] - t / 2 - 1.0, CB_H + 6.0),
                                                            z_dir=(0, 1, 0))))
    lab = [("8.10", -21.0, 14.0), ("8.15", -9.0, 14.0), ("8.25", 3.0, 14.0),
           ("6.0", -21.5, 0.6), ("6.1", -13.0, 0.6), ("6.2", -4.5, 0.6),
           ("4.0", 5.0, 0.6), ("4.1", 12.0, 0.6), ("4.2", 19.0, 0.6),
           ("H6.0", -22.0, -12.5), ("H6.1", -15.0, -12.5), ("H6.2", -8.0, -12.5),
           ("H4.0", 2.0, -12.5), ("H4.1", 9.0, -12.5), ("H4.2", 16.0, -12.5), ("EAR", 20.0, 13.5)]
    for txt, lx, ly in (lab if labels else []):
        tt = Text(txt, 2.0, font=FONT, font_style=FontStyle.BOLD).moved(Location((lx, ly, CB_H)))
        cuts.append(extrude(tt, -0.5))
    out = blk
    for c in cuts:
        out = out - c
    return out


def key_lug() -> Solid:
    """Spigot plate with one lug, in the coupon frame at its locked position (lug in the channel at
    the lock azimuth, aft face BAY_AX_CLEAR above the lip cone)."""
    c = _shell_consts()
    l, lw = BY["lock_theta_deg"][0], BY["lug_w_deg"]
    x_sp = _bx(c["R_SPIG_LUG"])                     # spigot surface
    y_l = sorted([_bay_y(l - lw / 2, c["R_LUG"]), _bay_y(l + lw / 2, c["R_LUG"])])
    yc = 0.5 * (y_l[0] + y_l[1])
    plate = _prism_xz([(x_sp, 0.0), (x_sp + 3.0, 0.0), (x_sp + 3.0, _bz(c["X_TIP"])), (x_sp, _bz(c["X_TIP"]))],
                      yc - 5.5, yc + 5.5)
    xo = _bx(c["R_LUG"])
    lug = [(x_sp + 0.3, _bz(c["X_TIP"] + c["LUG_T"]) - 0.3), (x_sp, _bz(c["X_TIP"] + c["LUG_T"])),
           (xo, _bz(c["X_TIP"] + c["LUG_T"] - c["LUG_H"])), (xo, _bz(c["X_TIP"] + c["lug_ch"])),
           (xo + c["lug_ch"], _bz(c["X_TIP"])), (x_sp + 0.3, _bz(c["X_TIP"]))]
    k = plate.fuse(_prism_xz(lug, *y_l)).clean()
    edges = [e for e in k.edges() if abs(e.center().X - (x_sp + 3.0)) < 1e-6]
    try:
        k = k.chamfer(0.4, None, edges)
    except Exception:  # noqa: BLE001
        pass
    return k


def key_tongue() -> Solid:
    """The fin tongue inset by CLEAR_SLIDE in the coupon groove, with a 6 mm stub (the fin root)."""
    tg = _prism_xy(groove_section_xy(P.CLEAR_SLIDE), 0.0, H)
    sec = groove_section_xy(P.CLEAR_SLIDE)
    w_lo, w_hi = sec[0][1], sec[1][1]
    x_root = -W / 2 - P.CLEAR_SLIDE
    stub = _prism_xy([(x_root - 6.0, w_lo), (x_root + 0.01, w_lo), (x_root + 0.01, w_hi), (x_root - 6.0, w_hi)], 0.0, H)
    k = tg.fuse(stub).clean()
    return k


def key_print_poses() -> tuple[Solid, Solid]:
    """Both keys moved into the window, in their print orientation, same frame as the block."""
    x0, x1, y0, y1 = WINDOW
    kl = key_lug()                                  # stands on z = 0 (spigot root down, tip up)
    bb = kl.bounding_box()
    kl = kl.moved(Location((x0 + 1.5 - bb.min.X, (y0 + y1) / 2 - bb.center().Y, -bb.min.Z)))
    kt = key_tongue().rotate(Axis.Y, 90.0).rotate(Axis.Z, 90.0)   # tongue bottom (+x face) down, length along y
    bb = kt.bounding_box()
    kt = kt.moved(Location((x1 - 1.5 - bb.max.X, (y0 + y1) / 2 - bb.center().Y, -bb.min.Z)))
    return kl, kt


_CACHE: dict = {}

OFFSET = Location((0, 160.0, 0))                    # world pose: parked beside the model (hidden)


def parts() -> list[dict]:
    if "parts" in _CACHE:
        return _CACHE["parts"]
    blk = coupon_block()
    kl, kt = key_print_poses()
    out = [
        dict(name="tolerance_coupon", shape=blk.moved(OFFSET), print=blk, material="stand", role="hidden",
             rotor=False, qty=1, bought=False, explode=[0, 0, 0],
             notes="Print FIRST, flat as modelled, with the two keys in its window (one 60 x 60 job), same "
                   "material and slicer profile as the real parts. Fit table in coupon.json: pick the 608 / 688 "
                   "seat, pin and insert diameters that give the right feel and adjust PRESS_FIT / hole "
                   "compensation; tap the M12; slide the keys in the dovetail and the bayonet channel."),
        dict(name="tolerance_coupon_key_lug", shape=kl.moved(OFFSET), print=kl, material="stand", role="hidden",
             rotor=False, qty=1, bought=False, explode=[0, 0, 0],
             notes="prints in the coupon window, tip up (as the aft shell's spigot)."),
        dict(name="tolerance_coupon_key_tongue", shape=kt.moved(OFFSET), print=kt, material="stand", role="hidden",
             rotor=False, qty=1, bought=False, explode=[0, 0, 0],
             notes="prints in the coupon window, tongue bottom down (as the fin)."),
    ]
    cb = coupon_b()
    out.append(dict(name="tolerance_coupon_b", shape=cb.moved(Location((0, 230.0, 0))), print=cb, material="stand",
                    role="hidden", rotor=False, qty=1, bought=False, explode=[0, 0, 0],
                    notes="Print with coupon A (flat as modelled, no supports): hub D-bore on a shaft offcut, the "
                          "6 x 2 x 2 block and 4 x 2 disc magnet seats vertical (hub, tail boss) and horizontal in the "
                          "-Y face (blade heel and shank), and a 3.7 mm upright ear with horizontal D2.9 / 3.0 / 3.1 "
                          "press holes like the hub ears."))
    _CACHE["parts"] = out
    return out


def preview_parts() -> list[dict]:
    return [dict(name=p["name"], shape=p["print"], material=m, role="body")
            for p, m in zip(parts(), ("shell", "fin", "blade"))]


def preview_test_fit() -> list[dict]:
    """Block with both keys at their test positions (lug locked in the channel, tongue in the groove)."""
    blk = coupon_block()
    return [dict(name="coupon", shape=blk, material="shell", role="body"),
            dict(name="key_lug", shape=key_lug(), material="blade", role="body"),
            dict(name="key_tongue", shape=key_tongue().moved(Location((0, 0, 0))), material="fin", role="body")]


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def _hole_d(solid: Solid, x: float, y: float, z: float, d: float) -> dict:
    """Probe a vertical hole: points at r = d/2 -/+ 0.01 must be outside / inside the solid."""
    ok = True
    for a in np.linspace(0, 2 * math.pi, 12, endpoint=False):
        pin = Vector(x + (d / 2 - 0.01) * math.cos(a), y + (d / 2 - 0.01) * math.sin(a), z)
        pout = Vector(x + (d / 2 + 0.01) * math.cos(a), y + (d / 2 + 0.01) * math.sin(a), z)
        ok &= (not solid.is_inside(pin)) and solid.is_inside(pout)
    return dict(xy=[x, y], z=z, d=d, ok=bool(ok))


def checks() -> dict:
    import trimesh
    ps = parts()
    blk, kl, kt = (p["print"] for p in ps[:3])
    cb = ps[3]["print"]
    res = dict(parts=[], features=[], interfaces=[])
    for p in ps:
        c = check(p["print"], p["name"])
        v, f = p["print"].tessellate(0.03, 0.2)
        m = trimesh.Trimesh(np.array([(q.X, q.Y, q.Z) for q in v]), np.array(f), process=True)
        n, zc = m.face_normals, m.triangles_center[:, 2]
        ovh = float(m.area_faces[(n[:, 2] < -math.cos(math.radians(45.0)) - 1e-3) & (zc > 0.05)].sum())
        c.update(watertight=bool(m.is_watertight), fits_bed=fits_bed(p["print"]), overhang_area_mm2_gt45=round(ovh, 3))
        res["parts"].append(c)
    F = res["features"]
    for x, y, d608, d688 in BORES:
        F.append(dict(feature="608 seat", **_hole_d(blk, x, y, 7.5, d608)))
        F.append(dict(feature="688 seat", **_hole_d(blk, x, y, 2.5, d688)))
    for x, y, d in PINS:
        F.append(dict(feature="pin hole (press)" if d < 3.15 else "pin hole (run)", **_hole_d(blk, x, y, 5.0, d)))
    F.append(dict(feature="M3 insert hole", **_hole_d(blk, INSERT[0], INSERT[1], H - 3.0, INSERT[2])))
    F.append(dict(feature="M12 tap drill", **_hole_d(blk, M12[0], M12[1], 5.0, M12[2])))
    for x, y, d in DBORES:
        okr = True
        for a in np.radians(np.arange(-170.0, -9.0, 20.0)):             # round half (the flat faces +y)
            okr &= (not cb.is_inside(Vector(x + (d / 2 - 0.01) * math.cos(a), y + (d / 2 - 0.01) * math.sin(a), 5.0))) \
                and cb.is_inside(Vector(x + (d / 2 + 0.01) * math.cos(a), y + (d / 2 + 0.01) * math.sin(a), 5.0))
        fr = d / 2 - DBORE_FLAT
        okf = (not cb.is_inside(Vector(x, y + fr - 0.01, 5.0))) and cb.is_inside(Vector(x, y + fr + 0.01, 5.0))
        F.append(dict(feature="B: hub D-bore D%.2f with a %.1f flat (round half and flat probed)" % (d, DBORE_FLAT),
                      xy=[x, y], ok=bool(okr and okf)))
    for x, y, d in VDISCS:
        F.append(dict(feature="B: disc seat vertical", **_hole_d(cb, x, y, CB_H - 1.0, d)))
    for x, y, (a, b) in VSEATS:
        ok = (not cb.is_inside(Vector(x + b / 2 - 0.02, y + a / 2 - 0.02, CB_H - 1.0))) and \
            cb.is_inside(Vector(x + b / 2 + 0.02, y, CB_H - 1.0)) if False else \
            (not cb.is_inside(Vector(x + a / 2 - 0.02, y + b / 2 - 0.02, CB_H - 1.0))) and \
            cb.is_inside(Vector(x + a / 2 + 0.02, y, CB_H - 1.0)) and cb.is_inside(Vector(x, y + b / 2 + 0.02, CB_H - 1.0))
        F.append(dict(feature="B: block seat vertical %.1f x %.1f" % (a, b), xy=[x, y], ok=bool(ok)))
    for x, d in EAR["holes"]:
        pin_ = Vector(x, EAR["y"], CB_H + 6.0 + d / 2 - 0.02)
        pout = Vector(x, EAR["y"], CB_H + 6.0 + d / 2 + 0.02)
        F.append(dict(feature="B: ear press hole horizontal D%.1f" % d, ok=bool((not cb.is_inside(pin_)) and cb.is_inside(pout))))
    I = res["interfaces"]

    def add(name, value, ok, **kw):
        I.append(dict(check=name, value=value, ok=bool(ok), **kw))

    bb = blk.bounding_box()
    add("coupon footprint incl. keys (mm)", [round(bb.size.X, 2), round(bb.size.Y, 2)],
        bb.size.X <= 60.0 + 1e-6 and bb.size.Y <= 60.0 + 1e-6)
    d1, d2 = min_distance(blk, kl), min_distance(blk, kt)
    add("keys clear of the coupon in the print layout (mm)", [round(d1, 3), round(d2, 3)], min(d1, d2) >= 0.8)
    add("keys clear of each other (mm)", round(min_distance(kl, kt), 3), min_distance(kl, kt) >= 0.8)
    blk0 = coupon_block(labels=False)
    lug = key_lug()
    ax_n = 0.2 / math.sqrt(2.0)                           # 0.2 axial play across the 45 deg lug / lip faces
    v = interference(blk0, lug)
    add("lug key locked in the channel: overlap (mm3)", round(v, 4), v < 1e-3)
    add("lug key locked: min gap to the socket (mm; 0.2 axial on 45 deg faces = 0.141 normal)",
        round(min_distance(blk0, lug), 3), min_distance(blk0, lug) >= ax_n - 1e-3)
    _, info = bayonet_cutters()
    c = _shell_consts()
    yc_lock = lug.bounding_box().center().Y
    gap_c = 0.5 * (info["y_gap"][0] + info["y_gap"][1])
    worst = 1e9
    for dz in np.linspace(-6.0, 0.0, 13):                 # lug coming up through the entry gap
        s = lug.moved(Location((0, gap_c - yc_lock, dz)))
        worst = min(worst, min_distance(blk0, s))
    add("lug key passes the entry gap (min gap over the stroke, mm)", round(worst, 3), worst >= 0.2)
    worst = 1e9
    for t in np.linspace(0.0, 1.0, 11):                   # twist (slide) from the gap to the lock end
        s = lug.moved(Location((0, (gap_c - yc_lock) * (1 - t), 0)))
        worst = min(worst, min_distance(blk0, s))
    add("lug key slides along the channel to the lock end (min gap, mm)", round(worst, 3), worst >= ax_n - 1e-3)
    tg = key_tongue()
    v = interference(blk0, tg)
    add("tongue key in the fin groove: overlap (mm3)", round(v, 4), v < 1e-3)
    add("tongue key to groove walls (mm, CLEAR_SLIDE)", round(min_distance(blk0, tg), 3),
        abs(min_distance(blk0, tg) - P.CLEAR_SLIDE) < 0.02)
    res["bayonet_map"] = dict(info, z_map="z = %.1f - x_shell, x = %.0f - (r - R_LIP)" % (BAY_Z0, W / 2),
                              note="unrolled at the feature radii (gap at R_LIP, channel at R_CH, lug at R_LUG)")
    res["rail_source"] = _rail()["source"]
    res["how_to_use"] = [
        "608 / 688: press each bearing with a vise; the right seat holds the bearing by hand but presses in "
        "without cracking. The tube's crush ribs (shells) target the same diameter (D 21.9 / 15.9 = nominal + "
        "PRESS_FIT per side); if the good seat is not 21.95 / 15.95, shift the rib tips by the difference.",
        "pins: press holes must hold a 3 mm dowel tight (hub ears, PRESS_FIT); running holes must let the "
        "blade knuckle turn freely (CLEAR_ROTATE); pick the pair and adjust hole compensation.",
        "insert: the CNC Kitchen insert must sink flush and hold an M3 screw.",
        "M12: tap M12x1.5 through and screw the gland in; 10 mm of thread.",
        "keys: the tongue slides in the groove without play you can feel; the lug passes the gap and slides "
        "to the lock end (0.25 radial, 0.2 axial play).",
        "coupon B, D-bore: an 8 mm g6 shaft offcut with a 0.5 flat must slide in by hand with no rock: pick the "
        "bore (8.10 / 8.15 / 8.25) and set ROTOR_HUB.bore_d.",
        "coupon B, magnet seats: a 6 x 2 x 2 block and a 4 x 2 disc must drop into the vertical and the "
        "horizontal (H) seats with a little play for the CA glue: pick the sizes for the hub (vertical), the "
        "blade (horizontal) and the tail boss (vertical, prints at 45 deg).",
        "coupon B, ear: the 3 mm dowel must press into the horizontal ear holes without splitting the 3.7 mm "
        "ear: pick 2.9 / 3.0 / 3.1 for the hub ears (they print horizontal, unlike the pin row of coupon A).",
    ]
    return res


def main() -> None:
    import time
    import trimesh
    t0 = time.time()
    ps = parts()
    res = checks()
    (OUT / "step").mkdir(parents=True, exist_ok=True)
    (OUT / "stl").mkdir(parents=True, exist_ok=True)
    files = []
    for p in ps:
        sp = OUT / "step" / f"{p['name']}.step"
        export_step(p["print"], str(sp))
        files.append(str(sp.relative_to(ROOT)))
        st = OUT / "stl" / f"{p['name']}.stl"
        export_stl(p["print"], str(st), tolerance=0.02, angular_tolerance=0.15)
        files.append(str(st.relative_to(ROOT)))
        for c in res["parts"]:
            if c["name"] == p["name"]:
                c["stl_watertight"] = bool(trimesh.load(st).is_watertight)
                c["tessellation_watertight"] = c.pop("watertight")
                c["watertight"] = c["stl_watertight"]          # the printed file is what counts
    res["files"] = files
    import sys
    if "--no-render" not in sys.argv:
        from .drivetrain import render_previews
        res["renders"] = render_previews("ovo.cad.coupon", [
            ("preview_parts", "coupon_print", "view=-30,55"),
            ("preview_test_fit", "coupon_bayonet_face", "view=90,12&zoom=1.6")])
    res["elapsed_s"] = round(time.time() - t0, 1)
    (OUT / "checks").mkdir(parents=True, exist_ok=True)
    (OUT / "checks" / "coupon.json").write_text(json.dumps(res, indent=1, default=str))
    for c in res["parts"]:
        print(c)
    for f in res["features"]:
        print(f)
    for i in res["interfaces"]:
        print(i)
    print("elapsed", res["elapsed_s"])


if __name__ == "__main__":
    main()
