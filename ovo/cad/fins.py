"""Roll fins of the Ovo: two identical Clark Y fins on recessed dovetail rails of the front shell,
and the two flush blanking strips that close the rails when the fins are off (FIN block).

    scripts/cadrun.sh .venv/bin/python -m ovo.cad.fins

Owner decision 2026-10-05: no proud pads. The fin root sits on the hull skin; the rail is a straight
half-dovetail groove sunk into a boss on the inside of the front shell wall (ovo.cad.shells:
FIN_RAIL, rail_frame); stowed, the fins come off and a blanking strip of the groove's section, its
outer face on the hull profile, keeps the egg at 75 mm.

Geometry (FIN = shells.FIN = params.FIN, written back from the fin study's FIN proposal in
data/layout_internal.json; FIN_MAGNET for the retention pair):

* Fin frame (shells.fin_frame / rail_frame, layout_internal.fin_frame): at the fin angle theta, u runs
  along the yawed root chord from the LE, w across it (+ toward +theta), a radial at the fin angle.
  The chord is yawed by the incidence (TE toward -theta) about the radial line through the yaw pivot.
  Every solid is built at theta 0 in WORLD coordinates and carried to theta by common.at_azimuth.
* Section: Clark Y from the classic ordinates, re-referenced to its chord (LE to mid-TE, the frame
  NeuralFoil and the layout use for the incidence), TE thickened symmetrically to MIN_TE; envelope
  centred on the rail (w = 0) as the layout's fold-clearance model assumes.
* Planform: the datum root line runs through the skin at the chord ends (FIN root_r_le / root_r_te,
  + root_gap); the tip line is the datum + span (straight), its corners rounded. The root FOLLOWS the
  skin: the airfoil is cut by the hull offset by root_gap (0.1 with the tongue seated on the floor).
* Root tongue (the keel): the shell groove inset by CLEAR_SLIDE on its walls and by 0.05 over the
  floor (the magnet pair pulls the tongue onto the floor), from ahead of the slot entry to the end stop
  (contact), cut by the same hull offset so it fills the slot up to the fin root; in the entry relief
  it fills the full groove width. The rail climbs inward going aft (slope -0.25): the fin slides aft
  and in onto the stop, drag pushes it home, and pulled forward every root point lifts off the skin.
* Retention: an N52 4 x 2 x 2 block (magnetised through 2) in a seat in the tongue bottom over a second
  one in the groove floor (shells): the cross pin of the proud pad cannot cross a flush groove (a
  tangential hole at the pin depth leaves the hull at a grazing angle about 12 mm each side of the
  slot). The pair holds normal to the floor; along the rail it holds by the shear of the pair
  (fins.json: both computed with the layout's calibrated charge model).
* Blanking strip: the same tongue section from the entry to the stop with its outer face ON the hull
  profile (flush when seated), the same magnet seat.
* Print: fin TIP DOWN (span vertical, the tongue on top: every face within 45 deg, the undercut flank
  of the tongue at 45 deg, the magnet seat opens upward); strip on its floor face (the seat opens on the
  bed). Both fins are the same part (theta FIN[0] and FIN[0] + 180 = the first turned 180 deg about X).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from build123d import (
    Axis,
    Compound,
    Edge,
    Face,
    Location,
    Plane,
    Solid,
    Vector,
    Wire,
    export_step,
    export_stl,
)

from .. import params as P
from .. import profile
from . import shells as SH
from .common import at_azimuth, check, fits_bed, interference, min_distance

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out" / "cad"
WORK = ROOT / "out" / "work" / "cad" / "fins"
F = SH.FIN
R = SH.FIN_RAIL

DESIGN = dict(
    root_gap=F.get("root_gap", 0.1),   # fin root face above the hull skin (tongue on the floor)
    clear_side=R["clear_side"],        # tongue / strip to the groove walls and flank (normal)
    clear_floor=R["clear_floor"],      # tongue / strip bottom over the floor (normal)
    tip_corner_r=(6.0, 2.5),           # planform rounding at the LE-tip / TE-tip corners
    n_section=90,                      # points per airfoil surface
    strip_tip_min=0.6,                 # the strip's entry wedge is trimmed where thinner than this (review P5)
    strip_tip_break=0.3,               # 45 deg break on the top edge of that trimmed end
    u_ahead=8.0,                       # tongue prism starts this far ahead of the slot entry (outside)
    root_flare_deg=43.0,               # review C4 / P3: root flare from the neck edges to the airfoil skin
                                       # (43 in its section: the root face's slope along the chord adds ~1 deg)
    tip_land_deg=42.0,                 # tip corner rounds stop where their face leans this from vertical
    nail_groove=dict(from_stop=2.6, w=1.2, depth=0.6),   # review C6: across the strip's outer face, steep
                                       # forward wall (the nail pulls on it), 45 deg aft flank
)

# --------------------------------------------------------------------------
# Clark Y coordinates generator
# --------------------------------------------------------------------------

# Classic Clark Y ordinates (percent chord, measured from the flat base line): station, upper, lower.
CLARKY_TABLE = [
    (0.0, 3.50, 3.50), (1.25, 5.45, 1.93), (2.5, 6.50, 1.47), (5.0, 7.90, 0.93), (7.5, 8.85, 0.63),
    (10.0, 9.60, 0.42), (15.0, 10.69, 0.15), (20.0, 11.36, 0.03), (30.0, 11.70, 0.0), (40.0, 11.40, 0.0),
    (50.0, 10.52, 0.0), (60.0, 9.15, 0.0), (70.0, 7.35, 0.0), (80.0, 5.22, 0.0), (90.0, 2.80, 0.0),
    (95.0, 1.49, 0.0), (100.0, 0.12, 0.0),
]


def clarky(n: int = 90, te: float = 0.0) -> tuple[np.ndarray, np.ndarray, dict]:
    """Clark Y for a unit chord in its chord frame: LE at (0, 0), TE mid-point at (1, 0).

    The classic ordinates are interpolated in s = sqrt(x) (a monotone cubic per surface keeps the
    round nose: dx/ds = 0 at the LE with a finite dy/ds), cosine spaced, then rotated about the
    LE so the chord (LE to the mid-point of the TE) is the x axis and scaled to unit length.
    te > 0 thickens the section symmetrically with an x^2 ramp so the TE is te thick (the mean
    line is unchanged). Returns (upper, lower) as (n, 2) arrays from the LE to the TE, and info."""
    from scipy.interpolate import PchipInterpolator

    t = np.array(CLARKY_TABLE, float) / 100.0
    s = np.sqrt(t[:, 0])
    fu, fl = PchipInterpolator(s, t[:, 1]), PchipInterpolator(s, t[:, 2])
    beta = np.linspace(0.0, math.pi, n)
    xs = (1 - np.cos(beta)) / 2
    up = np.column_stack([xs, fu(np.sqrt(xs))])
    lo = np.column_stack([xs, fl(np.sqrt(xs))])
    le = up[0].copy()
    te_mid = 0.5 * (up[-1] + lo[-1])
    ang = math.atan2(te_mid[1] - le[1], te_mid[0] - le[0])
    L = math.hypot(*(te_mid - le))
    rot = np.array([[math.cos(-ang), -math.sin(-ang)], [math.sin(-ang), math.cos(-ang)]])
    up = ((up - le) @ rot.T) / L
    lo = ((lo - le) @ rot.T) / L
    t_te0 = float(up[-1, 1] - lo[-1, 1])
    if te > t_te0:
        d = 0.5 * (te - t_te0) * np.clip(up[:, 0], 0, 1) ** 2
        up[:, 1] += d
        lo[:, 1] -= 0.5 * (te - t_te0) * np.clip(lo[:, 0], 0, 1) ** 2
    th = np.interp(np.linspace(0, 1, 401), up[:, 0], up[:, 1]) - np.interp(np.linspace(0, 1, 401), lo[:, 0], lo[:, 1])
    flat = lo[(lo[:, 0] > 0.3) & (lo[:, 0] < 0.98)]
    k = np.polyfit(flat[:, 0], flat[:, 1], 1)[0] if len(flat) > 3 else float("nan")
    info = dict(base_line_to_chord_deg=round(math.degrees(ang), 4), max_t_c=round(float(th.max()), 5),
                x_max_t=round(float(np.linspace(0, 1, 401)[int(np.argmax(th))]), 3),
                te_c=round(float(up[-1, 1] - lo[-1, 1]), 5), flat_bottom_to_chord_deg=round(math.degrees(math.atan(k)), 3),
                y_min=round(float(lo[:, 1].min()), 5), y_max=round(float(up[:, 1].max()), 5))
    return up, lo, info


def compare_uiuc(up: np.ndarray, lo: np.ndarray) -> dict:
    """Max deviation (fraction of chord) from the UIUC / aerosandbox 'clarky' coordinates the
    layout study fed to NeuralFoil (same chord frame), and the NeuralFoil alpha_L0 of both."""
    try:
        import aerosandbox as asb
    except Exception:  # noqa: BLE001
        return {}
    c = asb.Airfoil("clarky").coordinates
    i = int(np.argmin(c[:, 0]))
    u_ref, l_ref = c[: i + 1][::-1], c[i:]
    xq = np.linspace(0.002, 0.98, 300)
    du = np.interp(xq, up[:, 0], up[:, 1]) - np.interp(xq, u_ref[:, 0], u_ref[:, 1])
    dl = np.interp(xq, lo[:, 0], lo[:, 1]) - np.interp(xq, l_ref[:, 0], l_ref[:, 1])
    out = dict(max_dev_upper_c=round(float(np.abs(du).max()), 5), max_dev_lower_c=round(float(np.abs(dl).max()), 5))
    try:
        import neuralfoil as nf
        al = np.arange(-2.0, 5.01, 0.5)
        res = {}
        for nm, coords in (("uiuc", c), ("cad", np.vstack([up[::-1], lo[1:]]))):
            r = nf.get_aero_from_airfoil(asb.Airfoil(name=nm, coordinates=coords), alpha=al, Re=55000.0, n_crit=5.0,
                                         model_size="large")
            sl, ic = np.polyfit(np.radians(al), np.asarray(r["CL"]), 1)
            res[nm] = dict(alpha_L0_deg=round(-math.degrees(ic / sl), 3), a0_per_rad=round(float(sl), 3))
        out["neuralfoil_Re55k"] = res
    except Exception as ex:  # noqa: BLE001
        out["neuralfoil_error"] = str(ex)
    return out



# --------------------------------------------------------------------------
# Fin frame (theta 0, world coordinates)
# --------------------------------------------------------------------------


def frame() -> dict:
    """shells.rail_frame plus the fin's datum root line a_root(u) = root_r_le + k ci u. At theta 0 the
    fin-frame point (u, w, a) is o + u u_dir + w w_dir + a z (shells.fin_frame / rail_frame)."""
    RF = SH.rail_frame()
    Fr = RF["F"]
    k = (F["root_r_te"] - F["root_r_le"]) / (Fr["x_te"] - Fr["x_le"])

    def a_root(u):
        return F["root_r_le"] + k * Fr["ci"] * u

    def W(u, w, a):
        q = Fr["o"] + u * Fr["u"] + w * Fr["w"]
        return Vector(float(q[0]), float(q[1]), float(a))

    return dict(RF=RF, Fr=Fr, c=F["chord"], k=k, a_root=a_root, W=W, ci=Fr["ci"], si=Fr["si"])


def section_face(fr: dict, a0: float) -> tuple[Face, dict]:
    """Airfoil section in the horizontal plane a = a0 (theta 0), chord along u from the LE, envelope
    centred on w = 0."""
    c = fr["c"]
    up, lo, info = clarky(DESIGN["n_section"], te=P.MIN_TE / c)
    w_c = -0.5 * (info["y_min"] + info["y_max"]) * c        # chord line offset: envelope centred on the rail
    loop = np.vstack([up[::-1], lo[1:]])                     # TE upper -> LE -> TE lower
    pts = [fr["W"](c * x, w_c + c * y, a0) for x, y in loop]
    # Findings T19 / R1: one spline through the whole loop made a skin face whose GProp volume read 14 %
    # low (default) or 2.3 x (adaptive); the T19 face split fixed GProp but its seams tessellated into
    # folded slivers. Seven splines (three per surface, the nose piece round the LE from x/c ~0.07 upper to
    # ~0.07 lower), G1 at the joints (shared tangents from the loop's central differences), give GProp
    # within 0.3 % (default and adaptive) with no extra seams; probed variants (out/work/cad/polish_fins/
    # gprop_probe2.py): 3 or 5 pieces read 1.5 % low / 19 % high, and a joint AT the LE made the planform
    # boolean fail.
    i_le = len(up) - 1
    i_u = int(np.argmin(np.abs(loop[:i_le, 0] - 0.3)))
    i_l = i_le + int(np.argmin(np.abs(loop[i_le:, 0] - 0.3)))
    P_ = np.array([(p.X, p.Y, p.Z) for p in pts])

    def tan(i):
        d = P_[i + 1] - P_[i - 1]
        return Vector(*(d / np.linalg.norm(d)))
    n_ = len(pts) - 1
    cuts = [0, i_u // 2, i_u, (i_u + i_le) // 2, (i_le + i_l) // 2, i_l, (i_l + n_) // 2, n_]   # LE inside a piece
    edges = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        ta = tan(a) if a not in (0, len(pts) - 1) else None
        tb = tan(b) if b not in (0, len(pts) - 1) else None
        if ta is None and tb is None:
            edges.append(Edge.make_spline(pts[a:b + 1]))
        elif ta is None:
            edges.append(Edge.make_spline(pts[a:b + 1], tangents=[pts[a + 1] - pts[a], tb]))
        elif tb is None:
            edges.append(Edge.make_spline(pts[a:b + 1], tangents=[ta, pts[b] - pts[b - 1]]))
        else:
            edges.append(Edge.make_spline(pts[a:b + 1], tangents=[ta, tb]))
    edges.append(Edge.make_line(pts[-1], pts[0]))
    return Face(Wire(edges)), dict(info, chord_line_w=round(w_c, 4), section_splines=7)


def _hull_offset(gap: float) -> Solid:
    """Solid of revolution inside the hull skin offset by gap (radially), over the fin's x range."""
    x0, x1 = 12.0, 62.0
    return SH.revolve_nodes([("L", (x0, 0.0)), ("S", SH.curve(lambda x: SH.r_out(x) + gap, x0, x1, 140, cos=True)),
                             ("L", (x1, 0.0))])


def tongue_section(RF: dict, H: float = 20.0) -> tuple[list, list]:
    """(dovetail, relief) sections (w, h above the floor, given in the plane u = const) of the tongue / strip:
    the groove inset by clear_side normal to its walls and flank and by clear_floor normal to the floor.
    The rail climbs with slope s, so a wall shifted inside the section plane moves less than that in 3D
    unless it contains the w direction: the shifts below give the TRUE normal distances."""
    c, cf, s, k = DESIGN["clear_side"], DESIGN["clear_floor"], RF["s"], RF["k"]
    w0, w_n, w_b, hn = RF["w0"], RF["w_n"], RF["w_b"], R["uc_h"]
    w_lo, w_hi = w0 + c, w_n - c
    z_bot = cf * k                                          # floor normal (-s, 0, 1) / k
    fw0, fz0, fw1, fz1 = w_n, hn, w_b, 0.0                  # flank, neck end to floor end
    L = math.hypot(fw1 - fw0, fz1 - fz0)
    qw, qz = (fw1 - fw0) / L, (fz1 - fz0) / L
    slope = (fw1 - fw0) / (fz1 - fz0)                       # dw/dz along the flank (< 0)
    delta = -c * math.sqrt(1 + (s * qw) ** 2) / abs(qz)     # shift along w toward the groove interior
    z_n = fz0 + (w_hi - fw0 - delta) / slope
    w_bb = fw0 + (z_bot - fz0) * slope + delta
    dov = [(w_lo, z_bot), (w_bb, z_bot), (w_hi, z_n), (w_hi, H), (w_lo, H)]
    rel = [(w_lo, z_bot), (w_b - c, z_bot), (w_b - c, H), (w_lo, H)]
    return dov, rel


def _tongue_prism(fr: dict, u_end: float | None = None) -> Solid:
    RF = fr["RF"]
    dov, rel = tongue_section(RF)
    u_a = RF["u_entry"] - DESIGN["u_ahead"]
    u_end = RF["u_stop"] if u_end is None else u_end
    t = SH._prism_uwh(RF, dov, u_a, u_end).fuse(SH._prism_uwh(RF, rel, u_a, RF["u_relief"])).clean()
    sols = t.solids()
    return max(sols, key=lambda x: x.volume) if len(sols) > 1 else sols[0]


def _tongue_seat(fr: dict) -> Solid:
    """Magnet seat in the tongue / strip bottom (4.1 x 2.1 along the floor, 2.05 deep into the part)."""
    RF = fr["RF"]
    return SH._magnet_seat_local(RF, h_floor=DESIGN["clear_floor"] * RF["k"], up=True)


def _floor_frame(fr: dict, h_floor: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(forward-bottom corner of the 4.1 x 2.1 seat at the plane h_floor, t, w, n_f)."""
    RF = fr["RF"]
    M = R["magnet"]
    wm = 0.5 * (RF["w0"] + RF["w_n"])
    org = RF["Q"](RF["u_mag"][1], wm - M["seat"][1] / 2, h_floor) - M["seat"][0] * RF["t"]
    return org, RF["t"], RF["F"]["w"], RF["n_f"]


def magnet_block(fr: dict, side: str) -> Solid:
    """The 4 x 2 x 2 block in its seat: side 'floor' (shell: face recess below the floor) or 'tongue'
    (fin / strip: face recess above the tongue bottom), centred in the 4.1 x 2.1 seat."""
    M = R["magnet"]
    L, Wd, Hh = M["size"]
    if side == "floor":
        org, t, w, n = _floor_frame(fr, 0.0)
        n0 = -(M["recess"] + Hh)
    else:
        org, t, w, n = _floor_frame(fr, DESIGN["clear_floor"] * fr["RF"]["k"])
        n0 = M["recess"]
    o = org + 0.5 * (M["seat"][0] - L) * t + 0.5 * (M["seat"][1] - Wd) * w + n0 * n
    return Solid.make_box(L, Wd, Hh, Plane(origin=Vector(*o), x_dir=Vector(*t), z_dir=Vector(*n)))


# --------------------------------------------------------------------------
# Solids
# --------------------------------------------------------------------------


def _tip_corner(P_c: np.ndarray, side: np.ndarray, m_in: np.ndarray, n_in: np.ndarray, R: float) -> tuple:
    """Corner of the planform at P_c between a side line (unit direction `side`, from the corner toward the
    root; m_in its normal pointing into the planform) and the tip line (n_in its normal pointing into the
    planform): an arc of radius R tangent to the side line, down to the point where its outward normal leans
    45 deg from the print vertical (the fin prints tip down, so the tip line is the bed), then a straight land
    perpendicular to the tip line down to it (review C4: the plain round overhung up to 74 deg over its last
    1.85 mm). Returns (T_side, arc_mid, P45, F_tip): tangent point on the side line, a mid-arc point, the
    45 deg point and the foot of the land on the tip line."""
    a_ = R * (1.0 - m_in @ n_in) / (side @ n_in)
    Cc = P_c + a_ * side + R * m_in

    def pt(phi):                                                  # phi 90: the side face; 0: straight down
        q = math.radians(phi)
        return Cc + R * (math.cos(q) * (-n_in) + math.sin(q) * (-m_in))
    phi = 90.0 - DESIGN["tip_land_deg"]                          # overhang of the round's last face = 90 - phi
    P45 = pt(phi)
    F_tip = P45 + ((P_c - P45) @ n_in) * n_in                    # along -n_in down to the tip line
    return pt(90.0), pt(0.5 * (90.0 + phi)), P45, F_tip


def planform(fr: dict, a_lo: float) -> Solid:
    """Prism across w (+-10) of the (u, a) planform: datum root line - margin up to the tip line
    (datum + span), tip corners rounded (R6 at the LE, R2.5 at the TE) down to the 45 deg point of the
    round, then a short land perpendicular to the tip line (prints tip down: no corner face beyond 45 deg)."""
    c, span, ar = fr["c"], F["span"], fr["a_root"]
    r_le, r_te = DESIGN["tip_corner_r"]
    P_le = np.array([-0.5, ar(-0.5) + span])
    P_te = np.array([c + 0.5, ar(c + 0.5) + span])
    t_hat = (P_te - P_le) / np.linalg.norm(P_te - P_le)
    n_in = np.array([t_hat[1], -t_hat[0]])                        # toward the root (-a)
    if n_in[1] > 0:
        n_in = -n_in
    T1, M1, A1, F1 = _tip_corner(P_le, np.array([0.0, -1.0]), np.array([1.0, 0.0]), n_in, r_le + 0.5)
    T2, M2, A2, F2 = _tip_corner(P_te, np.array([0.0, -1.0]), np.array([-1.0, 0.0]), n_in, r_te + 0.5)
    Fr = fr["Fr"]
    o = Fr["o"] - 10.0 * Fr["w"]
    pl = Plane(origin=Vector(*o), x_dir=Vector(*Fr["u"]), z_dir=Vector(*(-Fr["w"])))   # local (X, Y) = (u, a)
    V = lambda q: pl.from_local_coords(Vector(float(q[0]), float(q[1]), 0.0))          # noqa: E731
    p0, p1 = np.array([-0.5, a_lo]), np.array([c + 0.5, a_lo])
    edges = [Edge.make_line(V(p0), V(p1)), Edge.make_line(V(p1), V(T2)),
             Edge.make_three_point_arc(V(T2), V(M2), V(A2)), Edge.make_line(V(A2), V(F2)),
             Edge.make_line(V(F2), V(F1)), Edge.make_line(V(F1), V(A1)),
             Edge.make_three_point_arc(V(A1), V(M1), V(T1)), Edge.make_line(V(T1), V(p0))]
    face = Face(Wire(edges))
    return Solid.extrude(face, Vector(*(20.0 * Fr["w"])))


def _root_flare(fr: dict, n_sec: int = 15) -> Solid:
    """Review C4 / P3: over the slot the tongue's neck (4.7 wide) is wider than the airfoil root above it, so
    the neck's top (the root face, root_gap over the skin) showed as 0.4 to 1.75 mm ledges beside the airfoil,
    facing the bed when the fin prints tip down. This flare fills them: at each chord station over the tongue,
    a 45 deg flank from each neck edge at the root face up and in to the airfoil skin (a root fillet, up to
    ~2 mm tall at the LE end), as a ruled loft of sections whose bottom dips 0.08 into the tongue (fused with
    it, never below the skin; a boolean cut with the hull offset returned nothing), from 0.2 behind the slot
    entry to 0.2 ahead of the end stop (the shell is under the root face beyond them)."""
    RF = fr["RF"]
    c = fr["c"]
    dov, _ = tongue_section(RF)
    w_lo, w_hi = dov[0][0], dov[2][0]
    up, lo, info = clarky(DESIGN["n_section"], te=P.MIN_TE / c)
    w_c = -0.5 * (info["y_min"] + info["y_max"]) * c
    g = DESIGN["root_gap"]
    k45 = 1.0 / math.tan(math.radians(DESIGN["root_flare_deg"]))
    u0, u1 = max(RF["u_entry"] + 0.2, 0.25), RF["u_stop"] - 0.2
    wires = []
    for u in np.linspace(u0, u1, n_sec):
        wu = w_c + c * float(np.interp(u / c, up[:, 0], up[:, 1]))
        wl = w_c + c * float(np.interp(u / c, lo[:, 0], lo[:, 1]))
        a_p, a_m = RF["a_hull"](u, w_hi, g), RF["a_hull"](u, w_lo, g)
        hu, hl = max(w_hi - wu, 0.05) * k45, max(wl - w_lo, 0.05) * k45
        e = 0.08                                       # flank foot just beyond the neck edge (inside the 0.1 root
        sec = [(w_lo + 0.05, a_m - 0.08), (w_hi - 0.05, a_p - 0.08), (w_hi + e, a_p), (wu, a_p + hu + e),  # gap):
               (wu - 0.05, a_p + hu + e + 0.3), (wl + 0.05, a_m + hl + e + 0.3), (wl, a_m + hl + e),       # no root-
               (w_lo - e, a_m)]                                                                          # face sliver
        wires.append(Wire.make_polygon([fr["W"](u, w, a) for w, a in sec], close=True))
    return Solid.make_loft(wires, ruled=True)


def fin0() -> tuple[Solid, dict]:
    """The theta-0 fin (world): airfoil cut by the hull offset (root on the skin), tongue filling the
    slot up to the root, a 45 deg root flare over the slot (review C4 / P3), magnet seat in the tongue
    bottom."""
    fr = frame()
    c, span = fr["c"], F["span"]
    a_lo = min(fr["a_root"](0.0), fr["a_root"](c)) - 14.0
    a_hi = max(fr["a_root"](0.0), fr["a_root"](c)) + span + 6.0
    face, info = section_face(fr, a_lo)
    body = Solid.extrude(face, Vector(0.0, 0.0, a_hi - a_lo))
    body = SH.inter(body, planform(fr, a_lo - 1.0))
    hull_g = _hull_offset(DESIGN["root_gap"])
    body = body.cut(hull_g)
    tg = SH.inter(_tongue_prism(fr), hull_g).cut(_relief_chamfer(fr))
    fin = body.fuse(tg).fuse(_root_flare(fr)).clean()
    fin = fin.cut(_tongue_seat(fr))
    sols = fin.solids()
    extra = sum(s_.volume for s_ in sols) - max(s_.volume for s_ in sols)
    fin = max(sols, key=lambda s_: s_.volume)
    info.update(n_solids_raw=len(sols), extra_mm3=round(extra, 4))
    return fin, info


def _relief_chamfer(fr: dict) -> Solid:
    """Fin only: in the slot's entry relief the tongue is wider than the airfoil nose above it, which left a
    2-2.6 mm ledge facing the bed when the fin prints tip down. Remove the relief part beyond the tongue's
    neck wall (w > w_hi) above a plane through the root face at w_hi that falls 1.2 per 1 toward +w (50 deg
    from the root face): the ledge becomes a 40 deg face; a V notch then shows beside the fin's LE root in
    the slot entry (the blanking strip keeps the full width)."""
    RF = fr["RF"]
    Fr = RF["F"]
    dov, rel = tongue_section(RF)
    w_hi = dov[2][0]
    g = DESIGN["root_gap"]
    ez = np.array([0.0, 0.0, 1.0])

    def P(u):
        q = Fr["o"] + u * Fr["u"] + w_hi * Fr["w"]
        return np.array([q[0], q[1], RF["a_hull"](u, w_hi, g)])
    p1, p2 = P(RF["u_entry"]), P(RF["u_relief"])
    d1 = p2 - p1
    d2 = Fr["w"] - 1.2 * ez
    nc = np.cross(d1, d2)
    if nc @ (Fr["w"] + ez) < 0:
        nc = -nc                                                # toward +w / +a: the side removed
    nc /= np.linalg.norm(nc)
    hs = SH.halfspace(tuple(p1), tuple(-nc))                    # the region on the +nc side
    zone = SH._prism_uwh(RF, [(w_hi, -1.0), (RF["w_b"] + 1.0, -1.0), (RF["w_b"] + 1.0, 20.0), (w_hi, 20.0)],
                         RF["u_entry"] - DESIGN["u_ahead"], RF["u_relief"] + 1.0)
    return SH.inter(zone, hs)


def strip_trim_u(fr: dict) -> float:
    """Chord station where the strip's entry wedge (hull over the floor) reaches strip_tip_min at its
    thinnest point across the strip (review P5: DESIGN strip_tip_min was declared but not applied)."""
    RF = fr["RF"]
    dov, rel = tongue_section(RF)
    z_bot = dov[0][1]
    ws = np.linspace(dov[0][0], rel[1][0], 9)
    for u in np.arange(RF["u_entry"] - 2.0, RF["u_stop"], 0.02):
        th = min(RF["a_hull"](u, w) - (RF["a_f"](u) + z_bot) for w in ws)
        if th >= DESIGN["strip_tip_min"]:
            return float(u)
    raise RuntimeError("strip never reaches strip_tip_min")


def _strip_end_cutters(fr: dict) -> list:
    """(1) the entry wedge ahead of strip_trim_u (a plane normal to the chord), with a strip_tip_break 45 deg
    break on its top edge; (2) the fingernail groove across the outer face near the aft end (review C6): its
    forward wall normal to the skin is what a nail pulls the strip forward by, its aft flank is 45 deg."""
    RF = fr["RF"]
    Fr = RF["F"]
    u_t = strip_trim_u(fr)
    b = DESIGN["strip_tip_break"]
    ez = np.array([0.0, 0.0, 1.0])
    p_t = fr["W"](u_t, 0.0, 0.0)
    cut_end = SH.halfspace((p_t.X, p_t.Y, p_t.Z), tuple(Fr["u"]))             # the -u side: u < u_t
    a_top = RF["a_hull"](u_t + b, 0.0)
    q = fr["W"](u_t + b, 0.0, a_top)
    nb = (-Fr["u"] + ez) / math.sqrt(2.0)
    cut_break = SH.halfspace((q.X, q.Y, q.Z), tuple(-nb))                      # ahead of / above the 45 deg plane
    G = DESIGN["nail_groove"]
    u_f = RF["u_stop"] - G["from_stop"] - G["w"] / 2                           # forward wall
    a_s = RF["a_hull"](u_f + G["w"] / 2, 0.0)
    d = G["depth"]
    prof = [(u_f, a_s + 3.0), (u_f, a_s - d), (u_f + G["w"] - d, a_s - d), (u_f + G["w"], a_s), (u_f + G["w"], a_s + 3.0)]
    o = -6.0
    face = Face(Wire.make_polygon([fr["W"](u, o, a) for u, a in prof], close=True))
    groove = Solid.extrude(face, Vector(*(12.0 * Fr["w"])))
    return [cut_end, cut_break, groove]


def strip0() -> Solid:
    """The theta-0 blanking strip (world): the tongue section from the slot entry to the stop, its outer
    face ON the hull (flush with the tongue bottom 0.05 over the floor), the magnet seat; its entry wedge
    trimmed where thinner than strip_tip_min with a 45 deg break (P5), a fingernail groove near its aft end
    to drag it forward out of the slot (C6)."""
    fr = frame()
    st = SH.inter(_tongue_prism(fr), _hull_offset(0.0))
    st = st.cut(_tongue_seat(fr))
    for cut in _strip_end_cutters(fr):
        st = st.cut(cut)
    sols = st.solids()
    return max(sols, key=lambda s_: s_.volume)


def print_fin(fin: Solid, fr: dict) -> Solid:
    """Tip down: the tip face (plane through the straight tip line and w) on the bed, span vertical."""
    Fr = fr["Fr"]
    d = Fr["u"] + fr["k"] * fr["ci"] * np.array([0.0, 0.0, 1.0])
    d /= np.linalg.norm(d)
    ez = np.array([0.0, 0.0, 1.0])
    n_tip = ez - (ez @ d) * d
    n_tip /= np.linalg.norm(n_tip)                          # outward (+a) normal of the tip face
    pl = Plane(origin=(0, 0, 0), x_dir=Vector(*d), z_dir=Vector(*(-n_tip)))
    return C_on_bed(pl.to_local_coords(fin))


def print_strip(st: Solid, fr: dict) -> Solid:
    """Floor face on the bed (the magnet seat opens on the bed), outer face up."""
    RF = fr["RF"]
    pl = Plane(origin=(0, 0, 0), x_dir=Vector(*RF["t"]), z_dir=Vector(*RF["n_f"]))
    return C_on_bed(pl.to_local_coords(st))


def C_on_bed(s):
    bb = s.bounding_box()
    return s.moved(Location((-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z)))


HEAL = dict(mode="none")


def _gprop_volumes(s: Solid) -> tuple[float, float]:
    """(default Gauss rule, adaptive 1e-6) B-rep volumes."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(s.wrapped, g)
    g2 = GProp_GProps()
    BRepGProp.VolumeProperties_s(s.wrapped, g2, 1e-6, False)
    return float(g.Mass()), float(g2.Mass())


_CACHE: dict = {}


def _built() -> dict:
    if "b" not in _CACHE:
        fr = frame()
        f0, info = fin0()
        s0 = strip0()
        _CACHE["b"] = dict(fr=fr, fin=f0, info=info, strip=s0, mag_floor=magnet_block(fr, "floor"),
                           mag_tongue=magnet_block(fr, "tongue"))
    return _CACHE["b"]


def _explode_vec(fr: dict, th: float, back: float, out: float) -> list:
    v = -back * fr["RF"]["t"] + out * fr["RF"]["n_f"]
    t = math.radians(th)
    y, z = v[1], v[2]                                       # rotate (y, z) by -th about X (at_azimuth)
    ct, st_ = math.cos(-t), math.sin(-t)
    return [float(v[0]), float(ct * y - st_ * z), float(st_ * y + ct * z)]


def parts() -> list[dict]:
    if "parts" in _CACHE:
        return _CACHE["parts"]
    B = _built()
    fr = B["fr"]
    pf, ps = print_fin(B["fin"], fr), print_strip(B["strip"], fr)
    out = []
    ths = list(F["theta_deg"])
    for k, th in enumerate(ths):
        out.append(dict(
            name=f"fin_{k}", shape=at_azimuth(B["fin"], th), print=pf, material="shell", role="body", rotor=False,
            qty=1, bought=False, explode=_explode_vec(fr, th, 22.0, 22.0), theta_deg=th, packed=False,
            notes="PETG, shell colour. Identical part x 2 (theta %g and %g). Prints TIP DOWN (span vertical, the "
                  "root tongue on top), no supports: the tongue's undercut flank and the root flare beside the "
                  "neck are 45 deg, the tip corners end on a short land at their 45 deg point; slender (73.7 mm on "
                  "a 3.6 mm wide first layer): 5 mm brim all round, outer walls <= 50 mm/s, 6-8 s minimum layer "
                  "time, and on a bed-slinger (A1) set the chord along the bed's motion axis (Y). Glue the "
                  "4 x 2 x 2 N52 block into the tongue-bottom seat (attracting the floor magnet). Fit: hook a "
                  "fingernail in the blanking strip's groove and drag it forward out of the slot, set the tongue "
                  "in at the slot entry and push the fin aft and in to the stop; it clicks on the magnet pair. "
                  "Yaw play in the groove about +-1.5 deg (0.25 mm side clearance over the 12.1 mm tongue): "
                  "irrelevant for the display build, a v2 item (couple ratio)." % tuple(ths)))
    for k, th in enumerate(ths):
        out.append(dict(
            name=f"fin_blank_{k}", shape=at_azimuth(B["strip"], th), print=ps, material="shell", role="body",
            rotor=False, qty=1, bought=False, explode=_explode_vec(fr, th, 16.0, 12.0), theta_deg=th, packed=True,
            notes="Blanking strip for the packed egg (fins off): the fin tongue's section with its outer face on "
                  "the hull profile, flush when seated; slides in like the fin and is held by the same magnet "
                  "pair; to remove it, hook a fingernail in the groove across its outer face near the aft end and "
                  "drag it forward. PLA/PETG in the shell colour; prints on its floor face (the magnet seat opens "
                  "on the bed), no supports; its flush face slopes ~29 deg, so print it with 0.08 to 0.12 mm "
                  "layers (or variable layer height) and no ironing to avoid terraces; the entry wedge is trimmed "
                  "at 0.6 mm."))
    mag_note = ("N52 block 4 x 2 x 2 magnetised through 2 (VERIFY), CA-glued in the "
                "4.1 x 2.1 x 2.05 seat, face %.2f mm %s; the floor / tongue pair attracts" % (R["magnet"]["recess"], "%s"))
    for k, th in enumerate(ths):
        out.append(dict(name=f"rail_magnet_{k}", shape=at_azimuth(B["mag_floor"], th), print=None, material="magnet",
                        role="body", rotor=False, qty=1, bought=True,
                        explode=_explode_vec(fr, th, 0.0, -6.0), theta_deg=th,
                        notes=mag_note % "below the groove floor (front shell)"))
        out.append(dict(name=f"fin_magnet_{k}", shape=at_azimuth(B["mag_tongue"], th), print=None, material="magnet",
                        role="body", rotor=False, qty=1, bought=True, packed=False,
                        explode=_explode_vec(fr, th, 22.0, 18.0), theta_deg=th,
                        notes=mag_note % "above the fin tongue's bottom"))
        out.append(dict(name=f"fin_blank_magnet_{k}", shape=at_azimuth(B["mag_tongue"], th), print=None,
                        material="magnet", role="body", rotor=False, qty=1, bought=True, packed=True,
                        explode=_explode_vec(fr, th, 16.0, 8.0), theta_deg=th,
                        notes=mag_note % "above the blanking strip's bottom"))
    _CACHE["parts"] = out
    return out


# --------------------------------------------------------------------------
# Magnet pair force (the layout's calibrated charge model, extended to rectangular poles)
# --------------------------------------------------------------------------


def _omega_rect(X, Y, z, a, b):
    """Solid angle of the rectangle [-a/2, a/2] x [-b/2, b/2] (plane 0) seen from (X, Y, z), z > 0."""
    def f(x, y):
        return np.arctan2(x * y, z * np.sqrt(x * x + y * y + z * z))
    x1, x2 = -a / 2 - X, a / 2 - X
    y1, y2 = -b / 2 - Y, b / 2 - Y
    return f(x2, y2) - f(x1, y2) - f(x2, y1) + f(x1, y1)


def _phi(z: float, a: float, b: float, dx: float = 0.0, n: int = 241) -> float:
    """Integral over the rectangle a x b (shifted by dx along a) of the solid angle of a second, coaxial
    a x b rectangle at distance z (mm^2 sr)."""
    z = max(z, 1e-4)
    xs = np.linspace(-a / 2, a / 2, n) + dx
    ys = np.linspace(-b / 2, b / 2, max(41, int(n * b / a) | 1))
    Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
    om = _omega_rect(Xg, Yg, z, a, b)
    return float(np.trapezoid(np.trapezoid(om, ys, axis=1), xs))


def block_pair_force(gap: float, dx: float = 0.0, size=(4.0, 2.0, 2.0), Br: float = 1.43) -> float:
    """Attraction (N, normal to the faces) between two identical blocks magnetised through size[2], pole
    faces size[0] x size[1] parallel at gap (mm), shifted dx along size[0]: surface-charge model
    F = Br^2 / (4 pi mu0) [Phi(g) - 2 Phi(g + h) + Phi(g + 2h)], times the layout's knock-down k
    (data/layout_rotor.json detents.calibration, 4 x 2 catalogue pull / model)."""
    a, b, h = size
    mu0 = 4e-7 * math.pi
    k = json.loads((ROOT / "data" / "layout_rotor.json").read_text())["detents"]["calibration"]["k_design"]
    s_ = _phi(gap, a, b, dx) - 2 * _phi(gap + h, a, b, dx) + _phi(gap + 2 * h, a, b, dx)
    return k * Br ** 2 / (4 * math.pi * mu0) * s_ * 1e-6          # mm^2 -> m^2


def block_pair_shear(gap: float, dx: float, size=(4.0, 2.0, 2.0)) -> float:
    """Restoring force (N) along size[0] at offset dx: -dU/dx with U(g, dx) = -int_g^inf F(g', dx) dg'."""
    gg = np.concatenate([np.linspace(gap, gap + 2.0, 21), np.geomspace(gap + 2.1, gap + 60.0, 25)])
    eps = 0.05

    def U(d):
        Fz = np.array([block_pair_force(g, d, size) for g in gg])
        return -float(np.trapezoid(Fz, gg))
    return -(U(dx + eps) - U(dx - eps)) / (2 * eps)


def disc_force_check() -> dict:
    """The rectangle model with a 4 x 4 square vs the layout's 4 x 2 disc table at 0.4 mm (sanity)."""
    t = json.loads((ROOT / "data" / "layout_rotor.json").read_text())["detents"]["force_vs_gap"]
    return dict(square_4x4x2_at_0p4_N=round(block_pair_force(0.4, 0.0, (4.0, 4.0, 2.0)), 3),
                disc_4x2_at_0p4_N=float(np.interp(0.4, t["gap_mm"], t["force_N"])),
                note="a 4 x 4 square has 1.27 x the disc's pole area")


def _root_stress(lay: dict) -> dict:
    """Bending stress at the fin root from the layout's tow-case root moment, on the Clark Y section
    (second moment about the chord-wise axis), against the PETG interlayer strength (ESTIMATE 18 MPa):
    printed tip down, this stress runs across the layers."""
    up, lo, _ = clarky(400, te=P.MIN_TE / F["chord"])
    c = F["chord"]
    x = np.linspace(0, 1, 400)
    yu = np.interp(x, up[:, 0], up[:, 1]) * c
    yl = np.interp(x, lo[:, 0], lo[:, 1]) * c
    dx = c / 399
    A = float(np.sum(yu - yl) * dx)
    yc = float(np.sum((yu ** 2 - yl ** 2) / 2) * dx / A)
    I = float(np.sum(((yu - yc) ** 3 - (yl - yc) ** 3) / 3) * dx)
    ymax = float(max(yu.max() - yc, yc - yl.min()))
    M = float(lay["root_moment_tow_Nm"]) * 1e3
    sig = M * ymax / I
    return dict(M_Nmm=round(M, 2), area_mm2=round(A, 2), I_mm4=round(I, 2), y_max_mm=round(ymax, 3),
                sigma_MPa=round(sig, 3), interlayer_MPa=18.0, utilisation=round(sig / 18.0, 4))


# --------------------------------------------------------------------------
# Previews
# --------------------------------------------------------------------------


def _pv(dicts):
    return [{k: v for k, v in p.items() if k in ("name", "shape", "material", "role")} for p in dicts]


def preview_parts() -> list[dict]:
    return _pv([p for p in parts() if p["name"].startswith(("fin_0", "fin_1"))])


def _front_shell():
    import os
    from build123d import import_brep
    f = os.environ.get("SHELLS_FRONT_BREP")
    return import_brep(f) if f else SH.build_all()["front_shell"]


def preview_on_shell() -> list[dict]:
    """Deployed look: the front shell with both fins and their magnets fitted."""
    out = [dict(name="front_shell", shape=_front_shell(), material="shell", role="body")]
    return out + _pv([p for p in parts() if p["name"] in ("fin_0", "fin_1")])


def preview_packed_on_shell() -> list[dict]:
    """Packed look: the front shell with the blanking strips in the slots (fins off)."""
    out = [dict(name="front_shell", shape=_front_shell(), material="shell", role="body")]
    return out + _pv([p for p in parts() if p["name"] in ("fin_blank_0", "fin_blank_1")])


def preview_rail_close() -> list[dict]:
    """Close-up of the theta FIN[0] rail: front shell clipped around it, fin_0 pulled 6 mm back along the
    rail (half inserted) with its magnet, the rail magnet in the floor."""
    B = _built()
    fr = B["fr"]
    th = F["theta_deg"][0]
    box = Solid.make_box(40.0, 30.0, 40.0).moved(Location((18.0, 16.0, -20.0)))
    sh = SH.inter(_front_shell(), box)
    back = at_azimuth(Solid.make_box(1, 1, 1).moved(Location(tuple(-6.0 * fr["RF"]["rv"]))), th).center() - \
        at_azimuth(Solid.make_box(1, 1, 1), th).center()
    fin = at_azimuth(B["fin"], th).moved(Location(back))
    mag = at_azimuth(B["mag_tongue"], th).moved(Location(back))
    return [dict(name="shell_clip", shape=sh, material="shell", role="body"),
            dict(name="fin_half_in", shape=fin, material="fin", role="body"),
            dict(name="fin_mag", shape=mag, material="magnet", role="body"),
            dict(name="rail_mag", shape=at_azimuth(B["mag_floor"], th), material="magnet", role="body")]


def preview_print() -> list[dict]:
    ps = parts()
    f = [p for p in ps if p["name"] == "fin_0"][0]["print"]
    s = [p for p in ps if p["name"] == "fin_blank_0"][0]["print"]
    return [dict(name="fin_print", shape=f, material="shell", role="body"),
            dict(name="strip_print", shape=s.moved(Location((45.0, 0.0, 0.0))), material="blade", role="body")]


def render_previews(views) -> list[str]:
    WORK.mkdir(parents=True, exist_ok=True)
    out = []
    for fn, stem, q in views:
        png = WORK / f"{stem}.png"
        r = subprocess.run([sys.executable, str(ROOT / "scripts/preview.py"), f"ovo.cad.fins:{fn}", "--out", str(png),
                            "--q", q], cwd=ROOT, capture_output=True, text=True)
        if r.returncode == 0:
            out.append(str(png.relative_to(ROOT)))
    return out


# --------------------------------------------------------------------------
# Checks (measured on the solids)
# --------------------------------------------------------------------------


def _mesh(shape, tol=0.03):
    import trimesh
    v, f = shape.tessellate(tol, 0.2)
    m = trimesh.Trimesh(np.array([(p.X, p.Y, p.Z) for p in v]), np.array(f), process=True)
    m.merge_vertices()
    return m


def _manifold(shape, tol=0.03):
    import manifold3d as m3
    from .common import closed_mesh
    m = closed_mesh(_mesh(shape, tol))
    return m3.Manifold(m3.Mesh(vert_properties=np.asarray(m.vertices, np.float32),
                               tri_verts=np.asarray(m.faces, np.uint32)))


def _overlap_gap(ma, mb, search=3.0) -> tuple[float, float]:
    ov = float((ma ^ mb).volume())
    return ov, (0.0 if ov > 1e-3 else float(ma.min_gap(mb, search)))


def _ray_walls(shape, origins, dirs) -> np.ndarray:
    return SH._occ_wall(shape, origins, dirs)


def seat_ligaments(part: Solid, fr: dict, h_plane: float) -> dict:
    """Material around a 4.1 x 2.1 x 2.05 tongue / strip magnet seat, MEASURED with OCC rays from inside
    the seat void: across w (both sides), along the floor (forward and aft walls) and up into the part
    (n_f), on a grid of starting points; the first material run along each ray is the ligament."""
    M = R["magnet"]
    org, t, w, n = _floor_frame(fr, h_plane)
    L, Wd, D = M["seat"][0], M["seat"][1], M["depth"]
    out = {}
    for key, d_, pts in (
            ("side_-w", -w, [org + a * t + 0.05 * w + h * n for a in np.linspace(0.3, L - 0.3, 7)
                            for h in np.linspace(0.2, D - 0.2, 5)]),
            ("side_+w", w, [org + a * t + (Wd - 0.05) * w + h * n for a in np.linspace(0.3, L - 0.3, 7)
                           for h in np.linspace(0.2, D - 0.2, 5)]),
            ("aft_wall", t, [org + (L - 0.05) * t + b * w + h * n for b in np.linspace(0.3, Wd - 0.3, 4)
                            for h in np.linspace(0.2, D - 0.2, 5)]),
            ("fwd_wall", -t, [org + 0.05 * t + b * w + h * n for b in np.linspace(0.3, Wd - 0.3, 4)
                             for h in np.linspace(0.2, D - 0.2, 5)]),
            ("ceiling", n, [org + a * t + b * w + (D - 0.05) * n for a in np.linspace(0.3, L - 0.3, 7)
                           for b in np.linspace(0.3, Wd - 0.3, 4)])):
        d = _ray_walls(part, pts, [d_] * len(pts))
        out[key] = round(float(np.min(d)), 3)
    out["min_mm"] = min(v for v in out.values())
    return out


def slide_sweep(m_part, m_shell, fr: dict, th: float, d_max: float, step: float = 0.5) -> dict:
    """Pull the part back along the rail (theta th) from seated to d_max and record the largest overlap
    with the front shell and the smallest gap: the insertion path is free if every overlap is 0."""
    rv = fr["RF"]["rv"]
    t = math.radians(th)
    ct, st_ = math.cos(-t), math.sin(-t)
    worst, gmin, at = 0.0, math.inf, None
    for d in np.arange(0.0, d_max + 1e-9, step):
        v = -d * rv
        vw = (float(v[0]), float(ct * v[1] - st_ * v[2]), float(st_ * v[1] + ct * v[2]))
        ov, gp = _overlap_gap(m_part.translate(vw), m_shell, 2.0)
        if ov > worst:
            worst, at = ov, float(d)
        gmin = min(gmin, gp)
    return dict(max_overlap_mm3=round(worst, 4), at_pull_mm=at, min_gap_mm=round(gmin, 4), d_max_mm=round(d_max, 2),
                step_mm=step)


def checks() -> dict:
    ps = parts()
    B = _built()
    fr, RF = B["fr"], B["fr"]["RF"]
    res = dict(section=B["info"], parts=[], interfaces=[], fin_block_vs_params=SH.FIN_VS_PARAMS,
               rail=dict(u_entry=RF["u_entry"], u_relief=RF["u_relief"], u_stop=RF["u_stop"], u_magnet=RF["u_mag"],
                         slope=RF["s"], design=DESIGN))
    up, lo, _ = clarky(DESIGN["n_section"], te=P.MIN_TE / F["chord"])
    cmp_ = compare_uiuc(up, lo)
    res["section_vs_uiuc"] = cmp_
    lay = json.loads((ROOT / "data" / "layout_internal.json").read_text())["fins"]
    try:
        AR = lay["span_eff_mm"] / F["chord"]
        a3 = lambda a0: a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))   # noqa: E731
        nfu, nfc = cmp_["neuralfoil_Re55k"]["uiuc"], cmp_["neuralfoil_Re55k"]["cad"]
        f = (a3(nfc["a0_per_rad"]) * (F["incidence_deg"] - nfc["alpha_L0_deg"])) / (
            a3(nfu["a0_per_rad"]) * (F["incidence_deg"] - nfu["alpha_L0_deg"]))
        r0 = lay["torque_basis"]["ratio"]
        res["couple_ratio_as_built"] = dict(layout_ratio=round(r0, 4), factor=round(f, 4), ratio=round(r0 * f, 4),
                                            band=[1.10, 1.155], ok=1.10 <= r0 * f <= 1.155,
                                            note="layout ratio (UIUC section) x NeuralFoil lift factor of the CAD "
                                                 "section (MIN_TE trailing edge) at Re 55 k, n_crit 5")
    except Exception as ex:  # noqa: BLE001
        res["couple_ratio_as_built"] = f"not computed: {ex}"
    I = res["interfaces"]

    def add(name, value, ok, **kw):
        I.append(dict(check=name, value=value, ok=bool(ok), **kw))
        print(("PASS " if ok else "FAIL ") + name + ": " + str(value), flush=True)

    for p in ps:
        c = check(p["shape"], p["name"])
        if p["print"] is not None:
            m = _mesh(p["print"])
            n = m.face_normals
            zc = m.triangles_center[:, 2]
            ovh = float(m.area_faces[(n[:, 2] < -math.cos(math.radians(45.0)) - 1e-6) & (zc > 0.05)].sum())
            c.update(print_size_mm=[round(v, 2) for v in p["print"].bounding_box().size],
                     tessellation_watertight=bool(m.is_watertight), fits_bed=fits_bed(p["print"]),
                     overhang_area_mm2_gt45=round(ovh, 3),
                     bed_contact_mm2=round(float(m.area_faces[(n[:, 2] < -0.999) & (zc < 0.05)].sum()), 1))
        res["parts"].append(c)
    fin = B["fin"]
    # geometry as built (theta 0 world): incidence, chord, span
    Fr = fr["Fr"]
    wc = B["info"]["chord_line_w"]
    inside = lambda u, w, a: fin.is_inside(fr["W"](u, w, a))            # noqa: E731
    a_mid = fr["a_root"](fr["c"] / 2) + F["span"] / 2
    le_ok = inside(0.05, wc, a_mid) and not inside(-0.05, wc, a_mid)
    te_ok = inside(fr["c"] - 0.05, wc, a_mid) and not inside(fr["c"] + 0.05, wc, a_mid)
    add("chord measured at mid span (LE / TE probes +-0.05 mm)", fr["c"], le_ok and te_ok)
    inc = math.degrees(math.atan2(-Fr["u"][1], Fr["u"][0]))
    add("fin incidence (chord vs flow, deg, TE toward -theta)", round(inc, 4), abs(inc - F["incidence_deg"]) < 1e-3)
    u_m = fr["c"] / 2
    a_tip = fr["a_root"](u_m) + F["span"]
    tip_ok = inside(u_m, wc + 0.5, a_tip - 0.1) and not inside(u_m, wc + 0.5, a_tip + 0.05)
    add("span measured at mid chord (datum root line to the tip, probes)", F["span"], tip_ok)
    # root on the skin: root face height over the hull along the chord (probes at w = chord line)
    gaps = []
    for u in np.linspace(1.0, fr["c"] - 1.0, 25):
        if RF["u_entry"] - 0.5 <= u <= RF["u_stop"] + 0.5:
            continue                                   # over the slot the keel continues down
        a_h = RF["a_hull"](u, wc)
        lo_, hi_ = a_h - 0.5, a_h + 2.0
        for _ in range(40):
            mid = 0.5 * (lo_ + hi_)
            lo_, hi_ = (lo_, mid) if inside(u, wc, mid) else (mid, hi_)
        gaps.append(hi_ - a_h)
    add("fin root face over the hull skin outside the slot (min / max, mm; design root_gap %.2f)" % DESIGN["root_gap"],
        [round(min(gaps), 3), round(max(gaps), 3)], abs(min(gaps) - DESIGN["root_gap"]) < 0.03
        and abs(max(gaps) - DESIGN["root_gap"]) < 0.03)
    # interface with the real front shell (cached BREP; manifold3d meshes, 0.03 mm)
    try:
        front = _front_shell()
        th = F["theta_deg"][0]
        m_sh = _manifold(front, 0.03)
        m_fin = _manifold(at_azimuth(fin, th), 0.02)
        m_st = _manifold(at_azimuth(B["strip"], th), 0.02)
        ov, gp = _overlap_gap(m_fin, m_sh)
        add("fin_0 vs front shell, seated: overlap mm3 / min gap mm (end stop contact expected)", [round(ov, 4), round(gp, 4)],
            ov < 0.05)
        ov2, gp2 = _overlap_gap(m_st, m_sh)
        add("blanking strip vs front shell, seated: overlap mm3 / min gap mm", [round(ov2, 4), round(gp2, 4)], ov2 < 0.05)
        d_max = RF["u_stop"] - RF["u_entry"] + DESIGN["u_ahead"]
        sw = slide_sweep(m_fin, m_sh, fr, th, d_max)
        add("fin_0 pulled out along the rail (0.5 mm steps to %.1f): largest overlap with the shell mm3" % d_max,
            sw, sw["max_overlap_mm3"] < 0.05)
        sw2 = slide_sweep(m_st, m_sh, fr, th, d_max)
        add("blanking strip pulled out along the rail: largest overlap with the shell mm3", sw2, sw2["max_overlap_mm3"] < 0.05)
        # strip flush: outer-face vertices vs the hull radius
        ms = _mesh(at_azimuth(B["strip"], th), 0.02)
        v = np.asarray(ms.vertices)
        rr = np.hypot(v[:, 1], v[:, 2]) - profile.radius(np.clip(v[:, 0], 0, P.HULL_L))
        top = rr > -0.3
        add("blanking strip outer face vs the hull profile (max proud / mean of the outer vertices), mm",
            [round(float(rr.max()), 4), round(float(rr[top].mean()), 4)], rr.max() <= 0.02)
        # magnets: seated pair gap, overlaps with their seats
        mf = _manifold(at_azimuth(B["mag_floor"], th), 0.01)
        mt = _manifold(at_azimuth(B["mag_tongue"], th), 0.01)
        add("rail magnet in its floor seat: overlap with the shell mm3", round(_overlap_gap(mf, m_sh)[0], 4),
            _overlap_gap(mf, m_sh)[0] < 0.01)
        add("fin magnet in its tongue seat: overlap with the fin mm3", round(_overlap_gap(mt, m_fin)[0], 4),
            _overlap_gap(mt, m_fin)[0] < 0.01)
        g_model = float(mf.min_gap(mt, 2.0))
        g_seat = g_model - DESIGN["clear_floor"]
        Fz_m, Fz_s = block_pair_force(g_model), block_pair_force(g_seat)
        sh_ = [(d, block_pair_shear(g_seat, d)) for d in (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)]
        peak = max(-v_ for _, v_ in sh_)
        mass_g = res["parts"][0]["volume_cm3"] * 1.27 * 0.6
        add("fin magnet pair (4 x 2 x 2 N52): face gap model / seated (tongue on the floor), mm; normal force N; peak "
            "shear (slide-out) N; fin weight N",
            dict(gap_model=round(g_model, 3), gap_seated=round(g_seat, 3), normal_N=[round(Fz_m, 2), round(Fz_s, 2)],
                 shear_peak_N=round(peak, 3), shear_vs_offset=[[d, round(-v_, 3)] for d, v_ in sh_],
                 slide_out_N_est=round(peak + 0.3 * Fz_s, 2), fin_weight_N=round(mass_g * 9.81e-3, 3)),
            Fz_s >= 2.0 and peak + 0.3 * Fz_s >= 10 * mass_g * 9.81e-3,
            note="charge model x layout k; slide-out = peak shear + friction 0.3 x normal (ESTIMATE); fin mass ~ "
                 "60 % of the solid PETG volume (ESTIMATE)")
        res["magnet_model_check"] = disc_force_check()
    except Exception as ex:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        add("rail interface with the front shell", f"not run: {type(ex).__name__}: {ex}", False)
    # ligaments around the magnet seats (MEASURED with rays on the solids)
    lig_f = seat_ligaments(fin, fr, DESIGN["clear_floor"] * RF["k"])
    add("fin tongue magnet seat ligaments (rays: -w / +w / aft / fwd / up), min mm", lig_f, lig_f["min_mm"] >= P.MIN_WALL - 0.02)
    lig_s = seat_ligaments(B["strip"], fr, DESIGN["clear_floor"] * RF["k"])
    add("blanking strip magnet seat ligaments (rays), min mm", lig_s, lig_s["min_mm"] >= P.MIN_WALL - 0.02)
    # GProp vs mesh volume (finding T19), no face split (R1)
    gd, ga = _gprop_volumes(fin)
    mv = float(_mesh(fin, 0.01).volume)
    add("fin B-rep volume: default / adaptive GProp / 0.01 mm mesh (mm3), agree within 0.5 %",
        [round(gd, 1), round(ga, 1), round(mv, 1)], abs(gd - mv) <= 0.005 * mv and abs(ga - mv) <= 0.005 * mv,
        heal=HEAL)
    f0, f1 = (p["shape"] for p in ps if p["name"] in ("fin_0", "fin_1"))
    f0r = f0.rotate(Axis.X, -180.0)
    dv = abs(f0r.volume - f1.volume) + (f0r.bounding_box().center() - f1.bounding_box().center()).length
    add("fin_1 = fin_0 turned 180 deg about X (volume + bbox centre deviation)", round(dv, 5), dv < 1e-3)
    rmax = float(np.hypot(_mesh(f0, 0.1).vertices[:, 1], _mesh(f0, 0.1).vertices[:, 2]).max())
    add("fin tip radius max (mm) vs layout r_tip_max %.2f" % lay["r_tip_max_mm"], round(rmax, 3),
        abs(rmax - lay["r_tip_max_mm"]) < 1.5)
    res["structure"] = dict(note="layout root moment (tow case, 2.5 m/s) on the Clark Y root section; printed tip "
                                 "down, the span-wise bending runs across layers", **_root_stress(lay))
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
    stale = [OUT / "step" / "fin_cross_pins.step"]
    for p_ in stale:
        if p_.exists():
            p_.unlink()
    for p in ps:
        sp = OUT / "step" / f"{p['name']}.step"
        export_step(p["shape"], str(sp))
        files.append(str(sp.relative_to(ROOT)))
        if p["print"] is not None:
            st = OUT / "stl" / f"{p['name']}.stl"
            export_stl(p["print"], str(st), tolerance=0.02, angular_tolerance=0.15)
            from .common import repair_stl
            fix = repair_stl(st)
            files.append(str(st.relative_to(ROOT)))
            m = trimesh.load(st)
            nrm = m.face_normals
            adj = m.face_adjacency
            fold = int((np.einsum("ij,ij->i", nrm[adj[:, 0]], nrm[adj[:, 1]]) < -0.985).sum())
            for c in res["parts"]:
                if c["name"] == p["name"]:
                    c.update(stl_watertight=bool(m.is_watertight), stl_repair=fix, stl_folded_pairs_gt170=fold,
                             stl_volume_cm3=round(float(m.volume) / 1000, 3))
    res["files"] = files
    res["design"] = DESIGN
    if "--no-render" not in sys.argv:
        res["renders"] = render_previews([
            ("preview_on_shell", "fins_on_shell", "view=150,18&zoom=1.1"),
            ("preview_packed_on_shell", "fins_packed_on_shell", "view=150,18&zoom=1.1"),
            ("preview_rail_close", "fins_rail_close", "view=165,25&zoom=1.2"),
            ("preview_print", "fins_print", "view=-30,15")])
    res["elapsed_s"] = round(time.time() - t0, 1)
    (OUT / "checks").mkdir(parents=True, exist_ok=True)
    (OUT / "checks" / "fins.json").write_text(json.dumps(res, indent=1, default=str))
    for c in res["parts"]:
        print({k: c.get(k) for k in ("name", "valid", "n_solids", "volume_cm3", "stl_watertight", "overhang_area_mm2_gt45",
                                     "print_size_mm", "stl_folded_pairs_gt170")})
    print("elapsed", res["elapsed_s"])


if __name__ == "__main__":
    main()
