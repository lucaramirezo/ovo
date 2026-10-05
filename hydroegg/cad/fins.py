"""Roll fins of the HydroEgg: two identical Clark Y fins on yawed dovetail rails of the front
shell (FIN block of hydroegg/params.py).

    scripts/cadrun.sh .venv/bin/python -m hydroegg.cad.fins

Geometry (all from params.FIN; the rail section and the cross pin from
hydroegg.cad.shells.FIN_RAIL, the shell side of the interface):

* Fin frame (same as layout_internal.fin_frame and shells.fin_frame): at the fin angle theta,
  u runs along the yawed root chord from the LE, w across it (+ toward +theta), a radially out.
  The chord is yawed by the incidence (12 deg, TE toward -theta) about the radial line through
  the yaw pivot (0.15 c, x 22.24), so the root chord spans x 16.37 to 55.5. The root line (pad
  top) runs from r 33.79 at the LE to r 44.60 at the TE; the tip line is the root line + span.
* Section: Clark Y from the classic ordinates (NACA, percent chord above the flat base line),
  re-referenced to its chord (LE to mid-TE, the frame NeuralFoil and the layout study use for
  the 12 deg incidence), TE thickened symmetrically to MIN_TE. The section envelope is centred
  on the rail (w = 0) as the layout's fold-clearance model assumed.
* Root: a half-dovetail tongue = the shell groove inset by CLEAR_SLIDE on every face, from the
  LE station to the end stop (u = chord, the TE: drag pushes the fin home), with the 2.3 mm
  cross-pin hole at u 32 (2 mm pin through the pad). Where the section is wider than the
  tongue's neck the root flares at 45 deg; the fin clears the pad top by CLEAR_SLIDE.
* Print: ROOT DOWN (tongue bottom on the bed, span up, 15 deg lean in the fin plane). A flat
  print is not possible with a straight rail along the chord: the Clark Y flat bottom is 2.1 deg
  off its chord, so the tongue (parallel to the chord, as the rail must be) cannot share the
  bed plane with the flat side; rotating the section to make them coplanar would raise the
  incidence to 14.1 deg (couple 1.27 x the load-line torque, out of the 1.10-1.155 band).
  Upright printing also traces the airfoil with perimeters (no stair steps across the chord)
  and keeps every face within 45 deg; the root bending stress then runs across layers, still
  at a small fraction of the PETG interlayer strength (fins.json).
* Both fins are the same part (theta 93 and 273 = the first turned 180 deg about X).
"""

from __future__ import annotations

import json
import math
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
    export_step,
    export_stl,
)

from .. import params as P
from .common import at_azimuth, check, fits_bed, interference, min_distance

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out" / "cad"
F = P.FIN

# Shell side of the rail (shells.FIN_RAIL); fallback copy if shells.py cannot be imported.
_RAIL_FALLBACK = dict(pad_w=F["rail_w"], groove_depth=F["groove_depth"], groove_w_top=4.5, groove_neck_h=2.5,
                      groove_w_bottom=6.0, groove_w_offset=-2.25, end_stop_u=F["chord"], pin_d=2.0,
                      pin_hole_d=2.0, pin_l=10.0, pin_u=32.0, pin_depth=2.0, tongue_pin_hole_d=2.3)


def rail() -> dict:
    try:
        from .shells import FIN_RAIL as R
        R = dict(R)
        R["source"] = "hydroegg.cad.shells.FIN_RAIL"
    except Exception:  # noqa: BLE001
        R = dict(_RAIL_FALLBACK, source="fins fallback copy")
    return R


DESIGN = dict(
    root_gap=P.CLEAR_SLIDE,        # fin body above the pad top
    tongue_clear=P.CLEAR_SLIDE,    # tongue inset from the groove on every face
    tongue_lead=0.6,               # 45 deg lead-in on the tongue's front end
    tongue_bed_chamfer=0.3,        # 45 deg chamfer on the tongue's flat-wall bed edge (finding T18)
    tip_corner_r=(8.0, 3.0),       # planform rounding at the LE-tip / TE-tip corners
    tip_edge_fillet=1.2,           # rounding of the tip edge (tried, then smaller radii)
    root_flare=1.1,                # dz/dw of the root flare beside the tongue neck: >= 1.036 keeps it
                                   # within 45 deg when printed root down (the root line climbs 15 deg)
    n_section=90,                  # points per surface
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
# Fin frame and solid
# --------------------------------------------------------------------------


def frame() -> dict:
    """Yawed fin frame numbers (identical to shells.fin_frame / layout_internal.fin_frame)."""
    i = math.radians(F["incidence_deg"])
    sense = -1.0 if F["te_side"] == "-theta" else 1.0
    ci, si = math.cos(i), sense * math.sin(i)
    c = F["chord"]
    x_te = F["x_te"]
    x_le = x_te - c * ci
    pc = F["yaw_pivot_chord_frac"] * c
    x_p = x_le + pc * ci
    k = (F["root_r_te"] - F["root_r_le"]) / (x_te - x_le)
    return dict(ci=ci, si=si, c=c, x_le=x_le, x_te=x_te, x_p=x_p, pc=pc, k=k, lat_le=-pc * si,
                a_root=lambda u: F["root_r_le"] + k * ci * u)


def plane(theta: float = 0.0) -> Plane:
    """Local frame of the fin at azimuth theta: local X = u, local Y = w, local Z = a."""
    fr = frame()
    t = math.radians(theta)
    e_r = Vector(0, math.sin(t), math.cos(t))
    e_t = Vector(0, math.cos(t), -math.sin(t))
    o = Vector(fr["x_le"], 0, 0) + e_t * fr["lat_le"]
    u = Vector(fr["ci"], 0, 0) + e_t * fr["si"]
    return Plane(origin=o, x_dir=u, z_dir=e_r)


def _poly_face(pts3) -> Face:
    return Face(Wire.make_polygon([Vector(*p) for p in pts3], close=True))


def _rail_prism(sec_wz, u0: float, u1: float) -> Solid:
    """Prism along the rail (the root line direction) of a (w, z) section, z measured from the root
    line at w = 0; local coordinates."""
    fr = frame()
    a0 = fr["a_root"](u0)
    face = _poly_face([(u0, w, a0 + z) for w, z in sec_wz])
    return Solid.extrude(face, Vector(u1 - u0, 0, fr["k"] * fr["ci"] * (u1 - u0)))


def tongue_section(R: dict | None = None) -> list[tuple[float, float]]:
    """The groove section (shells.FIN_RAIL) inset by CLEAR_SLIDE: flat wall, neck wall, 45 deg flank
    and bottom. (w, z) with z from the root line (pad top) at w = 0, + outward."""
    R = R or rail()
    c = DESIGN["tongue_clear"]
    w0, wt, wb = R["groove_w_offset"], R["groove_w_top"], R["groove_w_bottom"]
    hn, d = R["groove_neck_h"], R["groove_depth"]
    # the section lies in a radial (w, z) plane but the rail climbs with slope sr = dz/du, so a wall
    # shifted by c inside the section is closer than c in 3D unless it is vertical: shift each wall
    # so its true normal distance is c (planes containing the rail direction (1, 0, sr))
    fr = frame()
    sr = fr["k"] * fr["ci"]
    w_lo, w_hi = w0 + c, w0 + wt - c
    z_bot = -d + c * math.sqrt(1 + sr ** 2)
    fw0, fz0, fw1, fz1 = w0 + wt, -hn, w0 + wb, -d          # flank, from the neck to the bottom
    L = math.hypot(fw1 - fw0, fz1 - fz0)
    qw, qz = (fw1 - fw0) / L, (fz1 - fz0) / L
    slope = (fw1 - fw0) / (fz1 - fz0)                       # dw/dz along the flank (< 0)
    delta = -c * math.sqrt(1 + (sr * qw) ** 2) / abs(qz)    # shift along w toward the groove interior
    z_n = fz0 + (w_hi - fw0 - delta) / slope
    w_b = fw0 + (z_bot - fz0) * slope + delta
    top = DESIGN["root_gap"]                           # coplanar with the fin root face (clean fuse)
    ch = DESIGN["tongue_bed_chamfer"]                  # finding T18: the 90 deg bed edge (flat wall side)
    return [(w_lo, top), (w_hi, top), (w_hi, z_n), (w_b, z_bot), (w_lo + ch, z_bot), (w_lo, z_bot + ch)]


def section_face() -> tuple[Face, dict]:
    """Airfoil section face in the local (u, w) plane at a = a_bottom, chord along u, envelope
    centred on w = 0."""
    fr = frame()
    c = fr["c"]
    up, lo, info = clarky(DESIGN["n_section"], te=P.MIN_TE / c)
    w_c = -0.5 * (info["y_min"] + info["y_max"]) * c        # chord line offset: envelope centred on the rail
    a0 = F["root_r_le"] - 12.0
    pts_u = [Vector(c * x, w_c + c * y, a0) for x, y in up[::-1]]     # TE upper -> LE
    pts_l = [Vector(c * x, w_c + c * y, a0) for x, y in lo[1:]]       # LE -> TE lower
    spl = Edge.make_spline(pts_u + pts_l)
    te = Edge.make_line(pts_l[-1], pts_u[0])
    face = Face(Wire([spl, te]))
    info = dict(info, chord_line_w=round(w_c, 4), a0=a0)
    return face, info


def fin_local() -> tuple[Solid, dict]:
    """The fin + tongue in its local frame (u, w, a)."""
    fr = frame()
    R = rail()
    c, span = fr["c"], F["span"]
    face, info = section_face()
    body = Solid.extrude(face, Vector(0, 0, span + 40.0))
    # planform: root line - 6 .. tip line, rounded tip corners (u -0.5 .. c + 0.5 so LE / TE stay the airfoil's)
    r_le, r_te = DESIGN["tip_corner_r"]
    ar = fr["a_root"]
    pts = [(-0.5, ar(-0.5) - 8.0), (c + 0.5, ar(c + 0.5) - 8.0), (c + 0.5, ar(c + 0.5) + span),
           (-0.5, ar(-0.5) + span)]
    from build123d import Polygon, fillet as _fillet
    pl2 = Polygon(*pts, align=None)
    vs = sorted(pl2.vertices(), key=lambda v: v.Y)[-2:]           # the two tip corners
    v_le = min(vs, key=lambda v: v.X)
    v_te = max(vs, key=lambda v: v.X)
    pl2 = _fillet([v_le], r_le + 0.5)
    pl2 = _fillet([max(sorted(pl2.vertices(), key=lambda v: v.Y)[-2:], key=lambda v: v.X)], r_te + 0.5)
    pf = (Plane.XZ * pl2).faces()[0]                              # local (u, a) -> XZ plane at w = 0
    plan = Solid.extrude(pf.moved(Location((0, -10, 0))), Vector(0, 20, 0))
    body = body & plan
    # root: above the pad top by root_gap, 45 deg flare to the tongue neck
    t_sec = tongue_section(R)
    w_lo, w_hi = t_sec[0][0], t_sec[1][0]
    g = DESIGN["root_gap"]
    m = DESIGN["root_flare"]
    trough = [(w_lo - 40.0, g + 40.0 * m), (w_lo, g), (w_hi, g), (w_hi + 40.0, g + 40.0 * m), (w_hi + 40.0, 300.0),
              (w_lo - 40.0, 300.0)]
    body = body & _rail_prism(trough, -5.0, c + 5.0)
    # tongue from the LE station to the end stop, lead-in at the front
    u_end = R["end_stop_u"]
    tg = _rail_prism(t_sec, 0.0, u_end)
    lead = DESIGN["tongue_lead"]
    # lead-in on the front end's bottom and flank edges only (vertical side chamfers would leave the
    # fin's nose overhanging a notch when printed root down)
    front = [e for e in tg.edges() if abs(e.center().X) < 1e-6 and abs(e.tangent_at(0.5).Z) < 0.9
             and e.center().Z < fr["a_root"](0) + g - 1e-3]
    try:
        tg = tg.chamfer(lead, None, front)
    except Exception:  # noqa: BLE001
        pass
    fin = body.fuse(tg).clean()
    # tip edge rounding
    tip_faces = [f for f in fin.faces() if f.center().Z > fr["a_root"](f.center().X) + span - 1.0]
    tip_edges = [e for f in tip_faces for e in f.edges() if e.length > 2.0]
    info["tip_fillet"] = "none"
    for r in (DESIGN["tip_edge_fillet"], 0.8, 0.5, 0.3):
        try:
            cand = fin.fillet(r, tip_edges)
            if cand.is_valid and len(cand.solids()) == 1:
                fin, info["tip_fillet"] = cand, r
                break
        except Exception:  # noqa: BLE001
            continue
    # cross-pin hole: along w at u = pin_u, depth pin_depth below the root line; teardrop tip toward the
    # print 'up' (normal to the root line in the u-a plane), so it prints without support root down
    rp = R["tongue_pin_hole_d"] / 2
    phi = math.atan(fr["k"] * fr["ci"])
    up_p = Vector(-math.sin(phi), 0, math.cos(phi))
    e_r = Vector(math.cos(phi), 0, math.sin(phi))
    cpt = Vector(R["pin_u"], 10.0, fr["a_root"](R["pin_u"]) - R["pin_depth"])
    hole = Solid.make_cylinder(rp, 20.0, Plane(origin=cpt, x_dir=e_r, z_dir=(0, -1, 0)))
    s2 = rp
    sq = Solid.make_box(s2, s2, 20.0, Plane(origin=cpt, x_dir=(up_p + e_r).normalized(), z_dir=(0, -1, 0)))
    fin = fin - hole.fuse(sq)
    info.update(tongue_section_wz=[(round(a, 4), round(b, 4)) for a, b in t_sec], rail_source=R["source"],
                u_end=u_end, pin=(R["pin_u"], R["pin_depth"], R["tongue_pin_hole_d"]))
    v0 = _gprop_volumes(fin)
    fin = _heal(fin)
    info["gprop_volume_mm3_before_heal"] = {"default": round(v0[0], 2), "adaptive": round(v0[1], 2)}
    return fin, info


HEAL_MAX_FACE_AREA = 100.0       # finding T19: split faces above this area (see _heal)


def _gprop_volumes(s: Solid) -> tuple[float, float]:
    """(default Gauss rule, adaptive 1e-6) B-rep volumes."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(s.wrapped, g)
    g2 = GProp_GProps()
    BRepGProp.VolumeProperties_s(s.wrapped, g2, 1e-6, False)
    return float(g.Mass()), float(g2.Mass())


def _heal(s: Solid) -> Solid:
    """Finding T19: OCC's volume integration of the fin read 8.8 cm3 (default) or 16.2 cm3 (adaptive)
    against a 9.6 cm3 mesh; the big linear-extrusion skin face is to blame (NURBS conversion alone
    did not fix it). Splitting faces above HEAL_MAX_FACE_AREA plus a ShapeFix pass makes both rules
    agree with the mesh; checks() asserts it within 0.5 %."""
    from OCP.ShapeFix import ShapeFix_Shape
    from OCP.ShapeUpgrade import ShapeUpgrade_ShapeDivideArea
    from OCP.TopoDS import TopoDS
    d = ShapeUpgrade_ShapeDivideArea(s.wrapped)
    d.MaxArea = HEAL_MAX_FACE_AREA
    d.Perform()
    fx = ShapeFix_Shape(d.Result())
    fx.Perform()
    out = Solid(TopoDS.Solid(fx.Shape()))
    assert out.is_valid
    return out


def blank_local() -> Solid:
    """Blanking strip (finding T07 / T08): fills the pad's dovetail groove when the fins are packed.
    The tongue section (CLEAR_SLIDE in the groove) from the pad front to the end stop, its top flush
    with the pad top, the same cross-pin hole (teardrop toward the top face, which is up when it
    prints bottom down), a 0.6 mm lead-in at the front."""
    fr = frame()
    R = rail()
    sec = [(w, (0.0 if abs(z - DESIGN["root_gap"]) < 1e-9 else z)) for w, z in tongue_section(R)]
    st = _rail_prism(sec, 0.0, R["end_stop_u"])
    front = [e for e in st.edges() if abs(e.center().X) < 1e-6 and abs(e.tangent_at(0.5).Z) < 0.9
             and e.center().Z < fr["a_root"](0) - 1e-3]
    try:
        st = st.chamfer(DESIGN["tongue_lead"], None, front)
    except Exception:  # noqa: BLE001
        pass
    rp = R["tongue_pin_hole_d"] / 2
    phi = math.atan(fr["k"] * fr["ci"])
    up_p = Vector(-math.sin(phi), 0, math.cos(phi))           # toward the top face (up in print)
    e_r = Vector(math.cos(phi), 0, math.sin(phi))
    cpt = Vector(R["pin_u"], 10.0, fr["a_root"](R["pin_u"]) - R["pin_depth"])
    hole = Solid.make_cylinder(rp, 20.0, Plane(origin=cpt, x_dir=e_r, z_dir=(0, -1, 0)))
    sq = Solid.make_box(rp, rp, 20.0, Plane(origin=cpt, x_dir=(up_p + e_r).normalized(), z_dir=(0, -1, 0)))
    out = st - hole.fuse(sq)
    sols = out.solids()
    out = sols[0] if len(sols) == 1 else out
    assert out.is_valid
    return out


def blank_print(bl: Solid) -> Solid:
    """Bottom down, top face up (support-free: the 45 deg flank narrows going up)."""
    fr = frame()
    phi = math.degrees(math.atan(fr["k"] * fr["ci"]))
    s = bl.rotate(Axis.Y, phi)                                # root line horizontal, strip bottom on the bed
    bb = s.bounding_box()
    return s.moved(Location((-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z)))


_CACHE: dict = {}


def fin_solid(theta: float) -> Solid:
    if "local" not in _CACHE:
        _CACHE["local"] = fin_local()
    loc, _ = _CACHE["local"]
    return loc.moved(plane(theta).location)


def print_pose(fin_loc: Solid) -> Solid:
    """Root down: the tongue's bottom (parallel to the root line, in the local u-w plane tilted by
    the root slope) on the bed."""
    fr = frame()
    phi = math.degrees(math.atan(fr["k"] * fr["ci"]))
    s = fin_loc.rotate(Axis.Y, phi)          # local root line horizontal (rotate about w)
    bb = s.bounding_box()
    return s.moved(Location((-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z)))


def parts() -> list[dict]:
    if "parts" in _CACHE:
        return _CACHE["parts"]
    if "local" not in _CACHE:
        _CACHE["local"] = fin_local()
    loc, info = _CACHE["local"]
    pp = print_pose(loc)
    out = []
    for k, th in enumerate(F["theta_deg"]):
        t = math.radians(th)
        er = [0.0, 40.0 * math.sin(t), 40.0 * math.cos(t)]
        out.append(dict(
            name=f"fin_{k}", shape=fin_solid(th), print=pp, material="fin", role="body", rotor=False, qty=1,
            bought=False, explode=[-12.0, er[1], er[2]], theta_deg=th,
            notes="PETG. Identical part x 2 (theta 93 and 273). Prints ROOT DOWN: the tongue's bottom on the "
                  "bed, span up (leans 15 deg in its own plane), no supports; use brim ears at the two root ends "
                  "(not a full brim: the brim would bond to the tongue's sliding faces); the "
                  "perimeters trace the Clark Y (no stair steps across the chord). Slide aft onto the rail "
                  "to the end stop, then the 2 mm cross pin (ISO 8734 2 x 10) through the pad and the tongue."))
    for p_ in out:
        p_["packed"] = False                                   # removed for packing (assembly PACKED state)
    bl = blank_local()
    bp = blank_print(bl)
    for k, th in enumerate(F["theta_deg"]):
        out.append(dict(
            name=f"fin_blank_{k}", shape=bl.moved(plane(th).location), print=bp, material="shell", role="body",
            rotor=False, qty=1, bought=False, explode=[-12.0, 0.0, 0.0], theta_deg=th, packed=True,
            notes="Blanking strip for the packed egg (fins off): slides into the pad's dovetail groove like the "
                  "fin tongue, top flush with the pad, held by the same 2 mm cross pin. PLA/PETG in the shell "
                  "colour; prints bottom down (top face up), no supports."))
    pin_parts = cross_pins()
    out.append(dict(name="fin_cross_pins", shape=pin_parts, print=None, material="steel", role="body", rotor=False,
                    qty=2, bought=True, explode=[-12.0, 0.0, 0.0],
                    notes="2 x ISO 8734 2 x 10 stainless dowel (one per fin) through the pad and the fin tongue "
                          "(shells.FIN_RAIL), or through the blanking strip when the fins are packed; a 1.75 mm "
                          "filament stub works for display."))
    _CACHE["parts"] = out
    return out


def cross_pins():
    """The two 2 x 10 pins along w through the pads at u = pin_u (shell side position)."""
    from build123d import Compound
    R = rail()
    fr = frame()
    pins = []
    for th in F["theta_deg"]:
        pl = plane(th)
        c = Vector(R["pin_u"], -R["pin_l"] / 2, fr["a_root"](R["pin_u"]) - R["pin_depth"])
        pin = Solid.make_cylinder(R["pin_d"] / 2 - 0.005, R["pin_l"], Plane(origin=c, x_dir=(1, 0, 0), z_dir=(0, 1, 0)))
        ch = [e for e in pin.edges() if str(e.geom_type).endswith("CIRCLE")]
        try:
            pin = pin.chamfer(0.2, None, ch)
        except Exception:  # noqa: BLE001
            pass
        pins.append(pin.moved(pl.location))
    return Compound(pins)


def preview_parts() -> list[dict]:
    return [{k: v for k, v in p.items() if k in ("name", "shape", "material", "role")} for p in parts()]


def preview_with_pads() -> list[dict]:
    """Fins on the shell-side pads (shells._pad_local minus its groove), for a close look at the rail."""
    out = preview_parts()
    try:
        from . import shells as SH
        pad = SH._pad_local()
        g, pin = SH._groove_local()
        pad = pad - g - pin
        for th in F["theta_deg"]:
            out.append(dict(name=f"pad_{th:.0f}", shape=at_azimuth(pad, th), material="shell", role="body"))
    except Exception as ex:  # noqa: BLE001
        print("no shells pads:", ex)
    return out


def preview_print() -> list[dict]:
    p = parts()[0]["print"]
    return [dict(name="fin_print", shape=p, material="fin", role="body")]


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def _mesh(shape, tol=0.03):
    import trimesh
    v, f = shape.tessellate(tol, 0.2)
    m = trimesh.Trimesh(np.array([(p.X, p.Y, p.Z) for p in v]), np.array(f), process=True)
    m.merge_vertices()
    return m


def checks() -> dict:
    import trimesh
    ps = parts()
    loc, info = _CACHE["local"]
    fr = frame()
    res = dict(section=info, parts=[], interfaces=[])
    up, lo, _ = clarky(DESIGN["n_section"], te=P.MIN_TE / F["chord"])
    cmp_ = compare_uiuc(up, lo)
    res["section_vs_uiuc"] = cmp_
    try:
        # couple ratio of the section as built, scaled from the layout's 1.110 (UIUC section) with the
        # layout's lift model: CL = a3D (alpha_set - alpha_L0), a3D from a0 and AR (layout_internal._fin_eval)
        AR = F["span"] / F["chord"]
        a3 = lambda a0: a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))
        nfu, nfc = cmp_["neuralfoil_Re55k"]["uiuc"], cmp_["neuralfoil_Re55k"]["cad"]
        f = (a3(nfc["a0_per_rad"]) * (F["incidence_deg"] - nfc["alpha_L0_deg"])) / (
            a3(nfu["a0_per_rad"]) * (F["incidence_deg"] - nfu["alpha_L0_deg"]))
        res["couple_ratio_as_built"] = dict(layout_ratio=1.110, factor=round(f, 4), ratio=round(1.110 * f, 4),
                                            band=[1.10, 1.155], ok=1.10 <= 1.110 * f <= 1.155,
                                            note="Re 55 k, n_crit 5, NeuralFoil large; CAD section incl. the MIN_TE "
                                                 "trailing edge")
    except Exception as ex:  # noqa: BLE001
        res["couple_ratio_as_built"] = f"not computed: {ex}"
    I = res["interfaces"]

    def add(name, value, ok, **kw):
        I.append(dict(check=name, value=value, ok=bool(ok), **kw))

    for p in ps:
        c = check(p["shape"], p["name"])
        if p["print"] is not None:
            m = _mesh(p["print"])
            n = m.face_normals
            zc = m.triangles_center[:, 2]
            ovh = float(m.area_faces[(n[:, 2] < -math.cos(math.radians(45.0)) - 1e-6) & (zc > 0.05)].sum())
            c.update(print_size_mm=[round(v, 2) for v in p["print"].bounding_box().size], watertight=bool(m.is_watertight),
                     fits_bed=fits_bed(p["print"]), overhang_area_mm2_gt45=round(ovh, 3),
                     bed_contact_mm2=round(float(m.area_faces[(n[:, 2] < -0.999) & (zc < 0.05)].sum()), 1))
        res["parts"].append(c)
    # chord, incidence, LE / TE, span, measured on the world solid of fin_0
    th = F["theta_deg"][0]
    pl = plane(th)
    le = pl.from_local_coords(Vector(0, info["chord_line_w"], fr["a_root"](0)))
    tep = pl.from_local_coords(Vector(fr["c"], info["chord_line_w"], fr["a_root"](fr["c"])))
    t = math.radians(th)
    e_t = Vector(0, math.cos(t), -math.sin(t))
    d = tep - le
    inc = math.degrees(math.atan2(-d.dot(e_t), d.X))
    add("fin incidence (chord vs flow, deg, TE toward -theta)", round(inc, 4), abs(inc - F["incidence_deg"]) < 1e-3)
    add("root chord x_le / x_te at w = 0", [round(fr["x_le"], 4), round(fr["x_te"], 4)],
        abs(fr["x_le"] - F["x_le"]) < 1e-3 and abs(fr["x_te"] - F["x_te"]) < 1e-6)
    # measured on the local solid by probing: chord at mid span along the chord line, span at mid chord
    wc = info["chord_line_w"]
    a_mid = fr["a_root"](fr["c"] / 2) + F["span"] / 2
    inside = lambda u, w, a: loc.is_inside(Vector(u, w, a))            # noqa: E731
    le_ok = inside(0.05, wc, a_mid) and not inside(-0.05, wc, a_mid)
    te_ok = inside(fr["c"] - 0.05, wc, a_mid) and not inside(fr["c"] + 0.05, wc, a_mid)
    add("chord measured at mid span (LE / TE probes +-0.05 mm)", fr["c"], le_ok and te_ok)
    u_m = fr["c"] / 2
    a_tip = fr["a_root"](u_m) + F["span"]
    tip_ok = inside(u_m, wc + 1.0, a_tip - 0.35) and not inside(u_m, wc + 1.0, a_tip + 0.05)
    a_rt = fr["a_root"](u_m) + DESIGN["root_gap"]
    add("span measured at mid chord (root line to tip, probes; tip edge fillet %s)" % info.get("tip_fillet"),
        F["span"], tip_ok, tip_r_te_mm=round(fr["a_root"](fr["c"]) + F["span"], 3))
    add("section max t/c, x(t max), TE thickness (mm)", [info["max_t_c"], info["x_max_t"], round(info["te_c"] * fr["c"], 3)],
        info["te_c"] * fr["c"] >= P.MIN_TE - 1e-6)
    add("Clark Y flat bottom vs chord (deg): generator (no TE ramp) / UIUC", [clarky(120)[2]["flat_bottom_to_chord_deg"], 2.106],
        True, note="why the fin is not printed flat (see module docstring)")
    # rail interface with the shell's pad and groove (shells._pad_local / _groove_local at theta 0)
    try:
        from . import shells as SH
        pad = SH._pad_local()
        g, pin = SH._groove_local()
        pad_cut = pad - g - pin
        fin0 = loc.moved(plane(0.0).location)
        v = interference(fin0, pad_cut)
        add("fin vs shell pad (with its groove) overlap (mm3)", round(v, 4), v < 1e-3)
        dmin = min_distance(fin0, pad_cut)
        add("fin to pad min gap (end stop contact expected 0)", round(dmin, 4), dmin < 0.05)
        # tongue inside the groove: everything of the fin below the pad-top plane must be inside the groove
        below = fin0 & _below_pad_top()
        outside = below.volume - interference(below, g) if below.volume > 0 else 0.0
        add("tongue below the pad top lies inside the groove (volume outside, mm3)", round(outside, 4), outside < 1e-3)
        # clearance tongue <-> groove walls: distance from the tongue to the groove's complement within the pad
        # side / bottom clearance in the groove: slide the fin 1 mm forward along the rail (off the end
        # stop) and measure its distance to the pad
        rv = Vector(fr["ci"], fr["si"], fr["k"] * fr["ci"]).normalized()
        dside = min_distance(fin0.moved(Location(rv * -1.0)), pad_cut)
        add("tongue to groove walls, fin 1 mm off the end stop (mm, CLEAR_SLIDE 0.25)", round(dside, 4),
            abs(dside - P.CLEAR_SLIDE) < 0.02)
        cp = cross_pins()
        v2 = interference(cp, fin_solid(F["theta_deg"][0])) + interference(cp, fin_solid(F["theta_deg"][1]))
        add("cross pins vs fins overlap (mm3, pin D 2.0 in the 2.3 tongue hole)", round(v2, 4), v2 < 1e-3)
        bl0 = blank_local().moved(plane(0.0).location)
        vb = interference(bl0, pad_cut)
        add("blanking strip vs shell pad (with its groove) overlap (mm3)", round(vb, 4), vb < 1e-3)
        outside_b = bl0.volume - interference(bl0, g)
        add("blanking strip lies inside the groove (volume outside, mm3)", round(outside_b, 4), outside_b < 1e-3)
        R_ = rail()
        lig_pad = R_["pin_depth"] - R_["pin_hole_d"] / 2
        lig_tg = (R_["groove_depth"] - DESIGN["tongue_clear"]) - (R_["pin_depth"] + R_["tongue_pin_hole_d"] / 2)
        add("cross-pin ligaments: pad over the D2.0 hole / tongue under the D2.3 hole (mm, T18)",
            [round(lig_pad, 3), round(lig_tg, 3)], min(lig_pad, lig_tg) >= P.MIN_WALL - 0.06,
            note="groove depth %.1f, pin depth %.1f (PROPOSED FIN.groove_depth)" % (R_["groove_depth"], R_["pin_depth"]))
    except Exception as ex:  # noqa: BLE001
        add("rail check against shells._pad_local / _groove_local", f"not run: {ex}", False)
    gd, ga = _gprop_volumes(loc)
    mv = float(_mesh(loc, 0.01).volume)
    add("fin B-rep volume: default / adaptive GProp / 0.01 mm mesh (mm3), agree within 0.5 % (T19)",
        [round(gd, 1), round(ga, 1), round(mv, 1)], abs(gd - ga) <= 0.005 * ga and abs(ga - mv) <= 0.005 * mv,
        before_heal=info.get("gprop_volume_mm3_before_heal"))
    # fin pair: second fin = first turned 180 deg about X
    f0, f1 = fin_solid(F["theta_deg"][0]), fin_solid(F["theta_deg"][1])
    f0r = f0.rotate(Axis.X, -180.0)
    dv = abs(f0r.volume - f1.volume) + (f0r.bounding_box().center() - f1.bounding_box().center()).length
    add("fin_1 = fin_0 turned 180 deg about X (bbox centre deviation, mm)", round(dv, 5), dv < 1e-3)
    add("fin tip radius max (mm) vs layout r_tip_max 109.60", round(_max_radius(f0), 3),
        abs(_max_radius(f0) - (F["root_r_te"] + F["span"])) < 1.5)
    # stress at the root, printed upright: bending across layers (interlayer)
    res["structure"] = dict(note="layout root moment 0.1187 N m (tow case, 2.5 m/s) on the section modulus of "
                                 "the tongue neck + root", **_root_stress(info))
    return res


def _max_radius(s) -> float:
    m = _mesh(s, 0.1)
    v = m.vertices
    return float(np.hypot(v[:, 1], v[:, 2]).max())


def _root_stress(info) -> dict:
    """Bending stress at the fin root from the layout's tow-case root moment, on the Clark Y
    section (second moment about the chord-wise axis), compared with the PETG interlayer strength
    the layout uses (18 MPa): the fin prints upright, so this stress runs across layers."""
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
    M = 0.1187e3                                  # N mm, layout fins.root_moment_tow_Nm
    sig = M * ymax / I
    return dict(area_mm2=round(A, 2), I_mm4=round(I, 2), y_max_mm=round(ymax, 3), sigma_MPa=round(sig, 3),
                interlayer_MPa=18.0, utilisation=round(sig / 18.0, 4))


def _below_pad_top():
    return _pad_top_half(False)


def _pad_top_half(above: bool) -> Solid:
    """Half space above / below the pad-top plane of the theta-0 pad (shells.fin_frame: z = root_z(x),
    independent of y), clipped to the pad region."""
    fr = frame()
    phi = math.atan(fr["k"])
    n = Vector(-math.sin(phi), 0, math.cos(phi)) * (1 if above else -1)
    p0 = Vector(fr["x_le"], 0, F["root_r_le"])
    pl = Plane(origin=p0, z_dir=n)
    size = 200.0
    return Solid.make_box(size, size, 60.0, Plane(origin=p0 - pl.x_dir * (size / 2) - pl.y_dir * (size / 2),
                                                 x_dir=pl.x_dir, z_dir=n))


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
        export_step(p["shape"], str(sp))
        files.append(str(sp.relative_to(ROOT)))
        if p["print"] is not None:
            st = OUT / "stl" / f"{p['name']}.stl"
            export_stl(p["print"], str(st), tolerance=0.02, angular_tolerance=0.15)
            files.append(str(st.relative_to(ROOT)))
            for c in res["parts"]:
                if c["name"] == p["name"]:
                    c["stl_watertight"] = bool(trimesh.load(st).is_watertight)
    res["files"] = files
    res["design"] = DESIGN
    import sys
    if "--no-render" not in sys.argv:
        from .drivetrain import render_previews
        res["renders"] = render_previews("hydroegg.cad.fins", [
            ("preview_with_pads", "fins_on_pads", "view=top"),
            ("preview_print", "fins_print", "view=-30,15")])
    res["elapsed_s"] = round(time.time() - t0, 1)
    (OUT / "checks").mkdir(parents=True, exist_ok=True)
    (OUT / "checks" / "fins.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res["interfaces"], indent=1, default=str))
    for c in res["parts"]:
        print(c)
    print("elapsed", res["elapsed_s"])


if __name__ == "__main__":
    main()
