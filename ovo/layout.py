"""Rotor-side layout of Ovo: fold kinematics, stowage pockets, hub and hinge.

Run ``python -m ovo.layout`` to regenerate ``data/layout_rotor.json`` and
``docs/figures/layout_fold_side.png`` / ``layout_fold_end.png``.

Everything here is plain numpy (matplotlib only inside the figure function), so the
CAD can import the transforms, the blade surface and the pocket tables directly.

Frames and conventions
----------------------
World: X = body axis, flow travels +X, nose at x = 0, rotor at the tail, Z up.
Units: mm and deg for geometry, SI for physics.

Blade azimuth ``psi`` (deg) is measured about +X, right handed, from +Z:
``e_r(psi) = (0, -sin psi, cos psi)`` and ``e_t(psi) = (0, -cos psi, -sin psi)``,
so ``e_r x e_t = e_x`` and blade 0 points straight up. The rotor turns in +psi
(right handed about +X, counter-clockwise seen from aft); ``e_t`` is the direction
of motion, so the leading edge (LE) is on the +e_t side.

Hinge: a pin along ``e_t`` through ``H = (HINGE_X, 0, 0) + HINGE_R e_r``.

Blade frame (origin on the hinge axis at H, right handed, axes (s, c, n)):

* ``e_s``  span axis. Deployed (fold 0): ``e_r``, radially out.
* ``e_c``  chord axis, along the pin. Always ``e_t`` (toward the LE).
* ``e_n``  thickness axis, ``e_s x e_c``. Deployed: ``+e_x`` (aft). The suction
           (convex) side faces +n, the pressure side faces -n (upstream).

Fold angle ``fold`` (deg) rotates the blade about ``e_t`` by +fold:
``e_s = cos f e_r - sin f e_x`` and ``e_n = sin f e_r + cos f e_x``. fold = 0 is
deployed (radial), fold ~ 90 is stowed with the tip toward the nose; the aft axis
``e_n`` then points outward, so an aft pre-cone becomes an outward offset and the
pressure face lies against the hull.

Blade reference line (pitch axis at PITCH_AXIS chord) in the blade frame:
``P(r) = (s(r), lead(r), aft(r))``. ``aft`` is the pre-cone offset (deployed: aft of
the hinge plane), ``lead`` an in-plane offset toward the LE (tangential when
deployed and when stowed: the fold axis is e_t, so it is invariant). The deployed
perpendicular radius of P is ``r = sqrt((HINGE_R + s)^2 + lead^2)``.

Sections follow the bem.py convention: twist theta is the angle between the chord
and e_c measured in the section plane normal to the local pitch axis (plane
spanned by e_c and ``n_loc = unit(e_n - aft'(s) e_s)``), positive with the LE
upstream. LE direction ``L = cos th e_c - sin th n_loc``, suction normal
``U = sin th e_c + cos th n_loc``.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from . import params as P
from . import profile

ROOT = Path(__file__).resolve().parents[1]
DATA_JSON = ROOT / "data" / "layout_rotor.json"
BLADE_JSON = ROOT / "data" / "blade_design.json"
FIG_DIR = ROOT / "docs" / "figures"
WORK_DIR = ROOT / "out" / "work" / "layout_rotor"

PLA_DENSITY = 1.24e-3        # g/mm^3, solid PLA (estimate for masses)


# ==========================================================================
# Hull helpers
# ==========================================================================

def hull_r(x) -> np.ndarray:
    """Outer hull radius r(x) (mm), vectorised; 0 outside [0, HULL_L]."""
    return np.asarray(profile.radius(np.asarray(x, dtype=float)), dtype=float)


def hull_slope(x, h: float = 1e-3) -> np.ndarray:
    """dr/dx of the hull, vectorised central difference."""
    x = np.asarray(x, dtype=float)
    return (hull_r(x + h) - hull_r(x - h)) / (2 * h)


def height_above_hull(pts: np.ndarray) -> np.ndarray:
    """Signed distance of world points from the hull surface of revolution,
    measured along the local surface normal (first-order: radial gap times
    cos(slope angle)). Positive outside the hull."""
    x = pts[..., 0]
    rho = np.hypot(pts[..., 1], pts[..., 2])
    return (rho - hull_r(x)) / np.sqrt(1.0 + hull_slope(x) ** 2)


def cyl(pts: np.ndarray, psi_deg: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """World points -> (x, rho, azimuth_deg relative to psi, wrapped to +-180)."""
    x = pts[..., 0]
    rho = np.hypot(pts[..., 1], pts[..., 2])
    # azimuth about +X from +Z, right handed: az = atan2(-y, z)
    az = np.degrees(np.arctan2(-pts[..., 1], pts[..., 2])) - psi_deg
    az = (az + 180.0) % 360.0 - 180.0
    return x, rho, az


# ==========================================================================
# Blade definition (provisional, or data/blade_design.json when it exists)
# ==========================================================================

@dataclass(frozen=True)
class BladeDef:
    """Spanwise blade distributions vs deployed perpendicular radius r (mm)."""

    r_root: float = 30.0            # start of the aerodynamic span
    r_tip: float = P.ROTOR_R_TIP_TARGET
    chord_root: float = 26.0
    chord_tip: float = 14.0
    tc_root: float = 0.18
    tc_tip: float = 0.12
    twist_root: float = 22.0        # deg, chord to rotor plane, LE upstream
    twist_tip: float = 3.0
    pitch_axis: float = 0.30        # chord fraction of the stacking axis
    te_mm: float = P.MIN_TE
    camber_m: float = 0.04          # NACA 44xx, as bem.py
    camber_p: float = 0.40
    shank_chord: float = 20.0       # nominal, used only for reporting
    r_tab: tuple | None = None      # optional tabulated distributions (BEM output)
    chord_tab: tuple | None = None
    twist_tab: tuple | None = None
    tc_tab: tuple | None = None
    source: str = "provisional (task brief: linear chord/twist/t_c)"

    def _f(self, r):
        return (np.asarray(r, float) - self.r_root) / (self.r_tip - self.r_root)

    def chord(self, r):
        if self.r_tab is not None and self.chord_tab is not None:
            return np.interp(r, self.r_tab, self.chord_tab)
        return self.chord_root + (self.chord_tip - self.chord_root) * self._f(r)

    def twist(self, r):
        if self.r_tab is not None and self.twist_tab is not None:
            return np.interp(r, self.r_tab, self.twist_tab)
        return self.twist_root + (self.twist_tip - self.twist_root) * self._f(r)

    def t_c(self, r):
        if self.r_tab is not None and self.tc_tab is not None:
            return np.interp(r, self.r_tab, self.tc_tab)
        return self.tc_root + (self.tc_tip - self.tc_root) * np.clip(self._f(r), 0, 1)


def _find_dist(obj):
    """Depth-first search of a JSON tree for a dict holding r/chord/twist arrays."""
    if isinstance(obj, dict):
        keys = set(obj)
        rk = next((k for k in ("r_mm", "r") if k in keys), None)
        ck = next((k for k in ("chord_mm", "chord") if k in keys), None)
        tk = next((k for k in ("twist_deg", "twist") if k in keys), None)
        if rk and ck and tk and isinstance(obj[rk], list) and len(obj[rk]) > 2:
            return obj, rk, ck, tk
        for v in obj.values():
            hit = _find_dist(v)
            if hit:
                return hit
    elif isinstance(obj, list):
        for v in obj:
            hit = _find_dist(v)
            if hit:
                return hit
    return None


def load_blade(path: Path = BLADE_JSON) -> BladeDef:
    """BladeDef from data/blade_design.json if present and readable, else the
    provisional blade. Accepts either tabulated arrays (r_mm, chord_mm,
    twist_deg, optional t_c) anywhere in the tree, preferring a top-level
    "final"/"blade"/"design" entry, or linear root/tip scalars."""
    if not path.exists():
        return BladeDef()
    try:
        d = json.loads(path.read_text())
    except Exception:  # noqa: BLE001  (half-written file while the BEM study runs)
        return BladeDef(source=f"provisional ({path.name} unreadable)")
    # schema ovo.blade_design/1 (ovo.bem): stations + section + rotor blocks
    stn = d.get("stations") if isinstance(d, dict) else None
    if isinstance(stn, dict) and all(k in stn for k in ("radius_mm", "chord_mm", "twist_deg")):
        sec = d.get("section", {}) or {}
        rot = d.get("rotor", {}) or {}
        r = np.asarray(stn["radius_mm"], float)
        kw = dict(r_tab=tuple(r), chord_tab=tuple(np.asarray(stn["chord_mm"], float)),
                  twist_tab=tuple(np.asarray(stn["twist_deg"], float)),
                  r_root=float(rot.get("r_aero_root_mm", r.min())),
                  r_tip=float(rot.get("r_tip_mm", r.max())),
                  chord_root=float(stn["chord_mm"][0]), chord_tip=float(stn["chord_mm"][-1]),
                  twist_root=float(stn["twist_deg"][0]), twist_tip=float(stn["twist_deg"][-1]),
                  source=f"{path.name} ({d.get('schema', '?')}): stations radius/chord/twist/t_c")
        if "thickness_ratio" in stn:
            kw["tc_tab"] = tuple(np.asarray(stn["thickness_ratio"], float))
            kw["tc_root"], kw["tc_tip"] = float(stn["thickness_ratio"][0]), float(stn["thickness_ratio"][-1])
        for k_json, k_def in (("pitch_axis_chord_fraction", "pitch_axis"), ("te_thickness_mm", "te_mm"),
                              ("camber_m", "camber_m"), ("camber_p", "camber_p")):
            if k_json in sec:
                kw[k_def] = float(sec[k_json])
        sh = d.get("shank", {}) or {}
        if "lug_chord_mm" in sh:
            kw["shank_chord"] = float(sh["lug_chord_mm"])
        return replace(BladeDef(), **kw)
    for pref in ("final", "blade", "design", "geometry"):
        if isinstance(d, dict) and pref in d:
            hit = _find_dist(d[pref])
            if hit:
                break
    else:
        hit = _find_dist(d)
    base = BladeDef()
    if hit:
        node, rk, ck, tk = hit
        r = np.asarray(node[rk], float)
        kw = dict(r_tab=tuple(r), chord_tab=tuple(np.asarray(node[ck], float)),
                  twist_tab=tuple(np.asarray(node[tk], float)),
                  r_root=float(node.get("r_aero_root_mm", r.min())),
                  r_tip=float(node.get("r_tip_mm", r.max())),
                  source=f"{path.name}: tabulated '{rk}/{ck}/{tk}'")
        tck = next((k for k in ("t_c", "tc", "thickness_ratio") if k in node), None)
        if tck:
            kw["tc_tab"] = tuple(np.asarray(node[tck], float))
        for k_json, k_def in (("pitch_axis", "pitch_axis"), ("te_mm", "te_mm")):
            if k_json in node:
                kw[k_def] = float(node[k_json])
        return replace(base, **kw)
    flat = d.get("final", d) if isinstance(d, dict) else {}
    keys = ("chord_root_mm", "chord_tip_mm", "twist_root_deg", "twist_tip_deg")
    if isinstance(flat, dict) and all(k in flat for k in keys):
        return replace(base, chord_root=float(flat["chord_root_mm"]),
                       chord_tip=float(flat["chord_tip_mm"]),
                       twist_root=float(flat["twist_root_deg"]),
                       twist_tip=float(flat["twist_tip_deg"]),
                       source=f"{path.name}: linear root/tip scalars")
    return replace(base, source=f"provisional ({path.name} has no blade distribution)")


# ==========================================================================
# Layout parameters (the study's choices; proposals for params.py)
# ==========================================================================

@dataclass(frozen=True)
class Layout:
    """All rotor-side layout choices. Defaults are the values this study proposes;
    ``Layout.from_params()`` builds the one currently in params.py."""

    static_tail_x: float = P.STATIC_TAIL_X
    hinge_x: float = P.HINGE_X
    hinge_r: float = P.HINGE_R
    fold_deg: float = 90.0          # stowed fold angle (0 = deployed)
    hub_gap: float = P.HUB_GAP
    n_blades: int = P.N_BLADES
    psi0_deg: float = 0.0
    # hinge hardware
    pin_d: float = P.PIN_D
    pin_l: float = P.PIN_L
    knuckle_r: float = 4.0          # blade knuckle outer radius about the pin axis
    knuckle_l: float = 14.0         # knuckle length along the pin
    ear_t: float = 3.5              # hub ear thickness along the pin
    ear_clear: float = 0.3          # knuckle <-> ear axial clearance (each side)
    ear_cover: float = 2.0          # target material around the pin
    spinner_wall: float = 1.2       # spinner shell over the hub (MIN_WALL)
    spinner_gap: float = 0.2        # radial gap hub core <-> spinner shell
    run_clear: float = 0.3          # running clearance blade <-> hub (Prusa: >= 0.3)
    magnet_recess: float = 0.2      # magnet face below its contact face (each side)
    shaft_end_dx: float = 12.0      # shaft end (hub aft face) aft of the hinge plane
    # deployed stop: heel tab under the knuckle
    heel_reach: float = 7.0         # pin axis -> heel tip
    heel_t: float = 4.0             # heel thickness (in the s-n plane)
    heel_angle_deg: float = 0.0     # deployed heel direction from straight down, + aft
    heel_w: float = 10.0            # heel width along the pin (narrower than the knuckle)
    # shank
    shank_cone_deg: float = 25.0    # local cone of the shank spine at the knuckle
    shank_blend_start: float = 4.5  # s where the knuckle section starts to morph (> ear zone)
    # stowage
    proud: float = 0.0              # allowed protrusion of the stowed blade over the hull
    pocket_clear: float = 0.3       # pocket floor clearance under the aero span
    root_clear: float = 0.8         # min clearance stowed root <-> static tail
    side_clear: float = 0.8         # pocket side-wall clearance
    floor_overhang_deg: float = P.MAX_OVERHANG   # pocket floor slope limit for printing (see printable_floor)

    @classmethod
    def from_params(cls, **kw) -> "Layout":
        return cls(static_tail_x=P.STATIC_TAIL_X, hinge_x=P.HINGE_X, hinge_r=P.HINGE_R,
                   pin_l=P.PIN_L, **kw)

    def psi(self, i: int) -> float:
        return self.psi0_deg + 360.0 * i / self.n_blades


# ==========================================================================
# Frames and transforms
# ==========================================================================

def rotor_axes(psi_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(e_x, e_r, e_t) for azimuth psi (deg from +Z, right handed about +X)."""
    p = math.radians(psi_deg)
    return (np.array([1.0, 0.0, 0.0]),
            np.array([0.0, -math.sin(p), math.cos(p)]),
            np.array([0.0, -math.cos(p), -math.sin(p)]))


def hinge_point(lay: Layout, psi_deg: float) -> np.ndarray:
    _, e_r, _ = rotor_axes(psi_deg)
    return np.array([lay.hinge_x, 0.0, 0.0]) + lay.hinge_r * e_r


def blade_to_world(fold_deg: float, psi_deg: float, lay: Layout) -> np.ndarray:
    """Rigid 4x4 transform, blade frame (s, c, n) -> world, for a blade at azimuth
    psi folded by fold_deg (0 deployed, ~90 stowed toward the nose)."""
    e_x, e_r, e_t = rotor_axes(psi_deg)
    f = math.radians(fold_deg)
    e_s = math.cos(f) * e_r - math.sin(f) * e_x
    e_n = math.sin(f) * e_r + math.cos(f) * e_x
    T = np.eye(4)
    T[:3, 0], T[:3, 1], T[:3, 2] = e_s, e_t, e_n
    T[:3, 3] = hinge_point(lay, psi_deg)
    return T


def blade_transform(lay: Layout, i: int, fold_deg: float) -> np.ndarray:
    """Transform for blade index i (psi = psi0 + i * 360 / N)."""
    return blade_to_world(fold_deg, lay.psi(i), lay)


def apply(T: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Apply a 4x4 rigid transform to (..., 3) points."""
    return pts @ T[:3, :3].T + T[:3, 3]


# ==========================================================================
# Sections
# ==========================================================================

N_PER_SIDE = 41


def naca44(t_c: float, te_c: float, m: float = 0.04, p: float = 0.40,
           n: int = N_PER_SIDE) -> np.ndarray:
    """Unit-chord cambered NACA 4-digit section with a blunt TE, same construction
    as bem.naca44_section: (2n-1, 2) points, Selig order (upper TE -> LE -> lower
    TE), x from the LE, y toward the suction side."""
    b = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(b))
    yt = 5 * t_c * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                    + 0.2843 * x**3 - 0.1036 * x**4)
    yc = np.where(x <= p, m / p**2 * (2 * p * x - x**2),
                  m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))
    dyc = np.where(x <= p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th) + 0.5 * te_c * x
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th) - 0.5 * te_c * x
    return np.vstack([np.column_stack([xu[::-1], yu[::-1]]),
                      np.column_stack([xl[1:], yl[1:]])])


def rounded_rect(t_over_c: float, expo: float = 4.0, n: int = N_PER_SIDE) -> np.ndarray:
    """Unit-chord super-ellipse section (|2x-1|^e + |2y/t|^e = 1), same point count
    and order as naca44, for the shank/knuckle blend."""
    b = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(b))
    y = 0.5 * t_over_c * np.clip(1 - np.abs(2 * x - 1) ** expo, 0, None) ** (1 / expo)
    return np.vstack([np.column_stack([x[::-1], y[::-1]]),
                      np.column_stack([x[1:], -y[1:]])])


def section_frame(slope_aft: float, twist_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """(L, U) unit vectors in blade coordinates (s, c, n): LE direction and suction
    normal for a section whose pitch axis has d(aft)/ds = slope_aft."""
    slope_aft = np.asarray(slope_aft, float)
    th = np.radians(np.asarray(twist_deg, float))
    k = 1.0 / np.sqrt(1.0 + slope_aft**2)
    n_loc = np.stack([-slope_aft * k, np.zeros_like(k), k], axis=-1)
    e_c = np.stack([np.zeros_like(k), np.ones_like(k), np.zeros_like(k)], axis=-1)
    L = np.cos(th)[..., None] * e_c - np.sin(th)[..., None] * n_loc
    U = np.sin(th)[..., None] * e_c + np.cos(th)[..., None] * n_loc
    return L, U


def aero_sections(blade: BladeDef, r: np.ndarray, s: np.ndarray, lead: np.ndarray,
                  aft: np.ndarray, slope: np.ndarray) -> np.ndarray:
    """(Ns, Np, 3) blade-frame points of the aerodynamic sections at radii r whose
    pitch points are (s, lead, aft) with d(aft)/ds = slope."""
    r = np.atleast_1d(np.asarray(r, float))
    c = blade.chord(r)
    tc = blade.t_c(r)
    tw = blade.twist(r)
    L, U = section_frame(slope, tw)
    out = []
    for k in range(r.size):
        xy = naca44(float(tc[k]), blade.te_mm / float(c[k]), blade.camber_m, blade.camber_p)
        Pk = np.array([s[k], lead[k], aft[k]])
        out.append(Pk + np.outer((blade.pitch_axis - xy[:, 0]) * c[k], L[k])
                   + np.outer(xy[:, 1] * c[k], U[k]))
    return np.array(out)


class AeroStations:
    """Pre-computed unit sections at fixed radii; ``points`` places them anywhere."""

    def __init__(self, blade: BladeDef, r: np.ndarray, n_per_side: int = N_PER_SIDE):
        self.blade = blade
        self.r = np.asarray(r, float)
        self.c = np.asarray(blade.chord(self.r), float)
        self.tw = np.asarray(blade.twist(self.r), float)
        self.tc = np.asarray(blade.t_c(self.r), float)
        self.xy = np.stack([naca44(float(t), blade.te_mm / float(c), blade.camber_m,
                                   blade.camber_p, n_per_side)
                            for t, c in zip(self.tc, self.c)])

    def points(self, s, lead, aft, slope) -> np.ndarray:
        """(Ns, Np, 3) blade-frame section points for pitch points (s, lead, aft)."""
        L, U = section_frame(slope, self.tw)
        Pp = np.stack([s, lead, aft], axis=-1)[:, None, :]
        u = ((self.blade.pitch_axis - self.xy[..., 0]) * self.c[:, None])[..., None]
        v = (self.xy[..., 1] * self.c[:, None])[..., None]
        return Pp + u * L[:, None, :] + v * U[:, None, :]


def s_of_r(r, lead, hinge_r: float):
    """Span coordinate s of the pitch point whose deployed perpendicular radius is r."""
    return np.sqrt(np.asarray(r, float) ** 2 - np.asarray(lead, float) ** 2) - hinge_r


def stowed_heights(st: AeroStations, lay: Layout, s, lead, aft, slope,
                   fold_deg: float | None = None, psi_deg: float = 0.0) -> np.ndarray:
    """(Ns, Np) normal heights of the stowed section points above the hull."""
    fold = lay.fold_deg if fold_deg is None else fold_deg
    W = apply(blade_to_world(fold, psi_deg, lay), st.points(s, lead, aft, slope))
    return height_above_hull(W)


def place_flush(st: AeroStations, lay: Layout, lead, slope, fold_deg=None,
                proud: float | None = None, lo: float = -10.0, hi: float = 60.0,
                n_iter: int = 48) -> tuple[np.ndarray, np.ndarray]:
    """Aft offset per station such that the highest point of the stowed section sits
    ``proud`` mm above the hull (0 = flush). Vectorised bisection; returns (aft, s)."""
    proud = lay.proud if proud is None else proud
    lead = np.broadcast_to(np.asarray(lead, float), st.r.shape)
    slope = np.broadcast_to(np.asarray(slope, float), st.r.shape)
    s = s_of_r(st.r, lead, lay.hinge_r)
    a_lo = np.full(st.r.shape, lo)
    a_hi = np.full(st.r.shape, hi)
    for _ in range(n_iter):
        mid = 0.5 * (a_lo + a_hi)
        h = stowed_heights(st, lay, s, lead, mid, slope, fold_deg).max(axis=1)
        high = h > proud
        a_hi = np.where(high, mid, a_hi)
        a_lo = np.where(high, a_lo, mid)
    return 0.5 * (a_lo + a_hi), s


def flush_depth(st: AeroStations, lay: Layout, lead, slope, fold_deg=None, proud=None):
    """(depth, aft, s): depth = deepest point of the flush-placed stowed section
    below the hull (normal, mm)."""
    aft, s = place_flush(st, lay, lead, slope, fold_deg, proud)
    h = stowed_heights(st, lay, s, np.broadcast_to(lead, st.r.shape), aft,
                       np.broadcast_to(slope, st.r.shape), fold_deg)
    return -h.min(axis=1), aft, s


def optimise_lead(st: AeroStations, lay: Layout, slope, fold_deg=None, proud=None,
                  lo: float = -8.0, hi: float = 24.0, n_iter: int = 32) -> np.ndarray:
    """Per-station in-plane lead offset that minimises the flush pocket depth
    (vectorised golden-section search)."""
    g = (math.sqrt(5) - 1) / 2
    a = np.full(st.r.shape, lo)
    b = np.full(st.r.shape, hi)
    x1 = b - g * (b - a)
    x2 = a + g * (b - a)
    f1 = flush_depth(st, lay, x1, slope, fold_deg, proud)[0]
    f2 = flush_depth(st, lay, x2, slope, fold_deg, proud)[0]
    for _ in range(n_iter):
        left = f1 < f2
        b = np.where(left, x2, b)
        a = np.where(left, a, x1)
        x2n = np.where(left, x1, a + g * (b - a))
        x1n = np.where(left, b - g * (b - a), x2)
        x1, x2 = x1n, x2n
        fn = flush_depth(st, lay, np.where(left, x1, x2), slope, fold_deg, proud)[0]
        f1, f2 = np.where(left, fn, f2), np.where(left, f1, fn)
    return 0.5 * (a + b)


# ==========================================================================
# Blade reference line: pre-cone (aft) and in-plane lead, derived from the hull
# ==========================================================================

@dataclass
class BladeAxis:
    """Pitch-axis curve in the blade frame, tabulated vs deployed radius r.

    ``aft(r)``: pre-cone offset (deployed: aft of the hinge plane), ``lead(r)``:
    in-plane offset toward the LE, ``s(r)``: span coordinate, ``slope(r)`` =
    d(aft)/ds. Polynomial coefficients (in t = (r - r_root)/(r_tip - r_root))
    give a smooth closed form for the CAD."""

    r: np.ndarray
    s: np.ndarray
    lead: np.ndarray
    aft: np.ndarray
    slope: np.ndarray
    aft_coef: np.ndarray
    lead_coef: np.ndarray
    r_root: float
    r_tip: float
    lead_station_opt: np.ndarray | None = None

    def t(self, r):
        return (np.asarray(r, float) - self.r_root) / (self.r_tip - self.r_root)

    def lead_at(self, r):
        return np.polyval(self.lead_coef, self.t(r))

    def aft_at(self, r):
        return np.polyval(self.aft_coef, self.t(r))

    def slope_at(self, r, lead=None):
        """d(aft)/ds at r (ds = dr * r / sqrt(r^2 - lead^2) ... from s(r))."""
        r = np.asarray(r, float)
        dt = 1.0 / (self.r_tip - self.r_root)
        daft_dr = np.polyval(np.polyder(self.aft_coef), self.t(r)) * dt
        ld = self.lead_at(r) if lead is None else lead
        dld_dr = np.polyval(np.polyder(self.lead_coef), self.t(r)) * dt
        ds_dr = (r - ld * dld_dr) / np.sqrt(r**2 - ld**2)
        return daft_dr / ds_dr

    def cone_deg(self, r=None):
        """Local deployed cone of the pitch axis, deg (aft positive)."""
        r = self.r if r is None else r
        return np.degrees(np.arctan(self.slope_at(r)))


def solve_axis(blade: BladeDef, lay: Layout, n: int = 41, lead: str | float = "opt",
               lead_deg: int = 2, aft_deg: int = 5, iters: int = 3,
               proud: float | None = None) -> BladeAxis:
    """Derive the pre-cone (aft) and lead curves so that the stowed blade sits flush
    on the hull (highest point ``proud`` above the envelope) with the shallowest
    possible pocket.

    lead = "opt": per-station optimum, then a polynomial fit (degree lead_deg) so the
    planform stays smooth; a number: constant lead (0 = radial stacking).
    The section planes depend on d(aft)/ds, so the solve is iterated."""
    r = np.linspace(blade.r_root, blade.r_tip, n)
    t = (r - blade.r_root) / (blade.r_tip - blade.r_root)
    st = AeroStations(blade, r)
    slope = np.zeros_like(r)
    lead_opt = None
    for _ in range(iters):
        if lead == "opt":
            lead_opt = optimise_lead(st, lay, slope, proud=proud)
            lead_coef = np.polyfit(t, lead_opt, lead_deg)
        else:
            lead_coef = np.array([float(lead)])
        ld = np.polyval(lead_coef, t)
        aft, s = place_flush(st, lay, ld, slope, proud=proud)
        aft_coef = np.polyfit(t, aft, aft_deg)
        ax = BladeAxis(r, s, ld, np.polyval(aft_coef, t), slope, aft_coef, lead_coef,
                       blade.r_root, blade.r_tip, lead_opt)
        slope = ax.slope_at(r, ld)
    ax.slope = slope
    return ax


def axis_with_flush_residual(blade: BladeDef, lay: Layout, ax: BladeAxis) -> dict:
    """Per-station stowed metrics of a solved axis (heights from the fitted curve)."""
    st = AeroStations(blade, ax.r)
    h = stowed_heights(st, lay, ax.s, ax.lead, ax.aft, ax.slope)
    W = apply(blade_to_world(lay.fold_deg, 0.0, lay), st.points(ax.s, ax.lead, ax.aft, ax.slope))
    return {"h_max": h.max(axis=1), "depth": -h.min(axis=1), "world": W}


# ==========================================================================
# Root: shank, knuckle and heel (blade frame)
# ==========================================================================

def _hermite(s, s0, y0, m0, s1, y1, m1):
    """Cubic Hermite on [s0, s1]."""
    h = s1 - s0
    t = (np.asarray(s, float) - s0) / h
    h00, h10 = 2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t
    h01, h11 = -2 * t**3 + 3 * t**2, t**3 - t**2
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def _hermite_d(s, s0, y0, m0, s1, y1, m1):
    h = s1 - s0
    t = (np.asarray(s, float) - s0) / h
    return ((6 * t**2 - 6 * t) * y0 + (3 * t**2 - 4 * t + 1) * h * m0
            + (-6 * t**2 + 6 * t) * y1 + (3 * t**2 - 2 * t) * h * m1) / h


def shank_sections(blade: BladeDef, lay: Layout, ax: BladeAxis, n: int = 16,
                   n_per_side: int = N_PER_SIDE) -> tuple[np.ndarray, np.ndarray]:
    """(Ns, Np, 3) blade-frame sections of the root shank, from the knuckle (s = 0)
    to the aero root, and their span coordinates.

    The shank is a loft: a rounded rectangle (knuckle length x knuckle diameter,
    no twist) centred on the pin, morphing (smoothstep) into the root airfoil from
    s = shank_blend_start. Its chord-midpoint path leaves the knuckle with a local
    cone ``shank_cone_deg`` and meets the aero root tangentially."""
    r0 = blade.r_root
    s0, ld0, aft0, sl0 = float(ax.s[0]), float(ax.lead[0]), float(ax.aft[0]), float(ax.slope[0])
    c0, tw0 = float(blade.chord(r0)), float(blade.twist(r0))
    L0, _ = section_frame(sl0, tw0)
    M0 = np.array([s0, ld0, aft0]) + (blade.pitch_axis - 0.5) * c0 * L0   # root chord midpoint
    m_k = math.tan(math.radians(lay.shank_cone_deg))
    s = np.linspace(0.0, s0, n)
    n_mid = _hermite(s, 0.0, 0.0, m_k, s0, M0[2], sl0)
    sl = _hermite_d(s, 0.0, 0.0, m_k, s0, M0[2], sl0)
    sb = min(lay.shank_blend_start, 0.8 * s0)
    # stay within the knuckle width until clear of the pin ears, then shift to the root
    c_mid = np.where(s <= sb, 0.0, _hermite(np.maximum(s, sb), sb, 0.0, 0.0, s0, M0[1], 0.0))
    tt = np.clip((s - sb) / (s0 - sb), 0, 1)
    u = tt * tt * (3 - 2 * tt)
    tw = u * tw0
    foil = naca44(float(blade.t_c(r0)), blade.te_mm / c0, blade.camber_m, blade.camber_p, n_per_side)
    rect = rounded_rect(2 * lay.knuckle_r / lay.knuckle_l, 4.0, n_per_side)
    X = (1 - u)[:, None] * (0.5 - rect[:, 0])[None, :] * lay.knuckle_l \
        + u[:, None] * (0.5 - foil[:, 0])[None, :] * c0
    Y = (1 - u)[:, None] * rect[:, 1][None, :] * lay.knuckle_l \
        + u[:, None] * foil[:, 1][None, :] * c0
    L, U = section_frame(sl, tw)
    M = np.stack([s, c_mid, n_mid], axis=-1)
    pts = M[:, None, :] + X[..., None] * L[:, None, :] + Y[..., None] * U[:, None, :]
    return pts, s


def knuckle_points(lay: Layout, n_a: int = 48, n_t: int = 9) -> np.ndarray:
    """(n_t, n_a, 3) surface points of the knuckle cylinder (radius knuckle_r about
    the pin axis e_c, length knuckle_l), blade frame."""
    a = np.linspace(0, 2 * math.pi, n_a)
    tau = np.linspace(-lay.knuckle_l / 2, lay.knuckle_l / 2, n_t)
    S = lay.knuckle_r * np.cos(a)[None, :] * np.ones((n_t, 1))
    N = lay.knuckle_r * np.sin(a)[None, :] * np.ones((n_t, 1))
    C = tau[:, None] * np.ones((1, n_a))
    return np.stack([S, C, N], axis=-1)


def heel_dir(lay: Layout) -> tuple[np.ndarray, np.ndarray]:
    """(d, q): heel direction and its in-plane normal in blade (s, n) coordinates.
    Deployed, d points heel_angle_deg aft of straight down (-e_r)."""
    b = math.radians(lay.heel_angle_deg)
    d = np.array([-math.cos(b), math.sin(b)])
    q = np.array([math.sin(b), math.cos(b)])      # +q = toward aft (+n) side
    return d, q


def heel_points(lay: Layout, n_l: int = 12, n_w: int = 7, n_t: int = 9) -> np.ndarray:
    """Surface point cloud (M, 3) of the heel tab (stop lug under the knuckle):
    box from the pin axis to heel_reach along heel_dir, heel_t thick, heel_w long."""
    d, q = heel_dir(lay)
    lam = np.linspace(0.0, lay.heel_reach, n_l)
    w = np.linspace(-lay.heel_t / 2, lay.heel_t / 2, n_w)
    tau = np.linspace(-lay.heel_w / 2, lay.heel_w / 2, n_t)
    LL, WW, TT = np.meshgrid(lam, w, tau, indexing="ij")
    on = ((np.abs(WW) >= lay.heel_t / 2 - 1e-9) | (LL >= lay.heel_reach - 1e-9)
          | (np.abs(TT) >= lay.heel_w / 2 - 1e-9))
    S = LL * d[0] + WW * q[0]
    N = LL * d[1] + WW * q[1]
    pts = np.stack([S, TT, N], axis=-1)[on]
    return pts


@dataclass
class BladeSurface:
    """Blade-frame geometry: lofted grids (aero, shank) plus knuckle and heel."""

    aero: np.ndarray            # (Na, Np, 3)
    aero_r: np.ndarray
    shank: np.ndarray           # (Nk, Np, 3)
    shank_s: np.ndarray
    knuckle: np.ndarray         # (n_t, n_a, 3)
    heel: np.ndarray            # (M, 3)

    def clouds(self) -> dict[str, np.ndarray]:
        return {"aero": self.aero.reshape(-1, 3), "shank": self.shank.reshape(-1, 3),
                "knuckle": self.knuckle.reshape(-1, 3), "heel": self.heel}

    def all_points(self) -> np.ndarray:
        return np.vstack(list(self.clouds().values()))


def blade_surface(blade: BladeDef, lay: Layout, ax: BladeAxis, n_aero: int = 41,
                  n_shank: int = 16) -> BladeSurface:
    r = np.linspace(blade.r_root, blade.r_tip, n_aero)
    st = AeroStations(blade, r)
    ld = ax.lead_at(r)
    aero = st.points(s_of_r(r, ld, lay.hinge_r), ld, ax.aft_at(r), ax.slope_at(r, ld))
    sh, s_sh = shank_sections(blade, lay, ax, n_shank)
    return BladeSurface(aero, r, sh, s_sh, knuckle_points(lay), heel_points(lay))


# ==========================================================================
# Stowage pockets (static tail): raster of the stowed blade's inner surface
# ==========================================================================

def densify(G: np.ndarray, step: float = 0.25) -> np.ndarray:
    """Bilinear upsampling of a (Nu, Nv, 3) surface grid to <= step spacing."""
    Nu, Nv = G.shape[:2]
    du = np.linalg.norm(np.diff(G, axis=0), axis=-1).max() if Nu > 1 else 0.0
    dv = np.linalg.norm(np.diff(G, axis=1), axis=-1).max() if Nv > 1 else 0.0
    ku, kv = max(1, math.ceil(du / step)), max(1, math.ceil(dv / step))
    u = np.linspace(0, Nu - 1, (Nu - 1) * ku + 1)
    v = np.linspace(0, Nv - 1, (Nv - 1) * kv + 1)
    i0 = np.clip(np.floor(u).astype(int), 0, max(Nu - 2, 0))
    j0 = np.clip(np.floor(v).astype(int), 0, max(Nv - 2, 0))
    fu = (u - i0)[:, None, None]
    fv = (v - j0)[None, :, None]
    i1, j1 = np.minimum(i0 + 1, Nu - 1), np.minimum(j0 + 1, Nv - 1)
    a = G[i0][:, j0]
    b = G[i1][:, j0]
    c = G[i0][:, j1]
    d = G[i1][:, j1]
    return (1 - fu) * (1 - fv) * a + fu * (1 - fv) * b + (1 - fu) * fv * c + fu * fv * d


def _min_filter(A: np.ndarray, wi: int, wj: int) -> np.ndarray:
    """Min filter over a (2wi+1) x (2wj+1) window, edges padded (nan = empty)."""
    B = np.pad(np.where(np.isnan(A), np.inf, A), ((wi, wi), (wj, wj)), constant_values=np.inf)
    out = np.full(A.shape, np.inf)
    for di in range(2 * wi + 1):
        for dj in range(2 * wj + 1):
            out = np.minimum(out, B[di:di + A.shape[0], dj:dj + A.shape[1]])
    return np.where(np.isinf(out), np.nan, out)


@dataclass
class PocketMap:
    """Pocket floor raster for blade 0 (psi = psi0) on a (x, azimuth) grid.

    ``floor`` is the pocket floor radius (blade inner surface minus clearance,
    side clearance added by a min filter); nan where no pocket is needed. Other
    blades use the same map rotated by 360/N."""

    x: np.ndarray               # cell centres (mm)
    az: np.ndarray              # cell centres (deg, relative to the blade azimuth)
    floor: np.ndarray           # (nx, naz) pocket floor radius, nan = no pocket
    blade_floor: np.ndarray     # (nx, naz) raw inner-surface radius of the stowed blade
    blade_top: np.ndarray       # (nx, naz) raw outer-surface radius of the stowed blade

    def lookup(self, x, az):
        """Pocket floor radius at (x, az) (nearest cell), nan = no pocket."""
        i = np.clip(np.round((np.asarray(x) - self.x[0]) / (self.x[1] - self.x[0])).astype(int), 0, self.x.size - 1)
        j = np.clip(np.round((np.asarray(az) - self.az[0]) / (self.az[1] - self.az[0])).astype(int), 0, self.az.size - 1)
        inside = ((np.asarray(x) >= self.x[0] - 0.5 * (self.x[1] - self.x[0]))
                  & (np.asarray(x) <= self.x[-1] + 0.5 * (self.x[1] - self.x[0])))
        return np.where(inside, self.floor[i, j], np.nan)


def pocket_map(blade: BladeDef, lay: Layout, surf: BladeSurface, dx: float = 0.5,
               daz: float = 0.5, step: float = 0.25) -> PocketMap:
    """Rasterise the stowed blade 0 (aero + shank) onto an (x, az) grid over the
    static tail (x <= static_tail_x) and derive the pocket floor."""
    T = blade_to_world(lay.fold_deg, lay.psi0_deg, lay)
    parts = (("aero", surf.aero, lay.pocket_clear), ("shank", surf.shank, lay.root_clear))
    pts, clr = [], []
    for _, G, c in parts:
        W = apply(T, densify(G, step)).reshape(-1, 3)
        pts.append(W)
        clr.append(np.full(len(W), c))
    W = np.vstack(pts)
    clr = np.concatenate(clr)
    x, rho, az = cyl(W, lay.psi0_deg)
    x0 = math.floor(x.min()) - 2.0
    xs = np.arange(x0 + dx / 2, lay.static_tail_x + dx, dx)
    azs = np.arange(-90 + daz / 2, 90, daz)
    i = np.floor((x - x0) / dx).astype(int)
    j = np.floor((az + 90) / daz).astype(int)
    ok = (i >= 0) & (i < xs.size) & (j >= 0) & (j < azs.size)
    lin = i[ok] * azs.size + j[ok]
    bf = np.full(xs.size * azs.size, np.inf)
    ff = np.full(xs.size * azs.size, np.inf)
    bt = np.full(xs.size * azs.size, -np.inf)
    np.minimum.at(bf, lin, rho[ok])
    np.minimum.at(ff, lin, rho[ok] - clr[ok])
    np.maximum.at(bt, lin, rho[ok])
    bf, ff, bt = (a.reshape(xs.size, azs.size) for a in (bf, ff, bt))
    bf[np.isinf(bf)] = np.nan
    ff[np.isinf(ff)] = np.nan
    bt[np.isinf(bt)] = np.nan
    rh = hull_r(xs)[:, None]
    ff = np.where(ff < rh, ff, np.nan)                        # pocket only where below hull
    r_ref = np.maximum(rh, 5.0)
    wj = int(math.ceil(lay.side_clear / (np.radians(daz) * r_ref.min())))
    wi = int(math.ceil(lay.side_clear / dx))
    # side clearance: grow the footprint (min filter), but never beyond the hull
    fl = _min_filter(ff, wi, min(wj, 12))
    fl = np.where(fl < rh, fl, np.nan)
    fl = printable_floor(fl, rh[:, 0], dx, lay.floor_overhang_deg)
    return PocketMap(xs, azs, fl, bf, bt)


def printable_floor(fl: np.ndarray, rh: np.ndarray, dx: float, overhang_deg: float) -> np.ndarray:
    """Lower the pocket floor where it would overhang when printed. Both shells print with +x
    down (front shell on its equator face, aft shell on its tail end face), so going forward
    (up) the floor may rise by at most dx tan(overhang) per cell; faster rises (the root notch
    ramp, the pocket's forward end) become 45 deg ramps cut deeper into the shell. Scans each
    azimuth column from aft to forward; a floor can only get deeper (more blade clearance)."""
    k = dx * math.tan(math.radians(overhang_deg))
    out = fl.copy()
    for i in range(fl.shape[0] - 2, -1, -1):
        lim = out[i + 1] + k
        new = np.fmin(out[i], lim)
        out[i] = np.where(new < rh[i], new, np.nan)
    return out


def floor_overhang_max_deg(pm: "PocketMap") -> tuple[float, float]:
    """Steepest forward rise of the pocket floor (deg from vertical when printed +x down) and
    its x, over all azimuth columns where the floor exists in both cells."""
    d = pm.floor[:-1] - pm.floor[1:]                     # rise going forward by one cell
    dx = float(pm.x[1] - pm.x[0])
    ang = np.degrees(np.arctan(np.where(np.isfinite(d), np.maximum(d, 0.0), 0.0) / dx))
    i, _ = np.unravel_index(int(np.argmax(ang)), ang.shape)
    return float(ang.max()), float(pm.x[i])


def pocket_table(pm: PocketMap, lay: Layout, n_blades: int | None = None) -> dict:
    """Per-x pocket summary for the internal layout and the CAD (all radial, mm):
    depth (hull radius - floor), floor radius, keep-out radius (floor - MIN_WALL
    measured normal to the hull), centre azimuth and half-width (deg and mm of arc
    at the hull radius), land width between adjacent pockets."""
    n_blades = lay.n_blades if n_blades is None else n_blades
    rh = hull_r(pm.x)
    k = np.sqrt(1 + hull_slope(pm.x) ** 2)
    has = ~np.all(np.isnan(pm.floor), axis=1)
    out = {k_: np.full(pm.x.size, np.nan) for k_ in
           ("depth", "floor_r", "keepout_r", "az_lo", "az_hi", "az_c", "half_w_deg",
            "half_w_mm", "land_mm", "depth_normal")}
    for ii in np.nonzero(has)[0]:
        row = pm.floor[ii]
        m = ~np.isnan(row)
        fr = np.nanmin(row)
        lo, hi = pm.az[m].min() - 0.5 * (pm.az[1] - pm.az[0]), pm.az[m].max() + 0.5 * (pm.az[1] - pm.az[0])
        out["depth"][ii] = rh[ii] - fr
        out["depth_normal"][ii] = (rh[ii] - fr) / k[ii]
        out["floor_r"][ii] = fr
        out["keepout_r"][ii] = fr - P.MIN_WALL * k[ii]
        out["az_lo"][ii], out["az_hi"][ii] = lo, hi
        out["az_c"][ii] = 0.5 * (lo + hi)
        out["half_w_deg"][ii] = 0.5 * (hi - lo)
        out["half_w_mm"][ii] = np.radians(0.5 * (hi - lo)) * rh[ii]
        out["land_mm"][ii] = np.radians(360.0 / n_blades - (hi - lo)) * rh[ii]
    out["x"] = pm.x
    out["hull_r"] = rh
    return out


# ==========================================================================
# Hinge, hub and interference checks
# ==========================================================================

def hinge_coords(W: np.ndarray, lay: Layout, psi_deg: float | None = None):
    """World points -> hinge-plane coordinates (x, rho_e, tau) of the blade at psi:
    rho_e along e_r(psi), tau along e_t(psi) (the pin direction)."""
    psi = lay.psi0_deg if psi_deg is None else psi_deg
    _, e_r, e_t = rotor_axes(psi)
    return W[..., 0], W @ e_r, W @ e_t


def x_hub_front(lay: Layout) -> float:
    return lay.static_tail_x + lay.hub_gap


def slot_half_width(lay: Layout) -> float:
    return lay.knuckle_l / 2 + lay.ear_clear


def hinge_report(lay: Layout) -> dict:
    """Pin and ear envelope: how much hub material surrounds the pin and whether
    the ears fit under the hull (and under the spinner shell)."""
    pin_r = lay.pin_d / 2
    a = np.linspace(0, 2 * math.pi, 73)
    out = {}
    # pin surface, full length
    tau = np.linspace(-lay.pin_l / 2, lay.pin_l / 2, 41)
    A, Tt = np.meshgrid(a, tau, indexing="ij")
    _, e_r, e_t = rotor_axes(lay.psi0_deg)
    H = hinge_point(lay, lay.psi0_deg)
    def ring(rad):
        return (H + (rad * np.sin(A))[..., None] * np.array([1.0, 0, 0])
                + (rad * np.cos(A))[..., None] * e_r + Tt[..., None] * e_t)
    pin = ring(pin_r)
    out["pin_cover_to_hull_min_mm"] = float(-height_above_hull(pin).max())
    out["pin_cover_to_hub_core_min_mm"] = out["pin_cover_to_hull_min_mm"] - lay.spinner_wall - lay.spinner_gap
    out["pin_end_radius_max_mm"] = float(np.hypot(pin[..., 1], pin[..., 2]).max())
    # ear zone: pin_r + ear_cover around the axis, over the ear span (both ears)
    t_in = slot_half_width(lay)
    tau_e = np.concatenate([np.linspace(-lay.pin_l / 2, -t_in, 8), np.linspace(t_in, lay.pin_l / 2, 8)])
    A, Tt = np.meshgrid(a, tau_e, indexing="ij")
    ear = (H + ((pin_r + lay.ear_cover) * np.sin(A))[..., None] * np.array([1.0, 0, 0])
           + ((pin_r + lay.ear_cover) * np.cos(A))[..., None] * e_r + Tt[..., None] * e_t)
    h = height_above_hull(ear)
    out["ear_protrusion_over_hull_mm"] = float(max(h.max(), 0.0))
    out["ear_protrusion_over_hub_core_mm"] = float(max(h.max() + lay.spinner_wall, 0.0))
    out["ear_cover_available_min_mm"] = float(lay.ear_cover - max(h.max(), 0.0))
    out["pin_l_mm"] = lay.pin_l
    out["pin_l_needed_mm"] = lay.knuckle_l + 2 * lay.ear_clear + 2 * lay.ear_t
    out["hub_front_x_mm"] = x_hub_front(lay)
    out["pin_to_hub_front_mm"] = lay.hinge_x - x_hub_front(lay)
    out["hull_r_at_hinge_mm"] = float(hull_r(lay.hinge_x))
    return out


def _cloud_labels(surf: BladeSurface) -> tuple[np.ndarray, np.ndarray]:
    cl = surf.clouds()
    names = list(cl)
    pts = np.vstack([cl[k] for k in names])
    lab = np.concatenate([np.full(len(cl[k]), i) for i, k in enumerate(names)])
    return pts, lab, names


def _min_dist(A: np.ndarray, B: np.ndarray, chunk: int = 2000) -> float:
    best = np.inf
    for i in range(0, len(A), chunk):
        d = ((A[i:i + chunk, None, :] - B[None, :, :]) ** 2).sum(-1)
        best = min(best, float(d.min()))
    return math.sqrt(best)


def slot_profile(lay: Layout, surf: BladeSurface, folds: np.ndarray, cell: float = 0.2,
                 clear: float = 0.3) -> dict:
    """Swept profile, in the hinge plane (x, rho_e), of every blade point that lies
    between the ears (|tau| <= slot half width) over the fold range. Returns the
    occupancy grid (dilated by ``clear``) and derived hub numbers."""
    pts, lab, names = _cloud_labels(surf)
    xs = np.arange(lay.hinge_x - 30, lay.hinge_x + 20, cell)
    rs = np.arange(0.0, lay.hinge_r + 30, cell)
    occ = np.zeros((xs.size, rs.size), bool)
    hw = slot_half_width(lay)
    dense = np.vstack([densify(surf.shank, 0.4).reshape(-1, 3), densify(surf.knuckle, 0.4).reshape(-1, 3),
                       surf.heel, densify(surf.aero[:6], 0.4).reshape(-1, 3)])
    # fill the heel and knuckle interiors too (solid lugs)
    lam = np.linspace(0, lay.heel_reach, 30)
    w = np.linspace(-lay.heel_t / 2, lay.heel_t / 2, 12)
    d, q = heel_dir(lay)
    LL, WW = np.meshgrid(lam, w, indexing="ij")
    heel_fill = np.stack([LL * d[0] + WW * q[0], np.zeros_like(LL), LL * d[1] + WW * q[1]], -1).reshape(-1, 3)
    rr, aa = np.meshgrid(np.linspace(0, lay.knuckle_r, 12), np.linspace(0, 2 * math.pi, 60), indexing="ij")
    kn_fill = np.stack([rr * np.cos(aa), np.zeros_like(rr), rr * np.sin(aa)], -1).reshape(-1, 3)
    dense = np.vstack([dense, heel_fill, kn_fill])
    sel = np.abs(dense[:, 1]) <= hw
    dense = dense[sel]
    for f in folds:
        W = apply(blade_to_world(f, lay.psi0_deg, lay), dense)
        x, re, _ = hinge_coords(W, lay)
        i = np.floor((x - xs[0]) / cell).astype(int)
        j = np.floor(re / cell).astype(int)
        ok = (i >= 0) & (i < xs.size) & (j >= 0) & (j < rs.size)
        occ[i[ok], j[ok]] = True
    # dilate by the running clearance
    k = int(math.ceil(clear / cell))
    occd = occ.copy()
    for di in range(-k, k + 1):
        for dj in range(-k, k + 1):
            if di * di + dj * dj <= k * k:
                occd |= np.roll(np.roll(occ, di, 0), dj, 1)
    jj = np.nonzero(occd.any(axis=0))[0]
    floor = float(rs[jj.min()]) if jj.size else float("nan")
    return {"x": xs, "r": rs, "occ": occd, "slot_floor_r_mm": floor}


def stop_report(lay: Layout, prof: dict) -> dict:
    """Deployed stop: heel leading face against hub material. Measures the hub
    material ahead of the heel face along the heel's motion direction (thrust
    side), and the stop lever arm."""
    d, q = heel_dir(lay)
    lam_k = math.sqrt(max(lay.knuckle_r**2 - (lay.heel_t / 2) ** 2, 0.0))
    lam = np.linspace(lam_k, lay.heel_reach, 25)
    # leading face = -q side; under thrust (fold < 0) it moves along -q
    face = lam[:, None] * d[None, :] - (lay.heel_t / 2) * q[None, :]   # (s, n)
    x_face = lay.hinge_x + face[:, 1]
    r_face = lay.hinge_r + face[:, 0]
    dirn = np.array([-q[1], -q[0]])          # world (x, rho_e) direction of -q at fold 0
    xs, rs, occ = prof["x"], prof["r"], prof["occ"]
    cell = xs[1] - xs[0]
    xf = x_hub_front(lay)
    thick = []
    for x0, r0 in zip(x_face, r_face):
        t = 0.0
        while True:
            t += cell / 2
            xx, rr = x0 + t * dirn[0], r0 + t * dirn[1]
            if xx < xf or rr < 0 or rr > float(hull_r(xx)) - 1e-6:
                break
            i, j = int((xx - xs[0]) / cell), int(rr / cell)
            if 0 <= i < xs.size and 0 <= j < rs.size and occ[i, j] and t > 0.6:
                break
            if t > 30:
                break
        thick.append(t)
    thick = np.array(thick)
    good = thick >= 2.0
    lever = float(np.mean(lam[good])) if good.any() else float("nan")
    return {"stop_block_thickness_min_mm": float(thick.min()),
            "stop_block_thickness_max_mm": float(thick.max()),
            "stop_face_lambda_mm": [float(lam_k), float(lay.heel_reach)],
            "stop_contact_length_mm": float(np.ptp(lam[good])) if good.sum() > 1 else 0.0,
            "stop_lever_arm_mm": lever,
            "stop_face_x_mm": [float(x_face.min()), float(x_face.max())],
            "stop_face_r_mm": [float(r_face.min()), float(r_face.max())]}


def interference(blade: BladeDef, lay: Layout, surf: BladeSurface, pm: PocketMap,
                 d_fold: float = 1.0, nb_grid: int = 7) -> dict:
    """Swept-point checks over the fold range, for blade 0 (the others are copies):

    * static tail: every blade point forward of STATIC_TAIL_X must be outside the
      hull or above the pocket floor (floor already includes the design clearance,
      so >= 0 is the pass criterion; the stowed value equals the clearance);
    * ears: no blade point inside the pin ear zone;
    * hub envelope outside the slot: blade points that cut the hub/spinner beyond
      the knuckle width (the spinner slot must be widened there);
    * heel stowed and deployed inside the hub/spinner envelope;
    * neighbouring blades at any pair of fold angles; deployed blade vs hull."""
    folds = np.arange(0.0, lay.fold_deg + 1e-9, d_fold)
    pts, lab, names = _cloud_labels(surf)
    dense = {"aero": densify(surf.aero, 0.6).reshape(-1, 3),
             "shank": densify(surf.shank, 0.4).reshape(-1, 3),
             "knuckle": surf.knuckle.reshape(-1, 3), "heel": surf.heel}
    xf = x_hub_front(lay)
    hw = slot_half_width(lay)
    pin_zone = lay.pin_d / 2 + lay.ear_cover
    res = {"static_tail_clear_min_mm": np.inf, "static_tail_clear_fold_deg": None,
           "static_tail_clear_part": None, "ear_clear_min_mm": np.inf,
           }
    slot_tau = np.zeros(int(math.ceil(lay.hinge_x + 25 - xf)))
    for f in folds:
        T = blade_to_world(f, lay.psi0_deg, lay)
        for name, B in dense.items():
            W = apply(T, B)
            x, rho, az = cyl(W, lay.psi0_deg)
            m = x < lay.static_tail_x - 1e-6
            if m.any():
                fl = pm.lookup(x[m], az[m])
                allowed = np.where(np.isnan(fl), hull_r(x[m]), fl)
                c = rho[m] - allowed
                k = int(np.argmin(c))
                byp = res.setdefault("static_tail_clear_by_part_mm", {})
                byp[name] = min(byp.get(name, np.inf), float(c[k]))
                if c[k] < res["static_tail_clear_min_mm"]:
                    res.update(static_tail_clear_min_mm=float(c[k]), static_tail_clear_fold_deg=float(f),
                               static_tail_clear_part=name,
                               static_tail_clear_at=[float(x[m][k]), float(rho[m][k]), float(az[m][k])])
            xh, re, tau = hinge_coords(W, lay)
            # ear zone
            de = np.hypot(xh - lay.hinge_x, re - lay.hinge_r)
            ez = (np.abs(tau) > hw) & (np.abs(tau) <= lay.pin_l / 2)
            if ez.any():
                res["ear_clear_min_mm"] = min(res["ear_clear_min_mm"], float((de[ez] - pin_zone).min()))
            # hub envelope (x >= hub front, inside the hull): slot half-width needed per x
            hz = (xh >= xf) & (rho < hull_r(xh))
            if hz.any():
                ib = np.floor(xh[hz] - xf).astype(int)
                np.maximum.at(slot_tau, np.clip(ib, 0, slot_tau.size - 1), np.abs(tau[hz]))
    res["hub_slot_halfwidth_needed_mm"] = {"x_mm": (xf + 0.5 + np.arange(slot_tau.size)).tolist(),
                                           "tau_max_mm": np.where(slot_tau > 0, slot_tau, np.nan).tolist()}
    res["hub_slot_halfwidth_needed_max_mm"] = float(slot_tau.max())
    # heel envelope at both ends of travel
    for tag, f in (("deployed", 0.0), ("stowed", lay.fold_deg)):
        W = apply(blade_to_world(f, lay.psi0_deg, lay), surf.heel)
        x, rho, _ = cyl(W, lay.psi0_deg)
        res[f"heel_{tag}_below_spinner_inner_min_mm"] = float((hull_r(x) - lay.spinner_wall - rho).min())
        res[f"heel_{tag}_x_range_mm"] = [float(x.min()), float(x.max())]
        res[f"heel_{tag}_rho_range_mm"] = [float(rho.min()), float(rho.max())]
    # neighbours: root region only (aero is far apart), all fold pairs on a coarse grid
    root = np.vstack([surf.shank.reshape(-1, 3), surf.knuckle.reshape(-1, 3), surf.heel,
                      surf.aero[:10].reshape(-1, 3)])
    fgrid = np.linspace(0.0, lay.fold_deg, nb_grid)
    best = (np.inf, None)
    for fa in fgrid:
        A = apply(blade_to_world(fa, lay.psi(0), lay), root)
        for fb in fgrid:
            Bw = apply(blade_to_world(fb, lay.psi(1), lay), root)
            dmin = _min_dist(A, Bw)
            if dmin < best[0]:
                best = (dmin, (float(fa), float(fb)))
    res["neighbour_root_min_dist_mm"] = best[0]
    res["neighbour_root_min_folds_deg"] = best[1]
    Ad = apply(blade_to_world(0.0, lay.psi(0), lay), surf.aero.reshape(-1, 3))
    Bd = apply(blade_to_world(0.0, lay.psi(1), lay), surf.aero.reshape(-1, 3))
    res["deployed_neighbour_min_dist_mm"] = _min_dist(Ad, Bd)
    As = apply(blade_to_world(lay.fold_deg, lay.psi(0), lay), surf.aero.reshape(-1, 3))
    Bs = apply(blade_to_world(lay.fold_deg, lay.psi(1), lay), surf.aero.reshape(-1, 3))
    res["stowed_neighbour_min_dist_mm"] = _min_dist(As, Bs)
    res["deployed_aero_min_height_over_hull_mm"] = float(height_above_hull(Ad).min())
    return res


# ==========================================================================
# Heel (deployed stop lug) and hub slot
# ==========================================================================

SHAFT_BORE_CLEAR = 0.15     # diametral, "close running" (creative3dp ladder, data/hardware.json)
HUB_FLOOR_MIN = P.SHAFT_D / 2 + SHAFT_BORE_CLEAR / 2 + P.MIN_WALL   # lowest slot floor radius


def heel_reach_for(hinge_r: float, heel_t: float, run_clear: float = 0.3,
                   floor_min: float = HUB_FLOOR_MIN) -> float:
    """Longest heel (0.5 mm steps) whose leading tip corner, swept while folding,
    keeps the slot floor at ``floor_min`` (MIN_WALL of hub around the shaft bore)."""
    lim = hinge_r - run_clear - floor_min
    return math.floor(2.0 * math.sqrt(max(lim**2 - (heel_t / 2) ** 2, 0.0))) / 2.0


def heel_geometry(lay: Layout) -> dict:
    """Closed-form heel numbers (heel_angle 0: deployed it points straight down).
    The leading (thrust-side) face is flat from where it leaves the knuckle
    cylinder (lam_k) to the tip; its tip corner dips lowest at fold
    atan(t/2 / reach) and sets the slot floor."""
    t2 = lay.heel_t / 2
    lam_k = math.sqrt(max(lay.knuckle_r**2 - t2**2, 0.0))
    sweep_r = math.hypot(lay.heel_reach, t2)
    floor = lay.hinge_r - sweep_r - lay.run_clear
    return {"heel_reach_mm": lay.heel_reach, "heel_t_mm": lay.heel_t, "heel_w_mm": lay.heel_w,
            "flat_face_from_pin_mm": [lam_k, lay.heel_reach],
            "flat_face_length_mm": lay.heel_reach - lam_k,
            "swept_min_radius_mm": lay.hinge_r - sweep_r,
            "lowest_corner_fold_deg": math.degrees(math.atan2(t2, lay.heel_reach)),
            "slot_floor_r_mm": floor,
            "hub_wall_over_shaft_bore_mm": floor - P.SHAFT_D / 2 - SHAFT_BORE_CLEAR / 2}


def root_cloud(lay: Layout, surf: BladeSurface, n_aero_rows: int = 10, step: float = 0.4) -> np.ndarray:
    """Dense blade-frame cloud of the root: shank and first aero rows (surface),
    knuckle and heel (solid fill)."""
    lam = np.linspace(0, lay.heel_reach, 30)
    w = np.linspace(-lay.heel_t / 2, lay.heel_t / 2, 12)
    tw = np.linspace(-lay.heel_w / 2, lay.heel_w / 2, 25)
    d, q = heel_dir(lay)
    LL, WW, TT = np.meshgrid(lam, w, tw, indexing="ij")
    heel = np.stack([LL * d[0] + WW * q[0], TT, LL * d[1] + WW * q[1]], -1).reshape(-1, 3)
    rr, aa, tt = np.meshgrid(np.linspace(0, lay.knuckle_r, 10), np.linspace(0, 2 * math.pi, 60),
                             np.linspace(-lay.knuckle_l / 2, lay.knuckle_l / 2, 33), indexing="ij")
    kn = np.stack([rr * np.cos(aa), tt, rr * np.sin(aa)], -1).reshape(-1, 3)
    return np.vstack([densify(surf.shank, step).reshape(-1, 3),
                      densify(surf.aero[:n_aero_rows], step).reshape(-1, 3), heel, kn])


def hub_slot_table(lay: Layout, surf: BladeSurface, folds: np.ndarray, dx: float = 0.5) -> dict:
    """Swept envelope of blade 0's root inside the hub/spinner envelope (x >= hub
    front, rho < hull_r(x)) over the fold range, per x bin, in hinge-plane
    coordinates: lowest r_e (slot floor) and largest |tau| (slot half-width), with
    the running clearance added. This is the cut the hub and the spinner need for
    one blade; the other two are copies at +-120 deg."""
    cloud = root_cloud(lay, surf)
    xf = x_hub_front(lay)
    xs = np.arange(xf, lay.hinge_x + 30.0, dx)
    rmin = np.full(xs.size, np.inf)
    tmax = np.full(xs.size, -np.inf)
    for f in folds:
        W = apply(blade_to_world(f, lay.psi0_deg, lay), cloud)
        x, re, tau = hinge_coords(W, lay)
        rho = np.hypot(W[:, 1], W[:, 2])
        m = (x >= xf) & (rho < hull_r(x))
        i = np.clip(np.floor((x[m] - xf) / dx).astype(int), 0, xs.size - 1)
        np.minimum.at(rmin, i, re[m])
        np.maximum.at(tmax, i, np.abs(tau[m]))
    has = np.isfinite(rmin)
    return {"x_mm": (xs + dx / 2)[has], "slot_floor_re_mm": (rmin - lay.run_clear)[has],
            "slot_half_width_mm": (tmax + lay.run_clear)[has]}


def neighbour_check(lay: Layout, surf: BladeSurface, n: int = 10) -> dict:
    """Blade 0 vs blade 1 roots (shank, knuckle, heel, first aero rows) at every
    pair of fold angles on an n x n grid (blades fold independently)."""
    root = np.vstack([surf.shank.reshape(-1, 3), surf.knuckle.reshape(-1, 3), surf.heel,
                      surf.aero[:10].reshape(-1, 3)])
    g = np.linspace(0.0, lay.fold_deg, n)
    best = (np.inf, None)
    B = {fb: apply(blade_to_world(fb, lay.psi(1), lay), root) for fb in g}
    for fa in g:
        A = apply(blade_to_world(fa, lay.psi(0), lay), root)
        for fb in g:
            dmin = _min_dist(A, B[fb])
            if dmin < best[0]:
                best = (dmin, [float(fa), float(fb)])
    return {"neighbour_root_min_dist_mm": best[0], "at_folds_deg": best[1], "grid_n": n}


# ==========================================================================
# One-call evaluation (used by the sweep and the final run)
# ==========================================================================

def evaluate(blade: BladeDef, lay: Layout, ax: BladeAxis | None = None, d_fold: float = 3.0,
             n_axis: int = 21) -> dict:
    """Solve the axis (unless given) and return the decisive layout metrics."""
    ax = solve_axis(blade, lay, n=n_axis, iters=2) if ax is None else ax
    S = blade_surface(blade, lay, ax)
    pm = pocket_map(blade, lay, S)
    tb = pocket_table(pm, lay)
    hr = hinge_report(lay)
    hg = heel_geometry(lay)
    folds = np.arange(0.0, lay.fold_deg + 1e-9, d_fold)
    prof = slot_profile(lay, S, folds)
    sr = stop_report(lay, prof)
    it = interference(blade, lay, S, pm, d_fold=d_fold)
    x = tb["x"]
    tail = (x > lay.static_tail_x - 12) & ~np.isnan(tb["keepout_r"])
    aero = x < (lay.hinge_x - ax.s[0])
    byp = it.get("static_tail_clear_by_part_mm", {})
    rec = dict(HX=lay.hinge_x, HR=lay.hinge_r, front=lay.hinge_x - lay.static_tail_x,
               cone=lay.shank_cone_deg, fold=lay.fold_deg, XST=lay.static_tail_x,
               heel_t=lay.heel_t, heel_reach=lay.heel_reach,
               pin_cover_hull=hr["pin_cover_to_hull_min_mm"], pin_cover_core=hr["pin_cover_to_hub_core_min_mm"],
               depth_aero_max=float(np.nanmax(tb["depth"][aero])),
               keep_tail_min=float(np.nanmin(tb["keepout_r"][tail])),
               notch_floor=float(np.interp(lay.static_tail_x - 0.25, x, tb["floor_r"])),
               slot_floor=hg["slot_floor_r_mm"], hub_wall_over_bore=hg["hub_wall_over_shaft_bore_mm"],
               heel_flat=hg["flat_face_length_mm"], stop_block_max=sr["stop_block_thickness_max_mm"],
               stop_lever=sr["stop_lever_arm_mm"],
               clear_aero=byp.get("aero", np.nan), clear_shank=byp.get("shank", np.nan),
               ear_clear=it["ear_clear_min_mm"],
               heel_stow=it["heel_stowed_below_spinner_inner_min_mm"],
               heel_dep=it["heel_deployed_below_spinner_inner_min_mm"],
               nb_root=it["neighbour_root_min_dist_mm"],
               x_tip=float(x[~np.isnan(tb["depth"])].min()))
    return {"rec": rec, "axis": ax, "surf": S, "pm": pm, "table": tb, "hinge": hr, "heel": hg,
            "slot_profile": prof, "stop": sr, "interference": it}


FEASIBLE = dict(pin_cover_core=2.0, heel_stow=0.5, heel_flat=P.MAGNET["D"] + 0.6,
                ear_clear=0.5, nb_root=1.0, clear_shank=0.8, hub_wall_over_bore=P.MIN_WALL - 1e-6,
                stop_block_max=P.MAGNET["H"] + 0.2 + P.MIN_WALL + 0.6)


def feasible(rec: dict) -> bool:
    """Hard constraints of the hinge choice (see FEASIBLE): >= 2 mm of hub core
    around the pin, heel inside the spinner when stowed, room for a 4 mm magnet on
    the heel face and behind it in the stop block, ear and neighbour clearances,
    stowed root >= ROOT_CLEAR from the static tail while folding, MIN_WALL around
    the shaft bore."""
    return all(rec.get(k, -np.inf) >= v for k, v in FEASIBLE.items())


def hinge_sweep(grid, out: Path = WORK_DIR / "hinge_sweep.jsonl", blade: BladeDef | None = None,
                **lay_kw) -> list[dict]:
    """Evaluate (HX, HR, front) points; heel reach follows heel_reach_for(HR).
    Appends one JSON line per point to ``out`` and resumes from it."""
    blade = load_blade() if blade is None else blade
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = [json.loads(l) for l in out.read_text().splitlines()] if out.exists() else []
    done = {(r["HX"], r["HR"], r["front"], r["cone"]) for r in rows}
    base = chosen_layout(**lay_kw)
    for HX, HR, front in grid:
        if (HX, HR, front, base.shank_cone_deg) in done:
            continue
        lay = replace(base, hinge_x=HX, hinge_r=HR, static_tail_x=HX - front,
                      heel_reach=heel_reach_for(HR, base.heel_t, base.run_clear))
        ev = evaluate(blade, lay)
        rec = ev["rec"]
        rec["feasible"] = feasible(rec)
        rows.append(rec)
        with out.open("a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")
        print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
    return rows


# ==========================================================================
# The chosen layout (from the hinge sweep; see data/layout_rotor.json)
# ==========================================================================

CHOSEN = dict(hinge_x=116.0, hinge_r=14.0, static_tail_x=108.0, fold_deg=90.0,
              shank_cone_deg=55.0, knuckle_r=4.0, knuckle_l=14.0, ear_t=3.7, ear_clear=0.3,
              pin_l=22.0, heel_t=5.0, heel_reach=8.0, heel_w=10.0, heel_angle_deg=0.0)


def chosen_layout(**kw) -> Layout:
    """The Layout this study proposes (CHOSEN), with optional overrides."""
    d = dict(CHOSEN)
    d.update(kw)
    return Layout(**d)


# ==========================================================================
# Blade mass properties (estimate: solid print)
# ==========================================================================

def _planar_area_centroid(poly: np.ndarray) -> tuple[float, np.ndarray]:
    """Area and centroid of a planar 3D polygon (triangle fan)."""
    p0 = poly[0]
    a = poly[1:-1] - p0
    b = poly[2:] - p0
    cr = np.cross(a, b)
    n = cr.sum(axis=0)
    nn = np.linalg.norm(n)
    if nn == 0:
        return 0.0, poly.mean(axis=0)
    w = cr @ (n / nn) / 2.0
    cen = p0 + (a + b) / 3.0                       # triangle centroids
    A = float(w.sum())
    return abs(A), (w[:, None] * cen).sum(axis=0) / A


def blade_mass_props(lay: Layout, surf: BladeSurface, density_g_mm3: float = PLA_DENSITY) -> dict:
    """Volume, mass and centre of gravity (blade frame) of one blade from its
    section loft (aero + shank), knuckle cylinder (minus pin bore) and heel."""
    vols, cgs = [], []
    for G in (surf.aero, surf.shank):
        A, C = zip(*(_planar_area_centroid(sec) for sec in G))
        A, C = np.array(A), np.array(C)
        dl = np.linalg.norm(np.diff(C, axis=0), axis=1)
        Am = 0.5 * (A[1:] + A[:-1])
        Cm = 0.5 * (C[1:] + C[:-1])
        vols.append(float((Am * dl).sum()))
        cgs.append((Am * dl)[:, None] * Cm)
    v_k = math.pi * (lay.knuckle_r**2 - (lay.pin_d / 2) ** 2) * lay.knuckle_l
    lam_k = math.sqrt(max(lay.knuckle_r**2 - (lay.heel_t / 2) ** 2, 0.0))
    v_h = (lay.heel_reach - lam_k) * lay.heel_t * lay.heel_w
    d, _ = heel_dir(lay)
    c_h = 0.5 * (lay.heel_reach + lam_k) * np.array([d[0], 0.0, d[1]])
    V = sum(vols) + v_k + v_h
    cg = (cgs[0].sum(0) + cgs[1].sum(0) + v_h * c_h) / V
    m = V * density_g_mm3
    d_axis = float(math.hypot(cg[0], cg[2]))
    return {"volume_mm3": V, "volume_aero_mm3": vols[0], "volume_shank_mm3": vols[1],
            "volume_knuckle_mm3": v_k, "volume_heel_mm3": v_h, "density_g_cm3": density_g_mm3 * 1e3,
            "mass_g": m, "cg_blade_frame_mm": cg.tolist(), "cg_distance_from_pin_axis_mm": d_axis,
            "gravity_moment_max_Nmm": m * 1e-3 * 9.81 * d_axis,
            "note": "ESTIMATE: solid PLA (1.24 g/cm3, same assumption as the BEM study); "
                    "worst-case gravity moment about the pin = m g x cg distance from the pin axis"}


# ==========================================================================
# Magnet detents: charge-model force between two coaxial disc magnets
# ==========================================================================

MU0 = 4e-7 * math.pi


def _disc_coupling_m2(z_m: float, a_m: float) -> float:
    """I(z) = pi a^2 int_0^inf J1(ka)^2 exp(-kz) / k dk: the interaction integral of two
    coaxial uniformly charged discs of radius a at distance z (I(0) = pi a^2 / 2).
    Uses scipy (lazy import)."""
    if z_m <= 0:
        return math.pi * a_m**2 / 2
    from scipy.integrate import quad
    from scipy.special import j1
    zr = z_m / a_m
    f = lambda q: j1(q) ** 2 * math.exp(-q * zr) / q  # noqa: E731
    q_end = max(40.0, 40.0 / zr)
    edges = np.concatenate([[1e-12], np.arange(math.pi, q_end, math.pi), [q_end]])
    val = sum(quad(f, a, b, limit=100)[0] for a, b in zip(edges[:-1], edges[1:]))
    return math.pi * a_m**2 * val


def magnet_pair_force_N(gap_mm: float, D_mm: float = P.MAGNET["D"], H_mm: float = P.MAGNET["H"],
                        Br_T: float = 1.43) -> float:
    """Attraction (N) between two identical coaxial axially magnetised disc magnets
    face to face at ``gap_mm`` (surface-charge model, mu_r = 1). With the method of
    images, F(gap = 0) also models the pull of one magnet on a thick steel plate."""
    a, h, g = D_mm / 2000, H_mm / 1000, gap_mm / 1000
    I = lambda z: _disc_coupling_m2(z, a)  # noqa: E731
    return Br_T**2 / MU0 * (I(g) - 2 * I(g + h) + I(g + 2 * h))


def magnet_calibration() -> dict:
    """Knock-down factor k = catalogue pull to a steel plate / model pull at contact,
    from data/hardware.json (vendor pull figures) when present."""
    hw = {}
    try:
        hw = json.loads((ROOT / "data" / "hardware.json").read_text()).get("magnets", {})
    except Exception:  # noqa: BLE001
        pass
    pulls = []
    m4 = hw.get("n52_4x2", {}).get("pull", {})
    if isinstance(m4.get("value"), list):
        for v in m4["value"]:
            pulls.append(("4x2", 4.0, 2.0, float(v), f"{m4.get('source_url', '')}; {m4.get('note', '')}"))
    for key, D, H in (("n52_3x2_pull", 3.0, 2.0), ("n52_5x2_pull", 5.0, 2.0)):
        e = hw.get(key, {})
        if isinstance(e.get("value"), (int, float)):
            pulls.append((key[4:7], D, H, float(e["value"]), e.get("source_url", "")))
    if not pulls:   # fallback: jc-magnetics listing, 1.03 lb magnet to steel plate
        pulls = [("4x2", 4.0, 2.0, 4.58, "https://www.jc-magnetics.com/Magnet-N52-4mm-2mm-Discs")]
    rows = []
    for name, D, H, F, src in pulls:
        Fm = magnet_pair_force_N(0.0, D, H, 1.43)
        rows.append({"magnet": name, "D_mm": D, "H_mm": H, "catalogue_pull_N": F,
                     "model_pull_Br1p43_N": Fm, "k": F / Fm, "source": src})
    k4 = [r["k"] for r in rows if r["magnet"] == "4x2"]
    return {"rows": rows, "k_design": min(k4) if k4 else min(r["k"] for r in rows),
            "Br_T_model": 1.43,
            "note": "k_design = lowest 4x2 catalogue pull / model; applied to every model force"}


STOW_MAG_DX = 3.5      # stowed-detent magnet: this far forward of the static tail end face


def detents(lay: Layout, ax: BladeAxis, surf: BladeSurface, pm: PocketMap, tb: dict,
            mass: dict, cal: dict, stop: dict) -> dict:
    """Two magnet pairs per blade (4x2 N52, axially magnetised, faces recessed
    ``magnet_recess`` below the contact faces so plastic takes the contact):

    * deployed: blade magnet in the heel's leading face, hub magnet in the stop face;
      the pair closes when the heel lands on its stop and resists forward folding.
    * stowed: blade magnet in the root's inner (pressure) face just forward of the
      static tail end, static magnet in a boss on the notch floor; the boss is the
      stowed rest, and the pair also indexes the rotor in the bag.
    """
    k = cal["k_design"]
    D, H = P.MAGNET["D"], P.MAGNET["H"]
    gap0 = 2 * lay.magnet_recess
    gaps = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
    F_gap = [k * magnet_pair_force_N(g, D, H) for g in gaps]
    F0 = k * magnet_pair_force_N(gap0, D, H)
    hg = heel_geometry(lay)
    lam0, lam1 = hg["flat_face_from_pin_mm"]
    lam_m = 0.5 * (lam0 + lam1)
    dep = {"blade_magnet_blade_frame_mm": [-lam_m, 0.0, -lay.heel_t / 2 + lay.magnet_recess + H / 2],
           "blade_magnet_axis_blade_frame": [0.0, 0.0, -1.0],
           "hub_magnet_hinge_plane_x_re_tau_mm": [lay.hinge_x - lay.heel_t / 2 - lay.magnet_recess - H / 2,
                                                  lay.hinge_r - lam_m, 0.0],
           "hub_magnet_axis_world": "+x (faces the heel)",
           "edge_margin_on_heel_face_mm": (lam1 - lam0 - D) / 2,
           "gap_mm": gap0, "force_N": F0, "lever_mm": lam_m, "torque_Nmm": F0 * lam_m,
           "stop_face_x_mm": lay.hinge_x - lay.heel_t / 2}
    # stowed: inner face of the root at x_m on the notch centre azimuth
    x_m = lay.static_tail_x - STOW_MAG_DX
    i = int(np.argmin(np.abs(pm.x - x_m)))
    row = pm.blade_floor[i]
    j = int(np.nanargmin(row))
    r_in, az = float(row[j]), float(pm.az[j])
    t_here = float(pm.blade_top[i, j] - pm.blade_floor[i, j])
    psi = math.radians(lay.psi0_deg + az)
    p = np.array([pm.x[i], -r_in * math.sin(psi), r_in * math.cos(psi)])
    e_rho = np.array([0.0, -math.sin(psi), math.cos(psi)])
    T = blade_to_world(lay.fold_deg, lay.psi0_deg, lay)
    pb = (p - T[:3, 3]) @ T[:3, :3]                      # world -> blade frame
    _, _, e_t = rotor_axes(lay.psi0_deg)
    H0 = hinge_point(lay, lay.psi0_deg)
    lever = float(np.cross(p - H0, -e_rho) @ e_t)        # moment arm of a unit inward pull
    stow = {"x_mm": float(pm.x[i]), "az_deg": az, "blade_inner_r_mm": r_in,
            "blade_magnet_centre_blade_frame_mm": (pb + (lay.magnet_recess + H / 2) * (e_rho @ T[:3, :3])).tolist(),
            "blade_magnet_axis_blade_frame": (-(e_rho @ T[:3, :3])).tolist(),
            "blade_thickness_at_magnet_mm": t_here,
            "blade_skin_behind_magnet_mm": t_here - lay.magnet_recess - H,
            "boss_top_r_mm": r_in,
            "boss_height_over_pocket_floor_mm": (float(r_in - pm.floor[i, j]) if np.isfinite(pm.floor[i, j])
                                                 else lay.root_clear),
            "static_magnet_centre_r_mm": r_in - lay.magnet_recess - H / 2,
            "static_keepout_r_mm": r_in - lay.magnet_recess - H - P.MIN_WALL,
            "gap_mm": gap0, "force_N": F0, "lever_mm": lever, "torque_Nmm": F0 * lever}
    Mg = mass["gravity_moment_max_Nmm"]
    r_tip = float(ax.r_tip)
    return {"magnet": f"{D:g} x {H:g} mm N52 disc, axial (MAGNET in params.py)",
            "per_blade": "2 in the blade (heel face, root inner face), 1 in the hub stop face, "
                         "1 in the static tail notch floor; 12 in total",
            "polarity": "every pair attracting; all blades identical, all hub and all tail magnets identical",
            "calibration": cal, "force_vs_gap": {"gap_mm": gaps, "force_N": F_gap},
            "deployed": dep, "stowed": stow,
            "gravity_moment_max_Nmm": Mg,
            "safety_vs_gravity_deployed": dep["torque_Nmm"] / Mg,
            "safety_vs_gravity_stowed": stow["torque_Nmm"] / Mg,
            "finger_force_at_tip_to_fold_N": dep["torque_Nmm"] / (r_tip - lay.hinge_r),
            "finger_force_at_tip_to_deploy_N": stow["torque_Nmm"] / (r_tip - lay.hinge_r),
            "note": "Forces: surface-charge model scaled by k (catalogue pull / model). Contact faces "
                    "are plastic; magnets sit magnet_recess below them. A weaker stowed hold for "
                    "flow-opening in water: deepen the recesses (see force_vs_gap)."}


# ==========================================================================
# Hub concept and loads
# ==========================================================================

def bem_loads(path: Path = BLADE_JSON, lay: "Layout | None" = None) -> dict:
    """Flap moment and thrust per blade from data/blade_design.json (BEM study).

    Prefers the BEM's layout_axis_check (final blade flown on this study's hinge and solved
    curved pitch axis) when its hinge matches `lay`; otherwise falls back to the nominal-axis
    envelope (performance.envelope_up_to_runaway), which assumes the params.py hinge."""
    out = {"source": None}
    try:
        d = json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return out
    nb = d["rotor"]["n_blades"]
    env = d.get("performance", {}).get("envelope_up_to_runaway", [])
    se = d.get("structure_estimates", {})
    torque = {f"{e['flow_m_s']:g}": e["max_torque_Nm"] for e in env}
    lac = d.get("layout_axis_check") or {}
    match = (lay is not None and lac.get("rows")
             and abs(lac.get("hinge_x_mm", -1e9) - lay.hinge_x) < 1e-6
             and abs(lac.get("hinge_r_mm", -1e9) - lay.hinge_r) < 1e-6)
    rows = {}
    if match:
        for e in lac["rows"]:
            tag = f"{e['flow_m_s']:g}"
            rows[tag] = {"flap_moment_per_blade_Nm": e["max_flap_moment_hinge_Nm"],
                         "thrust_per_blade_N": e["max_thrust_N"] / nb,
                         "rotor_torque_Nm": torque.get(tag),
                         "flap_moment_nominal_axis_Nm": e["max_flap_moment_hinge_nominal_axis_Nm"],
                         "cp_peak": e["cp_peak"], "cp_peak_nominal_axis": e["cp_peak_nominal_axis"]}
        top = max(lac["rows"], key=lambda e: e["flow_m_s"])
        out.update(source=f"{path.name}: layout_axis_check (hinge {lac['hinge_x_mm']:g}, {lac['hinge_r_mm']:g}, "
                          "solved curved axis); torque from performance.envelope_up_to_runaway",
                   axis="curved (layout_axis_check)", max_flap_moment_Nm=top["max_flap_moment_hinge_Nm"],
                   max_flap_moment_nominal_axis_Nm=top["max_flap_moment_hinge_nominal_axis_Nm"],
                   max_thrust_N=top["max_thrust_N"], max_flow_m_s=top["flow_m_s"],
                   cp_peak=top["cp_peak"], cp_peak_nominal_axis=top["cp_peak_nominal_axis"])
    else:
        for e in env:
            rows[f"{e['flow_m_s']:g}"] = {"flap_moment_per_blade_Nm": e["max_flap_moment_hinge_per_blade_Nm"],
                                           "thrust_per_blade_N": e["max_thrust_N"] / nb,
                                           "rotor_torque_Nm": e["max_torque_Nm"]}
        out.update(source=f"{path.name}: performance.envelope_up_to_runaway, structure_estimates",
                   axis="nominal (params.py hinge, straight axis)",
                   max_flap_moment_Nm=se.get("max_flap_moment_2p5_m_s_Nm"),
                   max_thrust_N=se.get("max_thrust_2p5_m_s_N"))
    out.update(by_flow_m_s=rows, bem_hinge_r_mm=d["rotor"].get("hinge_r_mm"),
               bem_hinge_x_mm=d["rotor"].get("hinge_x_mm"), curved_axis_check_used=bool(match))
    return out


def hub_concept(lay: Layout, surf: BladeSurface, ev: dict, slot: dict, loads: dict) -> dict:
    """Hub, hinge, stop and spinner proposal with the numbers the CAD needs, plus
    first-order load checks (labelled estimates)."""
    hr, hg, sr = ev["hinge"], ev["heel"], ev["stop"]
    xf = x_hub_front(lay)
    x_end = lay.hinge_x + lay.shaft_end_dx
    psi_g = lay.psi0_deg + 180.0 / lay.n_blades           # between blade 0 and blade 1
    x_g = xf + 3.0
    core_r = lambda x: hull_r(x) - lay.spinner_wall - lay.spinner_gap  # noqa: E731
    # grub screw hole (M3, 2.5 tap drill) clearance to the swept roots of blades 0 and 1
    rr = np.linspace(P.SHAFT_D / 2, float(core_r(x_g)), 40)
    pg = np.array([0.0, -math.sin(math.radians(psi_g)), math.cos(math.radians(psi_g))])
    seg = np.array([x_g, 0.0, 0.0]) + rr[:, None] * pg
    cloud = root_cloud(lay, surf)
    dmin = np.inf
    for i in (0, 1):
        for f in np.arange(0.0, lay.fold_deg + 1e-9, 5.0):
            dmin = min(dmin, _min_dist(seg, apply(blade_to_world(f, lay.psi(i), lay), cloud)))
    # pin insertion channel along e_t from the pin end to the hub core surface
    H0 = hinge_point(lay, lay.psi0_deg)
    _, e_r, e_t = rotor_axes(lay.psi0_deg)
    tau_exit = math.sqrt(max(float(core_r(lay.hinge_x)) ** 2 - lay.hinge_r**2, 0.0))
    ch = H0 + np.linspace(lay.pin_l / 2, tau_exit, 20)[:, None] * e_t
    ch2 = H0 - np.linspace(lay.pin_l / 2, tau_exit, 20)[:, None] * e_t
    dch = np.inf
    for i in (1, 2):
        for f in np.arange(0.0, lay.fold_deg + 1e-9, 5.0):
            Wn = apply(blade_to_world(f, lay.psi(i), lay), cloud)
            dch = min(dch, _min_dist(ch, Wn), _min_dist(ch2, Wn))
    # loads: stop, pin, ears (first-order statics)
    lever = sr["stop_lever_arm_mm"]
    cases = {}
    for tag, row in loads.get("by_flow_m_s", {}).items():
        M, T = row["flap_moment_per_blade_Nm"], row["thrust_per_blade_N"]
        S = M / (lever * 1e-3)
        R = S + T
        cases[tag] = {"flap_moment_Nm": M, "thrust_per_blade_N": T, "stop_force_N": S,
                      **({"flap_moment_nominal_axis_Nm": row["flap_moment_nominal_axis_Nm"]}
                         if "flap_moment_nominal_axis_Nm" in row else {}),
                      "stop_contact_pressure_MPa": S / ((lay.heel_reach - hg["flat_face_from_pin_mm"][0]) * lay.heel_w),
                      "pin_reaction_N": R,
                      "pin_double_shear_MPa": R / (2 * math.pi * lay.pin_d**2 / 4),
                      "ear_bearing_MPa": R / (2 * lay.ear_t * lay.pin_d),
                      "knuckle_bearing_MPa": R / (lay.knuckle_l * lay.pin_d),
                      "knuckle_bending_MPa": M / ((lay.knuckle_l * (2 * lay.knuckle_r) ** 2 / 6
                                                   - lay.knuckle_l * lay.pin_d**2 / 6) * 1e-9) * 1e-6}
    return {
        "concept": [
            "Hub (printed, rotating): body of revolution from the hub front (static tail end + HUB_GAP) "
            "to the shaft end, outer surface = hull - spinner wall - gap. Bore 8 mm close-running with a "
            "D-flat; one M3 cup-point grub (ISO 4029) radial on the D-flat, between two blades.",
            "Clevis per blade: the blade root ends in a full-width knuckle (cylinder about the pin) that "
            "sits in a slot between two hub ears; the 3 mm pin runs TANGENTIALLY through ear, knuckle, ear.",
            "Pins: ISO 8734 parallel pins pressed into the ears (printed 3.0 hole), knuckle bore 3.0 + "
            "CLEAR_ROTATE. They go in along e_t through a 3 mm channel from the hub core surface; the "
            "spinner shell covers both channel ends and retains the pins (folding-prop practice).",
            "Deployed stop: a heel tab under the knuckle; under thrust its leading face lands on the hub "
            "front ring (stop face at x = hinge_x - heel_t/2), which takes the flapwise moment.",
            "Spinner (printed, rotating): 1.2 mm shell from the hub front to the blunt tip with three slots "
            "(slot table). Slides on from aft, keys on the slots, and is held by one axial M3 screw into a "
            "tapped hole in the shaft end; the hub's aft face bears on the spinner's inner aft wall, so "
            "rotor thrust goes hub -> spinner -> M3 screw -> shaft (the grub only locks rotation and "
            "holds the hub during assembly).",
            "Detents: N52 4x2 magnet pairs, heel face <-> stop face (deployed) and root inner face <-> "
            "boss on the static notch floor (stowed). No springs.",
        ],
        "hub_front_x_mm": xf, "hub_aft_face_x_mm": x_end, "shaft_end_x_mm": x_end,
        "shaft_engagement_mm": x_end - xf,
        "hub_core_r_at_hinge_mm": float(core_r(lay.hinge_x)),
        "bore_d_mm": P.SHAFT_D + SHAFT_BORE_CLEAR, "d_flat_depth_mm": 0.5,
        "grub": {"size": "M3 x 4 cup point (ISO 4029)", "x_mm": x_g, "azimuth_psi_deg": psi_g,
                 "hole_to_swept_roots_min_mm": dmin - 1.25,
                 "note": "radial M3 tapped hole (2.5 drill) at an azimuth between two blades; "
                         "clearance = hole axis to swept root surface minus the 1.25 mm drill radius"},
        "ears": {"thickness_along_pin_mm": lay.ear_t, "inner_face_tau_mm": slot_half_width(lay),
                 "outer_face_tau_mm": lay.pin_l / 2, "material_around_pin_in_core_min_mm": hr["pin_cover_to_hub_core_min_mm"],
                 "knuckle_side_clearance_mm": lay.ear_clear},
        "pin": {"d_mm": lay.pin_d, "l_mm": lay.pin_l, "standard": "ISO 8734 parallel pin 3 x 22 (VERIFY stock length)",
                "ear_hole_d_mm": lay.pin_d, "knuckle_bore_d_mm": lay.pin_d + P.CLEAR_ROTATE,
                "insertion_channel_tau_mm": [lay.pin_l / 2, tau_exit],
                "channel_to_neighbour_swept_roots_min_mm": dch - lay.pin_d / 2,
                "retention": "spinner shell covers both channel ends"},
        "knuckle": {"radius_mm": lay.knuckle_r, "length_mm": lay.knuckle_l},
        "heel": hg,
        "stop": {**sr, "stop_face_x_mm_nominal": lay.hinge_x - lay.heel_t / 2,
                 "stop_block_x_mm": [xf, lay.hinge_x - lay.heel_t / 2]},
        "spinner": {"wall_mm": lay.spinner_wall, "radial_gap_to_hub_mm": lay.spinner_gap,
                    "front_x_mm": xf, "tip_r_mm": P.SPINNER_TIP_R,
                    "retention": "axial M3 screw into a tapped hole in the shaft end (M3 tap, 8 deep)",
                    "slot": "per blade: hub_slot table (x, slot floor r_e, half-width tau) + run clearance"},
        "loads_source": loads.get("source"),
        "loads_axis": loads.get("axis"),
        "load_cases": cases,
        "load_note": "ESTIMATES, first-order statics: stop force = flap moment / stop lever; pin reaction = "
                     "stop force + blade thrust (same direction); stresses on projected areas. Compare "
                     "with PETG/PLA strength in the report; the 2.5 m/s case is the BEM worst case.",
    }


# ==========================================================================
# Fold-angle study and twist effect
# ==========================================================================

def cone_stats(ax: BladeAxis) -> dict:
    """Deployed cone implied by the solved axis."""
    s0, s1 = float(ax.s[0]), float(ax.s[-1])
    a0, a1 = float(ax.aft[0]), float(ax.aft[-1])
    c = ax.cone_deg()
    return {"cone_tip_local_deg": float(c[-1]), "cone_root_local_deg": float(c[0]),
            "cone_mean_hinge_to_tip_deg": math.degrees(math.atan2(a1, s1)),
            "cone_mean_aero_span_deg": math.degrees(math.atan2(a1 - a0, s1 - s0)),
            "aft_max_mm": float(ax.aft.max()), "aft_tip_mm": a1,
            "r_at_aft_max_mm": float(ax.r[int(np.argmax(ax.aft))])}


def fold_study(blade: BladeDef, lay: Layout, folds=(84.0, 87.0, 90.0, 93.0, 96.0)) -> list[dict]:
    """Stowed fit and deployed cone vs the stowed fold angle (axis re-solved each time)."""
    rows = []
    for f in folds:
        lf = replace(lay, fold_deg=float(f))
        ax = solve_axis(blade, lf, n=21, iters=2)
        S = blade_surface(blade, lf, ax)
        pm = pocket_map(blade, lf, S)
        tb = pocket_table(pm, lf)
        x = tb["x"]
        aero = x < (lf.hinge_x - ax.s[0])
        tail = (x > lf.static_tail_x - 12) & ~np.isnan(tb["keepout_r"])
        W = apply(blade_to_world(f, lf.psi0_deg, lf), S.heel)
        xh, rho, _ = cyl(W, lf.psi0_deg)
        rows.append({"fold_deg": float(f), **cone_stats(ax),
                     "pocket_depth_aero_max_mm": float(np.nanmax(tb["depth"][aero])),
                     "keepout_tail_min_mm": float(np.nanmin(tb["keepout_r"][tail])),
                     "stowed_tip_x_mm": float(x[~np.isnan(tb["depth"])].min()),
                     "heel_stowed_below_spinner_inner_mm": float((hull_r(xh) - lf.spinner_wall - rho).min())})
    return rows


def stowed_sections(blade: BladeDef, lay: Layout, ax: BladeAxis) -> dict:
    """Per aero station, stowed: pitch-point x, hull radius, outer surface height over
    the hull (max), deepest point below the hull (= flush pocket depth), section
    thickness, flat-chord sagitta on the round hull, twist rise c sin(twist)."""
    st = AeroStations(blade, ax.r)
    W = apply(blade_to_world(lay.fold_deg, lay.psi0_deg, lay), st.points(ax.s, ax.lead, ax.aft, ax.slope))
    h = height_above_hull(W)
    xp = apply(blade_to_world(lay.fold_deg, lay.psi0_deg, lay), np.stack([ax.s, ax.lead, ax.aft], -1))[:, 0]
    R = hull_r(xp)
    c = st.c
    up = h[:, :N_PER_SIDE]      # suction side = outer face when stowed
    return {"r_mm": ax.r, "x_stowed_mm": xp, "hull_r_mm": R, "chord_mm": c, "twist_deg": st.tw,
            "thickness_mm": st.tc * c, "top_over_hull_mm": h.max(axis=1), "depth_below_hull_mm": -h.min(axis=1),
            "outer_face_lowest_over_hull_mm": up.min(axis=1),
            "sagitta_flat_chord_mm": R - np.sqrt(np.maximum(R**2 - (c / 2) ** 2, 0)),
            "twist_rise_mm": c * np.sin(np.radians(st.tw)), "lead_mm": ax.lead, "aft_mm": ax.aft,
            "cone_local_deg": ax.cone_deg(),
            "yaw_vs_local_tangent_deg": np.degrees(np.arctan2(ax.lead, ax.r))}


def twist_effect(blade: BladeDef, lay: Layout, ax: BladeAxis) -> dict:
    """Flush pocket depth per station for the real blade, the same blade with radial
    stacking (lead = 0), and an untwisted blade (lead optimised): isolates what the
    twist costs and what the in-plane lead buys back."""
    base = stowed_sections(blade, lay, ax)
    ax0 = solve_axis(blade, lay, n=ax.r.size, lead=0.0, iters=2)
    no_lead = stowed_sections(blade, lay, ax0)
    flat = replace(blade, twist_tab=tuple(np.zeros(len(blade.twist_tab))) if blade.twist_tab else None,
                   twist_root=0.0, twist_tip=0.0)
    axf = solve_axis(flat, lay, n=ax.r.size, iters=2)
    no_tw = stowed_sections(flat, lay, axf)
    return {"r_mm": ax.r, "depth_mm": base["depth_below_hull_mm"],
            "depth_radial_stacking_mm": no_lead["depth_below_hull_mm"],
            "depth_untwisted_mm": no_tw["depth_below_hull_mm"],
            "max_depth_mm": float(base["depth_below_hull_mm"].max()),
            "max_depth_radial_stacking_mm": float(no_lead["depth_below_hull_mm"].max()),
            "max_depth_untwisted_mm": float(no_tw["depth_below_hull_mm"].max()),
            "cone_radial_stacking": cone_stats(ax0)}


# ==========================================================================
# Figures
# ==========================================================================

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
C_STOW, C_DEP, C_MID = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3
HULL_FILL, HUB_FILL = "#e6e5e1", "#f2f1ee"


def _blade_world(lay: Layout, surf: BladeSurface, fold: float, i: int = 0, aero_step: float = 0.5):
    """Dense world cloud of blade i (aero surface + solid-filled root) at a fold angle."""
    T = blade_to_world(fold, lay.psi(i), lay)
    return apply(T, np.vstack([densify(surf.aero, aero_step).reshape(-1, 3), root_cloud(lay, surf)]))


def _raster_outline(axp, u, v, cell, color, label=None, fill=True, alpha=0.9, lw=1.0, zorder=3):
    """Silhouette of a 2D point set: rasterise, close small gaps, fill holes, then
    draw it as a filled region or an outline (scipy.ndimage, lazy import)."""
    from scipy import ndimage
    u0, v0 = u.min() - 3 * cell, v.min() - 3 * cell
    nu = int((u.max() - u0) / cell) + 4
    nv = int((v.max() - v0) / cell) + 4
    img = np.zeros((nv, nu), bool)
    img[((v - v0) / cell).astype(int), ((u - u0) / cell).astype(int)] = True
    st = ndimage.generate_binary_structure(2, 1)
    img = ndimage.binary_closing(img, structure=st, iterations=2)
    img = ndimage.binary_fill_holes(img)
    U, V = np.meshgrid(u0 + (np.arange(nu) + 0.5) * cell, v0 + (np.arange(nv) + 0.5) * cell)
    if fill:
        axp.contourf(U, V, img.astype(float), levels=[0.5, 1.5], colors=[color], alpha=alpha, zorder=zorder)
        if label:
            axp.fill([], [], color=color, alpha=alpha, label=label)
    else:
        axp.contour(U, V, img.astype(float), levels=[0.5], colors=[color], linewidths=lw, zorder=zorder)
        if label:
            axp.plot([], [], color=color, lw=lw, label=label)


def figure_side(lay: Layout, surf: BladeSurface, tb: dict, slot: dict, det: dict, path: Path) -> None:
    """Side view in the meridional (x, rho) half-plane of blade 0: hull, static tail
    end, hub/spinner, pocket floor and keep-out, stowed, intermediate and deployed
    blade silhouettes, hinge, knuckle, heel and magnets. Full view plus hinge zoom."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    fig, axs = plt.subplots(2, 1, figsize=(10.5, 11.5), gridspec_kw={"height_ratios": [1.0, 1.05]})
    xx = np.linspace(0, P.HULL_L, 600)
    W_st = _blade_world(lay, surf, lay.fold_deg)
    W_dep = _blade_world(lay, surf, 0.0)
    mid = {f: _blade_world(lay, surf, f) for f in (30.0, 60.0)}
    xf = x_hub_front(lay)
    for k, axp in enumerate(axs):
        axp.fill_between(xx[xx <= lay.static_tail_x], 0, hull_r(xx[xx <= lay.static_tail_x]),
                         color=HULL_FILL, lw=0, label="static hull")
        xr = xx[xx >= xf]
        axp.fill_between(xr, 0, hull_r(xr), color=HUB_FILL, lw=0, hatch="////", edgecolor="#d8d7d2",
                         label="hub + spinner (rotating)")
        axp.plot(xx, hull_r(xx), color=INK2, lw=1.2)
        axp.plot(xx, profile.inner_radius(xx), color=MUTED, lw=0.7, ls=(0, (4, 3)))
        axp.axvline(lay.static_tail_x, color=INK2, lw=0.8, ls="--")
        m = ~np.isnan(tb["floor_r"])
        axp.plot(tb["x"][m], tb["floor_r"][m], color=INK, lw=1.2, ls=(0, (5, 2)), label="pocket floor (min r)")
        axp.plot(tb["x"][m], tb["keepout_r"][m], color=INK2, lw=0.9, ls=":", label="keep-out (floor - MIN_WALL)")
        axp.fill_between([60, slot["x_mm"].max()], 0, P.SHAFT_D / 2, color="#b9b8b2", lw=0, label="8 mm shaft")
        cell = 0.4 if k == 0 else 0.2
        for f in (30.0, 60.0):
            Wf = mid[f]
            _raster_outline(axp, Wf[:, 0], np.hypot(Wf[:, 1], Wf[:, 2]), cell, C_MID,
                            "fold 30 / 60 deg" if f == 30.0 else None, fill=False)
        _raster_outline(axp, W_dep[:, 0], np.hypot(W_dep[:, 1], W_dep[:, 2]), cell, C_DEP, "deployed (fold 0)")
        _raster_outline(axp, W_st[:, 0], np.hypot(W_st[:, 1], W_st[:, 2]), cell, C_STOW,
                        f"stowed (fold {lay.fold_deg:g} deg)")
        axp.plot(lay.hinge_x, lay.hinge_r, "o", ms=5, color=INK, zorder=5)
        axp.set_aspect("equal")
        axp.set_xlabel("x (mm), flow ->")
        axp.set_ylabel("radius from body axis (mm)")
        axp.grid(color="#e9e8e4", lw=0.5)
        axp.spines[["top", "right"]].set_visible(False)
        if k == 0:
            axp.set_xlim(0, 172)
            axp.set_ylim(0, 102)
            axp.set_title("Ovo rotor: fold kinematics, side view (meridional plane of one blade)",
                          loc="left", fontsize=11, color=INK)
            axp.annotate(f"hinge ({lay.hinge_x:g}, {lay.hinge_r:g})", (lay.hinge_x, lay.hinge_r),
                         (lay.hinge_x + 8, lay.hinge_r + 30), color=INK, fontsize=9,
                         arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
            axp.annotate(f"static tail end x = {lay.static_tail_x:g}", (lay.static_tail_x, 26),
                         (lay.static_tail_x - 60, 55), color=INK2, fontsize=9,
                         arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
            axp.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)
        else:
            axp.set_xlim(98, 142)
            axp.set_ylim(0, 33)
            axp.set_title("Hinge zone (zoom)", loc="left", fontsize=11, color=INK)
            axp.add_patch(Circle((lay.hinge_x, lay.hinge_r), lay.knuckle_r, fill=False, ec=INK, lw=0.9))
            axp.add_patch(Circle((lay.hinge_x, lay.hinge_r), lay.pin_d / 2, color=INK))
            ms = slot["x_mm"] <= lay.hinge_x + 14
            axp.plot(slot["x_mm"][ms], slot["slot_floor_re_mm"][ms], color=C_DEP, lw=1.0, ls="-.",
                     label="hub slot floor (swept root + 0.3)")
            xs_ = np.linspace(xf, slot["x_mm"].max(), 100)
            axp.plot(xs_, hull_r(xs_) - lay.spinner_wall, color=INK2, lw=0.7, ls=(0, (2, 2)),
                     label="spinner inner surface")
            # stop block and magnets
            xs0 = lay.hinge_x - lay.heel_t / 2
            lam0, lam1 = heel_geometry(lay)["flat_face_from_pin_mm"]
            axp.add_patch(Rectangle((xf, lay.hinge_r - lam1), xs0 - xf, lam1 - lam0,
                                    fill=False, ec=INK, lw=0.8, ls="--", zorder=6))
            axp.text(xf + 0.3, lay.hinge_r - lay.heel_reach - 1.6, "stop block", fontsize=8, color=INK)
            dx_, dr_ = det["deployed"]["hub_magnet_hinge_plane_x_re_tau_mm"][:2]
            H = P.MAGNET["H"]
            Dm = P.MAGNET["D"]
            mk = dict(color="#4a3aa7", zorder=7)
            axp.add_patch(Rectangle((dx_ - H / 2, dr_ - Dm / 2), H, Dm, label="N52 4x2 magnets", **mk))
            axp.add_patch(Rectangle((dx_ + H + 2 * lay.magnet_recess - H / 2, dr_ - Dm / 2), H, Dm, **mk))
            st = det["stowed"]
            axp.add_patch(Rectangle((st["x_mm"] - Dm / 2, st["static_magnet_centre_r_mm"] - H / 2), Dm, H, **mk))
            axp.add_patch(Rectangle((st["x_mm"] - Dm / 2, st["blade_inner_r_mm"] + lay.magnet_recess), Dm, H, **mk))
            axp.text(st["x_mm"] - 2.0, st["static_magnet_centre_r_mm"] - 3.2, "stowed\ndetent", fontsize=8,
                     color=INK, zorder=8)
            axp.text(dx_ - 1.0, dr_ + Dm / 2 + 0.6, "deployed\ndetent", fontsize=8, color=INK, zorder=8)
            axp.legend(loc="lower right", fontsize=8, frameon=True, framealpha=0.9, edgecolor="none")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_end(lay: Layout, surf: BladeSurface, path: Path) -> None:
    """View from aft (looking upstream, +Y to the right, Z up): three stowed blades on
    the tail (left) and three deployed blades (right)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(12, 6.3))
    t = np.linspace(0, 2 * math.pi, 361)
    Rmax = P.HULL_D / 2
    for k, axp in enumerate(axs):
        fold = lay.fold_deg if k == 0 else 0.0
        axp.fill(Rmax * np.cos(t), Rmax * np.sin(t), color=HULL_FILL, lw=0)
        axp.plot(Rmax * np.cos(t), Rmax * np.sin(t), color=INK2, lw=1.0, label=f"egg outline (D {P.HULL_D:g})")
        rs = float(hull_r(lay.static_tail_x))
        axp.plot(rs * np.cos(t), rs * np.sin(t), color=INK2, lw=0.7, ls="--",
                 label=f"static tail end (x {lay.static_tail_x:g})")
        rh = float(hull_r(lay.hinge_x))
        axp.fill(rh * np.cos(t), rh * np.sin(t), color=HUB_FILL, ec=INK2, lw=0.7, hatch="////",
                 label=f"hub + spinner at the hinge plane (x {lay.hinge_x:g})")
        axp.fill(P.SHAFT_D / 2 * np.cos(t), P.SHAFT_D / 2 * np.sin(t), color="#b9b8b2")
        for i in range(lay.n_blades):
            psi = lay.psi(i)
            W = _blade_world(lay, surf, fold, i)
            _raster_outline(axp, W[:, 1], W[:, 2], 0.25 if k == 0 else 0.5, C_STOW if k == 0 else C_DEP,
                            ("stowed blade" if k == 0 else "deployed blade") if i == 0 else None)
            _, e_r, e_t = rotor_axes(psi)
            Hh = lay.hinge_r * e_r
            pa = Hh - lay.pin_l / 2 * e_t
            pb = Hh + lay.pin_l / 2 * e_t
            axp.plot([pa[1], pb[1]], [pa[2], pb[2]], color=INK, lw=2.0,
                     label="hinge pin (tangential)" if i == 0 else None)
        if k == 1:
            rt = float(load_blade().r_tip)
            axp.plot(rt * np.cos(t), rt * np.sin(t), color=MUTED, lw=0.7, ls=":")
            axp.annotate("", xy=(-rt * 0.80 * math.sin(math.radians(35)), rt * 0.80 * math.cos(math.radians(35))),
                         xytext=(-rt * 0.80 * math.sin(math.radians(15)), rt * 0.80 * math.cos(math.radians(15))),
                         arrowprops=dict(arrowstyle="->", color=INK2, lw=1.2,
                                         connectionstyle="arc3,rad=0.15"))
            axp.text(-68, rt * 0.80, "rotation", fontsize=9, color=INK2)
            lim = rt + 6
            axp.set_title(f"Deployed (fold 0): rotor D {2 * rt:.0f} mm", loc="left", fontsize=11, color=INK)
        else:
            lim = Rmax + 6
            axp.set_title(f"Stowed (fold {lay.fold_deg:g} deg): blades inside the egg outline",
                          loc="left", fontsize=11, color=INK)
        axp.set_xlim(-lim, lim)
        axp.set_ylim(-lim, lim)
        axp.set_aspect("equal")
        axp.set_xlabel("y (mm)")
        axp.set_ylabel("z (mm)")
        axp.grid(color="#e9e8e4", lw=0.5)
        axp.spines[["top", "right"]].set_visible(False)
    hs, ls = {}, []
    for axp in axs:
        for h_, l_ in zip(*axp.get_legend_handles_labels()):
            if l_ not in hs:
                hs[l_] = h_
                ls.append(l_)
    fig.suptitle("Ovo rotor seen from aft (looking upstream)", x=0.01, y=0.985, ha="left",
                 fontsize=12, color=INK)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.88, bottom=0.2, wspace=0.2)
    fig.legend([hs[l_] for l_ in ls], ls, loc="lower center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 0.005))
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ==========================================================================
# JSON output and main
# ==========================================================================

X_AERO_ZONE = 100.0     # pocket summary: aero span zone (forward of the root/notch transition)

REC_KEYS = dict(HX="hinge_x_mm", HR="hinge_r_mm", front="hinge_x_minus_static_tail_x_mm",
                cone="shank_cone_deg", fold="fold_deg", XST="static_tail_x_mm", heel_t="heel_t_mm",
                heel_reach="heel_reach_mm", pin_cover_hull="pin_cover_to_hull_mm",
                pin_cover_core="pin_cover_in_hub_core_mm", depth_aero_max="pocket_depth_aero_max_mm",
                keep_tail_min="keepout_r_last_12mm_min_mm", notch_floor="notch_floor_r_mm",
                slot_floor="slot_floor_r_mm", hub_wall_over_bore="hub_wall_over_bore_mm",
                heel_flat="heel_flat_face_mm", stop_block_max="stop_block_thickness_mm",
                stop_lever="stop_lever_mm", clear_aero="stowed_clear_aero_mm",
                clear_shank="stowed_clear_root_mm", ear_clear="ear_clear_mm",
                heel_stow="heel_stowed_inside_spinner_mm", heel_dep="heel_deployed_inside_spinner_mm",
                nb_root="neighbour_root_min_mm", x_tip="stowed_tip_x_mm", feasible="feasible")


def _clean(o, nd: int = 5):
    """JSON-safe: numpy -> lists, floats to nd significant digits, nan/inf -> None."""
    if isinstance(o, dict):
        return {str(k): _clean(v, nd) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v, nd) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist(), nd)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return None if not math.isfinite(v) else float(f"{v:.{nd}g}")
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def margin_rank(rec: dict) -> tuple:
    """Selection key among feasible points: the smallest margin over the binding
    constraints (pin cover in the hub core, stowed heel, magnet room on the heel),
    to 0.01 mm; ties go to the smallest hinge-to-static-tail offset (longest static
    tail, shortest slotted hub ring)."""
    m = min(rec["pin_cover_core"] - FEASIBLE["pin_cover_core"], rec["heel_stow"] - FEASIBLE["heel_stow"],
            rec["heel_flat"] - FEASIBLE["heel_flat"])
    return (round(m, 2), -rec["front"])


def run(write: bool = True, figures: bool = True) -> dict:
    """Regenerate data/layout_rotor.json and the fold figures for the chosen layout."""
    import time
    t0 = time.time()
    blade = load_blade()
    lay = chosen_layout()
    ax = solve_axis(blade, lay, n=41, iters=3)
    ev = evaluate(blade, lay, ax=ax, d_fold=1.0)
    S, pm, tb = ev["surf"], ev["pm"], ev["table"]
    folds = np.arange(0.0, lay.fold_deg + 1e-9, 1.0)
    slot = hub_slot_table(lay, S, folds)
    nb = neighbour_check(lay, S, n=10)
    mass = blade_mass_props(lay, S)
    cal = magnet_calibration()
    det = detents(lay, ax, S, pm, tb, mass, cal, ev["stop"])
    loads = bem_loads(lay=lay)
    hub = hub_concept(lay, S, ev, slot, loads)
    print(f"core done {time.time() - t0:.0f} s", flush=True)
    fstudy = fold_study(blade, lay)
    tw = twist_effect(blade, lay, ax)
    secs = stowed_sections(blade, lay, ax)
    Tst = blade_to_world(lay.fold_deg, lay.psi0_deg, lay)
    proud = {f"{k}_mm": float(height_above_hull(apply(Tst, G)).max())
             for k, G in (("aero", densify(S.aero, 0.3).reshape(-1, 3)),
                          ("shank", densify(S.shank, 0.3).reshape(-1, 3)),
                          ("knuckle", S.knuckle.reshape(-1, 3)), ("heel", S.heel))}
    print(f"studies done {time.time() - t0:.0f} s", flush=True)
    # sweep record (out/work, regenerated with --sweep)
    sweep_rows = []
    sp = WORK_DIR / "hinge_sweep.jsonl"
    if sp.exists():
        sweep_rows = [json.loads(l) for l in sp.read_text().splitlines()]
    for r in sweep_rows:
        r["feasible"] = feasible(r)
    feas = [r for r in sweep_rows if r["feasible"]]
    pick = max(feas, key=margin_rank) if feas else None
    m = ~np.isnan(tb["depth"])
    xm = tb["x"][m]
    it = ev["interference"]
    it_slot = it.pop("hub_slot_halfwidth_needed_mm", None)  # superseded by hub_slot
    cs = cone_stats(ax)
    aero = xm <= X_AERO_ZONE
    out = {
        "schema": "ovo.layout_rotor/1",
        "generated_by": "python -m ovo.layout",
        "blade_source": blade.source,
        "conventions": {
            "frame": "X = body axis, flow +X, nose at x = 0, Z up, right handed; mm and deg",
            "azimuth": "psi about +X from +Z, right handed: e_r = (0, -sin psi, cos psi), "
                       "e_t = (0, -cos psi, -sin psi); blade i at psi = 120 i; rotor turns +psi "
                       "(clockwise seen from upstream, as data/blade_design.json)",
            "theta_internal_layout": "layout_internal uses theta from +Z toward +Y, i.e. theta = -psi",
            "blade_frame": "origin on the hinge axis at H = (hinge_x, 0, 0) + hinge_r e_r; e_s span "
                           "(deployed = e_r), e_c chord along the pin (= e_t, toward the LE), "
                           "e_n = e_s x e_c thickness (deployed = +X, suction side)",
            "fold": "rotation by +fold about e_t: e_s = cos f e_r - sin f e_x, e_n = sin f e_r + cos f e_x; "
                    "0 = deployed, fold_deg = stowed (tip forward over the static tail)",
            "transform": "blade_to_world(fold, psi): 4x4 with columns (e_s, e_t, e_n, H)",
            "pitch_axis": "pitch point P(r) = (s, lead, aft) in the blade frame at PITCH_AXIS chord; "
                          "r = sqrt((hinge_r + s)^2 + lead^2) is the deployed perpendicular radius",
            "twist": "as bem.py: chord vs plane of rotation in the section plane normal to the local "
                     "pitch axis, + = LE upstream",
            "pocket_azimuth": "az relative to the blade azimuth psi, + toward e_t (direction of rotation)",
        },
        "chosen": {
            "static_tail_x_mm": lay.static_tail_x, "hinge_x_mm": lay.hinge_x, "hinge_r_mm": lay.hinge_r,
            "fold_deg": lay.fold_deg, "hub_gap_mm": lay.hub_gap, "hub_front_x_mm": x_hub_front(lay),
            "pin_d_mm": lay.pin_d, "pin_l_mm": lay.pin_l, "knuckle_r_mm": lay.knuckle_r,
            "knuckle_l_mm": lay.knuckle_l, "ear_t_mm": lay.ear_t, "ear_clear_mm": lay.ear_clear,
            "heel_reach_mm": lay.heel_reach, "heel_t_mm": lay.heel_t, "heel_w_mm": lay.heel_w,
            "heel_angle_deg": lay.heel_angle_deg, "shank_cone_deg": lay.shank_cone_deg,
            "shank_blend_start_mm": lay.shank_blend_start, "spinner_wall_mm": lay.spinner_wall,
            "spinner_gap_mm": lay.spinner_gap, "run_clear_mm": lay.run_clear,
            "pocket_clear_mm": lay.pocket_clear, "root_clear_mm": lay.root_clear,
            "side_clear_mm": lay.side_clear, "proud_mm": lay.proud, "magnet_recess_mm": lay.magnet_recess,
            "shaft_end_x_mm": lay.hinge_x + lay.shaft_end_dx, "rotor_r_tip_mm": float(blade.r_tip),
            "r_aero_root_mm": float(blade.r_root),
        },
        "selection": {
            "rule": "feasible: pin cover in hub core >= 2, stowed heel >= 0.5 inside the spinner, heel face "
                    ">= magnet D + 0.6, stop block >= magnet H + recess + MIN_WALL + 0.6, ear clear >= 0.5, "
                    "neighbour >= 1, swept root >= ROOT_CLEAR from the static tail, hub wall >= MIN_WALL; "
                    "then max of the smallest margin (pin cover, heel, magnet room), ties to the smallest "
                    "hinge-to-static-tail offset (front 10 only adds keep-out the internal layout does not need)",
            "constraints": FEASIBLE,
            "chosen_matches_rule": None,
            "rule_pick": None if pick is None else {REC_KEYS[k]: v for k, v in pick.items() if k in REC_KEYS},
            "hinge_sweep": [{REC_KEYS[k]: v for k, v in r.items() if k in REC_KEYS} for r in sweep_rows],
            "sweep_source": "out/work/layout_rotor/hinge_sweep.jsonl (python -m ovo.layout --sweep), "
                            "fold 90, d_fold 3 deg",
        },
        "axis": {
            "r_mm": ax.r, "s_mm": ax.s, "lead_mm": ax.lead, "aft_mm": ax.aft,
            "slope_daft_ds": ax.slope, "cone_local_deg": ax.cone_deg(),
            "lead_station_optimum_mm": ax.lead_station_opt,
            "aft_poly_coef_mm": ax.aft_coef, "lead_poly_coef_mm": ax.lead_coef,
            "poly_variable": "t = (r - r_root) / (r_tip - r_root), numpy.polyval order (highest first)",
            "derivation": "for each station the aft offset is solved (bisection) so the highest point of the "
                          "stowed section is flush with the hull (proud 0); lead is the per-station "
                          "optimum (golden section) that minimises the flush pocket depth, then fitted by "
                          "polynomials and iterated with the section planes (they depend on d aft/ds)",
            **cs,
        },
        "fold_study": fstudy,
        "fold_choice": {"fold_deg": lay.fold_deg,
                        "reason": "stowed depth and keep-out barely change with fold (the axis re-solves); "
                                  "the fold angle mostly sets the deployed cone (+1 deg cone per +1 deg fold). "
                                  "90 deg gives a mean cone close to the BEM nominal 12 deg with the tip "
                                  "bending forward; see fold_study"},
        "stowed_sections": secs,
        "stowed_max_over_hull": {**proud, "note": "highest point of each stowed part above the hull "
                                 "envelope (normal); design value 0 (flush), residual = solver tolerance"},
        "twist_effect": tw,
        "pocket_depth": {
            "x_mm": tb["x"][m], "depth_mm": tb["depth"][m], "depth_normal_mm": tb["depth_normal"][m],
            "hull_r_mm": tb["hull_r"][m], "floor_r_mm": tb["floor_r"][m], "keepout_r_mm": tb["keepout_r"][m],
            "half_width_deg": tb["half_w_deg"][m], "half_width_mm": tb["half_w_mm"][m],
            "az_center_deg": tb["az_c"][m], "az_lo_deg": tb["az_lo"][m], "az_hi_deg": tb["az_hi"][m],
            "land_between_pockets_mm": tb["land_mm"][m],
            "note": "radial depth below the hull of the pocket floor (blade inner surface - clearance "
                    "0.3 under the aero span, 0.8 under the root, side clearance 0.8); keepout_r = floor "
                    "- MIN_WALL along the hull normal = the largest radius internal parts may use under a "
                    "pocket. Half-width is about az_center (relative to the blade azimuth).",
        },
        "pocket_summary": {
            "x_range_mm": [float(xm.min()), lay.static_tail_x],
            "aero_zone_x_max_mm": X_AERO_ZONE,
            "depth_aero_max_mm": float(np.nanmax(tb["depth"][m][aero])),
            "depth_aero_max_x_mm": float(xm[aero][int(np.nanargmax(tb["depth"][m][aero]))]),
            "depth_root_zone_max_mm": float(np.nanmax(tb["depth"][m][~aero])),
            "depth_max_mm": float(np.nanmax(tb["depth"][m])),
            "depth_max_x_mm": float(xm[int(np.nanargmax(tb["depth"][m]))]),
            "keepout_min_mm": float(np.nanmin(tb["keepout_r"][m])),
            "keepout_min_x_mm": float(xm[int(np.nanargmin(tb["keepout_r"][m]))]),
            "notch_at_static_tail_end": {"floor_r_mm": ev["rec"]["notch_floor"],
                                         "keepout_last_12mm_min_mm": ev["rec"]["keep_tail_min"]},
            "half_width_max_mm": float(np.nanmax(tb["half_w_mm"][m])),
            "land_min_mm": float(np.nanmin(tb["land_mm"][m])),
            "crosses_equator_x_mm": P.EQUATOR_X,
            "floor_print_rule": "both shells print +x down; the floor may rise at most tan(%g deg) per mm "
                                "going forward (printable_floor), so the root-notch ramp and the pocket's "
                                "forward end are cut as <= %g deg ramps, no supports" % (lay.floor_overhang_deg,
                                                                                          lay.floor_overhang_deg),
            "floor_overhang_max_deg": floor_overhang_max_deg(pm)[0],
            "floor_overhang_max_x_mm": floor_overhang_max_deg(pm)[1],
        },
        "pocket_width_mm": float(2 * np.nanmax(tb["half_w_mm"][m])),
        "pocket_theta_deg": [0.0, 120.0, 240.0],
        "shaft_aft_x_mm": lay.hinge_x + lay.shaft_end_dx,
        "hub_slot": {**slot, "note": "hinge-plane coords of blade 0 (r_e along e_r, tau along the pin), "
                                     "swept root over fold 0..stowed in 1 deg steps, + run clearance; "
                                     "cut the hub and the spinner to this per blade"},
        "hinge": ev["hinge"],
        "hub": hub,
        "detents": det,
        "blade_mass": mass,
        "interference": {**it, **nb, "pass": None},
        "proposed_params": {},
    }
    ok = (it["static_tail_clear_min_mm"] >= 0 and it["ear_clear_min_mm"] > 0
          and nb["neighbour_root_min_dist_mm"] >= 1.0 and it["deployed_neighbour_min_dist_mm"] > 5
          and it["heel_stowed_below_spinner_inner_min_mm"] > 0 and it["heel_deployed_below_spinner_inner_min_mm"] > 0
          and it["deployed_aero_min_height_over_hull_mm"] > 0)
    out["interference"]["pass"] = bool(ok)
    out["interference"]["criteria"] = (
        "static tail: swept blade never below the pocket floor (floor already includes 0.3/0.8 clearance); "
        "ears: no blade point in the pin + 2 mm zone; neighbours >= 1 mm at any fold pair; deployed blades "
        "apart; heel inside the spinner at both ends; deployed aero span above the hull")
    if pick is not None:
        out["selection"]["chosen_matches_rule"] = bool(
            pick["HX"] == lay.hinge_x and pick["HR"] == lay.hinge_r and pick["front"] == lay.hinge_x - lay.static_tail_x)
    out["proposed_params"] = {
        "STATIC_TAIL_X": lay.static_tail_x, "HINGE_X": lay.hinge_x, "HINGE_R": lay.hinge_r,
        "ROTOR_R_TIP_TARGET": float(blade.r_tip), "PIN_L": lay.pin_l, "FOLD_DEG": lay.fold_deg,
        "KNUCKLE_D": 2 * lay.knuckle_r, "KNUCKLE_L": lay.knuckle_l, "HUB_EAR_T": lay.ear_t,
        "SHANK_CONE_DEG": lay.shank_cone_deg, "HEEL_REACH": lay.heel_reach, "HEEL_T": lay.heel_t,
        "HEEL_W": lay.heel_w, "SHAFT_AFT_X": lay.hinge_x + lay.shaft_end_dx,
        "POCKET_CLEAR": lay.pocket_clear, "ROOT_CLEAR": lay.root_clear,
    }
    if loads.get("curved_axis_check_used"):
        bem_issue = ("BEM layout_axis_check flew the final blade on hinge (%g, %g) and the solved curved axis: "
                     "cp peak %.4f vs %.4f nominal, max flap moment %.3f vs %.3f Nm at %g m/s (+%.1f %%), used "
                     "for the hub load cases. blade_design.json rotor block still lists hinge (%g, %g) until "
                     "params.py adopts the proposed values; the in-plane lead is still ignored there."
                     % (lay.hinge_x, lay.hinge_r, loads["cp_peak"], loads["cp_peak_nominal_axis"],
                        loads["max_flap_moment_Nm"], loads["max_flap_moment_nominal_axis_Nm"],
                        loads["max_flow_m_s"],
                        100 * (loads["max_flap_moment_Nm"] / loads["max_flap_moment_nominal_axis_Nm"] - 1),
                        loads["bem_hinge_x_mm"], loads["bem_hinge_r_mm"]))
    else:
        bem_issue = ("BEM loads are for hinge (%s, %s) and a straight axis; re-run bem.py's layout axis check "
                     "with hinge (%g, %g) and the solved curved axis (axis block)."
                     % (loads.get("bem_hinge_x_mm"), loads.get("bem_hinge_r_mm"), lay.hinge_x, lay.hinge_r))
    out["open_issues"] = [
        bem_issue,
        "In-plane lead (pitch axis offset toward the LE, 4.9 to 11 mm) yaws the sections vs the local "
        "tangential flow by atan(lead/r) (stowed_sections.yaw_vs_local_tangent_deg); BEM assumes 0.",
        "Pockets cross the equator joint (x 57.5): the internal layout must thicken the shell there.",
        "Magnet force model calibrated to catalogue pull (k ~ 0.63); verify detent feel on a coupon.",
        "Pin length 22 mm: check stock (ISO 8734 3 x 22); fallback 3 x 20 with ear_t 2.7.",
    ]
    if write:
        DATA_JSON.parent.mkdir(parents=True, exist_ok=True)
        DATA_JSON.write_text(json.dumps(_clean(out), indent=1) + "\n")
        print(f"wrote {DATA_JSON.relative_to(ROOT)}", flush=True)
    if figures:
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        figure_side(lay, S, tb, slot, det, FIG_DIR / "layout_fold_side.png")
        figure_end(lay, S, FIG_DIR / "layout_fold_end.png")
        print("wrote docs/figures/layout_fold_side.png, layout_fold_end.png", flush=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)
    return out


SWEEP_GRID = ([(hx, hr, 8.0) for hx in (116.0, 117.0, 118.0) for hr in (13.5, 14.0, 14.5)]
              + [(116.0, 14.0, 6.0), (116.0, 14.0, 10.0)])


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", action="store_true", help="run/resume the hinge sweep first (~15 s/point)")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()
    if a.sweep:
        hinge_sweep(SWEEP_GRID)
    run(figures=not a.no_figures)
