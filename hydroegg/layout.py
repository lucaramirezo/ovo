"""Rotor-side layout of HydroEgg: fold kinematics, stowage pockets, hub and hinge.

Run ``python -m hydroegg.layout`` to regenerate ``data/layout_rotor.json`` and
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
    spinner_wall: float = 1.2       # spinner shell over the hub
    # deployed stop: heel tab under the knuckle
    heel_reach: float = 7.0         # pin axis -> heel tip
    heel_t: float = 4.0             # heel thickness (in the s-n plane)
    heel_angle_deg: float = 0.0     # deployed heel direction from straight down, + aft
    # shank
    shank_cone_deg: float = 25.0    # local cone of the shank spine at the knuckle
    shank_blend_start: float = 4.5  # s where the knuckle section starts to morph (> ear zone)
    # stowage
    proud: float = 0.0              # allowed protrusion of the stowed blade over the hull
    pocket_clear: float = 0.3       # pocket floor clearance under the aero span
    root_clear: float = 0.8         # min clearance stowed root <-> static tail
    side_clear: float = 0.8         # pocket side-wall clearance

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
    box from the pin axis to heel_reach along heel_dir, heel_t thick, knuckle_l long."""
    d, q = heel_dir(lay)
    lam = np.linspace(0.0, lay.heel_reach, n_l)
    w = np.linspace(-lay.heel_t / 2, lay.heel_t / 2, n_w)
    tau = np.linspace(-lay.knuckle_l / 2, lay.knuckle_l / 2, n_t)
    LL, WW, TT = np.meshgrid(lam, w, tau, indexing="ij")
    on = ((np.abs(WW) >= lay.heel_t / 2 - 1e-9) | (LL >= lay.heel_reach - 1e-9)
          | (np.abs(TT) >= lay.knuckle_l / 2 - 1e-9))
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
    return PocketMap(xs, azs, fl, bf, bt)


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
                 d_fold: float = 1.0) -> dict:
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
    fgrid = np.linspace(0.0, lay.fold_deg, 7)
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
