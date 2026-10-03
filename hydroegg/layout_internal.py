"""Internal layout study for HydroEgg: everything inside the static hull.

Run ``python -m hydroegg.layout_internal`` to regenerate ``data/layout_internal.json``
and ``docs/figures/layout_section*.png``. A colleague owns the rotor side (fold
kinematics, blade pockets, hub); this module owns the inside: drivetrain, generator
bay, equator joint, nose (tow eye + cable gland), roll fins, mass and buoyancy,
assembly sequence and part split.

Frame: X = body axis, flow travels +X, nose at x = 0, rotor at the tail. Z up.
theta = angle about +X measured from +Z toward +Y. Units: mm / deg in geometry,
SI in physics (function arguments say which).

Architecture chosen here (justified in the JSON ``decisions`` block):
  * Front shell F (nose -> equator): tow eye, M12 cable gland, motor carrier
    (twist-lock), generator with its bell facing aft, coupler half A. All
    electrics stay in F, no wire crosses the joint.
  * Aft shell A (equator -> static tail end): spigot with O-ring + 3 bayonet lugs,
    bearing tube with a 608 (front, locating, takes thrust) and a 688 (rear,
    floating), lip seal at the tail end, 8 mm shaft, coupler half B (also the
    thrust collar).
  * The coupler is a two-dog claw that engages blind when the shells close.
Inputs it reads when present: data/layout_rotor.json (pocket depth), data/
hardware.json (catalogue data), data/blade_design.json (cone, tip radius).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import params as P
from . import profile as prof

ROOT = Path(__file__).resolve().parents[1]
DATA_JSON = ROOT / "data" / "layout_internal.json"
ROTOR_JSON = ROOT / "data" / "layout_rotor.json"
HW_JSON = ROOT / "data" / "hardware.json"
BLADE_JSON = ROOT / "data" / "blade_design.json"
FIG_DIR = ROOT / "docs" / "figures"
WORK = ROOT / "out" / "work" / "layout_internal"
SOURCES_JSON = WORK / "sources" / "sources.json"

RHO = P.WATER_RHO
NU = P.WATER_NU
G = 9.81

# ==========================================================================
# Study-owned constants (choices, not catalogue data)
# ==========================================================================
T_MIN = 2.0                 # mm, structural wall of rings / pocket floors (5 perimeters)
POCKET_FLOOR_WALL = 2.0     # mm, shell left under a blade pocket floor
TUBE_WALL = 2.4             # mm, bearing tube wall (6 perimeters)
C_GEN = 0.5                 # mm, radial clearance generator envelope <-> anything
C_SEAL = P.CLEAR_SLIDE      # mm, radial gap spigot <-> seal bore (O-ring extrusion gap)

# Fallback pocket model (task brief) until data/layout_rotor.json exists
POCKET_FALLBACK = dict(x0=45.0, x1=118.0, depth=6.0, width=30.0, theta=(0.0, 120.0, 240.0))

# Equator joint (bayonet + radial O-ring gland on the aft-shell spigot)
OR_CS = P.ORING_EQUATOR["cs"]     # 2.0 mm
OR_SQUEEZE = 0.225                # nominal; +-0.15 mm radial print tolerance keeps 15-30 %
OR_GROOVE_W = 2.4                 # mm, groove bottom width; forward flank at 45 deg
OR_LAND_ROOT = 1.5                # mm, land between equator face and groove
OR_LAND_TIP = 1.5                 # mm, land forward of the groove
TOL_PRINT_R = 0.15                # mm, radial print tolerance used in the gland check
LUG_N = 3
LUG_H = 1.5                       # mm, radial height
LUG_T = 3.0                       # mm, axial thickness (1.5 straight + 1.5 at 45 deg)
LUG_W_DEG = 20.0                  # angular width
GAP_W_DEG = 24.0                  # entry gap in the socket lip
ENGAGE_DEG = 35.0                 # twist from entry to lock
LIP_T = 3.0                       # mm, socket lip axial thickness at its ID
CH_EXTRA = 0.4                    # mm, axial running clearance of the lug in its channel

# Drivetrain stations
SEAL_RECESS = 0.5                 # mm, lip seal outer face behind the static-tail end face
SEAL_SHOULDER_T = 1.5             # mm, ring between seal bore and rear bearing seat
SEAL_SHOULDER_ID = 12.5           # mm, supports 688 outer ring and seal back face
FRONT_SHOULDER_ID = 18.0          # mm, supports 608 outer ring, lets a 688 pass through
COUPLER_STACK = 16.0              # mm, bell face -> 608 front face (claw A 5, float 1, claw B 10)
COUPLER_OD = 22.0
SHAFT_AFT_FALLBACK = 146.0        # mm, shaft end inside hub/spinner until layout_rotor.json says
CARRIER_T = 4.0                   # mm, motor carrier plate
CARRIER_RING = 1.6                # mm, inward ring that holds the carrier tab slots
MOTOR_PCD_MIN = 25.0              # mm, smallest mount PCD the carrier must carry (params)

# Rotor loads (task brief)
CT_TOW = 0.8
V_TOW = 2.5                       # m/s
SF_TOW = 3.0
CP_DESIGN = 0.35
CD_HULL = 0.10                    # ESTIMATE, Myring body on frontal area (no source run)
RADIAL_LOAD_FRAC = 0.10           # ESTIMATE, radial rotor load as a fraction of thrust
THRUST_ASYM_FRAC = 0.20           # ESTIMATE, one blade off by 20 % of its share of thrust

# Tow eye (printed, part of the nose)
EYE_HOLE_D = 8.0
EYE_W = 5.0                       # mm, ring width in its plane
EYE_T = 6.0                       # mm, ring thickness
EYE_STANDOFF = 1.0                # mm, ring outer edge ahead of the nose tip... see nose()

# Cable gland (M12x1.5, threaded boss)
GLAND_X = 9.0                     # mm, boss centre on the surface, bottom (theta = 180)
GLAND_BOSS_D = 20.0
GLAND_ENGAGE = 9.0                # thread length of the gland (source: sunrom listing)

# Fins
FIN_SECTION = "clarky"            # flat lower surface: prints flat with no support
FIN_N_CRIT = 5.0                  # disturbed river inflow, same as the BEM study
FIN_CP_SPAN = 0.44                # ESTIMATE, spanwise centre of lift / span (rectangular)
FIN_E_OSWALD = 0.8                # ESTIMATE
FIN_CLMAX_3D_FACTOR = 0.9         # ESTIMATE, 3D CLmax / 2D CLmax
FIN_CL_RESERVE = 0.75             # design CL <= this x CLmax_3D
FIN_ROOT_GAP = 1.0                # mm, plate root above the max hull radius
FIN_TIP_CLEAR = 5.0               # mm, margin to the blade-tip sweep sphere
RAIL_H = 4.0                      # mm, dovetail rail height (raised) or slot depth (flush)
RAIL_X = (20.0, 43.0)             # mm, rail on the front shell (+-Y)

# Mass model
RHO_PLA = 1.24                    # g/cc, Prusament TDS
RHO_PETG = 1.27                   # g/cc, Prusament TDS
SKIN = 1.2                        # mm, 3 walls ~ 3 x 0.4 (also used for top/bottom shells)
INFILL = 0.15                     # gyroid
RHO_STEEL = 7.9                   # g/cc, ESTIMATE stainless 304 typical
RHO_LEAD = 11.34                  # g/cc, ESTIMATE
RHO_NBR = 1.3                     # g/cc, ESTIMATE
RHO_NDFEB = 7.5                   # g/cc, ESTIMATE

# ==========================================================================
# Inputs
# ==========================================================================


def _load(path: Path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def sources() -> dict:
    return _load(SOURCES_JSON) or {}


def _find(d, keys):
    """Return the first value found under any of the dotted key paths."""
    for k in keys:
        cur = d
        ok = True
        for part in k.split("."):
            if isinstance(cur, dict) and part in cur:
                cur = cur[part]
            else:
                ok = False
                break
        if ok:
            return cur
    return None


def rotor_inputs() -> dict:
    """Pocket model and rotor geometry. Uses data/layout_rotor.json when present."""
    out = dict(source="fallback (task brief: 6 mm deep, 30 mm wide, x 45-118, 3 x 120 deg)",
               **POCKET_FALLBACK, depth_x=None, depth_v=None,
               shaft_aft_x=SHAFT_AFT_FALLBACK, r_tip=P.ROTOR_R_TIP_TARGET, cone_deg=12.0)
    bd = _load(BLADE_JSON)
    if bd:
        rt = _find(bd, ["geometry.r_tip_mm", "rotor.r_tip_mm", "r_tip_mm", "proposed_params.ROTOR_R_TIP_TARGET"])
        cone = _find(bd, ["geometry.cone_deg", "rotor.cone_deg", "cone_deg", "design.cone_deg"])
        if isinstance(rt, (int, float)):
            out["r_tip"] = float(rt)
        if isinstance(cone, (int, float)):
            out["cone_deg"] = float(cone)
    else:
        out["cone_note"] = "cone 12 deg = CONE_NOMINAL in hydroegg/bem.py (blade_design.json absent)"
    rj = _load(ROTOR_JSON)
    if rj:
        pd = _find(rj, ["pocket_depth", "pockets.depth", "pocket.depth_profile", "pockets.pocket_depth"])
        xs = ds = None
        if isinstance(pd, dict):
            xs = _find(pd, ["x_mm", "x"])
            ds = _find(pd, ["depth_mm", "depth"])
        elif isinstance(pd, list) and pd and isinstance(pd[0], (list, tuple)):
            xs = [p[0] for p in pd]
            ds = [p[1] for p in pd]
        if xs and ds:
            out.update(depth_x=list(map(float, xs)), depth_v=list(map(float, ds)),
                       x0=float(min(xs)), x1=float(max(xs)), depth=float(max(ds)),
                       source="data/layout_rotor.json")
        w = _find(rj, ["pocket_width_mm", "pockets.width_mm", "pocket.width_mm"])
        th = _find(rj, ["pocket_theta_deg", "pockets.theta_deg", "pocket.theta_deg"])
        sx = _find(rj, ["shaft_aft_x_mm", "hub.shaft_end_x_mm", "shaft.aft_end_x_mm"])
        if isinstance(w, (int, float)):
            out["width"] = float(w)
        if isinstance(th, list) and len(th) == 3:
            out["theta"] = tuple(map(float, th))
        if isinstance(sx, (int, float)):
            out["shaft_aft_x"] = float(sx)
    return out


def hardware_inputs() -> dict:
    """Catalogue data. data/hardware.json wins when it carries the field."""
    src = sources()
    hw = dict(
        b608=dict(P.BEARING, C_N=3300.0, C0_N=1350.0, mass_g=11.9,
                  src="vxb listing (C, C0), SKF 608-2RSH page (mass), via web search 2026-10-03"),
        b688=dict(P.BEARING_ALT, C_N=1600.0, C0_N=710.0, mass_g=3.8,
                  src="vxb 688-2RS page (C, C0); mass from a web-search summary, VERIFY"),
        seal=dict(P.LIP_SEAL, mass_g=2.0, src="mass ESTIMATE"),
        motor=dict(name="GM4108H-120T (iPower)", D=47.0, H=22.7, H_enc=32.3, mass_g=124.0,
                   pcd=P.GENERATOR["mount_pcd"],
                   src="hardware study salvage quoting iflight pages; VERIFY"),
        gland=dict(name="M12x1.5 nylon gland, cable 3-6 mm", head_D=18.0, thread_L=9.0,
                   exposed_L=20.0, mass_g=5.0, src="sunrom listing via search; mass ESTIMATE"),
        source="defaults (out/work/layout_internal/sources/sources.json)",
    )
    hj = _load(HW_JSON)
    if hj:
        hw["source"] = "data/hardware.json + defaults for missing fields"
        for key, names in (("b608", ["bearings.608", "bearing_608", "608-2RS", "608"]),
                           ("b688", ["bearings.688", "bearing_688", "688-2RS", "688"])):
            v = _find(hj, names)
            if isinstance(v, dict):
                for f_out, f_in in (("C_N", ["C_N", "C_kN"]), ("C0_N", ["C0_N", "C0_kN"]),
                                    ("mass_g", ["mass_g"])):
                    val = _find(v, f_in)
                    if isinstance(val, (int, float)):
                        hw[key][f_out] = float(val) * (1000.0 if f_in[0].endswith("kN") and val < 50 else 1.0)
        m = _find(hj, ["generator", "motor", "generator_v2"])
        if isinstance(m, dict):
            for f_out, f_in in (("D", ["D_mm", "D"]), ("H", ["L_mm", "H_mm", "L"]), ("mass_g", ["mass_g"])):
                val = _find(m, f_in)
                if isinstance(val, (int, float)):
                    hw["motor"][f_out] = float(val)
    hw["_sources_file"] = bool(src)
    return hw


# ==========================================================================
# Hull geometry
# ==========================================================================


def r_out(x):
    return prof.radius(x)


def r_in(x):
    return prof.inner_radius(x)


def slope(x):
    x = np.asarray(x, dtype=float)
    h = 1e-3
    return (prof.radius(x + h) - prof.radius(x - h)) / (2 * h)


class Pockets:
    def __init__(self, ri: dict):
        self.ri = ri
        self.x0, self.x1 = ri["x0"], ri["x1"]
        self.width = ri["width"]
        self.theta = ri["theta"]

    def depth(self, x):
        x = np.asarray(x, dtype=float)
        if self.ri.get("depth_x"):
            d = np.interp(x, self.ri["depth_x"], self.ri["depth_v"], left=0.0, right=0.0)
        else:
            d = np.where((x >= self.x0) & (x <= self.x1), self.ri["depth"], 0.0)
        return d

    def floor_inner_r(self, x):
        """Radius of the shell's inner surface under a pocket (where the pocket exists)."""
        x = np.asarray(x, dtype=float)
        d = self.depth(x)
        rin_pk = r_out(x) - d - POCKET_FLOOR_WALL * np.sqrt(1 + slope(x) ** 2)
        return np.where(d > 0, np.minimum(rin_pk, r_in(x)), r_in(x))

    def half_angle_deg(self, x):
        return np.degrees(self.width / 2 / np.maximum(r_out(x), 1e-6))


# ==========================================================================
# Equator joint
# ==========================================================================


def joint(x_eq: float, pk: Pockets) -> dict:
    """Radial stack and axial stations of the bayonet + O-ring joint.

    Spigot on the aft shell (printed tail-down, the spigot is its top); socket in the
    front shell (printed equator-down). Seal zone at the spigot root, so the bayonet
    sits on the dry side; lugs at the spigot tip pass through the seal bore.
    """
    gd = OR_CS * (1 - OR_SQUEEZE)                       # gland depth (radial)
    groove_top = OR_GROOVE_W + gd                        # 45 deg forward flank
    L_seal = OR_LAND_ROOT + groove_top + OR_LAND_TIP
    xs = np.linspace(x_eq - L_seal, x_eq, 40)
    ro = r_out(xs)
    rb_lim = np.minimum(ro - P.HULL_WALL, ro - pk.depth(xs) - POCKET_FLOOR_WALL)
    r_b = float(np.floor(rb_lim.min() * 20) / 20)       # seal bore radius, 0.05 grid
    r_gr = r_b - gd
    r_s2 = r_b - C_SEAL
    r_lug = r_b - P.CLEAR_SLIDE                          # lugs pass through the seal bore
    r_s1 = r_lug - LUG_H
    r_lip = r_s1 + P.CLEAR_SLIDE
    r_ch = r_lug + P.CLEAR_SLIDE
    r_sp_in = min(r_gr - T_MIN, r_s1 - T_MIN)

    # stations (x decreasing = forward)
    x_seal0 = x_eq - L_seal
    x_lip0 = x_seal0 - LIP_T
    x_tip = x_lip0 - LUG_T                               # spigot tip
    x_ch0 = x_tip - CH_EXTRA                             # channel forward wall (at r_ch)
    x_ring0 = x_ch0 - (r_ch - r_lip)                     # 45 deg ceiling closes at r_lip
    groove = dict(x_aft=x_eq - OR_LAND_ROOT, x_fwd_bottom=x_eq - OR_LAND_ROOT - OR_GROOVE_W,
                  x_fwd_top=x_eq - OR_LAND_ROOT - groove_top)

    # O-ring choice: integer metric ID, stretch 1-5 % (globaloring: max 5 %)
    d_gr = 2 * r_gr
    oid = math.floor(d_gr / 1.01)
    stretch = d_gr / oid - 1
    a_or = math.pi * OR_CS ** 2 / 4

    def fill(g):
        return a_or / ((OR_GROOVE_W + g / 2) * g)

    gland = dict(
        cs_mm=OR_CS, oring_id_mm=oid, oring_name=f"{oid} x {OR_CS:g} NBR 70 (metric)",
        stretch_pct=100 * stretch, groove_root_d_mm=d_gr, bore_d_mm=2 * r_b,
        gland_depth_mm=gd, squeeze_nom_pct=100 * OR_SQUEEZE,
        squeeze_range_pct=[100 * (OR_CS - (gd + TOL_PRINT_R)) / OR_CS,
                           100 * (OR_CS - (gd - TOL_PRINT_R)) / OR_CS],
        groove_bottom_w_mm=OR_GROOVE_W, groove_top_w_mm=groove_top,
        fill_nom_pct=100 * fill(gd),
        fill_range_pct=[100 * fill(gd + TOL_PRINT_R), 100 * fill(gd - TOL_PRINT_R)],
        extrusion_gap_mm=C_SEAL,
        guidance="squeeze 15-30 %, fill 60-85 % (ficientdesign), ID stretch <= 5 % (globaloring)",
    )

    # angular layout: lugs centred between pockets
    th_p = np.array(pk.theta)
    lug_lock = (th_p + 60.0) % 360.0                     # locked lug centres
    lug_entry = (lug_lock - ENGAGE_DEG) % 360.0          # gap centres in the lip
    half_p = float(pk.half_angle_deg(x_eq))
    # angular span swept by a lug (entry -> lock), relative to the pocket edges
    sweep_lo = -ENGAGE_DEG - LUG_W_DEG / 2
    sweep_hi = LUG_W_DEG / 2
    margin_lo = 60.0 + sweep_lo - half_p
    margin_hi = 60.0 - sweep_hi - half_p
    seg_w = 360.0 / LUG_N - GAP_W_DEG
    overlap = min(LUG_W_DEG, seg_w - (ENGAGE_DEG - (GAP_W_DEG - LUG_W_DEG) / 2))
    return dict(
        x_eq=x_eq, L_seal=L_seal, L_spigot=x_eq - x_tip,
        r_bore=r_b, r_groove=r_gr, r_spigot_seal=r_s2, r_lug=r_lug, r_spigot_lug=r_s1,
        r_lip=r_lip, r_channel=r_ch, r_spigot_in=r_sp_in,
        x_seal0=x_seal0, x_lip0=x_lip0, x_tip=x_tip, x_ch0=x_ch0, x_ring0=x_ring0,
        groove=groove, gland=gland,
        walls=dict(socket_at_pocket_mm=float((ro - pk.depth(xs) - r_b).min()),
                   socket_plain_mm=float((ro - r_b).min()),
                   spigot_under_groove_mm=r_gr - r_sp_in, spigot_lug_zone_mm=r_s1 - r_sp_in,
                   channel_backing_mm=float(r_out(x_tip) - r_ch)),
        bayonet=dict(n_lugs=LUG_N, lug_h_mm=LUG_H, lug_t_mm=LUG_T, lug_w_deg=LUG_W_DEG,
                     lug_w_mm=math.radians(LUG_W_DEG) * r_s1,
                     gap_w_deg=GAP_W_DEG, engage_deg=ENGAGE_DEG, lip_t_mm=LIP_T,
                     lock_theta_deg=sorted(lug_lock.tolist()),
                     entry_theta_deg=sorted(lug_entry.tolist()),
                     lip_segment_deg=seg_w, lug_cover_at_lock_deg=overlap,
                     clear_to_pocket_edge_deg=[margin_lo, margin_hi],
                     load_face="45 deg cone on both parts (printable without support)",
                     detent="2 magnet pairs (4x2 N52) in the equator faces, attract at lock; "
                            "hard stop = channel end wall; pockets align at lock"),
    )


# ==========================================================================
# Front-shell cavity and generator bay
# ==========================================================================


def front_clear_r(x, J: dict, rail_flush: bool) -> np.ndarray:
    """Clear radius of the front-shell cavity (min over angle) for x < spigot tip."""
    x = np.asarray(x, dtype=float)
    r = np.array(r_in(x), dtype=float)
    # socket ring: 45 deg chamfer from r_lip at x_ring0 out to r_in going forward
    ring = J["r_lip"] + (J["x_ring0"] - x)
    r = np.where(x >= J["x_ring0"], np.minimum(r, J["r_lip"]), np.minimum(r, ring))
    r = np.where(x >= J["x_tip"], np.minimum(r, J["r_spigot_in"]), r)
    if rail_flush:
        intr = r_out(x) - RAIL_H - T_MIN
        r = np.where((x >= RAIL_X[0] - 2) & (x <= RAIL_X[1] + 2), np.minimum(r, intr), r)
    return r


def generator_bay(J: dict, rail_flush: bool, coupler_stack: float = COUPLER_STACK) -> dict:
    """Largest cylinder (diameter for each length) in the front shell, forward of the
    joint, with the motor carrier ahead of it. Returns the frontier and the pick."""
    x_aft_max = J["x_tip"] - 0.5                         # bell face never enters the spigot
    xs = np.linspace(2.0, J["x_eq"], 2000)
    rc = front_clear_r(xs, J, rail_flush)
    frontier = []
    for L in np.arange(18.0, 40.01, 1.0):
        best = (0.0, None)
        for x_aft in np.arange(x_aft_max, 20.0, -0.25):
            x_fwd = x_aft - L
            # carrier plate ahead of the motor must carry the mount PCD
            r_carrier = float(r_in(x_fwd - CARRIER_T)) - CARRIER_RING - P.CLEAR_SLIDE
            if r_carrier < MOTOR_PCD_MIN / 2 + 3.0:
                continue
            m = (xs >= x_fwd) & (xs <= x_aft)
            R = float(rc[m].min()) - C_GEN
            R = min(R, J["r_lip"] - C_GEN)               # must pass the socket lip to go in
            if 2 * R > best[0]:
                best = (2 * R, dict(x_fwd=x_fwd, x_aft=x_aft, r_carrier=r_carrier))
        frontier.append(dict(L_mm=float(L), D_max_mm=float(best[0]), **(best[1] or {})))
    # pick: maximise D^2 L (torque capacity of an outrunner scales ~ with rotor volume)
    cand = [f for f in frontier if 20 <= f["L_mm"] <= 35 and f["D_max_mm"] > 0]
    pick = max(cand, key=lambda f: f["D_max_mm"] ** 2 * f["L_mm"])
    D = math.floor(pick["D_max_mm"])
    L = pick["L_mm"]
    # place the chosen (rounded) envelope as far aft as possible: shortest shaft
    x_aft = pick["x_aft"]
    return dict(frontier=frontier, pick=dict(D_mm=float(D), L_mm=float(L),
                x_fwd_mm=x_aft - L, x_aft_mm=x_aft, r_carrier_mm=pick["r_carrier"],
                rule="max D^2 L for 20 <= L <= 35, D rounded down to 1 mm"),
                rail_flush=rail_flush)


def fits(bay_frontier, D, L) -> bool:
    for f in bay_frontier:
        if abs(f["L_mm"] - L) < 1e-6:
            return f["D_max_mm"] >= D
    return False


# ==========================================================================
# Drivetrain
# ==========================================================================


def drivetrain(J: dict, bay: dict, ri: dict, hw: dict, pk: Pockets) -> dict:
    x_end = P.STATIC_TAIL_X
    b6, b8, sl = hw["b608"], hw["b688"], hw["seal"]
    seal = (x_end - SEAL_RECESS - sl["B"], x_end - SEAL_RECESS)
    shoulder = (seal[0] - SEAL_SHOULDER_T, seal[0])
    rear = (shoulder[0] - b8["B"], shoulder[0])
    x_bell = bay["pick"]["x_aft_mm"]
    front = (x_bell + COUPLER_STACK, x_bell + COUPLER_STACK + b6["B"])
    cB = (front[0] - 10.0, front[0])                     # claw B = thrust collar
    cA = (x_bell, x_bell + 5.0)
    shaft = (cB[0], ri["shaft_aft_x"])
    xc_f = sum(front) / 2
    xc_r = sum(rear) / 2
    span = xc_r - xc_f
    # rotor load plane: hinge + thrust centroid (0.7 R) leaned aft by the cone
    r_c = 0.7 * ri["r_tip"]
    x_load = P.HINGE_X + (r_c - P.HINGE_R) * math.tan(math.radians(ri["cone_deg"]))
    a = x_load - xc_r

    # loads at the tow case (worst)
    A_rot = math.pi * (ri["r_tip"] / 1000) ** 2
    T_tow = 0.5 * RHO * A_rot * V_TOW ** 2 * CT_TOW
    Fr = RADIAL_LOAD_FRAC * T_tow
    M_asym = THRUST_ASYM_FRAC * T_tow / P.N_BLADES * r_c / 1000      # N m
    R_rear = Fr * (span + a) / span + M_asym / (span / 1000)
    R_front = Fr * a / span + M_asym / (span / 1000)
    # 608 takes all thrust (locating); P0 = 0.6 Fr + 0.5 Fa (ISO 76 for deep groove)
    P0_front = max(R_front, 0.6 * R_front + 0.5 * T_tow)
    s0_front = b6["C0_N"] / P0_front
    s0_rear = b8["C0_N"] / R_rear
    # tilt from seat slop and tip runout
    slop = 0.05
    tilt = 2 * slop / span
    tube_r_front = b6["D"] / 2 + TUBE_WALL
    tube_r_rear = b8["D"] / 2 + TUBE_WALL
    # pocket clearance of the tube
    xx = np.linspace(front[0], x_end, 200)
    tube_r = np.where(xx <= front[1] + 2, tube_r_front,
                      np.interp(xx, [front[1] + 2, rear[0] - 2, rear[0]],
                                [tube_r_front, tube_r_rear + 1.0, tube_r_rear]))
    gap = pk.floor_inner_r(xx) - tube_r
    return dict(
        stations_mm=dict(claw_A=list(cA), claw_B_thrust_collar=list(cB), bearing_front_608=list(front),
                         bearing_rear_688=list(rear), seal_shoulder=list(shoulder), lip_seal=list(seal),
                         static_tail_end=x_end, shaft=list(shaft), rotor_load_plane=x_load),
        bearings=dict(front=b6["name"], rear=b8["name"],
                      why="608 in the wide front where it can take the thrust (C0 1.35 kN); 688 in "
                          "the narrow tail; stepped bore so both press in from the front, and the "
                          "688 passes through the 608 seat",
                      spacing_mm=span, overhang_mm=a, spacing_over_overhang=span / a,
                      radial_load_N=Fr, asym_moment_Nm=M_asym, thrust_tow_N=T_tow,
                      reaction_front_N=R_front, reaction_rear_N=R_rear,
                      static_safety_front_s0=s0_front, static_safety_rear_s0=s0_rear,
                      seat_slop_assumed_mm=slop, tilt_mrad=1000 * tilt,
                      tip_axial_runout_mm=tilt * ri["r_tip"],
                      location="608 locating: claw B (set screw on a shaft flat) bears on its inner "
                               "race for +X thrust, outer race on the tube shoulder; 688 floating "
                               "(slip fit on the shaft). -X handling loads close the 0.6 mm hub gap."),
        tube=dict(r_out_front_mm=tube_r_front, r_out_rear_mm=tube_r_rear,
                  bore_front_mm=b6["D"], shoulder_front_id_mm=FRONT_SHOULDER_ID,
                  bore_rear_mm=b8["D"], shoulder_rear_id_mm=SEAL_SHOULDER_ID, seal_bore_mm=sl["D"],
                  x_from=front[0] - 1.0, x_to=x_end,
                  min_gap_to_pocket_floor_mm=float(gap.min()),
                  ribs="3 radial ribs 2.0 mm under the pocket floors (theta = pocket angles), "
                       "x from the tube front to the tail end wall; bottom sector left free for ballast"),
        seal=dict(name=sl["name"], x=list(seal), lip="main lip faces aft (water side)",
                  note="v1 display: pocket only, seal optional"),
        shaft=dict(d_mm=P.SHAFT_D, x_fwd=shaft[0], x_aft=shaft[1], length_mm=shaft[1] - shaft[0],
                   stock="8 mm h6 stainless ground rod, cut; chamfer ends 0.5 x 30 deg; two flats "
                         "(claw B set screw, hub key)",
                   mass_g=math.pi * (P.SHAFT_D / 20) ** 2 * (shaft[1] - shaft[0]) / 10 * RHO_STEEL),
        coupler=dict(type="two-dog claw coupling, V lead-in, 1.0 mm axial float, 0.3 mm radial "
                          "clearance per dog: engages blind as the shells close; tolerates "
                          "shell-to-shell misalignment so the motor bearings and the 608/688 pair "
                          "are not forced into line",
                     od_mm=COUPLER_OD, claw_A=list(cA), claw_B=list(cB),
                     fits_spigot_bore=COUPLER_OD / 2 < J["r_spigot_in"]),
    )


# ==========================================================================
# Nose: tow eye and cable gland
# ==========================================================================


def nose(ri: dict, fins: dict) -> dict:
    src = sources()
    A_rot = math.pi * (ri["r_tip"] / 1000) ** 2
    T = 0.5 * RHO * A_rot * V_TOW ** 2 * CT_TOW
    A_hull = math.pi * (P.HULL_D / 2000) ** 2
    D_hull = 0.5 * RHO * A_hull * V_TOW ** 2 * CD_HULL
    D_fins = fins.get("drag_tow_N", 0.0)
    F_work = T + D_hull + D_fins
    F_d = SF_TOW * F_work
    mats = {}
    for mat, key in (("PETG", "prusament_petg_tds"), ("PLA", "prusament_pla_tds")):
        s = src.get(key, {})
        adh = (s.get("interlayer_adhesion_MPa") or [18 if mat == "PETG" else 17])[0]
        sxy = (s.get("tensile_yield_xy_MPa") or [47 if mat == "PETG" else 51])[0]
        # eye legs: tension across layers (front shell printed nose-up)
        sig_leg = F_d / (2 * EYE_W * EYE_T)
        # crown: pin-loaded ring, M ~ F R / pi, bending stress in the layer plane
        Rm = EYE_HOLE_D / 2 + EYE_W / 2
        sig_crown = 6 * (F_d * Rm / math.pi) / (EYE_T * EYE_W ** 2)
        # bayonet lug: transverse shear at the root is interlaminar (aft shell printed tail-down)
        mats[mat] = dict(interlayer_MPa=adh, xy_yield_MPa=sxy,
                         eye_leg_stress_MPa=sig_leg, eye_leg_util=sig_leg / adh,
                         eye_crown_stress_MPa=sig_crown, eye_crown_util=sig_crown / sxy,
                         bearing_on_8mm_pin_MPa=F_d / (EYE_HOLE_D * EYE_T))
    return dict(
        thrust_2p5_m_s_N=T, hull_drag_est_N=D_hull, fin_drag_N=D_fins, working_load_N=F_work,
        safety_factor=SF_TOW, design_load_N=F_d,
        eye=dict(type="printed ring at the nose tip, on the axis, plane = XY (horizontal)",
                 hole_d_mm=EYE_HOLE_D, ring_w_mm=EYE_W, ring_t_mm=EYE_T,
                 hole_centre_x_mm=-(EYE_HOLE_D / 2 + EYE_W + EYE_STANDOFF) + EYE_W,
                 print="front shell nose-up: ring stands vertical, hole printed as a 45 deg teardrop",
                 materials=mats,
                 note="legs load the layers in tension (weakest direction) and still sit far "
                      "below the Prusament interlayer adhesion; use a stainless shackle or a "
                      "thimble, not a bare rope, to keep the bearing stress low"),
        gland=dict(thread="M12x1.5", boss_d_mm=GLAND_BOSS_D, boss_x_mm=GLAND_X, theta_deg=180.0,
                   axis="surface normal at the boss (forward-down)",
                   thread_engagement_mm=GLAND_ENGAGE, seal="flat NBR washer under the gland head",
                   boss_surface_r_mm=float(r_out(GLAND_X)),
                   note="printed or tapped M12x1.5 in a solid boss, no inner locknut needed; "
                        "cable runs aft through a hole in the motor carrier"),
    )


# ==========================================================================
# Roll stability: fins
# ==========================================================================


def _polar_cache() -> dict:
    p = WORK / "polars_fin.json"
    d = _load(p)
    return d if isinstance(d, dict) else {}


_POLAR_MEMO: dict = {}


def fin_polar(Re: float) -> dict:
    """NeuralFoil polar for the fin section at Re (cached on disk, Re grid then interp)."""
    key = round(float(Re), -1)
    if key in _POLAR_MEMO:
        return _POLAR_MEMO[key]
    _POLAR_MEMO[key] = _fin_polar(key)
    return _POLAR_MEMO[key]


def _fin_polar(Re: float) -> dict:
    grid = [20000, 30000, 40000, 55000, 75000, 100000, 150000, 250000]
    cache = _POLAR_MEMO.setdefault("_disk", _polar_cache())
    alphas = np.arange(-6.0, 16.01, 0.5)
    missing = [g for g in grid if str(g) not in cache]
    if missing:
        import aerosandbox as asb
        import neuralfoil as nf
        af = asb.Airfoil(FIN_SECTION)
        for g in missing:
            r = nf.get_aero_from_airfoil(af, alpha=alphas, Re=float(g), n_crit=FIN_N_CRIT,
                                         model_size="large")
            cache[str(g)] = dict(alpha=alphas.tolist(), CL=np.asarray(r["CL"]).tolist(),
                                 CD=np.asarray(r["CD"]).tolist(),
                                 conf=np.asarray(r["analysis_confidence"]).tolist())
        WORK.mkdir(parents=True, exist_ok=True)
        (WORK / "polars_fin.json").write_text(json.dumps(cache))
    Re = float(np.clip(Re, grid[0], grid[-1]))
    lo = max(g for g in grid if g <= Re)
    hi = min(g for g in grid if g >= Re)
    w = 0.0 if hi == lo else (math.log(Re) - math.log(lo)) / (math.log(hi) - math.log(lo))
    a = np.array(cache[str(lo)]["alpha"])
    CL = (1 - w) * np.array(cache[str(lo)]["CL"]) + w * np.array(cache[str(hi)]["CL"])
    CD = (1 - w) * np.array(cache[str(lo)]["CD"]) + w * np.array(cache[str(hi)]["CD"])
    conf = np.minimum(np.array(cache[str(lo)]["conf"]), np.array(cache[str(hi)]["conf"]))
    m = (a >= -2) & (a <= 5)
    a0, cl0 = np.polyfit(np.radians(a[m]), CL[m], 1)
    return dict(alpha=a, CL=CL, CD=CD, conf=conf, a0_per_rad=a0, alpha_L0_deg=-math.degrees(cl0 / a0),
                CLmax=float(CL.max()), alpha_CLmax=float(a[np.argmax(CL)]), Re=Re)


def roll_torque(v: float, ri: dict) -> dict:
    R = ri["r_tip"] / 1000
    P_ = CP_DESIGN * 0.5 * RHO * math.pi * R ** 2 * v ** 3
    om = P.DESIGN_TSR * v / R
    return dict(v_m_s=v, P_W=P_, omega_rad_s=om, rpm=om * 60 / (2 * math.pi), Q_Nm=P_ / om)


def upstream_factor(x_mm: float, x_rotor: float, R_mm: float) -> float:
    """1 - a_up: on-axis actuator-disc deceleration at x (conservative off axis)."""
    a = (1 - math.sqrt(1 - CT_TOW)) / 2
    xi = (x_mm - x_rotor) / R_mm
    return 1 - a * (1 + xi / math.sqrt(1 + xi * xi))


def fin_design(ri: dict, x_rotor: float) -> dict:
    x_te = fin_te_limit(ri)
    r_root = P.HULL_D / 2 + FIN_ROOT_GAP
    best = None
    rows = []
    for c in (30.0, 35.0, 40.0):
        x_le = x_te - c
        if x_le < 2.0:
            continue
        for s in np.arange(40.0, 130.01, 5.0):
            r_cp = r_root + FIN_CP_SPAN * s
            AR = s / c
            ok = True
            res = {}
            # incidence: smallest that gives torque ratio >= 1 at both speeds
            for alpha_set in np.arange(0.0, 14.01, 0.25):
                ok = True
                res = {}
                for v in (1.0, 1.5):
                    k = upstream_factor((x_le + x_te) / 2, x_rotor, ri["r_tip"])
                    vf = k * v
                    pol = fin_polar(vf * c / 1000 / NU)
                    a0 = pol["a0_per_rad"]
                    ar = a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))
                    CLmax3 = FIN_CLMAX_3D_FACTOR * pol["CLmax"]
                    CL = ar * math.radians(alpha_set - pol["alpha_L0_deg"])
                    q = 0.5 * RHO * vf ** 2
                    S = c * s / 1e6
                    tau = 2 * q * S * CL * r_cp / 1000
                    Q = roll_torque(v, ri)["Q_Nm"]
                    res[v] = dict(v_fin_m_s=vf, Re=pol["Re"], CL=CL, CLmax3D=CLmax3,
                                  CL_alpha_3D_per_rad=ar, alpha_L0_deg=pol["alpha_L0_deg"],
                                  torque_Nm=tau, Q_gen_Nm=Q, ratio=tau / Q,
                                  authority_ratio=2 * q * S * FIN_CL_RESERVE * CLmax3 * r_cp / 1000 / Q,
                                  stall_reserve=CLmax3 / CL if CL > 0 else float("inf"))
                    if tau < Q or CL > FIN_CL_RESERVE * CLmax3:
                        ok = False
                if ok:
                    break
            row = dict(chord_mm=c, span_mm=float(s), AR=AR, area_cm2=c * s / 100, r_cp_mm=r_cp,
                       alpha_set_deg=float(alpha_set) if ok else None, feasible=ok)
            rows.append(row)
            if ok and (best is None or c * s < best[0]["chord_mm"] * best[0]["span_mm"] - 1e-9):
                best = (row, res, x_le)
    row, res, x_le = best
    # drag at the tow speed for the nose load (chosen fin, chosen incidence)
    c, s = row["chord_mm"], row["span_mm"]
    k = upstream_factor(x_le + c / 2, x_rotor, ri["r_tip"])
    vt = k * V_TOW
    pol = fin_polar(vt * c / 1000 / NU)
    a0 = pol["a0_per_rad"]
    AR = s / c
    ar = a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))
    CL = ar * math.radians(row["alpha_set_deg"] - pol["alpha_L0_deg"])
    cd0 = float(np.interp(row["alpha_set_deg"], pol["alpha"], pol["CD"]))
    CD = cd0 + CL ** 2 / (math.pi * FIN_E_OSWALD * AR)
    D2 = 2 * 0.5 * RHO * vt ** 2 * (c * s / 1e6) * CD
    # root bending at the tow case (in the layer plane: fins print flat)
    L1 = 0.5 * RHO * vt ** 2 * (c * s / 1e6) * CL
    M_root = L1 * FIN_CP_SPAN * s / 1000
    t = 0.117 * c
    Z = 0.6 * c * t ** 2 / 6                              # ESTIMATE, airfoil ~ 0.6 x rectangle
    return dict(
        section=FIN_SECTION + " (flat lower surface, t/c 0.117), printed flat, no supports",
        n=2, layout="theta = 90 and 270 deg (+-Y), identical parts; rail yaw sets +alpha on +Y "
                    "and the 180 deg-rotated copy gives the opposite lift on -Y",
        chord_mm=c, span_mm=s, AR=AR, area_each_cm2=c * s / 100, r_root_mm=r_root,
        r_tip_mm=r_root + s, r_cp_mm=row["r_cp_mm"], x_le_mm=x_le, x_te_mm=x_te,
        incidence_deg=row["alpha_set_deg"], margins=res,
        drag_tow_N=D2, root_moment_tow_Nm=M_root, root_stress_tow_MPa=M_root * 1000 / Z,
        rail=dict(type="dovetail, 45 deg flanks, male tongue on the fin, groove on a raised pad of "
                       "the front shell", x_mm=list(RAIL_X), height_mm=RAIL_H,
                  inclination="rail line = chord of the hull meridian between its ends; yawed by "
                              "the incidence about the radial axis",
                  retention="fin slides aft onto an end stop (drag pushes it in); one M3 grub "
                            "screw in a heat-set insert for handling"),
        x_te_limit_mm=x_te,
        x_te_rule="fin TE root point stays FIN_TIP_CLEAR outside the sphere swept by the blade "
                  "tip about the hinge, for any rotor angle during deployment",
        sweep=rows,
    )


def fin_te_limit(ri: dict) -> float:
    r_root = P.HULL_D / 2 + FIN_ROOT_GAP
    Lb = ri["r_tip"] - P.HINGE_R
    need = Lb + FIN_TIP_CLEAR
    dz = r_root - P.HINGE_R
    return float(P.HINGE_X - math.sqrt(max(need ** 2 - dz ** 2, 0.0)))


def roll_summary(ri: dict, fins: dict, mb: dict) -> dict:
    out = {}
    for v in (1.0, 1.5, 2.5):
        out[f"{v:g}_m_s"] = roll_torque(v, ri)
    W = mb["v2"]["displacement_N"]
    BG = mb["v2"]["BG_mm"] / 1000
    Qd = out["1.5_m_s"]["Q_Nm"]
    return dict(
        torque=out, Q_scaling="Q = Cp 0.5 rho pi R^3 v^2 / TSR: grows with v^2 like fin lift, so a "
                              "fin set matched at one speed stays matched at constant TSR and Cp",
        pendulum_righting_max_Nm=W * BG, pendulum_vs_Qdesign=W * BG / Qd,
        heel_for_20pct_mismatch_deg=math.degrees(math.asin(min(1.0, 0.2 * Qd / (W * BG))))
        if W * BG > 0 else None,
        freewheel_note="with no electrical load Q ~ 0 and the fins spin the hull the other way; "
                       "roll damping of the fins limits it to p ~ alpha v / r_cp. v2 needs the "
                       "controller to keep a load (MPPT near TSR 4, dump load when the battery is "
                       "full) or a torque-carrying mount (pole/bracket) where fins are not needed",
    )


# ==========================================================================
# Mass and buoyancy
# ==========================================================================


def _rev_vol(x0, x1, rfun_out, rfun_in=None, n=800):
    x = np.linspace(x0, x1, n)
    ro = np.asarray(rfun_out(x), dtype=float)
    ri_ = np.zeros_like(ro) if rfun_in is None else np.asarray(rfun_in(x), dtype=float)
    a = math.pi * (ro ** 2 - np.minimum(ri_, ro) ** 2)
    v = np.trapezoid(a, x)
    xc = np.trapezoid(a * x, x) / v if v > 0 else (x0 + x1) / 2
    return float(v), float(xc)


def _solid_part(V, A_surf):
    """Printed volume of a 'solid' part: 3-wall skin + 15 % gyroid. mm^3."""
    skin = min(V, A_surf * SKIN)
    return skin + INFILL * (V - skin)


def mass_buoyancy(J, bay, DT, fins, ri, hw, pk: Pockets, nose_d) -> dict:
    x_eq = J["x_eq"]
    x_end = P.STATIC_TAIL_X
    parts = []

    def add(name, group, V_print_mm3, x_c, z_c=0.0, mat="print", V_disp_mm3=None, note=""):
        parts.append(dict(name=name, group=group, V_print_mm3=V_print_mm3, x_mm=x_c, z_mm=z_c,
                          mat=mat, V_disp_mm3=V_disp_mm3, note=note))

    # ---- front shell
    v, xc = _rev_vol(0, x_eq, r_out, r_in)
    add("front shell wall", "F", v, xc)
    v, xc = _rev_vol(J["x_ring0"] - 4, x_eq, r_in,
                     lambda x: np.where(x >= J["x_seal0"], J["r_bore"],
                                        np.where(x >= J["x_ring0"], J["r_lip"],
                                                 np.minimum(r_in(x), J["r_lip"] + (J["x_ring0"] - x)))))
    add("socket ring", "F", v, xc)
    Rm = EYE_HOLE_D / 2 + EYE_W / 2
    add("tow eye", "F", 2 * math.pi * Rm * EYE_W * EYE_T, -6.0)
    add("gland boss", "F", math.pi / 4 * (GLAND_BOSS_D ** 2 - 12 ** 2) * GLAND_ENGAGE, GLAND_X,
        -r_out(GLAND_X) + 4)
    xcar = bay["pick"]["x_fwd_mm"] - CARRIER_T / 2
    add("carrier ring", "F", 2 * math.pi * float(r_in(xcar)) * CARRIER_RING * 7.0, xcar)
    rail_len = RAIL_X[1] - RAIL_X[0]
    add("fin rail pads x2", "F", 2 * rail_len * 10.0 * RAIL_H, sum(RAIL_X) / 2)
    xp0 = max(pk.x0, 0)
    lp_F = max(0.0, x_eq - xp0)
    lp_A = x_end - max(x_eq, pk.x0)
    dpk = float(pk.ri["depth"])
    dV_pocket = pk.width * POCKET_FLOOR_WALL + 2 * dpk * T_MIN - pk.width * P.HULL_WALL
    add("pocket walls (F part)", "F", 3 * lp_F * dV_pocket, (xp0 + x_eq) / 2)
    # ---- aft shell
    v, xc = _rev_vol(x_eq, x_end, r_out, r_in)
    add("aft shell wall", "A", v, xc)
    add("pocket walls (A part)", "A", 3 * lp_A * dV_pocket, (x_eq + x_end) / 2)
    r_end = float(r_in(x_end - 2.5))
    add("tail end wall + seal boss", "A", math.pi * (r_end ** 2 - 4.2 ** 2) * 5.0, x_end - 2.5)
    vs = math.pi * ((J["r_spigot_seal"] ** 2 - J["r_spigot_in"] ** 2) * J["L_seal"]
                    + (J["r_spigot_lug"] ** 2 - J["r_spigot_in"] ** 2) * (LIP_T + LUG_T))
    vs += LUG_N * math.radians(LUG_W_DEG) * J["r_spigot_lug"] * LUG_H * LUG_T * 0.75
    add("spigot + lugs", "A", vs, x_eq - J["L_spigot"] / 2)
    t = DT["tube"]
    x0t, x1t = t["x_from"], t["x_to"]
    xm = DT["stations_mm"]["bearing_front_608"][1] + 2
    v1 = math.pi * (t["r_out_front_mm"] ** 2 - (t["bore_front_mm"] / 2) ** 2) * (xm - x0t)
    v2 = math.pi * ((t["r_out_front_mm"] + t["r_out_rear_mm"]) ** 2 / 4 - (FRONT_SHOULDER_ID / 2) ** 2) * (x1t - xm)
    add("bearing tube", "A", v1 + v2, (x0t + x1t) / 2)
    xr = np.linspace(x0t, x1t, 100)
    rib_len = np.trapezoid(np.maximum(pk.floor_inner_r(xr) - t["r_out_front_mm"] + 1, 0), xr)
    add("ribs x3", "A", 3 * 2.0 * rib_len, (x0t + x1t) / 2)
    # ---- internal printed parts (solid parts: skin + gyroid)
    rc = bay["pick"]["r_carrier_mm"]
    Vc = math.pi * rc ** 2 * CARRIER_T
    add("motor carrier", "int", _solid_part(Vc, 2 * math.pi * rc ** 2 + 2 * math.pi * rc * CARRIER_T),
        xcar, V_disp_mm3=0.0)
    for nm, x0_, x1_, r_ in (("claw A", *DT["coupler"]["claw_A"], COUPLER_OD / 2),
                             ("claw B", *DT["coupler"]["claw_B"], COUPLER_OD / 2)):
        V = math.pi * (r_ ** 2 - 4 ** 2) * (x1_ - x0_) * 0.7
        add(nm, "int", _solid_part(V, 2 * math.pi * r_ ** 2 + 2 * math.pi * r_ * (x1_ - x0_)),
            (x0_ + x1_) / 2, V_disp_mm3=0.0)
    Dg, Lg = bay["pick"]["D_mm"], bay["pick"]["L_mm"]
    Vg = math.pi * Dg ** 2 / 4 * Lg
    xg = (bay["pick"]["x_fwd_mm"] + bay["pick"]["x_aft_mm"]) / 2
    gen_replica = _solid_part(Vg, math.pi * Dg ** 2 / 2 + math.pi * Dg * Lg)
    # ---- fins (outside; displace their own volume)
    Vf = 2 * (0.074 * fins["chord_mm"] ** 2 * fins["span_mm"] + RAIL_H * 10 * rail_len)
    Af = 2 * (2 * fins["chord_mm"] * fins["span_mm"])
    add("fins x2", "fin", _solid_part(Vf, Af), (fins["x_le_mm"] + fins["x_te_mm"]) / 2,
        V_disp_mm3=Vf)
    # ---- rotor parts (ESTIMATES until the rotor study publishes masses)
    xh0 = x_end + P.HUB_GAP
    v_h, x_h = _rev_vol(xh0, P.HULL_L, r_out)
    sa_h = 2 * math.pi * np.trapezoid(r_out(np.linspace(xh0, P.HULL_L, 200)), np.linspace(xh0, P.HULL_L, 200))
    add("hub + spinner (ESTIMATE)", "rotor", _solid_part(v_h, float(sa_h) + math.pi * 20 ** 2),
        x_h, V_disp_mm3=v_h)
    rr = np.linspace(30.0, ri["r_tip"], 60)
    cc = np.interp(rr, [30, ri["r_tip"]], [26, 14])
    tc = np.interp(rr, [30, ri["r_tip"]], [0.18, 0.12])
    Vb = float(np.trapezoid(0.685 * tc * cc ** 2, rr)) + 22 * 8 * (30 - P.HINGE_R)
    Ab = float(np.trapezoid(2.05 * cc, rr))
    xb = P.HINGE_X + (0.5 * (30 + ri["r_tip"]) - P.HINGE_R) * math.tan(math.radians(ri["cone_deg"]))
    add("blades x3 (ESTIMATE)", "rotor", 3 * _solid_part(Vb, Ab), xb, V_disp_mm3=3 * Vb)

    # ---- hardware (grams directly)
    sh = DT["shaft"]
    hwl = [
        ("shaft 8 mm", sh["mass_g"], (sh["x_fwd"] + sh["x_aft"]) / 2, "computed, steel 7.9 g/cc ESTIMATE"),
        ("608-2RS", hw["b608"]["mass_g"], sum(DT["stations_mm"]["bearing_front_608"]) / 2, "catalogue"),
        ("688-2RS", hw["b688"]["mass_g"], sum(DT["stations_mm"]["bearing_rear_688"]) / 2, "catalogue, VERIFY"),
        ("lip seal 8x16x5", hw["seal"]["mass_g"], sum(DT["stations_mm"]["lip_seal"]) / 2, "ESTIMATE"),
        ("O-ring", RHO_NBR * math.pi ** 2 * (J["gland"]["oring_id_mm"] + OR_CS) * OR_CS ** 2 / 4 / 1000,
         x_eq - 3, "computed, NBR 1.3 g/cc ESTIMATE"),
        ("magnets x10 (6 blade detents, 4 joint)", 10 * RHO_NDFEB * math.pi * 4 * 2 / 1000, 90.0,
         "computed, NdFeB 7.5 g/cc ESTIMATE; x is a guess"),
        ("hinge pins 3x16 x3", 3 * RHO_STEEL * math.pi * 1.5 ** 2 * 16 / 1000, P.HINGE_X, "computed"),
        ("M3 screws, grub screws, inserts", 12 * 0.5, 40.0, "ESTIMATE"),
        ("M12 cable gland", hw["gland"]["mass_g"], GLAND_X, "ESTIMATE"),
    ]
    V_static = _rev_vol(0, x_end, r_out)[0]
    V_pockets = 3 * pk.width * float(np.trapezoid(pk.depth(np.linspace(pk.x0, pk.x1, 300)),
                                                  np.linspace(pk.x0, pk.x1, 300)))
    x_static = _rev_vol(0, x_end, r_out)[1]

    def build(mat_rho, generator: str, ballast_g: float = 0.0, fins_on: bool = True):
        m = 0.0
        mx = 0.0
        mz = 0.0
        rows = []
        for p in parts:
            if p["group"] == "fin" and not fins_on:
                continue
            g = p["V_print_mm3"] * mat_rho / 1000
            rows.append((p["name"], g, p["x_mm"]))
            m += g
            mx += g * p["x_mm"]
            mz += g * p["z_mm"]
        for nm, g, x, _ in hwl:
            rows.append((nm, g, x))
            m += g
            mx += g * x
        if generator == "replica":
            g = gen_replica * mat_rho / 1000
        else:
            g = hw["motor"]["mass_g"] + 10.0       # + 10 g leads (ESTIMATE)
        rows.append((f"generator ({generator})", g, xg))
        m += g
        mx += g * xg
        if ballast_g > 0:
            xb_ = 0.5 * (DT["stations_mm"]["bearing_front_608"][1] + DT["stations_mm"]["bearing_rear_688"][0])
            zb_ = -(float(pk.floor_inner_r(xb_)) - 4.0)      # bottom sector, ~4 mm above the wall
            rows.append(("lead ballast", ballast_g, xb_))
            m += ballast_g
            mx += ballast_g * xb_
            mz += ballast_g * zb_
        V_disp = V_static - V_pockets + sum(p["V_disp_mm3"] or 0.0 for p in parts
                                            if p["V_disp_mm3"] and (fins_on or p["group"] != "fin"))
        Vx = (V_static * x_static - V_pockets * (pk.x0 + pk.x1) / 2
              + sum((p["V_disp_mm3"] or 0.0) * p["x_mm"] for p in parts
                    if p["V_disp_mm3"] and (fins_on or p["group"] != "fin")))
        disp_g = V_disp * RHO / 1e6
        return dict(mass_g=m, x_cg_mm=mx / m, z_cg_mm=mz / m, displacement_cc=V_disp / 1000,
                    displacement_g=disp_g, x_cb_mm=Vx / V_disp, z_cb_mm=0.0,
                    net_buoyancy_N=(disp_g - m) * G / 1000, displacement_N=disp_g * G / 1000,
                    BG_mm=0.0 - mz / m, rows=rows)

    disp_pla = build(RHO_PLA, "replica", fins_on=False)
    v2_noball = build(RHO_PETG, "real", fins_on=True)
    # v2 target: net buoyancy +0.3 N with ballast in the bottom sector (keeps it rising
    # slowly in slack water and gives the pendulum BG)
    target_N = 0.3
    ballast = max(0.0, (v2_noball["net_buoyancy_N"] - target_N) * 1000 / G)
    v2 = build(RHO_PETG, "real", ballast_g=ballast, fins_on=True)
    xb_ = 0.5 * (DT["stations_mm"]["bearing_front_608"][1] + DT["stations_mm"]["bearing_rear_688"][0])
    lead_cc = ballast / RHO_LEAD
    return dict(
        model=dict(shells="6 perimeters = solid 2.4 mm wall", solid_parts=f"{SKIN} mm skin + "
                   f"{int(INFILL * 100)} % gyroid", rotor="ESTIMATE from BEM planform (chord 26->14, "
                   "t/c 0.18->0.12) and the hull profile for hub+spinner", sealed="displacement counts "
                   "the sealed static hull (minus wet pockets) plus the outer volume of rotor parts and "
                   "fins (closed perimeters trap the infill air)"),
        parts=[dict(name=p["name"], group=p["group"], V_print_cc=p["V_print_mm3"] / 1000,
                    m_PLA_g=p["V_print_mm3"] * RHO_PLA / 1000, m_PETG_g=p["V_print_mm3"] * RHO_PETG / 1000,
                    x_mm=p["x_mm"]) for p in parts] + [dict(name="generator replica", group="int",
                    V_print_cc=gen_replica / 1000, m_PLA_g=gen_replica * RHO_PLA / 1000,
                    m_PETG_g=gen_replica * RHO_PETG / 1000, x_mm=xg)],
        hardware=[dict(name=n, mass_g=g, x_mm=x, basis=b) for n, g, x, b in hwl],
        static_hull_outer_cc=V_static / 1000, pockets_cc=V_pockets / 1000,
        display_PLA=_strip(disp_pla), v2_PETG_no_ballast=_strip(v2_noball), v2=_strip(v2),
        v2_ballast=dict(lead_g=ballast, lead_cc=lead_cc, x_mm=xb_, target_net_buoyancy_N=target_N,
                        where="bottom sector of the aft shell between the ribs at 120 and 240 deg, "
                              "cast or shot-in-epoxy in a printed tray; trim by sliding it along x"),
    )


def _strip(d):
    return {k: v for k, v in d.items() if k != "rows"} | dict(
        breakdown=[dict(name=n, g=g, x_mm=x) for n, g, x in d["rows"]])


# ==========================================================================
# Assembly and part split
# ==========================================================================


def assembly(J, bay, DT) -> list:
    st = DT["stations_mm"]
    return [
        dict(step=1, sub="aft shell", action="Press the 688 down the tube from the open front, through "
             "the 608 seat (22 mm) and the 18 mm step, onto the rear shoulder (draw-bolt + washer, "
             "not a hammer).", check="688 OD 16 < step 18 mm: passes. OK"),
        dict(step=2, sub="aft shell", action="Press the lip seal into the tail-end bore from the aft "
             "face, main lip facing aft, onto the same shoulder (v2; optional on the display).",
             check="seal bore open aft, shoulder 1.5 mm. OK"),
        dict(step=3, sub="aft shell", action="Press the 608 into the front seat onto its 18 mm "
             "shoulder.", check="both bearings face shoulders that point forward: no retaining "
             "ring needed, tube printed tail-down so both shoulders print as floors. OK"),
        dict(step=4, sub="aft shell", action="Slide the shaft in from the front through 608, 688 and "
             "the seal (chamfered end, tape or a sleeve over the flat at the seal).",
             check="shaft passes the seal against its back face: needs the chamfer and a sleeve, flag"),
        dict(step=5, sub="aft shell", action="Fit claw B on the forward shaft end, push it against the "
             "608 inner race, tighten its grub screw onto the flat (hex key through the spigot bore).",
             check=f"spigot bore D {2 * J['r_spigot_in']:.1f} mm gives straight access. OK"),
        dict(step=6, sub="rotor", action="Rotor colleague: hub on the aft shaft end, blades on their "
             "pins, spinner retains the pins. Set axial play by sliding the hub before its key "
             "screw is tightened.", check="see layout_rotor"),
        dict(step=7, sub="aft shell", action="v2: drop the lead tray into the bottom sector between "
             "ribs and glue or screw it.", check="sector free of ribs. OK"),
        dict(step=8, sub="aft shell", action="Grease and fit the O-ring on the spigot groove; fit the "
             "joint magnets (glue).", check="O-ring stretch within 1-5 %. OK"),
        dict(step=9, sub="front shell", action="Screw the M12 gland into the nose boss with its "
             "washer; pass the cable.", check="done before the carrier blocks the nose cavity. OK"),
        dict(step=10, sub="front shell", action="Bench: bolt the motor stator to the carrier (M3 from "
             "the carrier front face), bolt claw A to the bell, solder the leads to the cable through "
             "the carrier hole.", check="all screws reachable on the bench. OK"),
        dict(step=11, sub="front shell", action="Push carrier + motor into the front shell through the "
             "socket, twist 30 deg into the carrier slots until it clicks.",
             check=f"motor D {bay['pick']['D_mm']:.0f} passes the socket lip D {2 * J['r_lip']:.1f}; "
                   "reaction torque drives the tabs into their hard stop. OK"),
        dict(step=12, sub="close", action="Align the entry marks, push the shells together (claws "
             "find each other on the V lead-in; turn the rotor by hand if they sit tip to tip), "
             f"twist {ENGAGE_DEG:.0f} deg to the stop; magnets click.",
             check="claw engagement is blind: V lead-in + free rotor makes it self-finding, flag for "
                   "the prototype"),
        dict(step=13, sub="fins", action="Slide each fin aft onto its rail to the stop, set the grub "
             "screw. Remove for transport.", check="fins only in tethered use. OK"),
        dict(step=14, sub="service", action="Disassembly is the reverse; bearings come out with a "
             "blind puller or by pushing the shaft (608) and a drift from the tail (688).",
             check="688 is behind the seal: remove the seal first. OK"),
    ]


def part_split(J, bay, fins) -> list:
    eye_out = EYE_HOLE_D / 2 + EYE_W + EYE_STANDOFF + EYE_HOLE_D / 2
    return [
        dict(part="front shell", material="PETG (v2) / PLA (display)",
             bbox_mm=[P.HULL_D + 2 * RAIL_H, P.HULL_D, J["x_eq"] + eye_out],
             orientation="equator face on the bed, nose up", supports="none: dome overhang is "
             "self-supporting with a 2.4 mm wall, socket and carrier-slot ceilings chamfered 45 deg, "
             "eye hole teardrop, rails vertical"),
        dict(part="aft shell", material="PETG / PLA",
             bbox_mm=[P.HULL_D, P.HULL_D, P.STATIC_TAIL_X - J["x_eq"] + J["L_spigot"]],
             orientation="tail end face on the bed, spigot up", supports="none: outer flank leans "
             f"<= {math.degrees(math.atan(abs(float(slope(P.STATIC_TAIL_X))))):.0f} deg from vertical, bearing "
             "shoulders print as floors, lug load faces and O-ring groove flank at 45 deg, seal "
             "shoulder 1.75 mm bridge"),
        dict(part="motor carrier", bbox_mm=[2 * bay["pick"]["r_carrier_mm"] + 2] * 2 + [CARRIER_T],
             orientation="flat", supports="none"),
        dict(part="claw A / claw B", bbox_mm=[COUPLER_OD, COUPLER_OD, 10.0],
             orientation="claw faces up", supports="none (V lead-in 45 deg)"),
        dict(part="generator replica (display)", bbox_mm=[bay["pick"]["D_mm"]] * 2 + [bay["pick"]["L_mm"]],
             orientation="base down", supports="none"),
        dict(part="fin x2", bbox_mm=[fins["span_mm"] + RAIL_H + 8, fins["chord_mm"] + 10, 0.117 * fins["chord_mm"]],
             orientation="flat side of the Clark Y on the bed", supports="none"),
        dict(part="ballast tray (v2)", bbox_mm=[30, 20, 40], orientation="open side up", supports="none"),
    ]


# ==========================================================================
# Figures
# ==========================================================================


def _style(plt):
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 9, "axes.titleweight": "bold", "savefig.dpi": 200})


def figures(J, bay, DT, fins, pk, nose_d, mb) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon

    _style(plt)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    C_HULL, C_WALL, C_PK = "#5B6670", "#AEB6BD", "#E8622C"
    C_GEN, C_STEEL, C_BRG, C_SEAL, C_RUB = "#3A6EA5", "#7D858C", "#2F3337", "#1A1A1A", "#111111"
    x = np.linspace(0, P.HULL_L, 1200)
    ro = r_out(x)
    rin = np.array(r_in(x))
    xs = x[x <= P.STATIC_TAIL_X]

    def hull(ax, top_pocket=True, xlim=None):
        ax.plot(x, ro, color=C_HULL, lw=1.4)
        ax.plot(x, -ro, color=C_HULL, lw=1.4)
        # walls (static part)
        rin_s = np.array(r_in(xs))
        ax.fill_between(xs, rin_s, r_out(xs), color=C_WALL, lw=0)
        ax.fill_between(xs, -r_out(xs), -rin_s, color=C_WALL, lw=0)
        if top_pocket:
            xp = np.linspace(pk.x0, pk.x1, 300)
            d = pk.depth(xp)
            ax.fill_between(xp, r_out(xp) - d, r_out(xp), color="white", lw=0, zorder=3)
            ax.plot(xp, r_out(xp) - d, color=C_PK, lw=1.0, zorder=4)
            ax.plot(xp, pk.floor_inner_r(xp), color=C_PK, lw=0.8, ls="--", zorder=4)
            ax.fill_between(xp, pk.floor_inner_r(xp), r_out(xp) - d, color=C_WALL, lw=0, zorder=3)
        # rotating part (hub + spinner), outline only
        xh = x[x >= P.STATIC_TAIL_X + P.HUB_GAP]
        ax.plot(xh, r_out(xh), color=C_HULL, lw=1.0, ls=(0, (3, 2)))
        ax.plot(xh, -r_out(xh), color=C_HULL, lw=1.0, ls=(0, (3, 2)))
        ax.axvline(J["x_eq"], color="#999", lw=0.6, ls=":")
        ax.axhline(0, color="#bbb", lw=0.5, ls="-.")

    def rect(ax, x0, x1, r0, r1, fc, ec=None, alpha=1.0, z=5, label=None, both=True):
        for s in ((1, -1) if both else (1,)):
            lo, hi = sorted((s * r0, s * r1))
            ax.add_patch(Rectangle((x0, lo), x1 - x0, hi - lo, fc=fc, ec=ec or fc, lw=0.6,
                                   alpha=alpha, zorder=z, label=label if s == 1 else None))

    def internals(ax):
        st = DT["stations_mm"]
        g = bay["pick"]
        # joint: socket ring (F) and spigot (A)
        xr = np.linspace(J["x_ring0"] - 6, J["x_eq"], 200)
        r_sock = np.where(xr >= J["x_seal0"], J["r_bore"],
                          np.where(xr >= J["x_ring0"], J["r_lip"],
                                   np.minimum(r_in(xr), J["r_lip"] + (J["x_ring0"] - xr))))
        for s in (1, -1):
            ax.fill_between(xr, s * r_sock, s * np.array(r_in(xr)), color="#C9B99A", lw=0, zorder=4)
        rect(ax, J["x_seal0"], J["x_eq"], J["r_spigot_in"], J["r_spigot_seal"], "#9C8B6E", z=5,
             label="spigot (aft shell) / socket ring (front shell)")
        rect(ax, J["x_tip"], J["x_seal0"], J["r_spigot_in"], J["r_spigot_lug"], "#9C8B6E", z=5)
        ax.add_patch(Rectangle((J["x_tip"], J["r_spigot_lug"]), LUG_T, LUG_H, fc="#6F5F45", lw=0, zorder=6,
                               label="bayonet lug (lock position)"))
        gr = J["groove"]
        rect(ax, gr["x_fwd_bottom"], gr["x_aft"], J["r_groove"], J["r_bore"], C_RUB, z=7, label="O-ring")
        # carrier + generator + claws
        rect(ax, g["x_fwd_mm"] - CARRIER_T, g["x_fwd_mm"], 0, g["r_carrier_mm"], "#7A9CC6", z=5,
             label="motor carrier (twist-lock)")
        rect(ax, g["x_fwd_mm"], g["x_aft_mm"], 0, g["D_mm"] / 2, C_GEN, alpha=0.85, z=5,
             label=f"generator envelope D{g['D_mm']:.0f} x {g['L_mm']:.0f}")
        rect(ax, *st["claw_A"], 0, COUPLER_OD / 2, "#8FA9C9", z=6, label="claw coupling A / B")
        rect(ax, *st["claw_B_thrust_collar"], 0, COUPLER_OD / 2, "#8FA9C9", z=6)
        # bearing tube
        t = DT["tube"]
        f0, f1 = st["bearing_front_608"]
        r0_, r1_ = st["bearing_rear_688"]
        xt = np.array([t["x_from"], f1 + 2, r0_ - 2, P.STATIC_TAIL_X])
        rt = np.array([t["r_out_front_mm"], t["r_out_front_mm"], t["r_out_rear_mm"] + 1, t["r_out_rear_mm"]])
        for s in (1, -1):
            ax.fill_between(xt, s * np.interp(xt, xt, rt), s * np.full_like(xt, 4.6), color="#D7D0C4",
                            lw=0, zorder=4)
        rect(ax, t["x_from"], P.STATIC_TAIL_X, 4.3, 0, "white", z=4.5)
        rect(ax, f0, f1, 4.0, t["bore_front_mm"] / 2, C_BRG, z=7, label="608-2RS (front) / 688-2RS (rear)")
        rect(ax, r0_, r1_, 4.0, t["bore_rear_mm"] / 2, C_BRG, z=7)
        rect(ax, *st["lip_seal"], 4.0, t["seal_bore_mm"] / 2, "#5A4A3A", z=7, label="lip seal 8x16x5 (v2)")
        # ribs (top, under the pocket at theta 0)
        xr2 = np.linspace(t["x_from"] + 1, P.STATIC_TAIL_X - 5, 100)
        ax.fill_between(xr2, np.interp(xr2, xt, rt), pk.floor_inner_r(xr2), color="#E4DED3", lw=0,
                        zorder=3.5, label="rib (under each pocket)")
        # shaft
        ax.add_patch(Rectangle((st["shaft"][0], -4), st["shaft"][1] - st["shaft"][0], 8, fc=C_STEEL,
                               lw=0, zorder=8, label=f"shaft 8 x {st['shaft'][1] - st['shaft'][0]:.0f}"))
        # ballast (bottom)
        b = mb["v2_ballast"]
        if b["lead_g"] > 0:
            L_b = 34.0
            h = b["lead_cc"] * 1000 / (L_b * 22.0)
            rb_ = float(pk.floor_inner_r(b["x_mm"]))
            rb_ = float(r_in(b["x_mm"])) - 0.5
            ax.add_patch(Rectangle((b["x_mm"] - L_b / 2, -rb_), L_b, h, fc="#888", ec="#444", lw=0.5,
                                   zorder=6, label=f"v2 lead ballast {b['lead_g']:.0f} g"))
        # nose: eye and gland
        e = nose_d["eye"]
        xc = e["hole_centre_x_mm"]
        th = np.linspace(0, 2 * np.pi, 80)
        Ro = EYE_HOLE_D / 2 + EYE_W
        ax.fill(xc + Ro * np.cos(th), Ro * np.sin(th) * 0 + (EYE_T / 2) * np.sign(np.sin(th)) * 0 +
                Ro * np.sin(th), color=C_WALL, zorder=5)
        ax.fill(xc + EYE_HOLE_D / 2 * np.cos(th), EYE_HOLE_D / 2 * np.sin(th), color="white", zorder=6)
        ax.plot(xc + Ro * np.cos(th), Ro * np.sin(th), color=C_HULL, lw=0.8, zorder=6)
        # gland on the bottom, along the surface normal
        xg_ = GLAND_X
        rg = float(r_out(xg_))
        sl = float(slope(xg_))
        n = np.array([-sl, 1.0]) / math.hypot(sl, 1.0)        # outward normal in (x, r)
        tdir = np.array([1.0, sl]) / math.hypot(1.0, sl)
        base = np.array([xg_, rg])
        pts = [base + tdir * 9 + n * 0, base - tdir * 9, base - tdir * 9 + n * 20, base + tdir * 9 + n * 20]
        poly = [(p_[0], -p_[1]) for p_ in pts]
        ax.add_patch(Polygon(poly, closed=True, fc="#444", ec="#222", lw=0.5, zorder=6,
                             label="M12 cable gland (theta 180)"))
        pts_in = [base + tdir * 6, base - tdir * 6, base - tdir * 6 - n * GLAND_ENGAGE,
                  base + tdir * 6 - n * GLAND_ENGAGE]
        ax.add_patch(Polygon([(p_[0], -p_[1]) for p_ in pts_in], closed=True, fc="#9A9A9A", lw=0, zorder=5))

    # ---- figure 1: side section
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    hull(ax)
    internals(ax)
    # fin silhouette (theta 90 is out of this plane: drawn as a ghost projected)
    fx0, fx1 = fins["x_le_mm"], fins["x_te_mm"]
    ax.add_patch(Rectangle((fx0, fins["r_root_mm"]), fx1 - fx0, fins["span_mm"], fc="none", ec="#2F3337",
                           lw=0.8, ls="--", zorder=2, label="fin at theta 90 (projected)"))
    ax.set_aspect("equal")
    ax.set_xlim(-22, 170)
    ax.set_ylim(-45, 45 + 0 * fins["span_mm"])
    ax.set_xlabel("x (mm), flow ->")
    ax.set_ylabel("z (mm)")
    ax.set_title("HydroEgg internal layout, side section (top half cut through a blade pocket at "
                 "theta 0, bottom half between pockets)")
    ax.annotate("equator joint", (J["x_eq"], 40), ha="center", fontsize=7, color="#555")
    ax.annotate("static tail end", (P.STATIC_TAIL_X, 40), ha="center", fontsize=7, color="#555")
    ax.annotate("pocket floor (6 mm)", (85, float(r_out(85)) - 7.5), fontsize=6.5, color=C_PK)
    h, lab = ax.get_legend_handles_labels()
    ax.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, fontsize=6.5, frameon=False)
    fig.tight_layout()
    p1 = FIG_DIR / "layout_section.png"
    fig.savefig(p1, bbox_inches="tight")
    plt.close(fig)
    out.append(p1)

    # ---- figure 2: joint detail + bay frontier
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.0), gridspec_kw=dict(width_ratios=[1.25, 1]))
    hull(a1)
    internals(a1)
    a1.set_aspect("equal")
    a1.set_xlim(J["x_ring0"] - 8, J["x_eq"] + 12)
    a1.set_ylim(J["r_spigot_in"] - 4, 39)
    a1.set_title(f"Equator joint at x = {J['x_eq']:.1f} (top: under a pocket)")
    a1.set_xlabel("x (mm)")
    a1.set_ylabel("r (mm)")
    for lbl, rr_ in (("seal bore", J["r_bore"]), ("lip", J["r_lip"]), ("spigot ID", J["r_spigot_in"])):
        a1.annotate(f"{lbl} D{2 * rr_:.1f}", (J["x_eq"] + 1, rr_), fontsize=6.5, va="center")
    fr = bay["frontier"]
    Ls = [f["L_mm"] for f in fr]
    a2.plot(Ls, [f["D_max_mm"] for f in fr], color=C_GEN, lw=1.6, label="raised fin rails (chosen)")
    if bay.get("frontier_flush"):
        a2.plot(Ls, [f["D_max_mm"] for f in bay["frontier_flush"]], color=C_GEN, lw=1.0, ls="--",
                label="flush fin slots")
    a2.add_patch(Rectangle((20, 40), 15, 16, fc="#ddd", ec="none", zorder=0, label="candidate box D40-56 x L20-35"))
    g = bay["pick"]
    a2.plot([g["L_mm"]], [g["D_mm"]], "o", color="#E8622C", label=f"chosen D{g['D_mm']:.0f} x {g['L_mm']:.0f}")
    a2.plot([22.7, 32.3], [47, 47], "s", color="#2F3337", ms=4, label="GM4108H (no enc / enc)")
    a2.set_xlabel("generator length L (mm)")
    a2.set_ylabel("max diameter D (mm)")
    a2.set_title("Generator bay frontier (front shell)")
    a2.legend(fontsize=6.5, frameon=False, loc="lower left")
    fig.tight_layout()
    p2 = FIG_DIR / "layout_section_joint_bay.png"
    fig.savefig(p2, bbox_inches="tight")
    plt.close(fig)
    out.append(p2)

    # ---- figure 3: cross sections
    st = DT["stations_mm"]
    cuts = [(30.0, "generator bay"), (J["x_seal0"] + 2, "joint seal zone"),
            (sum(st["bearing_front_608"]) / 2, "608"), (sum(st["bearing_rear_688"]) / 2, "688")]
    fig, axs = plt.subplots(1, 4, figsize=(10.5, 3.0))
    th = np.linspace(0, 2 * np.pi, 361)
    for ax, (xc, lbl) in zip(axs, cuts):
        R = float(r_out(xc))
        Ri = float(r_in(xc))
        d = float(pk.depth(xc))
        # outer with pockets
        rr = np.full_like(th, R)
        rr_in = np.full_like(th, Ri)
        for tp in pk.theta:
            dth = np.angle(np.exp(1j * (th - math.radians(tp))))
            hw_ = pk.width / 2 / R
            m = np.abs(dth) <= hw_
            rr[m] = R - d
            rr_in[m] = np.minimum(Ri, float(pk.floor_inner_r(xc)))
        Y = lambda r_: r_ * np.sin(th)  # noqa: E731
        Z = lambda r_: r_ * np.cos(th)  # noqa: E731
        ax.fill(np.r_[Y(rr), Y(rr_in)[::-1]], np.r_[Z(rr), Z(rr_in)[::-1]], color=C_WALL, lw=0)
        ax.plot(Y(rr), Z(rr), color=C_HULL, lw=1.0)
        if lbl == "generator bay":
            g = bay["pick"]
            ax.add_patch(plt.Circle((0, 0), g["D_mm"] / 2, color=C_GEN, alpha=0.85))
            for s_ in (1, -1):
                ax.add_patch(Rectangle((s_ * R - (RAIL_H if s_ < 0 else 0), -5), RAIL_H, 10, fc="#444"))
            ax.plot([R + RAIL_H, R + RAIL_H + 12], [0, 0], color="#2F3337", lw=2.5)
            ax.plot([-R - RAIL_H, -R - RAIL_H - 12], [0, 0], color="#2F3337", lw=2.5)
        elif lbl == "joint seal zone":
            ax.add_patch(plt.Circle((0, 0), J["r_bore"], fill=False, ec="#9C8B6E", lw=2.0))
            ax.add_patch(plt.Circle((0, 0), J["r_spigot_in"], fill=False, ec="#9C8B6E", lw=1.0))
            ax.add_patch(plt.Circle((0, 0), J["r_groove"] + 1, fill=False, ec="#111", lw=1.6))
            ax.add_patch(plt.Circle((0, 0), COUPLER_OD / 2, color="#8FA9C9"))
        else:
            t = DT["tube"]
            ro_t = t["r_out_front_mm"] if lbl == "608" else t["r_out_rear_mm"]
            bo = (t["bore_front_mm"] if lbl == "608" else t["bore_rear_mm"]) / 2
            ax.add_patch(plt.Circle((0, 0), ro_t, color="#D7D0C4"))
            ax.add_patch(plt.Circle((0, 0), bo, color=C_BRG))
            ax.add_patch(plt.Circle((0, 0), 4, color=C_STEEL))
            for tp in pk.theta:
                a_ = math.radians(tp)
                r_f = float(pk.floor_inner_r(xc))
                ax.plot([ro_t * math.sin(a_), r_f * math.sin(a_)], [ro_t * math.cos(a_), r_f * math.cos(a_)],
                        color="#BDB4A5", lw=2.0)
        ax.set_aspect("equal")
        ax.set_xlim(-42, 42)
        ax.set_ylim(-42, 42)
        ax.set_title(f"x = {xc:.1f}: {lbl}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle("Cross sections looking aft (+Z up, +Y right); pockets at theta 0/120/240", fontsize=8)
    fig.tight_layout()
    p3 = FIG_DIR / "layout_section_cross.png"
    fig.savefig(p3, bbox_inches="tight")
    plt.close(fig)
    out.append(p3)
    return out


# ==========================================================================
# Trade study on the equator position
# ==========================================================================


def equator_trade(pk: Pockets) -> list:
    rows = []
    for x_eq in (40.0, 55.0, 57.5, 60.0):
        J = joint(x_eq, pk)
        row = dict(x_eq_mm=x_eq, spigot_bore_d_mm=2 * J["r_spigot_in"], lip_d_mm=2 * J["r_lip"],
                   seal_bore_d_mm=2 * J["r_bore"])
        if x_eq >= 50:
            bay = generator_bay(J, rail_flush=False)
            row.update(gen_front_shell=bay["pick"], frontier_L26_D=next(f["D_max_mm"] for f in bay["frontier"]
                                                                      if f["L_mm"] == 26.0))
        else:
            row.update(note="joint forward of the pockets: the socket sits on the sloped nose; a "
                            "generator in the aft shell must pass the spigot bore")
        rows.append(row)
    sens = []
    for d in (3.0, 4.0, 5.0, 6.0):
        ri2 = dict(pk.ri, depth=d, depth_x=None, depth_v=None)
        J = joint(P.EQUATOR_X, Pockets(ri2))
        sens.append(dict(pocket_depth_at_equator_mm=d, seal_bore_d_mm=2 * J["r_bore"],
                         spigot_bore_d_mm=2 * J["r_spigot_in"], lip_d_mm=2 * J["r_lip"]))
    return rows, sens


# ==========================================================================
# Main
# ==========================================================================


def run(x_eq: float | None = None, write: bool = True, log=print) -> dict:
    ri = rotor_inputs()
    hw = hardware_inputs()
    pk = Pockets(ri)
    x_eq = P.EQUATOR_X if x_eq is None else x_eq
    J = joint(x_eq, pk)
    bay = generator_bay(J, rail_flush=False)
    bay_flush = generator_bay(J, rail_flush=True)
    bay["frontier_flush"] = bay_flush["frontier"]
    bay["pick_if_flush"] = bay_flush["pick"]
    DT = drivetrain(J, bay, ri, hw, pk)
    x_rotor = DT["stations_mm"]["rotor_load_plane"]
    fins = fin_design(ri, x_rotor)
    nose_d = nose(ri, fins)
    mb = mass_buoyancy(J, bay, DT, fins, ri, hw, pk, nose_d)
    roll = roll_summary(ri, fins, mb)
    trade, sens = equator_trade(pk)
    res = dict(inputs=dict(rotor=ri, hardware=hw), joint=J, generator_bay=bay, drivetrain=DT,
               nose=nose_d, fins=fins, roll=roll, mass=mb, equator_trade=trade,
               pocket_sensitivity=sens, assembly=assembly(J, bay, DT), parts=part_split(J, bay, fins))
    return res


def main() -> None:
    res = run()
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "last_run.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("joint",)}, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
