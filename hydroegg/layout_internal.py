"""Internal layout study for HydroEgg: everything inside the static hull.

Run ``python -m hydroegg.layout_internal`` to regenerate ``data/layout_internal.json`` and
``docs/figures/layout_section*.png``. A colleague owns the rotor side (fold kinematics, blade
pockets, hub: ``hydroegg/layout.py`` -> ``data/layout_rotor.json``); this module owns the
inside: drivetrain, generator bay, equator joint, nose (tow eye + cable gland), roll fins,
mass and buoyancy, assembly sequence and part split.

Frame: X = body axis, flow travels +X, nose at x = 0, rotor at the tail. Z up.
theta = angle about +X measured from +Z toward +Y. Units: mm / deg in geometry, SI in physics.

Architecture chosen here (justified in the JSON ``decisions`` block):
  * Front shell F (nose -> equator): tow eye on a bridged nose plug, M12 cable gland, motor
    carrier (twist-lock bulkhead) with a motor-specific adapter, generator with its rotor
    face aft. The bay straddles the joint: the motor's aft end slides into the spigot bore
    of the aft shell when the shells close. All electrics stay in F.
  * Aft shell A (equator -> static tail end): spigot with O-ring + 3 bayonet lugs, bearing
    tube on 3 ribs under the pockets with a 608 (front, locating, takes thrust) and a 688
    (rear, floating), lip seal at the tail end, 8 mm shaft.
  * A two-dog claw coupling engages blind when the shells close; claw B is the thrust collar.
Rotor inputs come from ``hydroegg.layout`` (the rotor half): pocket depth / floor / keep-out /
half-width and azimuth tables, the solved blade axis (load plane, fold sweep), hub stations,
detent magnets and blade mass. They are read from ``layout.DATA_JSON`` (data/layout_rotor.json,
written by ``python -m hydroegg.layout``); if that file is missing or does not match
``layout.CHOSEN`` the rotor study is re-run through ``layout.run()`` first. data/hardware.json
(catalogue data) and data/blade_design.json (thrust, cone nominal) are read as before.

The merged, consistent layout (rotor + internal + interface checks) is written to
``data/layout.json`` by the same run (see ``merged()``).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import layout as RL          # rotor half: single source for pockets, axis, hub
from . import params as P
from . import profile as prof

ROOT = Path(__file__).resolve().parents[1]
DATA_JSON = ROOT / "data" / "layout_internal.json"
MERGED_JSON = ROOT / "data" / "layout.json"
ROTOR_JSON = RL.DATA_JSON
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
C_SEAL = 0.10               # mm, radial gap spigot <-> seal bore (O-ring extrusion gap); the
                            # hardware study's value: 0.20 diametral > 0.13 table limit, accepted
                            # because the hull sees < 0.1 bar

# Fallback pocket model (task brief) until data/layout_rotor.json exists
POCKET_FALLBACK = dict(x0=45.0, x1=118.0, depth=6.0, width=30.0, theta=(0.0, 120.0, 240.0))

# Equator joint (bayonet + radial O-ring gland on the aft-shell spigot)
# O-ring: ISO 3601 / AS568 1xx series, cs 2.62 (the hardware study's choice: a 2.0 / 1.78
# section is too sensitive to +-0.1 mm print error). Free ID and ID tolerance (mm) from the
# Apple Rubber "ISO 3601 Metric Size O-Rings" chart, cached by the hardware study at
# out/work/hardware/src/iso3601_apple.pdf, parsed 2026-10-03.
AS568_1XX = {"-130": (40.94, 0.38), "-131": (42.52, 0.38), "-132": (44.12, 0.38),
             "-133": (45.69, 0.38), "-134": (47.29, 0.38), "-135": (48.90, 0.43),
             "-136": (50.47, 0.43), "-137": (52.07, 0.43), "-138": (53.64, 0.43),
             "-139": (55.25, 0.43), "-140": (56.82, 0.43), "-141": (58.42, 0.51),
             "-142": (59.99, 0.51), "-143": (61.60, 0.51), "-144": (63.17, 0.51),
             "-145": (64.77, 0.51)}
OR_CS = 2.62                      # mm (params.ORING_EQUATOR cs 2.0 -> proposed 2.62)
OR_CS_TOL = 0.08
# Static radial gland, cs 2.62, no back-up ring (same Apple Rubber chart, via data/hardware.json)
OR_DEPTH = 0.5 * (2.03 + 2.18)    # mm, radial gland depth (table 2.03-2.18)
OR_TABLE_W = 3.56                 # mm, rectangular groove width in the table
OR_TABLE_SQUEEZE = (14, 25)       # %
OR_TABLE_DCLR = 0.13              # mm, max diametral clearance (extrusion guard)
OR_STRETCH = (0.02, 0.04)         # nominal ID stretch window onto the groove bottom (max 5 %)
OR_GROOVE_W = 2.9                 # mm, groove BOTTOM width; forward flank at 45 deg so the
                                  # spigot prints tip-up without support (mean width ~3.95)
OR_LAND_ROOT = 1.5                # mm, land between equator face and groove
OR_LAND_TIP = 1.5                 # mm, land forward of the groove
TOL_PRINT_D = 0.10                # mm, +- on printed diameters after calibration (hardware study)
LUG_N = 3
LUG_H = 1.5                       # mm, radial height
LUG_T = 3.0                       # mm, axial thickness (1.5 straight + 1.5 at 45 deg)
LUG_W_DEG = 16.0                  # angular width (narrow enough that the twist channel clears
                                  # the pockets: 3 pockets and 3 lugs share 120 deg sectors)
GAP_W_DEG = 20.0                  # entry gap in the socket lip (lug + 2 deg each side)
ENGAGE_DEG = 35.0                 # twist from entry to lock
LIP_T = 3.0                       # mm, socket lip axial thickness at its ID
CH_EXTRA = 0.4                    # mm, axial running clearance of the lug in its channel

# Drivetrain stations
SEAL_RECESS = 0.5                 # mm, lip seal outer face behind the static-tail end face
SEAL_SHOULDER_T = 1.5             # mm, ring between seal bore and rear bearing seat
SEAL_SHOULDER_ID = 12.5           # mm, supports 688 outer ring and seal back face
FRONT_SHOULDER_ID = 18.0          # mm, supports 608 outer ring, lets a 688 pass through
CLAW_A_L = 5.0                    # mm, claw A (on the generator bell / hollow shaft)
CLAW_FLOAT = 1.0                  # mm, axial float between the claws
CLAW_B_L = 10.5                   # mm, 0.5 shim + claw B: printed sleeve pressed over a Ruland
                                  # MSC-8-SS set-screw collar (8 x 16 x 8, hardware.json) = thrust collar
COUPLER_STACK = CLAW_A_L + CLAW_FLOAT + CLAW_B_L   # bell face -> 608 front face
COUPLER_OD = 22.0
CLAW_B_STEP_L = 2.0               # mm, claw B steps down over its last 2 mm (ahead of the shim) so the
CLAW_B_STEP_OD = 20.0             # rotating claw clears the static tube mouth (608 seat bore 22) radially
SHIM_OD = 12.0                    # mm, shim washer on the 608 inner race (inner-race land, clears the seals)
TUBE_LEAD_IN = 0.5                # mm x 45 deg chamfer at the tube mouth (608 press-in lead)
SHAFT_AFT_FALLBACK = 146.0        # mm, shaft end inside hub/spinner until layout_rotor.json says
CARRIER_T = 4.0                   # mm, motor carrier plate
CARRIER_RING = 1.6                # mm, inward ring that holds the carrier tab slots
CARRIER_EDGE = 1.0                # mm, plate material outside the screw counterbores
# Shaft stability criteria (choices). Loads are tiny (bearing static safety > 15), so geometry
# governs: a minimum span, the blade-tip axial runout produced by printed-seat slop, the static
# safety of the bearings and the shaft deflection at the rotor. span / overhang is reported
# against the 1.25 rule of thumb for overhung rotors; with the static tail at x 108 and the
# coned rotor (load plane ~x 136) no bay position meets it, so the explicit checks replace it
# (merge decision, see docs/notes/layout.md).
SPAN_RATIO_REF = 1.25             # rule of thumb, reported only
SPAN_MIN = 4 * P.SHAFT_D          # mm
SEAT_SLOP = 0.05                  # mm radial, each printed bearing seat (ESTIMATE)
TIP_RUNOUT_MAX = 0.5              # mm at the blade tip
S0_MIN = 4.0                      # min static safety C0 / P0 (choice, generous for 2RS balls)
SHAFT_DEFL_MAX = 0.05             # mm at the rotor load plane (radial load + asymmetry moment)
E_STEEL = 200e3                   # MPa, ESTIMATE for 440C / 304 stainless
R_THRUST_FRAC = 0.7               # ESTIMATE, thrust centroid radius / tip radius
BAY_L_MARGIN = 1.0                # mm over the longer of the recommended motor and its alternative

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

# Cable gland (M12x1.5, threaded boss). Hummel HSK-K M12 / Lapp Skintop ST M12 (data/
# hardware.json): thread length 8, hex SW 15 (across corners 17.3), height 21, cable 3-6.5 mm.
# The boss sits on the nose beside the tow eye; its place and axis come from gland_place().
GLAND_THETA = 180.0               # deg, bottom of the nose
GLAND_SEAT_D = 18.0               # mm, flat seat (spot face) for the hex and washer: 17.3 across corners + 0.7
GLAND_SEAT_T = P.HULL_WALL        # mm, depth of the seat region (the shell wall under the spot face)
GLAND_BOSS_D = 12.0 + 2 * T_MIN   # mm, thread boss inside the shell: M12 + 2 x T_MIN wall
GLAND_ENGAGE = 8.0                # mm, thread length of the gland (fallback; hardware.json wins)
GLAND_HOLE = 12.5                 # mm, clearance hole when not threaded (hardware.json)
GLAND_CLEAR = 1.0                 # mm, boss to the carrier plate and to the carrier ring
GLAND_TOW_CLEAR = 2.0             # mm, gland body outside the hull to the tow eye and the shackle envelope
GLAND_X_SEARCH = (3.0, 12.0, 0.25)       # mm, seat centre x on the surface (from, to, step)
GLAND_TILT_SEARCH = (-10.0, 30.0, 2.5)   # deg, outward axis turned from the surface normal toward -X
CARRIER_RING_OVER = 1.2           # mm, the shell's carrier ring runs this far ahead of / behind the plate
CABLE_CLEAR = 1.0                 # mm, cable notch in the carrier: cable OD + this, both ways
# Tow eye orientation and shackle envelope. The eye ring lies in the XZ plane (hole axis Y), so a
# shackle pin runs along Y and the shackle body lies in the XY plane; the gland below (theta 180)
# stays out of it. Shackle: ESTIMATE for a 5 mm stainless D-shackle (pin 5, inner width 10):
# from the hole's aft edge forward, |pin axis| <= 12, |across| <= 5.
EYE_PLANE = "XZ"
SHACKLE_HALF_PIN = 12.0           # mm, ESTIMATE
SHACKLE_HALF_ACROSS = 5.0         # mm, ESTIMATE

# Fins
FIN_SECTION = "clarky"            # flat lower surface: prints flat with no support
FIN_N_CRIT = 5.0                  # disturbed river inflow, same as the BEM study
FIN_CP_SPAN = 0.44                # ESTIMATE, spanwise centre of lift / span (rectangular)
FIN_E_OSWALD = 0.8                # ESTIMATE
FIN_CLMAX_3D_FACTOR = 0.9         # ESTIMATE, 3D CLmax / 2D CLmax
FIN_CL_RESERVE = 0.75             # design CL <= this x CLmax_3D
FIN_DESIGN_V = 1.5                # m/s, the fins are matched at the design flow (params.DESIGN_FLOW_V)
FIN_SPEEDS = (1.0, 1.5, 2.5)      # m/s, flows at which the fin couple is evaluated (2.5 = tow case)
FIN_TORQUE_MARGIN = 0.10          # fin couple at the design flow >= (1 + margin) x the recommended
                                  # generator's load-line torque (data/blade_design.json): covers lift
                                  # that a printed Clark Y at Re ~5e4 may not reach (NeuralFoil is
                                  # smooth-surface) and the estimated GL40 drag torque; trimmed later
                                  # with fin variants (fin_design: trim)
FIN_RATIO_BAND = 0.05             # and <= (1 + margin)(1 + band): no fin that overshoots at zero incidence
FIN_ALPHA_STEP = 0.25             # deg, incidence grid
FIN_TRIM_DEG = 1.0                # deg, trim fin variant: root tongue turned by this about the pivot
FIN_TIP_CLEAR = 5.0               # mm, margin to the blade-tip sweep torus ("torus" rule)
FIN_FOLD_RULE = "index"           # merge decision: the blades fold only at the stowing index (the
                                  # pockets require it), checked against the swept blade surface;
                                  # "torus" = any rotor azimuth (TE outside the fold torus), kept as
                                  # the comparison case
FIN_TE_EQ_MARGIN = 2.0            # mm, rail end stop ahead of the equator (fins on the front shell)
FIN_THETA_SEARCH = (70.0, 110.0)  # deg, +side fin angle search; the other fin sits at +180
FIN_FOLD_CLEAR = 2.0              # mm, folding blade (stowing index) to the fin and its rail pad
FIN_POCKET_CLEAR = 2.0            # mm, rail pad edge to a stowed-blade pocket edge (arc)
FIN_CHORDS = (25.0, 30.0, 35.0, 40.0, 45.0)   # chords whose LE would fall ahead of x 4 are skipped
FIN_THETA_STEP = 0.5              # deg, fin angle search step
FIN_PIVOT_STEP = 0.05             # yaw pivot search step (chord fraction from the LE)
FIN_PLACE_MAX_EVAL = 400          # cap on fold-clearance evaluations per fin candidate
FIN_T_HALF = 0.5 * 0.117          # fin half-thickness envelope / chord (Clark Y max t/c)
GROOVE_W = 6.0                    # mm, mean width of the half-dovetail groove in the 10 mm pad (>= 2 mm
                                  # walls); the fin's root tongue fills it (ESTIMATE until the CAD)
RAIL_H = 4.0                      # mm, min pad height above the hull = dovetail depth
RAIL_W = 10.0                     # mm, pad width (tangential)
RAIL_X = (20.0, 43.0)             # mm, nominal rail extent (replaced by the fin root chord)
BRIDLE_A = 20.0                   # mm, v2 option: bridle eyes at +-Y on the nose
BRIDLE_L = 80.0                   # mm, bridle leg length

# Mass model
RHO_PLA = 1.24                    # g/cc, Prusament TDS
RHO_PETG = 1.27                   # g/cc, Prusament TDS
SKIN = 1.2                        # mm, 3 walls ~ 3 x 0.4 (also used for top/bottom shells)
INFILL = 0.15                     # gyroid
RHO_STEEL = 7.75                  # g/cc, ESTIMATE AISI 440C (hardware.json shaft) ~7.7-7.8
RHO_LEAD = 11.34                  # g/cc, ESTIMATE
RHO_NBR = 1.3                     # g/cc, ESTIMATE
RHO_NDFEB = 7.5                   # g/cc, ESTIMATE
RHO_FOAM = 0.03                   # g/cc, ESTIMATE closed-cell foam
SPINNER_ACCESS_D = 6.0            # mm, screw-access bore on the axis of the spinner (M3 socket head 5.5 + 0.5)
SPINNER_AFT_WALL = 3.0            # mm, spinner wall clamped by the M3 screw against the hub aft face
M3_HEAD_H = 3.0                   # mm, ISO 4762 M3 head height

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


def _num(v):
    """Scalar from a number or a {'value': ...} record (data/hardware.json style)."""
    if isinstance(v, dict) and "value" in v:
        v = v["value"]
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    return float(v)


def _arr(seq) -> np.ndarray:
    return np.array([float(v) if isinstance(v, (int, float)) and not isinstance(v, bool)
                     else np.nan for v in seq], dtype=float)


def _rotor_json(allow_live: bool = True) -> dict | None:
    """data/layout_rotor.json as written by hydroegg.layout. If it is missing or its chosen
    hinge / static tail / fold differ from layout.CHOSEN, re-run the rotor study through
    layout.run() (about 70 s) so both halves always describe the same rotor."""
    rj = _load(ROTOR_JSON)
    ch = (rj or {}).get("chosen", {})
    fresh = bool(rj) and all(
        abs(float(ch.get(f"{k}_mm" if k != "fold_deg" else k, math.nan)) - float(RL.CHOSEN[k])) < 1e-6
        for k in ("hinge_x", "hinge_r", "static_tail_x", "fold_deg"))
    if not fresh and allow_live:
        print("layout_rotor.json missing or stale: running hydroegg.layout.run()", flush=True)
        RL.run(write=True, figures=False)
        rj = _load(ROTOR_JSON)
    return rj


def _rotor_tables(rj: dict) -> dict:
    """The rotor half's tables this module uses, as numpy arrays (radial, mm; azimuths in the
    rotor's psi frame relative to each blade)."""
    pd = rj["pocket_depth"]
    tab = {k: _arr(pd[k]) for k in ("x_mm", "depth_mm", "floor_r_mm", "keepout_r_mm", "half_width_mm",
                                     "half_width_deg", "az_center_deg", "az_lo_deg", "az_hi_deg",
                                     "land_between_pockets_mm") if k in pd}
    ax = rj.get("axis", {})
    axis = {k: _arr(ax[k]) for k in ("r_mm", "s_mm", "lead_mm", "aft_mm") if k in ax}
    ss = rj.get("stowed_sections", {})
    sec = {k: _arr(ss[k]) for k in ("r_mm", "chord_mm", "twist_deg", "thickness_mm") if k in ss}
    return dict(pocket=tab, axis=axis, sections=sec)


def rotor_inputs() -> dict:
    """Pocket model and rotor geometry from the rotor half (hydroegg.layout), plus
    data/blade_design.json for the BEM thrust and the nominal cone."""
    out = dict(source="fallback (task brief: 6 mm deep, 30 mm wide, x 45-118, 3 x 120 deg)",
               **POCKET_FALLBACK, depth_x=None, depth_v=None, width_v=None, az_lo_v=None, az_hi_v=None,
               floor_r_v=None, keepout_r_v=None,
               static_tail_x=P.STATIC_TAIL_X, hinge_x=P.HINGE_X, hinge_r=P.HINGE_R,
               shaft_aft_x=None, r_tip=P.ROTOR_R_TIP_TARGET, cone_deg=12.0,
               thrust_bem_2p5_N=None, blade_mass_solid_PLA_g=None, omega_vector="+X",
               blade_source="params (data/blade_design.json absent)", x_load=None,
               sweep_len=None, rotor=None)
    bd = _load(BLADE_JSON)
    if bd:
        out["blade_source"] = "data/blade_design.json"
        for key, paths in (("r_tip", ["rotor.r_tip_mm"]), ("cone_deg", ["rotor.cone_deg"]),
                           ("thrust_bem_2p5_N", ["structure_estimates.max_thrust_2p5_m_s_N"]),
                           ("blade_mass_solid_PLA_g", ["structure_estimates.blade_mass_g"])):
            v = _num(_find(bd, paths))
            if v is not None:
                out[key] = v
        om = _find(bd, ["rotor.omega_vector"])
        if isinstance(om, str):
            out["omega_vector"] = om
    rj = _rotor_json()
    if rj and "pocket_depth" in rj:
        ch = rj["chosen"]
        T = _rotor_tables(rj)
        tp, ax, sec = T["pocket"], T["axis"], T["sections"]
        x = tp["x_mm"]
        out.update(
            source="hydroegg.layout -> data/layout_rotor.json (pocket_depth table, asymmetric about "
                   "each blade azimuth)",
            static_tail_x=float(ch["static_tail_x_mm"]), hinge_x=float(ch["hinge_x_mm"]),
            hinge_r=float(ch["hinge_r_mm"]), shaft_aft_x=float(rj.get("shaft_aft_x_mm", ch["shaft_end_x_mm"])),
            r_tip=float(ch["rotor_r_tip_mm"]),
            depth_x=x.tolist(), depth_v=np.nan_to_num(tp["depth_mm"]).tolist(),
            width_v=np.nan_to_num(2 * tp["half_width_mm"]).tolist(),
            az_lo_v=tp["az_lo_deg"].tolist(), az_hi_v=tp["az_hi_deg"].tolist(),
            floor_r_v=tp["floor_r_mm"].tolist(), keepout_r_v=tp["keepout_r_mm"].tolist(),
            x0=float(x.min()), x1=float(ch["static_tail_x_mm"]),
            depth=float(np.nanmax(tp["depth_mm"])), width=float(2 * np.nanmax(tp["half_width_mm"])),
            # blade i sits at psi = 120 i; theta = -psi, so the blade (and pocket) azimuths are
            # theta 0 / 240 / 120; the pocket spans theta_k + [-az_hi, -az_lo] at each x
            theta=(0.0, 120.0, 240.0),
            cone_deg=float(rj["axis"].get("cone_mean_hinge_to_tip_deg", out["cone_deg"])),
            floor_overhang_max_deg=rj["pocket_summary"].get("floor_overhang_max_deg"),
            floor_print_rule=rj["pocket_summary"].get("floor_print_rule"),
            blade_mass_solid_PLA_g=float(rj["blade_mass"]["mass_g"]),
            blade_cg_n_mm=float(rj["blade_mass"]["cg_blade_frame_mm"][2]),
            rotor=dict(
                hub_front_x=float(ch["hub_front_x_mm"]), pin_d=float(ch["pin_d_mm"]),
                pin_l=float(ch["pin_l_mm"]), hub_aft_x=float(rj["hub"]["hub_aft_face_x_mm"]),
                grub_x=float(rj["hub"]["grub"]["x_mm"]),
                stop_block_x=[float(v) for v in rj["hub"]["stop"]["stop_block_x_mm"]],
                stow_magnet_x=float(rj["detents"]["stowed"]["x_mm"]),
                stow_magnet_r=float(rj["detents"]["stowed"]["static_magnet_centre_r_mm"]),
                stow_magnet_keepout_r=float(rj["detents"]["stowed"]["static_keepout_r_mm"]),
                stow_boss_h=float(rj["detents"]["stowed"].get("boss_height_over_pocket_floor_mm", float("nan"))),
                n_detent_magnets=12, keepout_rule="floor - MIN_WALL x hull normal factor",
                notch_floor_r=float(rj["pocket_summary"]["notch_at_static_tail_end"]["floor_r_mm"]),
                keepout_tail_min_r=float(rj["pocket_summary"]["notch_at_static_tail_end"]
                                         ["keepout_last_12mm_min_mm"])),
        )
        if ax:
            # thrust load plane on the solved (curved) axis: deployed, the blade frame's aft
            # offset is +X, so a station at radius r sits at x = hinge_x + aft(r)
            r_c = R_THRUST_FRAC * out["r_tip"]
            out["x_load"] = out["hinge_x"] + float(np.interp(r_c, ax["r_mm"], ax["aft_mm"]))
            out["x_load_rule"] = ("hinge_x + aft(r) of the solved blade axis at r = %.2f R (thrust "
                                  "centroid radius, ESTIMATE)" % R_THRUST_FRAC)
            # fold sweep: largest distance of the blade from the pin axis (the pin runs along the
            # chord, so only span s and the thickness direction count): pitch-axis point plus half
            # the section thickness and the twisted chord's out-of-plane half extent
            if sec:
                t_half = 0.5 * np.interp(ax["r_mm"], sec["r_mm"], sec["thickness_mm"])
                c_out = 0.5 * np.interp(ax["r_mm"], sec["r_mm"], sec["chord_mm"]) * np.abs(
                    np.sin(np.radians(np.interp(ax["r_mm"], sec["r_mm"], sec["twist_deg"]))))
            else:
                t_half = c_out = 0.0
            out["sweep_len"] = float(np.max(np.hypot(ax["s_mm"], np.abs(ax["aft_mm"]) + t_half + c_out)))
    if out["shaft_aft_x"] is None:
        out["shaft_aft_x"] = out["static_tail_x"] + (SHAFT_AFT_FALLBACK - P.STATIC_TAIL_X)
        out["shaft_aft_rule"] = "fallback: 28 mm aft of the static tail end (hub + spinner)"
    return out


def hardware_inputs() -> dict:
    """Catalogue data. data/hardware.json wins when it carries the field; the defaults are
    the values this study sourced itself (out/work/layout_internal/sources/sources.json)."""
    hw = dict(
        b608=dict(P.BEARING, C_N=3300.0, C0_N=1350.0, mass_g=11.9,
                  src="vxb listing (C, C0), SKF 608-2RSH page (mass), via web search 2026-10-03"),
        b688=dict(P.BEARING_ALT, C_N=1600.0, C0_N=710.0, mass_g=3.8,
                  src="vxb 688-2RS page (C, C0); mass from a web-search summary, VERIFY"),
        seal=dict(P.LIP_SEAL, mass_g=2.0, src="mass ESTIMATE"),
        motor=dict(name="iPower GM4108H-120T", D=47.0, L=22.7, D_reserve=47.0, L_reserve=33.0,
                   mass_g=124.0, pcd=38.14, rotor_pcd=None, n_screw=4, screw="M3", cbore_d=6.6,
                   alt_name=None, alt_pcd=None, centre_bore=8.0, src="defaults"),
        gland=dict(name="Hummel HSK-K M12x1.5 (1.209.1200.50)", thread_L=8.0, SW=15.0,
                   height=21.0, clamp_mm=[3.0, 6.5], mass_g=5.0, src="defaults; mass ESTIMATE"),
        envelope=dict(D_max=56.0, L_max=35.0), candidates=[],
        source="defaults (out/work/layout_internal/sources/sources.json)",
    )
    hj = _load(HW_JSON)
    if not hj:
        return hw
    hw["source"] = "data/hardware.json (+ defaults for fields it does not carry)"
    s6 = _find(hj, ["bearings.b608_2rs.stainless_S608_2RS"]) or {}
    if _num(s6.get("C_dyn")) and _num(s6.get("C0")):
        hw["b608"].update(C_N=_num(s6["C_dyn"]), C0_N=_num(s6["C0"]),
                          src="data/hardware.json (stainless S608-2RS, bearingsdirect); mass SKF")
    b8 = _find(hj, ["bearings.b688_2rs"]) or {}
    if _num(b8.get("C_dyn")) and _num(b8.get("C0")):
        hw["b688"].update(C_N=_num(b8["C_dyn"]), C0_N=_num(b8["C0"]),
                          src="data/hardware.json (bearingsdirect 688-2RS); mass web summary, VERIFY")
    cad = _find(hj, ["generator.recommendation.cad_interface"]) or {}
    mo = hw["motor"]
    for k_out, k_in in (("D", "body_D_mm"), ("L", "body_L_mm"), ("D_reserve", "bay_D_reserve_mm"),
                        ("L_reserve", "bay_L_reserve_mm"), ("centre_bore", "shaft_bore_mm")):
        v = _num(cad.get(k_in))
        if v is not None:
            mo[k_out] = v
    sf = _find(cad, ["stator_face_holes.value"]) or _find(cad, ["base_holes.value"]) or {}
    rf = _find(cad, ["rotor_face_holes.value"]) or {}
    if _num(sf.get("pcd")):
        mo["pcd"] = _num(sf["pcd"])
        mo["n_screw"] = int(_num(sf.get("n")) or 4)
    if _num(rf.get("pcd")):
        mo["rotor_pcd"] = _num(rf["pcd"])
    choice = _find(hj, ["generator.recommendation.choice"])
    if isinstance(choice, str):
        mo["name"] = choice
    alt = _find(hj, ["generator.recommendation.alternative"])
    if isinstance(alt, str):
        mo["alt_name"] = alt.split(" (")[0]
        # the alternative's base pattern (4 x M3 on PCD 38.14, community STEP, VERIFY)
        mo["alt_pcd"] = 38.14 if "GM4108" in alt else None
    mo["src"] = "data/hardware.json generator.recommendation.cad_interface"
    env = _find(hj, ["generator.envelope"]) or {}
    hw["envelope"] = dict(D_max=_num(env.get("D_max_mm")) or 56.0, L_max=_num(env.get("L_max_mm")) or 35.0)
    for c in _find(hj, ["generator.candidates"]) or []:
        D, L, m = _num(c.get("D")), _num(c.get("L")), _num(c.get("mass"))
        hw["candidates"].append(dict(id=c.get("id"), name=c.get("name"), D_mm=D, L_mm=L, mass_g=m))
        if c.get("name") == mo["name"] and m:
            mo["mass_g"] = m
    gl = _find(hj, ["cable_gland.hummel_hsk_k_m12.value"]) or {}
    if gl:
        hw["gland"].update(thread_L=_num(gl.get("thread_length")) or 8.0, SW=_num(gl.get("SW")) or 15.0,
                           height=_num(gl.get("height")) or 21.0, clamp_mm=gl.get("clamping", [3.0, 6.5]),
                           src="data/hardware.json cable_gland.hummel_hsk_k_m12; mass ESTIMATE")
    # equator O-ring: the hardware study's final pick (size, groove form and diameters)
    orj = hj.get("oring_equator") or {}
    rec = (orj.get("recommended") or {}).get("value") or {}
    gd = (orj.get("gland_design") or {}).get("value") or {}
    if rec and gd:
        hw["oring"] = dict(
            designation=rec["designation"], dash="-" + rec["designation"].split("-")[1].split()[0],
            id=_num(rec["ID"]), id_tol=_num(rec["ID_tol"]), cs=_num(rec["cs"]), cs_tol=_num(rec["cs_tol"]),
            material=rec.get("material"), groove_form=gd.get("groove_form"), flank=gd.get("flank"),
            groove_bottom_d=_num(gd["spigot_groove_bottom_D"]), groove_bottom_w=_num(gd["groove_bottom_width"]),
            groove_top_w=_num(gd["groove_top_width"]), gland_depth=_num(gd["gland_depth"]),
            groove_depth_in_spigot=_num(gd["groove_depth_in_spigot"]), bore_d=_num(gd["bore_D"]),
            spigot_od=_num(gd["spigot_OD"]), radial_clearance=_num(gd["radial_clearance"]),
            performance=(orj.get("performance") or {}).get("value") or {},
            src="data/hardware.json oring_equator (recommended, gland_design, performance)")
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
    """Blade pockets in the static hull: depth(x) (radial, deepest point over azimuth, from the
    hull surface), width(x) (arc at the hull radius), the blade azimuths theta and, with the
    rotor table, the pocket's angular extent about each blade azimuth (asymmetric: the stowed
    blade's lead and twist push it toward the direction of rotation)."""

    def __init__(self, ri: dict):
        self.ri = ri
        self.x0, self.x1 = ri["x0"], ri["x1"]
        self.width = ri["width"]
        self.theta = tuple(ri["theta"])

    def depth(self, x):
        x = np.asarray(x, dtype=float)
        if self.ri.get("depth_x"):
            return np.interp(x, self.ri["depth_x"], self.ri["depth_v"], left=0.0, right=0.0)
        return np.where((x >= self.x0) & (x <= self.x1), self.ri["depth"], 0.0)

    def width_at(self, x):
        x = np.asarray(x, dtype=float)
        if self.ri.get("width_v"):
            return np.interp(x, self.ri["depth_x"], self.ri["width_v"], left=0.0, right=0.0)
        return np.where(self.depth(x) > 0, self.width, 0.0)

    def theta_span(self, x):
        """(lo, hi) theta offsets (deg) of the pocket relative to its blade azimuth at x; nan
        where there is no pocket. theta = -psi, so the rotor's [az_lo, az_hi] maps to
        [-az_hi, -az_lo]."""
        x = np.asarray(x, dtype=float)
        has = self.depth(x) > 0
        if self.ri.get("az_lo_v") and self.ri.get("depth_x"):
            xs = self.ri["depth_x"]
            lo = -np.interp(x, xs, self.ri["az_hi_v"])
            hi = -np.interp(x, xs, self.ri["az_lo_v"])
        else:
            h = self.half_angle_deg(x)
            lo, hi = -h, h
        return np.where(has, lo, np.nan), np.where(has, hi, np.nan)

    def edges_deg(self, x):
        """List of (lo, hi) absolute theta intervals, one per pocket, at x (scalar)."""
        lo, hi = self.theta_span(float(x))
        if not np.isfinite(lo):
            return []
        return [((tp + float(lo)) % 360.0, (tp + float(hi)) % 360.0) for tp in self.theta]

    def keepout_r(self, x):
        """The rotor study's keep-out radius under a pocket (floor - MIN_WALL), nan elsewhere."""
        x = np.asarray(x, dtype=float)
        if not (self.ri.get("keepout_r_v") and self.ri.get("depth_x")):
            return np.where(self.depth(x) > 0, r_out(x) - self.depth(x) - P.MIN_WALL, np.nan)
        k = np.interp(x, self.ri["depth_x"], self.ri["keepout_r_v"])
        return np.where(self.depth(x) > 0, k, np.nan)

    def floor_inner_r(self, x):
        """Radius of the shell's inner surface under a pocket (where the pocket exists)."""
        x = np.asarray(x, dtype=float)
        d = self.depth(x)
        rin_pk = r_out(x) - d - POCKET_FLOOR_WALL * np.sqrt(1 + slope(x) ** 2)
        return np.where(d > 0, np.minimum(rin_pk, r_in(x)), r_in(x))

    def half_angle_deg(self, x):
        return np.degrees(self.width_at(x) / 2 / np.maximum(r_out(x), 1e-6))


# ==========================================================================
# Equator joint
# ==========================================================================


def _stretched_cs(cs: float, id_free: float, d_groove: float) -> float:
    """Volume-conserving cross-section after stretching the ring onto d_groove (same model
    as the hardware study)."""
    dm0 = id_free + cs
    cs1 = cs
    for _ in range(30):
        cs1 = cs * math.sqrt(dm0 / (d_groove + cs1))
    return cs1


def _gland_area(w_bot: float, local_depth: float, groove_depth: float) -> float:
    """Gland cross-section between the groove bottom and the bore over a local radial depth, for a
    groove with one 45 deg (tip-side) flank: the width grows 1:1 with height inside the groove (up
    to groove_depth, the depth cut into the spigot), then stays at the top width in the radial
    clearance above it (same model as the hardware study, hw_calc.gland_area)."""
    h_in = min(local_depth, groove_depth)
    return w_bot * local_depth + 0.5 * h_in * h_in + h_in * max(local_depth - groove_depth, 0.0)


def oring_pick(d_bore_max: float, hw_or: dict | None = None) -> dict:
    """Equator O-ring and spigot groove. With hw_or (data/hardware.json oring_equator) the
    hardware study's final pick is adopted as is (size, groove bottom and bore diameters, bottom
    width, 45 deg tip-side flank, radial clearance) and only checked here: the bore must fit under
    d_bore_max (POCKET_FLOOR_WALL under the pockets over the seal zone). Without it, the largest
    AS568 1xx ring whose bore fits is picked at OR_STRETCH. Either way the gland is re-evaluated with
    the hardware study's worst-case model: print +-TOL_PRINT_D on diameters, full eccentricity in
    the radial clearance, ring cs and ID tolerances in the adverse direction (thinnest section =
    smallest cs and smallest free ID, i.e. most stretch), groove bottom width -0.13."""
    if hw_or:
        dash, idf, id_tol = hw_or["dash"], hw_or["id"], hw_or["id_tol"]
        cs, cs_tol = hw_or["cs"], hw_or["cs_tol"]
        d_gr, d_bore, c_rad = hw_or["groove_bottom_d"], hw_or["bore_d"], hw_or["radial_clearance"]
        w_bot, depth = hw_or["groove_bottom_w"], hw_or["gland_depth"]
        source = hw_or["src"]
    else:
        best = None
        for dash_, (idf_, _) in AS568_1XX.items():
            st = min(OR_STRETCH[1] - 0.005, (d_bore_max - 2 * OR_DEPTH) / idf_ - 1)   # 3.5 % cap
            if st >= OR_STRETCH[0] and (best is None or idf_ * (1 + st) > best[1]):
                best = (dash_, idf_ * (1 + st))
        dash, d_gr = best
        idf, id_tol = AS568_1XX[dash]
        cs, cs_tol = OR_CS, OR_CS_TOL
        d_gr = math.floor(d_gr * 20) / 20
        depth, c_rad, w_bot = OR_DEPTH, C_SEAL, OR_GROOVE_W
        d_bore = d_gr + 2 * depth
        source = "layout_internal fallback pick (data/hardware.json has no oring_equator)"
    g = depth - c_rad                                     # groove depth cut into the spigot
    cs_n = _stretched_cs(cs, idf, d_gr)
    L_min = depth - TOL_PRINT_D - c_rad
    L_max = depth + TOL_PRINT_D + c_rad
    cs_min = _stretched_cs(cs - cs_tol, idf - id_tol, d_gr + TOL_PRINT_D)
    cs_max = _stretched_cs(cs + cs_tol, idf + id_tol, d_gr - TOL_PRINT_D)
    a_max = math.pi / 4 * cs_max ** 2
    out = dict(
        designation=f"ISO 3601 / AS568 {dash}", material=(hw_or or {}).get("material") or "EPDM 70 (water) or NBR 70",
        source=source, cs_mm=cs, id_free_mm=idf, id_tol_mm=id_tol,
        groove_bottom_d_mm=d_gr, bore_d_mm=d_bore, spigot_seal_od_mm=d_bore - 2 * c_rad,
        gland_depth_mm=depth, groove_depth_in_spigot_mm=g, groove_bottom_w_mm=w_bot,
        groove_top_w_mm=w_bot + g, groove_forward_flank_deg=45.0, groove_aft_flank_deg=90.0,
        groove_form="aft flank radial; tip-side (forward) flank 45 deg from the groove bottom up to the spigot "
                    "OD, so the spigot prints tip-up without support",
        radial_clearance_mm=c_rad, table_diametral_clearance_max_mm=OR_TABLE_DCLR,
        table_groove_w_rect_mm=OR_TABLE_W, table_squeeze_pct=list(OR_TABLE_SQUEEZE),
        bore_limit_d_mm=d_bore_max, bore_fits=bool(d_bore <= d_bore_max + 1e-9),
        stretch_nom_pct=100 * (d_gr / idf - 1),
        stretch_range_pct=[100 * ((d_gr - TOL_PRINT_D) / (idf + id_tol) - 1),
                           100 * ((d_gr + TOL_PRINT_D) / (idf - id_tol) - 1)],
        squeeze_nom_pct=100 * (cs_n - depth) / cs_n,
        squeeze_worst_pct=[100 * (cs_min - L_max) / cs_min, 100 * (cs_max - L_min) / cs_max],
        fill_nom_pct=100 * (math.pi / 4 * cs_n ** 2) / _gland_area(w_bot, depth, g),
        fill_worst_no_ecc_pct=100 * a_max / _gland_area(w_bot - 0.13, depth - TOL_PRINT_D, g - TOL_PRINT_D),
        fill_worst_pct=100 * a_max / _gland_area(w_bot - 0.13, L_min, L_min),
        guidance="Apple Rubber static radial table for cs 2.62: depth 2.03-2.18, width 3.56, "
                 "squeeze 14-25 %, diametral clearance <= 0.13; globaloring: fill 60-85 %, "
                 "stretch <= 5 %. Worst cases include +-0.10 mm print error on each diameter "
                 "and full eccentricity. 0.20 mm diametral clearance exceeds the 0.13 table "
                 "limit, which guards extrusion under pressure: acceptable at < 0.1 bar",
    )
    perf = (hw_or or {}).get("performance") or {}
    if perf:
        out["hardware_study_performance"] = dict(
            squeeze_nom_pct=100 * perf["squeeze_nominal"], squeeze_worst_pct=[100 * v for v in perf["squeeze_worst"]],
            fill_nom_pct=100 * perf["fill_nominal"], fill_worst_pct=100 * perf["fill_worst"],
            note="data/hardware.json oring_equator.performance; the values above are recomputed here with "
                 "the same model and must agree")
    return out


def joint(x_eq: float, pk: Pockets, hw_or: dict | None = None) -> dict:
    """Radial stack and axial stations of the bayonet + O-ring joint.

    Spigot on the aft shell (printed tail-down, the spigot is its top); socket in the
    front shell (printed equator-down). Seal zone at the spigot root, so the bayonet
    sits on the dry side; lugs at the spigot tip pass through the seal bore. hw_or: the
    hardware study's O-ring and groove (adopted when its bore fits the walls, else the
    fallback pick of oring_pick, e.g. in the equator trade).
    """
    g_sp = (hw_or["groove_depth_in_spigot"] if hw_or else OR_DEPTH - C_SEAL)
    w_bot = hw_or["groove_bottom_w"] if hw_or else OR_GROOVE_W
    groove_top = w_bot + g_sp                            # 45 deg forward flank over the spigot groove depth
    L_seal = OR_LAND_ROOT + groove_top + OR_LAND_TIP
    x_seal0 = x_eq - L_seal
    x_lip0 = x_seal0 - LIP_T
    x_tip = x_lip0 - LUG_T                               # spigot tip
    xs = np.linspace(x_tip - 2.0, x_eq, 60)
    ro = r_out(xs)
    rb_lim = np.minimum(ro - P.HULL_WALL, ro - pk.depth(xs) - POCKET_FLOOR_WALL)
    xs_seal = xs >= x_seal0
    d_lim = 2 * float(rb_lim[xs_seal].min())
    oring = oring_pick(d_lim, hw_or)
    if hw_or and not oring["bore_fits"]:
        return joint(x_eq, pk, None)
    r_b = oring["bore_d_mm"] / 2                        # seal bore radius (socket)
    r_gr = oring["groove_bottom_d_mm"] / 2
    r_s2 = r_b - oring["radial_clearance_mm"]
    r_lug = r_b - P.CLEAR_SLIDE                          # lugs pass through the seal bore
    r_s1 = r_lug - LUG_H
    r_lip = r_s1 + P.CLEAR_SLIDE
    r_ch = r_lug + P.CLEAR_SLIDE
    r_sp_in = min(r_gr - T_MIN, r_s1 - T_MIN)
    x_ch0 = x_tip - CH_EXTRA                             # channel forward wall (at r_ch)
    x_ring0 = x_ch0 - (r_ch - r_lip)                     # 45 deg ceiling closes at r_lip
    groove = dict(x_aft=x_eq - OR_LAND_ROOT, x_fwd_bottom=x_eq - OR_LAND_ROOT - w_bot,
                  x_fwd_top=x_eq - OR_LAND_ROOT - groove_top)

    # angular layout: lugs centred between pockets; the channel swept by a lug from entry to
    # lock must stay FIN_POCKET_CLEAR (arc) clear of the pocket edges wherever both exist
    th_p = np.array(pk.theta)
    xc = np.linspace(x_ch0, x_lip0, 12)
    has_pk = pk.depth(xc) > 0
    margin_deg = math.degrees(FIN_POCKET_CLEAR / r_ch)
    engage_min = LUG_W_DEG / 2 + GAP_W_DEG / 2 + 2.0     # lug fully behind the lip at lock (+2)
    if has_pk.any():
        # land between two neighbouring pockets over the channel's x range (pocket extents are
        # asymmetric about the blade azimuth): from the max 'hi' edge of one pocket to the min
        # 'lo' edge of the next, FIN_POCKET_CLEAR (arc at the channel radius) inside each edge
        lo, hi = pk.theta_span(xc[has_pk])
        p_lo, p_hi = float(np.nanmin(lo)), float(np.nanmax(hi))
        land0 = p_hi + margin_deg
        land1 = 360.0 / LUG_N + p_lo - margin_deg
        engage_max = (land1 - land0) - LUG_W_DEG
        engage = float(min(ENGAGE_DEG, math.floor(engage_max)))
        slack = (land1 - land0) - LUG_W_DEG - engage
        lock_off = land1 - slack / 2 - LUG_W_DEG / 2       # lug centre at lock, rel. to the blade
        half_p = 0.5 * (p_hi - p_lo)
        edge_clear = slack / 2 + margin_deg
    else:
        p_lo = p_hi = 0.0
        engage_max = ENGAGE_DEG
        engage = ENGAGE_DEG
        lock_off = 60.0
        half_p = 0.0
        edge_clear = None
    lug_lock = (th_p + lock_off) % 360.0
    lug_entry = (lug_lock - engage) % 360.0
    seg_w = 360.0 / LUG_N - GAP_W_DEG
    root_ring_L = float(r_in(x_eq)) - r_sp_in             # 45 deg underside from the spigot bore
    return dict(
        x_eq=x_eq, L_seal=L_seal, L_spigot=x_eq - x_tip, root_ring_L=root_ring_L,
        r_bore=r_b, r_groove=r_gr, r_spigot_seal=r_s2, r_lug=r_lug, r_spigot_lug=r_s1,
        r_lip=r_lip, r_channel=r_ch, r_spigot_in=r_sp_in,
        x_seal0=x_seal0, x_lip0=x_lip0, x_tip=x_tip, x_ch0=x_ch0, x_ring0=x_ring0,
        groove=groove, oring=oring, bore_limit_d_mm=2 * float(rb_lim[xs_seal].min()),
        walls=dict(socket_at_pocket_mm=float((ro - pk.depth(xs) - r_b)[xs_seal].min()),
                   socket_plain_mm=float((ro - r_b)[xs_seal].min()),
                   spigot_under_groove_mm=r_gr - r_sp_in, spigot_lug_zone_mm=r_s1 - r_sp_in,
                   channel_backing_between_pockets_mm=float(r_out(x_tip) - r_ch)),
        bayonet=dict(n_lugs=LUG_N, lug_h_mm=LUG_H, lug_t_mm=LUG_T, lug_w_deg=LUG_W_DEG,
                     lug_w_mm=math.radians(LUG_W_DEG) * r_s1,
                     gap_w_deg=GAP_W_DEG, engage_deg=engage, engage_window_deg=[engage_min, engage_max],
                     engage_ok=engage_min <= engage <= engage_max + 1e-9,
                     lip_t_mm=LIP_T, lip_segment_deg=seg_w,
                     lock_theta_deg=sorted(lug_lock.tolist()),
                     entry_theta_deg=sorted(lug_entry.tolist()),
                     pocket_half_angle_at_channel_deg=half_p,
                     pocket_span_at_channel_rel_blade_deg=[p_lo, p_hi],
                     channel_to_pocket_edge_deg=edge_clear,
                     channel_to_pocket_edge_mm=(math.radians(edge_clear) * r_ch) if edge_clear is not None else None,
                     load_face="aft face of each lug and forward face of the lip: 45 deg cone on "
                               "both (the lug underside is an overhang when the aft shell prints "
                               "tip-up)",
                     lip_underside="45 deg chamfer from the seal bore in to the lip (overhang when "
                                   "the front shell prints equator-down)",
                     detent="2 magnet pairs (N52 4x2) in the equator faces, in two of the lands "
                            "between pockets; they attract at lock. Hard stop = the channel end "
                            "wall, reached when the pockets of both shells align"),
    )


def joint_strength(J: dict, F_d: float) -> dict:
    """Bayonet under the tow design load (rotor thrust + drags, x SF): the whole tow load
    passes from the nose eye (front shell) through the lugs into the aft shell.
    Both shells print with the body axis vertical, so the lug-root and lip-root shear
    planes carry shear across the layers: compared with the Prusament interlayer
    adhesion (tensile) as a proxy, with util < 0.5 required (shear strength not sourced)."""
    src = sources()
    adh = (src.get("prusament_petg_tds", {}).get("interlayer_adhesion_MPa") or [18])[0]
    w_lug = math.radians(LUG_W_DEG) * J["r_spigot_lug"]
    w_lip = math.radians(LUG_W_DEG) * J["r_channel"]
    tau_lug = F_d / (LUG_N * w_lug * LUG_T)
    tau_lip = F_d / (LUG_N * w_lip * LIP_T)
    # bearing on the 45 deg face: contact band = lug height x width / cos 45
    sig_b = F_d / (LUG_N * w_lug * LUG_H / math.cos(math.radians(45)))
    return dict(design_load_N=F_d, lug_w_mm=w_lug, tau_lug_root_MPa=tau_lug,
                tau_lip_root_MPa=tau_lip, contact_pressure_45deg_MPa=sig_b,
                interlayer_adhesion_PETG_MPa=adh, util_lug=tau_lug / adh, util_lip=tau_lip / adh,
                ok=max(tau_lug, tau_lip) / adh < 0.5)


# ==========================================================================
# Front-shell cavity and generator bay
# ==========================================================================


def cavity_clear_r(x, J: dict, pk: Pockets, rail_flush: bool = False) -> np.ndarray:
    """Clear radius (min over angle) available to the generator at x, with both shells
    closed. Front shell: hull inner wall and the socket ring (45 deg ceiling). From the
    spigot tip aft: the spigot bore of the aft shell (the generator, mounted in the front
    shell, slides into it when the shells close); further aft the pocket floors."""
    x = np.asarray(x, dtype=float)
    r = np.array(pk.floor_inner_r(x), dtype=float)
    ring = J["r_lip"] + (J["x_ring0"] - x)
    r = np.where(x >= J["x_ring0"], np.minimum(r, J["r_lip"]), np.minimum(r, ring))
    r = np.where(x >= J["x_tip"], np.minimum(np.array(pk.floor_inner_r(x)), J["r_spigot_in"]), r)
    if rail_flush:
        intr = r_out(x) - RAIL_H - T_MIN
        r = np.where((x >= RAIL_X[0] - 2) & (x <= RAIL_X[1] + 2), np.minimum(r, intr), r)
    return r


def ring_lip_r(x_fwd: float) -> float:
    """Inner radius of the shell's carrier ring where it overhangs the bay: the ring runs
    CARRIER_RING_OVER behind the carrier plate (x_fwd .. x_fwd + CARRIER_RING_OVER), CARRIER_RING
    deep from the shell inner wall."""
    xq = np.linspace(x_fwd, x_fwd + CARRIER_RING_OVER, 7)
    return float(np.min(np.asarray(r_in(xq)))) - CARRIER_RING


def carrier_r_needed(hw: dict) -> float:
    """The carrier (bulkhead) takes a motor-specific adapter plate; it is sized for the largest
    screw pattern among the recommended motor and its alternative."""
    m = hw["motor"]
    pcd = max(m["pcd"], m.get("alt_pcd") or 0.0)
    return pcd / 2 + m["cbore_d"] / 2 + CARRIER_EDGE


def bearing_layout(x_bell: float, ri: dict, hw: dict) -> dict:
    """Axial stations of coupler, bearings and seal for a generator bell face at x_bell, and
    the shaft-stability numbers (tow load case: C_T 0.8 at 2.5 m/s)."""
    x_end = ri["static_tail_x"]
    b6, b8, sl = hw["b608"], hw["b688"], hw["seal"]
    seal = (x_end - SEAL_RECESS - sl["B"], x_end - SEAL_RECESS)
    shoulder = (seal[0] - SEAL_SHOULDER_T, seal[0])
    rear = (shoulder[0] - b8["B"], shoulder[0])
    cA = (x_bell, x_bell + CLAW_A_L)
    cB = (x_bell + CLAW_A_L + CLAW_FLOAT, x_bell + COUPLER_STACK)
    front = (cB[1], cB[1] + b6["B"])
    xc_f, xc_r = sum(front) / 2, sum(rear) / 2
    span = xc_r - xc_f
    r_c = R_THRUST_FRAC * ri["r_tip"]                     # thrust centroid radius (ESTIMATE)
    x_load = ri.get("x_load") or (ri["hinge_x"] + (r_c - ri["hinge_r"])
                                  * math.tan(math.radians(ri["cone_deg"])))
    a = x_load - xc_r
    A_rot = math.pi * (ri["r_tip"] / 1000) ** 2
    T_tow = 0.5 * RHO * A_rot * V_TOW ** 2 * CT_TOW
    Fr = RADIAL_LOAD_FRAC * T_tow
    M_asym = THRUST_ASYM_FRAC * T_tow / P.N_BLADES * r_c / 1000      # N m
    R_rear = Fr * (span + a) / span + M_asym / (span / 1000)
    R_front = Fr * a / span + M_asym / (span / 1000)
    P0_front = max(R_front, 0.6 * R_front + 0.5 * T_tow)               # ISO 76, X0 0.6 Y0 0.5
    tilt = 2 * SEAT_SLOP / span
    # shaft (steel, rigid seats) as a beam on two supports with an overhang a: deflection at the
    # load plane from the radial load F (F a^2 (L + a) / 3EI) and the moment M (M a (2L + 3a) / 6EI)
    I = math.pi * P.SHAFT_D ** 4 / 64
    defl = (Fr * a ** 2 * (span + a) / (3 * E_STEEL * I)
            + M_asym * 1000 * a * (2 * span + 3 * a) / (6 * E_STEEL * I))
    s0f, s0r = b6["C0_N"] / P0_front, b8["C0_N"] / R_rear
    checks = dict(span_min=span >= SPAN_MIN, tip_runout=tilt * ri["r_tip"] <= TIP_RUNOUT_MAX,
                  static_safety=min(s0f, s0r) >= S0_MIN, shaft_deflection=defl <= SHAFT_DEFL_MAX)
    ok = all(checks.values())
    return dict(claw_A=cA, claw_B=cB, front=front, rear=rear, shoulder=shoulder, seal=seal,
                x_end=x_end, xc_f=xc_f, xc_r=xc_r, span=span, overhang=a, ratio=span / a,
                x_load=x_load, T_tow=T_tow, Fr=Fr, M_asym=M_asym, R_front=R_front, R_rear=R_rear,
                s0_front=s0f, s0_rear=s0r, tilt=tilt, shaft_defl=defl, checks=checks,
                tip_runout=tilt * ri["r_tip"], hub_face_runout=tilt * 20.0, ok=ok)


def x_bell_limit(ri: dict, hw: dict) -> float:
    """Aft-most generator face that still leaves a stable shaft (bearing_layout checks: SPAN_MIN,
    TIP_RUNOUT_MAX)."""
    xb = 20.0
    for x in np.arange(20.0, ri["static_tail_x"] - 40.0, 0.25):
        if bearing_layout(float(x), ri, hw)["ok"]:
            xb = float(x)
        else:
            break
    return xb


def alt_motor_L(hw: dict) -> float:
    """Body length of the hardware study's alternative motor (0 if unknown)."""
    mo = hw["motor"]
    return next((c["L_mm"] for c in hw["candidates"]
                 if mo.get("alt_name") and c["name"] == mo["alt_name"] and c["L_mm"]), 0.0)


def bay_length(hw: dict) -> float:
    """Chosen bay length: the longer of the recommended motor and its alternative's base body,
    + BAY_L_MARGIN, rounded up to 1 mm (merge decision; see generator_bay)."""
    return float(math.ceil(max(hw["motor"]["L"], alt_motor_L(hw)) + BAY_L_MARGIN))


def generator_bay(J: dict, pk: Pockets, ri: dict, hw: dict, rail_flush: bool = False,
                  x_bell_max: float | None = None) -> dict:
    """Largest cylinder (max diameter for each length) for the generator, with the motor
    carrier (bulkhead) ahead of it in the front shell and the bell face no further aft than
    the shaft-stability limit. Returns the frontier, the pick and the candidate fits."""
    xbm = x_bell_limit(ri, hw) if x_bell_max is None else x_bell_max
    xs = np.linspace(2.0, ri["static_tail_x"] - 30.0, 2400)
    rc = cavity_clear_r(xs, J, pk, rail_flush)
    r_car_min = carrier_r_needed(hw)
    frontier = []
    for L in np.arange(18.0, 40.01, 1.0):
        best = (0.0, None)
        for x_aft in np.arange(xbm, 20.0, -0.25):
            x_fwd = x_aft - L
            r_carrier = float(r_in(x_fwd - CARRIER_T)) - CARRIER_RING - P.CLEAR_SLIDE
            if r_carrier < r_car_min:
                continue
            m = (xs >= x_fwd) & (xs <= x_aft)
            R = min(float(rc[m].min()), J["r_lip"], ring_lip_r(x_fwd)) - C_GEN   # socket lip, carrier ring
            if 2 * R > best[0] + 1e-9:
                best = (2 * R, dict(x_fwd=x_fwd, x_aft=x_aft, r_carrier=r_carrier))
        frontier.append(dict(L_mm=float(L), D_max_mm=float(best[0]), **(best[1] or {})))
    env = hw["envelope"]
    cand = [f for f in frontier if 20 <= f["L_mm"] <= env["L_max"] and f["D_max_mm"] > 0]
    if not cand:                                         # nothing >= 20 long fits: report L 18-19
        cand = [f for f in frontier if f["D_max_mm"] > 0] or [dict(frontier[0], D_max_mm=0.0,
                                                                   x_fwd=float("nan"), x_aft=float("nan"))]
    best_vol = max(cand, key=lambda f: min(f["D_max_mm"], env["D_max"]) ** 2 * f["L_mm"])
    largest = dict(D_mm=float(math.floor(min(best_vol["D_max_mm"], env["D_max"]) * 2) / 2),
                   L_mm=best_vol["L_mm"], x_fwd_mm=best_vol["x_fwd"], x_aft_mm=best_vol["x_aft"],
                   rule="max D^2 L inside the candidate box (D <= %.0f, 20 <= L <= %.0f), D "
                        "rounded down to 0.5 mm" % (env["D_max"], env["L_max"]))

    def dmax_at(L):
        Ls = [f["L_mm"] for f in frontier]
        return float(np.interp(L, Ls, [f["D_max_mm"] for f in frontier]))

    def place(D, L, x_aft_max=None):
        """Forward-most placement of a D x L cylinder: the rotor face sets the 608 position
        (coupler stack), so forward-most gives the longest bearing span."""
        x_hi = xbm if x_aft_max is None else x_aft_max
        for x_aft in np.arange(20.0, x_hi + 1e-9, 0.25):
            x_fwd = x_aft - L
            r_carrier = float(r_in(x_fwd - CARRIER_T)) - CARRIER_RING - P.CLEAR_SLIDE
            if r_carrier < r_car_min:
                continue
            m = (xs >= x_fwd) & (xs <= x_aft)
            if min(float(rc[m].min()), J["r_lip"], ring_lip_r(x_fwd)) - C_GEN >= D / 2:
                return dict(x_fwd_mm=float(x_fwd), x_aft_mm=float(x_aft), r_carrier_mm=r_carrier)
        return None

    mo = hw["motor"]
    # chosen envelope: the largest one that fits (task rule), i.e. max D^2 L under the bearing-span limit
    # on the rotor face, placed forward-most; candidate motors are tested against it
    L_alt = alt_motor_L(hw)
    D_sel, L_sel = largest["D_mm"], largest["L_mm"]
    rule = (f"largest envelope that fits: D{D_sel:g} x L{L_sel:g} (max D^2 L with the rotor face at or ahead of "
            f"x {xbm:.2f} for the bearing span and the carrier ahead of the bay), placed forward-most")
    pl = place(D_sel, L_sel)
    if pl is None:
        pl = dict(x_fwd_mm=largest["x_fwd_mm"], x_aft_mm=largest["x_aft_mm"], r_carrier_mm=float("nan"))
    fits_tab = []
    for c in hw["candidates"]:
        if c["D_mm"] and c["L_mm"]:
            fits_tab.append(dict(name=c["name"], D_mm=c["D_mm"], L_mm=c["L_mm"], D_max_at_L_mm=dmax_at(c["L_mm"]),
                                 fits_bay=bool(c["D_mm"] <= D_sel and c["L_mm"] <= L_sel),
                                 radial_room_mm=(D_sel - c["D_mm"]) / 2, length_room_mm=L_sel - c["L_mm"]))
    pick = dict(D_mm=D_sel, L_mm=L_sel, **pl, rule=rule)
    pick["motor"] = mo["name"]
    pick["hardware_reserve_D_mm"] = mo["D_reserve"]
    pick["alt_motor_L_mm"] = L_alt
    pick["motor_body_x_mm"] = [pick["x_aft_mm"] - mo["L"], pick["x_aft_mm"]]
    pick["adapter_spacer_L_mm"] = L_sel - mo["L"]
    pick["clear_to_spigot_bore_mm"] = J["r_spigot_in"] - D_sel / 2
    # bay length trade at the hardware study's reserve diameter: forward-most placement with no span
    # limit, and what each length does to the bearing span
    trade = []
    for L_t in sorted({L_sel, float(mo["L_reserve"]), float(mo["L"]), 28.0}):
        p_ = place(mo["D_reserve"], L_t, x_aft_max=ri["static_tail_x"] - 45.0)
        if p_ is None:
            trade.append(dict(L_mm=L_t, fits=False))
            continue
        b_ = bearing_layout(p_["x_aft_mm"], ri, hw)
        trade.append(dict(L_mm=L_t, fits=True, x_fwd_mm=p_["x_fwd_mm"], x_rotor_face_mm=p_["x_aft_mm"],
                          span_mm=b_["span"], overhang_mm=b_["overhang"], span_over_overhang=b_["ratio"],
                          tip_runout_mm=b_["tip_runout"], shaft_defl_mm=b_["shaft_defl"],
                          s0_min=min(b_["s0_front"], b_["s0_rear"]), checks=b_["checks"], ok=b_["ok"]))
    return dict(frontier=frontier, largest=largest, pick=pick, candidate_fits=fits_tab,
                x_bell_max_mm=xbm, carrier_r_min_mm=r_car_min, rail_flush=rail_flush,
                length_trade=trade)


# ==========================================================================
# Drivetrain
# ==========================================================================


def drivetrain(J: dict, bay: dict, ri: dict, hw: dict, pk: Pockets) -> dict:
    b6, b8, sl = hw["b608"], hw["b688"], hw["seal"]
    x_bell = bay["pick"]["x_aft_mm"]
    BL = bearing_layout(x_bell, ri, hw)
    x_end = BL["x_end"]
    front, rear, seal, shoulder = BL["front"], BL["rear"], BL["seal"], BL["shoulder"]
    cA, cB = BL["claw_A"], BL["claw_B"]
    shaft = (cB[0], ri["shaft_aft_x"])
    tube_r_front = b6["D"] / 2 + TUBE_WALL
    tube_r_rear = b8["D"] / 2 + TUBE_WALL
    xx = np.linspace(front[0] - 1.0, x_end, 300)
    tube_r = np.where(xx <= front[1] + 2, tube_r_front,
                      np.interp(xx, [front[1] + 2, rear[0] - 2, rear[0]],
                                [tube_r_front, tube_r_rear + 1.0, tube_r_rear]))
    gap = pk.floor_inner_r(xx) - tube_r
    k_ro = pk.keepout_r(xx)
    gap_ko = np.where(np.isfinite(k_ro), k_ro - tube_r, np.inf)
    trade = []
    for xb in np.arange(36.0, BL["x_end"] - 45.0, 2.0):
        b = bearing_layout(float(xb), ri, hw)
        trade.append(dict(x_bell_mm=float(xb), span_mm=b["span"], overhang_mm=b["overhang"],
                          span_over_overhang=b["ratio"], tip_runout_mm=b["tip_runout"],
                          shaft_defl_mm=b["shaft_defl"], R_rear_N=b["R_rear"], ok=b["ok"]))
    # rotating claw B vs the static tube mouth (x_from .. 608 front face): radial clearance of each
    # claw B segment that overlaps the mouth axially
    x_mouth = (front[0] - 1.0, front[0])
    segs = [("claw B body OD %g" % COUPLER_OD, cB[0], cB[1] - 0.5 - CLAW_B_STEP_L, COUPLER_OD),
            ("claw B step OD %g" % CLAW_B_STEP_OD, cB[1] - 0.5 - CLAW_B_STEP_L, cB[1] - 0.5, CLAW_B_STEP_OD),
            ("shim OD %g" % SHIM_OD, cB[1] - 0.5, cB[1], SHIM_OD)]
    mouth = []
    for nm, x0_, x1_, od in segs:
        if x1_ > x_mouth[0] + 1e-9 and x0_ < x_mouth[1] - 1e-9:
            mouth.append(dict(part=nm, x_mm=[max(x0_, x_mouth[0]), min(x1_, x_mouth[1])],
                              radial_clear_mm=(b6["D"] + 2 * P.PRESS_FIT) / 2 - od / 2))
    mouth_clear = min((m["radial_clear_mm"] for m in mouth), default=math.inf)
    return dict(
        stations_mm=dict(claw_A=list(cA), claw_B_thrust_collar=list(cB), bearing_front_608=list(front),
                         bearing_rear_688=list(rear), seal_shoulder=list(shoulder), lip_seal=list(seal),
                         static_tail_end=x_end, shaft=list(shaft), rotor_load_plane=BL["x_load"]),
        bearings=dict(front=b6["name"], rear=b8["name"],
                      why="608 in the wide front where it takes the thrust; 688 in the narrow tail "
                          "under the deepest pockets; stepped bore (22 -> 18 -> 16) so both press in "
                          "from the front and the 688 passes through the 608 seat; seal bore 16 = "
                          "688 bore (one reamed diameter)",
                      ratings=dict(b608=dict(C_N=b6["C_N"], C0_N=b6["C0_N"], src=b6["src"]),
                                   b688=dict(C_N=b8["C_N"], C0_N=b8["C0_N"], src=b8["src"])),
                      spacing_mm=BL["span"], overhang_mm=BL["overhang"], spacing_over_overhang=BL["ratio"],
                      criteria=dict(span_min_mm=SPAN_MIN, tip_runout_max_mm=TIP_RUNOUT_MAX,
                                    seat_slop_mm=SEAT_SLOP, s0_min=S0_MIN,
                                    shaft_deflection_max_mm=SHAFT_DEFL_MAX,
                                    span_over_overhang_reference=SPAN_RATIO_REF,
                                    note="span / overhang is reported against the 1.25 rule of thumb "
                                         "only; it cannot be met with the static tail at x %.0f and the "
                                         "coned rotor (load plane x %.1f)" % (x_end, BL["x_load"])),
                      checks=BL["checks"], ok=BL["ok"], shaft_deflection_at_load_plane_mm=BL["shaft_defl"],
                      radial_load_N=BL["Fr"], asym_moment_Nm=BL["M_asym"], thrust_tow_N=BL["T_tow"],
                      reaction_front_N=BL["R_front"], reaction_rear_N=BL["R_rear"],
                      static_safety_front_s0=BL["s0_front"], static_safety_rear_s0=BL["s0_rear"],
                      tilt_mrad=1000 * BL["tilt"], tip_axial_runout_mm=BL["tip_runout"],
                      hub_face_runout_at_r20_mm=BL["hub_face_runout"], hub_gap_mm=P.HUB_GAP,
                      load_basis="tow case: T = 0.5 rho A v^2 C_T (C_T 0.8, 2.5 m/s); radial load 10 % "
                                 "of T and one blade 20 % off its share of T (both ESTIMATES)",
                      location="608 locating: claw B (thrust collar, 2 M3 grub screws on a shaft flat "
                               "+ a retaining-ring groove, VERIFY with the shaft supplier) bears on the "
                               "608 inner race for +X thrust; outer race on the 18 mm tube shoulder. 688 "
                               "floating (slip fit on the shaft). -X handling loads close the hub gap"),
        span_trade=trade,
        tube=dict(r_out_front_mm=tube_r_front, r_out_rear_mm=tube_r_rear,
                  bore_front_mm=b6["D"], shoulder_front_id_mm=FRONT_SHOULDER_ID,
                  bore_rear_mm=b8["D"], shoulder_rear_id_mm=SEAL_SHOULDER_ID, seal_bore_mm=sl["D"],
                  x_from=front[0] - 1.0, x_to=x_end,
                  min_gap_to_pocket_floor_mm=float(gap.min()),
                  min_gap_to_pocket_floor_x_mm=float(xx[int(np.argmin(gap))]),
                  min_gap_to_rotor_keepout_mm=float(gap_ko.min()),
                  pocket_floor_wall_mm=POCKET_FLOOR_WALL,
                  ribs="3 radial ribs 2.0 mm thick under the pocket floors (theta = pocket angles), "
                       "from the tube front to the tail end wall; the sectors between them stay free "
                       "(ballast tray in the bottom one)"),
        seal=dict(name=sl["name"], x=list(seal), lip="main lip faces aft (water side)",
                  bore_mm=sl["D"], note="v1 display: pocket only, seal optional"),
        shaft=dict(d_mm=P.SHAFT_D, x_fwd=shaft[0], x_aft=shaft[1], length_mm=shaft[1] - shaft[0],
                   aft_end_rule=ri.get("shaft_aft_rule", "data/layout_rotor.json"),
                   stock="8 mm g6/h6 stainless ground rod (hardware study: hardened 440C for the "
                         "lip seal); two flats (claw B; hub D-flat 0.5 deep under the M3 grub at x "
                         "%.1f), one retaining-ring groove ahead of claw B and an M3 x 8 tapped hole "
                         "in the aft end (spinner screw, carries the rotor thrust), ordered machined; "
                         "chamfer both ends" % ((ri.get("rotor") or {}).get("grub_x", float("nan"))),
                   hub_engagement_mm=shaft[1] - ((ri.get("rotor") or {}).get("hub_front_x",
                                                                              x_end + P.HUB_GAP)),
                   mass_g=math.pi * (P.SHAFT_D / 20) ** 2 * (shaft[1] - shaft[0]) / 10 * RHO_STEEL),
        coupler=dict(type="two-dog claw coupling, V lead-in, 1.0 mm axial float, 0.3 mm radial "
                          "clearance per dog: engages blind as the shells close and tolerates "
                          "shell-to-shell misalignment, so the motor's own bearings and the 608/688 "
                          "pair are never forced into line (no 4-bearing shaft)",
                     od_mm=COUPLER_OD, claw_A=list(cA), claw_B=list(cB),
                     claw_B_segments=[dict(part=nm, x_mm=[x0_, x1_], od_mm=od) for nm, x0_, x1_, od in segs],
                     tube_mouth=dict(x_mm=list(x_mouth), bore_d_mm=b6["D"] + 2 * P.PRESS_FIT,
                                     lead_in="%.1f x 45 deg chamfer at x %.1f" % (TUBE_LEAD_IN, x_mouth[0]),
                                     overlaps=mouth, min_radial_clear_mm=mouth_clear,
                                     note="the rotating claw B steps down to OD %g over its last %g mm, so it "
                                          "runs %.2f mm clear of the static 608 seat bore where it enters the "
                                          "tube mouth" % (CLAW_B_STEP_OD, CLAW_B_STEP_L, mouth_clear)),
                     claw_A_mount="v2: on the bell face holes if the motor has them (unconfirmed for the "
                                  "GM4108H), else on an 8 mm steel stub bonded into the motor's hollow "
                                  "shaft (ID 8) with retaining compound; display: printed replica",
                     fits_spigot_bore=COUPLER_OD / 2 < J["r_spigot_in"]),
    )


# ==========================================================================
# Nose: tow eye and cable gland
# ==========================================================================


def nose_plug() -> dict:
    """Front shell prints nose-up: forward of x_45 the inner wall leans more than 45 deg from
    vertical (|dr/dx| > 1) and no 45 deg cone fits inside the shell there (the elliptic nose
    closes faster than any printable ceiling). The tip forward of x_45 is therefore a plug
    (perimeters + infill) whose aft face is a flat bridge across the cavity; the plug also
    roots the tow eye."""
    xs = np.linspace(0.5, P.HULL_NOSE_A, 4000)
    ri_ = np.array(r_in(xs))
    d = np.gradient(ri_, xs)
    i = int(np.argmax((d <= 1.0) & (ri_ > 1.0)))          # first x where the inner wall is <= 45 deg
    x45, r45 = float(xs[i]), float(ri_[i])
    V_cav, _ = _rev_vol(0.0, x45, r_in)                   # cavity volume the plug fills
    return dict(x_45_mm=x45, r_inner_at_45_mm=r45, bridge_d_mm=2 * r45, plug_cavity_cc=V_cav / 1000,
                rule="inner wall |dr/dx| = 1 is the 45 deg overhang limit when printed nose-up; "
                     "forward of it the cavity is filled (infill) and closed by a %.1f mm bridge"
                     % (2 * r45),
                alternative="leave the dome open: no supports either, but the hidden inner ceiling "
                            "forward of x_45 sags (overhang up to 90 deg at the apex)")


def _cyl_samples(S, A, u, v, R, t0, t1, nt=25, nphi=73, nr=9):
    """Points filling a cylinder of radius R on the axis S + t A, t in [t0, t1] (u, v span the
    cross-section)."""
    t = np.linspace(t0, t1, nt)[:, None, None]
    ph = np.linspace(0.0, 2 * math.pi, nphi)[None, :, None]
    rr = np.linspace(0.0, R, nr)[None, None, :]
    P_ = S + t[..., None] * A + rr[..., None] * (np.cos(ph)[..., None] * u + np.sin(ph)[..., None] * v)
    return P_.reshape(-1, 3)


def gland_geom(x_s: float, tilt_deg: float, hw: dict, theta_deg: float = GLAND_THETA) -> dict:
    """Gland on the nose: seat centre on the hull at x_s (azimuth theta), outward axis turned by
    tilt from the surface normal toward -X. Returns the axis, the boss samples (seat D18 over the
    shell wall, thread boss D16 to the thread length), the seat rim height above / below the hull
    and the external gland body samples (hex across corners as a cylinder, height - thread)."""
    r_s = float(r_out(x_s))
    sl = float(slope(x_s))
    a = math.atan2(1.0, -sl) + math.radians(tilt_deg)        # outward axis angle from +X, meridian plane
    t = math.radians(theta_deg)
    er = np.array([0.0, math.sin(t), math.cos(t)])
    A = math.cos(a) * np.array([1.0, 0.0, 0.0]) + math.sin(a) * er
    S = np.array([x_s, 0.0, 0.0]) + r_s * er
    u = np.array([0.0, math.cos(t), -math.sin(t)])
    v = np.cross(A, u)
    v /= np.linalg.norm(v)
    L = hw["gland"]["thread_L"]
    boss = np.vstack([_cyl_samples(S, A, u, v, GLAND_SEAT_D / 2, -GLAND_SEAT_T, 0.0, nt=7),
                      _cyl_samples(S, A, u, v, GLAND_BOSS_D / 2, -L, -GLAND_SEAT_T)])
    ph = np.linspace(0.0, 2 * math.pi, 181)
    rim = S + (GLAND_SEAT_D / 2) * (np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * v)
    rr = np.hypot(rim[:, 1], rim[:, 2])
    hr = np.where(rim[:, 0] < 0, 0.0, r_out(np.clip(rim[:, 0], 0.0, P.HULL_L)))
    r_body = hw["gland"]["SW"] / math.sqrt(3.0)              # hex across corners / 2
    body = _cyl_samples(S, A, u, v, r_body, 0.0, hw["gland"]["height"] - L, nr=2)
    return dict(x_s=x_s, r_s=r_s, tilt_deg=tilt_deg, axis_angle_from_x_deg=math.degrees(a), S=S, A=A,
                inner_end=S - L * A, boss=boss, body=body, r_body=r_body,
                seat_raise_mm=float((rr - hr).max()), seat_recess_mm=float((hr - rr).max()))


def _gland_clearances(g: dict, car_front: float, r_carrier: float, eye_c: float) -> dict:
    """Clearances of a gland placement: boss to the carrier plate front face and to the carrier
    ring (inward ring on the shell, from car_front - CARRIER_RING_OVER), external body to the eye
    ring (XZ plane) and to the shackle envelope."""
    B = g["boss"]
    xb, rb = B[:, 0], np.hypot(B[:, 1], B[:, 2])
    near_plate = rb <= r_carrier + GLAND_CLEAR
    c_plate = float((car_front - xb[near_plate]).min()) if near_plate.any() else math.inf
    x_ring0 = car_front - CARRIER_RING_OVER
    ring_r = np.asarray(r_in(np.clip(xb, 0.0, P.HULL_L))) - CARRIER_RING
    near_ring = rb >= ring_r - GLAND_CLEAR - 3.0
    c_ring = float((x_ring0 - xb[near_ring]).min()) if near_ring.any() else math.inf
    Bd = g["body"]
    # eye ring in the XZ plane: annulus about (eye_c, 0, 0) in (x, z), |y| <= EYE_T / 2
    rp = np.hypot(Bd[:, 0] - eye_c, Bd[:, 2])
    dr = np.maximum(np.maximum(EYE_HOLE_D / 2 - rp, rp - (EYE_HOLE_D / 2 + EYE_W)), 0.0)
    dy = np.maximum(np.abs(Bd[:, 1]) - EYE_T / 2, 0.0)
    c_eye = float(np.sqrt(dr ** 2 + dy ** 2).min())
    x_hole_aft = eye_c + EYE_HOLE_D / 2
    dx = np.maximum(Bd[:, 0] - x_hole_aft, 0.0)
    dpin = np.maximum(np.abs(Bd[:, 1]) - SHACKLE_HALF_PIN, 0.0)
    dacr = np.maximum(np.abs(Bd[:, 2]) - SHACKLE_HALF_ACROSS, 0.0)
    c_sh = float(np.sqrt(dx ** 2 + dpin ** 2 + dacr ** 2).min())
    return dict(boss_to_carrier_plate_mm=c_plate, boss_to_carrier_ring_mm=c_ring, body_to_eye_ring_mm=c_eye,
                body_to_shackle_envelope_mm=c_sh, boss_x_max_mm=float(xb.max()),
                ok=(c_plate >= GLAND_CLEAR and c_ring >= GLAND_CLEAR and c_eye >= GLAND_TOW_CLEAR
                    and c_sh >= GLAND_TOW_CLEAR))


def gland_place(hw: dict, bay: dict, eye_c: float) -> dict:
    """Gland position and axis on the nose (theta GLAND_THETA): among seat positions x and axis
    tilts (GLAND_X_SEARCH, GLAND_TILT_SEARCH) whose whole boss volume stays GLAND_CLEAR ahead of
    the motor carrier and its ring, and whose external body stays GLAND_TOW_CLEAR from the eye and
    the shackle envelope, take the one with the lowest seat rim above the hull (smallest raised
    spot face), ties to the larger carrier clearance."""
    car_front = bay["pick"]["x_fwd_mm"] - CARRIER_T
    r_car = bay["pick"]["r_carrier_mm"]
    best, n_ok, n_all = None, 0, 0
    for tilt in np.arange(GLAND_TILT_SEARCH[0], GLAND_TILT_SEARCH[1] + 1e-9, GLAND_TILT_SEARCH[2]):
        for xs in np.arange(GLAND_X_SEARCH[0], GLAND_X_SEARCH[1] + 1e-9, GLAND_X_SEARCH[2]):
            g = gland_geom(float(xs), float(tilt), hw)
            c = _gland_clearances(g, car_front, r_car, eye_c)
            n_all += 1
            if not c["ok"]:
                continue
            n_ok += 1
            key = (round(g["seat_raise_mm"], 3), -min(c["boss_to_carrier_plate_mm"], c["boss_to_carrier_ring_mm"]))
            if best is None or key < best[0]:
                best = (key, g, c)
    if best is None:
        raise ValueError("no gland placement clears the carrier and the tow eye")
    _, g, c = best
    # the original placement (x 12 on the surface normal) for the record
    g12 = gland_geom(12.0, 0.0, hw)
    c12 = _gland_clearances(g12, car_front, r_car, eye_c)
    return dict(g=g, c=c, n_feasible=n_ok, n_searched=n_all, car_front=car_front, r_carrier=r_car,
                previous_x12=dict(boss_x_max_mm=c12["boss_x_max_mm"],
                                  boss_to_carrier_plate_mm=c12["boss_to_carrier_plate_mm"]))


def cable_notch(hw: dict, bay: dict) -> dict:
    """Notch in the carrier rim at theta 180 for the gland cable: cable OD + CABLE_CLEAR wide and
    deep, so the cable passes inside the shell's carrier ring (no ring notch). The cable is soldered
    to the motor leads on the bench and turns with the carrier through its 30 deg twist (slack loop
    in the nose cavity), so the notch needs no twist arc. Checked against the motor screw
    counterbores (pattern clocked 45 deg) and the centre bore."""
    d_c = float(max(hw["gland"]["clamp_mm"]))
    w = d_c + CABLE_CLEAR
    r_car = bay["pick"]["r_carrier_mm"]
    r_bot = r_car - w
    x_c = bay["pick"]["x_fwd_mm"] - CARRIER_T / 2
    ring_in = float(r_in(x_c)) - CARRIER_RING
    mo = hw["motor"]
    clear_screws = []
    for pcd in {mo["pcd"], mo.get("alt_pcd") or mo["pcd"]}:
        r_s = pcd / 2
        for th_s in (135.0, 225.0):                         # pattern clocked 45 deg: nearest screws to 180
            # nearest distance from the counterbore circle to the notch (a w-wide radial slot at 180)
            p = np.array([r_s * math.sin(math.radians(th_s)), r_s * math.cos(math.radians(th_s))])
            # notch: |y| <= w/2, z from -r_car to -r_bot (theta 180 is -Z)
            dy = max(abs(p[0]) - w / 2, 0.0)
            dz = max(-r_car - p[1], p[1] - (-r_bot), 0.0)
            clear_screws.append(math.hypot(dy, dz) - mo["cbore_d"] / 2)
    return dict(theta_deg=180.0, width_mm=w, depth_from_rim_mm=w, bottom_r_mm=r_bot, cable_d_max_mm=d_c,
                cable_inside_ring=bool(r_car <= ring_in), carrier_ring_inner_r_mm=ring_in,
                web_to_centre_bore_mm=r_bot - 8.0, clear_to_screw_counterbores_mm=float(min(clear_screws)),
                ok=bool(min(clear_screws) >= T_MIN and r_bot - 8.0 >= T_MIN),
                note="the cable lies in the notch, inside the shell's carrier ring, and turns with the carrier "
                     "through the 30 deg lock twist; leave a slack loop (~%d mm, ESTIMATE) in the nose cavity"
                     % 120)


def nose(ri: dict, fins: dict, bay: dict, hw: dict) -> dict:
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
        sig_leg = F_d / (2 * EYE_W * EYE_T)               # legs: tension across the layers
        Rm = EYE_HOLE_D / 2 + EYE_W / 2
        sig_crown = 6 * (F_d * Rm / math.pi) / (EYE_T * EYE_W ** 2)   # pin-loaded ring, in-layer
        mats[mat] = dict(interlayer_adhesion_MPa=adh, xy_tensile_yield_MPa=sxy,
                         eye_leg_stress_MPa=sig_leg, eye_leg_util=sig_leg / adh,
                         eye_crown_stress_MPa=sig_crown, eye_crown_util=sig_crown / sxy,
                         bearing_on_8mm_pin_MPa=F_d / (EYE_HOLE_D * EYE_T))
    cone = nose_plug()
    eye_c = -(EYE_HOLE_D / 2 + EYE_STANDOFF)
    gp = gland_place(hw, bay, eye_c)
    g, c = gp["g"], gp["c"]
    notch = cable_notch(hw, bay)
    return dict(
        loads=dict(thrust_2p5_m_s_N=T, thrust_bem_2p5_m_s_N=ri.get("thrust_bem_2p5_N"),
                   hull_drag_est_N=D_hull, hull_cd_est=CD_HULL, fin_drag_N=D_fins,
                   working_load_N=F_work, safety_factor=SF_TOW, design_load_N=F_d,
                   basis="T = 0.5 rho pi R^2 v^2 C_T with C_T 0.8 at 2.5 m/s (task rule); BEM "
                         "cross-check from data/blade_design.json; hull drag C_D 0.10 on the "
                         "frontal area (ESTIMATE)"),
        eye=dict(type="printed ring at the nose tip, on the axis, ring plane = %s (hole axis Y)" % EYE_PLANE,
                 hole_d_mm=EYE_HOLE_D, ring_w_mm=EYE_W, ring_t_mm=EYE_T, hole_centre_x_mm=eye_c,
                 plane_why="with the hole axis along Y the shackle pin runs along Y and the shackle body lies "
                           "in the XY plane, so the gland under the nose (theta 180) stays clear of it",
                 shackle_envelope=dict(model="ESTIMATE, 5 mm stainless D-shackle: from the hole's aft edge "
                                             "forward, |y| <= %g (pin), |z| <= %g" % (SHACKLE_HALF_PIN,
                                                                                    SHACKLE_HALF_ACROSS)),
                 print="front shell prints nose-up: the ring stands vertical on the nose plug, "
                       "hole printed as a 45 deg teardrop",
                 materials=mats,
                 note="legs load the layers in tension (weakest direction) and still sit far below "
                      "the Prusament interlayer adhesion; use a stainless shackle or a thimble, not "
                      "a bare rope, to keep the bearing stress low. v2 pull test >= design load"),
        nose_plug=cone,
        gland=dict(part=hw["gland"]["name"], thread="M12x1.5", thread_engagement_mm=hw["gland"]["thread_L"],
                   hex_SW_mm=hw["gland"]["SW"], height_mm=hw["gland"]["height"], cable_mm=hw["gland"]["clamp_mm"],
                   theta_deg=GLAND_THETA, seat_x_mm=g["x_s"], seat_surface_r_mm=g["r_s"],
                   seat_d_mm=GLAND_SEAT_D, boss_d_mm=GLAND_BOSS_D,
                   boss_model="seat D%g over the %.1f mm shell wall + thread boss D%g (M12 + 2 x %.1f) to the "
                              "%.0f mm thread length, all along the gland axis" % (GLAND_SEAT_D, GLAND_SEAT_T,
                                                                                 GLAND_BOSS_D, T_MIN,
                                                                                 hw["gland"]["thread_L"]),
                   axis_tilt_from_normal_deg=g["tilt_deg"], axis_angle_from_x_deg=g["axis_angle_from_x_deg"],
                   axis="outer surface normal turned %g deg toward -X (points forward-down)" % g["tilt_deg"],
                   inner_end_x_mm=float(g["inner_end"][0]),
                   inner_end_r_mm=float(np.hypot(g["inner_end"][1], g["inner_end"][2])),
                   seat_raise_mm=g["seat_raise_mm"],
                   seat_note="the D%g spot face is flat on a curved nose: its rim stands up to %.1f mm proud of "
                             "the hull (a raised boss, prints nose-up without support: its flank leans %.0f deg "
                             "from vertical)" % (GLAND_SEAT_D, g["seat_raise_mm"], 180.0 - g["axis_angle_from_x_deg"]),
                   carrier_front_x_mm=gp["car_front"], boss_x_max_mm=c["boss_x_max_mm"],
                   boss_to_carrier_plate_mm=c["boss_to_carrier_plate_mm"],
                   boss_to_carrier_ring_mm=c["boss_to_carrier_ring_mm"],
                   body_to_eye_ring_mm=c["body_to_eye_ring_mm"],
                   body_to_shackle_envelope_mm=c["body_to_shackle_envelope_mm"], clear_ok=bool(c["ok"]),
                   placement_rule="search x %g-%g mm (step %g) and tilt %g to %g deg (step %g): whole boss "
                                  ">= %.1f mm ahead of the carrier plate and its ring, gland body >= %.1f mm from "
                                  "the eye ring and the shackle envelope; pick the lowest seat rim (%d of %d "
                                  "placements feasible)" % (*GLAND_X_SEARCH, *GLAND_TILT_SEARCH, GLAND_CLEAR,
                                                            GLAND_TOW_CLEAR, gp["n_feasible"], gp["n_searched"]),
                   previous_x12=dict(gp["previous_x12"], note="the earlier boss at x 12 on the surface normal "
                                                              "ran through the carrier plate into the bay"),
                   seal="flat NBR washer under the gland head (IP68 with the extra O-ring, hardware.json)",
                   clearance_hole_if_unthreaded_mm=GLAND_HOLE,
                   thread_form="tapped M12x1.5 in the solid boss (or printed thread), no inner locknut; the "
                               "bore runs on through the nose plug into the cavity",
                   cable_route="gland -> nose cavity (slack loop) -> notch in the carrier rim at theta 180 "
                               "-> motor leads behind the carrier", cable_notch=notch),
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


_GEN_LL: dict = {}


def gen_load_line() -> dict:
    """Operating points of the recommended generator on its load line (final blade, BEM study:
    data/blade_design.json generator_coupling_estimate.load_line_by_flow, the 'capped' point =
    best net power with the cavitation cap), by flow speed; {} if the BEM file lacks it."""
    if "points" in _GEN_LL or "none" in _GEN_LL:
        return _GEN_LL
    bd = _load(BLADE_JSON) or {}
    gc = bd.get("generator_coupling_estimate") or {}
    rec = next((g for g in (gc.get("load_line_by_flow") or {}).get("generators", [])
                if g.get("role") == "recommended"), None)
    if not rec:
        _GEN_LL["none"] = True
        return _GEN_LL
    best = next((r_ for r_ in gc.get("rows", []) if r_.get("role") == "recommended" and r_.get("blade") == "final"),
                None)
    pts = {}
    for p_ in rec["points"]:
        q = dict(p_.get("capped") or p_["unconstrained"])
        q["cap_binding"] = bool(p_.get("cap_binding"))
        pts[round(float(p_["flow_m_s"]), 2)] = q
    _GEN_LL.update(generator=rec["generator"], short=rec.get("short"), kv=rec.get("kv"), r_ll_ohm=rec.get("r_ll_ohm"),
                   points=pts, best_net_1p5=(best or {}).get("best_net"),
                   source="data/blade_design.json generator_coupling_estimate.load_line_by_flow (recommended "
                          "generator, final blade, 'capped' point = best net power with the cavitation cap)",
                   model=gc.get("model"))
    return _GEN_LL


def roll_torque(v: float, ri: dict) -> dict:
    """Shaft torque the hull must react at flow v: the recommended generator's load-line point
    (rotor torque = generator electromagnetic + drag torque, both react on the stator, i.e. the
    hull). The old layout rule (Cp 0.35 at TSR 4) is kept as Q_cp_rule_Nm for reference."""
    R = ri["r_tip"] / 1000
    P_ = CP_DESIGN * 0.5 * RHO * math.pi * R ** 2 * v ** 3
    om = P.DESIGN_TSR * v / R
    ll = gen_load_line()
    pt = (ll.get("points") or {}).get(round(float(v), 2))
    if pt:
        return dict(v_m_s=v, Q_Nm=pt["torque_Nm"], tsr=pt["tsr"], rpm=pt["rpm"], shaft_power_W=pt["shaft_power_W"],
                    net_power_W=pt["net_power_W"], cap_binding=pt["cap_binding"],
                    basis="%s load line (BEM, capped point)" % ll["generator"], Q_cp_rule_Nm=P_ / om)
    return dict(v_m_s=v, P_W=P_, omega_rad_s=om, rpm=om * 60 / (2 * math.pi), Q_Nm=P_ / om, tsr=P.DESIGN_TSR,
                basis="layout rule Cp %.2f at TSR %g (no generator load line in data/blade_design.json)"
                      % (CP_DESIGN, P.DESIGN_TSR), Q_cp_rule_Nm=P_ / om)


def upstream_factor(x_mm: float, x_rotor: float, R_mm: float) -> float:
    """1 - a_up: on-axis actuator-disc deceleration at x (conservative off axis)."""
    a = (1 - math.sqrt(1 - CT_TOW)) / 2
    xi = (x_mm - x_rotor) / R_mm
    return 1 - a * (1 + xi / math.sqrt(1 + xi * xi))


def fin_root_line(x_le: float, x_te: float) -> dict:
    """Straight dovetail rail on a raised pad at theta 90 / 270 under the fin root chord: the
    line through the hull points at the chord ends, lifted so the pad is >= RAIL_H high
    everywhere (the hull is convex in between). The fin root edge sits on this line."""
    xs = np.linspace(x_le, x_te, 80)
    ro = np.array(r_out(xs))
    k = (ro[-1] - ro[0]) / (x_te - x_le)
    line0 = ro[0] + k * (xs - x_le)
    lift = float((ro - line0).max()) + RAIL_H
    pad = line0 + lift - ro
    return dict(x_le=x_le, x_te=x_te, y_le=float(ro[0] + lift), y_te=float(ro[-1] + lift),
                y_mean=float(np.mean(line0 + lift)), incl_deg=math.degrees(math.atan(k)),
                pad_min_mm=float(pad.min()), pad_max_mm=float(pad.max()))


def blade_sweep_len(ri: dict) -> float:
    """Largest distance of the blade from its pin axis (radius of the torus it sweeps while
    folding): from the rotor half's solved axis and sections when available, else a straight
    coned blade."""
    if ri.get("sweep_len"):
        return ri["sweep_len"]
    return (ri["r_tip"] - ri["hinge_r"]) / math.cos(math.radians(ri["cone_deg"]))


def fin_te_limit(ri: dict, c: float, rule: str = FIN_FOLD_RULE,
                 x_eq: float = P.EQUATOR_X) -> tuple[float, dict]:
    """Aft-most fin TE. 'torus': the TE root point stays FIN_TIP_CLEAR outside the torus swept
    by the blades about the hinge circle (any fold angle, any rotor azimuth). 'index': the fins
    live on the front shell (TE at the rail end stop ahead of the equator); clearance to the
    blades folding at the stowing index is checked separately (fin_fold_clearance)."""
    if rule == "index":
        x_te = x_eq - FIN_TE_EQ_MARGIN
        return x_te, fin_root_line(x_te - c, x_te)
    need = blade_sweep_len(ri) + FIN_TIP_CLEAR
    x_te = ri["hinge_x"] - need
    for _ in range(6):
        rl = fin_root_line(x_te - c, x_te)
        dz = rl["y_te"] - ri["hinge_r"]
        x_te = float(ri["hinge_x"] - math.sqrt(max(need ** 2 - dz ** 2, 0.0)))
    return x_te, fin_root_line(x_te - c, x_te)


def fin_te_sense(ri: dict) -> float:
    """Tangential side (+1 = +theta) toward which the fin TE is yawed. The fins must give a -X
    couple against the +X generator reaction (rotor omega +X): the +Y fin lifts toward -Z, so its
    TE points toward +Z, which is -theta at theta 90; the second fin is the first turned 180 deg
    about X, so both TEs point toward -theta."""
    return -1.0 if ri.get("omega_vector", "+X") == "+X" else 1.0


def fin_frame(x_te: float, c: float, inc_deg: float, pivot: float, sense: float,
              x_eq: float = P.EQUATOR_X) -> dict:
    """Yawed fin root and rail pad. The root chord (length c) and its pad are yawed by the
    incidence about the radial line through the pivot (chord fraction from the LE), so the rail
    runs along the fin chord and the fin prints flat with its root block. The TE sits on the rail
    end stop x_te. Local coordinates at a fin angle: u chordwise from the LE, w across the chord
    (+ toward +theta at zero yaw), 'along' radial. The pad (RAIL_W wide, top = the root line)
    runs from the LE aft and ends with a 45 deg chamfer from the TE if that fits ahead of the
    equator, else it runs down to the equator face (the front shell prints on that face, so the
    pad's aft end stands on the bed: no overhang)."""
    i = math.radians(inc_deg)
    ci, si = math.cos(i), sense * math.sin(i)
    x_le = x_te - c * ci
    x_p = x_le + pivot * c * ci
    rl = fin_root_line(x_le, x_te)
    k = (rl["y_te"] - rl["y_le"]) / (x_te - x_le)
    h_te = rl["y_te"] - float(r_out(x_te))
    chamfer = x_te + h_te <= x_eq
    return dict(c=c, inc_deg=inc_deg, pivot=pivot, sense=sense, ci=ci, si=si, x_le=x_le, x_te=x_te,
                x_p=x_p, rl=rl, k=k, h_te=h_te, chamfer=chamfer, x_end=x_te + h_te if chamfer else x_eq)


def fin_offset(fg: dict, x) -> np.ndarray:
    """Tangential offset (mm, + = +theta) of the yawed chord line from the fin angle line at x."""
    return (np.asarray(x, float) - fg["x_p"]) * fg["si"] / fg["ci"]


def _root_y(fg: dict, x):
    return fg["rl"]["y_le"] + fg["k"] * (np.asarray(x, float) - fg["x_le"])


def _pad_top(fg: dict, x):
    x = np.asarray(x, float)
    y = _root_y(fg, x)
    if fg["chamfer"]:
        y = np.where(x > fg["x_te"], np.minimum(y, r_out(x) + fg["h_te"] - (x - fg["x_te"])), y)
    return y


def _to_fin(fg: dict, x, lat):
    dx = np.asarray(x, float) - fg["x_p"]
    u = dx * fg["ci"] + lat * fg["si"] + fg["pivot"] * fg["c"]
    w = -dx * fg["si"] + lat * fg["ci"]
    return u, w


def pad_footprint(fg: dict, du: float = 0.25, nw: int = 11):
    """Sampled pad footprint on the hull as (x, lat): u from the LE aft, |w| <= RAIL_W / 2, cut at
    the pad end x."""
    u_max = (fg["x_end"] - fg["x_p"]) / fg["ci"] + fg["pivot"] * fg["c"] + RAIL_W
    uu = np.arange(0.0, u_max + 1e-9, du)[:, None]
    ww = np.linspace(-RAIL_W / 2, RAIL_W / 2, nw)[None, :]
    up = uu - fg["pivot"] * fg["c"]
    x = fg["x_p"] + up * fg["ci"] - ww * fg["si"]
    lat = up * fg["si"] + ww * fg["ci"]
    keep = x <= fg["x_end"]
    return x[keep], lat[keep]


def fin_pocket_clearance(fg: dict, pk: Pockets, thetas) -> np.ndarray:
    """Arc clearance (mm, at the hull radius) between the yawed rail pads of the fin pair
    (theta, theta + 180) and the nearest stowed-blade pocket edge, for each theta; negative =
    overlap; inf where no pad point has a pocket at its x."""
    thetas = np.atleast_1d(np.asarray(thetas, float))
    x, lat = pad_footprint(fg)
    has = pk.depth(x) > 0
    if not has.any():
        return np.full(thetas.shape, np.inf)
    x, lat = x[has], lat[has]
    R = np.asarray(r_out(x), float)
    d = np.degrees(lat / R)
    lo, hi = pk.theta_span(x)
    out = np.full(thetas.shape, np.inf)
    for off in (0.0, 180.0):
        th = thetas[:, None] + off + d[None, :]
        for tp in pk.theta:
            a, b = tp + lo, tp + hi
            mid = 0.5 * (a + b)
            thw = mid + ((th - mid + 180.0) % 360.0 - 180.0)
            gap = np.where(thw < a, a - thw, np.where(thw > b, thw - b, -np.minimum(thw - a, b - thw)))
            out = np.minimum(out, (np.radians(gap) * R).min(axis=1))
    return out


_FOLD_CLOUD: dict = {}
_FOLD_NEAR: dict = {}


def blade_fold_cloud() -> dict:
    """World points of the three blades folding from deployed to stowed at the stowing index
    (blade i at psi = 120 i), from the rotor half's surface (layout.solve_axis +
    blade_surface, about 4 s). Returns x, rho, theta (theta from +Z toward +Y)."""
    if not _FOLD_CLOUD:
        blade = RL.load_blade()
        lay = RL.chosen_layout()
        ax = RL.solve_axis(blade, lay, n=41, iters=3)
        S = RL.blade_surface(blade, lay, ax)
        G = np.vstack([RL.densify(S.aero, 1.0).reshape(-1, 3), RL.densify(S.shank, 1.0).reshape(-1, 3),
                       S.knuckle.reshape(-1, 3), S.heel])
        pts = []
        x_keep = P.EQUATOR_X + 10.0          # the fins live on the front shell: keep x <= x_eq + 10
        for psi in (lay.psi(i) for i in range(lay.n_blades)):
            for f in np.arange(0.0, lay.fold_deg + 1e-9, 1.0):
                Wf = RL.apply(RL.blade_to_world(float(f), psi, lay), G)
                pts.append(Wf[Wf[:, 0] <= x_keep].astype(np.float32))
        W = np.vstack(pts).astype(float)
        _FOLD_CLOUD.update(x=W[:, 0], rho=np.hypot(W[:, 1], W[:, 2]),
                           theta=np.degrees(np.arctan2(W[:, 1], W[:, 2])) % 360.0,
                           n_points=int(W.shape[0]), fold_step_deg=1.0, fold_deg=lay.fold_deg,
                           x_keep_mm=x_keep, sampling="aero and shank surfaces densified to 1 mm, "
                                                      "fold in 1 deg steps (clearances +- ~0.5 mm)")
    return _FOLD_CLOUD


def _fold_cloud_near() -> dict:
    """The fold cloud restricted to where a fin pair in FIN_THETA_SEARCH can be (theta within 50
    deg of the search window or its opposite, x <= x_eq + 2)."""
    if not _FOLD_NEAR:
        C = blade_fold_cloud()
        lo, hi = FIN_THETA_SEARCH[0] - 50.0, FIN_THETA_SEARCH[1] + 50.0
        th = C["theta"]
        m = (C["x"] <= P.EQUATOR_X + 2.0) & (((th >= lo) & (th <= hi)) | ((th >= lo + 180) & (th <= hi + 180)))
        _FOLD_NEAR.update(x=C["x"][m], rho=C["rho"][m], theta=th[m], n_points=int(m.sum()),
                          fold_deg=C["fold_deg"], fold_step_deg=C["fold_step_deg"], sampling=C["sampling"])
    return _FOLD_NEAR


def fin_fold_clearance(fg: dict, span: float, fin_th) -> dict:
    """Smallest 3D distance between the blades folding at the stowing index and (a) the yawed rail
    pads, (b) the yawed fins (Clark Y envelope |w| <= 0.0585 c, from the root line to the tip),
    for the fin angles fin_th. Distances are taken to the pad / fin boxes in the yawed frame."""
    C = _fold_cloud_near()
    t_half = FIN_T_HALF * fg["c"]
    pad, fin = [], []
    for th_f in fin_th:
        dth = np.radians((C["theta"] - th_f + 180.0) % 360.0 - 180.0)
        near = np.abs(dth) < math.pi / 2
        x, rho, d_ = C["x"][near], C["rho"][near], dth[near]
        lat, al = rho * np.sin(d_), rho * np.cos(d_)
        u, w = _to_fin(fg, x, lat)
        top = _pad_top(fg, x)
        dp = np.sqrt(np.maximum(0.0, -u) ** 2 + np.maximum(0.0, np.abs(w) - RAIL_W / 2) ** 2
                     + np.maximum(0.0, x - fg["x_end"]) ** 2 + np.maximum(0.0, al - top) ** 2)
        y_r = _root_y(fg, x)
        df = np.sqrt(np.maximum(0.0, np.maximum(-u, u - fg["c"])) ** 2 + np.maximum(0.0, np.abs(w) - t_half) ** 2
                     + np.maximum(0.0, np.maximum(y_r - al, al - y_r - span)) ** 2)
        pad.append(float(dp.min()) if dp.size else math.inf)
        fin.append(float(df.min()) if df.size else math.inf)
    return dict(pad_clear_min_mm=min(pad), fin_clear_min_mm=min(fin),
                per_fin=[dict(theta_deg=float(t), pad_clear_mm=p, fin_clear_mm=f_)
                         for t, p, f_ in zip(fin_th, pad, fin)])


def fin_place(ri: dict, pk: Pockets, c: float, inc_deg: float, span: float, x_te: float,
              fold_check: bool = True, x_eq: float = P.EQUATOR_X) -> dict:
    """Fin angle theta (the pair at theta, theta + 180) and yaw pivot that maximise the smallest of
    the pad-to-pocket arc clearance and (fold_check) the pad / fin clearances to the blades
    folding at the stowing index. Exhaustive over the pocket grid (theta in FIN_THETA_SEARCH,
    pivot 0..1 in FIN_PIVOT_STEP), the fold check pruned by the pocket clearance (an upper bound
    of the minimum)."""
    sense = fin_te_sense(ri)
    thetas = np.arange(FIN_THETA_SEARCH[0], FIN_THETA_SEARCH[1] + 1e-9, FIN_THETA_STEP)
    cand = []
    for f in np.arange(0.0, 1.0 + 1e-9, FIN_PIVOT_STEP):
        fg = fin_frame(x_te, c, inc_deg, float(f), sense, x_eq)
        for th, v in zip(thetas, fin_pocket_clearance(fg, pk, thetas)):
            cand.append((float(v), round(float(f), 3), float(th)))
    cand.sort(key=lambda t: (-t[0], abs(t[2] - 90.0), abs(t[1] - 0.5)))
    best = None
    for n, (v, f, th) in enumerate(cand):
        if best is not None and (v <= best["min_clear_mm"] or n >= FIN_PLACE_MAX_EVAL):
            break
        fg = fin_frame(x_te, c, inc_deg, f, sense, x_eq)
        fc = fin_fold_clearance(fg, span, (th, th + 180.0)) if fold_check else None
        m = min(v, fc["pad_clear_min_mm"], fc["fin_clear_min_mm"]) if fc else v
        if best is None or m > best["min_clear_mm"] + 1e-9:
            best = dict(min_clear_mm=m, theta=th, pivot=f, pocket_clear_mm=v, fold=fc, fg=fg, n_eval=n + 1)
        if not fold_check:
            break
    best["ok_pocket"] = best["pocket_clear_mm"] >= FIN_POCKET_CLEAR
    best["ok_fold"] = (not fold_check) or min(best["fold"]["pad_clear_min_mm"],
                                              best["fold"]["fin_clear_min_mm"]) >= FIN_FOLD_CLEAR
    best["ok"] = best["ok_pocket"] and best["ok_fold"]
    return best


def _clarky_area_ratio() -> float:
    try:
        import aerosandbox as asb
        return float(asb.Airfoil(FIN_SECTION).area())
    except Exception:
        return 0.078                                       # ESTIMATE if aerosandbox is missing


def _fin_eval(c: float, s: float, r_cp: float, alpha_set: float, x_mid: float, x_rotor: float,
              ri: dict) -> dict:
    res = {}
    for v in FIN_SPEEDS:
        vf = upstream_factor(x_mid, x_rotor, ri["r_tip"]) * v
        pol = fin_polar(vf * c / 1000 / NU)
        a0 = pol["a0_per_rad"]
        AR = s / c
        ar = a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))
        CLmax3 = FIN_CLMAX_3D_FACTOR * pol["CLmax"]
        CL = ar * math.radians(alpha_set - pol["alpha_L0_deg"])
        q = 0.5 * RHO * vf ** 2
        S = c * s / 1e6
        tau = 2 * q * S * CL * r_cp / 1000
        Q = roll_torque(v, ri)["Q_Nm"]
        res[v] = dict(v_fin_m_s=vf, Re=pol["Re"], CL=CL, CLmax3D=CLmax3, CL_alpha_3D_per_rad=ar,
                      alpha_L0_deg=pol["alpha_L0_deg"], torque_Nm=tau, Q_gen_Nm=Q, ratio=tau / Q,
                      ratio_per_deg=2 * q * S * ar * math.radians(1.0) * r_cp / 1000 / Q,
                      authority_ratio=2 * q * S * FIN_CL_RESERVE * CLmax3 * r_cp / 1000 / Q,
                      stall_reserve=CLmax3 / CL if CL > 0 else float("inf"),
                      roll_damping_Nm_s=2 * q * S * ar * (r_cp / 1000) ** 2 / vf,
                      nf_confidence_min=float(pol["conf"].min()))
    return res


def fin_design(ri: dict, x_rotor: float, pk: Pockets, rule: str = FIN_FOLD_RULE,
               x_eq: float = P.EQUATOR_X, compare: bool = True) -> dict:
    """Fin pair with the smallest tip-to-tip size whose fixed incidence gives, at the design flow
    FIN_DESIGN_V, a couple of (1 + FIN_TORQUE_MARGIN) to (1 + FIN_TORQUE_MARGIN)(1 + FIN_RATIO_BAND)
    times the recommended generator's load-line torque (roll_torque), with CL <= FIN_CL_RESERVE x
    CLmax_3D at every flow in FIN_SPEEDS, TE at the rule's limit
    (fin_te_limit), LE on the nose, and whose YAWED rail pads (fin_frame: pad yawed with the fin
    about a pivot) clear the asymmetric pockets by FIN_POCKET_CLEAR and (rule 'index') whose pads
    and fins clear the blades folding at the stowing index by FIN_FOLD_CLEAR; theta and the pivot
    are chosen per fin by fin_place. The root line (and r_cp) use the yawed chord's x extent."""
    rows, res_of = [], {}
    for c in FIN_CHORDS:
        x_te, _ = fin_te_limit(ri, c, rule, x_eq)
        if x_te - c < 4.0:
            continue
        for s in np.arange(40.0, 130.01, 5.0):
            ok, alpha_ok, res, rl = False, None, {}, None
            target = 1.0 + FIN_TORQUE_MARGIN
            for alpha_set in np.arange(0.0, 14.0 + 1e-9, FIN_ALPHA_STEP):
                rl = fin_root_line(x_te - c * math.cos(math.radians(alpha_set)), x_te)
                r_cp = rl["y_mean"] + FIN_CP_SPAN * s
                res = _fin_eval(c, s, r_cp, alpha_set, x_te - c / 2, x_rotor, ri)
                if res[FIN_DESIGN_V]["ratio"] >= target:
                    ok = (all(m["CL"] <= FIN_CL_RESERVE * m["CLmax3D"] for m in res.values())
                          and res[FIN_DESIGN_V]["ratio"] <= target * (1.0 + FIN_RATIO_BAND))
                    alpha_ok = float(alpha_set)
                    break
            row = dict(chord_mm=c, span_mm=float(s), AR=s / c, area_cm2=c * s / 100,
                       r_cp_mm=rl["y_mean"] + FIN_CP_SPAN * s, x_le_mm=rl["x_le"], x_te_mm=x_te,
                       alpha_set_deg=alpha_ok, torque_ok=bool(ok), tip_to_tip_mm=2 * (rl["y_te"] + s),
                       theta_deg=None, pivot_chord_frac=None, clear_min_mm=None, feasible=None)
            rows.append(row)
            res_of[(c, float(s))] = (res, rl)

    placed = {}

    def place(row):
        key = (row["chord_mm"], row["span_mm"])
        if key not in placed:
            b = fin_place(ri, pk, row["chord_mm"], row["alpha_set_deg"], row["span_mm"], row["x_te_mm"],
                          fold_check=(rule == "index"), x_eq=x_eq)
            placed[key] = b
            row.update(theta_deg=b["theta"], pivot_chord_frac=b["pivot"], clear_min_mm=b["min_clear_mm"],
                       feasible=bool(b["ok"]))
        return placed[key]

    cand = [r_ for r_ in rows if r_["torque_ok"]]
    pick = next((r_ for r_ in sorted(cand, key=lambda r_: (round(r_["tip_to_tip_mm"], 3), r_["area_cm2"]))
                 if place(r_)["ok"]), None)
    if pick is None:
        raise ValueError("no feasible fin for rule %s" % rule)
    min_area = next((r_ for r_ in sorted(cand, key=lambda r_: (r_["area_cm2"], r_["tip_to_tip_mm"]))
                     if place(r_)["ok"]), None)
    b = placed[(pick["chord_mm"], pick["span_mm"])]
    fg = b["fg"]
    res, rl = res_of[(pick["chord_mm"], pick["span_mm"])]
    phi = b["theta"]
    alt = None
    if compare:
        o = "torus" if rule == "index" else "index"
        try:
            a_ = fin_design(ri, x_rotor, pk, rule=o, x_eq=x_eq, compare=False)
            alt = dict(rule=o, chord_mm=a_["chord_mm"], span_mm=a_["span_mm"], x_le_mm=a_["x_le_mm"],
                       x_te_mm=a_["x_te_mm"], theta_deg=a_["theta_deg"], incidence_deg=a_["incidence_deg"],
                       pivot_chord_frac=a_["yaw"]["pivot_chord_frac"], tip_to_tip_mm=a_["tip_to_tip_mm"],
                       area_each_cm2=a_["area_each_cm2"], drag_tow_N=a_["drag_tow_N"])
        except (TypeError, ValueError):
            alt = dict(rule=o, feasible=False)
    c, s, x_le, x_te = pick["chord_mm"], pick["span_mm"], fg["x_le"], fg["x_te"]
    inc = pick["alpha_set_deg"]
    # tow case (2.5 m/s): drag for the nose load and root bending
    vt = upstream_factor(x_te - c / 2, x_rotor, ri["r_tip"]) * V_TOW
    pol = fin_polar(vt * c / 1000 / NU)
    a0 = pol["a0_per_rad"]
    AR = s / c
    ar = a0 / (math.sqrt(1 + (a0 / (math.pi * AR)) ** 2) + a0 / (math.pi * AR))
    CL = ar * math.radians(inc - pol["alpha_L0_deg"])
    cd0 = float(np.interp(inc, pol["alpha"], pol["CD"]))
    CD = cd0 + CL ** 2 / (math.pi * FIN_E_OSWALD * AR)
    D2 = 2 * 0.5 * RHO * vt ** 2 * (c * s / 1e6) * CD
    L1 = 0.5 * RHO * vt ** 2 * (c * s / 1e6) * CL
    M_root = L1 * FIN_CP_SPAN * s / 1000
    t = 0.117 * c
    Z = 0.6 * c * t ** 2 / 6                              # ESTIMATE, airfoil ~ 0.6 x rectangle
    src = sources()
    sxy = (src.get("prusament_petg_tds", {}).get("tensile_yield_xy_MPa") or [47])[0]
    x_cl = ri["hinge_x"] - math.sqrt(max((blade_sweep_len(ri) + FIN_TIP_CLEAR) ** 2
                                         - (rl["y_te"] - ri["hinge_r"]) ** 2, 0.0))
    # pad footprint and volumes (per pad): outer volume above the hull, the dovetail groove (filled
    # by the fin's root block when the fin is on)
    xq = np.linspace(x_le, fg["x_end"], 200)
    h_q = np.maximum(np.asarray(_pad_top(fg, xq)) - np.asarray(r_out(xq)), 0.0)
    pad_V = RAIL_W / fg["ci"] * float(np.trapezoid(h_q, xq))
    groove_V = GROOVE_W * RAIL_H * c
    xs_fp, lat_fp = pad_footprint(fg)
    pocket_cl = b["pocket_clear_mm"]
    fold_c = b["fold"]
    off_te, off_le, off_end = (float(fin_offset(fg, v)) for v in (x_te, x_le, fg["x_end"]))
    qd = roll_torque(FIN_DESIGN_V, ri)
    md = res[FIN_DESIGN_V]
    far = max(b["pivot"], 1.0 - b["pivot"]) * c                # chord length from the pivot to the far end
    torque_basis = dict(
        generator=gen_load_line().get("generator"), basis=qd["basis"], design_flow_m_s=FIN_DESIGN_V,
        Q_design_Nm=qd["Q_Nm"], tsr=qd.get("tsr"), rpm=qd.get("rpm"), margin_pct=100 * FIN_TORQUE_MARGIN,
        target_ratio=1.0 + FIN_TORQUE_MARGIN, target_ratio_max=(1.0 + FIN_TORQUE_MARGIN) * (1.0 + FIN_RATIO_BAND),
        fin_torque_Nm=md["torque_Nm"], ratio=md["ratio"], Q_cp_rule_Nm=qd["Q_cp_rule_Nm"],
        source=gen_load_line().get("source"),
        margin_reason="+%.0f %% over the load-line torque: NeuralFoil polars are smooth-surface, and a printed "
                      "Clark Y at Re ~%.0f k tends to give less lift; the GL40 drag torque (%.4f N m in the BEM "
                      "model) is an estimate. Lift shortfall is the likelier error, so the bias is deliberate; "
                      "if the prototype over-balances, a trim fin with less incidence corrects it"
                      % (100 * FIN_TORQUE_MARGIN, md["Re"] / 1000,
                         next((r_["drag_Nm"] for r_ in ((_load(BLADE_JSON) or {}).get("generator_coupling_estimate")
                                                        or {}).get("rows", []) if r_.get("role") == "recommended"),
                              float("nan"))),
        old_rule_note="the merge sized the fins for Q = Cp %.2f at TSR %g (%.3f N m at %.1f m/s), %.2f x the "
                      "load-line torque" % (CP_DESIGN, P.DESIGN_TSR, qd["Q_cp_rule_Nm"], FIN_DESIGN_V,
                                            qd["Q_cp_rule_Nm"] / qd["Q_Nm"]))
    trim = dict(
        variant_deg=FIN_TRIM_DEG, ratio_change_per_deg_at_design_flow=md["ratio_per_deg"],
        torque_change_per_deg_Nm=md["ratio_per_deg"] * md["Q_gen_Nm"],
        fin_box_shift_max_mm=far * math.sin(math.radians(FIN_TRIM_DEG)),
        how="print fin variants whose root tongue is turned +-%g deg about the yaw pivot relative to the chord: "
            "the shell's yawed pad and groove stay, the incidence changes by +-%g deg (%.0f %% of the "
            "torque per degree); the fin box moves by at most %.2f mm at its far end, inside the fold "
            "clearance" % (FIN_TRIM_DEG, FIN_TRIM_DEG, 100 * md["ratio_per_deg"],
                           far * math.sin(math.radians(FIN_TRIM_DEG))))
    return dict(
        section=FIN_SECTION + " (flat lower surface, t/c 0.117), printed flat, no supports",
        polars="NeuralFoil large, n_crit 5, cached in out/work/layout_internal/polars_fin.json",
        n=2, theta_deg=[phi, phi + 180.0], fold_rule=rule,
        layout="identical parts at theta %g and %g (diametrically opposite, %g deg off +-Y so the "
               "yawed rail pads clear the asymmetric pockets); each fin and its rail pad are yawed by the "
               "incidence about the radial line through the pivot at %.2f c, TE toward -theta; the second "
               "fin is the first turned 180 deg about X, so the lifts are opposite and form a pure couple"
               % (phi, phi + 180.0, abs(90.0 - phi), b["pivot"]),
        torque_sense=("rotor omega %s: the generator reaction turns the hull the same way (+X), "
                      "so the +Y fin lifts toward -Z and the -Y fin toward +Z (fin couple -X)"
                      % ri.get("omega_vector", "+X")),
        chord_mm=c, span_mm=s, AR=AR, area_each_cm2=c * s / 100,
        root=rl, r_root_mean_mm=rl["y_mean"], r_tip_max_mm=rl["y_te"] + s, r_cp_mm=pick["r_cp_mm"],
        tip_to_tip_mm=2 * (rl["y_te"] + s),
        x_le_mm=x_le, x_te_mm=x_te, incidence_deg=inc, margins=res, torque_basis=torque_basis, trim=trim,
        yaw=dict(pivot_chord_frac=b["pivot"], pivot_x_mm=fg["x_p"], te_side="-theta" if fg["sense"] < 0 else "+theta",
                 offset_per_mm_x=fg["si"] / fg["ci"], pad_half_width_at_constant_x_mm=RAIL_W / 2 / fg["ci"],
                 offset_at_le_mm=off_le, offset_at_te_mm=off_te, offset_at_pad_end_mm=off_end,
                 offset_sign="+ = toward +theta, at the hull, from the fin angle line",
                 x_extent_note="the yawed root chord spans x %.2f-%.2f (c cos i)" % (x_le, x_te)),
        drag_tow_N=D2, lift_tow_each_N=L1, root_moment_tow_Nm=M_root,
        root_stress_tow_MPa=M_root * 1000 / Z, root_util_xy=M_root * 1000 / Z / sxy,
        root_note="lift bends the fin about its chord; printed flat, the span-wise bending "
                  "stress lies in the layer plane (XY strength)",
        rail=dict(type="half dovetail (flat on the fin's flat side, 45 deg flank on the other): male tongue "
                       "on the fin root, female groove (mean width %g, depth %g) in a %g mm raised pad on the front "
                       "shell; the pad is yawed with the fin, so the rail runs along the fin chord"
                       % (GROOVE_W, RAIL_H, RAIL_W),
                  x_mm=[x_le, x_te], pad_x_mm=[float(xs_fp.min()), fg["x_end"]], pad_w_mm=RAIL_W,
                  pad_h_mm=[rl["pad_min_mm"], float(_pad_top(fg, fg["x_end"]) - r_out(fg["x_end"]))
                            if not fg["chamfer"] else rl["pad_max_mm"]],
                  pad_aft_end=("45 deg chamfer from the TE to x %.1f" % fg["x_end"]) if fg["chamfer"] else
                              ("runs down to the equator face (x %.1f): the front shell prints on that face, so "
                               "the pad's aft end stands on the bed (a 45 deg chamfer from the rail stop would "
                               "need %.1f mm, past the equator)" % (fg["x_end"], fg["h_te"])),
                  inclination_deg=rl["incl_deg"], pad_outer_volume_each_cc=pad_V / 1000,
                  groove_volume_each_cc=groove_V / 1000,
                  retention="the fin slides aft onto the end stop at x %.1f (drag pushes it home); one M3 grub "
                            "screw in a heat-set insert for handling" % x_te),
        clearance=dict(blade_sweep_len_mm=blade_sweep_len(ri), torus_te_limit_x_mm=x_cl,
                       tip_clear_mm=FIN_TIP_CLEAR,
                       stowed_pockets=dict(min_clear_mm=pocket_cl, ok=bool(b["ok_pocket"]), limit_mm=FIN_POCKET_CLEAR,
                                           model="yawed pad footprint (pad_footprint, 0.25 mm x 11 samples) vs "
                                                 "each pocket's theta span at the sample x, arc at the hull radius"),
                       fold_at_index=(dict(rule="blades fold at the stowing index (blade i at psi 120 i), fold 0 to "
                                                "%.0f deg in %.0f deg steps, rotor surface from hydroegg.layout"
                                                % (_fold_cloud_near()["fold_deg"], _fold_cloud_near()["fold_step_deg"]),
                                           sampling=_fold_cloud_near()["sampling"],
                                           n_points_near_fins=_fold_cloud_near()["n_points"],
                                           pad_clear_min_mm=fold_c["pad_clear_min_mm"],
                                           fin_clear_min_mm=fold_c["fin_clear_min_mm"], per_fin=fold_c["per_fin"],
                                           fin_half_thickness_mm=FIN_T_HALF * c, limit_mm=FIN_FOLD_CLEAR,
                                           model="3D distance to the yawed pad box (|w| <= %g, up to the pad top) "
                                                 "and fin box (|w| <= t/2, root line to tip) in the yawed frame"
                                                 % (RAIL_W / 2),
                                           ok=bool(b["ok_fold"])) if fold_c else None),
                       search=dict(theta_range_deg=list(FIN_THETA_SEARCH), theta_step_deg=FIN_THETA_STEP,
                                   pivot_step=FIN_PIVOT_STEP, evaluated_fold_checks=b["n_eval"]),
                       te_rule=("index: TE at the rail end stop %.1f mm ahead of the equator; the blades "
                                "fold only at the stowing index (the pockets need it) and clear the fins "
                                "by >= %.1f mm there (fold_at_index); folding at another rotor azimuth "
                                "with the fins on can touch a fin (fins come off for transport)"
                                % (FIN_TE_EQ_MARGIN, FIN_FOLD_CLEAR)) if rule == "index" else
                               "torus: TE root point >= sweep length + clearance from the hinge circle",
                       rule="yawed rail pad >= %.1f mm (arc) from any pocket edge" % FIN_POCKET_CLEAR),
        alternative_rule=alt,
        selection="smallest fin-tip diameter (handling, stowage) among fins whose couple at %.1f m/s is "
                  "%.2f-%.3f x the %s load-line torque (CL <= %.2f CLmax_3D at %s m/s) AND whose yawed pads / "
                  "fins clear the pockets and the blades folding at the stowing index; candidates evaluated in "
                  "that order; the least-area feasible fin is listed for comparison"
                  % (FIN_DESIGN_V, 1 + FIN_TORQUE_MARGIN, (1 + FIN_TORQUE_MARGIN) * (1 + FIN_RATIO_BAND),
                     gen_load_line().get("short") or "generator", FIN_CL_RESERVE,
                     "/".join("%g" % v for v in FIN_SPEEDS)),
        least_area_alternative=min_area,
        sweep=rows,
    )


def roll_summary(ri: dict, fins: dict, mb: dict) -> dict:
    out = {}
    for v in (1.0, 1.5, 2.5):
        out[f"{v:g}_m_s"] = roll_torque(v, ri)
    W = mb["v2"]["displacement_N"]
    BG = mb["v2"]["BG_mm"] / 1000
    Qd = out["1.5_m_s"]["Q_Nm"]
    bal = {}
    for v, m in fins["margins"].items():
        dT = m["torque_Nm"] - m["Q_gen_Nm"]
        p = dT / m["roll_damping_Nm_s"]
        T_v = 0.5 * RHO * math.pi * (ri["r_tip"] / 1000) ** 2 * v ** 2 * CT_TOW
        k_br = T_v * (BRIDLE_A / 1000) ** 2 / (BRIDLE_L / 1000)
        bal[f"{v:g}_m_s"] = dict(excess_torque_Nm=dT, excess_pct=100 * dT / m["Q_gen_Nm"],
                                 free_roll_rate_rpm=p * 60 / (2 * math.pi),
                                 pendulum_heel_deg=(math.degrees(math.asin(min(1.0, abs(dT) / (W * BG))))
                                                    if W * BG > 0 else 90.0),
                                 bridle_stiffness_Nm_rad=k_br,
                                 bridle_heel_deg=math.degrees(abs(dT) / k_br))
    q10, q15, q25 = (out[k]["Q_Nm"] for k in ("1_m_s", "1.5_m_s", "2.5_m_s"))
    m15 = fins["margins"][FIN_DESIGN_V]
    p_free = m15["torque_Nm"] / m15["roll_damping_Nm_s"]
    return dict(
        torque=out,
        Q_scaling="on the %s load line the operating TSR rises with flow (%s), so Q does not grow with v^2 "
                  "like the fin lift: Q(1.0) / Q(1.5) = %.3f against (1.0 / 1.5)^2 = %.3f, Q(2.5) / Q(1.5) = %.3f "
                  "against %.3f (the 2.5 m/s point is held by the cavitation cap). A fixed fin set matched at "
                  "1.5 m/s therefore under-balances at 1.0 m/s and over-balances at 2.5 m/s (balance below)"
                  % (gen_load_line().get("short") or "generator",
                     ", ".join("TSR %.2f at %g m/s" % (out[k].get("tsr") or float("nan"), out[k]["v_m_s"])
                               for k in ("1_m_s", "1.5_m_s", "2.5_m_s")),
                     q10 / q15, (1.0 / 1.5) ** 2, q25 / q15, (2.5 / 1.5) ** 2),
        balance=bal,
        pendulum_righting_max_Nm=W * BG, pendulum_vs_Qdesign=W * BG / Qd,
        bridle=dict(eyes="two eyes at +-Y on the nose", a_mm=BRIDLE_A, leg_mm=BRIDLE_L,
                    model="k = T a^2 / L (legs carry T/2 each)"),
        conclusion="fins carry the mean torque; nothing on an axial tow line restores roll, and the "
                   "pendulum (near-neutral hull) is negligible, so any mismatch turns into a slow "
                   "hull roll that twists the cable. v2 should add a roll reference: a two-leg "
                   "bridle on two nose eyes (numbers above) or a cable swivel/slip ring",
        freewheel_roll_rate_1p5_m_s_rpm=p_free * 60 / (2 * math.pi),
        freewheel_note="with no electrical load Q ~ 0 and the fins roll the hull the other way at "
                       "p ~ tau_fin / damping (%.0f rpm at 1.5 m/s); v2 must keep a load on the generator "
                       "(the load line above, a dump load when the battery is full)" % (p_free * 60 / (2 * math.pi)),
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


def _thick_feature(V, A_exposed):
    """Thick ring inside a shell: the shell's 6-wall skin on its exposed faces, infill inside."""
    skin = min(V, A_exposed * P.HULL_WALL)
    return skin + INFILL * (V - skin)


def hub_wet(ri: dict) -> dict:
    """Volumes of the rotating tail (hub + spinner, x from the hub front to the tip) that flood,
    so the displacement does not count them: (1) the three blade slots through hub and spinner,
    from the rotor half's hub_slot table (per x: slot floor r_e along e_r and half-width along the
    pin), cut straight out through the hull envelope; (2) the screw-access bore on the axis from
    the spinner screw head to the tip (the spinner aft of the hub face is printed solid around it,
    perimeters + infill, closed); (3) the gap between the spinner shell and the hub core. The
    blade roots inside the slots are counted with the blades."""
    rj = _load(ROTOR_JSON) or {}
    hub, slot, ch = rj.get("hub", {}), rj.get("hub_slot", {}), rj.get("chosen", {})
    x_front = float(ch.get("hub_front_x_mm", ri["static_tail_x"] + P.HUB_GAP))
    x_aft = float(hub.get("hub_aft_face_x_mm", ri["shaft_aft_x"]))
    V_slot, A_wall = 0.0, 0.0
    if slot.get("x_mm"):
        xs = np.asarray(slot["x_mm"], float)
        fl = np.asarray(slot["slot_floor_re_mm"], float)
        hwid = np.asarray(slot["slot_half_width_mm"], float)
        dx = float(np.median(np.diff(xs))) if xs.size > 1 else 0.5
        for x_, f_, h_ in zip(xs, fl, hwid):
            R = float(r_out(x_))
            tt = np.linspace(-h_, h_, 41)
            hh = np.maximum(np.sqrt(np.maximum(R ** 2 - tt ** 2, 0.0)) - f_, 0.0)
            V_slot += float(np.trapezoid(hh, tt)) * dx
            A_wall += 2 * max(R - f_, 0.0) * dx
    sp = hub.get("spinner", {})
    gap = float(sp.get("radial_gap_to_hub_mm", 0.2))
    wall = float(sp.get("wall_mm", P.MIN_WALL))
    x_acc0 = x_aft + SPINNER_AFT_WALL + M3_HEAD_H
    V_acc = math.pi * (SPINNER_ACCESS_D / 2) ** 2 * max(P.HULL_L - x_acc0, 0.0)
    xg = np.linspace(x_front, x_aft, 200)
    V_gap = float(np.trapezoid(2 * math.pi * (np.asarray(r_out(xg)) - wall - gap / 2) * gap, xg))
    return dict(x_front_mm=x_front, x_aft_mm=x_aft, slots_cc=3 * V_slot / 1000, slot_each_cc=V_slot / 1000,
                slot_wall_area_mm2=A_wall, access_bore_cc=V_acc / 1000, access_bore_x_mm=[x_acc0, P.HULL_L],
                access_bore_d_mm=SPINNER_ACCESS_D, gap_film_cc=V_gap / 1000,
                spinner_tip="aft of the hub face the spinner is printed solid (closed perimeters + %d %% "
                            "gyroid) around a D%g screw-access bore on the axis, which floods (a push-in plug "
                            "would add %.1f cc of displacement)" % (int(INFILL * 100), SPINNER_ACCESS_D,
                                                                     V_acc / 1000),
                source="data/layout_rotor.json hub_slot + hub.spinner")


def mass_buoyancy(J, bay, DT, fins, ri, hw, pk: Pockets, nose_d) -> dict:
    x_eq = J["x_eq"]
    x_end = ri["static_tail_x"]
    parts = []

    def add(name, group, V_print_mm3, x_c, z_c=0.0, mat="print", V_disp_mm3=None, note=""):
        parts.append(dict(name=name, group=group, V_print_mm3=V_print_mm3, x_mm=x_c, z_mm=z_c,
                          mat=mat, V_disp_mm3=V_disp_mm3, note=note))

    # ---- front shell
    v, xc = _rev_vol(0, x_eq, r_out, r_in)
    add("front shell wall", "F", v, xc)
    x_r0 = J["x_ring0"] - (float(r_in(J["x_ring0"])) - J["r_lip"])
    v, xc = _rev_vol(x_r0, x_eq, r_in,
                     lambda x: np.where(x >= J["x_seal0"], J["r_bore"],
                                        np.where(x >= J["x_ring0"], J["r_lip"],
                                                 np.minimum(r_in(x), J["r_lip"] + (J["x_ring0"] - x)))))
    A_exp = 2 * math.pi * J["r_lip"] * (x_eq - x_r0) + math.pi * (float(r_in(x_eq)) ** 2 - J["r_bore"] ** 2)
    add("socket ring (6-wall skin + infill)", "F", _thick_feature(v, A_exp), xc)
    Rm = EYE_HOLE_D / 2 + EYE_W / 2
    add("tow eye", "F", 2 * math.pi * Rm * EYE_W * EYE_T, nose_d["eye"]["hole_centre_x_mm"])
    npl = nose_d["nose_plug"]
    Vp = npl["plug_cavity_cc"] * 1000
    add("nose plug (infill, bridged)", "F", _solid_part(Vp, math.pi * npl["r_inner_at_45_mm"] ** 2),
        0.6 * npl["x_45_mm"])
    gl = nose_d["gland"]
    add("gland boss (seat + thread boss)", "F", math.pi / 4 * ((GLAND_SEAT_D ** 2 - 12 ** 2) * GLAND_SEAT_T
                                                              + (GLAND_BOSS_D ** 2 - 12 ** 2)
                                                              * (gl["thread_engagement_mm"] - GLAND_SEAT_T)),
        gl["seat_x_mm"], -gl["seat_surface_r_mm"] + 4)
    xcar = bay["pick"]["x_fwd_mm"] - CARRIER_T / 2
    add("carrier ring", "F", 2 * math.pi * float(r_in(xcar)) * CARRIER_RING * 7.0, xcar)
    rail = fins["rail"]
    pad_net = (rail["pad_outer_volume_each_cc"] - rail["groove_volume_each_cc"]) * 1000
    add("fin rail pads x2 (outside the hull, minus the wet groove)", "F", 2 * pad_net,
        0.5 * sum(rail["pad_x_mm"]), V_disp_mm3=2 * pad_net)
    xp0 = max(pk.x0, 0)

    def pocket_wall_vol(xa, xb):
        """Extra print volume of the 3 pocket troughs (floor + side walls - removed shell)."""
        if xb <= xa:
            return 0.0, xa
        xq = np.linspace(xa, xb, 200)
        dq = pk.depth(xq)
        wq = pk.width_at(xq)
        a = np.where(dq > 0, wq * POCKET_FLOOR_WALL + 2 * dq * T_MIN - wq * P.HULL_WALL, 0.0)
        v = 3 * float(np.trapezoid(a, xq))
        return v, (float(np.trapezoid(a * xq, xq) / np.trapezoid(a, xq)) if abs(v) > 1e-9 else xa)

    add("pocket walls (F part)", "F", *pocket_wall_vol(xp0, x_eq))
    # ---- aft shell
    v, xc = _rev_vol(x_eq, x_end, r_out, r_in)
    add("aft shell wall", "A", v, xc)
    add("pocket walls (A part)", "A", *pocket_wall_vol(max(x_eq, pk.x0), min(x_end, pk.x1)))
    # tail end block: solid from the 688 seat's aft face to the end face (seal shoulder, seal
    # bore, the 3 notch floors and the 3 stowed-detent magnet seats of the rotor half)
    xb0 = DT["stations_mm"]["bearing_rear_688"][1]
    r_end = float(r_in(0.5 * (xb0 + x_end)))
    V_end = math.pi * (r_end ** 2 - (P.LIP_SEAL["D"] / 2) ** 2) * (x_end - xb0)
    add("tail end block (seal boss, notch floors, magnet seats)", "A",
        _thick_feature(V_end, 2 * math.pi * (r_end ** 2 - (P.LIP_SEAL["D"] / 2) ** 2)), 0.5 * (xb0 + x_end))
    vs = math.pi * ((J["r_spigot_seal"] ** 2 - J["r_spigot_in"] ** 2) * J["L_seal"]
                    + (J["r_spigot_lug"] ** 2 - J["r_spigot_in"] ** 2) * (LIP_T + LUG_T))
    vs += LUG_N * math.radians(LUG_W_DEG) * J["r_spigot_lug"] * LUG_H * LUG_T * 0.75
    add("spigot + lugs", "A", vs, x_eq - J["L_spigot"] / 2)
    v, xc = _rev_vol(x_eq, x_eq + J["root_ring_L"], r_in,
                     lambda x: np.minimum(J["r_spigot_in"] + (x - x_eq), r_in(x)))
    A_exp = (math.pi * (float(r_in(x_eq)) + J["r_spigot_in"]) * J["root_ring_L"] * math.sqrt(2)
             + math.pi * (float(r_in(x_eq)) ** 2 - J["r_spigot_seal"] ** 2))
    add("spigot root ring (6-wall skin + infill)", "A", _thick_feature(v, A_exp), xc)
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
    mo = hw["motor"]
    Dg, Lg = mo["D"], mo["L"]
    Vg = math.pi * Dg ** 2 / 4 * Lg
    xg = 0.5 * sum(bay["pick"]["motor_body_x_mm"])
    gen_replica = _solid_part(Vg, math.pi * Dg ** 2 / 2 + math.pi * Dg * Lg)
    Ls = bay["pick"]["adapter_spacer_L_mm"]
    if Ls > 0.5:
        r_o = min(mo["D"], bay["pick"]["D_mm"]) / 2 - 1.0
        add("motor adapter / spacer ring", "int", math.pi * (r_o ** 2 - (r_o - 2.4) ** 2) * Ls
            + math.pi * (r_o ** 2 - 4.5 ** 2) * 3.0, bay["pick"]["x_fwd_mm"] + Ls / 2, V_disp_mm3=0.0)
    # ---- fins (outside; displace their own volume)
    Vf = 2 * (_clarky_area_ratio() * fins["chord_mm"] ** 2 * fins["span_mm"]
              + rail["groove_volume_each_cc"] * 1000)       # blade + root tongue (fills the groove)
    Af = 2 * (2 * fins["chord_mm"] * fins["span_mm"])
    add("fins x2", "fin", _solid_part(Vf, Af), (fins["x_le_mm"] + fins["x_te_mm"]) / 2,
        V_disp_mm3=Vf)
    # ---- rotor parts: hub + spinner envelope minus what floods (hub_wet)
    hw_ = hub_wet(ri)
    xh0 = hw_["x_front_mm"]
    v_h, x_h = _rev_vol(xh0, P.HULL_L, r_out)
    sa_h = 2 * math.pi * np.trapezoid(r_out(np.linspace(xh0, P.HULL_L, 200)), np.linspace(xh0, P.HULL_L, 200))
    v_solid = v_h - (hw_["slots_cc"] + hw_["access_bore_cc"] + hw_["gap_film_cc"]) * 1000
    add("hub + spinner (envelope minus 3 blade slots, screw access bore, shell gap)", "rotor",
        _solid_part(v_solid, float(sa_h) + math.pi * 20 ** 2 + 3 * hw_["slot_wall_area_mm2"]), x_h,
        V_disp_mm3=v_solid)
    if ri.get("blade_cg_n_mm") is not None:          # rotor half: CG in the blade frame, n = +X deployed
        xb = ri["hinge_x"] + ri["blade_cg_n_mm"]
    else:
        xb = ri["hinge_x"] + (0.5 * (30 + ri["r_tip"]) - ri["hinge_r"]) * math.tan(math.radians(ri["cone_deg"]))
    if ri.get("blade_mass_solid_PLA_g"):
        Vb = ri["blade_mass_solid_PLA_g"] / RHO_PLA * 1000          # mm^3, printed solid
        src_b = "layout_rotor.json blade_mass" if ri.get("rotor") else "blade_design.json"
        add(f"blades x3 ({src_b}, solid incl. knuckle and heel)", "rotor", 3 * Vb, xb, V_disp_mm3=3 * Vb)
    else:
        rr = np.linspace(30.0, ri["r_tip"], 60)
        cc = np.interp(rr, [30, ri["r_tip"]], [26, 14])
        tc = np.interp(rr, [30, ri["r_tip"]], [0.18, 0.12])
        Vb = float(np.trapezoid(0.685 * tc * cc ** 2, rr)) + 22 * 8 * (30 - ri["hinge_r"])
        add("blades x3 (ESTIMATE)", "rotor", 3 * Vb, xb, V_disp_mm3=3 * Vb)

    # ---- hardware (grams directly)
    sh = DT["shaft"]
    rr_ = ri.get("rotor") or {}
    pin_l = rr_.get("pin_l", P.PIN_L)
    m_mag = RHO_NDFEB * math.pi * (P.MAGNET["D"] / 2) ** 2 * P.MAGNET["H"] / 1000
    hwl = [
        ("shaft 8 mm", sh["mass_g"], (sh["x_fwd"] + sh["x_aft"]) / 2, "computed, steel 7.75 g/cc ESTIMATE"),
        ("608-2RS", hw["b608"]["mass_g"], sum(DT["stations_mm"]["bearing_front_608"]) / 2, "catalogue"),
        ("688-2RS", hw["b688"]["mass_g"], sum(DT["stations_mm"]["bearing_rear_688"]) / 2, "catalogue, VERIFY"),
        ("lip seal 8x16x5", hw["seal"]["mass_g"], sum(DT["stations_mm"]["lip_seal"]) / 2, "ESTIMATE"),
        ("O-ring", RHO_NBR * math.pi ** 2 * (J["oring"]["id_free_mm"] + OR_CS) * OR_CS ** 2 / 4 / 1000,
         x_eq - 3, "computed, NBR 1.3 g/cc ESTIMATE"),
        ("magnets x12 (rotor detents: 6 blade, 3 hub stop, 3 tail notch)", 12 * m_mag,
         (6 * ri["hinge_x"] + 3 * (rr_.get("stop_block_x") or [ri["hinge_x"] - 3.7] * 2)[1]
          + 3 * rr_.get("stow_magnet_x", x_end - 3.5)) / 12, "computed, NdFeB 7.5 g/cc ESTIMATE"),
        ("magnets x4 (equator joint, 2 pairs)", 4 * m_mag, x_eq, "computed, NdFeB 7.5 g/cc ESTIMATE"),
        ("hinge pins 3x%g x3" % pin_l, 3 * RHO_STEEL * math.pi * 1.5 ** 2 * pin_l / 1000, ri["hinge_x"],
         "computed"),
        ("M3 screws, grub screws, inserts", 12 * 0.5, 40.0, "ESTIMATE"),
        ("M12 cable gland", hw["gland"]["mass_g"], nose_d["gland"]["seat_x_mm"], "ESTIMATE"),
    ]
    V_static = _rev_vol(0, x_end, r_out)[0]
    xq = np.linspace(pk.x0, min(pk.x1, x_end), 300)
    aq = pk.width_at(xq) * pk.depth(xq)
    V_pockets = 3 * float(np.trapezoid(aq, xq))
    x_pockets = float(np.trapezoid(aq * xq, xq) / np.trapezoid(aq, xq)) if V_pockets > 0 else 0.0
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
            g = hw["motor"]["mass_g"] + 10.0       # + 10 g leads and solder (ESTIMATE)
        rows.append((f"generator ({generator})", g, xg))
        m += g
        mx += g * xg
        if ballast_g > 0:
            xb_ = 0.5 * (DT["stations_mm"]["bearing_front_608"][1] + DT["stations_mm"]["bearing_rear_688"][0])
            zb_ = -(float(r_in(xb_)) - 4.0)                  # bottom sector (theta 180), no pocket
            rows.append(("lead ballast", ballast_g, xb_))
            m += ballast_g
            mx += ballast_g * xb_
            mz += ballast_g * zb_
        V_disp = V_static - V_pockets + sum(p["V_disp_mm3"] or 0.0 for p in parts
                                            if p["V_disp_mm3"] and (fins_on or p["group"] != "fin"))
        Vx = (V_static * x_static - V_pockets * x_pockets
              + sum((p["V_disp_mm3"] or 0.0) * p["x_mm"] for p in parts
                    if p["V_disp_mm3"] and (fins_on or p["group"] != "fin")))
        disp_g = V_disp * RHO / 1e6
        return dict(mass_g=m, x_cg_mm=mx / m, z_cg_mm=mz / m, displacement_cc=V_disp / 1000,
                    displacement_g=disp_g, x_cb_mm=Vx / V_disp, z_cb_mm=0.0,
                    net_buoyancy_N=(disp_g - m) * G / 1000, displacement_N=disp_g * G / 1000,
                    BG_mm=0.0 - mz / m, rows=rows)

    disp_pla = build(RHO_PLA, "replica", fins_on=False)
    v2_noball = build(RHO_PETG, "real", fins_on=True)
    # v2 target: net buoyancy +0.3 N (rises slowly in slack water, never sits on the bed)
    target_N = 0.3
    ballast = max(0.0, (v2_noball["net_buoyancy_N"] - target_N) * 1000 / G)
    v2 = build(RHO_PETG, "real", ballast_g=ballast, fins_on=True)
    xb_ = 0.5 * (DT["stations_mm"]["bearing_front_608"][1] + DT["stations_mm"]["bearing_rear_688"][0])
    lead_cc = ballast / RHO_LEAD
    # what v2 needs
    deficit_g = max(0.0, (target_N - v2_noball["net_buoyancy_N"]) * 1000 / G)
    shell_cc = sum(p["V_print_mm3"] for p in parts if p["name"] in ("front shell wall", "aft shell wall")) / 1000
    save_4p_g = shell_cc * (1 - 1.6 / P.HULL_WALL) * RHO_PETG
    foam_cc = deficit_g / (1 - RHO_FOAM)
    W = v2["displacement_N"]
    dx_cb = (v2["x_cg_mm"] - v2["x_cb_mm"]) / 1000
    T15 = 0.5 * RHO * math.pi * (ri["r_tip"] / 1000) ** 2 * 1.5 ** 2 * CT_TOW
    k_pitch = T15 * (DT["stations_mm"]["rotor_load_plane"] - 0.0) / 1000      # N m / rad, eye at x 0
    test_ballast_g = 20.0
    v2_test = build(RHO_PETG, "real", ballast_g=test_ballast_g, fins_on=True)
    needs = dict(
        target_net_buoyancy_N=target_N, net_buoyancy_no_ballast_N=v2_noball["net_buoyancy_N"],
        mass_to_shed_or_float_g=deficit_g,
        option_4_perimeter_shells_saves_g=save_4p_g,
        option_external_foam_cc=foam_cc, foam_density_g_cc=RHO_FOAM,
        option_combined=dict(after_4_perimeter_g=max(0.0, deficit_g - save_4p_g),
                             foam_after_4_perimeter_cc=max(0.0, deficit_g - save_4p_g) / (1 - RHO_FOAM)),
        option_note="4-perimeter (1.6 mm) shells first, then a small foam float on the tether for the rest "
                    "(option_combined); ballast only after that, because the hull is already near neutral",
        trim=dict(x_cg_minus_x_cb_mm=1000 * dx_cb, pitch_moment_Nm=W * dx_cb,
                  tether_pitch_stiffness_1p5_m_s_Nm_rad=k_pitch,
                  trim_angle_1p5_m_s_deg=math.degrees(W * dx_cb / k_pitch),
                  slack_water=("near-neutral, CG %.1f mm %s of CB: in still water it settles %s; in "
                               "current the rotor thrust aligns it with the flow (angle above)"
                               % (abs(1000 * dx_cb), "ahead" if dx_cb < 0 else "aft",
                                  "nose-down" if dx_cb < 0 else "tail-down"))),
        roll_pendulum=dict(BG_no_ballast_mm=v2_noball["BG_mm"],
                           BG_with_20g_lead_mm=v2_test["BG_mm"],
                           righting_max_with_20g_lead_Nm=v2_test["displacement_N"] * v2_test["BG_mm"] / 1000,
                           note="compare with the generator torque in roll.torque: the pendulum is "
                                "two orders of magnitude too weak, the fins must carry the torque"),
    )
    return dict(
        model=dict(shells="6 perimeters = solid 2.4 mm wall", solid_parts=f"{SKIN} mm skin + "
                   f"{int(INFILL * 100)} % gyroid", rotor="blades from layout_rotor.json blade_mass (solid, "
                   "incl. knuckle and heel); hub + spinner = hull envelope from the hub front to the tip minus "
                   "what floods (hub_wet: 3 blade slots %.1f cc, screw-access bore %.1f cc, spinner gap %.1f cc)"
                   % (hw_["slots_cc"], hw_["access_bore_cc"], hw_["gap_film_cc"]),
                   sealed="displacement counts the sealed static hull (minus the wet pockets), the fin rail "
                   "pads outside the hull (minus their wet grooves), and the outer volume of the rotor parts "
                   "and fins (closed perimeters trap the infill air)"),
        rotor_wet=hw_,
        parts=[dict(name=p["name"], group=p["group"], V_print_cc=p["V_print_mm3"] / 1000,
                    m_PLA_g=p["V_print_mm3"] * RHO_PLA / 1000, m_PETG_g=p["V_print_mm3"] * RHO_PETG / 1000,
                    x_mm=p["x_mm"]) for p in parts] + [dict(name="generator replica", group="int",
                    V_print_cc=gen_replica / 1000, m_PLA_g=gen_replica * RHO_PLA / 1000,
                    m_PETG_g=gen_replica * RHO_PETG / 1000, x_mm=xg)],
        hardware=[dict(name=n, mass_g=g, x_mm=x, basis=b) for n, g, x, b in hwl],
        static_hull_outer_cc=V_static / 1000, pockets_cc=V_pockets / 1000,
        display_PLA=_strip(disp_pla), v2_PETG_no_ballast=_strip(v2_noball), v2=_strip(v2),
        v2_ballast=dict(lead_g=ballast, lead_cc=lead_cc, x_mm=xb_, target_net_buoyancy_N=target_N,
                        where="bottom sector of the aft shell (theta 180) between the ribs at 120 "
                              "and 240 deg, shot in epoxy in a printed tray; trim by sliding it along x"),
        v2_needs=needs,
    )


def _strip(d):
    return {k: v for k, v in d.items() if k != "rows"} | dict(
        breakdown=[dict(name=n, g=g, x_mm=x) for n, g, x in d["rows"]])


# ==========================================================================
# Assembly and part split
# ==========================================================================


def assembly(J, bay, DT, nose_d, fins, hw, ri_rot: dict | None = None) -> list:
    st = DT["stations_mm"]
    g = bay["pick"]
    o = J["oring"]
    b = J["bayonet"]
    steps = [
        dict(step=1, sub="aft shell", action="Press the 688 down the tube from the open equator end, "
             "through the 608 seat (22) and the 18 mm bore, onto the seal shoulder (ID %.1f). Use a "
             "draw bolt (M6 rod, washers) through the tail end, not a hammer." % SEAL_SHOULDER_ID,
             check="688 OD 16 < 18 mm bore: passes. Seat depth %.0f mm below the tube mouth: draw "
                   "bolt needed. OK" % (st["bearing_rear_688"][0] - DT["tube"]["x_from"]), flag=False),
        dict(step=2, sub="aft shell", action="v2: press the lip seal into the 16 mm bore from the tail "
             "end face, main lip facing aft, onto the same shoulder.",
             check="seal bore open aft, recessed %.1f mm. OK" % SEAL_RECESS, flag=False),
        dict(step=3, sub="aft shell", action="Press the 608 into the front seat onto its 18 mm shoulder "
             "(same draw bolt).",
             check="both shoulders face forward, so no retaining rings; the tube prints tail-down and "
                   "both shoulders print as floors. OK", flag=False),
        dict(step=4, sub="aft shell", action="Slide the shaft in from the front through 608, 688 and "
             "the seal (chamfered end; a thin sleeve over the flat while it passes the seal lip).",
             check="the flat and the ring groove pass the seal lip: needs the chamfer and a sleeve. "
                   "FLAG for the prototype", flag=True),
        dict(step=5, sub="aft shell", action="Fit claw B on the forward shaft end against the 608 "
             "inner race; retaining ring in its groove; tighten the 2 grub screws on the flat (1.5 mm "
             "hex key).",
             check="claw B sits x %.1f-%.1f: its front %.1f mm stand ahead of the equator face, inside the "
                   "spigot (tip x %.1f); the spigot bore D%.1f gives radial access for a short 1.5 mm key. OK"
                   % (st["claw_B_thrust_collar"][0], st["claw_B_thrust_collar"][1],
                      max(0.0, J["x_eq"] - st["claw_B_thrust_collar"][0]), J["x_tip"], 2 * J["r_spigot_in"]),
             flag=False),
        dict(step=6, sub="rotor", action="Rotor (layout_rotor): blades into the hub clevis slots, 3 x 22 "
             "pins pressed into the ears through their tangential channels; hub onto the D-flat with a %.1f "
             "mm feeler in the hub gap, M3 grub tightened; spinner slid on from aft (keys on the slots, "
             "covers the pin channels), M3 screw into the tapped shaft end." % P.HUB_GAP,
             check="hub front x %.1f, shaft end x %.1f; rotor thrust: hub -> spinner -> M3 screw -> shaft "
                   "-> claw B -> 608. OK" % ((ri_rot or {}).get("hub_front_x", float("nan")),
                                              DT["shaft"]["x_aft"]), flag=False),
        dict(step=7, sub="aft shell", action="v2: bond the ballast tray (if any) in the bottom sector "
             "between the ribs.", check="sector free of ribs. OK", flag=False),
        dict(step=8, sub="aft shell", action="Grease the %s and roll it into the spigot groove; glue the "
             "2 joint magnets." % o["designation"],
             check="stretch %.1f %% nominal (%.1f-%.1f %% worst) in the groove; on the way the ring "
                   "rolls over the lugs and the forward land (D%.2f), a momentary %.0f %% stretch from its "
                   "%.2f ID (elastomer elongation limit not sourced; a smooth 45 deg lead-in on the lugs "
                   "helps). OK" % (o["stretch_nom_pct"], *o["stretch_range_pct"], o["spigot_seal_od_mm"],
                                   100 * (o["spigot_seal_od_mm"] / o["id_free_mm"] - 1), o["id_free_mm"]),
             flag=False),
        dict(step=9, sub="front shell", action="Screw the M12 gland into the nose boss (seat x %.1f, bottom) "
             "with its washer; pass the cable in and out through the open equator end, leaving ~120 mm "
             "(ESTIMATE) inside for the slack loop." % nose_d["gland"]["seat_x_mm"],
             check="the whole boss ends %.2f mm ahead of the carrier plate and %.2f mm ahead of its ring, so "
                   "the carrier never meets the gland; wrench SW %g on the outside. OK"
                   % (nose_d["gland"]["boss_to_carrier_plate_mm"], nose_d["gland"]["boss_to_carrier_ring_mm"],
                      hw["gland"]["SW"]),
             flag=not nose_d["gland"]["clear_ok"]),
        dict(step=10, sub="front shell", action="Bench: bolt the motor base to the carrier (%d x %s on "
             "PCD %.2f, heads in counterbores on the carrier front face), fit claw A, lay the cable in the "
             "%.1f mm rim notch at theta 180 and solder it to the motor leads." % (
                 hw["motor"]["n_screw"], hw["motor"]["screw"], hw["motor"]["pcd"],
                 nose_d["gland"]["cable_notch"]["width_mm"]),
             check="all screws reachable on the bench; notch %.1f mm from the nearest counterbore. OK"
                   % nose_d["gland"]["cable_notch"]["clear_to_screw_counterbores_mm"], flag=False),
        dict(step=11, sub="front shell", action="Push carrier + motor into the front shell through the "
             "socket (the cable slack folds into the nose cavity), tabs through the ring notches, twist 30 "
             "deg until it clicks.",
             check="motor D%.1f passes the socket lip D%.1f; the cable lies in the notch inside the carrier ring "
                   "and turns with the carrier; the generator reaction torque drives the tabs into their hard "
                   "stop. OK" % (hw["motor"]["D"], 2 * J["r_lip"]), flag=False),
        dict(step=12, sub="close", action="Align the entry marks, push the shells together (the claws "
             "find each other on the V lead-in; turn the rotor by hand if the dogs sit tip to tip), "
             "twist %.0f deg to the stop; the magnets click." % b["engage_deg"],
             check="the spigot (bore D%.1f) slides over the aft %.1f mm of the bay envelope (D%.0f): %.1f mm "
                   "radial clearance. Claw engagement is blind: self-finding with the V lead-in and a "
                   "free rotor. FLAG for the prototype"
                   % (2 * J["r_spigot_in"], max(0.0, g["x_aft_mm"] - J["x_tip"]), g["D_mm"],
                      J["r_spigot_in"] - g["D_mm"] / 2), flag=True),
        dict(step=13, sub="fins", action="Slide each fin aft onto its rail to the stop, set the grub "
             "screw. Remove for transport. To stow the blades, turn the rotor to its index (blades over "
             "their pockets) first, then fold.", check="at the index the folding blades clear the fins "
             "by >= %.1f mm; at other azimuths they can touch. FLAG for the user guide" % FIN_FOLD_CLEAR,
             flag=True),
        dict(step=14, sub="service", action="Reverse order. The 608 comes out by pushing the shaft "
             "from the tail (claw B drives it), the 688 with a drift from the tail after the seal.",
             check="688 is behind the seal: remove the seal first (it is destroyed; spare). OK", flag=False),
    ]
    return steps


def part_split(J, bay, fins, nose_d, ri) -> list:
    eye_fwd = -nose_d["eye"]["hole_centre_x_mm"] + EYE_HOLE_D / 2 + EYE_W
    pad = fins["root"]["pad_max_mm"]
    bed = P.BED[0]
    parts = [
        dict(part="front shell", material="PETG (v2) / PLA (display)",
             bbox_mm=[P.HULL_D + 2 * (pad + 2), P.HULL_D, J["x_eq"] + eye_fwd],
             orientation="equator face on the bed, nose up",
             supports="none: outer dome self-supporting; forward of x %.1f the tip is an infilled "
                      "plug whose aft face bridges %.1f mm; socket ring, lip underside and carrier-slot "
                      "ceilings chamfered 45 deg; eye hole teardrop; rail grooves run near-vertical; the fin "
                      "pads run down to the equator face, so their aft ends stand on the bed (no ledge); the "
                      "gland seat boss leans %.0f deg from vertical; pocket floors rise <= %.0f deg"
                      % (nose_d["nose_plug"]["x_45_mm"], nose_d["nose_plug"]["bridge_d_mm"],
                         180.0 - nose_d["gland"]["axis_angle_from_x_deg"],
                         ri.get("floor_overhang_max_deg") or float("nan"))),
        dict(part="aft shell (with bearing tube, ribs, spigot)", material="PETG / PLA",
             bbox_mm=[P.HULL_D, P.HULL_D, ri["static_tail_x"] - J["x_eq"] + J["L_spigot"]],
             orientation="tail end face on the bed, spigot up",
             supports="none: outer flank leans <= %.0f deg from vertical; bearing shoulders print "
                      "as floors; lug undersides and the O-ring groove forward flank at 45 deg; seal "
                      "shoulder is a %.2f mm ring overhang; pocket floors (incl. the root-notch ramp) rise "
                      "<= %.0f deg from vertical going up (rotor half printable_floor: the floor is cut deeper "
                      "where the blade would allow a steeper ramp); the stowed-detent boss on the notch floor is "
                      "a 4 mm wide, %.1f mm high bump (bridges)"
                      % (math.degrees(math.atan(abs(float(slope(ri["static_tail_x"]))))),
                         (P.LIP_SEAL["D"] - SEAL_SHOULDER_ID) / 2, ri.get("floor_overhang_max_deg") or float("nan"),
                         (ri.get("rotor") or {}).get("stow_boss_h", float("nan")))),
        dict(part="motor carrier (bulkhead)", material="PETG",
             bbox_mm=[2 * bay["pick"]["r_carrier_mm"] + 4] * 2 + [CARRIER_T],
             orientation="flat, front face down", supports="none (counterbores open downward are "
                                                            "holes, not overhangs)"),
        dict(part="claw A / claw B", material="PETG", bbox_mm=[COUPLER_OD, COUPLER_OD, CLAW_B_L],
             orientation="dog faces up", supports="none (45 deg V lead-in)"),
        dict(part="generator replica (display)", material="PLA", bbox_mm=[bay["pick"]["D_mm"]] * 2 +
             [bay["pick"]["L_mm"]], orientation="base down", supports="none"),
        dict(part="fin x2", material="PETG", bbox_mm=[round(fins["span_mm"] + RAIL_H + 2, 1),
             round(fins["chord_mm"] / math.cos(math.radians(fins["root"]["incl_deg"])) + 2, 1), 8.0],
             orientation="flat side of the Clark Y on the bed (the yawed root tongue lies in the fin "
                         "plane)", supports="none"),
        dict(part="ballast tray (v2)", material="PETG", bbox_mm=[30, 20, 40],
             orientation="open side up", supports="none"),
    ]
    for p in parts:
        p["fits_bed"] = max(p["bbox_mm"]) <= bed
    return parts


# ==========================================================================
# Figures
# ==========================================================================


def _style(plt):
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 9, "axes.titleweight": "bold", "savefig.dpi": 200,
                         "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e",
                         "font.family": "DejaVu Sans", "hatch.linewidth": 0.4})


# palette: neutrals for structure, reference categorical slots for the functional parts
COL_INK, COL_INK2 = "#0b0b0b", "#52514e"
COL_HULL, COL_WALL, COL_WALL2 = "#5B6670", "#C4C9CE", "#DADDE0"
COL_PK = "#eb6834"          # pockets (slot 2)
COL_GEN = "#2a78d6"         # generator (slot 1)
COL_CPL = "#1baf7a"         # coupler (slot 3)
COL_CAR = "#4a3aa7"         # carrier / adapter (slot 7)
COL_JNT, COL_JNT2 = "#B9A886", "#8E7B58"
COL_STEEL, COL_BRG, COL_RUB = "#8A9299", "#2F3337", "#111111"


def figures(res: dict) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Polygon, Rectangle, Wedge

    J, bay, DT, fins = res["joint"], res["generator_bay"], res["drivetrain"], res["fins"]
    nose_d, mb, pk, ri = res["nose"], res["mass"], res["_pk"], res["inputs"]["rotor"]
    hw = res["inputs"]["hardware"]
    _style(plt)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    x_end = ri["static_tail_x"]
    st = DT["stations_mm"]
    g = bay["pick"]
    t = DT["tube"]
    x = np.linspace(0, P.HULL_L, 1600)
    xs = x[x <= x_end]

    def sock_r(xr, lug=False):
        """Inner radius of the front-shell socket ring. lug=True: section through a lug, where
        the twist channel (with its 45 deg ceiling) is cut ahead of the lip."""
        base = np.where(xr >= J["x_seal0"], J["r_bore"],
                        np.where(xr >= J["x_lip0"], J["r_lip"],
                                 np.where(xr >= J["x_ring0"], J["r_lip"],
                                          np.minimum(r_in(xr), J["r_lip"] + (J["x_ring0"] - xr)))))
        if lug:
            ch = np.where((xr >= J["x_ch0"]) & (xr < J["x_lip0"]), J["r_channel"],
                          np.where((xr >= J["x_ring0"]) & (xr < J["x_ch0"]),
                                   J["r_lip"] + (xr - J["x_ring0"]), 0.0))
            base = np.maximum(base, ch)
        return base

    def hull(ax, sides=(1, -1), pocket_top=True):
        for s_ in sides:
            ax.plot(x, s_ * r_out(x), color=COL_HULL, lw=1.3, zorder=6)
            rin_s = np.array(r_in(xs))
            ax.fill_between(xs, s_ * rin_s, s_ * r_out(xs), color=COL_WALL, lw=0, zorder=2)
        if pocket_top and 1 in sides:
            xp = np.linspace(pk.x0, min(pk.x1, x_end), 400)
            d = pk.depth(xp)
            ax.fill_between(xp, r_out(xp) - d, r_out(xp) + 0.05, color="white", lw=0, zorder=3)
            ax.fill_between(xp, pk.floor_inner_r(xp), r_out(xp) - d, color=COL_WALL, lw=0, zorder=3)
            ax.plot(xp, r_out(xp) - d, color=COL_PK, lw=1.1, zorder=6, label="blade pocket floor (deepest over azimuth)")
            ax.plot(xp, r_out(xp), color=COL_PK, lw=0.6, ls=(0, (2, 2)), zorder=6)
        xh = x[x >= x_end + P.HUB_GAP]
        for k, s_ in enumerate(sides):
            ax.plot(xh, s_ * r_out(xh), color="#9AA3AB", lw=1.0, ls=(0, (4, 2)), zorder=6,
                    label="hub + spinner outline (rotor side)" if k == 0 else None)
        ax.axhline(0, color="#9a9a9a", lw=0.5, ls=(0, (8, 3, 1, 3)), zorder=1)

    def rect(ax, x0, x1, r0, r1, sides=(1, -1), **kw):
        lab = kw.pop("label", None)
        for k, s_ in enumerate(sides):
            lo, hi = sorted((s_ * r0, s_ * r1))
            ax.add_patch(Rectangle((x0, lo), x1 - x0, hi - lo, label=lab if k == 0 else None, **kw))

    def joint_parts(ax, sides=(1, -1), lug_side=-1):
        xr = np.linspace(J["x_ring0"] - 8, J["x_eq"], 300)
        for s_ in sides:
            # socket ring (front shell): from the socket bore out to the hull inner wall, or to
            # the pocket floor on the pocket side (theta 0, top)
            r_s = np.array(sock_r(xr, lug=(s_ == lug_side)))
            top = np.array(r_in(xr))
            if s_ > 0:
                d_ = pk.depth(xr)
                top = np.where(d_ > 0, np.minimum(top, np.array(r_out(xr)) - d_), top)
            ax.fill_between(xr, s_ * r_s, s_ * top, color=COL_JNT, lw=0, zorder=4,
                            label="socket ring (front shell)" if s_ == sides[0] else None)
        # spigot (aft shell) and its root ring (45 deg underside, prints tail-down)
        xq = np.linspace(J["x_eq"], J["x_eq"] + J["root_ring_L"], 60)
        for s_ in sides:
            top = np.array(r_in(xq))
            if s_ > 0:
                top = np.minimum(top, pk.floor_inner_r(xq))
            inner = np.minimum(J["r_spigot_in"] + (xq - J["x_eq"]), top)
            ax.fill_between(xq, s_ * inner, s_ * top, color=COL_JNT2, lw=0, zorder=4.5)
        rect(ax, J["x_seal0"], J["x_eq"], J["r_spigot_in"], J["r_spigot_seal"], sides, fc=COL_JNT2, lw=0,
             zorder=5, label="spigot + root ring (aft shell)")
        rect(ax, J["x_tip"], J["x_seal0"], J["r_spigot_in"], J["r_spigot_lug"], sides, fc=COL_JNT2, lw=0, zorder=5)
        gr = J["groove"]
        for s_ in sides:
            # groove cut (bottom width + 45 deg forward flank) and the O-ring section
            poly = [(gr["x_aft"], s_ * J["r_spigot_seal"]), (gr["x_aft"], s_ * J["r_groove"]),
                    (gr["x_fwd_bottom"], s_ * J["r_groove"]), (gr["x_fwd_top"], s_ * J["r_spigot_seal"])]
            ax.add_patch(Polygon(poly, closed=True, fc="white", lw=0, zorder=5.5))
            cx = 0.5 * (gr["x_aft"] + gr["x_fwd_bottom"])
            ax.add_patch(Circle((cx, s_ * (J["r_groove"] + 0.5 * (J["r_bore"] - J["r_groove"]))),
                                0.5 * (J["r_bore"] - J["r_groove"]) * 1.05, fc=COL_RUB, lw=0, zorder=6,
                                label=J["oring"]["designation"].split(" / ")[-1] + " O-ring"
                                if s_ == sides[0] else None))
        if lug_side in sides:
            ax.add_patch(Polygon([(J["x_tip"], lug_side * J["r_spigot_lug"]),
                                  (J["x_tip"], lug_side * J["r_lug"]),
                                  (J["x_tip"] + LUG_T - LUG_H, lug_side * J["r_lug"]),
                                  (J["x_tip"] + LUG_T, lug_side * J["r_spigot_lug"])],
                                 closed=True, fc="#5E4E33", lw=0, zorder=6, label="bayonet lug (locked)"))

    def internals(ax, sides=(1, -1)):
        # nose plug + bridge
        npl = nose_d["nose_plug"]
        xq = np.linspace(0.3, npl["x_45_mm"], 80)
        for s_ in sides:
            ax.fill_between(xq, 0, s_ * np.array(r_in(xq)), fc="none", ec=COL_HULL, hatch="////", lw=0,
                            zorder=2.5, label="nose plug (infill, bridged)" if s_ == sides[0] else None)
            ax.plot([npl["x_45_mm"]] * 2, [0, s_ * npl["r_inner_at_45_mm"]], color=COL_HULL, lw=0.8, zorder=6)
        # tow eye (ring in the XZ plane: the side section shows it face-on)
        e = nose_d["eye"]
        xc = e["hole_centre_x_mm"]
        Ro = EYE_HOLE_D / 2 + EYE_W
        ax.add_patch(Wedge((xc, 0.0), Ro, 0, 360, width=EYE_W, fc=COL_WALL, ec=COL_HULL, lw=0.7, zorder=6,
                           label="tow eye (ring in the XZ plane, hole axis Y)"))
        # gland at theta 180 (bottom): seat D18 + thread boss D16 inside, gland body outside
        if -1 in sides:
            gl = nose_d["gland"]
            a_ = math.radians(gl["axis_angle_from_x_deg"])
            n = np.array([math.cos(a_), math.sin(a_)])          # outward axis in (x, r)
            tdir = np.array([math.sin(a_), -math.cos(a_)])
            base = np.array([gl["seat_x_mm"], gl["seat_surface_r_mm"]])
            Lt = hw["gland"]["thread_L"]
            hgt = hw["gland"]["height"] - Lt
            rb = hw["gland"]["SW"] / math.sqrt(3.0)

            def quad(w0, w1, d0, d1):
                return [(p_[0], -p_[1]) for p_ in (base + tdir * w0 + n * d0, base - tdir * w0 + n * d0,
                                                    base - tdir * w1 + n * d1, base + tdir * w1 + n * d1)]
            ax.add_patch(Polygon(quad(GLAND_BOSS_D / 2, GLAND_BOSS_D / 2, -GLAND_SEAT_T, -Lt), closed=True,
                                 fc=COL_WALL, lw=0, zorder=4.6))
            ax.add_patch(Polygon(quad(GLAND_SEAT_D / 2, GLAND_SEAT_D / 2, 0.0, -GLAND_SEAT_T), closed=True,
                                 fc=COL_WALL, lw=0, zorder=4.6))
            ax.add_patch(Polygon(quad(6.0, 6.0, 0.0, -Lt), closed=True, fc="#9A9A9A", lw=0, zorder=5))
            ax.add_patch(Polygon(quad(rb, rb * 0.8, 0.0, hgt), closed=True, fc="#444", ec="#222", lw=0.5,
                                 zorder=6, label="M12 cable gland (theta 180, x %.1f)" % gl["seat_x_mm"]))
        # carrier ring on the shell + carrier plate + adapter spacer
        xcar0 = g["x_fwd_mm"] - CARRIER_T
        for s_ in sides:
            rr = float(r_in(xcar0 + CARRIER_T / 2))
            ax.add_patch(Rectangle((xcar0 - 1.2, s_ * rr if s_ < 0 else rr - CARRIER_RING), CARRIER_T + 2.4,
                                   CARRIER_RING, fc=COL_WALL, lw=0, zorder=4.5))
        rect(ax, xcar0, g["x_fwd_mm"], 0, g["r_carrier_mm"], sides, fc=COL_CAR, lw=0, zorder=5,
             label="motor carrier (bulkhead, twist-lock)")
        x_body0 = g["motor_body_x_mm"][0]
        if x_body0 - g["x_fwd_mm"] > 0.5:
            r_o = min(hw["motor"]["D"], g["D_mm"]) / 2 - 1.0
            rect(ax, g["x_fwd_mm"], x_body0, r_o - 2.4, r_o, sides, fc="#9C93D6", lw=0, zorder=5,
                 label="adapter / spacer ring")
            rect(ax, x_body0 - 3.0, x_body0, 4.5, r_o, sides, fc="#9C93D6", lw=0, zorder=5)
        # generator: reserved envelope (dashed) and the recommended motor body
        rect(ax, g["x_fwd_mm"], g["x_aft_mm"], 0, g["D_mm"] / 2, sides, fc="none", ec=COL_GEN, lw=0.9,
             ls=(0, (3, 2)), zorder=7, label=f"bay reserve D{g['D_mm']:g} x {g['L_mm']:g}")
        rect(ax, x_body0, g["motor_body_x_mm"][1], 4.0, hw["motor"]["D"] / 2, sides, fc=COL_GEN, alpha=0.9,
             lw=0, zorder=5.5, label=f"{hw['motor']['name']} D{hw['motor']['D']:g} x {hw['motor']['L']:g}")
        # coupler
        rect(ax, *st["claw_A"], 4.0, COUPLER_OD / 2, sides, fc=COL_CPL, lw=0, zorder=6,
             label="claw coupling A / B (B = thrust collar)")
        for seg in DT["coupler"]["claw_B_segments"]:
            rect(ax, *seg["x_mm"], 4.0, seg["od_mm"] / 2, sides, fc=COL_CPL, lw=0, zorder=6)
        rect(ax, st["claw_B_thrust_collar"][0] + 1.0, st["claw_B_thrust_collar"][1] - 1.0, 4.0, 8.0, sides,
             fc="#0E7A55", lw=0, zorder=6.2)
        # bearing tube, ribs, bearings, seal
        f0, f1 = st["bearing_front_608"]
        r0_, r1_ = st["bearing_rear_688"]
        xt = np.array([t["x_from"], f1 + 2, r0_ - 2, x_end])
        rt = np.array([t["r_out_front_mm"], t["r_out_front_mm"], t["r_out_rear_mm"] + 1, t["r_out_rear_mm"]])
        rb = np.array([t["bore_front_mm"] / 2, t["shoulder_front_id_mm"] / 2, t["shoulder_front_id_mm"] / 2,
                       t["seal_bore_mm"] / 2])
        for k, s_ in enumerate(sides):
            ax.fill_between(xt, s_ * rb, s_ * rt, color=COL_WALL2, lw=0, zorder=4,
                            label="bearing tube (aft shell)" if k == 0 else None)
        xr2 = np.linspace(t["x_from"] + 1, x_end - 2, 120)
        if 1 in sides:
            ax.fill_between(xr2, np.interp(xr2, xt, rt), pk.floor_inner_r(xr2), color="#E6E8EA", lw=0,
                            zorder=3.5, label="rib under each pocket (theta 0/120/240)")
            ax.plot(xr2, np.interp(xr2, xt, rt), color="#B5B9BD", lw=0.4, zorder=3.6)
        xb0 = st["bearing_rear_688"][1]
        for s_ in sides:
            xq_ = np.linspace(xb0, x_end, 30)
            top_ = np.array(r_in(xq_))
            if s_ > 0:
                top_ = np.minimum(top_, np.array(r_out(xq_)) - pk.depth(xq_))
            ax.fill_between(xq_, s_ * t["seal_bore_mm"] / 2, s_ * top_, color=COL_WALL2, lw=0, zorder=4)
        rect(ax, f0, f1, 4.0, t["bore_front_mm"] / 2, sides, fc=COL_BRG, lw=0, zorder=7,
             label="608-2RS (front) / 688-2RS (rear)")
        rect(ax, r0_, r1_, 4.0, t["bore_rear_mm"] / 2, sides, fc=COL_BRG, lw=0, zorder=7)
        rect(ax, *st["lip_seal"], 4.0, t["seal_bore_mm"] / 2, sides, fc="#6B5A47", lw=0, zorder=7,
             label="lip seal TC 8x16x5 (v2)")
        # shaft
        ax.add_patch(Rectangle((st["shaft"][0], -4), st["shaft"][1] - st["shaft"][0], 8, fc=COL_STEEL,
                               lw=0, zorder=8, label=f"shaft 8 x {st['shaft'][1] - st['shaft'][0]:.0f}"))
        # ballast tray position (bottom sector)
        if -1 in sides:
            b = mb["v2_ballast"]
            xb_ = b["x_mm"]
            rb_ = float(r_in(xb_))
            ax.add_patch(Rectangle((xb_ - 15, -rb_ + 0.3), 30, 6, fc="none", ec=COL_INK2, lw=0.7, ls=(0, (2, 1.5)),
                                   zorder=6, label="v2 ballast tray (bottom sector)"))

    rj = _load(ROTOR_JSON) or {}

    def rotor_parts(ax):
        """Rotor half (data/layout_rotor.json): hub core with blade 0's slot (top half), spinner
        shell, pin, knuckle and deployed heel of blade 0, the stowed blade 0 envelope and the
        rotor's keep-out line (floor - MIN_WALL) under the pocket."""
        if not rj or "hub" not in rj:
            return
        ch, hub = rj["chosen"], rj["hub"]
        xh = np.linspace(ch["hub_front_x_mm"], hub["hub_aft_face_x_mm"], 200)
        r_core = np.array(r_out(xh)) - ch["spinner_wall_mm"] - ch["spinner_gap_mm"]
        r_bore = hub["bore_d_mm"] / 2
        sl = rj["hub_slot"]
        floor = np.interp(xh, sl["x_mm"], sl["slot_floor_re_mm"], left=np.inf, right=np.inf)
        top = np.minimum(r_core, floor)
        ax.fill_between(xh, r_bore, top, color="#D5DCE6", lw=0, zorder=4, label="hub core (rotor half)")
        ax.fill_between(xh, -r_bore, -r_core, color="#D5DCE6", lw=0, zorder=4)
        xs_ = np.linspace(ch["hub_front_x_mm"], P.HULL_L, 300)
        for s_ in (1, -1):
            ax.fill_between(xs_, s_ * (np.array(r_out(xs_)) - ch["spinner_wall_mm"]), s_ * np.array(r_out(xs_)),
                            color="#AEB9C7", lw=0, zorder=4.2, label="spinner shell" if s_ > 0 else None)
        hx, hr = ch["hinge_x_mm"], ch["hinge_r_mm"]
        st_ = hub["stop"]
        ax.add_patch(Rectangle((st_["stop_block_x_mm"][0], st_["stop_face_r_mm"][0]),
                               st_["stop_block_x_mm"][1] - st_["stop_block_x_mm"][0],
                               st_["stop_face_r_mm"][1] - st_["stop_face_r_mm"][0], fc="#8FA0B5", lw=0, zorder=4.6,
                               label="stop block (hub front ring)"))
        ht, hrch = ch["heel_t_mm"], ch["heel_reach_mm"]
        ax.add_patch(Rectangle((hx - ht / 2, hr - hrch), ht, hrch, fc=COL_PK, alpha=0.55, lw=0, zorder=5,
                               label="blade 0 knuckle + heel (deployed)"))
        ax.add_patch(Circle((hx, hr), ch["knuckle_r_mm"], fc=COL_PK, alpha=0.55, lw=0, zorder=5))
        ax.add_patch(Circle((hx, hr), ch["pin_d_mm"] / 2, fc=COL_STEEL, ec=COL_INK, lw=0.4, zorder=8,
                            label="hinge pin %g x %g (tangential)" % (ch["pin_d_mm"], ch["pin_l_mm"])))
        ss = rj["stowed_sections"]
        xs_s = np.array(ss["x_stowed_mm"])
        hr_s = np.array(ss["hull_r_mm"])
        o_ = np.argsort(xs_s)
        ax.fill_between(xs_s[o_], (hr_s - np.array(ss["depth_below_hull_mm"]))[o_],
                        (hr_s + np.array(ss["top_over_hull_mm"]))[o_], color=COL_PK, alpha=0.30, lw=0, zorder=5,
                        label="stowed blade 0, aero span (deepest point per station)")
        xp_ = np.array(rj["pocket_depth"]["x_mm"])
        ax.plot(xp_, rj["pocket_depth"]["keepout_r_mm"], color=COL_PK, lw=0.7, ls=(0, (1, 1.5)), zorder=6,
                label="rotor keep-out (pocket floor - MIN_WALL)")

    # ---------------------------------------------------------------- figure 1: side section
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    hull(ax)
    joint_parts(ax)
    internals(ax)
    rotor_parts(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-16, 170)
    ax.set_ylim(-47, 47)
    ax.set_xlabel("x (mm), flow ->")
    ax.set_ylabel("z (mm)")
    ax.set_title("HydroEgg merged layout: side section (top half: blade 0 plane, theta 0, pocket floor = "
                 "deepest point over azimuth; bottom half: theta 180, between pockets)", loc="left")
    for xv, lab in ((J["x_eq"], "equator %.1f" % J["x_eq"]), (x_end, "static tail end %.0f" % x_end)):
        ax.axvline(xv, color="#9a9a9a", lw=0.6, ls=":", zorder=1)
        ax.annotate(lab, (xv, 44.5), ha="center", fontsize=6.5, color=COL_INK2)
    ax.annotate("rotor half: hub, spinner,\nblade 0 root (layout_rotor.json)", (x_end + 22, 27), fontsize=6.5,
                color=COL_INK2)
    ax.annotate("pockets: " + ("layout_rotor.json, x %.1f-%.0f, max %.1f mm deep" % (ri["x0"], ri["x1"], ri["depth"])
                               if ri["depth_x"] else "fallback 6 mm x 30 mm, x 45-118"),
                (40, 39.5), fontsize=6.5, color=COL_PK)
    ax.annotate("rotor face plane x %.2f" % g["x_aft_mm"], (g["x_aft_mm"], -44.5), ha="center", fontsize=6.5,
                color=COL_INK2)
    h, lab = ax.get_legend_handles_labels()
    ax.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=6.5, frameon=False)
    fig.tight_layout()
    p1 = FIG_DIR / "layout_section.png"
    fig.savefig(p1, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    out.append(p1)

    # ---------------------------------------------------------------- figure 2: joint detail + bay frontier
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw=dict(width_ratios=[1.3, 1]))
    hull(a1, sides=(-1,), pocket_top=False)
    joint_parts(a1, sides=(-1,), lug_side=-1)
    internals(a1, sides=(-1,))
    xp = np.linspace(J["x_ring0"] - 8, J["x_eq"] + 12, 200)
    a1.plot(xp, -(np.array(r_out(xp)) - pk.depth(xp)), color=COL_PK, lw=0.9, ls=(0, (3, 2)),
            label="pocket floor at theta 0 (other half-plane)")
    a1.set_aspect("equal")
    a1.set_xlim(J["x_ring0"] - 6, J["x_eq"] + 9)
    a1.set_ylim(-38.5, -21.5)
    a1.set_title("Equator joint, x = %.1f (theta 180 half-plane, lug drawn in its channel)" % J["x_eq"], loc="left")
    a1.set_xlabel("x (mm)")
    a1.set_ylabel("z (mm)")
    o = J["oring"]
    a1.text(0.99, 0.97, "seal bore D%.2f, spigot OD %.2f\nlip D%.2f, lug OD %.2f\nspigot bore D%.2f\n"
            "O-ring %s: stretch %.1f %%, squeeze %.0f-%.0f %% (worst)\nbayonet: 3 lugs %.0f deg, twist %.0f deg"
            % (2 * J["r_bore"], 2 * J["r_spigot_seal"], 2 * J["r_lip"], 2 * J["r_lug"], 2 * J["r_spigot_in"],
               o["designation"].split(" / ")[-1], o["stretch_nom_pct"], *o["squeeze_worst_pct"],
               LUG_W_DEG, J["bayonet"]["engage_deg"]),
            transform=a1.transAxes, ha="right", va="top", fontsize=6.3, color=COL_INK,
            bbox=dict(fc="white", ec="#cccccc", lw=0.5, pad=3))
    h, lab = a1.get_legend_handles_labels()
    keep = [i for i, l_ in enumerate(lab) if l_.startswith(("socket", "spigot", "AS568", "bayonet", "pocket",
                                                             "bay reserve", hw["motor"]["name"]))]
    a1.legend([h[i] for i in keep], [lab[i] for i in keep], fontsize=6, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=3)
    fr = bay["frontier"]
    Ls = np.array([f_["L_mm"] for f_ in fr])
    Ds = np.array([f_["D_max_mm"] for f_ in fr])
    m = Ds > 0
    a2.add_patch(Rectangle((20, 40), hw["envelope"]["L_max"] - 20, hw["envelope"]["D_max"] - 40, fc="#EEEEEC",
                           ec="none", zorder=0, label="candidate box D 40-56 x L 20-35"))
    a2.step(Ls[m], Ds[m], where="mid", color=COL_GEN, lw=2.0, label="max D for each L (bay frontier)")
    for c_ in bay["candidate_fits"]:
        mk = "o" if c_["fits_bay"] else "x"
        a2.plot([c_["L_mm"]], [c_["D_mm"]], mk, color=COL_INK if c_["fits_bay"] else "#e34948", ms=5, mew=1.2,
                zorder=5)
        short = c_["name"].split(" (")[0].replace("CubeMars ", "").replace("T-motor ", "").replace("iPower ", "")
        short = short.replace("Rctimer ", "")
        off = {"GL40 II KV82.5": (4, -9), "GBM4114HS 100T": (-8, -10), "GL40 KV70": (-14, -10),
               "GM4108H-120T": (4, 3), "GM5208-24": (4, -8)}.get(short, (4, 2))
        a2.annotate(short, (c_["L_mm"], c_["D_mm"]), xytext=off, textcoords="offset points", fontsize=5.8,
                    color=COL_INK2)
    a2.plot([], [], "o", color=COL_INK, ms=5, label="catalogue motor that fits the chosen bay")
    a2.plot([], [], "x", color="#e34948", ms=5, mew=1.2, label="does not fit the chosen bay")
    a2.add_patch(Rectangle((0, 0), g["L_mm"], g["D_mm"], fc="none", ec=COL_PK, lw=1.4, zorder=4,
                           label=f"chosen bay reserve D{g['D_mm']:g} x L{g['L_mm']:g}"))
    a2.set_xlim(17, 37)
    a2.set_ylim(38, 63)
    a2.annotate("D limited by the nose wall and the carrier ring at the forward end of the bay",
                (17.3, Ds[0] + 0.8), fontsize=6, color=COL_GEN)
    if m.any():
        a2.plot([Ls[m][-1] + 0.5] * 2, [38, Ds[m][-1]], color=COL_GEN, lw=1.0, ls=(0, (2, 2)))
        a2.annotate("no longer bay: span >= %g mm\n(rotor face <= x %.2f)" % (SPAN_MIN, bay["x_bell_max_mm"]),
                    (Ls[m][-1] + 0.8, 51.0), fontsize=6, color=COL_GEN)

    a2.set_xlabel("generator length L (mm)")
    a2.set_ylabel("generator diameter D (mm)")
    a2.set_title("Generator bay: largest cylinder vs length", loc="left")
    a2.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    fig.tight_layout()
    p2 = FIG_DIR / "layout_section_joint_bay.png"
    fig.savefig(p2, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    out.append(p2)

    # ---------------------------------------------------------------- figure 3: cross sections
    x_fin_cut = min(max(40.0, fins["x_le_mm"] + 2.0), fins["x_te_mm"] - 2.0)
    x_fin_cut = min(x_fin_cut, g["x_aft_mm"] - 1.0) if g["x_aft_mm"] - 1.0 > fins["x_le_mm"] else x_fin_cut
    cuts = [(x_fin_cut, "generator + fins"), (J["x_tip"] + 2.0, "bayonet lugs"),
            (sum(st["bearing_front_608"]) / 2, "front bearing 608"), (sum(st["bearing_rear_688"]) / 2,
                                                                       "rear bearing 688")]
    fig, axs = plt.subplots(1, 4, figsize=(11.5, 3.3))
    th = np.linspace(0, 2 * np.pi, 721)
    for ax, (xc_, lbl) in zip(axs, cuts):
        R = float(r_out(xc_))
        Ri = float(r_in(xc_))
        d = float(pk.depth(xc_))
        rr = np.full_like(th, R)
        rr_in = np.full_like(th, Ri)
        lo_, hi_ = (float(v) for v in pk.theta_span(xc_))
        for tp in pk.theta:
            dth = np.degrees(np.angle(np.exp(1j * (th - math.radians(tp)))))
            m_ = (dth >= lo_) & (dth <= hi_) if d > 0 else np.zeros_like(th, dtype=bool)
            rr[m_] = R - d
            rr_in[m_] = np.minimum(Ri, float(pk.floor_inner_r(xc_)))
        Y = rr * np.sin(th)
        Z = rr * np.cos(th)
        ax.fill(np.r_[Y, (rr_in * np.sin(th))[::-1]], np.r_[Z, (rr_in * np.cos(th))[::-1]], color=COL_WALL, lw=0)
        ax.plot(Y, Z, color=COL_HULL, lw=1.0)
        if d > 0:
            for tp in pk.theta:
                a_ = math.radians(tp)
                ax.annotate("", (0.9 * R * math.sin(a_), 0.9 * R * math.cos(a_)), fontsize=5)
        if "fins" in lbl:
            ax.add_patch(Circle((0, 0), hw["motor"]["D"] / 2, color=COL_GEN, alpha=0.9))
            ax.add_patch(Circle((0, 0), g["D_mm"] / 2, fill=False, ec=COL_GEN, ls=(0, (3, 2)), lw=0.9))
            rl = fins["root"]
            y_root = rl["y_le"] + (xc_ - rl["x_le"]) / (rl["x_te"] - rl["x_le"]) * (rl["y_te"] - rl["y_le"])
            tl = 0.117 * fins["chord_mm"]

            yw = fins["yaw"]
            off_c = (xc_ - yw["pivot_x_mm"]) * yw["offset_per_mm_x"]   # yawed chord line at this x

            def radial_bar(th_deg, r0, r1, half_w, off=0.0, **kw):
                a_ = math.radians(th_deg)
                er = np.array([math.sin(a_), math.cos(a_)])          # (Y, Z) of the radial direction
                et = np.array([math.cos(a_), -math.sin(a_)])         # + theta
                ax.add_patch(Polygon([r0 * er + (off - half_w) * et, r1 * er + (off - half_w) * et,
                                      r1 * er + (off + half_w) * et, r0 * er + (off + half_w) * et],
                                     closed=True, **kw))
            for th_f in fins["theta_deg"]:
                radial_bar(th_f, R - 0.5, y_root, yw["pad_half_width_at_constant_x_mm"], off_c, fc="#7a7a7a", lw=0)
                radial_bar(th_f, y_root, y_root + 60, tl / 2, off_c, fc="#2F3337", lw=0)
            ax.annotate("fins at theta %g / %g (span %.0f, clipped)" % (*fins["theta_deg"], fins["span_mm"]),
                        (0, -R - 6), ha="center", fontsize=6, color=COL_INK2)
            lim = 44
        elif "lugs" in lbl:
            ax.add_patch(Circle((0, 0), J["r_lip"], fill=False, ec=COL_JNT, lw=2.4))
            ax.add_patch(Circle((0, 0), J["r_spigot_in"], fill=False, ec=COL_JNT2, lw=1.0))
            ax.add_patch(Circle((0, 0), (J["r_spigot_lug"] + J["r_spigot_in"]) / 2, fill=False, ec=COL_JNT2,
                                lw=(J["r_spigot_lug"] - J["r_spigot_in"]) * 2.2))
            for tl_ in J["bayonet"]["lock_theta_deg"]:
                ax.add_patch(Wedge((0, 0), J["r_lug"], 90 - tl_ - LUG_W_DEG / 2, 90 - tl_ + LUG_W_DEG / 2,
                                   width=LUG_H, fc="#5E4E33", lw=0))
                ax.add_patch(Wedge((0, 0), J["r_channel"] + 0.6, 90 - tl_ - LUG_W_DEG / 2,
                                   90 - tl_ + LUG_W_DEG / 2 + J["bayonet"]["engage_deg"], width=0.6,
                                   fc=COL_INK2, lw=0))
            ax.add_patch(Circle((0, 0), hw["motor"]["D"] / 2, color=COL_GEN, alpha=0.9))
            ax.annotate("lugs locked at %s deg,\ntwist %.0f deg (arc above each)"
                        % ("/".join("%.0f" % v for v in J["bayonet"]["lock_theta_deg"]), J["bayonet"]["engage_deg"]),
                        (0, -R - 7), ha="center", fontsize=6, color=COL_INK2)
            lim = 44
        else:
            ro_t = t["r_out_front_mm"] if "608" in lbl else t["r_out_rear_mm"]
            bo = (t["bore_front_mm"] if "608" in lbl else t["bore_rear_mm"]) / 2
            for tp in pk.theta:
                a_ = math.radians(tp)
                r_f = float(pk.floor_inner_r(xc_))
                ax.plot([ro_t * math.sin(a_), r_f * math.sin(a_)], [ro_t * math.cos(a_), r_f * math.cos(a_)],
                        color="#AFB4B8", lw=2.2)
            ax.add_patch(Circle((0, 0), ro_t, color=COL_WALL2))
            ax.add_patch(Circle((0, 0), bo, color=COL_BRG))
            ax.add_patch(Circle((0, 0), 4, color=COL_STEEL))
            ax.annotate("tube + 3 ribs under the pockets", (0, -R - 6), ha="center", fontsize=6, color=COL_INK2)
            lim = 44
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(f"x = {xc_:.1f}: {lbl}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle("Cross sections seen from aft, looking upstream (+Z up, +Y right); pockets at theta %s deg; "
                 "fins and pads drawn with their yaw offset at the cut" %
                 "/".join(f"{v:g}" for v in pk.theta), fontsize=8)
    fig.tight_layout()
    p3 = FIG_DIR / "layout_section_cross.png"
    fig.savefig(p3, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    out.append(p3)

    # ---------------------------------------------------------------- figure 4: fins, plan view
    fig, (ax, a2) = plt.subplots(1, 2, figsize=(12.5, 6.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    xx = np.linspace(0, P.HULL_L, 600)
    ax.plot(xx, r_out(xx), color=COL_HULL, lw=1.3)
    ax.plot(xx, -np.array(r_out(xx)), color=COL_HULL, lw=1.3)
    rl = fins["root"]
    s_ = fins["span_mm"]
    pad_x1 = fins["rail"]["pad_x_mm"][1]
    k_rl = (rl["y_te"] - rl["y_le"]) / (rl["x_te"] - rl["x_le"])
    for sg in (1, -1):
        poly = [(rl["x_le"], sg * rl["y_le"]), (rl["x_te"], sg * rl["y_te"]),
                (rl["x_te"], sg * (rl["y_te"] + s_)), (rl["x_le"], sg * (rl["y_le"] + s_))]
        ax.add_patch(Polygon(poly, closed=True, fc="#D9DCE0", ec=COL_BRG, lw=1.0,
                             label="fin (Clark Y, %.0f x %.0f, incidence %.2f deg)" %
                                   (fins["chord_mm"], s_, fins["incidence_deg"]) if sg > 0 else None))
        xq = np.linspace(rl["x_le"], pad_x1, 60)
        ax.fill_between(xq, sg * np.array(r_out(xq)), sg * (rl["y_le"] + (xq - rl["x_le"]) * k_rl),
                        color="#8a8a8a", lw=0,
                        label="rail pad + dovetail (runs to the equator face)" if sg > 0 else None)
        ax.plot([rl["x_le"], rl["x_te"]], [sg * (rl["y_le"] + FIN_CP_SPAN * s_ + 0 * s_),
                                           sg * (rl["y_te"] + FIN_CP_SPAN * s_)], color=COL_INK2, lw=0.6,
                ls=(0, (2, 2)))
    Lb = fins["clearance"]["blade_sweep_len_mm"]
    for sg in (1, -1):
        hc = (ri["hinge_x"], sg * ri["hinge_r"])
        ang = np.linspace(math.pi / 2, math.pi * 1.0, 100) if sg > 0 else np.linspace(math.pi, 1.5 * math.pi, 100)
        ax.plot(hc[0] + Lb * np.cos(ang), hc[1] + Lb * np.sin(ang), color=COL_PK, lw=1.2,
                label="blade fold torus, any azimuth (reference: blades fold at the index)" if sg > 0 else None)
        ax.plot(hc[0] + (Lb + FIN_TIP_CLEAR) * np.cos(ang), hc[1] + (Lb + FIN_TIP_CLEAR) * np.sin(ang),
                color=COL_PK, lw=0.7, ls=(0, (3, 2)), label="+ %.0f mm clearance" % FIN_TIP_CLEAR if sg > 0 else None)
        ax.plot([hc[0]], [hc[1]], "o", color=COL_PK, ms=3)
    x_tip_dep = ri["hinge_x"] + float((rj.get("axis") or {}).get("aft_tip_mm", 0.0))
    ax.plot([x_tip_dep] * 2, [-ri["r_tip"], ri["r_tip"]], color=COL_PK, lw=0.8, ls=":",
            label="deployed blade tips (R %.0f, x %.1f)" % (ri["r_tip"], x_tip_dep))
    ax.set_aspect("equal")
    ax.set_xlim(-10, 175)
    ax.set_ylim(-fins["r_tip_max_mm"] - 8, fins["r_tip_max_mm"] + 8)
    ax.set_xlabel("x (mm), flow ->")
    ax.set_ylabel("y (mm)")
    ax.set_title("Roll fins, plan view from above (fins on the front shell)", loc="left")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    # end view: blades folding at the stowing index, points inside the fin x range
    fc = fins["clearance"].get("fold_at_index") or {}
    if fins.get("fold_rule") == "index":
        C = blade_fold_cloud()
        m = (C["x"] >= fins["x_le_mm"]) & (C["x"] <= fins["rail"]["pad_x_mm"][1])
        idx = np.nonzero(m)[0]
        if idx.size > 25000:
            idx = idx[:: idx.size // 25000 + 1]
        th_c = np.radians(C["theta"][idx])
        a2.scatter(C["rho"][idx] * np.sin(th_c), C["rho"][idx] * np.cos(th_c), s=0.3, color=COL_PK, alpha=0.25,
                   lw=0, label="blades folding 0-90 deg at the index")
        xc_ = fins["x_te_mm"]
        tt = np.linspace(0, 2 * np.pi, 721)
        R_ = float(r_out(xc_))
        a2.plot(R_ * np.sin(tt), R_ * np.cos(tt), color=COL_HULL, lw=1.0)
        rl2 = fins["root"]
        yw = fins["yaw"]
        for th_f in fins["theta_deg"]:
            a_ = math.radians(th_f)
            er = np.array([math.sin(a_), math.cos(a_)])
            et = np.array([math.cos(a_), -math.sin(a_)])
            # pad and fin drawn at the TE (largest yaw offset over the fin), half-widths at constant x
            off = (xc_ - yw["pivot_x_mm"]) * yw["offset_per_mm_x"]
            for r0, r1, w, col in ((R_, rl2["y_te"], yw["pad_half_width_at_constant_x_mm"], "#8a8a8a"),
                                   (rl2["y_te"], rl2["y_te"] + fins["span_mm"], fc.get("fin_half_thickness_mm", 2.0),
                                    "#2F3337")):
                a2.add_patch(Polygon([r0 * er + (off - w) * et, r1 * er + (off - w) * et, r1 * er + (off + w) * et,
                                      r0 * er + (off + w) * et], closed=True, fc=col, lw=0, alpha=0.85))
        a2.set_aspect("equal")
        lim_ = rl2["y_te"] + fins["span_mm"] + 5
        a2.set_xlim(-lim_, lim_)
        a2.set_ylim(-lim_, lim_)
        a2.set_xticks([])
        a2.set_yticks([])
        for sp in a2.spines.values():
            sp.set_visible(False)
        a2.set_title("Seen from aft, looking upstream (+Y right), x %.1f-%.1f: yawed pads and fins (at the TE)\n"
                     "vs blades folding 0-90 deg at the stowing index; 3D min clearance pad %.1f / fin %.1f mm"
                     % (fins["x_le_mm"], fins["rail"]["pad_x_mm"][1], fc.get("pad_clear_min_mm") or float("nan"),
                        fc.get("fin_clear_min_mm") or float("nan")), fontsize=8, loc="left")
    else:
        a2.set_visible(False)
    fig.tight_layout()
    p4 = FIG_DIR / "layout_section_fins.png"
    fig.savefig(p4, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    out.append(p4)
    return out


# ==========================================================================
# Trade study on the equator position
# ==========================================================================


def equator_trade(pk: Pockets, ri: dict, hw: dict) -> tuple[list, list]:
    rows = []
    mo = hw["motor"]
    for x_eq in (52.5, 55.0, 57.5, 60.0):
        J = joint(x_eq, pk, hw.get("oring"))
        bay = generator_bay(J, pk, ri, hw)
        fr = bay["frontier"]
        rows.append(dict(x_eq_mm=x_eq, oring=J["oring"]["designation"], seal_bore_d_mm=2 * J["r_bore"],
                         spigot_bore_d_mm=2 * J["r_spigot_in"], lip_d_mm=2 * J["r_lip"],
                         bayonet_engage_deg=J["bayonet"]["engage_deg"],
                         largest_envelope=bay["largest"],
                         D_max_at_bay_L_mm=float(np.interp(bay_length(hw), [f["L_mm"] for f in fr],
                                                           [f["D_max_mm"] for f in fr]))))
    sens = []
    for d in (2.0, 4.0, 6.0, 8.0, 10.0):
        ri2 = dict(pk.ri, depth=d, depth_x=None, depth_v=None, width_v=None)
        pk2 = Pockets(ri2)
        J = joint(P.EQUATOR_X, pk2, hw.get("oring"))
        bay = generator_bay(J, pk2, ri, hw)
        fr = bay["frontier"]
        sens.append(dict(pocket_depth_mm=d, oring=J["oring"]["designation"],
                         seal_bore_d_mm=2 * J["r_bore"], spigot_bore_d_mm=2 * J["r_spigot_in"],
                         D_max_at_bay_L_mm=float(np.interp(bay_length(hw), [f["L_mm"] for f in fr],
                                                           [f["D_max_mm"] for f in fr])),
                         bayonet_engage_deg=J["bayonet"]["engage_deg"],
                         bayonet_ok=J["bayonet"]["engage_ok"]))
    return rows, sens


# ==========================================================================
# Main
# ==========================================================================


def run(x_eq: float | None = None, trades: bool = True) -> dict:
    ri = rotor_inputs()
    hw = hardware_inputs()
    pk = Pockets(ri)
    x_eq = P.EQUATOR_X if x_eq is None else x_eq
    J = joint(x_eq, pk, hw.get("oring"))
    bay = generator_bay(J, pk, ri, hw)
    DT = drivetrain(J, bay, ri, hw, pk)
    fins = fin_design(ri, DT["stations_mm"]["rotor_load_plane"], pk)
    nose_d = nose(ri, fins, bay, hw)
    JS = joint_strength(J, nose_d["loads"]["design_load_N"])
    mb = mass_buoyancy(J, bay, DT, fins, ri, hw, pk, nose_d)
    roll = roll_summary(ri, fins, mb)
    trade, sens = equator_trade(pk, ri, hw) if trades else ([], [])
    return dict(inputs=dict(rotor=ri, hardware=hw), joint=J, joint_strength=JS, generator_bay=bay,
                drivetrain=DT, nose=nose_d, fins=fins, roll=roll, mass=mb, equator_trade=trade,
                pocket_sensitivity=sens, assembly=assembly(J, bay, DT, nose_d, fins, hw, ri.get("rotor")),
                parts=part_split(J, bay, fins, nose_d, ri), _pk=pk)


def _clean(o):
    """JSON-safe copy: numpy -> python, tuples -> lists, NaN -> None, floats rounded."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_clean(v) for v in o.tolist()]
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        if not math.isfinite(f):
            return None
        return round(f, 4)
    return o


def pocket_depth_limit(ri: dict, hw: dict, D_need: float, L: float | None = None) -> float:
    """Largest uniform pocket depth at the equator for which the bay still holds D_need at length L
    (default: the reserved length; bisection on the fallback pocket model)."""
    L = bay_length(hw) if L is None else L

    def dmax(d):
        pk2 = Pockets(dict(ri, depth=d, depth_x=None, depth_v=None, width_v=None,
                           x0=POCKET_FALLBACK["x0"], x1=POCKET_FALLBACK["x1"]))
        J = joint(P.EQUATOR_X, pk2, hw.get("oring"))
        fr = generator_bay(J, pk2, ri, hw)["frontier"]
        return float(np.interp(L, [f["L_mm"] for f in fr], [f["D_max_mm"] for f in fr]))

    lo, hi = 2.0, 12.0
    if dmax(lo) < D_need:
        return 0.0
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dmax(mid) >= D_need else (lo, mid)
    return lo


def summarize(res: dict) -> dict:
    """Curated, unit-suffixed output for data/layout_internal.json."""
    J, JS, bay, DT = res["joint"], res["joint_strength"], res["generator_bay"], res["drivetrain"]
    nose_d, fins, roll, mb = res["nose"], res["fins"], res["roll"], res["mass"]
    ri, hw = res["inputs"]["rotor"], res["inputs"]["hardware"]
    o, b, g, st, br = J["oring"], J["bayonet"], bay["pick"], DT["stations_mm"], DT["bearings"]
    mo = hw["motor"]
    d_lim = pocket_depth_limit(ri, hw, g["D_mm"], g["L_mm"])
    trade_L = bay["length_trade"]
    t33 = next((t for t in trade_L if abs(t["L_mm"] - mo["L_reserve"]) < 1e-6), None)
    fm = fins["margins"]
    tube, cn, gl_, rl_ = DT["tube"], nose_d["gland"]["cable_notch"], nose_d["gland"], fins["root"]
    x_car0 = g["x_fwd_mm"] - CARRIER_T
    pk_ = res["_pk"]
    proposed = dict(
        EQUATOR_X=J["x_eq"],
        BEARING_FRONT=dict(P.BEARING, x0=st["bearing_front_608"][0], role="front, locating"),
        BEARING_REAR=dict(P.BEARING_ALT, x0=st["bearing_rear_688"][0], role="rear, floating"),
        LIP_SEAL=dict(P.LIP_SEAL, x0=st["lip_seal"][0], recess=SEAL_RECESS),
        BEARING_TUBE=dict(
            x0=tube["x_from"], x1=tube["x_to"], r_out_front=tube["r_out_front_mm"], r_out_rear=tube["r_out_rear_mm"],
            profile_x=[tube["x_from"], st["bearing_front_608"][1] + 2, st["bearing_rear_688"][0] - 2,
                       st["bearing_rear_688"][0], tube["x_to"]],
            profile_r=[tube["r_out_front_mm"], tube["r_out_front_mm"], tube["r_out_rear_mm"] + 1.0,
                       tube["r_out_rear_mm"], tube["r_out_rear_mm"]],
            seat_front_d=P.BEARING["D"], shoulder_front_id=FRONT_SHOULDER_ID,
            shoulder_front_x=[st["bearing_front_608"][1], st["bearing_rear_688"][0]],
            seat_rear_d=P.BEARING_ALT["D"], seal_shoulder_x=st["seal_shoulder"], seal_shoulder_id=SEAL_SHOULDER_ID,
            seal_bore_d=P.LIP_SEAL["D"], mouth_lead_in=TUBE_LEAD_IN,
            seat_note="printed seats = nominal D + 2 x PRESS_FIT; the 18 mm bore between the seats lets the 688 "
                      "pass through the 608 seat",
            ribs=dict(n=3, t=2.0, theta_deg=sorted(pk_.theta), x=[tube["x_from"], tube["x_to"]],
                      rule="radial ribs under the pocket floors, tube to shell")),
        TAIL_END_BLOCK=dict(x0=st["bearing_rear_688"][1], x1=ri["static_tail_x"],
                            note="solid ring: seal shoulder, seal bore, the 3 notch floors and the stowed-detent "
                                 "magnet seats (layout_rotor.json detents.stowed)"),
        GENERATOR=dict(name=mo["name"], D=mo["D"], L=mo["L"], mount_pcd=mo["pcd"], mount_screw="M3",
                       n_mount=mo["n_screw"], rotor_face_pcd=mo.get("rotor_pcd"),
                       bay_D=g["D_mm"], bay_L=g["L_mm"], x_bay_fwd=g["x_fwd_mm"], x_rotor_face=g["x_aft_mm"],
                       motor_body_x=g["motor_body_x_mm"], adapter_spacer_L=g["adapter_spacer_L_mm"],
                       alternative=mo.get("alt_name"), alt_mount_pcd=mo.get("alt_pcd")),
        CARRIER=dict(x0=x_car0, t=CARRIER_T, r=g["r_carrier_mm"], twist_deg=30.0, centre_bore_d=16.0,
                     screw_pattern_clock_deg=45.0,
                     ring=dict(x=[x_car0 - CARRIER_RING_OVER, g["x_fwd_mm"] + CARRIER_RING_OVER], t_radial=CARRIER_RING,
                               r_in_at_plate=cn["carrier_ring_inner_r_mm"]),
                     cable_notch=dict(theta_deg=cn["theta_deg"], width=cn["width_mm"], depth=cn["depth_from_rim_mm"],
                                      bottom_r=cn["bottom_r_mm"])),
        COUPLER=dict(od=COUPLER_OD, claw_a_L=CLAW_A_L, float=CLAW_FLOAT, claw_b_L=CLAW_B_L,
                     claw_a_x=st["claw_A"], claw_b_x=st["claw_B_thrust_collar"],
                     claw_b_step_od=CLAW_B_STEP_OD, claw_b_step_L=CLAW_B_STEP_L, shim_od=SHIM_OD, shim_t=0.5,
                     thrust_collar="Ruland MSC-8-SS 8x16x8"),
        SHAFT=dict(d=P.SHAFT_D, x0=st["shaft"][0], x1=st["shaft"][1], L=st["shaft"][1] - st["shaft"][0],
                   hub_flat_x=(ri.get("rotor") or {}).get("grub_x"), hub_flat_depth=0.5, aft_tap="M3 x 8"),
        JOINT=dict(x_eq=J["x_eq"], x_spigot_tip=J["x_tip"], x_lip=[J["x_lip0"], J["x_seal0"]],
                   x_seal_zone=[J["x_seal0"], J["x_eq"]], x_channel=[J["x_ch0"], J["x_lip0"]],
                   x_ring_ceiling_start=J["x_ring0"], groove_x=J["groove"],
                   seal_bore_d=2 * J["r_bore"], spigot_seal_od=2 * J["r_spigot_seal"],
                   groove_bottom_d=2 * J["r_groove"], lug_od=2 * J["r_lug"], spigot_lug_zone_od=2 * J["r_spigot_lug"],
                   lip_d=2 * J["r_lip"], channel_d=2 * J["r_channel"], spigot_bore_d=2 * J["r_spigot_in"],
                   root_ring_L=J["root_ring_L"], land_root=OR_LAND_ROOT, land_tip=OR_LAND_TIP),
        ORING_EQUATOR=dict(cs=o["cs_mm"], name=o["designation"], id=o["id_free_mm"],
                           groove_d=o["groove_bottom_d_mm"], bore_d=o["bore_d_mm"], spigot_od=o["spigot_seal_od_mm"],
                           groove_bottom_w=o["groove_bottom_w_mm"], groove_top_w=o["groove_top_w_mm"],
                           groove_depth=o["groove_depth_in_spigot_mm"], gland_depth=o["gland_depth_mm"],
                           flank="aft flank radial, tip-side flank 45 deg", radial_clearance=o["radial_clearance_mm"]),
        BAYONET=dict(n_lugs=LUG_N, lug_w_deg=LUG_W_DEG, lug_h=LUG_H, lug_t=LUG_T, twist_deg=b["engage_deg"],
                     lock_theta_deg=b["lock_theta_deg"], entry_theta_deg=b["entry_theta_deg"], gap_w_deg=GAP_W_DEG,
                     lip_t=LIP_T, channel_extra=CH_EXTRA),
        FIN=dict(section=FIN_SECTION, chord=fins["chord_mm"], span=fins["span_mm"],
                 incidence_deg=fins["incidence_deg"], x_le=fins["x_le_mm"], x_te=fins["x_te_mm"],
                 theta_deg=fins["theta_deg"], yaw_pivot_chord_frac=fins["yaw"]["pivot_chord_frac"],
                 yaw_pivot_x=fins["yaw"]["pivot_x_mm"], te_side=fins["yaw"]["te_side"],
                 root_r_le=rl_["y_le"], root_r_te=rl_["y_te"], rail_w=RAIL_W, rail_h_min=RAIL_H,
                 groove_w_mean=GROOVE_W, groove_depth=RAIL_H, groove_form="half dovetail, 45 deg flank",
                 pad_x=fins["rail"]["pad_x_mm"], end_stop_x=fins["x_te_mm"],
                 retention="M3 grub in a heat-set insert", trim_variant_deg=FIN_TRIM_DEG),
        GLAND=dict(thread="M12x1.5", x=gl_["seat_x_mm"], theta_deg=GLAND_THETA,
                   tilt_from_normal_deg=gl_["axis_tilt_from_normal_deg"],
                   axis_angle_from_x_deg=gl_["axis_angle_from_x_deg"],
                   seat_d=GLAND_SEAT_D, boss_d=GLAND_BOSS_D, seat_raise=gl_["seat_raise_mm"],
                   inner_end_x=gl_["inner_end_x_mm"], inner_end_r=gl_["inner_end_r_mm"]),
        TOW_EYE=dict(hole_d=EYE_HOLE_D, ring_w=EYE_W, ring_t=EYE_T, plane=EYE_PLANE,
                     hole_centre_x=nose_d["eye"]["hole_centre_x_mm"]),
        NOSE_PLUG=dict(x_45=nose_d["nose_plug"]["x_45_mm"], bridge_d=nose_d["nose_plug"]["bridge_d_mm"]),
        SPINNER_ACCESS=dict(d=SPINNER_ACCESS_D, x=mb["rotor_wet"]["access_bore_x_mm"], aft_wall=SPINNER_AFT_WALL),
        BALLAST_V2=dict(lead_g=mb["v2_ballast"]["lead_g"], x=mb["v2_ballast"]["x_mm"], theta_deg=180.0),
        POCKET_FLOOR_WALL=POCKET_FLOOR_WALL,
    )
    fc = fins["clearance"]
    alt_f = fins.get("alternative_rule") or {}
    reasons = dict(
        EQUATOR_X_reason="kept: the bay is limited by the bearing span and the nose (carrier), not by "
                         "x_eq; 57.5 keeps the seal in the parallel band",
        BEARING_reason="608 front at x %.1f (takes the thrust, wide hull), 688 rear at x %.1f under the tail "
                       "notch (passes through the 608 seat so both press in from the front); span %.1f mm"
                       % (st["bearing_front_608"][0], st["bearing_rear_688"][0], br["spacing_mm"]),
        GENERATOR_reason="bay = the largest envelope that fits, D%g x L%g at x %.2f-%.2f (the nose wall and the "
                         "carrier ring set D, the bearing span sets L); it takes %s. GL40 (hardware-study pick) "
                         "stays the recommended motor. The hardware study's L%g reserve would cut the span to "
                         "%.1f mm (< %g)" % (g["D_mm"], g["L_mm"], g["x_fwd_mm"], g["x_aft_mm"],
                                             ", ".join(c_["name"] for c_ in bay["candidate_fits"] if c_["fits_bay"]),
                                             mo["L_reserve"],
                                             t33["span_mm"] if t33 and t33.get("fits") else float("nan"), SPAN_MIN),
        ORING_EQUATOR_reason="the hardware study's final pick, adopted as is: %s, groove bottom D%.2f (stretch %.1f %%), "
                             "bore D%.2f, spigot OD %.2f, bottom width %.2f with a 45 deg tip-side flank (top %.3f); "
                             "the bore fits the wall limit D%.2f (%.2f mm of shell under the %.1f mm deep pockets at "
                             "the joint); recomputed here: squeeze %.1f %% nominal, %.1f-%.1f %% worst, fill %.1f %% "
                             "nominal, %.1f %% worst"
                             % (o["designation"], o["groove_bottom_d_mm"], o["stretch_nom_pct"], o["bore_d_mm"],
                                o["spigot_seal_od_mm"], o["groove_bottom_w_mm"], o["groove_top_w_mm"],
                                o["bore_limit_d_mm"], J["walls"]["socket_at_pocket_mm"],
                                float(res["_pk"].depth(J["x_eq"])), o["squeeze_nom_pct"], *o["squeeze_worst_pct"],
                                o["fill_nom_pct"], o["fill_worst_pct"]),
        BAYONET_reason="3 x %g deg lugs locked at theta %s, %g deg twist; the twist channel stays %.1f mm clear of "
                       "the asymmetric pockets; lug shear util %.2f at the tow design load"
                       % (LUG_W_DEG, "/".join("%.1f" % v for v in b["lock_theta_deg"]), b["engage_deg"],
                          b["channel_to_pocket_edge_mm"] or float("nan"), JS["util_lug"]),
        FIN_reason="smallest tip-to-tip fin pair whose couple at %.1f m/s is %.3f x the %s load-line torque "
                   "(%.4f N m at TSR %.2f; %.0f %% margin) with CL <= %.2f CLmax_3D, and whose YAWED rail pads (pad "
                   "yawed with the fin about the pivot at %.2f c, TE toward %s) clear the pockets; fins on the front "
                   "shell (TE x %.1f) at theta %g/%g, pad-to-pocket %.2f mm; blades folding at the stowing index "
                   "clear pads and fins by %.2f mm. The any-azimuth fold-torus rule would force %g x %g fins at x "
                   "%.1f-%.1f, %.0f mm tip to tip"
                   % (FIN_DESIGN_V, fins["torque_basis"]["ratio"], fins["torque_basis"]["generator"],
                      fins["torque_basis"]["Q_design_Nm"], fins["torque_basis"]["tsr"] or float("nan"),
                      fins["torque_basis"]["margin_pct"], FIN_CL_RESERVE,
                      fins["yaw"]["pivot_chord_frac"], fins["yaw"]["te_side"],
                      fins["x_te_mm"], *fins["theta_deg"], fc["stowed_pockets"]["min_clear_mm"] or float("nan"),
                      min(v for v in ((fc.get("fold_at_index") or {}).get("pad_clear_min_mm"),
                                      (fc.get("fold_at_index") or {}).get("fin_clear_min_mm")) if v is not None)
                      if fc.get("fold_at_index") else float("nan"),
                      alt_f.get("chord_mm", float("nan")), alt_f.get("span_mm", float("nan")),
                      alt_f.get("x_le_mm", float("nan")), alt_f.get("x_te_mm", float("nan")),
                      alt_f.get("tip_to_tip_mm", float("nan"))),
        BEARING_TUBE_reason="aft-shell tube: 608 seat at the mouth, 18 mm bore, 688 seat, seal shoulder and seal "
                            "bore, one stepped bore pressed from the front; 3 ribs under the pocket floors",
        TAIL_END_BLOCK_reason="solid end of the aft shell behind the 688 (seal boss, notch floors, magnet seats)",
        JOINT_reason="radial stack and stations of the bayonet + O-ring joint from joint(); spigot on the aft "
                     "shell, socket in the front shell",
        NOSE_PLUG_reason="forward of x %.1f the nose inner wall leans past 45 deg (printed nose-up): infilled "
                         "plug closed by a %.1f mm bridge" % (nose_d["nose_plug"]["x_45_mm"],
                                                               nose_d["nose_plug"]["bridge_d_mm"]),
        SPINNER_ACCESS_reason="the spinner aft of the hub face is printed solid (closed perimeters + infill) "
                              "around a D%g screw-access bore that floods; the buoyancy counts it wet" % SPINNER_ACCESS_D,
        BALLAST_V2_reason="v2 only, to reach +0.3 N net buoyancy if needed (0 g means the hull is already below "
                          "the target and needs float, see mass_buoyancy.v2_needs)",
        CARRIER_reason="twist-lock bulkhead at the forward limit set by the PCD 38.14 adapter; screw pattern "
                       "clocked 45 deg so the cable notch at theta 180 misses the counterbores",
        COUPLER_reason="blind self-finding claw across the joint; claw B is the thrust collar and steps to OD "
                       "%g over its last %g mm (%.2f mm radial clearance in the static tube mouth)"
                       % (CLAW_B_STEP_OD, CLAW_B_STEP_L, DT["coupler"]["tube_mouth"]["min_radial_clear_mm"]),
    )
    fr_L = [f_ for f_ in bay["frontier"] if f_["D_max_mm"] > 0]
    decisions = [
        dict(topic="generator bay", decision="motor in the front shell on a twist-lock carrier, rotor face "
             "aft at x %.2f; bay = the largest envelope that fits, D%g x L%g from x %.2f (frontier: D %.2f at "
             "L %g, longest bay %g mm because the rotor face may not pass x %.2f); motors that fit it: %s"
             % (g["x_aft_mm"], g["D_mm"], g["L_mm"], g["x_fwd_mm"], fr_L[-1]["D_max_mm"], fr_L[-1]["L_mm"],
                fr_L[-1]["L_mm"], bay["x_bell_max_mm"],
                ", ".join(c_["name"] for c_ in bay["candidate_fits"] if c_["fits_bay"])),
             why="uses the empty nose volume, keeps every wire in the front shell (no connector across the "
                 "joint) and leaves the tail to the bearings. With the static tail at x %.0f the generator "
                 "length sets the bearing span: forward-most placement maximises it" % ri["static_tail_x"]),
        dict(topic="drivetrain", decision="608-2RS at x %.2f-%.2f (locating), 688-2RS at %.1f-%.1f "
             "(floating), lip seal %.1f-%.1f, shaft 8 x %.2f from x %.2f to %.1f (hub on x %.1f-%.1f)"
             % (*st["bearing_front_608"], *st["bearing_rear_688"], *st["lip_seal"],
                st["shaft"][1] - st["shaft"][0], *st["shaft"], ri["rotor"]["hub_front_x"] if ri.get("rotor")
                else ri["static_tail_x"] + P.HUB_GAP, st["shaft"][1]),
             why="span %.1f mm vs overhang %.1f mm (span/overhang %.2f, below the %.2f rule of thumb, which no "
                 "position holding the GL40 meets); explicit checks pass: tip runout %.2f mm from 0.05 mm "
                 "seat slop, shaft deflection %.3f mm, static safety s0 %.0f / %.0f under the tow load"
                 % (br["spacing_mm"], br["overhang_mm"], br["spacing_over_overhang"], SPAN_RATIO_REF,
                    br["tip_axial_runout_mm"], br["shaft_deflection_at_load_plane_mm"],
                    br["static_safety_front_s0"], br["static_safety_rear_s0"])),
        dict(topic="coupling", decision="two-dog claw coupling across the joint; claw B is a printed "
             "sleeve over a Ruland MSC-8-SS collar and is the thrust collar on the 608 inner race",
             why="engages blind when the shells close, tolerates shell misalignment (no 4-bearing "
                 "shaft), no machining of a hardened shaft"),
        dict(topic="equator joint", decision="%s O-ring on the aft-shell spigot, 3 bayonet lugs %.0f deg "
             "wide locked at theta %s, %.0f deg twist, 45 deg load faces, 2 magnet pairs as detent"
             % (o["designation"], LUG_W_DEG, "/".join("%.1f" % v for v in b["lock_theta_deg"]), b["engage_deg"]),
             why="seal on the wet side of the lugs, everything prints without supports, lugs in the lands "
                 "between the (asymmetric) pockets so the channel never cuts a pocket floor"),
        dict(topic="nose", decision="printed tow eye (8 mm hole, 5 x 6 ring, ring plane %s) on the axis; M12 "
             "gland beside it on the bottom, seat at x %.2f, axis tilted %g deg from the surface normal toward the "
             "nose; tip plug forward of x %.1f; cable to the motor through a %.1f mm notch in the carrier rim"
             % (EYE_PLANE, nose_d["gland"]["seat_x_mm"], nose_d["gland"]["axis_tilt_from_normal_deg"],
                nose_d["nose_plug"]["x_45_mm"], nose_d["gland"]["cable_notch"]["width_mm"]),
             why="eye legs at %.2f of the PETG interlayer adhesion under 3 x the tow load; the whole gland boss "
                 "stays %.2f mm ahead of the carrier plate and %.2f mm ahead of its ring, the gland body %.1f mm "
                 "from the eye and %.1f mm from the shackle envelope (the earlier boss at x 12 reached x %.1f, "
                 "through the carrier into the bay)"
                 % (nose_d["eye"]["materials"]["PETG"]["eye_leg_util"], nose_d["gland"]["boss_to_carrier_plate_mm"],
                    nose_d["gland"]["boss_to_carrier_ring_mm"], nose_d["gland"]["body_to_eye_ring_mm"],
                    nose_d["gland"]["body_to_shackle_envelope_mm"],
                    nose_d["gland"]["previous_x12"]["boss_x_max_mm"])),
        dict(topic="roll", decision="two identical slide-in Clark Y fins %.0f x %.0f mm at theta %g / %g, incidence "
             "%.2f deg, on raised dovetail pads yawed with the fins (pivot %.2f c), fin x %.1f-%.1f, pads to x %.1f"
             % (fins["chord_mm"], fins["span_mm"], *fins["theta_deg"], fins["incidence_deg"],
                fins["yaw"]["pivot_chord_frac"], fins["x_le_mm"], fins["x_te_mm"], fins["rail"]["pad_x_mm"][1]),
             why="sized on the %s load line: at %.1f m/s the couple is %.3f x its %.4f N m (%.0f %% margin, see "
                 "fins.torque_basis); off design the ratio is %.3f at 1.0 m/s and %.3f at 2.5 m/s because the load "
                 "line's TSR moves with the flow; +-%g deg trim fins change it by +-%.0f %%. Blades fold at the "
                 "stowing index, so the fins only need to clear that sweep"
                 % (fins["torque_basis"]["generator"], FIN_DESIGN_V, fm[1.5]["ratio"],
                    fins["torque_basis"]["Q_design_Nm"], fins["torque_basis"]["margin_pct"], fm[1.0]["ratio"],
                    fm[2.5]["ratio"], FIN_TRIM_DEG, 100 * fins["trim"]["ratio_change_per_deg_at_design_flow"])),
        dict(topic="part split", decision="front shell (nose up) and aft shell with integral tube, ribs "
             "and spigot (tail down): no supports on either shell", why="both fit the 256 mm bed; the fin pads run "
             "down to the equator face (no ledge) and the pocket floors rise at most %.0f deg from vertical "
             "(rotor half: printable_floor)" % (ri.get("floor_overhang_max_deg") or float("nan"))),
    ]
    open_issues = [
        "Bearing span %.1f mm against a %.1f mm overhang to the rotor thrust centroid (x %.1f on the solved "
        "blade axis at 0.7 R): span/overhang %.2f, below the %.2f rule of thumb. Explicit checks pass (tip "
        "runout %.2f mm, shaft deflection %.3f mm, s0 >= %.0f), but measure rotor wobble on the prototype. "
        "Only the rotor side could restore the ratio (a longer static tail), and its hinge sweep rejects "
        "STATIC_TAIL_X 110 (root clearance 0.77 < 0.8 mm)."
        % (br["spacing_mm"], br["overhang_mm"], br["rotor_load_plane"] if "rotor_load_plane" in br else
           st["rotor_load_plane"], br["spacing_over_overhang"], SPAN_RATIO_REF, br["tip_axial_runout_mm"],
           br["shaft_deflection_at_load_plane_mm"], min(br["static_safety_front_s0"], br["static_safety_rear_s0"])),
        "Bay length: the bay is D%g x L%g (largest that fits). The hardware study's L%g reserve (GM4108H "
        "encoder variant, GL40 II, GBM4114HS) does not fit with a %g mm span; it would leave %.1f mm."
        % (g["D_mm"], g["L_mm"], mo["L_reserve"], SPAN_MIN, t33["span_mm"] if t33 and t33.get("fits") else float("nan")),
        "Roll: the fins carry the mean torque but nothing restores roll angle (pendulum %.4f N m vs Q %.3f N m "
        "at 1.5 m/s). If the models are exact, the fin couple is %+.0f %% of Q at 1.5 m/s (hull rolls %.1f rpm %s), "
        "%+.0f %% at 1.0 m/s (%.1f rpm %s) and %+.0f %% at 2.5 m/s (%.1f rpm %s). v2 needs a bridle on two nose eyes "
        "(heel %.1f deg at 1.5 m/s) or a cable swivel; trim fins (+-%g deg) set the balance after the first wet test."
        % (roll["pendulum_righting_max_Nm"], roll["torque"]["1.5_m_s"]["Q_Nm"],
           *[w for k in ("1.5_m_s", "1_m_s", "2.5_m_s")
             for w in (roll["balance"][k]["excess_pct"], abs(roll["balance"][k]["free_roll_rate_rpm"]),
                       "against the rotor" if roll["balance"][k]["excess_pct"] > 0 else "with the rotor")],
           roll["balance"]["1.5_m_s"]["bridle_heel_deg"], FIN_TRIM_DEG),
        "Fins: %.0f mm tip to tip (rotor %.0f mm), removable for transport. They clear the blades only when "
        "the blades fold at the stowing index; folding at another rotor azimuth with the fins on can touch "
        "a fin (assembly guide: index the rotor, or remove the fins, before folding). Margins: yawed pad to "
        "pocket %.2f mm, folding blades to pad %.2f mm and to fin %.2f mm (limits %.1f); the pad footprint, the "
        "pocket raster and the fold sweep (1 deg steps) are sampled, so confirm them on the CAD."
        % (fins["tip_to_tip_mm"], 2 * ri["r_tip"], fins["clearance"]["stowed_pockets"]["min_clear_mm"],
           (fins["clearance"].get("fold_at_index") or {}).get("pad_clear_min_mm") or float("nan"),
           (fins["clearance"].get("fold_at_index") or {}).get("fin_clear_min_mm") or float("nan"), FIN_POCKET_CLEAR),
        "v2 buoyancy: net %+.2f N without ballast (blade slots, spinner access bore and gap counted wet); "
        "reaching +%.1f N needs %.0f g less mass (4-perimeter shells save %.0f g) or ~%.0f cc of external "
        "float." % (mb["v2_needs"]["net_buoyancy_no_ballast_N"], mb["v2_needs"]["target_net_buoyancy_N"],
                    mb["v2_needs"]["mass_to_shed_or_float_g"], mb["v2_needs"]["option_4_perimeter_shells_saves_g"],
                    mb["v2_needs"]["option_external_foam_cc"]),
        "O-ring worst-case squeeze %.0f %% is below the 14 %% table minimum (print error + full "
        "eccentricity); fine at < 0.1 bar, verify with a printed test ring." % o["squeeze_worst_pct"][0],
        "Assembly flags: shaft flat passing the lip seal (sleeve), blind claw engagement on closing.",
        "Lug and lip shear use the interlayer tensile adhesion as a proxy (shear strength across layers "
        "not sourced); utilisation %.2f." % max(JS["util_lug"], JS["util_lip"]),
        "Hull drag C_D 0.10, radial rotor load 10 % of thrust, one blade 20 % off its share, thrust centroid "
        "at 0.7 R: ESTIMATES.",
    ]

    def fins_margin(v):
        m = fm.get(v) if v in fm else fm.get(str(v))
        return dict(v_fin_m_s=m["v_fin_m_s"], Re=m["Re"], CL=m["CL"], CLmax_3D=m["CLmax3D"],
                    fin_torque_Nm=m["torque_Nm"], generator_torque_Nm=m["Q_gen_Nm"], torque_ratio=m["ratio"],
                    torque_ratio_per_deg_incidence=m["ratio_per_deg"],
                    authority_ratio_at_CL_reserve=m["authority_ratio"], stall_reserve=m["stall_reserve"],
                    roll_damping_Nm_s=m["roll_damping_Nm_s"], neuralfoil_confidence_min=m["nf_confidence_min"])

    out = dict(
        schema="hydroegg.layout_internal/1",
        generated_by="python -m hydroegg.layout_internal",
        conventions=dict(frame="X body axis, flow +X, nose at x = 0, Z up, right-handed",
                         theta="angle about +X from +Z toward +Y; theta = -psi of layout_rotor (blade i at psi "
                               "120 i, i.e. theta 0/240/120); pockets span theta_blade + [-az_hi, -az_lo]",
                         units="mm and deg in geometry, SI in physics; keys carry the unit",
                         estimates="any value marked ESTIMATE is a modelling assumption, not a sourced number"),
        inputs=dict(pockets_source=ri["source"],
                    pocket_model=dict(x0_mm=ri["x0"], x1_mm=ri["x1"], depth_max_mm=ri["depth"],
                                      width_mm=ri["width"], theta_deg=list(ri["theta"])),
                    blade_source=ri["blade_source"], hardware_source=hw["source"],
                    r_tip_mm=ri["r_tip"], cone_deg=ri["cone_deg"], static_tail_x_mm=ri["static_tail_x"],
                    hinge_x_mm=ri["hinge_x"], hinge_r_mm=ri["hinge_r"], rotor_omega=ri["omega_vector"]),
        decisions=decisions,
        joint=dict(
            x_eq_mm=J["x_eq"], seal_bore_d_mm=2 * J["r_bore"], groove_bottom_d_mm=2 * J["r_groove"],
            spigot_seal_od_mm=2 * J["r_spigot_seal"], lug_od_mm=2 * J["r_lug"], spigot_lug_zone_od_mm=2 * J["r_spigot_lug"],
            lip_d_mm=2 * J["r_lip"], channel_d_mm=2 * J["r_channel"], spigot_bore_d_mm=2 * J["r_spigot_in"],
            stations_x_mm=dict(spigot_tip=J["x_tip"], lip=[J["x_lip0"], J["x_seal0"]],
                               seal_zone=[J["x_seal0"], J["x_eq"]], channel=[J["x_ch0"], J["x_lip0"]],
                               ring_ceiling_start=J["x_ring0"], groove=J["groove"],
                               root_ring=[J["x_eq"], J["x_eq"] + J["root_ring_L"]]),
            spigot_length_mm=J["L_spigot"], bore_limit_d_mm=J["bore_limit_d_mm"],
            oring=o, bayonet=b, walls_mm=J["walls"],
            strength=dict(design_load_N=JS["design_load_N"], lug_root_shear_MPa=JS["tau_lug_root_MPa"],
                          lip_root_shear_MPa=JS["tau_lip_root_MPa"], contact_pressure_MPa=JS["contact_pressure_45deg_MPa"],
                          interlayer_adhesion_PETG_MPa=JS["interlayer_adhesion_PETG_MPa"],
                          util_lug=JS["util_lug"], util_lip=JS["util_lip"], ok=JS["ok"],
                          basis="tow design load passes eye -> front shell -> lugs -> aft shell; "
                                "interlayer adhesion used as a shear proxy, util < 0.5 required"),
        ),
        generator_bay=dict(
            pick=dict(bay_D_mm=g["D_mm"], bay_L_mm=g["L_mm"], x_fwd_mm=g["x_fwd_mm"], x_rotor_face_mm=g["x_aft_mm"],
                      motor=g["motor"], motor_D_mm=mo["D"], motor_L_mm=mo["L"], motor_mass_g=mo["mass_g"],
                      motor_body_x_mm=g["motor_body_x_mm"], adapter_spacer_L_mm=g["adapter_spacer_L_mm"],
                      clear_to_spigot_bore_mm=g["clear_to_spigot_bore_mm"], rule=g["rule"]),
            largest_envelope=dict(D_mm=bay["largest"]["D_mm"], L_mm=bay["largest"]["L_mm"],
                                  x_fwd_mm=bay["largest"]["x_fwd_mm"], x_aft_mm=bay["largest"]["x_aft_mm"],
                                  rule=bay["largest"]["rule"]),
            frontier=[dict(L_mm=f_["L_mm"], D_max_mm=f_["D_max_mm"], x_fwd_mm=f_.get("x_fwd"),
                           x_aft_mm=f_.get("x_aft")) for f_ in bay["frontier"] if f_["D_max_mm"] > 0],
            frontier_note="the bay runs from the carrier (needs r >= %.2f mm for the PCD %.2f adapter, so x_fwd "
                          ">= %.2f) to the rotor face, which may not pass x %.2f (bearing span >= %g mm): hence L "
                          "<= %g. D_max is set by the nose wall and the carrier ring (it overhangs the bay by %.1f mm) at "
                          "the bay's forward end (socket lip D %.2f is not binding); the spigot bore (D %.2f) only matters aft of the spigot tip (x %.2f)"
                          % (bay["carrier_r_min_mm"], max(mo["pcd"], mo.get("alt_pcd") or 0), g["x_fwd_mm"],
                             bay["x_bell_max_mm"], SPAN_MIN,
                             max(f_["L_mm"] for f_ in bay["frontier"] if f_["D_max_mm"] > 0), CARRIER_RING_OVER,
                             2 * J["r_lip"], 2 * J["r_spigot_in"], J["x_tip"]),
            length_trade=bay["length_trade"],
            rotor_face_x_max_mm=bay["x_bell_max_mm"],
            candidate_fits=bay["candidate_fits"],
            carrier=dict(type="twist-lock bulkhead plate (30 deg, tabs under an inward ring of the front "
                              "shell), motor-specific adapter / spacer between it and the stator face",
                         x_mm=[g["x_fwd_mm"] - CARRIER_T, g["x_fwd_mm"]], r_mm=g["r_carrier_mm"],
                         r_required_mm=bay["carrier_r_min_mm"],
                         holes="motor pattern clocked 45 deg so no screw sits at theta 180 (cable); "
                               "cable notch in the rim at theta 180 (%.1f wide, %.1f deep, bottom r %.2f, %.1f mm "
                               "from the nearest screw counterbore); centre bore 16 mm"
                               % (nose_d["gland"]["cable_notch"]["width_mm"],
                                  nose_d["gland"]["cable_notch"]["depth_from_rim_mm"],
                                  nose_d["gland"]["cable_notch"]["bottom_r_mm"],
                                  nose_d["gland"]["cable_notch"]["clear_to_screw_counterbores_mm"]),
                         cable_notch=nose_d["gland"]["cable_notch"]),
        ),
        drivetrain=dict(
            stations_x_mm=dict(st, rotor_face=g["x_aft_mm"]),
            bearings=dict({k: v for k, v in br.items() if k not in ("ratings",)},
                          ratings=br["ratings"]),
            span_trade=DT["span_trade"], tube=DT["tube"], seal=DT["seal"],
            shaft=dict(d_mm=DT["shaft"]["d_mm"], x_fwd_mm=DT["shaft"]["x_fwd"], x_aft_mm=DT["shaft"]["x_aft"],
                       length_mm=DT["shaft"]["length_mm"], mass_g=DT["shaft"]["mass_g"],
                       aft_end_rule=DT["shaft"]["aft_end_rule"], stock=DT["shaft"]["stock"]),
            coupler=dict(DT["coupler"], thrust_collar="Ruland MSC-8-SS (bore 8, OD 16, width 8, 303 SS; "
                                                      "hardware.json shaft_retention) inside claw B, 0.5 mm "
                                                      "shim washer to the 608 inner race"),
        ),
        nose=nose_d,
        fins=dict({k: v for k, v in fins.items() if k not in ("margins", "sweep")},
                  margins={f"{v:.1f}_m_s": fins_margin(v) for v in FIN_SPEEDS},
                  sweep_feasible=[r_ for r_ in fins["sweep"] if r_["feasible"]]),
        roll=roll,
        mass_buoyancy=dict(model=mb["model"], display_PLA=mb["display_PLA"], v2_PETG=mb["v2"],
                           v2_needs=mb["v2_needs"], v2_ballast=mb["v2_ballast"],
                           static_hull_outer_cc=mb["static_hull_outer_cc"], pockets_cc=mb["pockets_cc"],
                           printed_parts=mb["parts"], hardware=mb["hardware"]),
        assembly=res["assembly"],
        part_split=res["parts"],
        trades=dict(equator_position=res["equator_trade"], pocket_depth_sensitivity=res["pocket_sensitivity"],
                    max_pocket_depth_at_equator_for_bay_mm=d_lim, bay_length=trade_L,
                    fin_fold_rule_alternative=fins.get("alternative_rule")),
        proposed_param_updates=dict(proposed, **reasons),
        open_issues=open_issues,
    )
    return out


# ==========================================================================
# Merge: one consistent layout (rotor half + internal half + interface checks)
# ==========================================================================

ROTOR_PARAM_REASONS = dict(
    STATIC_TAIL_X="rotor hinge sweep: 8 mm ahead of the hinge leaves a 4.9 mm stop block, keeps the "
                  "folding root >= 1.36 mm clear (front 6 gives 0.77) and a tail keep-out radius >= 14.2 mm",
    HINGE_X="rotor hinge sweep: only feasible HX with comfortable margins (117: 0.05 mm pin cover, 118 fails)",
    HINGE_R="rotor hinge sweep: best smallest margin (pin cover 2.46 mm in the hub core, stowed heel 1.35 mm "
            "inside the spinner, 4.9 mm heel face for a 4 mm magnet)",
    ROTOR_R_TIP_TARGET="unchanged (BEM owns it); stows with the tip at x 31.75, flush with the hull",
    PIN_L="14 mm knuckle + 2 x 0.3 clearance + 2 x 3.7 mm ears; ISO 8734 3 x 22, VERIFY stock (fallback "
          "3 x 20 with 2.7 mm ears)",
    FOLD_DEG="new: stowed fold angle; stowed fit is nearly fold-independent, 90 gives a 12.6 deg mean cone",
    KNUCKLE_D="new: blade root knuckle diameter (BEM lug thickness 8)",
    KNUCKLE_L="new: knuckle length along the pin (clevis slot 14.6)",
    HUB_EAR_T="new: hub ear thickness along the pin (ear bearing 4.0 MPa at 1.5 m/s, ESTIMATE)",
    SHANK_CONE_DEG="new: shank leaves the knuckle at 55 deg (raises the stowed root and notch floor ~1.5 mm)",
    HEEL_REACH="new: longest heel that keeps MIN_WALL of hub over the 8 mm bore",
    HEEL_T="new: heel thickness; 4.9 mm flat stop face holds a 4 x 2 magnet",
    HEEL_W="new: heel width along the pin",
    SHAFT_AFT_X="new: shaft end = hub aft face = HINGE_X + 12 (19.4 mm hub engagement); the internal "
                "layout uses it (no fallback)",
    POCKET_CLEAR="new: pocket floor clearance under the aero span",
    ROOT_CLEAR="new: stowed root to static tail clearance",
)


def _interval(name, x0, x1, owner, r=None, note=""):
    d = dict(part=name, x0_mm=float(x0), x1_mm=float(x1), owner=owner)
    if r is not None:
        d["r_mm"] = r
    if note:
        d["note"] = note
    return d


def merged(res: dict, summ: dict) -> dict:
    """data/layout.json: the rotor half (data/layout_rotor.json) and the internal half (this
    module) in one file, with the axial station list, the interface checks between the halves
    and one consolidated set of proposed params.py updates."""
    rj = _load(ROTOR_JSON)
    ri, hw = res["inputs"]["rotor"], res["inputs"]["hardware"]
    J, bay, DT, fins, nose_d = res["joint"], res["generator_bay"], res["drivetrain"], res["fins"], res["nose"]
    pk = res["_pk"]
    st, g = DT["stations_mm"], bay["pick"]
    ch, hub, hinge = rj["chosen"], rj["hub"], rj["hinge"]
    br = DT["bearings"]
    x_end = ri["static_tail_x"]
    mo = hw["motor"]
    eye_c = nose_d["eye"]["hole_centre_x_mm"]
    ring_R = EYE_HOLE_D / 2 + EYE_W

    stations = [
        _interval("tow eye (ring, hole centre at x %.1f)" % eye_c, eye_c - ring_R, eye_c + ring_R, "internal"),
        _interval("nose plug (infill, bridged)", 0.0, nose_d["nose_plug"]["x_45_mm"], "internal"),
        _interval("M12 gland boss (theta 180, seat x %.2f)" % nose_d["gland"]["seat_x_mm"],
                  float(np.min(gland_geom(nose_d["gland"]["seat_x_mm"], nose_d["gland"]["axis_tilt_from_normal_deg"],
                                          hw)["boss"][:, 0])), nose_d["gland"]["boss_x_max_mm"], "internal"),
        _interval("fins (theta %g / %g)" % tuple(fins["theta_deg"]), fins["x_le_mm"], fins["x_te_mm"], "internal"),
        _interval("fin rail pads (yawed, to the equator face)", *fins["rail"]["pad_x_mm"], "internal"),
        _interval("motor carrier", g["x_fwd_mm"] - CARRIER_T, g["x_fwd_mm"], "internal", r=g["r_carrier_mm"]),
        _interval("generator bay D%g x L%g" % (g["D_mm"], g["L_mm"]), g["x_fwd_mm"], g["x_aft_mm"], "internal",
                  r=g["D_mm"] / 2),
        _interval("%s body" % mo["name"], *g["motor_body_x_mm"], "internal", r=mo["D"] / 2),
        _interval("blade pockets (static hull, 3x)", ri["x0"], x_end, "rotor",
                  note="pocket starts with a 45 deg ramp at x %.2f (printable floor)" % ri["x0"]),
        _interval("claw A (motor side)", *st["claw_A"], "internal", r=COUPLER_OD / 2),
        _interval("claw B + thrust collar", *st["claw_B_thrust_collar"], "internal", r=COUPLER_OD / 2),
        _interval("bayonet spigot (aft shell)", J["x_tip"], J["x_eq"], "internal", r=J["r_spigot_seal"]),
        _interval("socket ring (front shell)", J["x_ring0"], J["x_eq"], "internal", r=J["r_bore"]),
        _interval("equator joint plane", J["x_eq"], J["x_eq"], "internal"),
        _interval("bearing tube (aft shell)", DT["tube"]["x_from"], DT["tube"]["x_to"], "internal",
                  r=DT["tube"]["r_out_front_mm"]),
        _interval("608-2RS (locating)", *st["bearing_front_608"], "internal", r=P.BEARING["D"] / 2),
        _interval("688-2RS (floating)", *st["bearing_rear_688"], "internal", r=P.BEARING_ALT["D"] / 2),
        _interval("seal shoulder", *st["seal_shoulder"], "internal"),
        _interval("lip seal TC 8x16x5 (v2)", *st["lip_seal"], "internal", r=P.LIP_SEAL["D"] / 2),
        _interval("tail end block (seal boss, notch floors)", st["bearing_rear_688"][1], x_end, "internal"),
        _interval("stowed-detent magnet seats (3x, notch floor)", ri["rotor"]["stow_magnet_x"] - P.MAGNET["D"] / 2,
                  ri["rotor"]["stow_magnet_x"] + P.MAGNET["D"] / 2, "rotor", r=ri["rotor"]["stow_magnet_r"]),
        _interval("static tail end face", x_end, x_end, "both"),
        _interval("hub gap", x_end, ch["hub_front_x_mm"], "both"),
        _interval("hub (rotating)", ch["hub_front_x_mm"], hub["hub_aft_face_x_mm"], "rotor"),
        _interval("stop block (hub front ring)", *hub["stop"]["stop_block_x_mm"], "rotor"),
        _interval("M3 grub on the shaft D-flat", hub["grub"]["x_mm"], hub["grub"]["x_mm"], "rotor"),
        _interval("hinge pins (tangential, r %g)" % ch["hinge_r_mm"], ch["hinge_x_mm"] - ch["knuckle_r_mm"],
                  ch["hinge_x_mm"] + ch["knuckle_r_mm"], "rotor", r=ch["hinge_r_mm"],
                  note="pin axis at x %g" % ch["hinge_x_mm"]),
        _interval("heel, deployed", *rj["interference"]["heel_deployed_x_range_mm"], "rotor"),
        _interval("shaft (8 mm)", *st["shaft"], "internal", r=P.SHAFT_D / 2),
        _interval("spinner shell", ch["hub_front_x_mm"], P.HULL_L, "rotor"),
        _interval("rotor thrust centroid (0.7 R, ESTIMATE)", st["rotor_load_plane"], st["rotor_load_plane"], "both"),
    ]
    stations.sort(key=lambda d: (d["x0_mm"], d["x1_mm"]))

    # ---- interface checks between the halves
    xs = np.linspace(g["x_fwd_mm"], g["x_aft_mm"], 200)
    bay_floor_gap = float((pk.floor_inner_r(xs) - g["D_mm"] / 2).min())
    xr = np.linspace(DT["tube"]["x_from"], DT["tube"]["x_to"], 300)
    mag_gap = ri["rotor"]["stow_magnet_keepout_r"] - P.LIP_SEAL["D"] / 2
    fold_c = fins["clearance"].get("fold_at_index") or {}
    fold_min = min(v for v in (fold_c.get("pad_clear_min_mm"), fold_c.get("fin_clear_min_mm"), 1e9) if v is not None)
    checks = [
        dict(check="static tail end: rotor half = internal half", value_mm=x_end, limit_mm=ch["static_tail_x_mm"],
             passed=abs(x_end - ch["static_tail_x_mm"]) < 1e-9),
        dict(check="hinge (x, r): rotor half = internal half", value_mm=[ri["hinge_x"], ri["hinge_r"]],
             limit_mm=[ch["hinge_x_mm"], ch["hinge_r_mm"]],
             passed=ri["hinge_x"] == ch["hinge_x_mm"] and ri["hinge_r"] == ch["hinge_r_mm"]),
        dict(check="shaft aft end = hub aft face (spinner M3 screw thread)", value_mm=st["shaft"][1],
             limit_mm=hub["hub_aft_face_x_mm"], passed=abs(st["shaft"][1] - hub["hub_aft_face_x_mm"]) < 1e-9),
        dict(check="hub gap = HUB_GAP", value_mm=ch["hub_front_x_mm"] - x_end, limit_mm=P.HUB_GAP,
             passed=abs(ch["hub_front_x_mm"] - x_end - P.HUB_GAP) < 1e-9),
        dict(check="lip seal outer face ahead of the hub front (axial gap)", value_mm=ch["hub_front_x_mm"] - st["lip_seal"][1],
             limit_mm=P.HUB_GAP, passed=ch["hub_front_x_mm"] - st["lip_seal"][1] >= P.HUB_GAP),
        dict(check="no axial overlap: 688 / shoulder / seal inside the static tail", value_mm=st["lip_seal"][1],
             limit_mm=x_end, passed=st["bearing_rear_688"][1] <= st["seal_shoulder"][0] + 1e-9
             and st["seal_shoulder"][1] <= st["lip_seal"][0] + 1e-9 and st["lip_seal"][1] <= x_end),
        dict(check="M3 grub (shaft D-flat) inside the hub on the shaft", value_mm=hub["grub"]["x_mm"],
             limit_mm=[ch["hub_front_x_mm"], st["shaft"][1]],
             passed=ch["hub_front_x_mm"] < hub["grub"]["x_mm"] < st["shaft"][1]),
        dict(check="bearing tube to pocket floor (internal 2.0 mm floor wall)", value_mm=DT["tube"]["min_gap_to_pocket_floor_mm"],
             limit_mm=0.0, passed=DT["tube"]["min_gap_to_pocket_floor_mm"] >= 0,
             note="at x %.1f (tail notch)" % DT["tube"]["min_gap_to_pocket_floor_x_mm"]),
        dict(check="bearing tube to rotor keep-out (floor - MIN_WALL)", value_mm=DT["tube"]["min_gap_to_rotor_keepout_mm"],
             limit_mm=0.0, passed=DT["tube"]["min_gap_to_rotor_keepout_mm"] >= 0),
        dict(check="688 tube radius under the tail notch keep-out", value_mm=DT["tube"]["r_out_rear_mm"],
             limit_mm=ri["rotor"]["keepout_tail_min_r"], passed=DT["tube"]["r_out_rear_mm"] < ri["rotor"]["keepout_tail_min_r"]),
        dict(check="stowed-detent magnet keep-out to the seal bore (end block material)", value_mm=mag_gap,
             limit_mm=T_MIN, passed=mag_gap >= T_MIN),
        dict(check="socket wall under the pocket floor at the joint", value_mm=J["walls"]["socket_at_pocket_mm"],
             limit_mm=POCKET_FLOOR_WALL, passed=J["walls"]["socket_at_pocket_mm"] >= POCKET_FLOOR_WALL - 1e-9),
        dict(check="bayonet twist channel to pocket edges (arc)", value_mm=J["bayonet"]["channel_to_pocket_edge_mm"],
             limit_mm=FIN_POCKET_CLEAR, passed=(J["bayonet"]["channel_to_pocket_edge_mm"] or 0) >= FIN_POCKET_CLEAR),
        dict(check="generator bay to the shell inner wall and pocket floors (radial)", value_mm=bay_floor_gap, limit_mm=C_GEN,
             passed=bay_floor_gap >= C_GEN),
        dict(check="yawed fin rail pads to pocket edges (arc)", value_mm=fins["clearance"]["stowed_pockets"]["min_clear_mm"],
             limit_mm=FIN_POCKET_CLEAR, passed=bool(fins["clearance"]["stowed_pockets"]["ok"]),
             note="pad yawed with the fin about the pivot at %.2f c, run to the equator face"
                  % fins["yaw"]["pivot_chord_frac"]),
        dict(check="yawed fins and pads to blades folding at the stowing index (3D)",
             value_mm=fold_min if fold_min < 1e9 else None,
             limit_mm=FIN_FOLD_CLEAR, passed=bool(fold_c.get("ok", True)),
             note="pad %.2f, fin %.2f" % (fold_c.get("pad_clear_min_mm") or float("nan"),
                                         fold_c.get("fin_clear_min_mm") or float("nan"))),
        dict(check="fin pad aft end printable (runs to the equator face or 45 deg chamfer)", value_mm=fins["rail"]["pad_x_mm"][1],
             limit_mm=J["x_eq"], passed=bool(fins["rail"]["pad_x_mm"][1] <= J["x_eq"] + 1e-9)),
        dict(check="pocket floor overhang (printed +x down, rotor half)", value_mm=None,
             limit_mm=None, passed=bool((ri.get("floor_overhang_max_deg") or 90.0) <= P.MAX_OVERHANG + 1e-6),
             note="max %.1f deg from vertical vs MAX_OVERHANG %.0f" % (ri.get("floor_overhang_max_deg") or float("nan"),
                                                                       P.MAX_OVERHANG)),
        dict(check="claw B (rotating) to the static tube mouth, radial", value_mm=DT["coupler"]["tube_mouth"]["min_radial_clear_mm"],
             limit_mm=P.CLEAR_ROTATE, passed=DT["coupler"]["tube_mouth"]["min_radial_clear_mm"] >= P.CLEAR_ROTATE),
        dict(check="bearing span rule (span >= 4 d, tip runout, s0, shaft deflection)", value_mm=br["spacing_mm"],
             limit_mm=SPAN_MIN, passed=bool(br["ok"]),
             note="span/overhang %.2f vs %.2f rule of thumb (reported, not enforced)"
                  % (br["spacing_over_overhang"], SPAN_RATIO_REF)),
        dict(check="hub wall over the shaft bore at the heel slots (rotor)", value_mm=hub["heel"]["hub_wall_over_shaft_bore_mm"],
             limit_mm=P.MIN_WALL, passed=hub["heel"]["hub_wall_over_shaft_bore_mm"] >= P.MIN_WALL - 1e-9),
        dict(check="rotor interference (folding, neighbours, heel, deployed span)", value_mm=None, limit_mm=None,
             passed=bool(rj["interference"]["pass"]), note=rj["interference"]["criteria"]),
        dict(check="gland boss volume (seat + thread boss) ahead of the carrier plate", value_mm=nose_d["gland"]["boss_to_carrier_plate_mm"],
             limit_mm=GLAND_CLEAR, passed=nose_d["gland"]["boss_to_carrier_plate_mm"] >= GLAND_CLEAR,
             note="the carrier only moves at x >= %.1f (push in from aft, 30 deg twist in place), so this also "
                  "clears its twist path; the bay starts behind it" % nose_d["gland"]["carrier_front_x_mm"]),
        dict(check="gland boss volume ahead of the carrier ring on the shell", value_mm=nose_d["gland"]["boss_to_carrier_ring_mm"],
             limit_mm=GLAND_CLEAR, passed=nose_d["gland"]["boss_to_carrier_ring_mm"] >= GLAND_CLEAR),
        dict(check="gland body to the tow eye and the shackle envelope", value_mm=min(nose_d["gland"]["body_to_eye_ring_mm"],
                                                                               nose_d["gland"]["body_to_shackle_envelope_mm"]),
             limit_mm=GLAND_TOW_CLEAR, passed=min(nose_d["gland"]["body_to_eye_ring_mm"],
                                                  nose_d["gland"]["body_to_shackle_envelope_mm"]) >= GLAND_TOW_CLEAR),
        dict(check="carrier cable notch clear of the motor screw counterbores and the centre bore",
             value_mm=min(nose_d["gland"]["cable_notch"]["clear_to_screw_counterbores_mm"],
                          nose_d["gland"]["cable_notch"]["web_to_centre_bore_mm"]),
             limit_mm=T_MIN, passed=bool(nose_d["gland"]["cable_notch"]["ok"])),
        dict(check="joint lug shear utilisation (interlayer proxy)", value_mm=None, limit_mm=0.5,
             passed=bool(res["joint_strength"]["ok"]), note="util %.2f" % max(res["joint_strength"]["util_lug"],
                                                                                res["joint_strength"]["util_lip"])),
    ]

    # ---- cross-study checks: roll torque basis vs the BEM generator load line; O-ring vs the hardware study
    fm15 = fins["margins"][FIN_DESIGN_V]
    tb = fins["torque_basis"]
    cross = [dict(
        check="fin couple at %.1f m/s vs the %s load-line torque (BEM)" % (FIN_DESIGN_V, tb["generator"]),
        fin_torque_Nm=fm15["torque_Nm"], load_line_torque_Nm=tb["Q_design_Nm"], tsr=tb["tsr"], rpm=tb["rpm"],
        ratio=fm15["ratio"], limit=[tb["target_ratio"], tb["target_ratio_max"]],
        passed=bool(tb["target_ratio"] - 1e-9 <= fm15["ratio"] <= tb["target_ratio_max"] + 1e-9),
        off_design={f"{v:.1f}_m_s": dict(ratio=fins["margins"][v]["ratio"], load_line_torque_Nm=fins["margins"][v]["Q_gen_Nm"])
                    for v in FIN_SPEEDS if v != FIN_DESIGN_V},
        note="fins sized on the load line with a %.0f %% margin (fins.torque_basis); the merge's Cp 0.35 / TSR 4 rule "
             "(%.3f N m) is retired" % (tb["margin_pct"], tb["Q_cp_rule_Nm"]))]
    o_ = J["oring"]
    hw_o = hw.get("oring") or {}
    if hw_o:
        dev = max(abs(o_["bore_d_mm"] - hw_o["bore_d"]), abs(o_["groove_bottom_d_mm"] - hw_o["groove_bottom_d"]),
                  abs(o_["spigot_seal_od_mm"] - hw_o["spigot_od"]), abs(o_["groove_top_w_mm"] - hw_o["groove_top_w"]))
        hp = o_.get("hardware_study_performance") or {}
        cross.append(dict(
            check="equator O-ring and groove: layout = hardware study (%s)" % hw_o["designation"],
            max_dimension_deviation_mm=dev, passed=bool(dev < 1e-6 and o_["bore_fits"]),
            layout=dict(bore_d_mm=o_["bore_d_mm"], groove_bottom_d_mm=o_["groove_bottom_d_mm"],
                        spigot_od_mm=o_["spigot_seal_od_mm"], groove_top_w_mm=o_["groove_top_w_mm"],
                        squeeze_worst_pct=o_["squeeze_worst_pct"], fill_worst_pct=o_["fill_worst_pct"]),
            hardware=dict(bore_d_mm=hw_o["bore_d"], groove_bottom_d_mm=hw_o["groove_bottom_d"],
                          spigot_od_mm=hw_o["spigot_od"], groove_top_w_mm=hw_o["groove_top_w"],
                          squeeze_worst_pct=hp.get("squeeze_worst_pct"), fill_worst_pct=hp.get("fill_worst_pct")),
            note="the layout adopts the hardware study's ring and groove and recomputes the gland with the same "
                 "worst-case model"))
    issues_cross = []

    # ---- consolidated params.py proposal
    rp = rj["proposed_params"]
    ip = summ["proposed_param_updates"]
    params = {k: rp[k] for k in rp}
    det = rj["detents"]
    params.update(
        ROTOR_HUB=dict(x_front=ch["hub_front_x_mm"], x_aft_face=hub["hub_aft_face_x_mm"], bore_d=hub["bore_d_mm"],
                       d_flat_depth=hub["d_flat_depth_mm"], grub=hub["grub"]["size"], grub_x=hub["grub"]["x_mm"],
                       grub_psi_deg=hub["grub"]["azimuth_psi_deg"], stop_face_x=det["deployed"]["stop_face_x_mm"],
                       ear_hole_d=hub["pin"]["ear_hole_d_mm"], knuckle_bore_d=hub["pin"]["knuckle_bore_d_mm"],
                       pin_channel_tau=hub["pin"]["insertion_channel_tau_mm"],
                       spinner_wall=ch["spinner_wall_mm"], spinner_gap=ch["spinner_gap_mm"],
                       slot_table="data/layout_rotor.json hub_slot (cut per blade, hinge-plane coords)"),
        DETENTS=dict(magnet=det["magnet"],
                     deployed_hub_magnet_x_re_tau=det["deployed"]["hub_magnet_hinge_plane_x_re_tau_mm"],
                     deployed_blade_magnet_blade_frame=det["deployed"]["blade_magnet_blade_frame_mm"],
                     stowed_x=det["stowed"]["x_mm"], stowed_az_deg=det["stowed"]["az_deg"],
                     stowed_static_magnet_r=det["stowed"]["static_magnet_centre_r_mm"],
                     stowed_boss_top_r=det["stowed"]["boss_top_r_mm"],
                     stowed_blade_magnet_blade_frame=det["stowed"]["blade_magnet_centre_blade_frame_mm"]),
        POCKETS=dict(x=[ri["x0"], x_end], theta_deg=list(pk.theta),
                     tables="data/layout_rotor.json pocket_depth (x, depth, floor_r, keepout_r, az_lo, az_hi per x; "
                            "theta = -psi)", floor_print_deg=ri.get("floor_overhang_max_deg")),
        BLADE_AXIS=dict(table="data/layout_rotor.json axis (solved curved pitch axis)",
                        lead_mm=[rj["axis"]["lead_mm"][0], rj["axis"]["lead_mm"][-1]]),
    )
    params.update({k: v for k, v in ip.items() if not k.endswith("_reason")})
    params["BEARING"] = params.pop("BEARING_FRONT")
    params["BEARING_ALT"] = params.pop("BEARING_REAR")
    reasons = {f"{k}_reason": v for k, v in ROTOR_PARAM_REASONS.items() if k in rp}
    reasons.update({k: v for k, v in ip.items() if k.endswith("_reason")})
    reasons["BEARING_ALT_reason"] = reasons.get("BEARING_reason", "")
    reasons["ROTOR_HUB_reason"] = "rotor half (layout.py): hub, grub, stop face, pin holes and spinner shell"
    reasons["DETENTS_reason"] = ("rotor half: N52 4x2 magnet pairs, deployed (heel <-> hub stop face) and stowed "
                                 "(blade root <-> boss on the static notch floor)")
    reasons["POCKETS_reason"] = ("rotor half: pocket tables, printable floor (<= %.0f deg) and asymmetric azimuth span"
                                 % P.MAX_OVERHANG)
    reasons["BLADE_AXIS_reason"] = "rotor half: solved hull-following pitch axis the blade CAD sweeps along"
    reasons["LIP_SEAL_reason"] = ("TC 8x16x5 in the 688's 16 mm bore, %.1f mm recessed from the static tail end, "
                                  "%.1f mm ahead of the hub front" % (SEAL_RECESS, ch["hub_front_x_mm"] - st["lip_seal"][1]))
    reasons["SHAFT_reason"] = ("8 x %.1f from claw B (x %.1f) to the hub aft face (x %.1f): hub D-flat under the "
                               "M3 grub at x %.1f, M3 x 8 tapped aft end for the spinner screw"
                               % (st["shaft"][1] - st["shaft"][0], st["shaft"][0], st["shaft"][1], hub["grub"]["x_mm"]))
    reasons["GLAND_reason"] = ("M12x1.5 beside the eye, bottom, seat x %.2f, axis %.0f deg from the surface normal "
                               "toward the nose: the whole boss (D%g seat + D%g thread boss) ends %.2f mm ahead of "
                               "the carrier (the boss at x 12 ran into it)"
                               % (nose_d["gland"]["seat_x_mm"], nose_d["gland"]["axis_tilt_from_normal_deg"],
                                  GLAND_SEAT_D, GLAND_BOSS_D, min(nose_d["gland"]["boss_to_carrier_plate_mm"],
                                                                  nose_d["gland"]["boss_to_carrier_ring_mm"])))
    reasons["TOW_EYE_reason"] = ("legs at %.2f of PETG interlayer adhesion under 3 x the tow load; ring in the XZ "
                                 "plane so the shackle (pin along Y) stays clear of the gland" % (
                                     nose_d["eye"]["materials"]["PETG"]["eye_leg_util"]))
    reasons["POCKET_FLOOR_WALL_reason"] = ("new: shell under every pocket floor used by the internal parts (2.0, "
                                           "stricter than the rotor keep-out at MIN_WALL)")

    # ---- open issues: rotor items still open + internal items
    bem_open = next((t for t in rj.get("open_issues", []) if t.startswith("BEM")), None)
    rotor_open = [
        (bem_open + " " if bem_open else "") + "The rotor load plane used by the bearing checks (x %.1f) is the "
        "solved axis at 0.7 R (ESTIMATE), not a BEM thrust centroid." % st["rotor_load_plane"],
        "In-plane lead (pitch axis toward the LE, %.1f mm at the root to %.1f mm at the tip) yaws the sections "
        "by up to %.0f deg against the local tangential flow (rotor.stowed_sections.yaw_vs_local_tangent_deg); "
        "BEM assumes 0." % (rj["axis"]["lead_mm"][0], rj["axis"]["lead_mm"][-1],
                            max(rj["stowed_sections"]["yaw_vs_local_tangent_deg"])),
        "Magnet forces come from a model scaled to catalogue pull (k %.2f); the stowed hold (%.1f N mm) may stop "
        "the flow from opening the blades. Verify both detents on a test print."
        % (rj["detents"]["calibration"]["k_design"], rj["detents"]["stowed"]["torque_Nmm"]),
        "Pin 3 x 22 (ISO 8734): check stock; fallback 3 x 20 with 2.7 mm ears.",
        "Hub wall over the shaft bore is exactly MIN_WALL (%.2f mm) in the three heel slots (x 114-122)."
        % hub["heel"]["hub_wall_over_shaft_bore_mm"],
        "Rotor loads (stop %.0f N, pin shear %.1f MPa, ear bearing %.1f MPa at 2.5 m/s) are first-order "
        "ESTIMATES; printed PETG allowables were not sourced for them."
        % (hub["load_cases"]["2.5"]["stop_force_N"], hub["load_cases"]["2.5"]["pin_double_shear_MPa"],
           hub["load_cases"]["2.5"]["ear_bearing_MPa"]),
        "Stowed, the section tops are flush but the outer LE/TE edges sit 1.0-1.6 mm below the hull: a visible "
        "recess, accepted.",
    ]
    issues = summ["open_issues"][:1] + issues_cross + summ["open_issues"][1:] + rotor_open
    resolved = [
        "Internal layout now reads the real rotor tables (pocket depth / floor / keep-out / azimuth extent, "
        "blade axis, hub stations) through hydroegg.layout; the fallback pockets are gone.",
        "STATIC_TAIL_X 118 -> 108: bearings and seal moved 10 mm forward; the bay moved forward and shrank to "
        "L %g so the bearing span stays >= %g mm." % (g["L_mm"], SPAN_MIN),
        "Pockets crossing the equator: %.1f mm deep there (not 6), so the joint grew (%s, seal bore D%.2f); the "
        "socket keeps %.2f mm under the pocket floor and the lugs sit in the lands (%.1f mm channel clearance)."
        % (float(pk.depth(J["x_eq"])), J["oring"]["designation"].split()[-1], 2 * J["r_bore"],
           J["walls"]["socket_at_pocket_mm"], J["bayonet"]["channel_to_pocket_edge_mm"] or float("nan")),
        "Hub / bearing / seal axial overlap: none (seal ends x %.1f, hub starts x %.1f); shaft end = hub aft face "
        "x %.1f; the tail end block carries the stowed-detent magnet seats." % (st["lip_seal"][1],
                                                                                 ch["hub_front_x_mm"], st["shaft"][1]),
        "Fins: with the new hinge (116, 14) the any-azimuth fold torus limits the fin TE to x %.0f; fins now follow the stowing-"
        "index rule and sit at theta %g / %g." % (fins["clearance"]["torus_te_limit_x_mm"], *fins["theta_deg"]),
        "Review fix round: the fin rail pads are modelled yawed with the fins (pivot %.2f c) in the selection, "
        "the pocket-clearance and the fold-clearance checks, and run to the equator face; the gland moved to "
        "seat x %.2f with the full boss volume checked; buoyancy counts the hub/spinner slots, the spinner "
        "access bore and the pads; claw B steps to OD %g at the tube mouth; the bay is the largest envelope "
        "(D%g x L%g, carrier ring included) and candidate fits are tested against it; the pocket floors are "
        "printable (<= 45 deg); the end-view figures say 'seen from aft, looking upstream'."
        % (fins["yaw"]["pivot_chord_frac"], nose_d["gland"]["seat_x_mm"], CLAW_B_STEP_OD, g["D_mm"], g["L_mm"]),
        "Fix round completion: the equator joint adopts the hardware study's %s (bore D%.2f, groove bottom D%.2f, "
        "spigot OD %.2f, 45 deg tip-side flank, top width %.3f), recomputed with the corrected worst case "
        "(squeeze %.1f-%.1f %%, fill %.1f %%); the fins are re-sized on the %s load line (%.4f N m at TSR %.2f, "
        "1.5 m/s) with a %.0f %% margin: %g x %g mm at %.2f deg, %.0f mm tip to tip (the Cp 0.35 / TSR 4 rule gave "
        "%.3f N m, %.2f x the load-line torque)."
        % (J["oring"]["designation"], J["oring"]["bore_d_mm"], J["oring"]["groove_bottom_d_mm"],
           J["oring"]["spigot_seal_od_mm"], J["oring"]["groove_top_w_mm"], *J["oring"]["squeeze_worst_pct"],
           J["oring"]["fill_worst_pct"], fins["torque_basis"]["generator"], fins["torque_basis"]["Q_design_Nm"],
           fins["torque_basis"]["tsr"] or float("nan"), fins["torque_basis"]["margin_pct"], fins["chord_mm"],
           fins["span_mm"], fins["incidence_deg"], fins["tip_to_tip_mm"], fins["torque_basis"]["Q_cp_rule_Nm"],
           fins["torque_basis"]["Q_cp_rule_Nm"] / fins["torque_basis"]["Q_design_Nm"]),
    ]
    rot_keep = {k: rj[k] for k in ("chosen", "axis", "fold_choice", "stowed_sections", "stowed_max_over_hull",
                                   "pocket_depth", "pocket_summary", "pocket_width_mm", "pocket_theta_deg",
                                   "hub_slot", "hinge", "hub", "detents", "blade_mass", "interference") if k in rj}
    rot_keep["selection"] = {k: v for k, v in rj["selection"].items() if k != "hinge_sweep"}
    int_keep = {k: summ[k] for k in ("joint", "nose", "roll", "assembly", "part_split") if k in summ}
    gb = summ["generator_bay"]
    int_keep["generator_bay"] = {k: gb[k] for k in ("pick", "carrier", "length_trade", "candidate_fits",
                                                    "rotor_face_x_max_mm", "frontier_note") if k in gb}
    dtj = summ["drivetrain"]
    int_keep["drivetrain"] = {k: dtj[k] for k in ("stations_x_mm", "bearings", "tube", "seal", "shaft", "coupler")}
    int_keep["fins"] = {k: v for k, v in summ["fins"].items() if k not in ("sweep_feasible",)}
    mbj = summ["mass_buoyancy"]
    int_keep["mass_buoyancy"] = {k: mbj[k] for k in ("model", "display_PLA", "v2_PETG", "v2_needs", "v2_ballast",
                                                     "static_hull_outer_cc", "pockets_cc", "hardware")}
    return dict(
        schema="hydroegg.layout/1",
        generated_by="python -m hydroegg.layout (rotor half -> data/layout_rotor.json), then python -m "
                     "hydroegg.layout_internal (internal half -> data/layout_internal.json, and this merge)",
        sources=dict(rotor="data/layout_rotor.json (hydroegg/layout.py)",
                     internal="data/layout_internal.json (hydroegg/layout_internal.py)",
                     blade="data/blade_design.json", hardware="data/hardware.json"),
        conventions=dict(
            frame="X = body axis, flow +X, nose at x = 0, Z up, right handed; mm and deg, SI in physics",
            psi="rotor half: azimuth about +X from +Z, e_r = (0, -sin psi, cos psi); blade i at psi = 120 i; "
                "pocket az relative to the blade azimuth, + toward e_t (direction of rotation)",
            theta="internal half: azimuth about +X from +Z toward +Y; theta = -psi (blades at theta 0/240/120); "
                  "a pocket spans theta_blade + [-az_hi, -az_lo]",
            radial="pocket depth, floor and keep-out are radial (from the body axis); depth = deepest point over "
                   "azimuth at each x",
            estimates="any value marked ESTIMATE is a modelling assumption, not a sourced number"),
        summary=dict(
            static_tail_x_mm=x_end, hub_front_x_mm=ch["hub_front_x_mm"], hinge_x_mm=ch["hinge_x_mm"],
            hinge_r_mm=ch["hinge_r_mm"], fold_deg=ch["fold_deg"], shaft_end_x_mm=st["shaft"][1],
            equator_x_mm=J["x_eq"], bay_D_mm=g["D_mm"], bay_L_mm=g["L_mm"], bay_x_mm=[g["x_fwd_mm"], g["x_aft_mm"]],
            bearing_608_x_mm=st["bearing_front_608"], bearing_688_x_mm=st["bearing_rear_688"],
            lip_seal_x_mm=st["lip_seal"], bearing_span_mm=br["spacing_mm"], rotor_overhang_mm=br["overhang_mm"],
            span_over_overhang=br["spacing_over_overhang"], tip_runout_mm=br["tip_axial_runout_mm"],
            oring=J["oring"]["designation"], seal_bore_d_mm=2 * J["r_bore"], spigot_bore_d_mm=2 * J["r_spigot_in"],
            bayonet_lock_theta_deg=J["bayonet"]["lock_theta_deg"], bayonet_twist_deg=J["bayonet"]["engage_deg"],
            pocket_x_mm=[ri["x0"], x_end], pocket_depth_at_equator_mm=float(pk.depth(J["x_eq"])),
            pocket_depth_aero_max_mm=rj["pocket_summary"]["depth_aero_max_mm"],
            pocket_depth_max_mm=rj["pocket_summary"]["depth_max_mm"],
            fin_chord_mm=fins["chord_mm"], fin_span_mm=fins["span_mm"], fin_theta_deg=fins["theta_deg"],
            fin_x_mm=[fins["x_le_mm"], fins["x_te_mm"]], fin_incidence_deg=fins["incidence_deg"],
            fin_tip_to_tip_mm=fins["tip_to_tip_mm"],
            display_PLA_mass_g=summ["mass_buoyancy"]["display_PLA"]["mass_g"],
            display_PLA_net_buoyancy_N=summ["mass_buoyancy"]["display_PLA"]["net_buoyancy_N"],
            v2_mass_g=summ["mass_buoyancy"]["v2_PETG"]["mass_g"],
            v2_net_buoyancy_no_ballast_N=summ["mass_buoyancy"]["v2_needs"]["net_buoyancy_no_ballast_N"],
            tow_design_load_N=nose_d["loads"]["design_load_N"],
            interface_checks_passed=sum(bool(c["passed"]) for c in checks), interface_checks=len(checks)),
        stations=stations,
        interface_checks=checks,
        cross_study_checks=cross,
        merge_decisions=[
            dict(topic="who moved", decision="the internal half adapted to the rotor half everywhere; no rotor "
                 "value changed", why="task rule: resolve conflicts on the internal side first"),
            dict(topic="bearing span criterion", decision="replace the span/overhang >= 1.25 rule of thumb by "
                 "explicit checks: span >= %g mm (4 d), tip runout <= %g mm from %.2f mm seat slop, s0 >= %g, "
                 "shaft deflection <= %g mm" % (SPAN_MIN, TIP_RUNOUT_MAX, SEAT_SLOP, S0_MIN, SHAFT_DEFL_MAX),
                 why="with the static tail at x %.0f and the thrust centroid at x %.1f (solved curved axis) the "
                     "overhang is %.1f mm; even the GL40 alone at its forward-most position gives span/overhang "
                     "%.2f. Loads are small (s0 > 20), so the rule of thumb is not the binding physics"
                     % (x_end, st["rotor_load_plane"], br["overhang_mm"],
                        next((t["span_over_overhang"] for t in bay["length_trade"]
                              if abs(t["L_mm"] - mo["L"]) < 1e-6 and t.get("fits")), float("nan")))),
            dict(topic="generator bay", decision="D%g x L%g = the largest envelope that fits (task rule), "
                 "forward-most at x %.2f-%.2f; fits %s"
                 % (g["D_mm"], g["L_mm"], g["x_fwd_mm"], g["x_aft_mm"],
                    ", ".join(c_["name"] for c_ in bay["candidate_fits"] if c_["fits_bay"])),
                 why="D is set by the nose wall and the carrier ring at the forward end, L by the bearing span; "
                     "the L%g reserve would cut the span to %.1f mm" % (
                     mo["L_reserve"], next((t["span_mm"] for t in bay["length_trade"]
                                            if abs(t["L_mm"] - mo["L_reserve"]) < 1e-6 and t.get("fits")), float("nan")))),
            dict(topic="fins", decision="stowing-index fold rule, fins %g x %g at theta %g/%g, x %.1f-%.1f, pads "
                 "yawed with the fins (pivot %.2f c) and run to the equator face; sized on the %s load line at "
                 "%.1f m/s with a %.0f %% margin"
                 % (fins["chord_mm"], fins["span_mm"], *fins["theta_deg"], fins["x_le_mm"], fins["x_te_mm"],
                    fins["yaw"]["pivot_chord_frac"], fins["torque_basis"]["generator"], FIN_DESIGN_V,
                    fins["torque_basis"]["margin_pct"]),
                 why="the any-azimuth torus rule gives %s" % (
                     "%g x %g fins, %.0f mm tip to tip" % (fins["alternative_rule"]["chord_mm"],
                                                          fins["alternative_rule"]["span_mm"],
                                                          fins["alternative_rule"]["tip_to_tip_mm"])
                     if (fins.get("alternative_rule") or {}).get("chord_mm") else "no feasible fin")),
        ],
        decisions_internal=summ["decisions"],
        rotor=rot_keep,
        internal=int_keep,
        proposed_param_updates=dict(params, **reasons),
        resolved_conflicts=resolved,
        open_issues=issues,
    )


def main() -> None:
    res = run()
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "last_run.json").write_text(json.dumps(_clean(res), indent=1))
    out = _clean(summarize(res))
    DATA_JSON.parent.mkdir(parents=True, exist_ok=True)
    DATA_JSON.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", DATA_JSON)
    mg = _clean(merged(res, out))
    MERGED_JSON.write_text(json.dumps(mg, indent=1) + "\n")
    print("wrote", MERGED_JSON, "(interface checks passed %d / %d)" % (
        mg["summary"]["interface_checks_passed"], mg["summary"]["interface_checks"]))
    for f in figures(res):
        print("wrote", f)


if __name__ == "__main__":
    main()
