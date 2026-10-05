"""Single source of truth for Ovo. Every module reads from here; nothing
hard-codes a dimension that lives in this file.

Units: mm and degrees for geometry, SI for physics.
Frame: X = body axis, flow travels +X (nose at x = 0, rotor at the tail).
       Z = up, Y = lateral, right-handed. The rotor spins about the X axis.
Azimuth theta is measured about +X from +Z toward +Y unless a block says otherwise.

History: values were seeded on 2026-10-03 and promoted on 2026-10-04 from the verified
design studies (docs/notes/blade-design.md, hardware-and-prior-art.md, layout.md). Tables
too big for this file stay in data/*.json and are named in the blocks that use them:
  data/blade_design.json   blade stations (chord, twist, t/c, NACA), performance
  data/layout_rotor.json   solved curved pitch axis, pocket_depth table, hub_slot table
  data/layout.json         merged layout incl. interface checks
  data/hardware.json       bought parts with sources
VERIFY marks catalogue data to check against the part actually bought.
"""

# --------------------------------------------------------------------------
# Hull: Myring body of revolution, reads as an elongated egg (L/D ~ 2.2)
# --------------------------------------------------------------------------
HULL_D = 75.0             # max diameter
HULL_NOSE_A = 55.0        # nose length (super-ellipse)
HULL_NOSE_N = 2.0         # nose exponent: 2 = elliptic, the "wide end of an egg"
HULL_MID_B = 5.0          # short parallel band at the equator
HULL_TAIL_C = 105.0       # Myring cubic tail, closes to a point at HULL_L
HULL_TAIL_THETA = 25.0    # tail half-angle at the tip (deg); 25 keeps the tail full
HULL_L = HULL_NOSE_A + HULL_MID_B + HULL_TAIL_C   # 165

HULL_WALL = 2.4           # shell wall: 6 perimeters at 0.4 mm nozzle

# The egg splits at its equator (middle of the parallel band), Kinder-egg style.
EQUATOR_X = HULL_NOSE_A + HULL_MID_B / 2          # 57.5

# The static nacelle ends here; aft of it the hub + spinner (rotating) continue the
# same profile to the tip. Chosen by the layout hinge sweep (root clearance >= 0.8 mm).
STATIC_TAIL_X = 108.0
HUB_GAP = 0.6             # axial running gap static tail <-> hub (hub front at 108.6)
SPINNER_TIP_R = 3.0       # blunt the Myring point so it prints and survives knocks

# --------------------------------------------------------------------------
# Rotor: 3 folding blades at the tail, downstream rotor
# --------------------------------------------------------------------------
N_BLADES = 3
# Blades hinge about a TANGENTIAL axis on the hub. Deployed = radial with an aft
# pre-cone (12.6 deg mean, from the solved curved axis); stowed = folded FORWARD over
# the static tail into shallow pockets. Thrust opens a folded blade and holds it on
# the deployed stop.
HINGE_X = 116.0
HINGE_R = 14.0
FOLD_DEG = 90.0
ROTOR_R_TIP_TARGET = 97.0 # deployed tip radius (rotor D 194)

PIN_D = 3.0               # ISO 8734 dowel, m6
PIN_L = 22.0              # ear + knuckle + ear; 3 x 22 is a standard length, VERIFY stock

KNUCKLE_D = 8.0           # blade root knuckle around the pin
KNUCKLE_L = 14.0          # knuckle length along the pin
HUB_EAR_T = 3.7           # hub ear thickness along the pin
SHANK_CONE_DEG = 55.0     # shank cone at the knuckle
HEEL_REACH = 8.0          # heel (deployed stop + magnet) reach from the pin
HEEL_T = 5.0
HEEL_W = 10.0
POCKET_CLEAR = 0.3        # pocket floor clearance under the stowed aero span
ROOT_CLEAR = 0.8          # stowed root to static tail clearance
POCKET_FLOOR_WALL = 2.0   # shell under every pocket floor where internal parts sit

# Blade design (BEM study owns the stations: data/blade_design.json)
DESIGN_FLOW_V = 1.5       # m/s, a decent river current; towed use is 2-2.5 m/s
DESIGN_TSR = 4.0          # best Cp of the printable blade (0.367 at TSR 3.8-4)
WATER_RHO = 1000.0        # kg/m^3
WATER_NU = 1.0e-6         # m^2/s at ~20 C
MIN_CHORD = 14.0          # printability floor at the tip
MIN_TE = 0.6              # trailing edge thickness floor (1.5 x nozzle line)
MAX_ROOT_TWIST = 25.0     # folding constraint (not binding: final root twist 17.2)
BLADE_R_AERO_ROOT = 30.0  # aerodynamic root radius; shank from the knuckle to here
BLADE_MAX_CHORD = 28.0    # three stowed blades side by side around the tail
BLADE_PITCH_AXIS = 0.30   # stacking axis as a chord fraction

# v2 generator load line (GL40, flow m/s -> TSR) with the cavitation cap from 2.2 m/s,
# and the flow-independent cap as (rpm, minimum shaft torque N m). Controller data.
GEN_TSR_SCHEDULE = [(1.0, 5.00), (1.5, 5.75), (2.0, 6.50), (2.2, 6.56), (2.3, 6.37), (2.4, 6.19), (2.5, 6.03)]
CAV_MIN_SHAFT_TORQUE = [(1295, 0.031), (1335, 0.106), (1377, 0.194), (1420, 0.296), (1463, 0.413), (1485, 0.477)]

# --------------------------------------------------------------------------
# Bought parts (data/hardware.json has sources and confidence)
# --------------------------------------------------------------------------
SHAFT_D = 8.0             # MISUMI PSSFJ8 class: SUS440C 56 HRC g6 (lip seal needs >= 45 HRC)
MAGNET = dict(D=4.0, H=2.0, grade="N52")   # 4.05-4.58 N pull to steel; 16 per build
INSERT_M3 = dict(L=5.7, hole_d=4.0, min_wall=1.6)   # CNC Kitchen M3 standard
SHAFT_STOP = dict(type="set-screw collar", d=8.0, OD=16.0, W=8.0,
                  alt="DIN 471 A8 circlip")         # forward of the 608: thrust pulls aft
BEARING = {'name': '608-2RS', 'd': 8, 'D': 22, 'B': 7, 'x0': 63, 'role': 'front, locating'}
BEARING_ALT = {'name': '688-2RS', 'd': 8, 'D': 16, 'B': 5, 'x0': 96, 'role': 'rear, floating'}
LIP_SEAL = {'name': 'TC 8x16x5', 'd': 8, 'D': 16, 'B': 5, 'x0': 102.5, 'recess': 0.5}
GENERATOR = {'name': 'CubeMars GL40 KV70',
 'D': 46.5,
 'L': 21.5,
 'mount_pcd': 24,
 'mount_screw': 'M3',
 'n_mount': 4,
 'rotor_face_pcd': 24,
 'bay_D': 50.5,
 'bay_L': 24,
 'x_bay_fwd': 22.5,
 'x_rotor_face': 46.5,
 'motor_body_x': [25, 46.5],
 'adapter_spacer_L': 2.5,
 'alternative': 'iPower GM4108H-120T',
 'alt_mount_pcd': 38.14}
ORING_EQUATOR = {'name': 'ISO 3601 / AS568 -141',
 'cs': 2.62,
 'id': 58.42,
 'groove_d': 59.88,
 'bore_d': 64.09,
 'spigot_od': 63.89,
 'groove_bottom_w': 2.9,
 'groove_top_w': 4.905,
 'groove_depth': 2.005,
 'gland_depth': 2.105,
 'flank': 'aft flank radial, tip-side flank 45 deg',
 'radial_clearance': 0.1}

# --------------------------------------------------------------------------
# Layout interfaces (layout study, 29/29 interface checks). x in mm along the axis.
# --------------------------------------------------------------------------
SHAFT = {'d': 8,
 'x0': 52.5,
 'x1': 128,
 'L': 75.5,
 'hub_flat_x': 111.6,
 'hub_flat_depth': 0.5,
 'aft_tap': 'M3 x 8'}

BEARING_TUBE = {'x0': 62,
 'x1': 108,
 'r_out_front': 13.4,
 'r_out_rear': 10.4,
 'profile_x': [62, 72, 94, 96, 108],
 'profile_r': [13.4, 13.4, 11.4, 10.4, 10.4],
 'seat_front_d': 22,
 'shoulder_front_id': 18,
 'shoulder_front_x': [70, 96],
 'seat_rear_d': 16,
 'seal_shoulder_x': [101, 102.5],
 'seal_shoulder_id': 12.5,
 'seal_bore_d': 16,
 'mouth_lead_in': 0.5,
 'ribs': {'n': 3, 't': 2, 'theta_deg': [0, 120, 240], 'x': [62, 108]}}

TAIL_END_BLOCK = {'x0': 101, 'x1': 108}

COUPLER = {'od': 22,
 'claw_a_L': 5,
 'float': 1,
 'claw_b_L': 10.5,
 'claw_a_x': [46.5, 51.5],
 'claw_b_x': [52.5, 63],
 'claw_b_step_od': 20,
 'claw_b_step_L': 2,
 'shim_od': 12,
 'shim_t': 0.5,
 'thrust_collar': 'Ruland MSC-8-SS 8x16x8'}

CARRIER = {'x0': 18.5,
 't': 4,
 'r': 23.3972,
 'twist_deg': 30,
 'centre_bore_d': 16,
 'screw_pattern_clock_deg': 45,
 'ring': {'x': [17.3, 23.7], 't_radial': 1.6, 'r_in_at_plate': 24.8669},
 'cable_notch': {'theta_deg': 180, 'width': 7.5, 'depth': 7.5, 'bottom_r': 15.8972}}

JOINT = {'x_eq': 57.5,
 'x_spigot_tip': 43.595,
 'x_lip': [46.595, 49.595],
 'x_seal_zone': [49.595, 57.5],
 'x_channel': [43.195, 46.595],
 'x_ring_ceiling_start': 41.695,
 'groove_x': {'x_aft': 56, 'x_fwd_bottom': 53.1, 'x_fwd_top': 51.095},
 'seal_bore_d': 64.09,
 'spigot_seal_od': 63.89,
 'groove_bottom_d': 59.88,
 'lug_od': 63.59,
 'spigot_lug_zone_od': 60.59,
 'lip_d': 61.09,
 'channel_d': 64.09,
 'spigot_bore_d': 55.88,
 'root_ring_L': 7.16,
 'land_root': 1.5,
 'land_tip': 1.5}

BAYONET = {'n_lugs': 3,
 'lug_w_deg': 16,
 'lug_h': 1.5,
 'lug_t': 3,
 'twist_deg': 35,
 'lock_theta_deg': [72, 192, 312],
 'entry_theta_deg': [37, 157, 277],
 'gap_w_deg': 20,
 'lip_t': 3,
 'channel_extra': 0.4}

ROTOR_HUB = {'x_front': 108.6,
 'x_aft_face': 128,
 'bore_d': 8.15,
 'd_flat_depth': 0.5,
 'grub': 'M3 x 4 cup point (ISO 4029)',
 'grub_x': 111.6,
 'grub_psi_deg': 60,
 'stop_face_x': 113.5,
 'ear_hole_d': 3,
 'knuckle_bore_d': 3.3,
 'pin_channel_tau': [11, 16.941],
 'spinner_wall': 1.2,
 'spinner_gap': 0.2,
 'slot_table': 'data/layout_rotor.json hub_slot'}

DETENTS = {'magnet': '4 x 2 mm N52 disc, axial',
 'deployed_hub_magnet_x_re_tau': [112.3, 8.4388, 0],
 'deployed_blade_magnet_blade_frame': [-5.5612, 0, -1.3],
 'stowed_x': 104.25,
 'stowed_az_deg': 18.75,
 'stowed_static_magnet_r': 19.996,
 'stowed_boss_top_r': 21.196,
 'stowed_blade_magnet_blade_frame': [11.75, 7.1989, 7.2073]}

POCKETS = {'x': [30.75, 108],
 'theta_deg': [0, 120, 240],
 'tables': 'data/layout_rotor.json pocket_depth (x, depth, floor_r, keepout_r, az_lo, az_hi; theta '
           '= -psi)',
 'floor_print_deg': 45}

BLADE_AXIS = {'table': 'data/layout_rotor.json axis (solved curved pitch axis)', 'lead_mm': [11.543, 4.752]}

FIN = {'section': 'clarky',
 'chord': 40,
 'span': 65,
 'incidence_deg': 12,
 'x_le': 16.3741,
 'x_te': 55.5,
 'theta_deg': [93, 273],
 'yaw_pivot_chord_frac': 0.15,
 'yaw_pivot_x': 22.243,
 'te_side': '-theta',
 'root_r_le': 33.7926,
 'root_r_te': 44.5966,
 'rail_w': 10,
 'rail_h_min': 4,
 'groove_w_mean': 6,
 'groove_depth': 4,
 'groove_form': 'half dovetail, 45 deg flank',
 'pad_x': [15.3345, 57.5],
 'end_stop_x': 55.5,
 'retention': 'M3 grub in a heat-set insert',
 'trim_variant_deg': 1}

GLAND = {'thread': 'M12x1.5',
 'x': 8.75,
 'theta_deg': 180,
 'tilt_from_normal_deg': 10,
 'axis_angle_from_x_deg': 146.6534,
 'seat_lift': 3,
 'spot_face_centre_x': 6.2439,
 'spot_face_centre_r': 21.9432,
 'seat_d': 18,
 'seat_raise': 6.7887,
 'thread_length': 8,
 'inner_end_x': 12.9268,
 'inner_end_r': 17.5456,
 'thread_engagement_min': 6.05,
 'boss_in_cavity': False}

TOW_EYE = {'hole_d': 8, 'ring_w': 5, 'ring_t': 6, 'plane': 'XZ', 'hole_centre_x': -5}

NOSE_PLUG = {'x_45': 11.0348, 'bridge_d': 38.5825}

SPINNER_ACCESS = {'d': 6, 'x': [134, 165], 'aft_wall': 3}

BALLAST_V2 = {'lead_g': 0, 'x': 83, 'theta_deg': 180}

# --------------------------------------------------------------------------
# Printing (Bambu Lab, 256 mm cube)
# --------------------------------------------------------------------------
BED = (256.0, 256.0, 256.0)
NOZZLE = 0.4
LAYER = 0.2
MIN_WALL = 1.2            # 3 perimeters
CLEAR_SLIDE = 0.25        # radial clearance, sliding fit (bayonet, fins)
CLEAR_ROTATE = 0.30       # radial clearance, free rotation (blade on pin)
PRESS_FIT = -0.05         # PER SIDE (-0.10 diametral) for pins and bearing seats, with 3
                          # crush ribs; 608 PLA pocket reference 21.95. Calibrated by the
                          # tolerance coupon before the final print.
MAX_OVERHANG = 45.0       # deg from vertical without supports
