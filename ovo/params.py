"""Single source of truth for Ovo. Every module reads from here; nothing
hard-codes a dimension that lives in this file.

Units: mm and degrees for geometry, SI for physics.
Frame: X = body axis, flow travels +X (nose at x = 0, rotor at the tail).
       Z = up, Y = lateral, right-handed. The rotor spins about the X axis.
Azimuth theta is measured about +X from +Z toward +Y unless a block says otherwise.

History: values were seeded on 2026-10-03 and promoted on 2026-10-04 from the verified
design studies (docs/notes/blade-design.md, hardware-and-prior-art.md, layout.md). On
2026-10-05 the as-built CAD deviations (review fix round T01-T31, polish round, owner
decisions; docs/notes/cad.md sections 9 and 10) were written back here, so the CAD and
this file agree again. Tables too big for this file stay in data/*.json and are named in
the blocks that use them:
  data/blade_design.json   blade stations (chord, twist, t/c, NACA), performance
  data/layout_rotor.json   solved curved pitch axis, pocket_depth table, hub_slot table
  data/layout.json         merged layout incl. interface checks
  data/hardware.json       bought parts with sources
The layout JSONs also carry each study's own proposed_param_updates: those are the study
outputs (pre-CAD); where they differ from this file, this file is the as-built truth.
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
HULL_L = HULL_NOSE_A + HULL_MID_B + HULL_TAIL_C   # 165: profile length (Myring point). The built
                          # egg ends at SPINNER_CAP x_tip_built 160.99 (R SPINNER_TIP_R round), L/D 2.15

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
ROTOR_R_TIP_TARGET = 97.0 # design tip radius of the pitch axis (BEM R, blade stations; rotor D 194)
ROTOR_R_ENVELOPE = 97.31  # MEASURED max radius of the deployed blade solids (97.313): the tip
                          # section's LE/TE corners lie off the radial line. Reporting only; the
                          # blade and BEM keep ROTOR_R_TIP_TARGET.

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
# CAD rule (blade.py, polish round R1): the skin's U parameter (around the section) is
# centripetal, averaged over the aero rows only; it replaces the T14 face split.

# Blade print cradle (finding T13, blade.py): sacrificial support rails in the blade's print
# coordinates, printed with the blade (slicer supports off).
BLADE_CRADLE = {'gap': 0.2,          # rail top below the blade (the slicer's support top-Z distance)
 'min_w': 1.6,                       # rail width floor (4 lines of 0.4); else the blade's lowest band + 0.4
 'station_step': 0.25,               # mm between rail stations along the span (review P6: at 1.0 the window
                                     # minimum over +-1 mm put the top 0.47-0.6 under the steep shank ramp)
 'min_float': 0.4,                   # rails only where the blade floats >= gap + one 0.2 layer
 'qty': 3}                           # one per blade

# v2 generator load line (GL40, flow m/s -> TSR) with the cavitation cap from 2.2 m/s,
# and the flow-independent cap as (rpm, minimum shaft torque N m). Controller data.
GEN_TSR_SCHEDULE = [(1.0, 5.00), (1.5, 5.75), (2.0, 6.50), (2.2, 6.56), (2.3, 6.37), (2.4, 6.19), (2.5, 6.03)]
CAV_MIN_SHAFT_TORQUE = [(1295, 0.031), (1335, 0.106), (1377, 0.194), (1420, 0.296), (1463, 0.413), (1485, 0.477)]

# --------------------------------------------------------------------------
# Bought parts (data/hardware.json has sources and confidence)
# --------------------------------------------------------------------------
SHAFT_D = 8.0             # MISUMI PSSFJ8 class: SUS440C 56 HRC g6 (lip seal needs >= 45 HRC)
# Stowed-detent disc (blade shank seat + aft-shell boss): 4.05-4.58 N pull to steel. The deployed
# detent uses 6 x 2 x 2 blocks (DETENTS deployed_magnet), the fin rails 4 x 2 x 2 blocks
# (FIN_MAGNET); 18 magnets per build in all. The layout's two equator pairs are not built.
MAGNET = dict(D=4.0, H=2.0, grade="N52", qty=6)
INSERT_M3 = dict(L=5.7, hole_d=4.0, min_wall=1.6,      # CNC Kitchen M3 standard
                 stator="none (T21): the replica stator takes the 4 carrier M3 x 8 screws self-tapping",
                 stator_pilot_d=2.5)                   # M3 self-tap pilot in the replica stator
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

# Fin-rail retention (recessed rail, owner decision 2026-10-05): one block in the groove floor
# (front shell) and one in the tongue bottom of the fin or the blanking strip.
FIN_MAGNET = {'size': [4.0, 2.0, 2.0],
 'magnetised': 'through 2 (pole faces 4 x 2)',
 'grade': 'N52',
 'seat': [4.1, 2.1],                 # glued seat (printed holes close up ~0.1)
 'seat_depth': 2.05,
 'recess': 0.05,                     # magnet face below the floor / tongue bottom
 'u_from_stop': 1.75,                # seat aft wall ahead of the rail end stop (along u)
 'qty': 6,                           # 2 shell floor, 2 fins, 2 blanking strips
 'part': 'N52 block 4 x 2 x 2 magnetised through 2 (pole faces 4 x 2), e.g. supermagnete '
         'Q-04-02-02-N52 (VERIFY)'}

# --------------------------------------------------------------------------
# Layout interfaces (layout study, 30/30 interface checks) and their as-built CAD values.
# x in mm along the axis.
# --------------------------------------------------------------------------
SHAFT = {'d': 8,
 'x0': 52.5,
 'x1': 128,
 'L': 75.5,
 'hub_flat_x': [108.4, 128.0],       # hub D-flat from just behind the hub front (108.6) out through the
                                     # aft end, so the hub's D-bore slides on; the grub sits at ROTOR_HUB grub_x
 'hub_flat_depth': 0.5,
 'collar_flat_x': [55.0, 62.0],      # flat under the thrust-collar set screw (x 58.5)
 'collar_flat_depth': 0.5,
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

# Interleaved claw coupling as built (T31 / T20): claw A on the GL40 rotor face (flange + 2 dogs
# carrying the screw counterbores, dogs out to the flange OD for a 1.2 mm skin), claw B pressed
# over the thrust collar (2 dogs forward). x ranges include the dogs; float 1 mm both ways.
COUPLER = {'od': 22,                 # claw B body OD
 'claw_a_L': 5.8,
 'float': 1,
 'claw_b_L': 13.0,
 'claw_a_x': [46.5, 52.3],
 'claw_b_x': [49.5, 62.5],
 'claw_a_od': 32.4,                  # flange / dog lobes: counterbore D6 at r 12 + 1.2 wall
 'claw_a_flange_t': 2.0,
 'claw_a_dogs_to_od': True,
 'claw_b_floor_t': 1.2,              # claw B floor ahead of the collar (MIN_WALL)
 'claw_b_cone_x': [60.5, 61.5],      # OD 22 -> claw_b_step_od cone (45 deg, prints aft face down)
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
 'cable_notch': {'theta_deg': 180, 'width': 7.5, 'depth': 7.5, 'bottom_r': 15.8972},
 # as built (T31): the GL40 solder-pad keep-out window at theta 180 is merged with the cable notch,
 # so the opening reaches r 12.0 (window r 12.0-20.5, half width 6.0, R2 corners)
 'pad_window': 'data/hardware.json solder_pad_keepout cad_window, merged with cable_notch',
 'engraving': 'chevron only (0.9 mm stroke, 0.5 deep, twist direction); no text (T23)'}

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
 'land_tip': 1.5,
 # O-ring lead-in at the equator (T05): the stretched -141 stands to r 32.528, so the mouth opens
 # beyond it at every azimuth; cut after the pocket liners (wall to the pocket floors >= 1.595)
 'seal_lead': [0.85, 18.0],          # radial mm, deg
 'seal_lead_mouth_d': 65.79}

BAYONET = {'n_lugs': 3,
 'lug_w_deg': 16,
 'lug_h': 1.5,
 'lug_t': 3,
 'twist_deg': 35,
 'lock_theta_deg': [72, 192, 312],
 'entry_theta_deg': [37, 157, 277],
 'gap_w_deg': 20,
 'lip_t': 3,
 'channel_extra': 0.4,
 # Click (T17): 3 curved flexure tines in the spigot wall between the lugs, each with a bump that
 # rides the socket ring ID (0.25 interference) and drops into an axial dimple at the lock angle.
 'click': {'theta_lock_deg': [120.0, 240.0, 0.0],   # bump azimuth when locked
           'tine_arc_deg': 28.0,
           'end_offset_deg': 3.0,
           'x_gap': 46.0,                            # tine root slot station
           # review 2026-10-05 (C1 / P1): the 45 deg underside slot is 1.15 along x (0.81 normal; it was
           # 0.6 / 0.42), and the print solid carries a breakaway tie (tie_h along x, full wall) across
           # each free-end slit at the tine's bottom, so the tine's first layers bridge between two
           # anchors instead of cantilevering 14 mm; cut the tie through the slit after printing
           'gap_w': 1.15,
           'tie_h': 0.8,
           'slit_w': 0.6,
           'bump_h': 0.5,
           'bump_top_hw': 0.6,
           'flank_deg': 30.0,
           'bump_x': [44.2, 46.0],
           'dimple_clear': 0.1,
           'dimple_x': [41.0, 46.3]}}

ROTOR_HUB = {'x_front': 108.6,
 'x_aft_face': 128,
 'bore_d': 8.15,
 'd_flat_depth': 0.5,
 'grub': 'M3 x 4 cup point (ISO 4029)',
 'grub_x': 111.6,
 'grub_psi_deg': 60,
 'stop_face_x': 113.5,
 'ear_hole_d': 2.9,                  # printed press bore = PIN_D + 2 x PRESS_FIT (3.0 is the nominal pin)
 'knuckle_bore_d': 3.3,              # running bore in the blade = PIN_D + 2 x CLEAR_ROTATE
 'pin_channel_tau': [11, 16.876],    # channel exit on the hub core (layout re-run with spinner_wall 1.25)
 'pin_filler': 'ISO 8734 3 x 4, loose in each pin channel, trapped by the spinner shell',
 'pin_filler_l': 4.0,
 'spinner_wall': 1.25,               # 1.2 measured 1.167 min on the mesh (T01); 1.25 keeps MIN_WALL
 'spinner_gap': 0.2,
 'slot_table': 'data/layout_rotor.json hub_slot'}

# Spinner tip cap (T02): printed push-in cap that continues profile.radius from the spinner's
# joint plane to an R SPINNER_TIP_R round.
SPINNER_CAP = {'x_joint': 151.725,   # spinner / cap joint plane (end of the 45 deg countersink)
 'spigot_d': 5.7,                    # in the D6 access bore (0.15 radial)
 'spigot_l': 12.0,
 'rib_tip_d': 6.1,                   # 3 crush ribs (0.05 radial interference)
 'rib_w': 0.6,
 'n_ribs': 3,
 'shoulder_deg': 45.0,               # shoulder seated in the spinner's countersink
 'groove_chamfer': 0.8,              # fingernail V-groove: 45 deg chamfers both sides of the joint
 'x_tip_built': 160.99}              # measured end of the egg (L/D 2.15)

DETENTS = {'magnet': 'deployed pair: N52 block 6 x 2 x 2 (deployed_magnet); stowed pair: N52 disc 4 x 2, '
                     'axial (MAGNET)',
 # deployed pair: blade heel face <-> hub stop face (T03: a D4 disc did not fit the 4.5 mm stop zone)
 'deployed_magnet': {'size_tau_re_x': [6.0, 2.0, 2.0],   # 6 along the pin, 2 radial, magnetised along x
                     'seat_tau_re': [6.1, 2.1],          # glued seats
                     'seat_depth': 2.2,
                     'r_e': 8.2,                         # block centre, hinge-plane r_e
                     'recess': 0.2,                      # each face below its stop face
                     'gap': 0.4,
                     'grade': 'N52',
                     'qty': 6,                           # 3 hub + 3 heel
                     'source': 'e.g. supermagnete Q-06-02-02-N52 (VERIFY)'},
 'deployed_hub_magnet_x_re_tau': [112.3, 8.2, 0.0],
 'deployed_blade_magnet_blade_frame': [-5.8, 0.0, -1.3],
 # stowed pair: blade shank inner face <-> aft-shell boss on the root-notch ramp
 'stowed_x': 104.25,
 'stowed_az_deg': 18.75,
 # layout anchor on the radial DETENTS line: blade.py centres the seat where this line meets the
 # blade's inner face and runs it along the face normal (the seat stays here, owner call T15)
 'stowed_blade_magnet_blade_frame': [11.75, 7.1989, 7.2073],
 'stowed_blade_magnet_seated_blade_frame': [10.899, 7.1625, 7.0023],   # as built
 'stowed_blade_magnet_axis_blade_frame': [0.7734, -0.2306, -0.5905],   # face normal, out of the face
 'stowed_blade_below_mouth': 0.1,    # blade magnet face below its seat mouth
 'stowed_boss_axis': 'the blade face normal at the shared mouth point (lean 49.03 deg from radial, '
                     'owner call T15)',
 'stowed_boss_lean_max_deg': 49.5,
 'stowed_boss_r': 3.25,
 'stowed_face_clear': 0.1,           # boss top below the blade's inner face
 'stowed_boss_recess': 0.05,         # static magnet face below the boss top
 'stowed_boss_top': 'contoured to the blade inner face (cubic fit) where it comes toward the boss, '
                    'flat plane elsewhere',
 'stowed_chin': {'len': 6.5, 'half_width': 1.25, 'clear_to_blade': 0.25},   # 45 deg start on the ramp
 'stowed_rest_zone_r': 8.25,         # boss + chin lie within this of the boss face centre (len + hw + 0.5)
 'stowed_static_magnet_centre_world_theta0': [103.259, -6.5972, 19.539],
 'stowed_static_magnet_axis_world_theta0': [0.7549, -0.2349, 0.6124],
 'stowed_pair_face_gap': 0.413,      # design (layout surface); 0.420 measured on the CAD blade, 1.66 deg
 # the layout study's radial pair (data/layout_rotor.json detents.stowed), superseded by the
 # leaned boss above; blade.py still reports against it
 'stowed_static_magnet_r': 19.996,
 'stowed_boss_top_r': 21.196}

POCKETS = {'x': [30.75, 108],
 'theta_deg': [0, 120, 240],
 'tables': 'data/layout_rotor.json pocket_depth (x, depth, floor_r, keepout_r, az_lo, az_hi; theta '
           '= -psi)',
 'floor_print_deg': 45,
 # CAD (shells.py): the ruled loft twists exact 45 deg ramps past 45, so the floor is re-limited
 # to 44.5 deg; on the root-notch ramp pocket_depth floor_r is therefore a MINIMUM depth (the CAD
 # floor is cut up to 0.156 mm deeper there)
 'cad_floor_print_deg': 44.5}

BLADE_AXIS = {'table': 'data/layout_rotor.json axis (solved curved pitch axis)', 'lead_mm': [11.543, 4.752]}

# Roll fins (owner decision 2026-10-05): recessed half-dovetail rail sunk into a thickened inner
# boss of the front shell, flush blanking strip when the fins are off (the stowed egg stays at
# HULL_D). Sizing and rail from the fin study (layout_internal.py, data/layout_internal.json FIN).
# Fin frame: u along the yawed root chord from the LE, w across (+ toward +theta), a radial.
FIN = {'section': 'clarky',
 'chord': 30.0,
 'span': 70.0,
 'incidence_deg': 12.5,
 'x_le': 26.2111,
 'x_te': 55.5,
 'theta_deg': [87.5, 267.5],
 'yaw_pivot_chord_frac': 0.05,
 'yaw_pivot_x': 27.6756,
 'te_side': '-theta',
 'root_r_le': 32.0525,              # root on the skin
 'root_r_te': 37.6,
 'root_gap': 0.1,                   # fin root face above the skin, tongue seated
 'rail': 'recessed (groove in an inner boss of the front shell, flush blanking strip when the fins are off)',
 'rail_w': 5.2,                     # groove neck = slot mouth
 'groove_w_floor': 6.5,             # 45 deg undercut to this at the floor
 'groove_w0': -2.6,                 # flat wall (the fin's flat Clark Y side faces -theta)
 'groove_undercut_h': 1.3,
 'groove_form': 'half dovetail, 45 deg flank on +w',
 'rail_slope': -0.25,               # da/du: the fin slides aft and inward onto the stop
 'rail_stop_x': 39.4,
 'rail_floor_a_stop': 29.6,
 'boss_floor_t': 1.6,
 'boss_side_t': 1.6,
 'boss_magnet_t': 3.45,
 'boss_cone_x': 43.1,
 'pad_x': [27.0045, 40.127],        # the fin study's planar footprint envelope of the slot (interval checks);
                                     # the built slot mouth at the skin spans x 27.80 to 39.95 (review C9)
 'end_stop_x': 39.4,
 'clear_floor': 0.05,               # tongue / strip bottom over the floor (the magnet pair pulls it down)
 'retention': 'N52 4 x 2 x 2 block pair (tongue bottom / groove floor), 1.8 mm ahead of the stop',
 'trim_variant_deg': 1.0,
 # Strip outline round (2026-10-05, review V8): the slot mouth reads as a deliberate inlay. The groove's aft
 # corners (end stop to the walls, the undercut flank included) are rounded slot_corner_r in the floor plane;
 # the tongue and the strip follow (R - their inset at each height, tangent to the same end stop). The entry
 # relief is no longer a full-width box with a 1.3 mm step: its +w wall is one straight taper ('ramp') from the
 # full width to the neck, sized by shells.rail_frame so the undercut lip left beside it stays >= MIN_WALL +
 # relief_lip_margin; the strip's +w edge follows it. The strip's across-face nail groove is replaced by a
 # half-round fingernail notch at its aft end (against the end stop).
 'slot_corner_r': 2.0,
 'relief_shape': 'ramp',            # +w edge: one straight taper from the full width at the entry to the neck at the
                                    # aft corner ('s': a short S-curve near the entry, read as a boot in review)
 'relief_lip_margin': 0.1,          # covers the ruled loft's chord between levels (<= 0.07) and the smooth max
 'strip_nail_notch': {'r': 1.5, 'depth': 0.8}}   # half disc on the strip centreline at its aft end, depth
                                    # below the skin (the nail pulls on its curved forward wall)

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
 'boss_in_cavity': False,
 'boss_root_fillet_r': 6.0,          # review V2: the seat boss blends into the nose (one smooth bump); R7 fails in OCC
 'tap_drill_d': 10.5,                # M12x1.5 tap drill (bore through block, wall and nose plug)
 # Owner decision 2026-10-05 (nose round): the DISPLAY front shell (front_shell) has a clean egg nose, no seat,
 # no bore and no plug (the blanking plug part is removed from the build and the BOM). The seat block and the
 # tap bore above exist only on the v2 front shell (front_shell_v2, role "v2"); the cutaway front shell is cut
 # from it, so the cutaway and the viewer's v2 mode show the seat with the gland in it.
 'display': 'none: clean nose (front_shell); the seat, the bore and the M12 gland are v2 only (front_shell_v2, '
            'the cutaway and the viewer v2 mode)',
 'v2_shell': 'front_shell_v2'}

TOW_EYE = {'hole_d': 8, 'ring_w': 5, 'ring_t': 6, 'plane': 'XZ', 'hole_centre_x': -5,
 # T06: outer profile a teardrop parallel to the hole's (45 deg roof toward -X), so the crown keeps
 # ring_w normal to the flanks (5.0 normal / 7.45 along the pull); the line bears on two flanks
 'outer_profile': 'teardrop parallel to the hole',
 'outer_apex_r': 1.5,
 'hole_apex_flat': 2.0,              # hole apex cut to a flat (a 2 mm bridge printed nose-up)
 'root_fillet_r': 2.0}

NOSE_PLUG = {'x_45': 11.0348, 'bridge_d': 38.5825}

# D6 screw-access bore aft of the hub face: from the M3 head seat (x_aft_face + aft_wall) to a
# 45 deg countersink that ends on the SPINNER_CAP joint plane
SPINNER_ACCESS = {'d': 6, 'x': [131.0, 148.354], 'aft_wall': 3, 'countersink_x': [148.354, 151.725]}

BALLAST_V2 = {'lead_g': 0, 'x': 83, 'theta_deg': 180}

# --------------------------------------------------------------------------
# Printing (Bambu Lab, 256 mm cube)
# --------------------------------------------------------------------------
BED = (256.0, 256.0, 256.0)
NOZZLE = 0.4
LAYER = 0.2
MIN_WALL = 1.2            # 3 perimeters
CLEAR_SLIDE = 0.25        # radial clearance per side, sliding fit (bayonet, fins)
CLEAR_ROTATE = 0.15       # RADIAL per side, rotating fit (blade knuckle bore 3.3 on the 3.0 pin)
RUN_CLEAR = 0.30          # radial per side, swept / running clearance between moving bodies
                          # (blade root to hub and spinner, knuckle to ears, replica bell journal)
PRESS_FIT = -0.05         # PER SIDE (-0.10 diametral) for pins and bearing seats, with 3
                          # crush ribs; 608 PLA pocket reference 21.95. Calibrated by the
                          # tolerance coupon before the final print.
MAX_OVERHANG = 45.0       # deg from vertical without supports
