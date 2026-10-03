"""Single source of truth for HydroEgg. Every module reads from here; nothing
hard-codes a dimension that lives in this file.

Units: mm and degrees for geometry, SI for physics.
Frame: X = body axis, flow travels +X (nose at x = 0, rotor at the tail).
       Z = up, Y = lateral, right-handed. The rotor spins about the X axis.

Values marked PROVISIONAL are starting points that a design study is expected to
refine; values marked VERIFY come from typical catalogue data and must be checked
against the datasheet of the part actually bought.
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
# same profile to the tip. PROVISIONAL: layout study may move it.
STATIC_TAIL_X = 118.0
HUB_GAP = 0.6             # axial running gap static tail <-> hub
SPINNER_TIP_R = 3.0       # blunt the Myring point so it prints and survives knocks

# --------------------------------------------------------------------------
# Rotor: 3 folding blades at the tail, downstream rotor
# --------------------------------------------------------------------------
N_BLADES = 3
# Blades hinge about a TANGENTIAL axis on the hub. Deployed = radial (with aft
# pre-cone); stowed = folded FORWARD over the static tail into shallow pockets.
# Flow thrust pushes a folded blade open and holds it against the deployed stop.
HINGE_X = 124.0           # PROVISIONAL
HINGE_R = 18.5            # PROVISIONAL, distance of the pin axis from the body axis
ROTOR_R_TIP_TARGET = 97.0 # PROVISIONAL, deployed tip radius (rotor ~ D 194)

PIN_D = 3.0               # steel dowel / ground pin, VERIFY length on purchase
PIN_L = 16.0

# Blade design inputs (BEM study owns the final numbers in data/blade_design.json)
DESIGN_FLOW_V = 1.5       # m/s, a decent river current; towed use is 2-2.5 m/s
DESIGN_TSR = 4.0          # PROVISIONAL tip speed ratio
WATER_RHO = 1000.0        # kg/m^3
WATER_NU = 1.0e-6         # m^2/s at ~20 C
MIN_CHORD = 14.0          # printability floor at the tip
MIN_TE = 0.6              # trailing edge thickness floor (1.5 x nozzle line)
MAX_ROOT_TWIST = 25.0     # folding constraint: big root twist pokes out when stowed

# --------------------------------------------------------------------------
# Drivetrain (display build has no electronics; v2 drops a real motor in)
# --------------------------------------------------------------------------
SHAFT_D = 8.0             # stainless rod, h6 if bought ground
BEARING = dict(name="608-2RS", d=8.0, D=22.0, B=7.0)    # VERIFY
BEARING_ALT = dict(name="688-2RS", d=8.0, D=16.0, B=5.0)  # VERIFY, compact option
LIP_SEAL = dict(name="TC 8x16x5", d=8.0, D=16.0, B=5.0)   # VERIFY, v2 only

# Generator envelope for v2: low-KV gimbal outrunner (low rpm, many poles).
# Display build prints a replica with the same envelope and mount.
GENERATOR = dict(name="GBM4108-class gimbal motor", D=46.0, L=26.0,  # VERIFY
                 mount_pcd=25.0, mount_screw="M3", n_mount=4)      # VERIFY

MAGNET = dict(D=4.0, H=2.0, grade="N52")   # blade detents, VERIFY
ORING_EQUATOR = dict(cs=2.0)               # cross-section; ID chosen by the groove

# --------------------------------------------------------------------------
# Printing (Bambu Lab, 256 mm cube)
# --------------------------------------------------------------------------
BED = (256.0, 256.0, 256.0)
NOZZLE = 0.4
LAYER = 0.2
MIN_WALL = 1.2            # 3 perimeters
CLEAR_SLIDE = 0.25        # radial clearance, sliding fit (bayonet, fins)
CLEAR_ROTATE = 0.30       # radial clearance, free rotation (blade on pin)
PRESS_FIT = -0.05         # interference for pins / bearing seats, PROVISIONAL
                          # (calibrated by the tolerance coupon before final print)
MAX_OVERHANG = 45.0       # deg from vertical without supports
