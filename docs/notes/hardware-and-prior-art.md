# Hardware selection, v2 generator trade study and prior art

Date: 2026-10-03, revision 3 (second fix round after the adversarial review, see section 10). Data: `data/hardware.json`
(every value with `source_url` and `confidence`) and `data/bom.json` are the record of this study. It is a sourced
desk study: its working scripts (`hw_calc.py`: losses, the coupled rotor + generator scan incl. a time-domain bridge
model, ratings, brake, gland, joint wall, hull fit, prior-art check; `bridge_td.py`: bridge-model arbitration;
`gl40_pads.py`: the GL40 drawing measured from its vector paths; `step_mount.py`: a community motor-mount STEP;
`build_outputs.py`, `write_notes.py`: this file) ran outside the repository against downloaded vendor datasheets
and are not published, so the study cannot be re-run from the repository. Where the text below says "computed in
`hw_calc.py`", the result is recorded in `data/hardware.json`. Inputs from the sibling studies: rotor curves from
`data/blade_design.json` (sha256 65c3c51185a4, design TSR 4.0, Cp peak at TSR 3.822);
blade-pocket depth from `data/layout_rotor.json` (sha256 04636f4af337); joint stations and
groove form from `data/layout_internal.json` (sha256 b346a3455bd4). The study would need a re-run if they change.
Confidence: high = manufacturer datasheet or standard; medium = dealer, community or derived; low = inference.

## 1. Decisions at a glance

| Item | Choice | Key dimensions | Confidence |
|---|---|---|---|
| v2 generator | CubeMars GL40 KV70 (alternative: iPower GM4108H-120T after a bench test) | D 46.5 x L 21.5, 8 mm through bore; stator 4-M3 x 3.5 on PCD 24 + 4-M2.5 x 3.5 on PCD 40, both at 0/90/180/270; pad block at 45 deg; rotor 4-M3 x 2.5 on PCD 24; bay reserve D 47 x L 33 | high (drawing) |
| Bearings | S608-2RS (440C), 688-2RS as compact alternative | 8x22x7; 8x16x5 | high |
| Shaft | 8 mm SUS440C hardened (MISUMI PSSFJ8), g6 | 56 HRC~, g6 = -0.005/-0.014 | medium |
| Shaft axial stop | set-screw collar Ruland MSC-8-SS forward of the forward bearing (or DIN 471 A8) | collar 8 x 16 x 8; A8: s 0.8, groove 7.6 x 0.9 | medium / high |
| Lip seal (v2) | TC 8x16x5 NBR | housing bore 16 (same as 688) | medium |
| Equator O-ring | ISO 3601-141 (AS568-141), EPDM 70 | ID 58.42 +-0.51 x CS 2.62 +-0.08; spigot groove bottom D 59.88, bottom 2.9 / top 4.905 (45 deg tip-side flank), gland depth 2.105; bore 64.09 | high (wall check medium: pockets not final) |
| Detent magnets | N52 disc 4x2 | pull 4.0-4.6 N to steel | medium |
| Hinge pins | ISO 8734 3 x 22 (A or C1), as the rotor layout needs | d 3 m6, c 0.5, l 21.5-22.5 | high |
| Inserts | CNC Kitchen M3 x 5.7 standard | hole 4.0, min wall 1.6 | high |
| Cable gland (v2) | Hummel HSK-K M12x1.5 + thread O-ring + PA locknut | clamp 3-6.5, GL 8, SW 15, H 21 | high |
| Tow eye | printed, integral to the nose | DIN 580 eye is 36 mm tall; test >= 211 N | high (DIN data) |
| Towing state (v2 electronics) | OPEN CONFLICT (superseded in part, see section 2.7): phases open was the first recommendation, but the BEM study shows free-wheel cavitation from 1.6 m/s and a capped load line above the GL40 rating from about 2.1 m/s | to be resolved in v2 (see the report's limitations) | medium |
| Printer | any Bambu A1 / P1S / X1C | 256 x 256 x 256 | high/medium |
| Material | PETG wetted parts; ASA only on enclosed printers | | high |

## 2. v2 generator trade study

### 2.1 Design point, rotor curve and conventions

Design band from the brief: 550-750 rpm and 0.20-0.35 Nm; nominal 650 rpm, 0.27 Nm =
18.38 W shaft. The BEM rotor at 1.5 m/s gives: Cp peak (TSR 3.822) 564.45 rpm /
0.3097 Nm / 18.309 W; TSR 4.0: 590.7 rpm / 0.2956 Nm;
TSR 5.0: 738.3 rpm / 0.2255 Nm. Runaway (zero torque) from the BEM study:
TSR 8.609 (1271.3 rpm) at 1.5 m/s, 1753.2 rpm at 2.0 m/s and
2236.1 rpm (TSR 9.086) towed at 2.5 m/s. Generator envelope: D <= 56 mm, L <= 35 mm.

KV is the hobby convention: no-load rpm per volt of DC bus, i.e. per volt of **line-to-line peak** back-EMF.
Define Kt0 = 60 / (2 pi KV) [Nm/A], so V_ll,peak = Kt0 w. With phase flux linkage lambda (V s/rad, mechanical) the
ll-peak EMF is sqrt3 lambda w, so Kt0 = sqrt3 lambda. Rectifier models:

- **Active (FOC) rectification**, sinusoidal currents of peak I_pk in phase with the EMF: tau = 1.5 lambda I_pk and
  P_cu = 1.5 I_pk^2 R_ph. Eliminating I_pk: `P_cu = 2 R_ph tau^2 / (3 lambda^2) = R_ll (tau / Kt0)^2`, with
  R_ll = 2 R_ph (star). So I = tau / Kt0 is a "DC-equivalent" current, losses use the **line-to-line** resistance, and
  the phase RMS current is I sqrt(2/3). Switch losses are ignored.
- **Passive 6-pulse diode bridge, block model** (used for the full candidate table): sinusoidal EMF, 120 deg current
  blocks, no inductance. The mean DC EMF is (3/pi) V_ll,peak, so tau = (3/pi) Kt0 I_dc and P_cu = I_dc^2 R_ll; two
  Schottky drops of 0.4 V (ASSUMED) conduct. At equal torque its copper loss is (pi/3)^2 = 1.10x the FOC value.
- **Passive bridge, time-domain model** (used wherever the vendor publishes L, and as the GL40 headline): star
  machine with phase R = R_ll/2 and L = L_ll/2, ideal diodes plus Vf, DC side a voltage sink (the input capacitor of a
  buck-boost charger); implicit Euler with an active-set diode solve; the DC voltage is bisected until the mean
  electromagnetic torque equals the rotor torque minus drag, then the rpm is scanned for maximum DC power.

Why the time-domain model. Commutation overlap is real for the GL40 (1.8 mH, 14 pole pairs: X_ph ~1.2 ohm against
R_ph 2.25 ohm near 900 rpm). The textbook overlap drop, (3/pi) w_e L_ph I_dc, assumes a stiff DC current. With a
capacitor link the inductance also smooths the current pulses, which lowers the RMS current and the copper loss.
`bridge_td.py` compares the models at 1.5 m/s:

| Case (1.5 m/s) | Time-domain best (W) | Block model, no L (W) | Overlap as a lossy drop (W) | Overlap, lossless textbook drop (W) | TD vs block model |
|---|---|---|---|---|---|
| GL40 KV70 | 7.65 (fine step 7.67) | 7.78 | 6.92 | 7.15 | -1.6% |
| GL40 KV70, R only (L = 0) | 7.89 (fine step 7.89) | 7.78 | 7.78 | 7.78 | +1.3% |
| GM4108H KV27 R_ll 11.1, R only | 11.23 (fine step 11.23) | 11.23 | 11.23 | 11.23 | -0.0% |

Without inductance the time-domain result is within 1.3% of the block model (it also matches the
review's independent R-only check, 7.89 W for the GL40). Adding the GL40's inductance costs
1.6% against the block model. Treating the overlap drop as dissipated (the review's estimate,
6.92 W) or as the lossless stiff-current drop (7.15 W) understates a capacitor-fed
bridge by 7-11%. The lossless value applies to a bridge with a stiff DC current (a DC choke);
the dissipative one is a lower bound. Convergence: a 4x finer time step over 16 cycles changes the GL40 result by
+0.2%.

Iron and friction loss is modelled as a constant drag torque Q_if = P_no-load / w_test from the published no-load input
(GM4108H: 0.07 A x 20 V at 540 rpm gives 0.0248 Nm; MN5008 and MN4006 from their idle currents).
Candidates with no published no-load data get the GM4108H value, labelled ASSUMED. A linear-in-rpm loss overstates
the iron loss at low rpm if part of it is eddy current, so this errs conservative.

Continuous ratings, where published (CubeMars GL40, GL40 II, RI50), cap the scan. The GL40's rated current is in the
torque-producing convention of this study (Kt 0.150 x 1.62 A = 0.243 Nm against a rated 0.25 Nm), so the caps are
I_dc <= 1.62 A (copper loss <= 11.81 W) and electromagnetic torque <= 0.25 Nm.

Where a vendor gives "internal resistance" without a convention, both readings are carried (value = R_ll, or value =
R_phase so R_ll = 2 x value). Km = Kt0 / sqrt(R_ll) is the figure of merit: P_cu / P_shaft = tau / (Km^2 w) (FOC).
Km needed for a copper-only efficiency target at the nominal point: 60%: 0.0996,
70%: 0.115, 80%: 0.1408, 90%: 0.1992 Nm/sqrt(W).

Holding-torque check (used in the bench test): a DC current I into one lead and out of another gives
tau(theta) = lambda I [sin theta - sin(theta - 120 deg)] = sqrt3 lambda I cos(theta - 60 deg), so the peak holding
torque is sqrt3 lambda I = **Kt0 x I**.

### 2.2 Candidates

| Candidate | Class | D x L (mm) | Shaft (mm) | Mount pattern | KV (rpm/V) | R published (ohm) | Mass (g) | Price | In envelope | Source |
|---|---|---|---|---|---|---|---|---|---|---|
| iPower GM4108H-120T | gimbal outrunner | 47.0 x 22.7 | hollow, OD 10 (-0.008/-0.012), ID 8 (+0.05/0) | base: 4 x M3 clearance holes on PCD 38.14 (adjacent 26.97), centre clearance 16 (one community source, possibly another motor revision); rotor face: unconfirmed (low) | 27.0/32.0 | 11.100 | 124 | USD 35.9, EUR 43.9 | yes | [src](https://shop.iflight.com/ipower-motor-gm4108h-120t-brushless-gimbal-motor-pro217) |
| Rctimer GBM4114HS 100T | gimbal outrunner | 46.0 x 31.0 | 5.0 hollow | not published (medium) | 32.0 | 15.200 | 150 | USD 34.99 | yes | [src](https://www.rctimer.com/4114-hollow-shaft-brushless-gimbal-motor-p0440.html) |
| CubeMars GL40 KV70 | gimbal outrunner | 46.5 x 21.5 | 8 through bore (shown on both faces) | rotor face: 4-M3 x 2.5 deep on PCD 24, centre bore 8; stator face: 4-M3 x 3.5 deep on PCD 24 + 4-M2.5 x 3.5 deep on PCD 40, both patterns on the same axes (0/90/180/270); phase solder-pad block on the 45 deg bisector between two holes, r 12.9-19.7 (high) | 63.5 | 4.500 | 107 | USD 72.99 | yes | [src](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf) |
| CubeMars GL40 II KV82.5 | gimbal outrunner (integrated driver) | 46.1 x 33.5 | 12.5 hollow | not extracted (low) | 82.5 | 3.000 | 125 | USD 133.99 | yes | [src](https://cubemars.com/goods-1182-GL40+II.html) |
| T-motor GB4106 KV53 | gimbal outrunner | 47.8 x 20.0 | hole size 5.6 | not extracted (low) | 53.0 | 15.300 | 70 | USD 45.9 | yes | [src](https://store.tmotor.com/product/gb4106-gimbal-type.html) |
| EMAX GB4008 66KV | gimbal outrunner | n/a x n/a | hollow | base centre pitch 19x19 or 25x25 (medium) | 66.0 | 7.500 | 101 | USD 23.49 | unknown | [src](https://shop.emaxmodel.com/products/emax-gb4008-brushless-gimbal-motor-hollow-shaft-66kv-90t) |
| CubeMars RI50 KV100 (frameless) | frameless inrunner torque motor | 54.0 x 27.0 | 22 hollow (frameless rotor) | frameless: no housing or bearings, user housing required (high) | 100.0 | 1.420 | 181 | USD 61.9 | yes | [src](https://www.cubemars.com/product/ri50-kv100-frameless-inrunner-torque-motor.html) |
| T-motor MN5008 KV170 | low-KV drone outrunner | 55.6 x 32.0 | 6 (inner) | not extracted (low) | 170.0 | 0.270 | 128 | USD 89.99 | yes | [src](https://store.tmotor.com/product/mn5008-kv170-motor-antigravity-type.html) |
| T-motor MN4006 KV380 | low-KV drone outrunner | 44.4 x 21.0 | not extracted | not extracted (low) | 380.0 | 0.194 | 68 | USD 149.9 | yes | [src](https://store.tmotor.com/product/mn4006-kv380-motor-antigravity-type.html) |
| maxon EC 45 flat 70 W 36 V (411815) | industrial flat outrunner | 45.0 x n/a | not extracted | not extracted (low) | 179.0 | 1.740 | 141 | EUR 226.94 | unknown | [src](https://www.maxongroup.com/maxon/view/product/motor/ecmotor/ecflat/ecflat45/411815) |
| iPower GM5208-24 (outside envelope) | gimbal outrunner | 59.5 x 31.2 | hollow OD 15, ID 12.6 | not published (medium) | 20.8 | 13.700 | 204 | USD 51.0 | no | [src](https://iflighteurope.com/products/ipower-gm5208-24-gimbal-motor) |
| T-motor GB54-2 KV26 (outside envelope) | gimbal outrunner | 60.7 x 25.0 | hollow, hole 12.7 | not extracted (low) | 26.0 | 15.000 | 156 | USD 74.9 | no | [src](https://store.tmotor.com/goods-445-T-Motor+GB54-2++Gimbal+Type+3-6S+UAV+Drone+Motor+26KV.html) |

Categories covered: gimbal outrunners (GM4108H, GBM4114HS, GL40, GL40 II, GB4106, GB4008), low-KV drone outrunners
(MN5008, MN4006), a frameless torque motor (RI50) and an industrial flat BLDC (EC 45 flat). For small hub or axial-flux
generators, no off-the-shelf unit with a datasheet exists within D 56 mm; search hits start at 50 W to 20 kW industrial
units ([search hit](https://www.permanent-magnetalternator.com/sale-13975594-low-rpm-axial-flux-coreless-permanent-magnet-generator-50w-20kw.html)). A custom PCB-stator axial-flux generator is a v3 idea, not a part.
Two out-of-envelope motors (GM5208-24, GB54-2) are kept as references: they show what a 60 mm hull would buy.
"Not extracted" means the vendor publishes a drawing image or file that was not parsed; "not published" means the
page has no hole pattern at all.

### 2.3 Losses at the nominal point (650 rpm, 0.27 Nm), iron loss excluded

| Candidate | KV | R_ll | Km (Nm/sqrtW) | Voc ll-pk / ll-rms (V) | I FOC / I_ph,rms (A) | P_cu FOC (W) | eta_cu FOC | I_dc bridge (A) | P_cu + P_diode bridge (W) | eta bridge | V_dc bridge (V) | eta_cu FOC over 550-750 rpm, 0.20-0.35 Nm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM4108H | 27.0 | 11.10 | 0.1062 | 24.1 / 17.0 | 0.76 / 0.62 | 6.5 | 0.65 | 0.80 | 7.1 + 0.64 | 0.58 | 13.3 | 0.46 to 0.77 |
| GM4108H | 27.0 | 22.20 | 0.0751 | 24.1 / 17.0 | 0.76 / 0.62 | 12.9 | 0.30 | 0.80 | 14.2 + 0.64 | 0.19 | 4.4 | -0.08 to 0.55 |
| GM4108H | 32.0 | 11.10 | 0.0896 | 20.3 / 14.4 | 0.91 / 0.74 | 9.1 | 0.51 | 0.95 | 10.0 + 0.76 | 0.42 | 8.1 | 0.24 to 0.68 |
| GM4108H | 32.0 | 22.20 | 0.0633 | 20.3 / 14.4 | 0.91 / 0.74 | 18.2 | 0.01 | 0.95 | 19.9 + 0.76 | -0.13 | -2.4 | -0.52 to 0.36 |
| GBM4114HS | 32.0 | 15.20 | 0.0765 | 20.3 / 14.4 | 0.91 / 0.74 | 12.4 | 0.32 | 0.95 | 13.7 + 0.76 | 0.22 | 4.2 | -0.04 to 0.56 |
| GBM4114HS | 32.0 | 30.40 | 0.0541 | 20.3 / 14.4 | 0.91 / 0.74 | 24.9 | -0.35 | 0.95 | 27.3 + 0.76 | -0.53 | -10.2 | -1.07 to 0.13 |
| GL40 KV70 | 63.5 | 4.50 | 0.0709 | 10.2 / 7.2 | 1.79 / 1.47 | 14.5 | 0.21 | 1.88 | 15.9 + 1.50 | 0.05 | 0.5 | -0.21 to 0.49 |
| GL40 II | 82.5 | 3.00 | 0.0668 | 7.9 / 5.6 | 2.33 / 1.91 | 16.3 | 0.11 | 2.44 | 17.9 + 1.95 | -0.08 | -0.6 | -0.36 to 0.43 |
| GL40 II | 86.8 | 3.00 | 0.0635 | 7.5 / 5.3 | 2.46 / 2.00 | 18.1 | 0.02 | 2.57 | 19.8 + 2.06 | -0.19 | -1.4 | -0.51 to 0.37 |
| GB4106 | 53.0 | 15.30 | 0.0461 | 12.3 / 8.7 | 1.50 / 1.22 | 34.4 | -0.87 | 1.57 | 37.7 + 1.26 | -1.12 | -13.1 | -1.86 to -0.20 |
| GB4008 | 66.0 | 7.50 | 0.0528 | 9.8 / 7.0 | 1.87 / 1.52 | 26.1 | -0.42 | 1.95 | 28.6 + 1.56 | -0.64 | -6.0 | -1.18 to 0.09 |
| RI50 | 100.0 | 1.42 | 0.0801 | 6.5 / 4.6 | 2.83 / 2.31 | 11.3 | 0.38 | 2.96 | 12.4 + 2.37 | 0.19 | 1.2 | 0.05 to 0.60 |
| RI50 | 79.6 | 1.42 | 0.1007 | 8.2 / 5.8 | 2.25 / 1.84 | 7.2 | 0.61 | 2.36 | 7.9 + 1.88 | 0.47 | 3.6 | 0.40 to 0.75 |
| MN5008 | 170.0 | 0.27 | 0.1081 | 3.8 / 2.7 | 4.81 / 3.92 | 6.2 | 0.66 | 5.03 | 6.8 + 4.03 | 0.41 | 1.5 | 0.48 to 0.78 |
| MN4006 | 380.0 | 0.19 | 0.0571 | 1.7 / 1.2 | 10.74 / 8.77 | 22.4 | -0.22 | 11.25 | 24.6 + 9.00 | -0.83 | -1.4 | -0.87 to 0.22 |
| EC 45 flat | 179.0 | 1.74 | 0.0404 | 3.6 / 2.6 | 5.06 / 4.13 | 44.6 | -1.43 | 5.30 | 48.9 + 4.24 | -1.89 | -6.5 | -2.71 to -0.56 |
| GM5208 (out) | 20.8 | 13.70 | 0.1240 | 31.2 / 22.1 | 0.59 / 0.48 | 4.7 | 0.74 | 0.62 | 5.2 + 0.49 | 0.69 | 20.6 | 0.60 to 0.83 |
| GM5208 (out) | 20.8 | 27.40 | 0.0877 | 31.2 / 22.1 | 0.59 / 0.48 | 9.5 | 0.48 | 0.62 | 10.4 + 0.49 | 0.41 | 12.2 | 0.21 to 0.67 |
| GB54-2 (out) | 26.0 | 15.00 | 0.0948 | 25.0 / 17.7 | 0.73 / 0.60 | 8.1 | 0.56 | 0.77 | 8.9 + 0.62 | 0.48 | 11.5 | 0.32 to 0.72 |

A fixed operating point penalises high-resistance and low-voltage machines: negative efficiencies mean the motor
cannot absorb 0.27 Nm at 650 rpm with positive output. The GL40 is an example: at 650 rpm its bridge output is only
0.51 V (block model). A real charger tracks the best point, so 2.4 is the comparison that matters.

### 2.4 Coupled rotor + generator (MPPT load)

Each variant loads the BEM rotor (`Cp_final_by_flow`, scanned from TSR 2.0 to runaway) and the scan finds the rpm
of maximum DC power, with iron + friction drag included. At 1.5 m/s, block model for every candidate:

| Candidate | KV | R_ll | Km | Iron drag (Nm) | Bridge: P_el DC (W) | rpm | TSR | V_dc (V) | I_dc (A) | eta | Active: P_el (W) | TSR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GM4108H | 27.0 | 11.10 | 0.1062 | 0.0248 (pub.) | **11.23** | 781 | 5.29 | 20.9 | 0.54 | 0.66 | 11.97 | 5.11 |
| GM4108H | 27.0 | 22.20 | 0.0751 | 0.0248 (pub.) | **8.67** | 859 | 5.82 | 19.9 | 0.43 | 0.56 | 9.40 | 5.76 |
| GM4108H | 32.0 | 11.10 | 0.0896 | 0.0248 (pub.) | **9.97** | 823 | 5.58 | 17.4 | 0.57 | 0.61 | 10.79 | 5.41 |
| GM4108H | 32.0 | 22.20 | 0.0633 | 0.0248 (pub.) | **7.26** | 923 | 6.25 | 17.6 | 0.41 | 0.53 | 7.98 | 6.04 |
| GBM4114HS | 32.0 | 15.20 | 0.0765 | 0.0248 (assumed) | **8.76** | 863 | 5.85 | 17.2 | 0.51 | 0.57 | 9.56 | 5.69 |
| GBM4114HS | 32.0 | 30.40 | 0.0541 | 0.0248 (assumed) | **6.07** | 969 | 6.56 | 17.7 | 0.34 | 0.49 | 6.71 | 6.37 |
| GL40 KV70 | 63.5 | 4.50 | 0.0709 | 0.0248 (assumed) | **7.78** | 894 | 6.06 | 8.6 | 0.91 | 0.54 | 8.92 | 5.79 |
| GL40 II | 82.5 | 3.00 | 0.0668 | 0.0248 (assumed) | **7.13** | 916 | 6.21 | 6.5 | 1.09 | 0.51 | 8.43 | 5.95 |
| GL40 II | 86.8 | 3.00 | 0.0635 | 0.0248 (assumed) | **6.72** | 923 | 6.25 | 6.0 | 1.12 | 0.49 | 8.00 | 6.10 |
| GB4106 | 53.0 | 15.30 | 0.0461 | 0.0248 (assumed) | **4.80** | 1011 | 6.85 | 10.3 | 0.47 | 0.43 | 5.50 | 6.71 |
| GB4008 | 66.0 | 7.50 | 0.0528 | 0.0248 (assumed) | **5.61** | 984 | 6.67 | 8.5 | 0.66 | 0.47 | 6.52 | 6.50 |
| RI50 | 100.0 | 1.42 | 0.0801 | 0.0248 (assumed) | **8.24** | 866 | 5.87 | 5.2 | 1.57 | 0.54 | 9.92 | 5.67 |
| RI50 | 79.6 | 1.42 | 0.1007 | 0.0248 (assumed) | **10.08** | 804 | 5.45 | 6.7 | 1.50 | 0.61 | 11.63 | 5.23 |
| MN5008 | 170.0 | 0.27 | 0.1081 | 0.0225 (pub.) | **9.24** | 811 | 5.49 | 2.9 | 3.19 | 0.56 | 12.19 | 5.20 |
| MN4006 | 380.0 | 0.19 | 0.0571 | 0.0075 (pub.) | **4.01** | 1043 | 7.07 | 1.1 | 3.51 | 0.40 | 7.61 | 6.56 |
| EC 45 flat | 179.0 | 1.74 | 0.0404 | 0.0248 (assumed) | **3.31** | 1056 | 7.15 | 2.7 | 1.21 | 0.34 | 4.61 | 6.95 |
| GM5208 (out) | 20.8 | 13.70 | 0.1240 | 0.0248 (assumed) | **12.28** | 740 | 5.01 | 26.9 | 0.46 | 0.70 | 12.90 | 5.00 |
| GM5208 (out) | 20.8 | 27.40 | 0.0877 | 0.0248 (assumed) | **9.98** | 819 | 5.55 | 26.5 | 0.38 | 0.61 | 10.63 | 5.51 |
| GB54-2 (out) | 26.0 | 15.00 | 0.0948 | 0.0248 (assumed) | **10.48** | 804 | 5.45 | 21.4 | 0.49 | 0.63 | 11.20 | 5.35 |

Time-domain bridge for the motors with a published inductance, and for the GM4108H with an ESTIMATED one (its L is
unpublished; the estimate gives it the GL40's L/R time constant of 0.40 ms):

| Generator reading | L_ll (mH) | L basis | R-only block model 1.5 m/s (W) | Time-domain 1.5 m/s: W | rpm | TSR | V_dc | I_dc | Change |
|---|---|---|---|---|---|---|---|---|---|
| GM4108H (KV 27.0, R_ll 11.10) | 4.44 | ESTIMATE (GL40 L/R) | 11.23 | **11.19** | 781 | 5.29 | 20.29 | 0.552 | -0.4% |
| GM4108H (KV 27.0, R_ll 22.20) | 8.88 | ESTIMATE (GL40 L/R) | 8.67 | **8.64** | 859 | 5.82 | 19.10 | 0.453 | -0.3% |
| GM4108H (KV 32.0, R_ll 11.10) | 4.44 | ESTIMATE (GL40 L/R) | 9.97 | **9.94** | 823 | 5.57 | 16.82 | 0.591 | -0.3% |
| GM4108H (KV 32.0, R_ll 22.20) | 8.88 | ESTIMATE (GL40 L/R) | 7.26 | **7.22** | 923 | 6.25 | 16.70 | 0.433 | -0.6% |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 1.80 | published | 7.78 | **7.66** | 894 | 6.05 | 7.96 | 0.962 | -1.5% |
| GL40 II (KV 82.5, R_ll 3.00) | 1.32 | published | 7.13 | **6.95** | 916 | 6.20 | 5.96 | 1.166 | -2.5% |
| GL40 II (KV 86.8, R_ll 3.00) | 1.32 | published | 6.72 | **6.53** | 943 | 6.39 | 5.90 | 1.107 | -2.8% |
| RI50 (KV 100.0, R_ll 1.42) | 1.50 | published | 8.24 | **8.00** | 870 | 5.89 | 4.82 | 1.661 | -2.9% |
| RI50 (KV 79.6, R_ll 1.42) | 1.50 | published | 10.08 | **9.91** | 812 | 5.50 | 6.41 | 1.545 | -1.7% |

GL40 across flows (time-domain bridge), against its continuous rating of 0.25 Nm / 1.62 A (11.81 W copper):

| Flow (m/s) | TD bridge W | rpm (TSR) | V_dc / I_dc | em torque (Nm) | P_cu (W) | Within rating | Capped at rating: W | rpm (TSR) | V_dc / I_dc | P_cu (W) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 2.12 | 499 (5.07) | 4.55 / 0.466 | 0.065 | 0.92 | yes | same (cap does not bind) | | | |
| 1.5 | 7.66 | 894 (6.05) | 7.96 / 0.962 | 0.131 | 3.78 | yes | same (cap does not bind) | | | |
| 2.0 | 17.05 | 1328 (6.74) | 11.48 / 1.484 | 0.195 | 8.75 | yes | same (cap does not bind) | | | |
| 2.5 | 30.30 | 1780 (7.23) | 14.52 / 2.086 | 0.263 | 16.84 | no | **29.47** | 1849 (7.51) | 17.19 / 1.715 | 11.62 |

Block model across flows for the main readings (V_dc at the optimum in brackets; "over rating" marks points above a
published continuous rating):

| Generator reading | 1.0 m/s bridge W (V_dc) | 1.5 m/s bridge W (V_dc) | 2.0 m/s bridge W (V_dc) | 2.5 m/s bridge W (V_dc) | 1.5 m/s active W |
|---|---|---|---|---|---|
| GM4108H (KV 27.0, R_ll 11.10) | 2.98 (11.5) | 11.23 (20.9) | 26.30 (28.5) | 49.33 (37.9) | 11.97 |
| GM4108H (KV 27.0, R_ll 22.20) | 2.43 (12.6) | 8.67 (19.9) | 19.55 (30.0) | 35.89 (39.0) | 9.40 |
| GM4108H (KV 32.0, R_ll 11.10) | 2.69 (10.5) | 9.97 (17.4) | 22.89 (24.7) | 42.62 (32.4) | 10.79 |
| GM4108H (KV 32.0, R_ll 22.20) | 2.10 (10.3) | 7.26 (17.6) | 16.21 (24.6) | 29.39 (32.2) | 7.98 |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 2.13 (5.1) | 7.78 (8.6) | 17.77 (12.2) | 32.77 (15.9) (over rating) | 8.92 |
| MN5008 (KV 170.0, R_ll 0.27) | 2.10 (1.6) | 9.24 (2.9) | 22.95 (4.3) | 44.74 (5.7) | 12.19 |

Findings:

1. With these generators the electrical optimum at 1.5 m/s is at TSR 5.29-6.25 for the GL40 and GM4108H
   (5.01-7.15 across all candidates), above the rotor's Cp peak
   (TSR 3.822): the MPPT trades rotor Cp for lower torque and copper loss. A rotor whose Cp peak sat
   nearer TSR 5-5.5 would recover part of that (BEM study's call; `tsr_sensitivity`: at equal shaft power a Km 0.07
   machine loses 81% to copper at TSR 4 and 36% at TSR 6).
2. In the coupled view the in-envelope candidates with Km ~0.07-0.11 land within a few watts of each other. What
   separates them is voltage, current and data quality. The GM4108H and GB54-2 deliver 17-21 V at
   0.4-0.6 A; the GL40 8.0 V at 0.96 A; the MN5008 and RI50 deliver 3-7 V at
   1.5-3.2 A, which needs an in-hull active rectifier to stay efficient.
3. The GL40 stays inside its continuous rating up to 2.0 m/s. Towed at 2.5 m/s the unconstrained optimum needs
   0.263 Nm and 16.8 W of copper loss; a charger current limit of 1.62 A moves the point to
   TSR 7.51 and costs 0.83 W (29.47 W instead of 30.30 W).

### 2.5 GM4108H: the datasheet contradicts itself

The iFlight page gives no KV. Its no-load line (513-567 rpm at 20 V) implies KV 27 (SVPWM reading) or 32 (sine-PWM
reading, also the value a SimpleFOC user runs), i.e. Kt0 0.2984-0.3537 Nm/A.
The same page lists "Load current 1.5A, Load volts 20V, Load torque(g.cm) 1200-1800", i.e. 0.1177-0.1765 Nm at
1.5 A, Kt 0.0785-0.1177 Nm/A, equivalent KV 81.1-121.6.
The KV readings predict 4564-5410 g.cm at 1.5 A,
2.54-4.51x the listed value. The two lines cannot both describe the
same back-EMF: with Kt from the load point the no-load speed at 20 V would be 1622-2432 rpm,
not ~540. If the load point is the right one, the coupled output collapses:

| KV (from load point) | R_ll | 1.5 m/s bridge P_el (W) | 1.5 m/s active P_el (W) |
|---|---|---|---|
| 121.6 | 11.1 | 1.51 | 1.95 |
| 121.6 | 22.2 | 0.80 | 1.04 |
| 81.1 | 11.1 | 3.14 | 3.78 |
| 81.1 | 22.2 | 1.77 | 2.16 |

Two further cross-checks:
- Size benchmark: the GL40 KV70 (46.5 x 21.5, KV and Kt agree on its datasheet) has Km 0.0709. The GM4108H readings
  span 0.0633-0.1062. The per-phase readings (R_ll 22.2) sit next to the GL40; the line-to-line
  readings would beat it by a quarter to a half. Inference, low confidence: Km also depends on laminations, magnets and
  slot fill.
- Mount provenance: the base pattern (4 x M3 on PCD 38.14) comes from one SimpleFOC user (csunwc), who posts
  "27 mm apart" in the forum and published the bracket STEP measured here. The same user reports a hollow shaft
  OD 15 / ID 8, whereas iFlight lists OD 10; the STEP's 16 mm centre hole suits the 15 mm shaft. It may be another
  revision. Confidence is low, and the rotor-face pattern is unconfirmed.

So the GM4108H numbers in 2.3 and 2.4 are conditional on the KV-from-no-load reading, and its block-model bridge
numbers are optimistic until its inductance is measured (time-domain with the ESTIMATED L: 7.22-11.19 W).

### 2.6 Recommendation: CubeMars GL40 KV70

Why:
- Only in-envelope candidate whose datasheet is complete and self-consistent: KV 63.5 and Kt 0.150 Nm/A agree (Kt0 = 0.1504), rated 0.25 Nm at 1.62 A matches, R and L are stated phase-to-phase, and a manufacturer drawing gives both mounting faces. Its trade-study numbers carry no convention spread.
- Coupled to the BEM rotor at 1.5 m/s it delivers 7.66 W DC through a passive diode bridge (7.96 V, 0.962 A, TSR 6.05), from a time-domain bridge model that includes its published 1.8 mH inductance (R-only block model: 7.78 W), and 8.92 W with active rectification; iron loss included (ASSUMED equal to the GM4108H drag torque). It runs at 52% of its rated torque there. That sits inside the GM4108H band (7.26-11.23 W, R-only; 7.22-11.19 W with an ESTIMATED inductance) and matches its readings with a size-plausible Km <= 0.08 (7.26-8.67 W).
- Low voltage at towing speed: freewheel at the BEM runaway speed is 35.2 V peak line-to-line at 2.5 m/s (GM4108H: 69.9-82.8 V), so a 60 V-rated input stage can simply leave the rotor open-circuit when towed; no brake is needed for voltage protection.
- D 46.5 x L 21.5 fits the hull under the layout_rotor blade pockets from x 19.0 to 93.0 (2 mm margin); in the equator joint the spigot bore of the -141 O-ring (55.88 with the layout study's lug-zone rule) leaves 4.44 mm radial clearance around the D 47 bay. 8 mm through bore matches SHAFT_D; the rotor face takes a bolted shaft hub on 4-M3 PCD 24.

Risks:
- Price: USD 72.99 against USD 35.90 for the GM4108H.
- No no-load loss is published for the GL40: its iron + friction drag is ASSUMED equal to the GM4108H value (0.0248 Nm). Measure the spin-down or no-load current on receipt.
- Continuous rating: 0.25 Nm / 1.62 A (rated copper loss 11.81 W). The unconstrained MPPT point towed at 2.5 m/s needs 0.2629 Nm and 16.84 W of copper loss, above the rating; capped at the rating the charger runs at TSR 7.51 and still delivers 29.47 W (time-domain model). At 2.0 m/s and below the optimum is inside the rating.
- Do not short the phases as a continuous tow brake: at the short-circuit equilibrium it dissipates 69.6 W in the windings at 2.5 m/s (5.89x the rated copper loss, 3.93 A against the 3.3 A peak rating) and already 16.7 W (1.41x) at 1.5 m/s, inside a sealed PETG hull. A short is acceptable only as a parking brake at <= 1.0 m/s (3.5 W, 0.3x).
- At 1.0 m/s the bridge output is 4.55 V at the optimum: a 5 V USB stage needs a buck-boost (or SEPIC) input, not a plain buck.
- Which part of the GL40 carries the 8 mm bore (rotating hollow shaft or stator hub) is inferred from the drawing (both faces show an 8 mm bore); torque transfer through the bolted rotor-face hub works either way. VERIFY on receipt.
- If the 8 mm shaft passes through the motor, the motor bearings become a third support; use a compliant coupling or let the motor float on its adapter plate.

Fixed-point totals at 650 rpm / 0.27 Nm with the iron estimate (passive bridge, block model; the coupled scan is the fairer view):

| Motor | KV | R_ll | Km | P_iron+friction (W) | P_el DC bridge (W) | eta total (est.) |
|---|---|---|---|---|---|---|
| GM4108H | 27.0 | 11.1 | 0.1062 | 1.69 | 8.96 | 0.49 |
| GM4108H | 27.0 | 22.2 | 0.0751 | 1.69 | 1.86 | 0.10 |
| GM4108H | 32.0 | 11.1 | 0.0896 | 1.69 | 5.97 | 0.33 |
| GM4108H | 32.0 | 22.2 | 0.0633 | 1.69 | -4.00 | -0.22 |
| GL40 KV70 | 63.5 | 4.5 | 0.0709 | 1.69 | -0.72 | -0.04 |

Alternative: iPower GM4108H-120T (same bay, own adapter plate; only after a bench test confirms Km >= 0.07 and measures its inductance). It is half the price and its high voltage suits a long tether, but every number
for it rests on the bench test:
- The base hole pattern is from a single community source whose motor may be another revision (shaft OD 15 vs iFlight's 10). The rotor face pattern is unconfirmed.
- Bench test that settles all of it in 10 minutes: lead-to-lead R (= R_ll); spin at a measured rpm and read the AC voltage between two leads (KV = rpm / (V_rms x sqrt2)); lead-to-lead inductance with an LCR meter (= L_ll); optionally pass a DC current I through two leads and read the peak holding torque on a lever arm, which equals Kt0 x I for a sinusoidal-EMF motor (derived in the notes); at 1.5 A the KV readings predict 0.448-0.531 Nm, the datasheet 0.1177-0.1765 Nm.

Runner-up outside the gimbal class: T-motor MN5008 KV170 only with in-hull active rectification and boost: at 1.5 m/s it gives 12.19 W active but 9.24 W through a diode bridge at 2.9 V and 3.188 A. D 55.6 fits the hull under the pockets only from x 27.5 to 79.5 and not through the equator joint, so it would sit wholly in one shell.

**Mount interface for CAD** (all mm; GL40 from the manufacturer drawing):

| Feature | Value | Confidence | Source |
|---|---|---|---|
| Body | D 46.5 (rotor bell; stator plate 45.5) x L 21.5 (body 21) | high | [GL40 drawing](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf) |
| Bay reserve | D 47.0 x L 33.0 (fits the GM4108H alternative too) | medium | `hw_calc.py` |
| Radial clearance motor to hull | 2.0 min | design margin | |
| Axial window (2 mm margin, layout_rotor pockets) | D 47: x 19.0 to 92.5 (73.5 long); GL40 D 46.5: x 19.0 to 93.0. Inside the joint (seal zone x 49.495-57.5) a bay that straddles the equator sits in the spigot bore: 55.88 with the layout study's lug-zone rule (57.48 at a 1.2 mm wall under the groove), i.e. 4.44 mm radial clearance around D 47 | medium | `hw_calc.py` |
| Stator face | 4 x M3 x 3.5 deep on PCD 24.0 and 4 x M2.5 x 3.5 deep on PCD 40.0, **both at 0/90/180/270 deg** (same axes) | high | [drawing](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf), `gl40_pads.py` |
| Solder-pad keep-out (stator face) | pad block on the 45 deg bisector between two adjacent holes: radial 12.93-19.05, tangential +-5.09 (r 12.93-19.72, +-21.5 deg); three ~2 x 2 pads at radial 16.6-18.6. CAD window in the adapter plate and replica: radial 12.0-20.5, tangential +-6.0 | high (window margin: design) | [drawing](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf), `gl40_pads.py` |
| Rotor face | 4 x M3 x 2.5 deep on PCD 24.0; centre bore 8.0 | high | [drawing](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf) |
| Screw length | clamped part + 3.0 (stator) / + 2.0 (rotor) max | medium | thread depths |
| Shaft | 8 mm g6 through the 8 mm bore | high | [drawing](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf) |
| Rotor coupling | hub bolted to the rotor face (4 x M3, PCD 24), M3 grub screw on a shaft flat; 8 mm axial zone aft of the rotor face | medium | design proposal |
| Adapter plate | printed, motor-specific, between bay bulkhead and stator face, with the pad window; clock it so the window faces the cable route; GM4108H version 4 x M3 on PCD 38.14 (low, VERIFY) | medium | design proposal |

The display replica copies D 46.5 x 21.5, both hole patterns and the pad block. Orient the stator face forward so the
phase leads run to the nose gland, with the rotor face and its shaft hub toward the bearings. If the layout clocks the
hole pattern at 45 deg (as the internal layout does, to keep screws off the cable path at theta 180), the pad block lands
on an axis, e.g. at theta 180 next to the cable notch.

Generator bench test before freezing v2 (10 minutes, multimeter + LCR meter + drill + phone tachometer):
1. R between any two leads at 20 C = R_ll; L between any two leads = L_ll.
2. Spin at a measured rpm and read the AC volts between two leads: KV = rpm / (V_rms x sqrt2).
3. Optional: hold a DC current I through two leads and read the peak holding torque on a lever arm: it should be Kt0 x I.
4. Recompute Km = 60 / (2 pi KV) / sqrt(R_ll) and re-run `hw_calc.py` with the measured values. For the GL40, measure the
   no-load current or a spin-down to replace the assumed iron drag.

### 2.7 Freewheel, brake and towing

> **Superseded in part (2026-10-06).** Electrically, the open-circuit freewheel below is safe. Hydrodynamically
> it is not: the blade study finds that a free-wheeling rotor cavitates from 1.6 m/s, and the layout study that
> the fins then roll the hull, so v2 must keep the rotor loaded when towed. The 10.02 ohm dump-load point
> (TSR 7.51 at 2.5 m/s) also lies in the cavitating region, and from about 2.1 m/s the cavitation-capped load
> line needs more than the GL40's continuous rating. The towing state is an open v2 issue: see the technical
> report (generator section and limitations).

An unloaded rotor runs to its zero-torque (runaway) speed, so the open-circuit voltage at runaway sizes the
electronics:

| Generator reading | Voc ll-pk at runaway 1.5 / 2.0 / 2.5 m/s (V) | L included in the brake |
|---|---|---|
| GM4108H (KV 27.0, R_ll 11.10) | 47.1 / 64.9 / 82.8 | no (L unpublished) |
| GM4108H (KV 27.0, R_ll 22.20) | 47.1 / 64.9 / 82.8 | no (L unpublished) |
| GM4108H (KV 32.0, R_ll 11.10) | 39.7 / 54.8 / 69.9 | no (L unpublished) |
| GM4108H (KV 32.0, R_ll 22.20) | 39.7 / 54.8 / 69.9 | no (L unpublished) |
| GBM4114HS (KV 32.0, R_ll 15.20) | 39.7 / 54.8 / 69.9 | no (L unpublished) |
| GBM4114HS (KV 32.0, R_ll 30.40) | 39.7 / 54.8 / 69.9 | no (L unpublished) |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 20.0 / 27.6 / 35.2 | yes |
| GL40 II (KV 82.5, R_ll 3.00) | 15.4 / 21.3 / 27.1 | yes |
| GL40 II (KV 86.8, R_ll 3.00) | 14.6 / 20.2 / 25.8 | yes |
| GB4106 (KV 53.0, R_ll 15.30) | 24.0 / 33.1 / 42.2 | no (L unpublished) |
| GB4008 (KV 66.0, R_ll 7.50) | 19.3 / 26.6 / 33.9 | no (L unpublished) |
| RI50 (KV 100.0, R_ll 1.42) | 12.7 / 17.5 / 22.4 | yes |
| RI50 (KV 79.6, R_ll 1.42) | 16.0 / 22.0 / 28.1 | yes |
| MN5008 (KV 170.0, R_ll 0.27) | 7.5 / 10.3 / 13.2 | no (L unpublished) |
| MN4006 (KV 380.0, R_ll 0.19) | 3.3 / 4.6 / 5.9 | no (L unpublished) |
| EC 45 flat (KV 179.0, R_ll 1.74) | 7.1 / 9.8 / 12.5 | no (L unpublished) |
| GM5208 (out) (KV 20.8, R_ll 13.70) | 61.1 / 84.3 / 107.5 | no (L unpublished) |
| GM5208 (out) (KV 20.8, R_ll 27.40) | 61.1 / 84.3 / 107.5 | no (L unpublished) |
| GB54-2 (out) (KV 26.0, R_ll 15.00) | 48.9 / 67.4 / 86.0 | no (L unpublished) |

GL40: 20.0 / 27.6 / 35.2 V at 1.5 / 2.0 / 2.5 m/s. A 60 V-rated input stage therefore
tolerates the open-circuit freewheel when towed, which made phases open the first recommended towing state (superseded in part, see the banner in section 2.7: the free-wheeling rotor cavitates from 1.6 m/s): phases open (or the
current-limited charger running), with only iron and friction loss in the motor (ASSUMED drag x runaway speed:
5.8 W at 2.5 m/s). The GM4108H reaches
69.9-82.8 V towed and needs >= 100 V parts or a clamp.

**A continuous phase short is not a tow brake.** With all three phases shorted, the rotor settles where the brake
torque Kt0^2 w R_ll / (R_ll^2 + X_ll^2) plus iron drag equals the rotor torque (scanning down from runaway), and the
windings then dissipate P_cu = Kt0^2 w^2 R_ll / (R_ll^2 + X_ll^2) (the brake power; DC-equivalent current
sqrt(P_cu / R_ll), phase RMS x sqrt(2/3)):

| Motor | Flow (m/s) | Equilibrium rpm (TSR) | Brake torque (Nm) | Winding copper loss (W) | Current DC-equiv / phase RMS (A) | vs rated copper loss | vs peak current |
|---|---|---|---|---|---|---|---|
| GM4108H (KV 27.0, R_ll 11.10) | 2.5 | 826 (3.36) | 1.0 | **84.3** | 2.76 / 2.25 | no rating | no rating |
| GM4108H (KV 27.0, R_ll 22.20) | 2.5 | 1163 (4.73) | 0.711 | **83.6** | 1.94 / 1.58 | no rating | no rating |
| GM4108H (KV 32.0, R_ll 11.10) | 2.5 | 984 (4.0) | 0.852 | **85.2** | 2.77 / 2.26 | no rating | no rating |
| GM4108H (KV 32.0, R_ll 22.20) | 2.5 | 1334 (5.42) | 0.585 | **78.3** | 1.88 / 1.53 | no rating | no rating |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 1.0 | 254 (2.58) | 0.156 | **3.5** | 0.88 / 0.72 | 0.3x | 0.27x |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 1.5 | 582 (3.94) | 0.299 | **16.7** | 1.93 / 1.57 | 1.41x | 0.58x |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 2.0 | 987 (5.01) | 0.414 | **40.2** | 2.99 / 2.44 | 3.4x | 0.91x |
| GL40 KV70 (KV 63.5, R_ll 4.50) | 2.5 | 1494 (6.07) | 0.47 | **69.6** | 3.93 / 3.21 | 5.89x | 1.19x |

At 2.5 m/s the GL40 would dissipate 69.6 W in its windings, 5.89x its rated copper loss, at
3.93 A (1.19x the 3.3 A peak rating), inside a sealed PETG hull. Even at 1.5 m/s it is
16.7 W (1.41x). A short is acceptable only as a parking brake at 1.0 m/s or less
(3.5 W, 0.3x). Pulsing it does not help at towing speed: a pulsed short lets the rotor run
between the short equilibrium and runaway, where the short-circuit current is 3.93-4.75 A, above the
3.3 A peak rating in every pulse. If the rotor must be slowed while towed (battery full), use a DC dump
load across the bridge output instead: 10.02 ohm holds the rotor at 1849 rpm (TSR 7.51) at 2.5 m/s with the motor at
its rated copper loss (11.62 W of 11.81 W) and puts 29.47 W into the resistor, outside the motor
(time-domain model; rate the resistor >= 2x, a design choice).

At 1.0 m/s the GL40's bridge optimum is 4.55 V, so the 5 V stage must be buck-boost. Phase current at
the GL40's 1.5 m/s optimum is 0.81 A RMS (FOC): a 3 x 0.5 mm2 tether is ample. Tether loss over an
ASSUMED 5 m of 0.5 mm2 (0.344 ohm loop) is 0.34 W at that point.

## 3. Bearings, shaft, axial stop and lip seal

| Part | d x D x B | Material | C / C0 (N) | Max rpm | Price | Source |
|---|---|---|---|---|---|---|
| S608-2RS | 8 x 22 x 7 | AISI 440C | 3296 / 1383 | 34,000 | USD 5.63 | [src](https://bearingsdirect.com/s608-2rs-stainless-steel-miniature-ball-bearing-8x22x7mm-sealed/) |
| 688-2RS | 8 x 16 x 5 | 52100 chrome (SS688-2RS stainless USD 5.63) | 1254 / 587 | 36,000 | USD 4.38 | [src](https://bearingsdirect.com/688-2rs-miniature-ball-bearing-8x16x5-sealed-mr688-2rs/) |
| 688-2RS (VXB) | 8 x 16 x 5 | 52100 | 1600 / 710 | 28,000 | USD 7.77 | [src](https://vxb.com/products/688-2rs-bearing-8x16x5-sealed-miniature) |

Load check with the BEM thrust (maximum over TSR up to runaway): 25.0 N at 1.5 m/s and 70.3 N towed at
2.5 m/s. Static ratings are 19.7x (S608) and 8.3x (688) the towed thrust, so the choice is set by space. Keep
**608-2RS** (cheap, robust, many printed-seat references). Use 688-2RS where the tail is tight: its 16 mm OD equals
the lip seal OD, so one 16 mm bore can seat both. SKF makes a stainless W 608-2RS1 ([src](https://www.skf.com/us/products/rolling-bearings/ball-bearings/deep-groove-ball-bearings/productid-W%20608-2RS1); page did not render).

**Axial stop.** The rotor is downstream, so thrust pulls the hub and shaft aft (+X). The shaft needs a stop forward of
the forward bearing, bearing on its inner race; the bearing outer races need shoulders on their aft side. Options:
- set-screw collar Ruland MSC-8-SS: bore 8, OD 16, width 8, 303 stainless ([listing](https://www.amazon.com/Ruland-MSC-8-SS-Collar-Stainless-Metric/dp/B0063KWIYC), medium: from the
  listing title, the vendor page refuses scripts). No machining; works on the hardened shaft.
- DIN 471 A8 circlip: s 0.8 (h11), free d3 7.4, groove d2 7.6, groove width min 0.9, a max 3.2; also in
  stainless X45CrMoV15 ([Ametric DIN 471](https://ametric.com/Images/document/RetainingRings-Metric.pdf)). Needs a groove: fine on plain rod, a grinding job on 56 HRC.
The opposite direction (handling, stowed) is taken by the hub or a spacer against the aft bearing's inner race.

Shaft: MISUMI PSSFJ8, SUS440C induction hardened to 56 HRC~ with g6 tolerance (-0.005/-0.014). It meets the lip-seal
hardness guideline; plain 303/304 rod does not, though it is fine for the display build. Medium confidence: from a search
extract of the MISUMI catalogue PDF ([src](https://us.misumi-ec.com/pdf/fa/2019/2019_US_0202_003.pdf)).

Lip seal, 8 mm shaft sizes found: TC 8x16x5 NBR ([src](https://www.amazon.com/8x16x5-Rubber-Covered-Double-Garter/dp/B00MW7GOD6)),
TC 8x16x8 ([src](https://avxseals.com/products/tc8x16x8-shaft-oil-seal-rubber-double-lip)), 8x16x7 SKF HMS4 ([src](https://intechbearing.com/products/8x16x7-mm-oil-seal)),
8x16x7 Freudenberg BA ([src](https://tinex-bell.si/en/product/15135/oil-seals-freudenberg-8x16x7-ba-nbr)), 8x14x4 NBR ([src](https://www.123bearing.com/seals-OA-8X14X4-NBR)).
Keep **TC 8x16x5**. Install per the vendor guideline ([src](https://www.oil-seal-stocks.com/all-oil-seals/oil-seal-tc-rotary-shaft-seal/)): shaft h11,
>= 45 HRC, Ra 0.2-0.8 um; housing H8, Ra 1.6-6.3 um. In a print, model the 16 mm bore like a bearing seat and calibrate
it on the tolerance coupon. Lip surface speed at 750 rpm is 0.314 m/s, so wear
is negligible. VERIFY the garter spring material, because a carbon spring rusts in water.

## 4. Equator O-ring and radial gland

Source table (Apple Rubber, ISO 3601 static radial gland, no back-up ring):
CS 2.62 -> radial gland depth 2.03-2.18, squeeze 14-25%, diametral clearance max 0.13, groove width 3.56 +-0.13,
groove radius 0.13-0.51, eccentricity max 0.08 ([src](https://www.applerubber.com/src/pdf/iso-3601-metric-size-o-rings.pdf)).
General guidance: squeeze 15-30%, fill 60-85%, stretch <= 5% ([src](https://www.globaloring.com/o-ring-groove-design/)).

The calc assumes 2.5% nominal stretch onto the groove bottom and printed diameters +-0.10 mm (ASSUMED, after
calibration). Worst cases combine: the O-ring CS tolerance; the ID tolerance in the right direction (the thinnest
section comes from the smallest free ID, which stretches most, on the largest groove bottom; the thickest from the
largest ID on the smallest groove); full eccentricity within the radial clearance, applied to the squeeze **and** to
the local fill; and the groove bottom width at its -0.13 tolerance.

**Joint geometry and the wall under the blade pockets.** The aft shell carries the spigot (male, groove on its OD);
the front shell carries the socket (female bore). The three blade pockets of `data/layout_rotor.json` run from x
31.75 to 108.25 and cross the equator, so over the socket's seal zone (x 49.495-57.5,
from `data/layout_internal.json`) the female wall under a pocket is hull radius - pocket depth - bore/2, and it must
stay >= MIN_WALL 1.2. The pocket is deepest at the aft end of the zone (3.02 mm at the equator face).
The layout study places the bayonet channel forward of the seal zone between the pockets, so it is not under them.
Each CS 2.62 size, with the layout study's groove form (45 deg flank, bottom 2.9, gland depth 2.105):

| O-ring | ID | bore D | spigot OD | female wall under pockets (min, x) | pocket reserve | wall with 6 mm placeholder | spigot ID (1.2 wall / layout lug rule) | D47 bay radial clearance | fill worst | passes |
|---|---|---|---|---|---|---|---|---|---|---|
| -136 | 50.47 | 55.94 | 55.74 | 6.51 (x 57.5) | 5.31 | 3.34 | 49.33 / 47.73 | 0.36 | 0.797 | no |
| -137 | 52.07 | 57.58 | 57.38 | 5.69 (x 57.5) | 4.49 | 2.52 | 50.97 / 49.37 | 1.18 | 0.797 | no |
| -138 | 53.64 | 59.19 | 58.99 | 4.89 (x 57.5) | 3.69 | 1.72 | 52.58 / 50.98 | 1.99 | 0.796 | no |
| -139 | 55.25 | 60.84 | 60.64 | 4.06 (x 57.5) | 2.86 | 0.89 | 54.23 / 52.63 | 2.82 | 0.796 | no |
| -140 | 56.82 | 62.45 | 62.25 | 3.26 (x 57.5) | 2.06 | 0.09 | 55.84 / 54.24 | 3.62 | 0.796 | yes |
| -141 | 58.42 | 64.09 | 63.89 | 2.44 (x 57.5) | 1.24 | -0.73 | 57.48 / 55.88 | 4.44 | 0.797 | **pick** |
| -142 | 59.99 | 65.7 | 65.5 | 1.63 (x 57.5) | 0.43 | -1.54 | 59.09 / 57.49 | 5.25 | 0.796 | no |
| -143 | 61.6 | 67.35 | 67.15 | 0.81 (x 57.5) | -0.39 | -2.36 | 60.74 / 59.14 | 6.07 | 0.796 | no |
| -144 | 63.17 | 68.96 | 68.76 | 0.00 (x 57.5) | -1.20 | -3.17 | 62.35 / 60.75 | 6.88 | 0.796 | no |
| -145 | 64.77 | 70.6 | 70.4 | -0.82 (x 57.5) | -2.02 | -3.99 | 63.99 / 62.39 | 7.70 | 0.796 | no |
| -146 | 66.34 | 72.21 | 72.01 | -1.62 (x 57.5) | -2.82 | -4.79 | 65.60 / 64.00 | 8.50 | 0.796 | no |

"Pocket reserve" is how much deeper the pockets could get before the wall drops below 1.2 mm. Rule: largest CS 2.62 size whose female socket wall over the seal zone x 49.495-57.5 stays >= MIN_WALL 1.2 under the layout_rotor pockets with 1.0 mm reserve on the pocket depth, spigot OD in 62-68, worst fill <= 0.85, worst squeeze <= 0.30, and >= 2.0 mm radial clearance from a D 47.0 bay to the spigot bore (layout_internal lug-zone rule).
The 1.0 mm reserve is a design choice: the rotor layout will re-derive the pockets after the BEM re-run with the new
hinge. -143 (the previous pick) leaves 0.81 mm of wall under the pockets and -2.36 mm under the
old 6 mm placeholder, so it is rejected; -142 passes today with only 0.43 mm reserve. Under the old 6 mm
placeholder only -138 and smaller would pass. Every size from -139 up leaves the D 47 bay >= 2 mm clear of the spigot bore.

**Groove printability.** The spigot prints tip-up, so its groove width runs along Z and a rectangular groove would
leave a 2 mm flat overhang on its tip-side flank (MAX_OVERHANG 45). The groove therefore has a radial aft flank and a
45 deg tip-side flank, as in the internal layout, and the fill is computed on the trapezoidal area (groove plus the
radial clearance band). Sweep for -141 at 0.10 radial clearance (rule: worst fill <= 0.85, worst squeeze <= 0.30, then
maximise the worst minimum squeeze, then the bottom width closest to the layout study's 2.9):

| bottom width | top width | gland depth | squeeze nom | squeeze worst | fill nom | fill worst (ecc) | passes rule |
|---|---|---|---|---|---|---|---|
| 2.6 | 4.53 | 2.03 | 0.216 | 0.107 to 0.318 | 0.718 | 0.912 | no |
| 2.6 | 4.565 | 2.065 | 0.203 | 0.093 to 0.305 | 0.703 | 0.890 | no |
| 2.6 | 4.605 | 2.105 | 0.187 | 0.077 to 0.290 | 0.686 | 0.866 | no |
| 2.6 | 4.64 | 2.14 | 0.174 | 0.063 to 0.277 | 0.671 | 0.847 | yes |
| 2.6 | 4.68 | 2.18 | 0.158 | 0.047 to 0.262 | 0.655 | 0.825 | yes |
| 2.8 | 4.73 | 2.03 | 0.216 | 0.107 to 0.318 | 0.681 | 0.861 | no |
| 2.8 | 4.765 | 2.065 | 0.203 | 0.093 to 0.305 | 0.666 | 0.841 | no |
| 2.8 | 4.805 | 2.105 | 0.187 | 0.077 to 0.290 | 0.650 | 0.819 | yes |
| 2.8 | 4.84 | 2.14 | 0.174 | 0.063 to 0.277 | 0.636 | 0.800 | yes |
| 2.8 | 4.88 | 2.18 | 0.158 | 0.047 to 0.262 | 0.622 | 0.780 | yes |
| 2.9 | 4.83 | 2.03 | 0.216 | 0.107 to 0.318 | 0.663 | 0.838 | no |
| 2.9 | 4.865 | 2.065 | 0.203 | 0.093 to 0.305 | 0.649 | 0.818 | no |
| 2.9 | 4.905 | 2.105 | 0.187 | 0.077 to 0.290 | 0.634 | 0.797 | **pick** |
| 2.9 | 4.94 | 2.14 | 0.174 | 0.063 to 0.277 | 0.620 | 0.779 | yes |
| 2.9 | 4.98 | 2.18 | 0.158 | 0.047 to 0.262 | 0.606 | 0.759 | yes |
| 3.0 | 4.93 | 2.03 | 0.216 | 0.107 to 0.318 | 0.647 | 0.816 | no |
| 3.0 | 4.965 | 2.065 | 0.203 | 0.093 to 0.305 | 0.633 | 0.797 | no |
| 3.0 | 5.005 | 2.105 | 0.187 | 0.077 to 0.290 | 0.618 | 0.776 | yes |
| 3.0 | 5.04 | 2.14 | 0.174 | 0.063 to 0.277 | 0.605 | 0.758 | yes |
| 3.0 | 5.08 | 2.18 | 0.158 | 0.047 to 0.262 | 0.591 | 0.739 | yes |
| 3.2 | 5.13 | 2.03 | 0.216 | 0.107 to 0.318 | 0.616 | 0.775 | no |
| 3.2 | 5.165 | 2.065 | 0.203 | 0.093 to 0.305 | 0.603 | 0.757 | no |
| 3.2 | 5.205 | 2.105 | 0.187 | 0.077 to 0.290 | 0.589 | 0.737 | yes |
| 3.2 | 5.24 | 2.14 | 0.174 | 0.063 to 0.277 | 0.577 | 0.721 | yes |
| 3.2 | 5.28 | 2.18 | 0.158 | 0.047 to 0.262 | 0.564 | 0.703 | yes |

Reference: the same size with the rectangular table groove (3.56) has a worst fill of 0.865; a slicer
rounding the 2.9 bottom down to 2.8 (14 layers) gives 0.819.

**Pick: ISO 3601-141 (58.42 x 2.62), EPDM 70** (Apple Rubber lists EPDM for water applications; use
silicone grease, not petroleum grease). The male groove sits on the spigot:

| Gland feature | Value (mm) |
|---|---|
| Groove bottom D (on spigot) | 59.88 |
| Groove widths | bottom 2.9, top 4.905; aft flank radial, tip-side flank 45 deg; corner radius 0.13-0.51 |
| Radial gland depth | 2.105 (table 2.03-2.18) |
| Groove depth in the spigot | 2.005 |
| Spigot OD | 63.89 (radial clearance 0.1) |
| Female bore D | 64.09 |
| Female wall: under the pockets (min) / plain | 2.44 at x 57.5 / 5.27 |
| Pocket reserve / largest uniform pocket depth keeping 1.2 mm | 1.24 / 4.07 |
| Spigot ID: 1.2 mm wall under the groove / layout lug-zone rule | 57.48 / 55.88 |
| D 47 bay radial clearance in the spigot bore (lug rule) | 4.44 |
| Squeeze nominal / worst | 0.187 / 0.077 to 0.290 |
| Fill nominal / worst without / with eccentricity | 0.634 / 0.747 / 0.797 |
| At 0.15 radial clearance: squeeze worst / fill worst | 0.057 to 0.308 / 0.824 |
| Stretch range | 0.0144 to 0.0357 |

The spigot OD of 63.89 lies in the 62-68 target. The worst minimum squeeze (7.7%) is the
price of printed tolerances plus full eccentricity: centre the halves with the bayonet and print the coupon. The
0.20 mm diametral clearance exceeds the 0.13 mm table limit, which exists to stop extrusion under pressure; at less
than 0.1 bar of submersion that does not apply, and a tighter printed fit would bind. The 1.78 section is too
sensitive to print error (worst squeeze 8% to 38%). Params currently say CS 2.0: change to 2.62.

At the Apple Rubber table groove (rectangular 3.56, 0.10 radial clearance), for reference:

| O-ring | radial clr | bore D | spigot OD | female wall under pockets | squeeze nom | squeeze worst | fill nom | fill worst (no ecc / ecc) |
|---|---|---|---|---|---|---|---|---|
| ISO 3601 / AS568 -037 (1.78 x 63.22) | 0.1 | 67.47 | 67.27 | 0.75 | 0.241 | 0.082 to 0.385 | 0.755 | 0.950 / 1.034 |
| ISO 3601 / AS568 -136 (2.62 x 50.47) | 0.1 | 55.94 | 55.74 | 6.51 | 0.187 | 0.078 to 0.290 | 0.703 | 0.822 / 0.865 |
| ISO 3601 / AS568 -137 (2.62 x 52.07) | 0.1 | 57.58 | 57.38 | 5.69 | 0.187 | 0.078 to 0.290 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -138 (2.62 x 53.64) | 0.1 | 59.19 | 58.99 | 4.89 | 0.187 | 0.078 to 0.290 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -139 (2.62 x 55.25) | 0.1 | 60.84 | 60.64 | 4.06 | 0.187 | 0.078 to 0.289 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -140 (2.62 x 56.82) | 0.1 | 62.45 | 62.25 | 3.26 | 0.187 | 0.078 to 0.289 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -141 (2.62 x 58.42) | 0.1 | 64.09 | 63.89 | 2.44 | 0.187 | 0.077 to 0.290 | 0.703 | 0.821 / 0.865 |
| ISO 3601 / AS568 -142 (2.62 x 59.99) | 0.1 | 65.7 | 65.5 | 1.63 | 0.187 | 0.078 to 0.290 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -143 (2.62 x 61.6) | 0.1 | 67.35 | 67.15 | 0.81 | 0.187 | 0.078 to 0.290 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -144 (2.62 x 63.17) | 0.1 | 68.96 | 68.76 | 0.00 | 0.187 | 0.078 to 0.289 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -145 (2.62 x 64.77) | 0.1 | 70.6 | 70.4 | -0.82 | 0.187 | 0.078 to 0.289 | 0.703 | 0.821 / 0.864 |
| ISO 3601 / AS568 -146 (2.62 x 66.34) | 0.1 | 72.21 | 72.01 | -1.62 | 0.187 | 0.078 to 0.289 | 0.703 | 0.820 / 0.863 |
| ISO 3601 / AS568 -230 (3.53 x 63.09) | 0.1 | 70.41 | 70.21 | -0.72 | 0.178 | 0.091 to 0.259 | 0.729 | 0.830 / 0.861 |

## 5. Small parts

| Part | Spec | Key data | Source | Conf. |
|---|---|---|---|---|
| Magnet N52 4x2 | 4 +-0.05 x 2, Ni-Cu-Ni, Br 1.43-1.48 T, 80 C max | pull 4.05-4.58 N to steel; USD 0.30 each, minimum order 50 | [src](https://www.jc-magnetics.com/Magnet-N52-4mm-2mm-Discs), [src](https://themagnetbaron.com/collections/4mm-disc-magnets) | medium |
| Magnet N52 3x2 / 5x2 | alternatives | 2.54 N / 5.52 N | [src](https://themagnetbaron.com/collections/3mm-disc-magnets), [src](https://themagnetbaron.com/collections/5mm-disc-magnets) | medium |
| Dowel pin 3 x 22 | ISO 8734, type A or C1 stainless | d 3 m6, chamfer c 0.5, l 21.5-22.5 (22 is in the Table 1 length series; d=3 commercial lengths 8-30); 3 x 16 if the hinge shortens | [ISO 8734](https://cdn.standards.iteh.ai/samples/20002/ba95529de6a64437ac7454e14dbd2ef6/ISO-8734-1997.pdf), [Table 1](https://cdn.standards.iteh.ai/samples/20002/ba95529de6a64437ac7454e14dbd2ef6/ISO-8734-1997.pdf) | high |
| Heat-set insert M3 | CNC Kitchen standard M3 x 5.7 | L 5.7, D1 4.6, D2 3.9, hole 4.0, min wall 1.6; EUR 12.09/100 | [src](https://www.3djake.com/cnc-kitchen/threaded-inserts-m3-standard) | high |
| Grub screw M3 | ISO 4029 cup point, A2/A4 | pitch 0.5, hex 1.5, cup d <= 1.4, socket depth >= 1.2; M3x3, M3x4 sold | [src](https://www.fasteners.eu/standards/iso/4029/) | medium |
| Shaft collar 8 mm | Ruland MSC-8-SS | bore 8, OD 16, width 8, 303 stainless | [listing](https://www.amazon.com/Ruland-MSC-8-SS-Collar-Stainless-Metric/dp/B0063KWIYC) | medium |
| Circlip 8 mm | DIN 471 A8 | s 0.8, d3 7.4, groove 7.6 x 0.9 min | [Ametric](https://ametric.com/Images/document/RetainingRings-Metric.pdf) | high |
| Cable gland | Hummel HSK-K M12x1.5 1.209.1200.50 | clamp 3-6.5, GL 8, SW 15, H 21, 1.5 Nm, IP68 10 bar/30 min with extra O-ring (IP54 without); no locknut supplied | [src](https://www.hummel.com/en/cable-glands/products/plastic-cable-glands/hsk-k/1209120050-hsk-k-grey-m12x1-5/) | high |
| Gland locknut | Hummel 1.261.1200.50 (PA, grey RAL 7035, M12x1.5) | EUR 0.06; alt. 1.262.1200.50: 5 thick, SW 17 (search extract, low) | [src](https://actif-shop.eu/en/welcome/5460-hummel-polyamide-cable-gland-nut-m12-x-15-grey-ral-7035-1261120050.html) | medium |
| Gland alternative | LAPP SKINTOP ST-M 12x1.5 (53111000) | clamp 3.5-7, GL 8, SW 15, L 30, IP68 5 bar | [src](https://products.lappgroup.com/online-catalogue/cable-glands/skintop-cable-glands-plastic-metric/standard/skintop-st-m.html) | high |
| Gland hole | M12x1.5 clearance | 12.5 drill, or print/tap the thread | [src](https://www.metric-bolt-torque.com/metric-thread/cable-gland-sizes/) | medium |
| Eye bolt | DIN 580 | smallest in standard M8 (an A4 M6 exists with the same eye): d3 36, h 36, l 13 | [src](https://www.fasteners.eu/standards/DIN/580/), [src](https://www.jwwinco.com/en-us/products/3.4-Installing-Lifting-Damping-with-Leveling-Feet-Lifting-Gear-and-Rubber-Bumpers/Lifting-eye-bolts-Lifting-eye-nuts/DIN-580-A4-Stainless-Steel-Lifting-Eye-Bolts) | high |

Notes:
- Magnet pull is quoted magnet-to-thick-steel at contact. It drops steeply behind a printed skin, so set detent
  strength on the coupon. The 3x2 and 5x2 sizes give a 0.55x to 1.2x range around the 4x2.
- ISO 8734 specifies m6. One vendor guide calls it h6, but the difference is micrometres and irrelevant against
  printed holes. Hinge: pin pressed into the hub (PRESS_FIT), blade turning on the pin (CLEAR_ROTATE).
- The gland thread is 8 mm long against a 2.4 mm hull wall: either a locknut inside the nose (fitted before the
  shells close) or a printed boss at least 8 mm deep carrying the thread.
- A DIN 580 eye is 36 mm tall on a 75 mm nose. Use a printed eye integral to the nose and pull-test it to >= 3 x the
  BEM maximum rotor thrust at 2.5 m/s (211 N) plus the hull drag, which this study does not compute.

## 6. Printing: build volume, fits, materials

Build volume 256 x 256 x 256 mm on all target printers: X1C ([spec PDF](https://public-cdn.bambulab.com/store/bambulab-X1-carbon-tech-specs.pdf)),
P1P ([spec PDF](https://public-cdn.bambulab.com/store/bambulab-P1P-tech-specs.pdf)), A1 and P1S (Bambu tech-spec pages via search extract, medium).
P1S and X1C are enclosed; X1C lists ASA as supported and P1P as "capable".

Fits:
- Prusa: moving parts need "at least 0.3 mm"; printer accuracy is about 0.2 mm; overhangs 45-60 deg
  ([src](https://help.prusa3d.com/article/modeling-with-3d-printing-in-mind_164135)).
- Creative3DP fit ladder, diametral (PLA/PETG/ABS, 5-25 mm): press -0.10, snug +0.05, close running +0.15, free +0.35.
  Printed holes come out 0.1-0.3 undersize. Press interference is 0.05-0.10 per side for PLA and up to 0.10 for PETG.
  A 608 pocket for PLA is modelled at 21.95 with standard hole compensation; PETG sliding fits need +0.3 total; add a
  0.5 x 45 deg lead-in chamfer ([src](https://tools.creative3dp.com/blog/press-fit-tolerances-3d-printing/)). Medium confidence:
  a tool vendor's blog that cites Prusa, Markforged and a 2023 Processes paper.
- Crush ribs (small deformable standoffs) make press fits tolerant of +-0.1 mm ([src](https://hackaday.com/2020/10/15/adding-crush-ribs-to-3d-printed-parts-for-a-better-press-fit/)).

Against params: CLEAR_SLIDE 0.25 radial (0.5 diametral) and CLEAR_ROTATE 0.30 radial sit at or above the sources, which is
sensible for PETG in silty water. PRESS_FIT -0.05 matches the sources only if read **per side** (-0.10 diametral). Make
that explicit and add 3 crush ribs to the bearing, seal and pin seats.

Materials ([Prusa PETG](https://help.prusa3d.com/article/petg_2059), [Prusa ASA](https://help.prusa3d.com/article/asa_1809)):
PETG is "water and humidity resistant" and suits most exterior use below 80 C. ASA is "suitable for outdoor use" for its
UV and temperature resistance (to 93 C), but needs an enclosure, a bed of >= 100 C, ventilation (styrene) and care with
warping. Use PETG for every wetted part; use ASA only for long sun exposure and only on P1S/X1C.

## 7. Prior art

| Product | Rotor D | Rated power | Flow speed | Mass | Size | Price | Stowage | Source |
|---|---|---|---|---|---|---|---|---|
| Enomad Estream | 210 mm (deployed) | 2.5-5 W | not published | 800 g | stowed 245 x 65 mm | USD 180 (Kickstarter, 2016) | folding blades; turbine and battery modules | [src](https://newatlas.com/enomad-estream-portable-water-turbine-charger/44897/), [src](https://www.greenoptimistic.com/estream-portable-water-turbine-20160816/) |
| WaterLily (discontinued) | n/a (unit 178 mm dia x 76 mm) | up to 15 W; 2-3 W under tow | min 1 km/h; peak 7.2 mph (3.2 m/s) | 1.3 kg | 7 x 3 in | USD 159 | rigid ducted wheel, tethered or towed | [src](https://get.waterlilyturbine.com/paddling-charger/), [src](https://www.popsci.com/gear/best-hydroelectric-generators/), [src](https://worldwaterreserve.com/waterlily-turbine-review/), [src](https://www.theoutbound.com/neonconflagrations/waterlily-turbine-portable-power-generator-review) |
| Blue Freedom Portable | 120 mm | max 5 W (15 V nominal) | not published by the maker; a translated review says 'at least 1 m/s' | 690 g | body 220 x 55 mm | EUR 299 | prop stored in body, flexible shaft | [src](https://blue-freedom.net/hydro/portable/), [review](https://gigazine.net/gsc_news/en/20160418-blue-freedom/) |
| Ampair Aquair UW | 318 mm (12.5 in) | 70 W @ 6 mph, 100 W @ 9 mph | 2.7 / 4.0 m/s | 13.6 kg | 14.5 in long | not for sale | fixed prop, pole or tow | [src](https://nooutage.com/aquair-uw.htm) |
| **Ovo (this project)** | 194 mm deployed | 18.3 W shaft @ 1.5 m/s (BEM, Cp peak); 7.66 W DC with the GL40 via a diode bridge (8.92 W active) | 1.5 m/s design | n/a | stowed 75 x 165 mm | display build priced parts EUR 35.48 to buy | blades fold forward into tail pockets | `hw_calc.py`, `data/blade_design.json` |

Power sanity check, P = 0.5 rho A V^3 Cp (computed in `hw_calc.py`; Cp_el = electrical output over kinetic flux):

| Product | Rotor D (m) | P (W) | V (m/s) | P_avail (W) | Cp implied | V needed @ Cp_el 0.20 | V needed @ Cp_el 0.30 |
|---|---|---|---|---|---|---|---|
| Enomad Estream | 0.21 | 5.0 | n/a | n/a | n/a | 1.13 | 0.99 |
| Blue Freedom Portable (5 W max; no rating speed published) | 0.12 | 5.0 | n/a | n/a | n/a | 1.64 | 1.43 |
| WaterLily (7 in housing as upper bound on rotor) | 0.1778 | 15.0 | 3.219 | 414.0 | 0.036 | 1.82 | 1.59 |
| Aquair UW @6 mph | 0.3175 | 70.0 | 2.682 | 763.9 | 0.092 | 2.07 | 1.81 |
| Aquair UW @9 mph | 0.3175 | 100.0 | 4.023 | 2578.2 | 0.039 | 2.33 | 2.03 |
| Ovo BEM shaft (Cp peak, 1.5 m/s) | 0.194 | 18.31 | 1.5 | 49.9 | 0.367 | 1.84 | 1.6 |

Reading:
- Blue Freedom states "max. 5 W" and no flow speed. A machine-translated review says at least 1 m/s is required,
  which reads as a minimum operating speed. At a plausible Cp_el of 0.20-0.30 a 120 mm rotor needs
  1.43-1.64 m/s for 5 W.
- The Estream (210 mm, 5 W) needs ~0.99-1.13 m/s at the same Cp_el; its rating speed is unpublished.
- WaterLily's 15 W peak sits at 3.2 m/s, i.e. towed-boat speeds.
- Ovo's 194 mm rotor gives 18.3 W shaft at 1.5 m/s (BEM Cp peak) and 7.66 W DC through a diode bridge with
  the GL40 (Cp_el 0.153). Per swept area it is in the Estream's class (power scales with V^3, and the Estream's
  rating speed is unpublished). Its distinctive points are stowage and openness: blades fold flush into a
  0.401 L egg, against 0.813 L for the Estream's 245 x 65 mm tube, and the design is fully
  parametric and documented.
- For reference, a full-scale 1.58 m river HAHT designed for 1 kW at 1.5 m/s measured Cp 0.4087 ([Chica et al. 2025](https://doi.org/10.3390/en18071657)).
  The constrained-BEM Cp of 0.367 is consistent with that.

## 8. BOM

| Item | Qty | Spec | EUR per-piece x qty | EUR to buy (pack / MOQ) | Purchase basis | Display | v2 only |
|---|---|---|---|---|---|---|---|
| Deep groove ball bearing S608-2RS | 2 | 8x22x7, AISI 440C stainless, 2RS | 10.03 | 10.03 | sold singly, USD 5.63 each, no minimum (bearingsdirect) | yes |  |
| Shaft 8 mm | 1 | MISUMI PSSFJ8 SUS440C hardened 56 HRC, g6; length set by the layout study | n/a | n/a | not priced | yes |  |
| Shaft collar 8 mm, set screw, stainless (Ruland MSC-8-SS) | 1 | bore 8, OD 16, width 8, 303 stainless; axial stop forward of the forward bearing (rotor thrust pulls aft) | n/a | n/a | not priced | yes |  |
| Dowel pin ISO 8734 3x22 | 3 | d 3 m6, l 22 (21.5-22.5), type A or C1 stainless; one per blade hinge; length from data/layout_rotor.json (ear-knuckle-ear) | n/a | n/a | not priced | yes |  |
| Magnet N52 disc 4x2 | 16 | 4 +-0.05 x 2 mm, Ni-Cu-Ni; 12 blade detents (layout_rotor: 4 per blade) + 4 bayonet detent (layout_internal: 2 pairs) | 4.28 | 13.36 | JC Magnetics minimum order 50 at USD 0.30 each (USD 0.17 from 100) | yes |  |
| Heat-set insert M3 x 5.7 (CNC Kitchen standard) | 8 | hole 4.0, min wall 1.6; qty ESTIMATE (spinner, fins, plate) | 0.97 | 12.09 | pack of 100 (3DJake EUR 12.09; cnckitchen.store EUR 9.40, ruthex EUR 8.99) | yes |  |
| Socket head cap screw M3 (ISO 4762), A2 | 8 | lengths set by CAD; qty ESTIMATE | n/a | n/a | not priced | yes |  |
| Grub screw M3 x 4 cup point (ISO 4029), A2/A4 | 2 | hex 1.5; hub on shaft flat, generator hub | n/a | n/a | not priced | yes |  |
| Filament PETG | 1 | all wetted parts; mass set by CAD | n/a | n/a | not priced | yes |  |
| Printed tow eye (nose) | 1 | integral to the nose shell; pull test >= 211 N (3 x BEM rotor thrust at 2.5 m/s) plus hull drag | 0.00 | 0.00 | printed | yes |  |
| Generator CubeMars GL40 KV70 | 1 | D46.5 x 21.5, 14 pole pairs, KV 63.5, 4.5 ohm ph-ph, 1.8 mH, 8 mm through bore; stator 4-M3 PCD 24 + 4-M2.5 PCD 40 (same axes), pad block at 45 deg; rotor 4-M3 PCD 24; rated 0.25 Nm / 1.62 A | 65.02 | 65.02 | sold singly, USD 72.99 |  | yes |
| M3 screws for the generator (ISO 4762 / ISO 7380), A2 | 8 | 4 stator (thread 3.5 deep) + 4 rotor hub (thread 2.5 deep); length = clamped part + 3.0 / 2.0 max | n/a | n/a | not priced |  | yes |
| Lip seal TC 8x16x5 NBR | 1 | double lip + garter spring; housing bore 16 | n/a | n/a | not priced |  | yes |
| O-ring ISO 3601-141 (AS568-141) | 2 | 58.42 x 2.62, EPDM 70 or NBR 70 (one spare) | n/a | n/a | not priced |  | yes |
| Cable gland Hummel HSK-K M12x1.5 grey + additional thread O-ring | 1 | 1.209.1200.50, clamp 3-6.5, IP68 with O-ring | 0.93 | 9.34 | pack of 10, USD 10.48 (tameson) |  | yes |
| Locknut M12x1.5 PA (Hummel 1.261.1200.50) | 1 | gland ships without a locknut; or a printed boss >= 8 mm deep | 0.06 | 0.06 | per piece, EUR 0.06 excl. tax, no minimum stated (actif-shop) |  | yes |
| 3-core cable 3 x 0.5 mm2 | 1 | OD 3-6.5 mm; length ASSUMED 5 m | n/a | n/a | not priced |  | yes |
| 3-phase rectifier + buck-boost to 5 V USB | 1 | input >= 60 V with the GL40 (open-circuit freewheel 35.2 V at 2.5 m/s, the towing state); >= 100 V if the GM4108H is used (82.8 V); buck-boost because the 1.0 m/s optimum is ~4.55 V; current limit so the winding copper loss stays at the rated 11.81 W (DC side 1.62 A block model / 1.715 A time-domain). No continuous phase-short brake (it would put 69.6 W into the windings at 2.5 m/s). Optional dump load 10.02 ohm across the DC link (29.47 W at 2.5 m/s; rate >= 2x) | n/a | n/a | not priced |  | yes |

Two totals, both for the priced lines only. Pro-rata per-piece cost: display build EUR 15.28, v2 additions
EUR 66.01. Cash outlay to buy them (minimum orders and packs: 50 magnets, 100 inserts, 10 glands):
display build EUR 35.48, v2 additions EUR 74.42.
11 items have no sourced price yet (null in `bom.json`). USD converted at the ECB reference rate of
1.1225 USD/EUR (2026-10-02).

## 9. Open issues

- GL40: no published no-load loss (iron drag ASSUMED); which part carries the 8 mm bore is inferred from the drawing
  (both faces show an 8 mm bore). Check both on receipt. Hole clocking and the pad block are now measured from the
  drawing's vector paths.
- GM4108H (alternative only): KV, resistance convention and the 1200-1800 g.cm load point are mutually inconsistent;
  inductance unpublished (time-domain numbers use an ESTIMATE); base pattern from one community source, possibly
  another revision; rotor face unconfirmed. Bench test first.
- The O-ring wall check uses the current `data/layout_rotor.json` pockets (sha256 04636f4af337); the pick keeps a
  1.0 mm reserve because those pockets will be re-derived after the BEM re-run. Re-run `hw_calc.py` when either layout file changes.
- The generator-optimal TSR at 1.5 m/s (5.29-6.25 for the GL40 and GM4108H) is above the rotor's Cp peak (3.822); the BEM study
  may want to evaluate a design TSR of 5-5.5.
- Hull drag is not in the tow-eye load; lip-seal garter spring material for water, and prices for seals, O-rings,
  shaft, pins, screws, collar and cable, still to be sourced.
- Ruland collar dimensions come from listing titles (the vendor pages refuse scripts); Hummel 1.262 locknut size from a
  search extract.
- A1/P1S build volume from a search extract (the Bambu tech-spec pages return 403 to scripts).

## 10. Changes in this revision

Second fix round (response to the second adversarial review):
1. Equator O-ring: the female socket wall is now checked under the real blade pockets (`data/layout_rotor.json`) over
   the seal zone. -143 failed (0.81 mm wall); the pick is -141 with 1.24 mm pocket reserve. Spigot bore and the
   bay clearance are re-derived from it (4.44 mm around a D 47 bay). The hull fit and the bay window use the
   layout_rotor pocket keep-out instead of the 6 mm placeholder (kept as a sensitivity).
2. GL40 stator face: the PCD 40 M2.5 holes are on the same axes as the PCD 24 M3 holes (0/90/180/270), measured from
   the drawing's vector paths; the solder-pad block at 45 deg gets a keep-out and a window in the adapter plate and replica.
3. Phase-short brake: its winding dissipation is now computed and it is no longer recommended for towing
   (69.6 W, 5.89x the rated copper loss at 2.5 m/s). Towing state: open conflict (superseded in part, section 2.7); phases open was the first recommendation, optional DC dump load.
4. Coupled scan capped at the published continuous rating; the GL40 2.5 m/s row is flagged and a capped value given.
5. Inductance: a time-domain bridge model (R + L, capacitor link) replaces the "reactance is small" assumption for every
   motor with a published L. GL40 at 1.5 m/s: 7.66 W (block model 7.78 W). The review's 6.9 W treats the overlap
   drop as dissipated with a stiff DC current; for a capacitor-fed bridge that is pessimistic (table in 2.1). GM4108H
   numbers are flagged optimistic until its L is measured.
6. Groove: 45 deg tip-side flank (printable tip-up), fill on the trapezoidal area, bottom width aligned with the
   internal layout (2.9).
7. BOM: purchase cost (minimum orders and pack sizes) next to the pro-rata cost; magnet count 16 (blade and bayonet
   detents); pins 3 x 22 as the rotor layout needs.

First fix round (summary): freewheel voltage at the BEM runaway speed; GM4108H load-point contradiction analysed and
its mount provenance downgraded; O-ring worst case corrected; iron loss and the 3/pi bridge factor added; Blue
Freedom Betz claim dropped; recommendation moved from the GM4108H to the GL40.
