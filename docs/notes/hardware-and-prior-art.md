# Hardware selection, v2 generator trade study and prior art

Date: 2026-10-03. Data: `data/hardware.json` (every value with `source_url` and `confidence`), `data/bom.json`.
Scripts (re-runnable): `out/work/hardware/hw_calc.py` (losses, gland, fit, prior-art check),
`out/work/hardware/step_mount.py` (community motor-mount STEP measurement), `out/work/hardware/build_outputs.py`,
`out/work/hardware/write_notes.py` (this file). Numbers below are generated from those outputs.
Confidence: high = manufacturer datasheet or standard; medium = dealer, community or derived; low = inference.

## 1. Decisions at a glance

| Item | Choice | Key dimensions | Confidence |
|---|---|---|---|
| v2 generator | iPower GM4108H-120T | D 47.0 x L 22.7 (bay reserve 33), hollow shaft OD 10 / ID 8, base 4 x M3 on PCD 38.14 | medium (interface), see 2.5 |
| Bearings | S608-2RS (440C), 688-2RS as compact alternative | 8x22x7; 8x16x5 | high |
| Shaft | 8 mm SUS440C hardened (MISUMI PSSFJ8), g6 | 56 HRC~, g6 = -0.005/-0.014 | medium |
| Lip seal (v2) | TC 8x16x5 NBR | housing bore 16 (same as 688) | medium |
| Equator O-ring | ISO 3601-143 (AS568-143), EPDM 70 | ID 61.60 +-0.51 x CS 2.62 +-0.08 | high |
| Detent magnets | N52 disc 4x2 | pull 4.0-4.6 N to steel | medium |
| Hinge pins | ISO 8734 3 x 16 (A or C1) | d 3 m6, c 0.5 | high |
| Inserts | CNC Kitchen M3 x 5.7 standard | hole 4.0, min wall 1.6 | high |
| Cable gland (v2) | Hummel HSK-K M12x1.5 + thread O-ring | clamp 3-6.5, GL 8, SW 15, H 21 | high |
| Tow eye | printed, integral to the nose | DIN 580 eye is 36 mm tall | high (DIN data) |
| Printer | any Bambu A1 / P1S / X1C | 256 x 256 x 256 | high/medium |
| Material | PETG wetted parts; ASA only on enclosed printers | | high |

## 2. v2 generator trade study

### 2.1 Design point and conventions

Design band 550-750 rpm and 0.20-0.35 Nm; nominal 650 rpm, 0.27 Nm = 18.38 W shaft.
The provisional BEM run (salvaged, 1.5 m/s) gives TSR 4.0: 591 rpm / 0.296 Nm; TSR 4.5: 665 / 0.259; TSR 5.0: 738 / 0.227,
inside the band. Low flow (0.8 m/s, TSR 4): 315.0 rpm. Towed freewheel (2.5 m/s, TSR ~6, ESTIMATE):
1476.7 rpm. Generator envelope: D <= 56 mm, L <= 35 mm.

KV is the hobby convention: no-load rpm per volt of DC bus, i.e. per volt of **line-to-line peak** back-EMF.
Define Kt0 = 60 / (2 pi KV) [Nm/A]. With phase flux linkage lambda (V s/rad, mechanical), the ll-peak EMF is
sqrt3 lambda w, so Kt0 = sqrt3 lambda. For sinusoidal (FOC) currents of peak I_pk in phase with the EMF:
tau = 1.5 lambda I_pk and P_cu = 1.5 I_pk^2 R_ph. Eliminating I_pk:

    P_cu = 2 R_ph tau^2 / (3 lambda^2) = R_ll (tau / Kt0)^2,   R_ll = 2 R_ph (star)

The same expression holds for six-step / diode-bridge operation, where the DC current I = tau / Kt0 flows through two
phases. So the study uses **I = tau / Kt0 (DC-equivalent) and line-to-line resistance**; the phase RMS current is
I sqrt(2/3). Km = Kt0 / sqrt(R_ll) is the figure of merit: P_cu / P_shaft = tau / (Km^2 w). Rectifier loss assumes two
Schottky drops of 0.4 V (ASSUMED). Iron and friction loss are estimated only where a no-load current is
published. Inductance is ignored. Where it is published, the line-to-line reactance at 650 rpm is 1.72 ohm (GL40 KV70),
1.26 ohm (GL40 II) and 0.71 ohm (RI50), small against R_ll.

Where a vendor gives "internal resistance" without a convention, both readings are carried (value = R_ll, or value = R_phase
so R_ll = 2 x value). KV for the GM4108H is not published. 27 rpm/V follows from the spec's no-load 513-567 rpm at 20 V,
and ~32 rpm/V is the value a SimpleFOC user runs. The two differ by the sine-PWM factor sqrt3/2.

Km needed for a copper-only efficiency target at the nominal point (`Km_required`): 60%: 0.0996,
70%: 0.115, 80%: 0.1408, 90%: 0.1992 Nm/sqrt(W).

### 2.2 Candidates

| Candidate | Class | D x L (mm) | KV (rpm/V) | R published (ohm) | Mass (g) | Price | In envelope | Source |
|---|---|---|---|---|---|---|---|---|
| iPower GM4108H-120T | gimbal outrunner | 47.0 x 22.7 | 27.0/32.0 | 11.100 | 124 | USD 35.9, EUR 43.9 | yes | [src](https://shop.iflight.com/ipower-motor-gm4108h-120t-brushless-gimbal-motor-pro217) |
| Rctimer GBM4114HS 100T | gimbal outrunner | 46.0 x 31.0 | 32.0 | 15.200 | 150 | USD 34.99 | yes | [src](https://www.rctimer.com/4114-hollow-shaft-brushless-gimbal-motor-p0440.html) |
| CubeMars GL40 KV70 | gimbal outrunner | 46.5 x 21.5 | 63.5 | 4.500 | 107 | USD 72.99 | yes | [src](https://www.cubemars.com/data/cms/202602/gl40-gimbal-motor-2d-drawing.pdf) |
| CubeMars GL40 II KV82.5 | gimbal outrunner (integrated driver) | 46.1 x 33.5 | 82.5 | 3.000 | 125 | USD 133.99 | yes | [src](https://cubemars.com/goods-1182-GL40+II.html) |
| T-motor GB4106 KV53 | gimbal outrunner | 47.8 x 20.0 | 53.0 | 15.300 | 70 | USD 45.9 | yes | [src](https://store.tmotor.com/product/gb4106-gimbal-type.html) |
| EMAX GB4008 66KV | gimbal outrunner | n/a x n/a | 66.0 | 7.500 | 101 | USD 23.49 | unknown | [src](https://shop.emaxmodel.com/products/emax-gb4008-brushless-gimbal-motor-hollow-shaft-66kv-90t) |
| CubeMars RI50 KV100 (frameless) | frameless inrunner torque motor | 54.0 x 27.0 | 100.0 | 1.420 | 181 | USD 61.9 | yes | [src](https://www.cubemars.com/product/ri50-kv100-frameless-inrunner-torque-motor.html) |
| T-motor MN5008 KV170 | low-KV drone outrunner | 55.6 x 32.0 | 170.0 | 0.270 | 128 | USD 89.99 | yes | [src](https://store.tmotor.com/product/mn5008-kv170-motor-antigravity-type.html) |
| T-motor MN4006 KV380 | low-KV drone outrunner | 44.4 x 21.0 | 380.0 | 0.194 | 68 | USD 149.9 | yes | [src](https://store.tmotor.com/product/mn4006-kv380-motor-antigravity-type.html) |
| maxon EC 45 flat 70 W 36 V (411815) | industrial flat outrunner | 45.0 x n/a | 179.0 | 1.740 | 141 | EUR 226.94 | unknown | [src](https://www.maxongroup.com/maxon/view/product/motor/ecmotor/ecflat/ecflat45/411815) |
| iPower GM5208-24 (outside envelope) | gimbal outrunner | 59.5 x 31.2 | 20.8 | 13.700 | 204 | USD 51.0 | no | [src](https://iflighteurope.com/products/ipower-gm5208-24-gimbal-motor) |
| T-motor GB54-2 KV26 (outside envelope) | gimbal outrunner | 60.7 x 25.0 | 26.0 | 15.000 | 156 | USD 74.9 | no | [src](https://store.tmotor.com/goods-445-T-Motor+GB54-2++Gimbal+Type+3-6S+UAV+Drone+Motor+26KV.html) |

Categories covered: gimbal outrunners (GM4108H, GBM4114HS, GL40, GL40 II, GB4106, GB4008), low-KV drone outrunners
(MN5008, MN4006), a frameless torque motor (RI50) and an industrial flat BLDC (EC 45 flat). For small hub or axial-flux
generators, no off-the-shelf unit with a datasheet exists within D 56 mm; search hits start at 50 W to 20 kW industrial
units ([search hit](https://www.permanent-magnetalternator.com/sale-13975594-low-rpm-axial-flux-coreless-permanent-magnet-generator-50w-20kw.html)). A custom PCB-stator axial-flux generator is a v3 idea, not a part.
Two out-of-envelope motors (GM5208-24, GB54-2) are kept as references: they show what a 60 mm hull would buy.

### 2.3 Losses at the nominal point (650 rpm, 0.27 Nm)

| Candidate | KV | R_ll | Km (Nm/sqrtW) | Voc ll-pk / ll-rms (V) | I_dc / I_ph,rms (A) | P_cu (W) | P_rect (W) | eta_cu | eta_cu+rect | eta_cu over 550-750 rpm, 0.20-0.35 Nm | Pmax matched (W) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| iPower GM4108H-120T (iFlight) | 27.0 | 11.10 | 0.1062 | 24.1 / 17.0 | 0.76 / 0.62 | 6.5 | 0.61 | 0.65 | 0.61 | 0.46 to 0.77 | 13.1 |
| iPower GM4108H-120T (iFlight) | 27.0 | 22.20 | 0.0751 | 24.1 / 17.0 | 0.76 / 0.62 | 12.9 | 0.61 | 0.30 | 0.26 | -0.08 to 0.55 | 6.5 |
| iPower GM4108H-120T (iFlight) | 32.0 | 11.10 | 0.0896 | 20.3 / 14.4 | 0.91 / 0.74 | 9.1 | 0.72 | 0.51 | 0.47 | 0.24 to 0.68 | 9.3 |
| iPower GM4108H-120T (iFlight) | 32.0 | 22.20 | 0.0633 | 20.3 / 14.4 | 0.91 / 0.74 | 18.2 | 0.72 | 0.01 | -0.03 | -0.52 to 0.36 | 4.7 |
| Rctimer GBM4114HS 100T | 32.0 | 15.20 | 0.0765 | 20.3 / 14.4 | 0.91 / 0.74 | 12.4 | 0.72 | 0.32 | 0.28 | -0.04 to 0.56 | 6.8 |
| Rctimer GBM4114HS 100T | 32.0 | 30.40 | 0.0541 | 20.3 / 14.4 | 0.91 / 0.74 | 24.9 | 0.72 | -0.35 | -0.39 | -1.07 to 0.13 | 3.4 |
| CubeMars GL40 KV70 | 63.5 | 4.50 | 0.0709 | 10.2 / 7.2 | 1.79 / 1.47 | 14.5 | 1.44 | 0.21 | 0.13 | -0.21 to 0.49 | 5.8 |
| CubeMars GL40 II KV82.5 (integrated driver) | 82.5 | 3.00 | 0.0668 | 7.9 / 5.6 | 2.33 / 1.91 | 16.3 | 1.87 | 0.11 | 0.01 | -0.36 to 0.43 | 5.2 |
| CubeMars GL40 II KV82.5 (integrated driver) | 86.8 | 3.00 | 0.0635 | 7.5 / 5.3 | 2.46 / 2.00 | 18.1 | 1.96 | 0.02 | -0.09 | -0.51 to 0.37 | 4.7 |
| T-motor GB4106 KV53 | 53.0 | 15.30 | 0.0461 | 12.3 / 8.7 | 1.50 / 1.22 | 34.4 | 1.20 | -0.87 | -0.94 | -1.86 to -0.20 | 2.5 |
| EMAX GB4008 66KV | 66.0 | 7.50 | 0.0528 | 9.8 / 7.0 | 1.87 / 1.52 | 26.1 | 1.49 | -0.42 | -0.50 | -1.18 to 0.09 | 3.2 |
| CubeMars RI50 KV100 (frameless inrunner) | 100.0 | 1.42 | 0.0801 | 6.5 / 4.6 | 2.83 / 2.31 | 11.3 | 2.26 | 0.38 | 0.26 | 0.05 to 0.60 | 7.4 |
| CubeMars RI50 KV100 (frameless inrunner) | 79.6 | 1.42 | 0.1007 | 8.2 / 5.8 | 2.25 / 1.84 | 7.2 | 1.80 | 0.61 | 0.51 | 0.40 to 0.75 | 11.8 |
| T-motor MN5008 KV170 | 170.0 | 0.27 | 0.1081 | 3.8 / 2.7 | 4.81 / 3.92 | 6.2 | 3.85 | 0.66 | 0.45 | 0.48 to 0.78 | 13.5 |
| T-motor MN4006 KV380 | 380.0 | 0.19 | 0.0571 | 1.7 / 1.2 | 10.74 / 8.77 | 22.4 | 8.60 | -0.22 | -0.69 | -0.87 to 0.22 | 3.8 |
| maxon EC 45 flat 70 W 36 V (411815) | 179.0 | 1.74 | 0.0404 | 3.6 / 2.6 | 5.06 / 4.13 | 44.6 | 4.05 | -1.43 | -1.65 | -2.71 to -0.56 | 1.9 |
| iPower GM5208-24 (out of envelope) | 20.8 | 13.70 | 0.1240 | 31.2 / 22.1 | 0.59 / 0.48 | 4.7 | 0.47 | 0.74 | 0.72 | 0.60 to 0.83 | 17.8 |
| iPower GM5208-24 (out of envelope) | 20.8 | 27.40 | 0.0877 | 31.2 / 22.1 | 0.59 / 0.48 | 9.5 | 0.47 | 0.48 | 0.46 | 0.21 to 0.67 | 8.9 |
| T-motor GB54-2 KV26 (out of envelope) | 26.0 | 15.00 | 0.0948 | 25.0 / 17.7 | 0.73 / 0.60 | 8.1 | 0.59 | 0.56 | 0.53 | 0.32 to 0.72 | 10.4 |

Read the copper-only efficiency as the generator's share of losses at a fixed rotor operating point. Negative values
mean the motor cannot absorb 0.27 Nm at 650 rpm with positive output: the design torque exceeds its short-circuit torque
or comes close to it. Pmax matched = Voc^2 / (4 R_ll) is the generator-limited DC power ceiling at 650 rpm.

### 2.4 Coupled rotor + generator at 1.5 m/s (MPPT load)

A fixed operating point penalises high-resistance machines unfairly, because a real charger (MPPT) lets the rotor speed
up until the generator's optimum is reached. Each variant is therefore coupled to the provisional BEM torque curve
(TSR 3.0-6.0 at 1.5 m/s). The scan runs over rpm, computing P_el = (Voc - I R_ll - 2 Vf) I with I = Q_rotor / Kt0.

| Candidate | KV | R_ll | rpm | TSR | tau (Nm) | P_shaft (W) | P_el DC (W) | V_dc (V) | I (A) | eta | at BEM edge |
|---|---|---|---|---|---|---|---|---|---|---|---|
| iPower GM4108H-120T (iFlight) | 27.0 | 11.10 | 822 | 5.57 | 0.193 | 16.6 | **12.9** | 23.6 | 0.55 | 0.78 | no |
| iPower GM4108H-120T (iFlight) | 27.0 | 22.20 | 886 | 6.00 | 0.168 | 15.6 | **10.2** | 21.5 | 0.47 | 0.65 | yes |
| iPower GM4108H-120T (iFlight) | 32.0 | 11.10 | 872 | 5.91 | 0.173 | 15.8 | **11.6** | 20.0 | 0.58 | 0.73 | no |
| iPower GM4108H-120T (iFlight) | 32.0 | 22.20 | 886 | 6.00 | 0.168 | 15.6 | **8.1** | 14.4 | 0.56 | 0.52 | yes |
| Rctimer GBM4114HS 100T | 32.0 | 15.20 | 886 | 6.00 | 0.168 | 15.6 | **10.3** | 18.3 | 0.56 | 0.66 | yes |
| Rctimer GBM4114HS 100T | 32.0 | 30.40 | 886 | 6.00 | 0.168 | 15.6 | **5.5** | 9.8 | 0.56 | 0.35 | yes |
| CubeMars GL40 KV70 | 63.5 | 4.50 | 886 | 6.00 | 0.168 | 15.6 | **9.1** | 8.1 | 1.12 | 0.58 | yes |
| CubeMars GL40 II KV82.5 (integrated driver) | 82.5 | 3.00 | 886 | 6.00 | 0.168 | 15.6 | **8.1** | 5.6 | 1.45 | 0.52 | yes |
| CubeMars GL40 II KV82.5 (integrated driver) | 86.8 | 3.00 | 886 | 6.00 | 0.168 | 15.6 | **7.4** | 4.8 | 1.53 | 0.47 | yes |
| T-motor GB4106 KV53 | 53.0 | 15.30 | 886 | 6.00 | 0.168 | 15.6 | **1.5** | 1.6 | 0.93 | 0.10 | yes |
| EMAX GB4008 66KV | 66.0 | 7.50 | 886 | 6.00 | 0.168 | 15.6 | **4.5** | 3.9 | 1.16 | 0.29 | yes |
| CubeMars RI50 KV100 (frameless inrunner) | 100.0 | 1.42 | 886 | 6.00 | 0.168 | 15.6 | **9.8** | 5.6 | 1.76 | 0.63 | yes |
| CubeMars RI50 KV100 (frameless inrunner) | 79.6 | 1.42 | 850 | 5.76 | 0.182 | 16.2 | **11.7** | 7.7 | 1.52 | 0.72 | no |
| T-motor MN5008 KV170 | 170.0 | 0.27 | 863 | 5.84 | 0.177 | 16.0 | **10.8** | 3.4 | 3.15 | 0.68 | no |
| T-motor MN4006 KV380 | 380.0 | 0.19 | 886 | 6.00 | 0.168 | 15.6 | **1.6** | 0.2 | 6.68 | 0.10 | yes |
| maxon EC 45 flat 70 W 36 V (411815) | 179.0 | 1.74 | no positive-voltage point | | | | | | | | |
| iPower GM5208-24 (out of envelope) | 20.8 | 13.70 | 790 | 5.35 | 0.206 | 17.0 | **13.9** | 31.0 | 0.45 | 0.82 | no |
| iPower GM5208-24 (out of envelope) | 20.8 | 27.40 | 870 | 5.89 | 0.174 | 15.9 | **11.6** | 30.6 | 0.38 | 0.73 | no |
| T-motor GB54-2 KV26 (out of envelope) | 26.0 | 15.00 | 850 | 5.76 | 0.182 | 16.2 | **12.1** | 24.4 | 0.50 | 0.75 | no |

"at BEM edge" = the optimum sits at the last BEM point (TSR 6.0), so the true optimum may lie at higher TSR, outside the
salvaged curve. Two conclusions:

1. With these generators the electrical optimum is at TSR ~5.5-6, above the rotor's Cp peak (TSR ~4). The BEM study
   should design the rotor for TSR 5-5.5 or extend the curve beyond TSR 6. Moving the design point there also cuts
   copper loss (`tsr_sensitivity`: at equal shaft power, Km 0.07 loses 81% at TSR 4 and 36% at TSR 6).
2. In the coupled view the spread between candidates narrows. What separates them is voltage and current: the
   GM4108H, GB54-2 and GM5208 deliver 14-31 V at about 0.4-0.6 A. The MN5008 and RI50 deliver 3-8 V at 1.5-3.2 A,
   which needs an in-hull active rectifier and boost to stay efficient.

### 2.5 Recommendation: iPower GM4108H-120T

Why:
- Its plausible Km band (0.063-0.106 Nm/sqrt(W)) reaches ~0.1, the value needed for ~60% copper-only efficiency at 650 rpm / 0.27 Nm (Km_required). Of the three in-envelope motors that can reach ~0.1 (GM4108H, MN5008, RI50) it is the only high-voltage, low-current, framed, bolt-in one.
- High-voltage low-current winding: 20-24 V peak line-line at 650 rpm and ~0.8-0.9 A, so diode and tether losses stay below 1 W; raw 3-phase AC can run ashore on a thin cable.
- Coupled to the provisional BEM curve with an MPPT load it delivers 8.1-12.9 W DC at 1.5 m/s across all four parameter readings; the best reading is the highest in the envelope and the worst is on par with the datasheet-grade CubeMars GL40 (9.1 W). Optima sit at TSR 5.6-6.0, above the Cp peak.
- D 47 fits the tail bay up to x = 89.5 with 5 mm blade pockets and 2 mm margin; hollow 8 mm bore matches SHAFT_D.
- USD 35.90 / EUR 43.90, mature SimpleFOC community support.

Total efficiency at the nominal point including the estimated iron + friction loss (1.69 W, from 0.07 A x 20 V no-load at
540 rpm scaled with rpm, low confidence):

| KV | R_ll | Km | P_el DC (W) | eta total (est.) |
|---|---|---|---|---|
| 27.0 | 11.1 | 0.1062 | 9.6 | 0.52 |
| 27.0 | 22.2 | 0.0751 | 3.1 | 0.17 |
| 32.0 | 11.1 | 0.0896 | 6.9 | 0.37 |
| 32.0 | 22.2 | 0.0633 | -2.2 | -0.12 |

Risks and how to retire them:
- KV and the resistance convention are not published: copper-only efficiency at the nominal point spans 1-65% across the four readings. Bench test before v2 CAD freeze: measure lead-to-lead R (= R_ll) and spin with a drill at known rpm to measure Voc line-line (peak = RMS x sqrt2) -> KV.
- Height listed as 22.7 (iFlight) and 32.3 (dealers); bay reserves 33 mm.
- Base hole pattern (4 x M3 on PCD 38.1) comes from a community bracket (two independent community sources agree); rotor-face pattern unconfirmed. Use a printed motor-specific adapter plate and measure on receipt.
- Unloaded towed freewheel could reach ~46-55 V peak (estimate, 2.5 m/s, TSR 6): electronics must tolerate >= 60 V.
- If the 8 mm shaft passes through the motor bore, the motor bearings become a third support; use a compliant coupling or let the stator float on a compliant mount.

Runner-up: T-motor MN5008 KV170 (Km 0.108 if 0.27 ohm is line-line) only with in-hull active rectification and boost: Voc 3.8 V at 650 rpm and 4.8 A make diode (3.9 W) and tether (8 W) losses prohibitive; D 55.6 fits only at x <= 67.

**Mount interface for CAD** (all mm):

| Feature | Value | Confidence | Source |
|---|---|---|---|
| Body diameter | 47.0 (+-0.05) | high | [src](https://shop.iflight.com/ipower-motor-gm4108h-120t-brushless-gimbal-motor-pro217) |
| Body length (replica) | 22.7 (+-0.2) | medium | [src](https://shop.iflight.com/ipower-motor-gm4108h-120t-brushless-gimbal-motor-pro217) |
| Axial bay reserve | 33.0 | medium | covers 32.3 dealer listing |
| Radial clearance bell to hull | 2.0 min | design margin | |
| Max aft position of a D47 body | x <= 89.5 (with 5 mm blade pockets) | medium | `hw_calc.py` |
| Hollow shaft (rotating, inferred) | OD 10 (-0.008/-0.012), ID 8 (+0.05/0) | high (dims) / low (rotates) | [src](https://shop.iflight.com/ipower-motor-gm4108h-120t-brushless-gimbal-motor-pro217) |
| Stator base holes | 4 x M3 on PCD 38.14, at 0/90/180/270 deg, adjacent spacing 26.97 | medium | [mount2.STEP](https://github.com/csunwc/simpleFOC-test) + [forum](https://community.simplefoc.com/t/ipower-gm4108h-120t-motor-mount-bracket-stl-file/2060) |
| Bracket clearance / counterbore | 3.6 / 6.6 | medium | same STEP |
| Centre clearance for base boss | 16 | medium | same STEP |
| Rotor coupling | through-shaft + printed coupler on the rotor face, or retaining compound (0.17 MPa needed) | low | design proposal |

Design the bay with a printed, motor-specific **adapter plate** between the hull bulkhead and the motor base. A wrong
hole pattern then costs one small reprint, not a shell. The display replica copies D 47 x 22.7, the base pattern and the
16 mm boss.

Bench test before freezing v2 (10 minutes, multimeter + drill + phone tachometer):
1. R between any two leads at 20 C = R_ll (expect 11.1 or 22.2 ohm).
2. Spin at a measured rpm, read AC volts between two leads: KV = rpm / (V_rms x sqrt2).
3. Recompute Km = 60 / (2 pi KV) / sqrt(R_ll); re-run `hw_calc.py` with the measured pair.

### 2.6 Freewheel and low-flow voltages (GM4108H)

Voc at 0.8 m/s (TSR 4, 315.0 rpm): 9.8-11.7 V peak line-to-line. Voc towed unloaded (ESTIMATE 2.5 m/s, TSR 6):
46-55 V. The v2 rectifier/buck must accept >= 60 V. Phase current at the nominal point is 0.62-0.74 A RMS: the stock #24 AWG
leads and a 3 x 0.5 mm2 tether are ample. Tether loss over an ASSUMED 5 m of 0.5 mm2 is 0.2-0.3 W.

## 3. Bearings, shaft and lip seal

| Part | d x D x B | Material | C / C0 (N) | Max rpm | Price | Source |
|---|---|---|---|---|---|---|
| S608-2RS | 8 x 22 x 7 | AISI 440C | 3296 / 1383 | 34,000 | USD 5.63 | [src](https://bearingsdirect.com/s608-2rs-stainless-steel-miniature-ball-bearing-8x22x7mm-sealed/) |
| 688-2RS | 8 x 16 x 5 | 52100 chrome (SS688-2RS stainless USD 5.63) | 1254 / 587 | 36,000 | USD 4.38 | [src](https://bearingsdirect.com/688-2rs-miniature-ball-bearing-8x16x5-sealed-mr688-2rs/) |
| 688-2RS (VXB) | 8 x 16 x 5 | 52100 | 1600 / 710 | 28,000 | USD 7.77 | [src](https://vxb.com/products/688-2rs-bearing-8x16x5-sealed-miniature) |

Load check (ESTIMATE, C_T 0.85): rotor thrust 28.3 N at 1.5 m/s and 78.5 N towed at 2.5 m/s.
Both sizes carry this with large margin, so the choice is set by space. Keep **608-2RS** (cheap, robust, many printed-seat
references). Use 688-2RS where the tail is tight: its 16 mm OD equals the lip seal OD, so one 16 mm bore can seat both.
SKF makes a stainless W 608-2RS1 ([src](https://www.skf.com/us/products/rolling-bearings/ball-bearings/deep-groove-ball-bearings/productid-W%20608-2RS1); page did not render, size from the listing).

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

Source table (Apple Rubber, ISO 3601 static radial gland, no back-up ring, verified):
CS 2.62 -> radial gland depth 2.03-2.18, squeeze 14-25%, diametral clearance max 0.13, groove width 3.56 +-0.13,
groove radius 0.13-0.51, eccentricity max 0.08 ([src](https://www.applerubber.com/src/pdf/iso-3601-metric-size-o-rings.pdf)).
General guidance: squeeze 15-30%, fill 60-85%, stretch <= 5% ([src](https://www.globaloring.com/o-ring-groove-design/)).

The calc assumes 2.5% nominal stretch onto the groove bottom and printed diameters +-0.10 mm (ASSUMED, after calibration).
Worst cases add full eccentricity within the radial clearance and the O-ring's own ID and CS tolerances.

| O-ring | radial clr | bore D | spigot OD | female wall | squeeze nom | squeeze worst | fill nom | fill worst |
|---|---|---|---|---|---|---|---|---|
| ISO 3601 / AS568 -037 (1.78 x 63.22) | 0.1 | 67.47 | 67.27 | 3.77 | 0.241 | 0.089 to 0.381 | 0.755 | 0.937 |
| ISO 3601 / AS568 -037 (1.78 x 63.22) | 0.15 | 67.47 | 67.17 | 3.77 | 0.241 | 0.059 to 0.408 | 0.755 | 0.937 |
| ISO 3601 / AS568 -141 (2.62 x 58.42) | 0.1 | 64.09 | 63.89 | 5.45 | 0.187 | 0.085 to 0.284 | 0.703 | 0.808 |
| ISO 3601 / AS568 -141 (2.62 x 58.42) | 0.15 | 64.09 | 63.79 | 5.45 | 0.187 | 0.065 to 0.303 | 0.703 | 0.808 |
| ISO 3601 / AS568 -142 (2.62 x 59.99) | 0.1 | 65.7 | 65.5 | 4.65 | 0.187 | 0.085 to 0.284 | 0.703 | 0.808 |
| ISO 3601 / AS568 -142 (2.62 x 59.99) | 0.15 | 65.7 | 65.4 | 4.65 | 0.187 | 0.065 to 0.303 | 0.703 | 0.808 |
| ISO 3601 / AS568 -143 (2.62 x 61.6) | 0.1 | 67.35 | 67.15 | 3.83 | 0.187 | 0.085 to 0.284 | 0.703 | 0.808 |
| ISO 3601 / AS568 -143 (2.62 x 61.6) | 0.15 | 67.35 | 67.05 | 3.83 | 0.187 | 0.065 to 0.303 | 0.703 | 0.808 |
| ISO 3601 / AS568 -144 (2.62 x 63.17) | 0.1 | 68.96 | 68.76 | 3.02 | 0.187 | 0.085 to 0.284 | 0.703 | 0.808 |
| ISO 3601 / AS568 -144 (2.62 x 63.17) | 0.15 | 68.96 | 68.66 | 3.02 | 0.187 | 0.065 to 0.303 | 0.703 | 0.808 |
| ISO 3601 / AS568 -145 (2.62 x 64.77) | 0.1 | 70.6 | 70.4 | 2.2 | 0.187 | 0.085 to 0.284 | 0.703 | 0.809 |
| ISO 3601 / AS568 -145 (2.62 x 64.77) | 0.15 | 70.6 | 70.3 | 2.2 | 0.187 | 0.065 to 0.303 | 0.703 | 0.809 |
| ISO 3601 / AS568 -230 (3.53 x 63.09) | 0.1 | 70.41 | 70.21 | 2.3 | 0.178 | 0.097 to 0.254 | 0.729 | 0.818 |
| ISO 3601 / AS568 -230 (3.53 x 63.09) | 0.15 | 70.41 | 70.11 | 2.3 | 0.178 | 0.083 to 0.268 | 0.729 | 0.818 |

Hull inner diameter at the equator: 70.2 mm.

**Pick: ISO 3601-143 (61.60 x 2.62), EPDM 70** (Apple Rubber lists EPDM for water applications; use silicone grease, not
petroleum grease). The male groove sits on the spigot:

| Gland feature | Value (mm) |
|---|---|
| Groove bottom D (on spigot) | 63.14 |
| Groove width | 3.56 (+-0.13), corner radius 0.13-0.51 |
| Radial gland depth | 2.105 (table 2.03-2.18) |
| Spigot OD | 67.15 (radial clearance 0.1) |
| Female bore D | 67.35 (female wall at equator 3.83) |
| Spigot ID max (1.2 mm wall under groove) | 60.74 |
| Squeeze nominal / worst | 0.187 / 0.085 to 0.284 |
| Fill nominal / worst | 0.703 / 0.808 |
| Stretch range | 0.0150 to 0.0352 |

The spigot OD of 67.15 lies in the 62-68 target. The 0.20 mm diametral clearance exceeds the 0.13 mm table limit,
which exists to stop extrusion under pressure. At less than 0.1 bar of submersion that does not apply, and a tighter
printed fit would bind. The 1.78 section is too sensitive to print error (worst squeeze 9% to 38%). The -145 and -230 sizes
leave less than 2.4 mm of female wall. Params currently say CS 2.0: change to 2.62.

## 5. Small parts

| Part | Spec | Key data | Source | Conf. |
|---|---|---|---|---|
| Magnet N52 4x2 | 4 +-0.05 x 2, Ni-Cu-Ni, Br 1.43-1.48 T, 80 C max | pull 4.05-4.58 N to steel; USD 0.30 each | [src](https://www.jc-magnetics.com/Magnet-N52-4mm-2mm-Discs), [src](https://themagnetbaron.com/collections/4mm-disc-magnets) | medium |
| Magnet N52 3x2 / 5x2 | alternatives | 2.54 N / 5.52 N | [src](https://themagnetbaron.com/collections/3mm-disc-magnets), [src](https://themagnetbaron.com/collections/5mm-disc-magnets) | medium |
| Dowel pin 3 x 16 | ISO 8734, type A or C1 stainless | d 3 m6, chamfer c 0.5, l 15.5-16.5; d=3 lengths 8-30 | [ISO 8734](https://cdn.standards.iteh.ai/samples/20002/ba95529de6a64437ac7454e14dbd2ef6/ISO-8734-1997.pdf) | high |
| Heat-set insert M3 | CNC Kitchen standard M3 x 5.7 | L 5.7, D1 4.6, D2 3.9, hole 4.0, min wall 1.6; EUR 12.09/100 | [src](https://www.3djake.com/cnc-kitchen/threaded-inserts-m3-standard) | high |
| Grub screw M3 | ISO 4029 cup point, A2/A4 | pitch 0.5, hex 1.5, cup d <= 1.4, socket depth >= 1.2; M3x3, M3x4 sold | [src](https://www.fasteners.eu/standards/iso/4029/) | medium |
| Cable gland | Hummel HSK-K M12x1.5 1.209.1200.50 | clamp 3-6.5, GL 8, SW 15, H 21, 1.5 Nm, IP68 10 bar/30 min with extra O-ring (IP54 without) | [src](https://www.hummel.com/en/cable-glands/products/plastic-cable-glands/hsk-k/1209120050-hsk-k-grey-m12x1-5/) | high |
| Gland alternative | LAPP SKINTOP ST-M 12x1.5 (53111000) | clamp 3.5-7, GL 8, SW 15, L 30, IP68 5 bar | [src](https://products.lappgroup.com/online-catalogue/cable-glands/skintop-cable-glands-plastic-metric/standard/skintop-st-m.html) | high |
| Gland hole | M12x1.5 clearance | 12.5 drill, or print/tap the thread | [src](https://www.metric-bolt-torque.com/metric-thread/cable-gland-sizes/) | medium |
| Eye bolt | DIN 580 | smallest in standard M8 (an A4 M6 exists with the same eye): d3 36, h 36, l 13 | [src](https://www.fasteners.eu/standards/DIN/580/), [src](https://www.jwwinco.com/en-us/products/3.4-Installing-Lifting-Damping-with-Leveling-Feet-Lifting-Gear-and-Rubber-Bumpers/Lifting-eye-bolts-Lifting-eye-nuts/DIN-580-A4-Stainless-Steel-Lifting-Eye-Bolts) | high |

Notes:
- Magnet pull is quoted magnet-to-thick-steel at contact. It drops steeply behind a printed skin, so set detent
  strength on the coupon. The 3x2 and 5x2 sizes give a 0.55x to 1.2x range around the 4x2.
- ISO 8734 specifies m6. One vendor guide calls it h6, but the difference is micrometres and irrelevant against
  printed holes. Hinge: pin pressed into the hub (PRESS_FIT), blade turning on the pin (CLEAR_ROTATE).
- A DIN 580 eye is 36 mm tall on a 75 mm nose. Use a printed eye integral to the nose and pull-test it to >= 3 x the
  towed thrust estimate (>= 236 N).

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
| Blue Freedom Portable | 120 mm | 5 W max (15 V nominal) | 'at least 1 m/s' | 690 g | body 220 x 55 mm | EUR 299 | prop stored in body, flexible shaft | [src](https://blue-freedom.net/hydro/portable/), [src](https://gigazine.net/gsc_news/en/20160418-blue-freedom/) |
| Ampair Aquair UW | 318 mm (12.5 in) | 70 W @ 6 mph, 100 W @ 9 mph | 2.7 / 4.0 m/s | 13.6 kg | 14.5 in long | not for sale | fixed prop, pole or tow | [src](https://nooutage.com/aquair-uw.htm) |
| **HydroEgg (this project)** | 194 mm deployed | 18.38 W shaft @ 1.5 m/s (BEM, provisional); 8.1-12.9 W DC with GM4108H (MPPT) | 1.5 m/s design | n/a | stowed 75 x 165 mm | display build parts EUR 14.25 (priced items) | blades fold forward into tail pockets | `out/work/hardware/hw_calc.py` |

Power sanity check, P = 0.5 rho A V^3 Cp (computed in `hw_calc.py`; Cp_el = electrical output over kinetic flux):

| Product | Rotor D (m) | P (W) | V (m/s) | P_avail (W) | Cp implied | V needed @ Cp_el 0.20 | V needed @ Cp_el 0.30 |
|---|---|---|---|---|---|---|---|
| Enomad Estream | 0.21 | 5.0 | n/a | n/a | n/a | 1.13 | 0.99 |
| Blue Freedom Portable | 0.12 | 5.0 | 1.0 | 5.7 | 0.884 | 1.64 | 1.43 |
| WaterLily (7 in housing as upper bound on rotor) | 0.1778 | 15.0 | 3.219 | 414.0 | 0.036 | 1.82 | 1.59 |
| Aquair UW @6 mph | 0.3175 | 70.0 | 2.682 | 763.9 | 0.092 | 2.07 | 1.81 |
| Aquair UW @9 mph | 0.3175 | 100.0 | 4.023 | 2578.2 | 0.039 | 2.33 | 2.03 |
| HydroEgg BEM shaft (provisional) | 0.194 | 18.38 | 1.5 | 49.9 | 0.368 | 1.84 | 1.61 |

Reading:
- Blue Freedom's "5 W, at least 1 m/s" with a 120 mm propeller implies Cp 0.88, above the Betz limit of 0.593. At a
  plausible Cp_el of 0.30 it needs ~1.43 m/s.
- The Estream (210 mm, 5 W) is plausible at ~1.0-1.1 m/s.
- WaterLily's 15 W peak sits at 3.2 m/s, i.e. towed-boat speeds.
- HydroEgg's 194 mm rotor gives a provisional 18.4 W shaft and 8-13 W DC at 1.5 m/s. Per swept area it is in the same
  class as the Estream (power scales with V^3, and the Estream's rating speed is unpublished). Its distinctive points are
  stowage and openness. Blades fold flush into a 0.401 L egg, against
  0.813 L for the Estream's 245 x 65 mm tube, and the design is fully parametric and
  documented. The electrical figure rests on the generator bench test above.
- For reference, a full-scale 1.58 m river HAHT designed for 1 kW at 1.5 m/s measured Cp 0.4087 ([Chica et al. 2025](https://doi.org/10.3390/en18071657)).
  The provisional constrained-BEM Cp of 0.367 is consistent with that.

## 8. BOM

| Item | Qty | Spec | Est. EUR | Display | v2 only |
|---|---|---|---|---|---|
| Deep groove ball bearing S608-2RS | 2 | 8x22x7, AISI 440C stainless, 2RS | 10.04 | yes |  |
| Shaft 8 mm | 1 | MISUMI PSSFJ8 SUS440C hardened 56 HRC, g6; length set by the layout study | n/a | yes |  |
| Dowel pin ISO 8734 3x16 | 3 | d 3 m6, l 16, type A or C1 stainless; one per blade hinge | n/a | yes |  |
| Magnet N52 disc 4x2 | 12 | 4 +-0.05 x 2 mm, Ni-Cu-Ni; qty ESTIMATE: 3 blades x 2 detents x 2 | 3.24 | yes |  |
| Heat-set insert M3 x 5.7 (CNC Kitchen standard) | 8 | hole 4.0, min wall 1.6; qty ESTIMATE (spinner, fins, plate) | 0.97 | yes |  |
| Socket head cap screw M3 (ISO 4762), A2 | 8 | lengths set by CAD; qty ESTIMATE | n/a | yes |  |
| Grub screw M3 x 4 cup point (ISO 4029), A2/A4 | 2 | hex 1.5; hub on shaft flat, v2 coupler | n/a | yes |  |
| Filament PETG | 1 | all wetted parts; mass set by CAD | n/a | yes |  |
| Printed tow eye (nose) | 1 | integral to the nose shell; pull test >= 236 N | 0.00 | yes |  |
| Generator iPower GM4108H-120T | 1 | D47 x 22.7, 24N22P, 11.1 ohm, hollow 8 mm bore | 43.90 |  | yes |
| Lip seal TC 8x16x5 NBR | 1 | double lip + garter spring; housing bore 16 | n/a |  | yes |
| O-ring ISO 3601-143 (AS568-143) | 2 | 61.60 x 2.62, EPDM 70 or NBR 70 (one spare) | n/a |  | yes |
| Cable gland Hummel HSK-K M12x1.5 grey + additional thread O-ring | 1 | 1.209.1200.50, clamp 3-6.5, IP68 with O-ring | 0.93 |  | yes |
| 3-core cable 3 x 0.5 mm2 | 1 | OD 3-6.5 mm; length ASSUMED 5 m | n/a |  | yes |
| 3-phase rectifier + buck to 5 V USB | 1 | input >= 60 V (towed freewheel estimate); not specified here | n/a |  | yes |

Priced items: display build EUR 14.25, v2 additions EUR 44.83.
9 items have no sourced price yet (null in `bom.json`). USD converted at the ECB reference rate of
1.1225 USD/EUR (2026-10-02).

## 9. Open issues

- GM4108H KV and resistance convention unverified. Bench test (section 2.5) before freezing v2.
- GM4108H rotor-face hole pattern unconfirmed; base pattern from community sources. Use an adapter plate.
- Whether the GM4108H hollow shaft rotates is inferred from the encoder variant. Confirm on receipt.
- The BEM curve used for the coupled study is a salvaged provisional run. Several optima sit at its TSR 6 edge, so
  extend it to TSR 7 and re-run.
- Lip-seal garter spring material for water, and prices for seals, O-rings, shaft, pins and screws, still to be sourced.
- A1/P1S build volume from a search extract (the Bambu tech-spec pages return 403 to scripts).
