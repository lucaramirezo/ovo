# Layout: rotor and internal, merged

Status: v1 layout, 2026-10-03. Sources: `hydroegg/layout.py` (rotor half) and
`hydroegg/layout_internal.py` (internal half and merge). Outputs: `data/layout.json`
(schema `hydroegg.layout/1`, the single consistent layout), plus the two half files it is built
from, `data/layout_rotor.json` and `data/layout_internal.json`. Every number below comes from
`data/layout.json`; values marked ESTIMATE are modelling assumptions, not sourced numbers.

Regenerate in this order:

```
python -m hydroegg.layout            # rotor half, about 60 s, writes data/layout_rotor.json
python -m hydroegg.layout_internal   # internal half + merge, about 17 s
```

`layout_internal` reads the rotor tables through `hydroegg.layout` (`layout.DATA_JSON`). If that
file is missing or its hinge, static tail or fold angle differ from `layout.CHOSEN`, it re-runs
`layout.run()` first, so the two halves cannot drift apart.

## Summary

| result | value |
|---|---|
| static tail end / hub front / shaft end | x 108 / 108.6 / 128 |
| hinge (tangential pin, 3 x 22) | x 116, r 14, fold 90 deg |
| blade pockets in the static hull | x 31.75 to 108; 3.0 mm deep at the equator, 6.42 max in the aero zone, 11.52 in the root notch |
| generator bay | D47 x L24 at x 22.5 to 46.5 (CubeMars GL40 KV70, GM4108H base fits) |
| bearings | 608-2RS x 63 to 70 (locating), 688-2RS x 96 to 101 (floating), lip seal x 102.5 to 107.5 |
| bearing span / overhang to the thrust centroid | 32.0 mm / 37.8 mm (ratio 0.85); tip runout 0.30 mm, shaft deflection 0.015 mm, s0 28.5 / 22.1 |
| equator joint | AS568-141 O-ring, seal bore D64.66, spigot bore D56.45; 3 lugs locked at theta 72.1 / 192.1 / 312.1, 35 deg twist |
| roll fins | Clark Y 40 x 75 at theta 84.5 / 264.5, x 15.5 to 55.5, incidence 13.75 deg, 240 mm tip to tip |
| mass, display PLA / v2 PETG | 275 g (+1.17 N net) / 402 g (+0.15 N net, no ballast) |
| tow design load (3 x working) | 234 N |
| interface checks between the halves | 21 of 21 pass |

The rotor half was taken as fixed and the internal half adapted to it. Three internal choices
changed in the merge: the bearing-span criterion, the generator bay length, and the fin
placement rule (section 5).

## 1. Frame and conventions

X is the body axis, flow travels +X, nose at x = 0, Z up. The two halves use different azimuths:

* rotor half: `psi` about +X from +Z, `e_r = (0, -sin psi, cos psi)`; blade i sits at psi = 120 i.
  Pocket azimuths (`az_lo`, `az_hi`) are relative to the blade, + toward `e_t` (direction of
  rotation).
* internal half: `theta` about +X from +Z toward +Y, so `theta = -psi`. The blades sit at theta
  0 / 240 / 120 and a pocket spans `theta_blade + [-az_hi, -az_lo]`.

Pocket depth, floor and keep-out are radial. The depth at each x is the deepest point over
azimuth. The pockets are asymmetric: the stowed blade's lead and twist push them toward the
direction of rotation, by up to 44.5 deg on one side and 10 to 20 deg on the other.

## 2. Axial stations

| part | x from (mm) | x to (mm) | r (mm) | owner |
|---|---|---|---|---|
| tow eye (ring, hole centre at x -5.0) | -14 | 4 | - | internal |
| nose plug (infill, bridged) | 0 | 11.03 | - | internal |
| M12 gland boss (theta 180) | 2 | 22 | - | internal |
| fins + rail pads (theta 84.5 / 264.5) | 15.5 | 55.5 | - | internal |
| motor carrier | 18.5 | 22.5 | 23.4 | internal |
| generator bay D47 x L24 | 22.5 | 46.5 | 23.5 | internal |
| CubeMars GL40 KV70 body | 25 | 46.5 | 23.25 | internal |
| blade pockets (static hull, 3x) | 31.75 | 108 | - | rotor |
| socket ring (front shell) | 41.59 | 57.5 | 32.33 | internal |
| bayonet spigot (aft shell) | 43.49 | 57.5 | 32.23 | internal |
| claw A (motor side) | 46.5 | 51.5 | 11 | internal |
| claw B + thrust collar | 52.5 | 63 | 11 | internal |
| shaft (8 mm) | 52.5 | 128 | 4 | internal |
| equator joint plane | 57.5 | 57.5 | - | internal |
| bearing tube (aft shell) | 62 | 108 | 13.4 | internal |
| 608-2RS (locating) | 63 | 70 | 11 | internal |
| 688-2RS (floating) | 96 | 101 | 8 | internal |
| seal shoulder | 101 | 102.5 | - | internal |
| tail end block (seal boss, notch floors) | 101 | 108 | - | internal |
| stowed-detent magnet seats (3x, notch floor) | 102.25 | 106.25 | 20 | rotor |
| lip seal TC 8x16x5 (v2) | 102.5 | 107.5 | 8 | internal |
| static tail end face | 108 | 108 | - | both |
| hub gap | 108 | 108.6 | - | both |
| stop block (hub front ring) | 108.6 | 113.5 | - | rotor |
| hub (rotating) | 108.6 | 128 | - | rotor |
| spinner shell | 108.6 | 165 | - | rotor |
| M3 grub on the shaft D-flat | 111.6 | 111.6 | - | rotor |
| hinge pins (tangential, r 14) | 112 | 120 | 14 | rotor |
| heel, deployed | 113.5 | 118.5 | - | rotor |
| rotor thrust centroid (0.7 R, ESTIMATE) | 136.32 | 136.32 | - | both |

![Merged side section](../figures/layout_section.png)

`layout_section.png`: top half in the blade 0 plane (theta 0) with the deepest pocket floor, the
stowed aero span of blade 0, the rotor keep-out line, the hub core with blade 0's clevis slot,
the pin, knuckle, heel and stop block; bottom half at theta 180, between pockets.

## 3. Rotor half (hydroegg/layout.py)

Summary of `data/layout_rotor.json`. The full tables (blade axis, pocket table, hub slot, detents)
are copied into `layout.json` under `rotor`.

* **Hinge.** A hinge sweep over x 116 to 118, r 13.5 to 14.5 and a static tail 6, 8 or 10 mm ahead
  of the hinge (11 points) picked x 116, r 14, static tail x 108: the largest smallest margin
  (pin cover 2.46 mm in the hub core, stowed heel 1.35 mm inside the spinner, 4.9 mm heel face for a
  4 mm magnet). A static tail at x 110 fails the 0.8 mm root clearance (0.77).
* **Blade root.** An 8 x 14 mm knuckle between two 3.7 mm hub ears on a 3 x 22 tangential pin. A
  5 x 8 x 10 mm heel bears on the hub front ring (stop face x 113.5) as the deployed stop. The
  shank leaves the knuckle at a 55 deg cone.
* **Pre-cone and stowage.** The aft offset is solved per station so the stowed section top is
  flush with the hull (max +0.065 mm). An in-plane lead (11.5 mm at the root to 4.8 mm at the tip)
  cuts the deepest station pocket from 12.6 mm (radial stacking) to 4.9 mm. The deployed cone is
  12.6 deg mean hinge to tip, +25 deg local at the root, -16 deg at the tip.
* **Pockets.** x 31.75 to 108, at most 6.42 mm deep in the aero zone (x <= 100) and 11.52 mm in the
  root notch, widest 29.9 mm, land between pockets at least 29.2 mm. The rotor keep-out (floor minus
  MIN_WALL along the hull normal) is at least 14.19 mm in the last 12 mm of the static tail.
* **Hub.** D-flat bore with an M3 grub at x 111.6, pins pressed into the ears, the 1.2 mm spinner
  retains them, one M3 screw into the tapped shaft end carries the thrust. Hub wall over the bore
  is 1.24 mm (MIN_WALL) in the heel slots.
* **Detents.** 12 magnets (4 x 2 N52): 2.32 N per pair at 0.4 mm (model, calibration k 0.63),
  12.9 N mm holding deployed and 25.8 N mm holding stowed, 6.3x and 12.6x the worst gravity moment
  of the 8.0 g blade (ESTIMATE, solid PLA).
* **Loads at 2.5 m/s (ESTIMATES).** Stop force 224 N, stop contact 4.6 MPa, pin shear 17.5 MPa,
  ear bearing 11.2 MPa.

![Fold, side view](../figures/layout_fold_side.png)
![Fold, end view](../figures/layout_fold_end.png)

## 4. Internal half (hydroegg/layout_internal.py)

* **Generator bay.** The motor sits in the front shell on a twist-lock carrier with its rotor face
  aft. The carrier needs r >= 23.37 mm for the PCD 38.14 adapter of the alternative motor, which
  puts the bay front at x >= 22.5. The bay is D47 x L24 at x 22.5 to 46.5. Its rotor face reaches
  3 mm into the spigot bore (4.7 mm radial clearance). GL40, GM4108H, GB4106 and MN4006 fit; the
  L 27 to 33.5 motors (GBM4114HS, GL40 II, RI50, MN5008) do not.
* **Drivetrain.** Claw coupling (claw A on the motor at x 46.5 to 51.5, 1 mm float, claw B over a
  Ruland MSC-8-SS collar at x 52.5 to 63, the thrust collar on the 608 inner race). 608 locating at
  x 63 to 70, 688 floating at x 96 to 101, seal shoulder, TC 8x16x5 lip seal at x 102.5 to 107.5,
  0.5 mm behind the static tail face. The tube is 13.4 / 10.4 mm in outer radius on 3 ribs under
  the pocket floors. The shaft is 8 x 75.5 mm (x 52.5 to 128): a flat for claw B, the hub D-flat,
  a retaining-ring groove and an M3 x 8 tapped aft end.
* **Equator joint.** The pockets are only 3.0 mm deep at x 57.5 (the internal study's fallback
  assumed 6), so the joint grew to an AS568-141 ring (stretch 3.5 %, squeeze 18.4 % nominal and
  8.1 to 28.1 % worst case, fill 62.8 %). The seal bore is D64.66, the lip D61.66 and the spigot
  bore D56.45. The socket keeps 2.15 mm of wall under the pocket floor. The three 16 deg lugs lock
  at theta 72.1 / 192.1 / 312.1 after a 35 deg twist. Their channel stays 10.8 mm clear of the
  pocket edges, and lug shear is 0.17 of the PETG interlayer strength under the 234 N tow design load.
* **Nose.** Printed tow eye (8 mm hole, 5 x 6 ring; legs at 0.22 of the interlayer strength).
  M12 gland boss at x 12, bottom, with its thread end 1.3 mm ahead of the carrier. Forward of x 11.0
  the nose is filled solid, with a 38.6 mm bridge as its back face.
* **Fins.** Two Clark Y fins, 40 x 75 mm at theta 84.5 / 264.5 and x 15.5 to 55.5, incidence
  13.75 deg. Fin torque divided by generator torque is 1.016 at 1.0 m/s and 1.062 at 1.5 m/s, with
  CL at most 0.64 of CLmax_3D. Pad clearance to the pockets is 3.0 mm, and the blades folding at the
  stowing index clear the fins by 4.1 mm (pad) and 5.2 mm (fin). Tow drag is 2.8 N, and root
  stress is 0.05 of the PETG XY strength.
* **Mass (ESTIMATES for print fill and hardware).** The display PLA model is 275 g and floats with
  +1.17 N. The v2 (PETG, real GL40) is 402 g with +0.15 N net and its CG 6.2 mm ahead of the CB. The
  target of +0.3 N needs 16 g less mass (4-perimeter shells save 22 g) or about 16 cc of float.
* **Print split.** Both shells print without supports and every part fits the 256 mm bed (front
  shell bbox 93.6 x 75 x 71.5 mm with the fin pads).

![Joint and bay](../figures/layout_section_joint_bay.png)
![Cross sections](../figures/layout_section_cross.png)
![Fins](../figures/layout_section_fins.png)

## 5. Merge decisions

**Who moved.** No rotor value changed. Every conflict was resolved on the internal side.

**Static tail 118 to 108 and the bearing span.** The bearings and seal follow the static tail end,
so they moved 10 mm forward. The internal half also now takes the rotor load plane from the solved
curved blade axis: the thrust centroid (0.7 R, ESTIMATE) sits at x 136.3, against 127.5 for a
straight 12 deg cone from the new hinge. The overhang from the 688 centre is therefore 37.8 mm. With
the GL40 alone at its forward-most position the span is 34.5 mm (ratio 0.91). No bay that holds the
recommended motor meets the 1.25 span/overhang rule of thumb the internal study used. The loads are
small (s0 above 20), so the rule was replaced by explicit checks:

| criterion | limit | value |
|---|---|---|
| span | >= 4 d = 32 mm | 32.0 mm |
| blade-tip axial runout from 0.05 mm seat slop | <= 0.5 mm | 0.30 mm |
| static safety C0 / P0 (tow case) | >= 4 | 28.5 (608), 22.1 (688) |
| shaft deflection at the load plane (radial load + asymmetry moment) | <= 0.05 mm | 0.015 mm |
| span / overhang (reported only) | 1.25 rule of thumb | 0.85 |

**Generator bay length.** The span criterion caps the rotor face at x 46.5. The bay is now D47 x
L24: the longer of the GL40 (L 21.5) and the GM4108H base motor (L 22.7), plus 1 mm. The hardware
study's L33 reserve covered the GM4108H encoder variant. At its forward-most position it would cut
the span to 23.0 mm (tip runout 0.42 mm), so it was dropped.

| bay L (mm) | rotor face x | span (mm) | span / overhang | tip runout (mm) | passes |
|---|---|---|---|---|---|
| 21.5 | 44.0 | 34.5 | 0.91 | 0.28 | yes |
| 24 (chosen) | 46.5 | 32.0 | 0.85 | 0.30 | yes |
| 28 | 50.5 | 28.0 | 0.74 | 0.35 | no (span) |
| 33 | 55.5 | 23.0 | 0.61 | 0.42 | no (span) |

**Pockets at the joint.** The real pockets are shallower at the equator than the fallback (3.0 mm
instead of 6), but they extend forward to x 31.75 and are asymmetric. The joint keeps 2.15 mm of
socket wall under the floor, so the local shell thickening the rotor study asked for is the socket
and spigot rings themselves. The lugs were moved from the symmetric 60 deg offset into the land
centres (72.1 deg), and the twist went from 24 to 35 deg.

**Hub, bearing and seal axial stack.** There is no overlap. The seal ends at x 107.5 (0.5 mm inside
the static tail face) and the hub starts at x 108.6, after the 0.6 mm hub gap. The shaft end equals
the hub aft face (x 128) for the spinner screw, and the grub's D-flat (x 111.6) lies inside the
hub. The tail end block (x 101 to 108) is solid from the seal bore out to the notch floors. It
carries the three stowed-detent magnet seats (centre r 20; 9.8 mm of material between their keep-out and the seal bore).

**Keep-out definitions.** The rotor keep-out is the floor minus MIN_WALL (1.2 mm). The internal half
keeps 2.0 mm (`POCKET_FLOOR_WALL`) under every floor, so it is 0.8 mm stricter. The bearing tube
clears the floor shell by 2.93 mm and the rotor keep-out by 3.79 mm, both at the tail notch.

**Fins.** With the hinge at (116, 14), the any-azimuth fold torus limits the fin TE to x 31. That
rule gives 25 x 110 fins at x 4.7 to 29.7, 301 mm tip to tip. The blades must sit over their
pockets before they fold, so the fins now only clear the sweep at the stowing index. This is
checked on the rotor half's swept blade surface, fold 0 to 90 deg in 1 deg steps, all three blades.
The fins can then use the whole front shell (TE at x 55.5, 2 mm ahead of the joint) and come out at
40 x 75 mm, 240 mm tip to tip. They are rotated 5.5 deg off +-Y to clear the asymmetric pockets.

## 6. Interface checks

| check | value (mm) | limit (mm) | pass |
|---|---|---|---|
| static tail end: rotor half = internal half | 108 | 108 | yes |
| hinge (x, r): rotor half = internal half | 116 / 14 | 116 / 14 | yes |
| shaft aft end = hub aft face (spinner M3 screw thread) | 128 | 128 | yes |
| hub gap = HUB_GAP | 0.6 | 0.6 | yes |
| lip seal outer face ahead of the hub front (axial gap) | 1.1 | 0.6 | yes |
| no axial overlap: 688 / shoulder / seal inside the static tail | 107.5 | 108 | yes |
| M3 grub (shaft D-flat) inside the hub on the shaft | 111.6 | 108.6 / 128 | yes |
| bearing tube to pocket floor (internal 2.0 mm floor wall) | 2.93 | 0 | yes |
| bearing tube to rotor keep-out (floor - MIN_WALL) | 3.79 | 0 | yes |
| 688 tube radius under the tail notch keep-out | 10.4 | 14.19 | yes |
| stowed-detent magnet keep-out to the seal bore (end block material) | 9.8 | 2 | yes |
| socket wall under the pocket floor at the joint | 2.15 | 2 | yes |
| bayonet twist channel to pocket edges (arc) | 10.79 | 2 | yes |
| generator bay to the shell inner wall and pocket floors (radial) | 4.07 | 0.5 | yes |
| fin rail pads to pocket edges (arc) | 3.01 | 2 | yes |
| fins to blades folding at the stowing index | 4.14 | 2 | yes |
| bearing span rule (span >= 4 d, tip runout, s0, shaft deflection) | 32 | 32 | yes |
| hub wall over the shaft bore at the heel slots (rotor) | 1.24 | 1.2 | yes |
| rotor interference (folding, neighbours, heel, deployed span) | - | - | yes |
| gland thread inner end ahead of the carrier | 1.3 | 0 | yes |
| joint lug shear utilisation (interlayer proxy, ratio) | 0.17 | 0.5 | yes |

**Cross-study check (fails, open).** The fins balance Q = 0.282 N m at 1.5 m/s (Cp 0.35 at TSR 4,
the layout rule). The BEM's best-net operating point for the GL40 is TSR 5.75, 849 rpm and
0.176 N m (`data/blade_design.json`, `generator_coupling_estimate`). There, fixed fins give 1.70x
the generator torque and roll the hull the other way.

## 7. Proposed params.py updates (consolidated)

| param | value |
|---|---|
| `STATIC_TAIL_X` | 108 |
| `HINGE_X` / `HINGE_R` | 116 / 14 |
| `ROTOR_R_TIP_TARGET` | 97 (unchanged) |
| `PIN_L` | 22 (ISO 8734 3 x 22, VERIFY stock) |
| new `FOLD_DEG`, `KNUCKLE_D`, `KNUCKLE_L`, `HUB_EAR_T`, `SHANK_CONE_DEG` | 90, 8, 14, 3.7, 55 |
| new `HEEL_REACH`, `HEEL_T`, `HEEL_W` | 8, 5, 10 |
| new `SHAFT_AFT_X`, `POCKET_CLEAR`, `ROOT_CLEAR`, `POCKET_FLOOR_WALL` | 128, 0.3, 0.8, 2.0 |
| `EQUATOR_X` | 57.5 (kept) |
| `BEARING` (608-2RS) / `BEARING_ALT` (688-2RS) | x0 63 (locating) / x0 96 (floating) |
| `LIP_SEAL` | TC 8x16x5, x0 102.5 |
| `GENERATOR` | CubeMars GL40 KV70 D46.5 x 21.5, mount PCD 24 (4 x M3); bay D47 x L24, x 22.5 to 46.5; alternative GM4108H (PCD 38.14) |
| `ORING_EQUATOR` | cs 2.62, AS568-141 (ID 58.42), groove D60.45, bore D64.66 |
| new `BAYONET` | 3 lugs 16 deg, h 1.5, t 3, twist 35, lock theta 72.1 / 192.1 / 312.1 |
| new `SHAFT` | d 8, x 52.5 to 128 (L 75.5) |
| new `COUPLER`, `CARRIER` | claw OD 22 (5 + 1 + 10.5), Ruland MSC-8-SS; carrier x0 18.5, t 4, r 23.4 |
| new `FIN` | Clark Y, chord 40, span 75, incidence 13.75, x 15.5 to 55.5, theta 84.5 / 264.5 |
| new `GLAND`, `TOW_EYE` | M12x1.5 at x 12, theta 180, boss D20; eye hole 8, ring 5 x 6 |

The reasons for each value are in `data/layout.json`, `proposed_param_updates.*_reason`.

## 8. Open issues

1. **Bearing span.** 32.0 mm against a 37.8 mm overhang (ratio 0.85). The explicit checks pass,
   but measure rotor wobble on the prototype. Only a longer static tail would restore the ratio,
   and the rotor sweep rejects x 110.
2. **Roll torque basis.** The fins are sized for 0.282 N m at 1.5 m/s, while the BEM best-net load
   line for the GL40 takes 0.176 N m (1.70x over-balance). Size the fins, or an adjustable
   incidence, against the chosen generator load line.
3. **Bay length.** The L33 reserve (GM4108H encoder variant, GL40 II) no longer fits.
4. **Roll angle.** Nothing restores roll (pendulum 0.0004 N m). A 6.2 % mismatch rolls the hull at
   3.1 rpm, so v2 needs a two-leg bridle (7.5 deg heel) or a cable swivel.
5. **Fold procedure.** Fins come off for transport. With the fins on, fold the blades only at the
   stowing index.
6. **BEM inputs.** BEM ran with hinge (124, 18.5) and a straight 12 deg axis. Re-run it with hinge
   (116, 14) and the solved curved axis. It should also check the lead offset, which yaws sections
   by up to 21 deg.
7. **Magnet detents.** Magnet forces are modelled. Check the detent feel on a test print; the
   stowed hold may keep the flow from opening the blades.
8. **O-ring squeeze.** Worst-case squeeze is 8 %, under the 14 % table minimum (acceptable below
   0.1 bar). Verify with a printed test ring.
9. **Assembly flags.** The shaft flat must pass the lip seal (needs a sleeve), and the claws engage
   blind.
10. **Pin stock and hub wall.** Pin 3 x 22: check stock. The hub wall over the bore is exactly
    MIN_WALL in the heel slots.
11. **Unsourced values.** Rotor loads, the joint shear proxy, hull C_D 0.10, the radial load
    fraction and the 0.7 R thrust centroid are ESTIMATES. Printed PETG shear allowables were not
    sourced.
