# CAD: parts, assembly and checks

Status: v1 CAD, integrated 2026-10-04; design review findings T01 to T31 applied 2026-10-05
(section 9). Sources: `hydroegg/cad/*.py` (one module per subassembly, see the table below) and
`hydroegg/cad/assembly.py` (integration). Every dimension comes from `hydroegg/params.py` and the
layout tables in `data/`, except where the review round changed the CAD: those values are listed
in section 10 as proposed params changes (params.py, profile.py, layout*.py and bem.py were not
edited). The CAD numbers below are measured on the real solids and come from
`out/cad/checks/*.json`.

Regenerate everything with one command (about 900 s measured on 2026-10-05, most of it the
modules' own checks; logs in `out/work/build/`):

```
.venv/bin/python scripts/build.py            # modules, assembly, renders, summary table
.venv/bin/python scripts/build.py --fresh    # also rebuild the cached shell and part BREPs
.venv/bin/python scripts/build.py --summary  # only print the table from the existing checks
```

`build.py` runs the six part modules through `scripts/cadrun.sh` (at most 2 OCC jobs on the
machine), then `python -m hydroegg.cad.assembly --force` (rebuilds the part cache under
`out/work/cad/integrate/cache/`), then the seven renders through the viewer, and prints a summary
table (parts, interference, key clearances, spin sweep, viewer check, print plan, outputs).

## 1. Outputs

| output | what |
|---|---|
| `out/cad/step/<part>.step` | every part in its assembled (world) pose, except `blade.step` and the coupon parts, which are in print pose |
| `out/cad/stl/<part>.stl` | every printed part in print orientation (on z = 0, mm); the three blades share `blade.stl`, the two fins and the two fin blanks are identical pairs |
| `out/cad/step/hydroegg_assembly_stowed.step` | the machine as carried: PACKED state (blades stowed, fins off, blanking strips in), 54 named parts, 28.8 MB |
| `out/cad/step/hydroegg_assembly_deployed.step` | the machine running: blades deployed, fins fitted, 54 named parts, 34.7 MB |
| `out/web/hydroegg.glb`, `out/web/manifest.json` | viewer scene: 61 meshes, 276k triangles, 9.7 MB; roles, rotor flags, `packed` flags, explode vectors, blade pivots, cut-face section meshes, material table |
| `out/web/packed_kit.glb`, `out/web/packed_kit.json` | a separate light scene: the packed egg with the two fins laid out in front of the nose (24 meshes) |
| `out/renders/*.png` | hero_deployed, stowed, packed_kit, cutaway, exploded, side, front (1600 x 1000) |
| `out/cad/checks/<module>.json` | each module's own checks |
| `out/cad/checks/assembly.json` | part audit, STL audit, interference matrix, key clearances, spin sweep, viewer fold check, print plan, BOM |

## 2. Frame, states and conventions

* X is the body axis, flow toward +X, nose at x = 0, Z up. Azimuth theta is measured about +X from
  +Z toward +Y. The rotor half uses psi = -theta; blade i sits at psi = 120 i (theta 0, 240, 120).
* **Stowing index**: rotor spin 0. Each blade then lies over its own pocket (theta = -psi), which
  is the only angle at which the blades can fold home; the pockets index the rotor in the bag.
* **Fold**: 0 = deployed, 90 = stowed, as `blade.blade_world(i, fold)`. Every part dict carries its
  world pose with the blades STOWED. A blade (and the two magnets glued into it) carries
  `pivot = {origin, axis, deployed_deg}` from `blade.pivot(i)`: rotating the stowed shape by
  `deployed_deg` about `axis` through `origin` gives the deployed blade. The integrator checks this
  on the matrices (`R(pivot) T_stowed - T_deployed` = 1.1e-14) and in the browser (section 6).
* **States** used by the checks and exports: STOWED (fold 90), FOLD_60, FOLD_30, DEPLOYED (fold 0),
  CUTAWAY (cutaway shells and internals, blades deployed) and PACKED (fold 90, fins off, the two
  blanking strips `fin_blank_0/1` in the pad grooves: the carry configuration). Fins are fitted in
  every state except PACKED. A part dict may carry `packed`: `False` = removed for packing (the
  fins), `True` = only in the packed egg (the blanks); the viewer follows the same flag (`packed=1`).

## 3. Part list

63 parts: 53 from the modules plus 9 detent magnets that no builder modelled (added by
`assembly.magnets()` from the blade's seat geometry and the shells' boss geometry) and the
sector-cut O-ring for the cutaway (`assembly.oring_cutaway()`). The coupon parts are `hidden` (not
in the assembly; the brass seal sleeve is a BOM item only); the cutaway shells and the cut O-ring
replace the full ones in CUTAWAY only. Every part is a valid solid; the compounds are intended
(several identical bought items in one part, or a bearing's two seals). Volumes are B-rep volumes.

| part | module | role | rotor | kind | cm3 | notes |
|---|---|---|---|---|---|---|
| `front_shell` | shells | shell | | printed | 47.24 | nose, teardrop tow eye, gland seat, carrier ring and twist-lock groove, bayonet socket with the 0.85 x 18 deg O-ring lead-in, 2 yawed fin pads (R3 front, R2 root blend, chamfered dovetail groove mouth) |
| `aft_shell` | shells | shell | | printed | 57.68 | spigot with O-ring groove and 3 bayonet lugs, bearing tube on 3 ribs, 3 blade pockets (smooth rim, 0.6 mm chamfer), 3 stowed-detent bosses, static tail end x 108 |
| `front_shell_cutaway`, `aft_shell_cutaway` | shells | shell_cutaway | | printed (optional) | 39.53, 43.55 | sector theta -73 to 17 removed (faces the hero / iso cameras) |
| `hub` | rotor | body | yes | printed | 14.42 | D-bore on the shaft flat, 3 slotted ear pairs, stop faces with glued block-magnet seats, flat front face |
| `spinner` | rotor | body | yes | printed | 14.47 | 1.2 mm shell over the hub, one smooth window per blade, M3 screw seat, D6 access bore with a 45 deg countersink at x 148.35 to 151.73 |
| `spinner_cap` | rotor | body | yes | printed | 1.09 | push-in tip cap: D5.7 x 12 spigot with 3 crush ribs, 45 deg shoulder, R3 Myring tip at x 160.99, fingernail groove at the joint |
| `hinge_pin_0..2` | rotor | internal | yes | bought | 0.155 each | ISO 8734 3 x 22, pressed into the ears after each blade is in its slot |
| `pin_filler_0a..2b` | rotor | internal | yes | bought | 0.027 each | ISO 8734 3 x 4, loose in the pin channels (pin retention) |
| `hub_magnet_0..2` | rotor | internal | yes | bought | 0.024 each | N52 block 6 x 2 x 2, glued in the stop-face seats (6.1 x 2.1 x 2.2) |
| `hub_grub`, `spinner_screw` | rotor | internal | yes | bought | 0.024, 0.119 | M3 x 4 cup point on the D-flat; M3 x 8 into the shaft end |
| `blade_0..2` | blade | body | yes | printed | 5.47 each | 3 identical blades, curved pitch axis, knuckle with D3.3 bore, heel stop with a block seat, shank disc seat |
| `shaft` | drivetrain | internal | yes | bought (machined) | 3.69 | 8 x 75.5, x 52.5 to 128, hub D-flat from x 108.4 out through the end, collar flat, M3 tapped end |
| `shaft_collar`, `shim_8x12x0p5` | drivetrain | internal | yes | bought | 1.18, 0.03 | thrust stop on the 608 inner race |
| `claw_a`, `claw_b` | drivetrain | internal | yes | printed | 2.68, 2.20 | dog coupling, 2.8 mm engagement, 0.3 mm flank gap, 1 mm float; claw A dogs run to the flange OD |
| `screws_claw_a_m3x4` | drivetrain | internal | yes | bought | 0.18 (2) | claw A to the generator rotor face |
| `bearing_608` (+ `_seals`) | drivetrain | internal | | bought | 2.09 | 608-2RS stainless, x 63 to 70, locating |
| `bearing_688` (+ `_seals`) | drivetrain | internal | | bought | 0.67 | 688-2RS, x 96 to 101, floating |
| `lip_seal_tc8x16x5` (+ `_spring`) | drivetrain | internal | | bought | 0.54 | v2 shaft seal, x 102.5 to 107.5 |
| `oring_141` | drivetrain | shell | | bought | 1.02 | AS568 -141, modelled installed (squeezed section); shell role so it hides with the shells |
| `oring_141_cutaway` | assembly | shell_cutaway | | bought | 0.77 | the same ring with the cutaway sector removed (CUTAWAY only) |
| `motor_carrier`, `screws_carrier_m3x8` | drivetrain | internal | | printed / bought (4) | 5.25, 0.31 | twist-lock bulkhead for the generator; the screws self-tap into the replica stator |
| `adapter_spacer` | drivetrain | internal | | printed | 3.16 | GL40-specific spacer |
| `generator_stator`, `generator_bell` | drivetrain | internal | / yes | printed | 15.16, 8.35 | display stand-in for the CubeMars GL40 (real envelope); the bell turns on a journal; the stator has D2.5 pilots, no inserts |
| `inserts_m3_bell` | drivetrain | internal | yes | bought (2) | 0.08 | heat-set inserts in the replica bell (display only) |
| `cable_gland_m12` | drivetrain | body | | bought | 2.99 | Hummel HSK-K M12x1.5 in the nose seat |
| `fin_0`, `fin_1` | fins | body | | printed | 9.62 each | Clark Y 40 x 65 at theta 93 / 273, 12 deg incidence, dovetail tongue 4.8 deep; `packed: False` |
| `fin_blank_0`, `fin_blank_1` | fins | body | | printed | 0.73 each | blanking strips for the packed egg, top flush with the pad, held by the cross pins; `packed: True` |
| `fin_cross_pins` | fins | body | | bought | 0.06 (2) | ISO 8734 2 x 10, one per fin (or blank) |
| `tolerance_coupon` (+ 2 keys) | coupon | hidden | | printed | 23.17 | print first: seats, pins, insert, M12 tap drill, bayonet and dovetail keys |
| `tolerance_coupon_b` | coupon | hidden | | printed | 19.08 | print with coupon A: hub D-bore, block and disc magnet seats (vertical and horizontal), upright ear with D2.9 / 3.0 / 3.1 holes |
| `blade_magnet_heel_0..2` | assembly | body | yes | bought | 0.024 each | N52 block 6 x 2 x 2 glued in the heel stop face, pivots with its blade |
| `blade_magnet_stow_0..2` | assembly | body | yes | bought | 0.025 each | N52 disc 4 x 2 in the shank inner face (stowed detent), 0.1 mm below the seat mouth, pivots with its blade |
| `tail_magnet_0..2` | assembly | body (`_0`: shell) | | bought | 0.025 each | N52 disc 4 x 2 in the aft shell's stowed-detent bosses, 0.05 mm below the boss face; `tail_magnet_0` takes the shell role because its boss is in the cutaway sector |

## 4. How the parts fit

**Shells and joint.** The egg splits at the equator x 57.5. The aft shell's spigot (seal zone OD
63.89) enters the front shell's socket (bore D64.09, 0.1 mm radial) and carries the -141 O-ring in
its groove. The socket mouth is a 0.85 mm x 18 deg lead-in (mouth r 32.879 at every azimuth, the
stretched ring stands to r 32.528). Three lugs pass the socket lip through entry gaps and lock
after a 35 deg twist; the channels end exactly at the locked lugs, so the hard stop aligns the
pockets across the joint (measured: 0 mm3 overlap, lugs on the stop at 0.000 mm).

**Generator bay (front shell).** The motor carrier enters through the socket with its three tabs
at theta 120 / 240 / 0, passes the ring notches and twists 30 deg toward -theta to lock at 90 / 210 /
330; the generator reaction torque drives the tabs into the stop. Carrier, adapter spacer and the
replica stator stack on 4 M3 x 8 that self-tap into D2.5 pilots (1.65 mm from the bell journal
groove; the real GL40 has M3 threads); the bell turns on the stator journal (0.298 mm min gap over
a full turn) and carries claw A on 2 M3 x 4. The cable leaves through the carrier's theta 180 notch
to the M12 gland in the nose seat (gland to carrier 3.61 mm).

**Coupling and shaft (aft shell).** The 608 (locating) and the 688 (floating) sit on three crush
ribs each in the bearing tube. The shaft runs from x 52.5 to the hub at x 128. Claw B is a sleeve
pressed over the Ruland collar; collar, shim and the 608 inner race carry the rotor thrust. Claw A
and claw B interleave with 2.8 mm engagement, 0.300 mm flank gap and 1 mm axial float, so the
shells can be closed without aligning the dogs exactly. Claw B runs 1.262 mm clear of the tube
mouth. The lip seal (v2) goes on last from the tail end over a brass sleeve.

**Rotor.** The hub slides onto the shaft's D-flat from the aft end (0.074 mm radial clearance)
until its aft face is flush with the shaft end: the shaft end is the one axial datum, and the
0.6 mm hub gap to the static tail follows from the shaft setting (checked with a feeler, held by
the M3 grub). Each blade is dropped radially into its slot and its 3 x 22 pin is pressed through
the D3.5 channel into both ears; the knuckle runs free on it (0.149 mm bore clearance). A loose
3 x 4 dowel in each pin channel and the spinner keep a loosened pin in both ears (0.42 mm free
travel, 3.28 mm far-ear engagement left). The spinner slides on over the hub at any fold
(0.675 mm minimum, deployed; folding the blades 15 deg or more gives 1.7 mm and is easier), an M3 x 8
clamps spinner and hub to the shaft end, and the tip cap pushes into the access bore on its crush
ribs. Deployed, each blade's heel lands flat on its hub stop face (0.000 mm, 0 mm3), where the heel
and hub block magnets face each other across 0.400 mm (about 2.32 N).

**Stowed blades.** Folded 90 deg at the stowing index, each blade lies in its pocket: the aero span
0.331 mm above the pocket floors, the root 0.901 mm clear of the static tail, the shank's inner
face 0.129 mm above its detent boss, and 1.10 mm from the pocket rim at the skin (the rim is the
blade's swept footprint over the whole fold plus 1.0 mm). The fins (theta 93 / 273) stay
6.521 mm and their pads 3.066 mm from the blades over the whole fold.

**Fins.** Each fin's tongue slides aft in the pad's dovetail groove to the end stop and is held by
a 2 mm cross pin (0.155 mm clearance in the tongue's D2.3 hole; 1.2 mm ligaments over the pad hole
and under the tongue hole). For packing the fins come off and the blanking strips slide into the
grooves on the same pins, flush with the pad tops.

### Assembly order

From `drivetrain.json` (assembly, 12 steps) and the rotor module:

1. Press the 2 M3 inserts into the replica bell (the stator has D2.5 pilots, no inserts).
2. Stack carrier, adapter spacer and stator with 4 ISO 7380 M3 x 8 from the carrier counterbores,
   self-tapping into the pilots.
3. Slide the bell onto the stator journal; bolt claw A to the bell with 2 ISO 4762 M3 x 4.
4. Screw the M12 gland into the nose seat, pass the cable, lay it in the carrier notch.
5. Push carrier + generator + claw A into the front shell with the tabs at theta 120 / 240 / 0,
   then twist 30 deg toward -theta to lock.
6. Press the 688, then the 608, into the aft shell's bearing tube (draw bolt).
7. Slide the shaft in from the front until it stands 20.0 mm proud of the static tail face.
8. Shim against the 608 inner race, claw B (collar pressed in) on the shaft front; tighten the
   collar's set screw through claw B's access hole.
9. (v2) Lip seal from the tail end over the brass sleeve.
10. Rotor: glue the magnets (3 hub blocks into the stop-face seats; in each blade one block in the
    heel and one disc in the shank, polarity attracting the hub and tail magnets; 3 tail discs into
    the aft shell bosses). Slide the hub onto the D-flat until its aft face is flush with the shaft
    end (straightedge), check the 0.6 mm gap with the feeler, tighten the grub. Drop each blade
    radially into its slot and press its pin through the channel into both ears; the 3 x 4 filler
    dowels; the spinner (blades folded 15 deg or more); M3 x 8 into the shaft end; push the tip cap
    in (the fingernail groove pulls it out again for the screw).
11. Grease the O-ring and roll it into the spigot groove.
12. Unfold the blades (fold 80 deg or less: the stowed blades lie in pockets of both shells, 35 deg
    apart at the bayonet entry). Push the shells together at the bayonet entry marks (the claws find
    each other on their lead-ins) and twist 35 deg to the hard stop; only then fold the blades home.
    Nothing but the hub on the static tail stops the shaft toward -X: after stowing, pull the
    spinner aft so the hub runs free.

Packing: pull the cross pins, slide the fins off, slide the blanking strips in and pin them.
Before stowing, always twist the joint to its hard stop (there is no joint detent, section 8).

## 5. Interference and clearance results

Method: manifold3d booleans and `min_gap` on meshes of the real solids (0.02 mm chord, shells
0.03 mm; every shape's stored triangulation is dropped first so the tolerance is honoured). Exact
OCC booleans against the full shells take minutes each. Overlaps up to 0.5 mm3 on curved contacts
are tessellation facets; flat contacts tessellate exactly. Every pair of parts present in a state
is tested whose boxes come within 3 mm.

| state | fold deg | parts | pairs within 3 mm | intended overlaps | zero-gap contacts | defects (> 0.5 mm3, unexpected) |
|---|---|---|---|---|---|---|
| stowed | 90.0 | 54 | 264 | 12 | 32 | 0 |
| deployed | 0.0 | 54 | 57 | 0 | 9 | 0 |
| fold_30 | 30.0 | 54 | 69 | 0 | 6 | 0 |
| fold_60 | 60.0 | 54 | 98 | 0 | 6 | 0 |
| cutaway | 0.0 | 53 | 105 | 2 | 15 | 0 |
| packed | 90.0 | 54 | 11 | 0 | 2 | 0 |

(Static pairs are computed once in STOWED; the other states add the pairs that involve a folding
part, CUTAWAY the pairs that involve a cutaway part, PACKED the pairs that involve a blank. CUTAWAY
has 53 parts: `tail_magnet_0`'s boss lies in the removed sector, so the magnet goes with the full
shell.)

**Intended overlaps** (listed in `assembly.EXPECTED`, reported separately):

| pair | kind | overlap mm3 | state | why |
|---|---|---|---|---|
| aft_shell / bearing_608 | press | 0.9424 | stowed | 608 outer race on 3 crush ribs (tips D21.9) |
| aft_shell / bearing_688 | press | 0.5709 | stowed | 688 outer race on 3 crush ribs (tips D15.9) |
| hub / hinge_pin_0 | press | 3.1727 | stowed | 3 mm pin pressed into the D2.90 printed ears |
| hub / hinge_pin_1 | press | 3.1727 | stowed | 3 mm pin pressed into the D2.90 printed ears |
| hub / hinge_pin_2 | press | 3.1727 | stowed | 3 mm pin pressed into the D2.90 printed ears |
| hub / hub_grub | thread | 7.9601 | stowed | M3 grub tapped into the PLA hub |
| spinner / spinner_cap | press | 1.0093 | stowed | tip cap: 3 crush ribs (D6.1 over the tips) in the D6 access bore |
| shaft / bearing_608 | contact | 0.0196 | stowed | inner race on the shaft (zero clearance) |
| shaft / bearing_688 | contact | 0.02 | stowed | inner race on the shaft (slip fit) |
| shaft_collar / claw_b | press | 0.2651 | stowed | Ruland collar pressed into claw B's 3 crush ribs (D15.9) |
| screws_claw_a_m3x4 / inserts_m3_bell | thread | 6.783 | stowed | M3 screws in the inserts |
| screws_carrier_m3x8 / generator_stator | thread | 19.6872 | stowed | M3 x 8 self-tapped into the D2.5 stator pilots (T21) |
| aft_shell_cutaway / bearing_608 | press | 0.6283 | cutaway | 608 outer race on 3 crush ribs (tips D21.9) |
| aft_shell_cutaway / bearing_688 | press | 0.3806 | cutaway | 688 outer race on 3 crush ribs (tips D15.9) |

The O-ring, the lip seal (0.004 mm to the shaft, 0 to its bore), the glued magnets and the fin
cross pins touch or sit in their seats without overlap.

**Zero-gap contacts** (designed): shells at the lug stop; carrier tabs on their stop; fin tongues
and blanks on the groove end stops; hub, spinner, shaft end and spinner screw clamped together;
grub on the D-flat; the stacked carrier / spacer / stator and claw A / bell / inserts faces; collar
/ shim / 608 inner race; magnets on their seat floors; pin fillers on the pin ends; deployed heels
on the stop faces.

**Key clearances on the real solids** (compared with `data/layout.json` and the design rules):

|  | check (real solids) | measured | layout / design | limit |
|---|---|---|---|---|
| PASS | hub gap: static tail end (aft shell) to hub front | 0.600 mm | 0.6 | >= 0.6 |
| PASS | spinner to aft shell | 0.600 mm |  | >= 0.6 |
| PASS | lip seal outer face to hub front (axial) | 1.100 mm | 1.1 | >= 0.6 |
| PASS | shaft to hub D-bore (radial) | 0.074 mm, 0 mm3 | 0.075 | > 0 |
| PASS | stowed aero span (x < 100) to the pocket floors and walls, 3 blades | 0.331 mm (blades 0.331, 0.331, 0.331) | 0.3 | >= 0.3 (-0.05 tessellation) |
| PASS | stowed root (x > 100, outside the detent boss) to the static tail | 0.901 mm (blades 0.901, 0.901, 0.902) | 0.8 | >= 0.8 |
| PASS | stowed blade to its detent boss face | 0.129 mm (blades 0.129, 0.130, 0.129) |  | > 0 (shells design 0.1 to the layout surface) |
| FAIL | stowed detent magnet pair (blade shank / tail boss): face-centre gap, min gap, axis angle, lateral | 0.701 / 0.513 mm, 5.7 deg, 0.00 mm; about 1.73 N (layout 2.32 N at 0.4 mm) | 0.4 | <= 0.5 face gap (design 0.4) |
| PASS | deployed detent magnet pair (blade heel / hub stop face): min gap | 0.400 mm, about 2.32 N | 0.4 | 0.4 +- 0.15 (seat depth H + 0.1: 0.3 is slightly stronger) |
| PASS | deployed heel on the hub stop face: gap / overlap (3 blades) | 0.000 mm / 0.000 mm3 | 0.0 | contact, no overlap |
| PASS | blades folding 0..90 deg at the stowing index to the fins | 6.521 mm at fold 90 | 8.44 | >= 2.0 |
| PASS | blades folding 0..90 deg at the stowing index to the fin pads | 3.066 mm at fold 90 | 3.0948 | >= 2.0 |
| PASS | blades folding 6..90 deg to the hub (run clearance; below 6 deg the heel lifts off its stop) | 0.296 mm at fold 6 |  | >= 0.3 (-0.01 tessellation) |
| PASS | blades folding 0..90 deg to the spinner | 0.664 mm at fold 0 |  | >= 0.3 |
| PASS | blades folding 0..90 deg to the hinge pins (knuckle bore) | 0.149 mm at fold 60 | 0.15 | > 0.1 (bore 3.3 on a 3.0 pin) |
| PASS | blades folding 0..90 deg to the shells (minimum is the stowed boss) | front 0.331 mm at fold 90, aft 0.129 mm at fold 90 |  | > 0 |
| PASS | neighbour blades folding together | >= 3 mm at every fold |  | >= 1.0 (layout rule) |
| PASS | claw B (rotating) in the static tube mouth | 1.262 mm | 0.95 | >= 0.3 |
| PASS | claw A to claw B dog flanks | 0.300 mm | 0.3 | >= 0.25 |
| PASS | cable gland to motor carrier | 3.609 mm | 2.275 | >= 1.0 |
| PASS | generator replica to the shells | 4.681 mm | 0.637 | >= 0.5 (bay to shell) |
| PASS | front / aft shell locked: overlap, gap | 0.0000 mm3, 0.000 mm (lugs on the hard stop) | 0.0 | no overlap |

The fold sweep runs in 3 deg steps; hub, spinner, pins and shells are 3-fold symmetric, so blade 0
stands for all three there, while all three blades are checked against the fins and pads. Below
6 deg the heel lifts off its stop face and its gap grows from zero, so the hub run clearance is
read from 6 deg (as in `rotor.py`). The one FAIL is the stowed detent (section 8, item 1).

**Rotor spin sweep** (deployed, rotating parts turned 0..350 deg in 10 deg steps against every
static part within 3 mm): 23 pairs, none rubbing. Closest running gaps: generator bell to stator
0.298 mm, shim to the 608 seals 0.389 mm, collar and claw B to the 608 0.5 mm, hub and spinner to
the static tail 0.6 mm. The only zero gaps are the bearing inner races and the seal lip on the
shaft, which turn with it.

**Module checks** (each module's JSON; all pass unless listed in section 8): rotor (hub seat walls
1.27 / 1.26 / 3.08 mm, hub min wall 1.23 mm, spinner window bridge 15 mm with a 2.14 mm ring of
cone under the seat ceiling, cap spigot 2.2 mm clear of the M3 head, slide-on at every fold);
blade (block seat walls 1.15 heel tip / 1.95 flanks, far points classified OUT in every pose,
TE 0.600 mm, STL watertight); shells (rim 1.10 mm at the skin, 0 fold-sweep points outside the
pocket, flush 0.065, walls under the floors 2.135 / 2.069, seal lead-in wall 1.595, eye crown
5.0 / 7.45, pad root blend R2, pads 2.28 mm from the pocket edges); fins (ligaments 1.2 / 1.2,
B-rep volume default / adaptive / mesh 9616.7 / 9621.9 / 9607.8 mm3); drivetrain (claw skin 1.2 over the counterbores, stator pilots
1.65 mm from the journal groove, claw B bridge on a 0.2 mm membrane); coupon (all key and probe
checks).

## 6. Web scene and viewer check

`assembly.write_web()` tessellates each part for the web (0.05 mm chord, 0.25 rad; blades 0.08 mm
/ 0.3 rad, sub-pixel at 1600 px) and writes the GLB and manifest through
`export_web.write_scene`. Roles: shells and the O-ring `shell`, cutaway shells and the cut O-ring
`shell_cutaway`, everything inside `internal` (shown in cutaway and exploded views), visible
bought parts and the detent magnets `body`, coupon parts `hidden` (left out of the GLB). The cut
faces of the cutaway shells are separate meshes (`front_shell_cutaway_section`,
`aft_shell_cutaway_section`, colour #C98B5E, lifted 0.03 mm) so the section reads against the skin.
The manifest carries the `packed` flag (fins hidden and blanks shown with `packed=1`).

Explode vectors (`assembly.EXPLODE_GROUPS`) open the egg in two levels: the shells move apart
(front shell, fins and gland -52 mm, aft shell +62 mm), the drivetrain opens in assembly order
between them (carrier screws -22, carrier -12, spacer -7, stator -3, bell +4, claw A +9, its screws
+15, claw B +17, collar +22, shaft with shim and bearings +24, O-ring +47), and the rotor goes aft
(lip seal +70, hub group +78, spinner +100, cap +112, spinner screw +122); blades move 40 mm
radially, pins and fillers 15 mm along the pin axis, fins 26 mm outward, the tail magnets out of
their bosses.

Browser check (`assembly.verify_viewer_fold`, headless Chrome through `scripts/render.py`'s
server): the viewer's read-only `window.partBoxes()` returns each part's posed vertex box in the
CAD frame. At fold = 1 every blade box matches `blade.blade_world(i, 0)` within 0.063 mm (web
tessellation), and the blade magnets match their CAD pose within 0.002 mm; at fold = 0 the boxes
match the stowed solids.

## 7. Renders

`out/renders/` (1600 x 1000, light theme):

| render | viewer URL parameters | shows |
|---|---|---|
| `hero_deployed.png` | view hero, fold 1 | running: blades deployed on their stops, fins fitted |
| `stowed.png` | view hero, fold 0, packed 1 | carried: blades in their pockets, fins off, blanking strips flush in the pads |
| `packed_kit.png` | packed_kit scene | the packed egg with the two fins beside it |
| `cutaway.png` | cutaway, fold 1 | generator, carrier, coupling, shaft, bearings, seal; copper section faces |
| `exploded.png` | view 55 / 24, fold 1, explode 1 | every part in assembly order, all three blades visible |
| `side.png`, `front.png` | side, front | profile and frontal area |

All seven were inspected after the final build (2026-10-05): the deployed blades land on their hub
stops and the stowed blades lie flush in their pockets; the spinner shows one smooth window per
blade and ends in the tip cap; the tow eye reads as a teardrop lug; the fin pads carry the R2 root
blend and the chamfered groove mouth; in `stowed.png` the fins are off and the blanking strips sit
flush in the pads; the cutaway shows the generator, carrier, coupling, shaft, bearings and seal with
copper section faces and no O-ring arc across the opening; the exploded view keeps the parts
in assembly order with all three blades in sight.

## 8. Open issues

1. **Stowed detent weaker than designed (blade + shells, T15).** The stowed blade's inner face
   leans about 49 deg from radial, the boss axis is held to 45 deg by the print rule. With the
   blade magnet 0.1 mm below its mouth and the tail magnet 0.05 mm below the boss face, the pair
   sits 0.701 mm apart face to face (0.513 mm at the closest edge), 5.7 deg apart, about 1.73 N
   instead of 2.32 N (it was 1.045 mm and 1.27 N). Scaled from the earlier 14.1 N mm at 1.27 N,
   that holds about 19 N mm against the 2.05 N mm gravity moment of a stowed blade (about 9 x;
   the layout had 12.6 x): enough for the display model. Also the boss's lowest tip starts about 1.1 mm
   (print height) above the root-notch floor ramp and cannot be supported without hitting the
   stowed root: a small first-layer island that prints in mid-air. Fixes for Luca: lean the boss
   to the blade normal (49 deg, a 4 deg overhang on a 4.1 mm hole), a 5 x 2 magnet, or move the
   detent forward (DETENTS.stowed_x about 102.5, over a higher floor). The blade module's
   `stowed_detent_interface` still compares with the layout's radial DETENTS pair (superseded).
2. **Spinner windows (rotor, T01).** The 0.4 mm outer edge chamfer could not be built (OCC chamfer
   and fillet fail on every B-spline window edge); break the edges by hand. A 2 mm trim of the hub
   band under the windows was built and dropped (it left a 1.08 mm wedge and a 12 mm2 ledge at the
   pin channel), so the hub's band steps are visible through the window. Spinner wall: 1.167 mm
   minimum on the 1.2 mm shell mesh (p01 1.185); proposal spinner_wall 1.25 (section 10).
3. **Fin pads stand 9.2 mm proud (shells, T07).** The pads are rounded (R3 front, R2 root blend,
   chamfered groove mouth) but their height comes from the fin root standoff in the layout
   (root_r_le / root_r_te): a layout decision.
4. **Fins in the packed egg (T08).** Fins come off for packing (blanking strips in). A fin
   redesign that stays on inside the rotor disc changes the couple ratio band: left for Luca.
5. **Blade trailing edge on supports (blade, T13).** The TE floats 1.7 to 3.7 mm above the bed on
   tree supports (top Z distance 0.2, dense interface, 5 mm brim on the knuckle); the sacrificial
   fence suggested by the review was not modelled.
6. **Nose-plug bridge (shells, T16).** The 38.6 mm flat bridge inside the nose is the only bridge;
   a 45 deg cone from r 19.2 would reach the carrier (x 22.7). Good bridge settings, expect a
   little sag on that internal face.
7. **No joint detent (T17).** A 0.3 mm ramp bump before the lock (mirrored in the coupon) was not
   modelled; twist to the hard stop before stowing (it aligns the pockets across the joint).
8. **Colours (T27).** Spinner, hub and fins stay charcoal: a colour decision for Luca.
9. **CLEAR_ROTATE (T31).** params.py calls it radial; layout.py uses it as diametral for the
   knuckle bore (0.15 radial), rotor.py as a 0.3 radial running clearance. Proposal in section 10.
10. **Small residual faces beyond 45 deg** (all intended or below a layer): aft shell 3 x 0.43 mm2
    at 59 deg on pocket loft facets and the 29 mm2 seal stop ring; front shell the nose bridge and
    the 12 mm2 eye-hole flat; spinner the 3 x 63.9 mm2 window-top bridges and the D22.6 seat
    ceiling (354.7 mm2, on a complete cone ring); hub small arched-hole islands (largest 8 mm2,
    anchored); `tolerance_coupon_b` its horizontal test holes (63.5 mm2, the point of the test).
11. **Front shell STL.** OCC's tessellation of the front shell leaves 2 small open loops
    (perimeter up to 4.47 mm) at sliver faces of the B-rep; `common.repair_stl` closes them and the
    written STL is watertight with consistent winding (the same for the cutaway shell). The B-rep
    itself is valid.
12. **Cutaway sector.** The carrier's theta 330 lock tab (r 27.2) lies in the removed sector beyond
    r 25: kept, it shows the twist lock.
13. **File sizes.** The blade heal (faces above 200 mm2 split, T14) grows `blade.step` to 2.0 MB and
    the assembly STEPs to 28.8 / 34.7 MB (stowed / deployed; 23.7 MB before).
14. Bought items to verify before ordering: 3 x 22 and 3 x 4 pins, the Ruland collar's set screw,
    the 8 x 12 x 0.5 shim (DIN 988 lists 8 x 14), N52 6 x 2 x 2 blocks magnetised through a 2 mm
    face. The layout's two equator magnet pairs at the joint are not modelled (no room in the
    5.05 mm equator face); the joint relies on O-ring friction and the lug stop.

## 9. Review fix round (2026-10-05)

The design review (visual, mechanical, print and conformance reviewers, triaged into T01 to T31)
was applied in one pass across the CAD modules. params.py, profile.py, layout*.py and bem.py were
not edited: every value that now differs from params is listed in section 10 as a proposal.

| id | finding | what was done | measured after |
|---|---|---|---|
| T01 | spinner windows: stair-stepped outline, cantilevered treads | one outline per blade (convex hull of the slot bands + 1.5 mm, arc joins) with a 45 deg gable top; full depth only where the heel and knuckle sweep (it closes with 45 deg flanks and a 15 mm bridge at x 125.9, leaving a complete ring of the cone under the D22.6 seat ceiling); elsewhere a recess over a 1.2 mm skin whose floor rises at 45 deg out through the surface. The 0.4 mm outer chamfer could not be built (OCC chamfer fails on every B-spline window edge). A 2 mm hub trim was built and dropped (section 8) | no treads; window tops: 3 x 63.9 mm2 bridges at print z 17.3 anchored on both flanks; seat ceiling on a complete annulus; blade to spinner 0.664 mm (fold 0, was 0.635); spinner slides on at every fold, even deployed (min 0.675 mm) |
| T02 | blunt open tip, egg 151.7 mm long | push-in tip cap (`spinner_cap`): D5.7 x 12 spigot with 3 crush ribs to D6.1, 45 deg shoulder in a countersink, profile.radius continued to an R3 round, fingernail V-groove at the joint | egg ends at x 160.99 (L/D 2.15); groove 0.30 deep; spigot 2.2 mm clear of the M3 head; cap prints support-free |
| T03 | D4 deployed magnet does not fit the stop zone | N52 blocks 6 x 2 x 2 (6 along the pin, magnetised along x) in glued 6.1 x 2.1 seats at r_e 8.2, both faces 0.2 below their stop faces | hub seat walls 1.27 (knuckle arc) / 1.26 (mouth rim) / 3.08; blade heel tip 1.15, flanks 1.95; pair gap 0.400 mm, about 2.32 N |
| T04 | pocket rim serrated and oversized | rim at the skin = blade swept footprint (fold 0..90) in the lip band + 1.0 mm; 45 deg in-slice flare to the kept depth-wise bounds; 0.6 x 45 deg rim chamfer; smooth sewn lofts | rim to blade 1.10 mm at the skin (was 1.4 to 6.6); fold sweep 0 points outside; flush 0.065; floor 0.331; walls under floors 2.135 / 2.069 |
| T05 | O-ring lead-in smaller than the ring | 0.85 mm x 18 deg lead-in (mouth D65.79), cut after the pocket liners | mouth r 32.879 at all 72 azimuths vs ring r 32.528; wall to the pocket floors over the lead-in 1.595 (>= MIN_WALL) |
| T06 | eye crown 3.3 mm, line bears on a point | parallel teardrop outer profile, hole apex cut to a 2 mm flat, R2 root fillet | crown 5.0 normal / 7.45 along the pull; 19.1 MPa, utilisation 0.41 PETG / 0.37 PLA; eye to gland body 7.08 |
| T07 | fin pads: square rails, boxy | rounded front (R3 in plan), R2 root blend all round (a mixed R4 / R1 blend failed on the tangent edge chain), 0.5 x 45 deg groove-mouth chamfer, blanking strips (`fin_blank_0/1`) for the packed egg | R2 blend built; blades to pads 3.07 mm at fold 90 (was 3.89, limit 2.0); pads still 9.2 mm proud: lowering them needs a lower fin root standoff (layout, Luca) |
| T08 | fins fitted in every state | PACKED state (fins off, blanks in): stowed.png, the stowed STEP and the viewer's `packed=1`; packed_kit.png shows the fins beside the egg | stowed width = hull 75 + pads; fin redesign left for Luca |
| T09 | O-ring arc across the cutaway | `oring_141` takes the shell role, `oring_141_cutaway` (sector removed) the shell_cutaway role | only the carrier's theta 330 lock tab still lies in the sector beyond r 25 (r 27.2, informative, kept) |
| T10 | guide pressed the pins before the blades | assembly steps rewritten: pins pressed through the channel after each blade is dropped in | rotor.py notes, drivetrain assembly step 10, sections 4 and print plan |
| T11 | hub front-face V-grooves | removed | grub thread to the front face 1.5 mm |
| T12 | two axial datums on the hub | one datum, the shaft end: hub slid flush with it (straightedge), gap checked with the feeler | docs |
| T13 | blade TE on tree-support interface | documented the support settings (top Z 0.2, dense interface, 5 mm knuckle brim); TE height corrected to 1.7 to 3.7 mm (1.7 to 2.2 along the aero span) | the sacrificial fence was not modelled (section 8) |
| T14 | blade B-rep classified far points IN | faces above 200 mm2 split + ShapeFix (orientation fixes and sewing alone did not help) | far points OUT in the print, stowed and deployed poses; volume unchanged |
| T15 | stowed pair 1.045 mm apart; boss tip prints in mid-air | blade magnet 0.1 below its mouth, tail magnet 0.05 below the boss face | 0.701 mm face gap, about 1.73 N (was 1.27 N); the boss tip could not be supported (section 8) |
| T16 | 38.6 mm nose-plug bridge | documented: good bridge settings, expect some sag on that internal face | the 45 deg cone from r 19.2 would reach the carrier (x 22.7) |
| T17 | no joint detent | documented "twist to the hard stop before stowing"; MAGNET count proposed (section 10) | bump not modelled (section 8) |
| T18 | cross-pin ligaments 1.0 / 0.55 mm | groove and tongue 0.8 deeper (4.8, 45 deg undercut kept), pin depth 2.2; 0.3 x 45 deg tongue bed chamfer; brim ears at the root ends | ligaments 1.2 / 1.2 |
| T19 | fin GProp unstable | faces above 100 mm2 split + ShapeFix | default 9616.7 / adaptive 9621.9 / mesh 9607.8 mm3 (0.05 %) |
| T20 | claw A dog skin 0.69 mm | dogs to the flange OD (16.2), flange edge chamfer removed under them | skin 1.2; no overhang |
| T21 | stator insert holes 0.9 mm from the journal groove | no stator inserts: carrier screws self-tap into D2.5 pilots (45 deg lead-in) | 1.65 mm to the groove |
| T22 | claw B collar-pocket bridge ends on the bore | 0.2 mm sacrificial membrane over the bore in the print solid (drill it out) | full D16.2 bridge anchored all round |
| T23 | "LOCK" engraving too fine | text removed, the 0.9 mm chevron stays; the 0.9 mm counterbore web documented as non-structural | |
| T24 | nothing stops the shaft toward -X | documented: after stowing, pull the spinner aft | |
| T25 | collar set screw stub short of the flat | modelled down to the flat (r 3.5); v2 keying of claw B documented | contact, no overlap |
| T26 | coupon misses several fits | `tolerance_coupon_b`: D-bore 8.10 / 8.15 / 8.25 with the 0.5 flat, block seats 6.0 x 2.0 / 6.1 x 2.1 / 6.2 x 2.2 and disc seats D4.0 / 4.1 / 4.2 vertical and horizontal, an upright 3.7 mm ear with horizontal D2.9 / 3.0 / 3.1 holes | all probes pass |
| T27 | charcoal spinner, hub and fins | not changed: a colour decision for Luca (section 8) | |
| T28 | cut faces share the skin colour | `front_shell_cutaway_section` / `aft_shell_cutaway_section` meshes (cap colour #C98B5E) lifted 0.03 mm | cutaway.png |
| T29 | exploded view 500 mm long, one blade hidden | two-level explode (shells apart, drivetrain opened between them, rotor aft), view 55 / 24 | exploded.png shows all three blades |
| T30 | closing the shells with stowed blades | step 12: unfold the blades (fold <= 80 deg) before joining, stow after the twist | docs |
| T31 | deviations not written back | section 10 and the structured proposal list | |

## 10. Proposed params changes

The CAD departs from params.py / the layout in these places (`proposed_param_updates` of the
review round, also in `out/work/cad/fix/proposed_params.json`). params.py is read-only for the CAD
modules; Luca applies or rejects each.

| params key | proposed | was | why |
|---|---|---|---|
| `ROTOR_HUB.spinner_wall` | 1.25 | 1.2 | the 1.2 shell measures 1.167 min on the mesh |
| `ROTOR_HUB.ear_hole_d` | 2.90 (printed press bore) | 3.0 | 3.0 is the nominal pin |
| `DETENTS` deployed pair | N52 block 6 x 2 x 2, hub centre (x, r_e, tau) [112.3, 8.2, 0], blade [-5.8, 0, -1.3] (blade frame), gap 0.4 | D4 x 2 disc at r_e 8.44 / [-5.56, 0, -1.3], gap 0.3 | a D4 disc does not fit the 4.5 mm stop zone (T03) |
| `DETENTS` stowed pair | static centre [103.1233, -6.5048, 19.2891], axis [0.7069, -0.2514, 0.6612] (45 deg boss), boss face centre [103.8655, -6.7688, 19.9833]; blade seated centre [10.899, 7.1625, 7.0023], axis [0.7734, -0.2306, -0.5905], face 0.1 below the mouth; gap 0.701 | radial pair, static r 19.996, boss top r 21.196, gap 0.4 | the blade face leans 49 deg (T15); consider stowed_x 102.5 |
| `MAGNET` | 6 x N52 block 6 x 2 x 2 + 6 x N52 disc 4 x 2 | 16 x D4 x 2 | 12 detent magnets built; the equator pairs are not |
| `SPINNER_ACCESS` | x [131.0, 148.354], 45 deg countersink to x 151.725, push-in tip cap | x [134, 165] | T02 |
| `HULL_L` | keep 165; built tip x 160.99 (R3 round), L/D 2.15 | | the cap restores the Myring tip |
| `JOINT` seal lead-in | 0.85 radial x 18 deg (mouth D65.79), cut after the pocket liners | 0.4 x 20 deg | T05 |
| `TOW_EYE` | teardrop outer profile parallel to the hole, hole apex 2 mm flat, R2 root fillet | | T06: crown 5.0 / 7.45 |
| `FIN` groove | groove_depth 4.8, groove_w_mean 4.73 (top 4.5, bottom 6.0, neck 3.3) | 4, 6 | T18 ligaments 1.2 / 1.2 |
| `FIN` retention | ISO 8734 2 x 10 cross pin at u 32, depth 2.2 | M3 grub in a heat-set insert | the insert does not fit the 10 mm pad |
| `FIN` root standoff | decision: the pad top stands about 7.1 mm off the hull, pads 9.2 mm proud | | T07 |
| `SHAFT.hub_flat_x` | [108.4, 128.0] | 111.6 | the D-flat runs out through the aft end so the hub slides on |
| `COUPLER` | claw A OD 32.4, x [46.5, 52.3]; claw B x [49.5, 62.5], cone x [60.5, 61.5]; claw A dogs to the OD | OD 22, x [46.5, 51.5] / [52.5, 63] | as built |
| `CARRIER` | cable window open to r 12.0 over theta 171 to 189 (reviewer value); chevron only, no text | notch bottom r 15.9 | as built, T23 |
| `ROTOR_R_TIP_TARGET` | 97.31 | 97.0 | measured on the deployed blades |
| `CLEAR_ROTATE` | 0.30 radial (slots, running) plus a new KNUCKLE_BORE_CLEAR 0.15 radial | ambiguous | decision |
| `INSERT_M3` | no stator inserts: M3 x 8 self-tap into D2.5 pilots | inserts | T21 |
| `part_split` | carrier aft face down, fins root down, fin blanks bottom down, cap spigot down, claw B membrane; see print-plan.md | carrier front face down, fins flat | as built |
| layout notes | pocket_depth.floor_r is a minimum on the root-notch ramp (CAD floor re-limited to 44.5 deg, up to 0.156 mm deeper) | | |
| `BAYONET` detent | open: a 0.3 mm ramp bump before the lock | | T17 |
