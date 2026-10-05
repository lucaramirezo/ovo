# CAD: parts, assembly and checks

Status: v1 CAD, integrated 2026-10-04; design review findings T01 to T31 applied 2026-10-05
(section 9); polish round 2026-10-05 (section 9.1: recessed fin rails, stowed detent, bayonet click,
blade skin parameterisation, mesh quality, all-ivory colours, display blanking plug, exploded view,
viewer controls, blade print cradle, closing-path check, output pruning); params write-back the same
day (section 10); review 2 round the same day (section 9.2: bayonet-click breakaway ties, one detent
boss top, tight skin-normal spinner windows, blended gland boss, pocket tips that follow the blade tip,
fin root flare, strip groove, finer blade cradle, coupons, exploded view, studio backdrop); fresh
integration build 2026-10-05 (`build.py --fresh`, then the fins, the shells' exports and `--only-assembly`
once more; every number below).
Sources: `ovo/cad/*.py` (one module per subassembly, see the table below) and
`ovo/cad/assembly.py` (integration). Every dimension comes from `ovo/params.py` and the layout
tables in `data/`: the as-built deviations of the review and polish rounds were written back into
params.py (section 10), so params and the CAD agree. The CAD numbers below are measured on the
real solids and come from `out/cad/checks/*.json`.

Regenerate everything with one command (logs in `out/work/build/`). The integration build
`--fresh` of 2026-10-05 (review 2) took 2554 s: shells 1242, blade 256, rotor 223, drivetrain 20,
fins 56, coupon 89, assembly checks and exports 581, renders 85; `--only-assembly` takes about 630 s.

```
.venv/bin/python scripts/build.py                  # modules, assembly, renders, summary table
.venv/bin/python scripts/build.py --fresh          # also rebuild the cached shell and part BREPs and
                                                   # prune outputs the current parts no longer produce
.venv/bin/python scripts/build.py --only-assembly  # reuse the modules' outputs, rerun the assembly
.venv/bin/python scripts/build.py --summary        # only print the table from the existing checks
```

`--fresh` ends by deleting every file in `out/cad/step`, `out/cad/stl`, `out/cad/checks`, `out/web`
and `out/renders` that the build did not produce (`build.prune_outputs`: removed parts such as
`inserts_m3_stator.step`, test artefacts such as `out/web/test.glb` and `manifest_test.json`, and
the `hydroegg_*` files from before the rename). It never prunes without a part cache.

`build.py` runs the six part modules through `scripts/cadrun.sh` (at most 2 OCC jobs on the
machine), then `python -m ovo.cad.assembly --force` (rebuilds the part cache under
`out/work/cad/integrate/cache/`), then the seven renders through the viewer, and prints a summary
table (parts, interference, key clearances, spin sweep, viewer check, print plan, outputs).

## 1. Outputs

| output | what |
|---|---|
| `out/cad/step/<part>.step` | every part in its assembled (world) pose, except `blade.step` and the coupon parts, which are in print pose |
| `out/cad/stl/<part>.stl` | every printed part in print orientation (on z = 0, mm); the three blades share `blade.stl`, the two fins and the two fin blanks are identical pairs; every STL is watertight with 0 folded sliver pairs (the assembly's STL audit rewrites any that is not, section 5). The aft shells' STLs carry the bayonet-click breakaway ties (print only, review C1 / P1) |
| `out/cad/stl/blade_cradle.stl` | the blade's sacrificial print cradle, in the coordinates of `blade.stl` (print one per blade, T13) |
| `out/cad/step/ovo_assembly_stowed.step` | the machine as carried: PACKED state (blades stowed, fins off, blanking strips in, blanking plug), 57 named parts, 30.1 MB |
| `out/cad/step/ovo_assembly_deployed.step` | the machine running: blades deployed, fins fitted, blanking plug, 57 named parts, 32.0 MB |
| `out/web/ovo.glb`, `out/web/manifest.json` | viewer scene: 67 meshes, 277k triangles, 9.3 MB, 0 folded pairs; roles, rotor flags, `packed` and `v2` flags, explode vectors, blade pivots, cut-face section meshes, material table |
| `out/web/packed_kit.glb`, `out/web/packed_kit.json` | a separate light scene: the packed egg with the two fins (and their magnets) lying flat on the ground in front of the nose (29 meshes) |
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
  blanking strips `fin_blank_0/1` in the recessed rails: the carry configuration). Fins are fitted in
  every state except PACKED. A part dict may carry `packed`: `False` = removed for packing (the
  fins), `True` = only in the packed egg (the blanks); the viewer follows the same flag (`packed=1`).
* **Variants**: a part dict may carry `v2`. `True` marks v2 / cutaway hardware (the M12 cable gland):
  only in CUTAWAY and the viewer's v2 mode. `False` marks its display-build stand-in (the blanking
  plug): everywhere else. The STEP assemblies carry the plug.

## 3. Part list

72 parts: 62 from the modules (one of them the hidden brass seal sleeve, a BOM item only, so the
part audit lists 71) plus 9 detent magnets that no builder modelled (added by
`assembly.magnets()` from the blade's seat geometry and the shells' boss geometry) and the
sector-cut O-ring for the cutaway (`assembly.oring_cutaway()`). The coupon parts are `hidden` (not
in the assembly; the brass seal sleeve is a BOM item only); the cutaway shells and the cut O-ring
replace the full ones in CUTAWAY only. Every part is a valid solid; the compounds are intended
(several identical bought items in one part, or a bearing's two seals). Volumes are B-rep volumes.

| part | module | role | rotor | kind | cm3 | notes |
|---|---|---|---|---|---|---|
| `front_shell` | shells | shell | | printed | 45.50 | nose, teardrop tow eye, gland seat (D18 boss without rim chamfer, blended into the nose by an R6 root fillet), carrier ring and twist-lock groove, bayonet socket with the 0.85 x 18 deg O-ring lead-in, 2 recessed half-dovetail fin rails sunk into bosses on the inner wall (no proud pads; each boss widened toward its nearer pocket liner), groove-floor magnet seats, pocket tips that follow the blade tip |
| `aft_shell` | shells | shell | | printed | 57.60 | spigot with O-ring groove and 3 bayonet lugs, bearing tube on 3 ribs, 3 blade pockets (smooth rim, 0.6 mm chamfer), 3 stowed-detent bosses leaned 49 deg, each started by a chin on the root-notch ramp (boss and chin share one top), 3 bayonet click tines in the spigot wall (0.81 mm underside slots; the STL carries 3 breakaway ties), static tail end x 108 |
| `front_shell_cutaway`, `aft_shell_cutaway` | shells | shell_cutaway | | printed (optional) | 37.76, 43.49 | sector theta -73 to 17 removed (faces the hero / iso cameras) |
| `hub` | rotor | body | yes | printed | 14.38 | D-bore on the shaft flat, 3 slotted ear pairs, stop faces with glued block-magnet seats, flat front face; ivory |
| `spinner` | rotor | body | yes | printed | 15.14 | 1.25 mm shell over the hub, one tight window per blade cut normal to the skin (review V1), M3 screw seat, D6 access bore with a 45 deg countersink at x 148.35 to 151.73; ivory |
| `spinner_cap` | rotor | body | yes | printed | 1.09 | push-in tip cap: D5.7 x 12 spigot with 3 crush ribs, 45 deg shoulder, R3 Myring tip at x 160.99, fingernail groove at the joint; ivory |
| `hinge_pin_0..2` | rotor | internal | yes | bought | 0.155 each | ISO 8734 3 x 22, pressed into the ears after each blade is in its slot |
| `pin_filler_0a..2b` | rotor | internal | yes | bought | 0.027 each | ISO 8734 3 x 4, loose in the pin channels (pin retention) |
| `hub_magnet_0..2` | rotor | internal | yes | bought | 0.024 each | N52 block 6 x 2 x 2, glued in the stop-face seats (6.1 x 2.1 x 2.2) |
| `hub_grub`, `spinner_screw` | rotor | internal | yes | bought | 0.024, 0.119 | M3 x 4 cup point on the D-flat; M3 x 8 into the shaft end |
| `blade_0..2` | blade | body | yes | printed | 5.47 each | 3 identical blades, curved pitch axis, centripetal U parameter (no face split), knuckle with D3.3 bore, heel stop with a block seat, shank disc seat; orange |
| `blade_cradle` | blade | hidden | | printed (x 3) | 0.26 (2 rails) | sacrificial print cradle in the blade's print coordinates (T13), stations every 0.25 mm (P6); `print_with: blade` |
| `shaft` | drivetrain | internal | yes | bought (machined) | 3.69 | 8 x 75.5, x 52.5 to 128, hub D-flat from x 108.4 out through the end, collar flat, M3 tapped end |
| `shaft_collar`, `shim_8x12x0p5` | drivetrain | internal | yes | bought | 1.18, 0.03 | thrust stop on the 608 inner race |
| `claw_a`, `claw_b` | drivetrain | internal | yes | printed | 2.68, 2.20 | dog coupling, 2.8 mm engagement, 0.3 mm flank gap, 1 mm float; claw A dogs run to the flange OD |
| `screws_claw_a_m3x4` | drivetrain | internal | yes | bought | 0.18 (2) | claw A to the generator rotor face |
| `bearing_608` (+ `_seals`) | drivetrain | internal | | bought | 2.09 | 608-2RS stainless, x 63 to 70, locating |
| `bearing_688` (+ `_seals`) | drivetrain | internal | | bought | 0.67 | 688-2RS, x 96 to 101, floating |
| `lip_seal_tc8x16x5` (+ `_spring`) | drivetrain | internal | | bought | 0.54 | v2 shaft seal, x 102.5 to 107.5 |
| `oring_141` | drivetrain | shell | | bought | 1.02 | AS568 -141, modelled installed (squeezed section); shell role so it hides with the shells |
| `oring_141_cutaway` | assembly | shell_cutaway | | bought | 0.77 | the same ring with the cutaway sector removed (CUTAWAY only) |
| `motor_carrier`, `screws_carrier_m3x8` | drivetrain | internal | | printed / bought (4) | 5.27, 0.31 | twist-lock bulkhead for the generator, D6.0 counterbores; the screws self-tap into the replica stator |
| `adapter_spacer` | drivetrain | internal | | printed | 3.16 | GL40-specific spacer |
| `generator_stator`, `generator_bell` | drivetrain | internal | / yes | printed | 15.21, 8.35 | display stand-in for the CubeMars GL40 (real envelope); the bell turns on a journal; the stator has D2.5 pilots, no inserts, and D1.6 dimples on PCD 39 for the real M2.5 holes |
| `inserts_m3_bell` | drivetrain | internal | yes | bought (2) | 0.08 | heat-set inserts in the replica bell (display only) |
| `cable_gland_m12` | drivetrain | body, `v2: True` | | bought (v2) | 2.98 | Hummel HSK-K M12x1.5 in the nose seat; CUTAWAY and the viewer's v2 mode only |
| `gland_plug` | drivetrain | body, `v2: False` | | printed | 1.52 | display M12 blanking plug: low dome covering the whole D18 seat flush with its edge, push-in crush ribs, coin slot; ivory |
| `fin_0`, `fin_1` | fins | body | | printed | 5.60 each | Clark Y fins on recessed half-dovetail rails (no proud pad), 45 deg root flare beside the neck, tip corners on a short land, magnet retention; ivory; `packed: False` |
| `fin_blank_0`, `fin_blank_1` | fins | body | | printed | 0.168 each | blanking strips, outer face on the hull profile (flush), fingernail groove near the aft end, entry wedge trimmed at 0.6 mm, same magnet pair; ivory; `packed: True` |
| `rail_magnet_0/1`, `fin_magnet_0/1`, `fin_blank_magnet_0/1` | fins | body | | bought | 0.016 each | N52 4 x 2 x 2 blocks: groove floor, fin tongue, strip (`packed` as their host) |
| `tolerance_coupon` (+ 3 keys) | coupon | hidden | | printed | 23.10 (+ 0.26, 0.52, 0.70) | print first: seats, pins, insert, M12 tap drill, bayonet and dovetail keys; the standing key `tolerance_coupon_key_tongue_up` prints as the fin |
| `tolerance_coupon_b` | coupon | hidden | | printed | 18.40 | print with coupon A: hub D-bore, block and disc magnet seats (vertical, horizontal, and tilted 49 deg on the sloped +X end), upright ear with D2.9 / 3.0 / 3.1 holes |
| `blade_magnet_heel_0..2` | assembly | body | yes | bought | 0.024 each | N52 block 6 x 2 x 2 glued in the heel stop face, pivots with its blade |
| `blade_magnet_stow_0..2` | assembly | body | yes | bought | 0.025 each | N52 disc 4 x 2 in the shank inner face (stowed detent), 0.1 mm below the seat mouth, pivots with its blade |
| `tail_magnet_0..2` | assembly | body (`_0`: shell) | | bought | 0.025 each | N52 disc 4 x 2 in the aft shell's stowed-detent bosses on the blade face normal (49 deg from radial); `tail_magnet_0` takes the shell role because its boss is in the cutaway sector |

Colours (owner decision 2026-10-05): every part on the outside of the egg is ivory (shells, hub,
spinner, spinner cap, fins, blanking strips, blanking plug); the blades are orange; the generator
replica is blue; the printed drivetrain internals are shown charcoal (`drive`).

## 4. How the parts fit

**Shells and joint.** The egg splits at the equator x 57.5. The aft shell's spigot (seal zone OD
63.89) enters the front shell's socket (bore D64.09, 0.1 mm radial) and carries the -141 O-ring in
its groove. The socket mouth is a 0.85 mm x 18 deg lead-in (mouth r 32.879 at every azimuth, the
stretched ring stands to r 32.528). Three lugs pass the socket lip through entry gaps and lock
after a 35 deg twist; the channels end exactly at the locked lugs, so the hard stop aligns the
pockets across the joint (measured: 0 mm3 overlap, lugs on the stop at 0.000 mm). The lock clicks
(T17, `BAYONET.click`): three curved flexure tines in the spigot wall between the lugs (13.2 x 2.4 x
2.36 mm) each carry a 0.5 mm bump that rides 0.25 mm under the socket ring during the twist and drops
into a dimple at the lock angles 120 / 240 / 0 (measured: 0 mm3 locked, 1.82 mm3 interference along
the twist; ESTIMATE, PLA E 3.5 GPa: 2.97 N per tine, 17.7 MPa at the tine root). Each tine's 45 deg
underside slot is 0.81 mm normal (1.15 along x); the aft shell's print solid carries a 0.8 mm
breakaway tie across each tine's free-end slit at the tine's bottom (review C1 / P1), so the tine's
first layers bridge from its root to the wall beyond the slit. Cut the ties through the slits after
printing; the world solid and the checks have no tie.

**Generator bay (front shell).** The motor carrier enters through the socket with its three tabs
at theta 120 / 240 / 0, passes the ring notches and twists 30 deg toward -theta to lock at 90 / 210 /
330; the generator reaction torque drives the tabs into the stop. Carrier, adapter spacer and the
replica stator stack on 4 M3 x 8 that self-tap into D2.5 pilots (1.65 mm from the bell journal
groove; the real GL40 has M3 threads; its M2.5 holes on PCD 40 are D1.6 dimples on PCD 39 on the replica); the bell turns on the stator journal (0.298 mm min gap over
a full turn) and carries claw A on 2 M3 x 4. In v2 the cable leaves through the carrier's theta 180
notch to the M12 gland in the nose seat (gland to carrier 3.61 mm). The display build closes the
seat with the printed blanking plug (`gland_plug`): a low ivory dome that covers the whole D18 seat
flush with its edge (the seat has no rim chamfer and the boss blends into the nose through an R6
fillet, review V2), pushed into the D10.5 tap-drill bore on three crush ribs (0.192 mm3 of rib press,
seat face 0.000 mm), with a coin slot to twist it out.

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
3 x 4 dowel in each pin channel and the spinner keep a loosened pin in both ears (0.34 mm free
travel, 3.36 mm far-ear engagement left). The spinner slides on over the hub with the blades stowed
(1.07 mm minimum; also from fold 30 down to stowed; its tight windows block the deployed roots,
review V1), an M3 x 8 clamps spinner and hub to the shaft end, and the tip cap pushes into the access
bore on its crush ribs. The windows are the blade root's swept footprint in the spinner's material +
RUN_CLEAR + 0.7, cut normal to the skin, open at the front edge, with a rounded top at x 127.0
(1.0 mm of the cone ring left under the seat ceiling). Deployed, each blade's heel lands flat on its hub stop face (0.000 mm, 0 mm3), where the heel
and hub block magnets face each other across 0.400 mm (about 2.32 N).

**Stowed blades.** Folded 90 deg at the stowing index, each blade lies in its pocket: the aero span
0.337 mm above the pocket floors, the root 0.899 mm clear of the static tail outside the 8.25 mm
detent rest zone, the shank's inner face 0.171 mm above its leaned (49 deg) detent boss, and
1.10 mm from the pocket rim at the skin (the rim is the blade's swept footprint over the whole fold
plus 1.0 mm; over the tip zone the pocket follows the CAD blade's rounded tip + 2 mm, review V4). The
stowed magnet pair sits 0.420 mm apart face to face (0.365 mm at the closest edge), 1.66 deg apart:
about 2.27 N against the 2.32 N designed, within the owner's 0.45 mm / 2 deg / 2.0 N. The fins stay
5.635 mm from the blades over the whole fold.

**Fins** (`fins.py`, recessed rails, owner decision). Each Clark Y fin's root tongue slides aft
and inward along a half-dovetail groove. The groove is sunk into a boss on the inside of the front
shell wall, so no pad stands proud of the hull (measured: 0 mm3 of front shell more than 0.05 mm
outside the hull at the rails). The fin is held by an N52 4 x 2 x 2 block in the tongue over a
second one in the groove floor. For packing, pull each fin forward out of its slot and slide the
blanking strip in: the strip is the tongue section with its outer face on the hull profile, flush
when seated, held by the same magnet pair; to take it out, hook a fingernail in the groove across its
outer face near the aft end and drag it forward (review C6). Stowed, the egg keeps its 75 mm diameter.
Over the slot the fin's root carries a 45 deg flare from the neck edges up to the airfoil (review C4 /
P3), so the tip-down print has no ledge. Details and numbers: `out/cad/checks/fins.json` and
`shells.json`.

### Assembly order

From `drivetrain.json` (assembly, 12 steps) and the rotor module:

1. Press the 2 M3 inserts into the replica bell (the stator has D2.5 pilots, no inserts).
2. Stack carrier, adapter spacer and stator with 4 ISO 7380 M3 x 8 from the carrier counterbores,
   self-tapping into the pilots.
3. Slide the bell onto the stator journal; bolt claw A to the bell with 2 ISO 4762 M3 x 4.
4. Display build: push the blanking plug into the nose seat (coin slot outward). v2: tap the seat
   M12x1.5, screw the gland in, pass the cable and lay it in the carrier notch.
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
    dowels; the spinner (blades stowed; face its D22.6 clamp ring flat first, review P8); M3 x 8 into
    the shaft end; push the tip cap in (the fingernail groove pulls it out again for the screw).
11. Grease the O-ring and roll it into the spigot groove.
12. Take the fins off and unfold every blade by at least 10 deg (fold 80 deg or less; the sweep
    closes cleanly from fold 85 down and collides at 90); with the fins left on, unfold by at least
    45 deg (fold 45 or less; the sweep is clean from fold 50 down and overlaps 2.0 mm3 at 60).
    Then join the shells: the stowed blades lie in pockets of both shells, 35 deg apart at the
    bayonet entry. The limits come from `assembly.closing_sweep` on the real solids (entry, push,
    twist, both senses; section 5), with one sweep step of margin. Push the shells together at the
    bayonet entry marks (the claws find each other on their lead-ins) and twist 35 deg until the
    bayonet clicks on the hard stop; only then fold the blades home and refit the fins or the
    blanking strips. Nothing but the hub on the static tail stops the shaft toward -X: after
    stowing, pull the spinner aft so the hub runs free. `assembly.json` (`assembly_steps`)
    restates this step with the measured numbers.

Packing: pull each fin forward out of its rail (the magnet pair lets go) and slide its blanking
strip in to the stop (a fingernail in its groove drags it out again). Before stowing, twist the joint until the bayonet clicks into its lock (the
hard stop aligns the pockets across the joint).

## 5. Interference and clearance results

Method: manifold3d booleans and `min_gap` on meshes of the real solids (0.02 mm chord, shells
0.03 mm; every shape's stored triangulation is dropped first so the tolerance is honoured). Exact
OCC booleans against the full shells take minutes each. Overlaps up to 0.5 mm3 on curved contacts
are tessellation facets; flat contacts tessellate exactly. Every pair of parts present in a state
is tested whose boxes come within 3 mm.

| state | fold deg | parts | pairs within 3 mm | intended overlaps | zero-gap contacts | defects (> 0.5 mm3, unexpected) |
|---|---|---|---|---|---|---|
| stowed | 90.0 | 57 | 262 | 13 | 36 | 0 |
| deployed | 0.0 | 57 | 57 | 0 | 9 | 0 |
| fold_30 | 30.0 | 57 | 69 | 0 | 6 | 0 |
| fold_60 | 60.0 | 57 | 98 | 0 | 6 | 0 |
| cutaway | 0.0 | 56 | 107 | 2 | 17 | 0 |
| packed | 90.0 | 57 | 14 | 0 | 4 | 0 |

(Static pairs are computed once in STOWED; the other states add the pairs that involve a folding
part, CUTAWAY the pairs that involve a cutaway part or the gland (v2), PACKED the pairs that involve
a blanking strip. In CUTAWAY the gland replaces the plug; it has one part fewer because
`tail_magnet_0`'s boss lies in the removed sector, so that magnet goes with the full shell.)

**Intended overlaps** (listed in `assembly.EXPECTED`, reported separately):

| pair | kind | overlap mm3 | state | why |
|---|---|---|---|---|
| front_shell / gland_plug | press | 0.1925 | stowed | display blanking plug: 3 crush ribs (tips D10.6) in the D10.5 bore |
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

The O-ring, the lip seal, the glued magnets and the fin magnets touch or sit in their seats without
overlap.

**Zero-gap contacts** (designed): shells at the lug stop; carrier tabs on their stop; fin tongues
and blanking strips on the rail floors and end stops; the plug's band on the gland seat; hub,
spinner, shaft end and spinner screw clamped together; grub on the D-flat; the stacked carrier /
spacer / stator and claw A / bell / inserts faces; collar / shim / 608 inner race; magnets on their
seat floors; pin fillers on the pin ends; deployed heels on the stop faces.

**Key clearances on the real solids** (compared with `data/layout.json` and the design rules):

|  | check (real solids) | measured | layout / design | limit |
|---|---|---|---|---|
| PASS | hub gap: static tail end (aft shell) to hub front | 0.600 mm | 0.6 | >= 0.6 |
| PASS | spinner to aft shell | 0.600 mm |  | >= 0.6 |
| PASS | lip seal outer face to hub front (axial) | 1.100 mm | 1.1 | >= 0.6 |
| PASS | shaft to hub D-bore (radial) | 0.074 mm, 0 mm3 | 0.075 | > 0 |
| PASS | stowed aero span (x < 100) to the pocket floors and walls, 3 blades | 0.337 mm (blades 0.337, 0.337, 0.337) | 0.3 | >= 0.3 (-0.05 tessellation) |
| PASS | stowed root (x > 100, more than 8.25 mm from the detent boss face centre) to the static tail | 0.899 mm (blades 0.899, 0.899, 0.899) | 0.8 | >= 0.8 |
| PASS | stowed blade to its detent boss (blade within 8.25 mm of the boss face centre) | 0.171 mm (blades 0.171, 0.171, 0.171) |  | > 0 (shells design 0.1 to the layout surface) |
| PASS | stowed detent magnet pair (blade shank / tail boss): face-centre gap, min gap, axis angle, lateral | 0.420 / 0.365 mm, 1.66 deg, 0.00 mm; about 2.27 N (layout 2.32 N at 0.4 mm) | 0.4 | owner call T15: face gap <= 0.45 mm, axis angle <= 2.0 deg, force >= 2.0 N |
| PASS | deployed detent magnet pair (blade heel / hub stop face): min gap | 0.400 mm, about 2.32 N | 0.4 | 0.4 +- 0.15 (seat depth H + 0.1: 0.3 is slightly stronger) |
| PASS | deployed heel on the hub stop face: gap / overlap (3 blades) | 0.000 mm / 0.000 mm3 | 0.0 | contact, no overlap |
| PASS | blades folding 0..90 deg at the stowing index to the fins | 5.635 mm at fold 90 | 8.44 | >= 2.0 |
| PASS | front shell material proud of the hull at the fin rails (recessed rail, owner decision) | 0.000 mm3 more than 0.05 mm outside the hull |  | flush egg (75 mm stowed) |
| PASS | blades folding 6..90 deg to the hub (run clearance; below 6 deg the heel lifts off its stop) | 0.296 mm at fold 6 |  | >= RUN_CLEAR 0.3 (-0.01 tessellation) |
| PASS | blades folding 0..90 deg to the spinner | 1.043 mm at fold 0 |  | >= RUN_CLEAR 0.3 |
| PASS | blades folding 0..90 deg to the hinge pins (knuckle bore) | 0.149 mm at fold 60 | 0.15 | > 0.1 (bore 3.3 on a 3.0 pin) |
| PASS | blades folding 0..90 deg to the shells (minimum is the stowed boss) | front 0.337 mm at fold 90, aft 0.171 mm at fold 90 |  | > 0 |
| PASS | neighbour blades folding together | >= 3 mm at every fold |  | >= 1.0 (layout rule) |
| PASS | claw B (rotating) in the static tube mouth | 1.262 mm | 0.95 | >= 0.3 |
| PASS | claw A to claw B dog flanks | 0.300 mm | 0.3 | >= 0.25 |
| PASS | display blanking plug in the gland seat: rib press overlap, seat face | 0.192 mm3, 0.000 mm |  | crush ribs only (< 3 mm3) |
| PASS | cable gland to motor carrier | 3.609 mm | 2.275 | >= 1.0 |
| PASS | generator replica to the shells | 3.122 mm | 0.637 | >= 0.5 (bay to shell) |
| PASS | front / aft shell locked: overlap, gap | 0.0000 mm3, 0.000 mm (lugs on the hard stop) | 0.0 | no overlap |
| PASS | closing path (bayonet entry, push, twist) vs the blades | fins off: closes at fold <= 85 deg (unfold >= 5 deg); fins on: fold <= 50 deg (unfold >= 40 deg) |  | 0 mm3 overlap |

spin sweep: 23 pairs, 0 rubbing; [(0.2981, 'generator_bell', 'generator_stator'), (0.389, 'shim_8x12x0p5', 'bearing_608_seals'), (0.5, 'claw_b', 'bearing_608'), (0.5, 'shaft_collar', 'bearing_608'), (0.6, 'hub', 'aft_shell'), (0.6, 'spinner', 'aft_shell')]

The fold sweep runs in 3 deg steps; hub, spinner, pins and shells are 3-fold symmetric, so blade 0
stands for all three there, while all three blades are checked against the fins. Below 6 deg the
heel lifts off its stop face and its gap grows from zero, so the hub run clearance is read from 6 deg
(as in `rotor.py`). The detent rest zone is `DETENTS.stowed_rest_zone_r` = 8.25 mm around the boss
face centre (chin length + half width + 0.5): the boss and chin share one top that follows the blade
face 0.1 mm under it (layout surface) and never falls toward the bed (review C2 / P2), so that zone
is the boss check (0.171 mm on the CAD blade), not root clearance.

The closing path (recheck R3, `assembly.closing_sweep`) runs as follows. The front shell enters at
the bayonet entry, in both twist senses. It is pushed on from 40 mm at 8 stations, then twisted to
the lock in 2.5 deg steps. At every step it is checked against the three blades at folds 0 to 90.

* Fins off: 0 mm3 at every fold up to 85 deg, and 1 938 mm3 at 90.
* Fins on, measured on the front shell with both fins and their magnets: 0 mm3 up to 50 deg, then
  2.0 mm3 at 60, 15.3 at 65, 36.7 at 70, 75.6 at 78 and 107 at 85.

Step 12 of the assembly order follows from this with one sweep step of margin
(`assembly.assembly_steps`): fins off, fold 80 deg or less (unfold at least 10 deg); fins on, fold
45 deg or less (unfold at least 45 deg).

**Rotor spin sweep** (deployed, rotating parts turned 0..350 deg in 10 deg steps against every
static part within 3 mm): 23 pairs, none rubbing. Closest running gaps: generator bell to stator
0.298 mm, shim to the 608 seals 0.389 mm, collar and claw B to the 608 0.5 mm, hub and spinner to
the static tail 0.6 mm. The only zero gaps are the bearing inner races and the seal lip on the
shaft, which turn with it.

**Mesh quality** (recheck R1): every STL is watertight with consistent winding and 0 folded sliver
pairs (adjacent normals more than 170 deg apart; `stl[].folded_pairs_gt170`). The modules'
`aft_shell.stl` (15) and `aft_shell_cutaway.stl` (10) had some, so the audit rewrote them through
`export_web.write_clean_stl` (`stl[].regenerated` records it). The GLB has 0 in all 67 meshes,
measured after reloading the float32 file.

**Module checks** (each module's JSON): rotor (hub seat walls 1.27 / 1.26 / 3.08 mm, hub min wall
1.220 mm, spinner min wall 1.239 mm (p01 1.243) on the 1.25 shell; windows (V1): aft end x 127.0, a 12.2 mm
bridge at the top, 1.0 mm of cone ring left under the seat ceiling; blade roots to the spinner 0.98 mm
over the fold (was 0.69); slide-on over the stowed blades 1.07 mm (also from fold 30; fold 0 and 15 are
blocked by the tight windows); cap spigot 2.2 mm clear of the M3 head; pins: 0.34 mm free travel,
3.36 mm far-ear engagement left).

Blade:

* block seat walls 1.15 (heel tip) / 1.95 (flanks)
* trailing edge 0.600 mm
* skin end check 55.9 deg worst, 0 non-monotone row ends
* classification 0 / 1000 far and 0 / 1500 near-surface points wrong
* STL 12 042 triangles, watertight, 0 folded pairs
* cradle (P6): 328 stations every 0.25 mm in 2 rails, up to 2.77 mm tall, 0.200 mm to the blade,
  0 mm3 overlap; it carries 71 % of the 126 mm2 of blade faces beyond 45 deg (was 57 %); the top sits
  0.20 to 0.33 mm under the blade along the span (`te_gap.py`, every 0.5 mm)

Drivetrain: claw skin 1.2 over the counterbores, stator pilots 1.65 mm from the journal groove, stator
dimples 1.19 mm from the chamfered rim in the first layer (P9), carrier web 1.01 mm between the D6.0
counterbores and the D16 bore in every layer (P9, was 0.56), claw B bridge on a 0.2 mm membrane, plug
print overhang 12.6 mm2 (the slot-roof bridge).

Shells (42 checks, all pass):

* stowed blades: aero span 0.308 mm over the CAD pocket floors (layout surface; 0.337 on the CAD
  blade), side walls 0.829, rim 1.1, flush 0.065; pocket forward end: -theta bound at the front row
  -8.4 deg (the square table tip had -18.0), 19 tip rows follow the CAD tip (V4)
* stowed detent: boss axis on the blade face normal (lean 49.03 deg, 0.0 deg mismatch), boss + chin
  top 0.130 mm under the layout blade surface, its faces lean at most 49.01 deg (grid), the static
  magnet face stands up to 0.20 mm proud of the boss top at one side of the hole; pair 0.413 mm /
  2.29 N on the layout surface
* walls under the pocket floors 2.135 (front) / 2.069 (aft); generator bay to the shells 0.637
* fin rail groove walls (OCC rays, both fins): floor 1.552, magnet seat bottom 1.297, flat wall
  1.445, neck 5.342 (theta 87.5) / 2.9 (theta 267.5, its boss now widens on -w, C7), undercut lip
  2.355, end stop 2.634; slot to the pocket rims 5.809 (arc), grooves to the pocket voids 4.482; front
  shell proud of the hull at the fins 0.0 mm
* joint: 0 mm3 locked, lugs on the stop; O-ring lead-in mouth r 32.879 (ring 32.528); bayonet
  click 0 mm3 locked, 1.818 mm3 along the twist, tine 2.97 N / 17.7 MPa (ESTIMATE); underside slots
  0.813 mm normal (`tine_gap.py`: 0.800 to 0.818), 3 breakaway ties 4.41 mm3 in the print solid
* nose: gland boss root fillet R6 built; tow eye crown 5.0 / 7.45 mm, utilisation 0.41
  PETG; gland thread engagement 6.06 mm
* layer scans of the four shell STLs at the 0.1, 0.15 and 0.05 phases (0.2 mm layers): 0 islands

Fins (all pass): chord 30.0, span 70.0, incidence 12.5 deg (probes); root 0.101 to 0.102 mm over the
skin; fin seated 0.013 mm3 (end-stop contact facets) and strip 0 mm3, both pulled out along the rail in
0.5 mm steps over 20.1 mm with nothing larger; strip proud of the hull 0.0 mm; magnet pair 2.16 / 2.36 N
normal (model / seated), about 1.34 N to slide out; tongue and strip seat ligaments 1.30 mm; GProp
5606.0 / 5607.7 / 5599.2 mm3; couple ratio as built 1.1422 (band 1.10 to 1.155; but see the
yaw play, section 8); print: 0.66 mm2 of the fin beyond 45.5 deg (the first 0.2 mm of the slot entry;
was 25 mm2), the strip's only one its magnet seat's 2.1 mm ceiling bridge.

## 6. Web scene and viewer check

`assembly.write_web()` tessellates each part for the web through `export_web.to_trimesh` and writes
the GLB and manifest through `export_web.write_scene`.

* **Deflection**: absolute chord deflection, 0.05 mm / 0.25 rad, blades 0.03 mm / 0.2 rad
  (`export_web.occ_mesh`). build123d's own `tessellate` scales the deflection with each edge's
  length, so on long faces only the angular limit applied; the blade's trailing edge showed
  facets.
* **Mesh clean-up** (recheck R1): every mesh goes through `export_web.clean_mesh`, which welds in
  float64, collapses needles and flips caps until no adjacent pair of faces folds more than
  170 deg. The mesh stays watertight, and the result is kept only when it is better than the
  input.
* **Normals** are smooth across edges flatter than 28 deg and AREA-weighted. trimesh's default
  angle weighting let the 180 deg corner of a sliver drag its vertex normal, which caused the
  streaks.
* **GLB check**: `assembly.glb_quality` reloads the written GLB (float32 positions) and counts
  folded pairs per mesh. The counts are in `assembly.json` under `web.mesh_quality`.

Roles: shells and the O-ring `shell`, cutaway shells and the cut O-ring `shell_cutaway`, everything
inside `internal` (shown in cutaway and exploded views), visible bought parts and the detent
magnets `body`, coupon parts and the blade cradle `hidden` (left out of the GLB). The cut faces of
the cutaway shells are separate meshes (`*_cutaway_section`, colour #C98B5E, lifted 0.03 mm).

The manifest carries two variant flags:

* `packed`: `false` for the fins and their magnets, `true` for the blanking strips and theirs.
* `v2`: `true` for the M12 cable gland, which shows with `cutaway=1` or `v2=1`; `false` for the
  display blanking plug, which shows otherwise.

Materials (owner decision): shells, hub, spinner, spinner cap, fins, blanking strips and the
blanking plug are all ivory #ECE7DD (`shell`; the `hub` and `fin` keys are copies of it). The
blades are orange #E8622C. The generator replica stays blue. The printed drivetrain internals
(claws, carrier, spacer) take the new `drive` key (charcoal), so the cutaway keeps its contrast.

**Viewer** (`viewer/index.html`). URL parameters are unchanged, with these additions:

* `v2=0|1`: show the gland instead of the plug.
* `fit=tight&fill=<0..1>`: project the visible vertices for the chosen view, re-centre on their
  screen box and pull the camera in until the parts span `fill` of the frame (4 passes for
  perspective). The renders use it.
* `shadow=1` (with `fit=tight`): also frame the ground shadow, i.e. each fitted vertex projected
  onto the ground along the key light. The level side and front views use it; without it their
  shadow ran off the bottom of the frame.
* `ground=0`: no ground shadow. The exploded view uses it: a cloud of floating parts casts a
  scattered shadow that reads as clutter.
* `ui=1` adds **packed** and **v2 gland** checkboxes beside fold, explode, cutaway and spin. The
  controls start from the URL state.
* `theme=studio` (review V7): a mid warm-grey backdrop (#d9d4cb to #b8b2a8, ground shadow 0.22) for
  the portfolio renders; the ivory egg measured 1.10:1 against the light backdrop. `light` stays the
  default for the docs.
* The key light's shadow carries `normalBias = 0.5` (mm) next to its bias (review V6): it removes
  the self-shadow stipple on the orange faces.

The controls were driven headless (`out/work/cad/polish_present/ui_test.py`):

* packed on: the fins hide and the blanking strips show.
* v2 on: the plug hides and the gland shows.
* cutaway: the gland shows, with the cutaway shells.
* fold slider to 1: blade_0 deploys (box max z 97.3).

**Exploded view** (`assembly.EXPLODE_STATIONS`). One axial line in assembly order, in the deployed
pose:

1. front shell, with the fins sliding out along their rail path (the fins module's `explode`
   vector) and the plug or gland 28 mm out of the nose seat along its axis (16 before V3)
2. carrier screws, carrier, adapter spacer, stator, bell, claw A and its screws
3. claw B, collar, shim, 608, 688, O-ring
4. aft shell, with the shaft left in it: its aft end shows behind the tail, carrying the lip seal
5. hub, with the blades 60 mm radially out (34 before review V3: the deployed blades were drawn
   across the aft shell rim), the pins 13 mm out of their ears along the pin axis,
   the hub magnets out through the open slot and the grub out of its tap
6. spinner, spinner screw, cap

Each station is placed behind the previous one, its front face a set gap behind the previous
station's aft face, with the extents measured on the real parts. A geometry change can therefore
never make two stations collide.

Measured: exploded length 406.5 mm against the egg's 178.1 mm, so 2.28 x (`assembly.exploded_extent`,
target about 2.2 x). The render (review V3) uses view 20 / 18 with the rotor turned 30 deg (`spin=30`):
one horizontal line in build order with one blade straight up, one straight down and the third toward
the viewer, so the frame's lower half carries a blade instead of staying empty.

Browser check (`assembly.verify_viewer_fold`, headless Chrome through `scripts/render.py`'s
server): the viewer's read-only `window.partBoxes()` returns each part's posed vertex box in the
CAD frame. At fold = 1 every blade box matches `blade.blade_world(i, 0)` within 0.0105 mm (web
tessellation); at fold = 0 the boxes match the stowed solids.

## 7. Renders

`out/renders/` (1600 x 1000, `theme=studio`, `fit=tight`):

| render | viewer URL parameters | shows |
|---|---|---|
| `hero_deployed.png` | view -40 / 22, fold 1, fill 0.86 | running: blades deployed on their stops, fins fitted, blanking plug in the nose seat; the two fins read as a level pair across the egg (review V5) |
| `stowed.png` | view hero, fold 0, packed 1, fill 0.8 | carried: blades in their pockets, fins off, blanking strips in: the clean ivory egg with three orange stripes |
| `packed_kit.png` | packed_kit scene, fill 0.86 | the packed egg with the two loose fins lying flat on the ground in front of the nose (`write_kit_scene`: each fin rests on the largest flat facet of its convex hull, long side along the egg, lowest point on the egg's ground plane) |
| `cutaway.png` | view -40 / 22, cutaway, fold 1, fill 0.86 | generator, carrier, coupling, shaft, bearings, seal, the M12 gland (v2); copper section faces |
| `exploded.png` | view 20 / 18, spin 30, fold 1, explode 1, fill 0.92, ground 0 | every part on one horizontal line in assembly order, blades 60 mm radially out, up and down; no ground shadow |
| `side.png`, `front.png` | side / front, fold 1, fill 0.88 / 0.9, shadow 1 | profile and frontal area, with the whole ground shadow in frame |

All seven were inspected after the review 2 build (2026-10-05). On the studio backdrop the ivory
egg, its strips and the fins separate clearly. In `stowed.png` each stripe now runs into a smooth,
rounded spinner window cut normal to the skin (no saw-tooth outline; crop
`out/work/cad/fix3/zoom_stowed_tail_new.png`), the pocket's forward end follows the rounded blade tip
(no white triangle ahead of the stripe, `zoom_stowed_stripe_front_new.png`), and the gland boss and
plug read as one smooth bump blended into the nose (`zoom_stowed_nose_new.png`); it still stands 9 mm
off the skin (section 8). The hero and cutaway cameras show both fins as a level pair; the exploded
view fills the frame on one line with a blade above and below. Residual: the blanking strips keep
their sharp corners and entry-relief step (V8, section 8); through each window the hub's slot bands
still show as steps.

Blade close-ups (recheck R1, `out/work/cad/polish_present/blade_front_close.png`,
`blade_back_abs.png`, made with `blade.preview_parts` / `blade.preview_web`): the pressure and
suction faces shade smoothly, with no streak or triangular mark, and the trailing edge reads as
a straight crisp edge.

## 8. Open issues

1. **Stowed detent (blade + shells, T15, owner call built).** The aft-shell boss leans 49 deg on the
   blade face normal and starts from the root-notch ramp on a chin; boss and chin share one top
   (review C2 / P2), so no top face leans beyond the 49 deg boss plane and the layer scans find no
   island. The pair: 0.420 mm face to face (0.365 mm at the closest edge), 1.66 deg apart, about 2.27 N (layout 2.32 N), within the 0.45 mm / 2 deg / 2.0 N targets. The static magnet keeps its blade-driven depth, so its face stands up to
   0.20 mm proud of the boss top on one side of the hole (glued 1.8 to 2.0 mm deep; harmless). The blade
   module's `stowed_detent_interface` still compares with the layout's radial DETENTS pair (superseded),
   and `DETENTS.stowed_chin.clear_to_blade` (0.25) is no longer used (proposal).
2. **Rotor band behind the stripes (review V1, partly fixed).** The spinner windows are now tight
   and cut normal to the skin (section 9.2). Two of the three proposed changes were not built:
   * raising the hub's outer surface under the windows: with the tight outline, 88 to 90 % of each
     window is the blade root's swept slot (at r = R - 0.3), so a raised hub would fill about
     10 to 12 % of it; the slot's band steps stay visible in the window;
   * a shank fairing on r(x) in the stowed pose: it would sweep through the spinner band while
     folding (growing the hub slot and the window by its thickness) and over x 105 to 108 it would
     enter the aft-shell pocket rim, so it needs a coordinated blade / rotor / shell change.
   The reviewer's target (under 15 % open skin in x 108-128) is not reachable with folding blades:
   the deployed root crosses the spinner at x 118 to 127, and that window must stay open when the
   blades are stowed. Measured (packed display, `openarea.py`): x 100-108 21.7 %, x 108-118 46.1 % (was 57.1), x 118-128 36.7 % (was 50.5); the shells measure 7.5 to 9.9 %.
3. **Fins (owner decisions, v2 items).** The fins run on recessed rails with magnet retention and
   come off for packing; flush blanking strips keep the 75 mm egg. Not changed:
   * yaw play: the tongue has 0.25 mm side clearance over its 12.1 mm length, so the fin can yaw
     about +1.5 / -1.75 deg about the stop corners (review C5, `fin_yaw.py`); with a linear lift slope this is
     about +-10 % lift, so the couple ratio as built (1.142) is really a range of about 1.01 to 1.25 (ESTIMATE, linear lift slope, alpha_L0 -3.07 deg)
     against the 1.10 to 1.155 band. No effect on the display build; for v2, two or three crush
     ribs (0.05 interference) on the tongue's flat side near its forward end would locate it
     (proposal).
   * the slot mouth and strip outline keep their sharp corners and the entry-relief step (review
     V8, cosmetic): rounding them needs matching corners on the tongue and the strip and an S-blend
     of the relief that the fin's relief chamfer and the slide checks would have to follow.
   * the planform stays the constant-chord Clark Y plank (review V5 option 1, a tapered planform of
     equal area, changes the couple ratio: owner decision); the hero and cutaway cameras moved instead.
4. **Nose gland seat (review V2).** The display build keeps the raised v2 seat (spot face 6.79 above
   the nose, thread engagement 6.06 mm) and blends it into the nose with an R6 fillet, the plug flush
   on top; the bump stands 9.1 (crown; the boss 5.3) mm off the skin. Sinking the seat (a counterbored spot face in a
   thickened inner wall, GLAND boss_in_cavity True) would leave only the plug's dome proud, but it moves
   the v2 gland interface toward the carrier: owner decision (proposal).
5. **Blade print support (T13).** The blade prints on `blade_cradle.stl` (stations every 0.25 mm since
   review P6): it carries 71 % of the 126 mm2 of blade faces beyond 45 deg, its top 0.20 to 0.33 mm under the blade along the span; the rest are the bridged ceilings of the heel slot and the shank seat (inside the blade) and faces within one layer of the bed next to the knuckle. Tree supports (top Z 0.2) stay as the fallback.
6. **Nose-plug bridge (shells, T16).** The 38.6 mm flat bridge inside the nose is the only bridge; a
   45 deg cone from r 19.2 would reach the carrier (x 22.7). Good bridge settings, expect a little sag on
   that internal face.
7. **Joint click (T17).** The tine force and stress are an ESTIMATE from a cantilever model (2.97 N per
   tine, 17.7 MPa in PLA). The print solid's breakaway ties make the tines' first layers a bridge; they
   are cut through the slits after printing (a manual step, in the print plan). The coupon carries a
   bayonet key but no tine: judge the click on the first aft-shell print.
8. **Spinner seat ceiling (review P8).** The D22.6 clamp ceiling is a flat bridge pierced by the screw
   hole and the three peg holes (short cantilevers up to ~8 mm next to them). It is the hub clamp and the
   axial datum: face it flat on sandpaper over glass before assembly (print plan, step 10). A raised
   clamp ring with a 45 deg cone inside would remove the cantilevers but changes the T12 datum design.
9. **Small residual faces beyond 45 deg** (all intended, measured on the B-reps):
   * aft shell: the 29 mm2 seal stop ring; the detent boss and chin tops at the accepted 49.01 deg
     (3 x 32.3 mm2, owner call T15)
   * front shell: the nose bridge (1114 mm2) and the 12 mm2 eye-hole flat; 3 x 2.1 mm2 at 54 to 57 deg
     on the pocket walls at x 35.6 to 36.4, where the smooth pocket loft overshoots the tip chamfer's
     45 deg narrowing (review V4; a ruled tip loft was tried and was worse)
   * spinner: the D22.6 seat ceiling (item 8) and the three window tops: rounded arches whose last
     0.6 mm overhang as crescents up to 0.9 mm wide before a ~6 mm wide bridge closes them
   * hub: small arched-hole islands (largest 8 mm2, anchored)
   * fin: 0.66 mm2 of root ledge over the first 0.2 mm of the slot entry (was 25 mm2); blanking strip:
     the 2.1 mm wide ceiling of its magnet seat
   * blanking plug: its 12.6 mm2 slot-roof bridge
   * `tolerance_coupon_b`: its horizontal test holes (the point of the test)
10. **Mesh repair.** OCC's tessellation of some shells leaves small gaps or folded sliver pairs at B-rep
    sliver faces; the assembly's STL audit rewrites such files through `export_web.write_clean_stl` and
    keeps the module's file when the rewrite is not better. The B-reps are valid.
11. **Cutaway sector.** The carrier's theta 330 lock tab (r 27.2) lies in the removed sector beyond r 25:
    kept, it shows the twist lock.
12. **File sizes.** `blade.step` 0.75 MB; assembly STEPs 30.1 / 32.0 MB (stowed / deployed); GLB 9.3 MB, 277k
    triangles (Draco or meshopt compression would help the web page).
13. Bought items to verify before ordering: 3 x 22 and 3 x 4 pins, the Ruland collar's set screw, the
    8 x 12 x 0.5 shim (DIN 988 lists 8 x 14), N52 6 x 2 x 2 and 4 x 2 x 2 blocks magnetised through a
    2 mm face. The layout's two equator magnet pairs at the joint are not modelled (no room in the
    5.05 mm equator face); the joint relies on O-ring friction, the lug stop and the click.

## 9. Review fix round (2026-10-05)

The design review (visual, mechanical, print and conformance reviewers, triaged into T01 to T31)
was applied in one pass across the CAD modules. params.py, profile.py, layout*.py and bem.py were
not edited in that round: every value that differed from params was listed as a proposal, and the
params write-back has since applied them all (section 10).

| id | finding | what was done | measured after |
|---|---|---|---|
| T01 | spinner windows: stair-stepped outline, cantilevered treads | one outline per blade (convex hull of the slot bands + 1.5 mm, arc joins) with a 45 deg gable top; full depth only where the heel and knuckle sweep (it closes with 45 deg flanks and a 15 mm bridge at x 125.9, leaving a complete ring of the cone under the D22.6 seat ceiling); elsewhere a recess over a 1.2 mm skin whose floor rises at 45 deg out through the surface. The 0.4 mm outer chamfer could not be built (OCC chamfer fails on every B-spline window edge). A 2 mm hub trim was built and dropped (section 8) | no treads; window tops: 3 x 63.9 mm2 bridges at print z 17.3 anchored on both flanks; seat ceiling on a complete annulus; blade to spinner 0.664 mm (fold 0, was 0.635); spinner slides on at every fold, even deployed (min 0.675 mm); superseded by the tight radial windows of 9.2 (V1) |
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
| T13 | blade TE on tree-support interface | documented the support settings (top Z 0.2, dense interface, 5 mm knuckle brim); TE height corrected to 1.7 to 3.7 mm (1.7 to 2.2 along the aero span) | superseded by the print cradle (9.1) |
| T14 | blade B-rep classified far points IN | faces above 200 mm2 split + ShapeFix (orientation fixes and sewing alone did not help) | superseded (9.1, R1): the split caused the mesh streaks; the root cause was a parameterisation hook, now removed |
| T15 | stowed pair 1.045 mm apart; boss tip prints in mid-air | blade magnet 0.1 below its mouth, tail magnet 0.05 below the boss face | 0.701 mm face gap, about 1.73 N (was 1.27 N); superseded by the owner call (9.1): leaned 49 deg boss with a chin, 0.420 mm, 2.27 N |
| T16 | 38.6 mm nose-plug bridge | documented: good bridge settings, expect some sag on that internal face | the 45 deg cone from r 19.2 would reach the carrier (x 22.7) |
| T17 | no joint detent | documented "twist to the hard stop before stowing"; MAGNET count proposed (section 10) | superseded (9.1): three flexure tines click into dimples at the lock |
| T18 | cross-pin ligaments 1.0 / 0.55 mm | groove and tongue 0.8 deeper (4.8, 45 deg undercut kept), pin depth 2.2; 0.3 x 45 deg tongue bed chamfer; brim ears at the root ends | ligaments 1.2 / 1.2 computed (recheck R2: 0.9 / 1.1 measured); superseded (9.1): recessed rail with magnet retention, no cross pin |
| T19 | fin GProp unstable | faces above 100 mm2 split + ShapeFix | default 9616.7 / adaptive 9621.9 / mesh 9607.8 mm3 (0.05 %); superseded (9.1): skin of 7 G1 splines, no face split |
| T20 | claw A dog skin 0.69 mm | dogs to the flange OD (16.2), flange edge chamfer removed under them | skin 1.2; no overhang |
| T21 | stator insert holes 0.9 mm from the journal groove | no stator inserts: carrier screws self-tap into D2.5 pilots (45 deg lead-in) | 1.65 mm to the groove |
| T22 | claw B collar-pocket bridge ends on the bore | 0.2 mm sacrificial membrane over the bore in the print solid (drill it out) | full D16.2 bridge anchored all round |
| T23 | "LOCK" engraving too fine | text removed, the 0.9 mm chevron stays; the 0.9 mm counterbore web documented as non-structural | |
| T24 | nothing stops the shaft toward -X | documented: after stowing, pull the spinner aft | |
| T25 | collar set screw stub short of the flat | modelled down to the flat (r 3.5); v2 keying of claw B documented | contact, no overlap |
| T26 | coupon misses several fits | `tolerance_coupon_b`: D-bore 8.10 / 8.15 / 8.25 with the 0.5 flat, block seats 6.0 x 2.0 / 6.1 x 2.1 / 6.2 x 2.2 and disc seats D4.0 / 4.1 / 4.2 vertical and horizontal, an upright 3.7 mm ear with horizontal D2.9 / 3.0 / 3.1 holes | all probes pass |
| T27 | charcoal spinner, hub and fins | owner decision: all-ivory egg (9.1) | |
| T28 | cut faces share the skin colour | `front_shell_cutaway_section` / `aft_shell_cutaway_section` meshes (cap colour #C98B5E) lifted 0.03 mm | cutaway.png |
| T29 | exploded view 500 mm long, one blade hidden | two-level explode (shells apart, drivetrain opened between them, rotor aft), view 55 / 24 | replaced by the stacked axial explode (9.1, 2.24 x the egg) |
| T30 | closing the shells with stowed blades | step 12: unfold the blades before joining, stow after the twist | measured since by the closing sweep (9.1, R3) |
| T31 | deviations not written back | section 10 and the structured proposal list | applied by the params write-back (section 10) |

### 9.1 Polish round (2026-10-05): presentation and recheck items

| item | what was done | measured after |
|---|---|---|
| R1 blade mesh seams (major) | **Root cause**: the skin's U parameter (around the section) was the chord-length average over every row. The shank rows pulled it far enough that each U iso-curve left the trailing-edge corner backward and turned round within about 0.01 mm, a micro hook along the span. The hook made BRepClass3d call far points IN (the old T14), and made BRepMesh put about 38 600 nodes into the hook band: 7 400 folded sliver pairs, the streaks. **Fix**: centripetal U parameter averaged over the aero rows (`blade._u_params`). The T14 face split is removed. All meshes now go through `export_web.clean_mesh` with area-weighted normals (section 6). | no non-monotone row end of 122, worst end tangent 55.9 deg (was 180); classification on the unhealed solid: 0 / 1000 far points and 0 / 1500 near-surface points wrong (print, stowed, deployed, fold 45); blade.stl 12 042 triangles, 0 folded pairs, watertight; blade.step 0.75 MB (was 2.0); surface vs layout between stations 0.0028 mm, at stations 1e-13; GLB 0 folded pairs in all 67 meshes |
| R3 closing the shells | `assembly.closing_sweep`: front shell at the bayonet entry (both senses), pushed on from 40 mm, twisted to the lock in 2.5 deg steps, against the three blades at 13 folds, fins off and on. The result is a key-clearance row; `assembly_steps` restates step 12 with it. | fins off: 0 mm3 from fold 85 down (1 938 mm3 at 90); fins on: 0 mm3 from fold 50 down (2.0 at 60, 36.7 at 70, 107 at 85). Step 12: fins off and unfold by at least 10 deg, or by 45 deg with the fins on (one sweep step of margin on both, integration) |
| R5 presentation | packed and v2 viewer controls; fins lie flat on the ground in packed_kit; `fit=tight` framing for every render | sections 6 and 7 |
| R6 stale outputs | `build.py --fresh` prunes every output the current parts do not produce | run 2026-10-05 against the polish build's part cache: removed `hydroegg_assembly_{stowed,deployed}.step`, `inserts_m3_stator.step`, `hydroegg.glb`, `test.glb`, `manifest_test.json`; every current part's files kept |
| T13 blade on supports | `blade_cradle.stl`: a sacrificial support cradle in the blade's print coordinates. Rails under the heel and under the shank ramp, trailing edge and tip, top 0.2 mm below the blade (the support top-Z distance), 1.6 mm wide or the width of the blade's lowest 0.5 mm band + 0.4. Print it with the blade, slicer supports off. | 2 rails, 81 stations, up to 2.31 mm tall, 0.255 cm3; 0.200 mm to the blade, 0 mm3 overlap; it carries 57 % of the 126 mm2 of blade faces beyond 45 deg. The rest are the bridged ceilings of the two magnet seats and faces within 1.5 mm of the bed at the knuckle (section 8). Refined in 9.2 (P6): stations every 0.25 mm, 71 % |
| colours (owner) | all-ivory egg: `hub` and `fin` materials = `shell`; printed drivetrain internals moved to a new charcoal `drive` key | renders |
| cable gland (owner) | `gland_plug` (drivetrain): low domed M12 blanking plug, push-in. D10.0 spigot with 3 half-round crush ribs to D10.6 in the D10.5 tap-drill bore; it also grips a tapped thread on its crests. 0.4 mm band on the D18 seat flat (r 8.5), a shoulder arc from vertical to 45 deg, a D13.5 flat face with a 10 x 1.2 x 1.0 coin slot. Printed flat face down, spigot up. `cable_gland_m12` gets `v2: True`, the plug `v2: False`. | valid, 1.40 cm3; print overhang 12.6 mm2, all of it the slot-roof bridge; rib press overlap 0.192 mm3, seat face 0.000 mm |
| CLEAR_ROTATE (engineering call) | radial per side for rotating fits (the knuckle on its pin, 0.15). Running and swept clearances read `RUN_CLEAR` (0.30, `getattr(P, "RUN_CLEAR", 0.30)` until params carries it) in `assembly.py` and in the drivetrain's bell journal. | hub and spinner fold sweep against `RUN_CLEAR`; geometry unchanged |
| fins: recessed rail (owner decision; T07, T08, R2) | the proud pads are gone. The fin study (`layout_internal.py`) was re-run with the root on the skin: 30 x 70 Clark Y at 12.5 deg, theta 87.5 / 267.5, yaw pivot 0.05 c, x 26.2 to 55.5, 215 mm tip to tip (was 219). The rail is a straight half dovetail along the yawed chord (5.2 neck, 45 deg undercut to 6.5 at the floor) sunk into an inner boss of the front shell; it climbs inward going aft (slope -0.25), so the fin slides aft and in onto the end stop at x 39.4 and drag pushes it home. Entry relief where the undercut lip would be under 1.2 mm. Retention: an N52 4 x 2 x 2 block in the tongue over one in the groove floor (a cross pin cannot cross a flush groove). Flush blanking strips for the packed egg. Fin skin of 7 G1 splines (replaces the T19 face split) | shell proud of the hull at the fins 0.0 mm; couple 1.103 x the GL40 load-line torque (layout), 1.142 x as built (band 1.10 to 1.155); groove walls (OCC rays) floor 1.552, seat bottom 1.297, flat wall 1.445, neck 5.342, undercut lip 2.355, end stop 2.634; slot to pocket rims 5.348 (arc), grooves to pocket voids 4.444; bay to shells 0.637; fin root 0.101 to 0.102 over the skin; tongue 0.25 on the walls, 0.05 over the floor; fin and strip pulled out in 0.5 mm steps over 20.1 mm: 0 mm3 overlap; magnets 2.36 N seated, about 1.3 N to slide out (fin 0.04 N); fin prints tip down with no supports; GProp default 5588.2 / adaptive 5587.8 / mesh 5581.2 mm3, 0 folded pairs |
| T15 stowed detent (owner call) | the aft-shell boss axis follows the blade face normal (lean 49.03 deg), its top follows the blade's inner face 0.1 mm off it (cubic fit), and a 2.5 mm wide chin swept down the face plane into the root-notch ramp starts the boss (kept 0.25 mm off the blade's lower envelope); the blade seat stays where it is | pair face-centre gap 0.420 mm, axis angle 1.66 deg, about 2.27 N; boss to blade 0.088 mm, 0 mm3; targets 0.45 mm / 2 deg / 2.0 N met; boss and chin tops unified in 9.2 (C2 / P2): boss to blade 0.171 mm |
| T17 bayonet click | three curved flexure tines in the spigot wall between the lugs, each with a 0.5 mm bump that snaps 0.25 mm under the socket ring into dimples at the lock angles; tine slits vertical in print, tine undersides 45 deg slots | 0 mm3 locked, 1.82 mm3 along the twist; ESTIMATE 2.97 N per tine, 17.7 MPa; slots widened and print-only ties added in 9.2 (C1 / P1) |
| coupon | the fin sample is the recessed groove with the floor magnet seat (45 deg top wall) and a matching tongue key (0.05 floor, 0.25 walls) | side play verified by +-0.23 / 0.27 mm shifts; all coupon features pass |
| shells check fixes | pocket-floor wall rays skip start points inside the chin; the proud-of-hull check ignores the tow eye; the bayonet channel cutters' front face below R_LIP leans 45 deg (no flat ledge facing the bed where a rail boss reaches past x 41.7) | shells checks all pass |
| params write-back | every deviation above written into params.py (section 10); layout studies re-run | `data/layout.json` 30/30; `fin_block_vs_params` empty |
| integration (fresh) | `build.py --fresh` from params; renders: side and front frame their ground shadow (`shadow=1`), the exploded view drops the ground shadow (`ground=0`) | section 5 and 7 numbers |

### 9.2 Review 2 round (2026-10-05): changed-code, visual and print reviewers

Every finding was first reproduced on the built outputs (the reviewers' scripts in
`out/work/cad/review2_*` and probes in `out/work/cad/fix3`); the confirmed ones were fixed in one pass,
then the whole model was rebuilt (`build.py --fresh`). Refuted or deferred items are listed in section 8.

| id | finding (reproduced) | what was done | measured after |
|---|---|---|---|
| C1 / P1 click tines (major) | each tine's 45 deg underside started as a knife edge on the bore, cantilevered ~14 mm over a 0.42 mm (normal) slot: the first lines would sag and fuse the tine | the underside slot is 1.15 along x (0.81 normal, `BAYONET.click gap_w`); the PRINT solid (`aft_shell_print`, `aft_shell_cutaway_print`) carries a 0.8 mm breakaway tie across each free-end slit at the tine's bottom (`tie_h`), so the underside's first lines bridge between the tine root and the wall beyond the slit; cut the 3 ties through the slits after printing. The world solid has no tie. The axial-finger alternative was rejected: the slits must stop before the O-ring seal zone (x 49.6), so a finger is at most ~6 mm long and the 0.25 mm snap would strain it ~1.3 % across the layer lines (33 to 44 MPa in PLA) | slot 0.800 to 0.818 mm normal (`tine_gap.py`), ties 4.41 mm3 in the print solid; the tie section joins the tine's first lines to the wall beyond the slit (`out/work/cad/fix3/tie_sec.png`); twist overlap 1.818 mm3 and 0 mm3 locked as before |
| C2 / P2 detent boss and chin | the contoured boss top and the chin top leaned up to 57.6 / 58.3 deg; the boss rim hung 0.25 mm above a lower chin top, so a 0.64 mm2 island started at print z 1.9 | boss and chin share one top surface (`shells.stow_top_grid`): the stowed blade's inner face from the layout surface on a 0.25 mm grid, face_clear below it, never falling going print-down (a running minimum), built as bilinear patches; the magnet keeps its blade-driven depth | faces at most 49.01 deg (3 x 32.3 mm2 at the boss plane), 0 islands at the 0.05 / 0.1 / 0.15 phases (`chin.py`: no free component from z 1.6 to 2.2); boss to the CAD blade 0.171 mm (was 0.088); the static magnet face stands up to 0.20 mm proud of the top at one side of the hole |
| C3 stowed-pair gate | the assembly gate checked only face gap <= 0.5 | gated on the owner's three T15 targets (gap <= 0.45 mm, axis <= 2.0 deg, force >= 2.0 N); the boss and tail-magnet notes give the 49 deg lean | 0.420 / 0.365 mm, 1.66 deg, about 2.27 N: PASS on all three targets |
| C4 / P3 fin ledges | the 4.7 mm neck stood 0.4 to 1.75 mm beyond the 3.35 mm airfoil root: flat ledges facing the bed in the tip-down print; the tip-corner rounds overhung up to 74 deg near the bed | a root flare from each neck edge up to the airfoil skin at 43 deg (in section; a root fillet up to ~2 mm tall at the LE end, `fins._root_flare`); the tip-corner rounds stop where their face leans 42 deg and end on a short land perpendicular to the tip line (`fins._tip_corner`) | fin_0.stl: 0.66 mm2 beyond 45.5 deg (the first 0.2 mm of the slot entry), was 25 mm2 up to 88 deg; seated in the shell 0.013 mm3 (end-stop facets), slide-out path clean |
| C6 strip grip | the flush strip had nothing to pull it forward with | a fingernail groove across the strip's outer face 2.6 mm ahead of its aft end: 1.2 wide, 0.6 deep, steep forward wall (the nail pulls on it), 45 deg aft flank | strip seated and pulled out: 0 mm3; flush (max proud 0.0) |
| C7 rail boss of fin 1 | the boss was turned 180 deg for fin 1, so its 2.5 mm extension pointed at open wall and a 0.55 to 1.6 mm slit stayed against the theta-240 liner | each boss widens toward its nearer pocket (`shells.boss_extra_side`: +w for theta 87.5, -w for 267.5) | `crevice2.py`: no air run between the theta-267.5 boss and the 240 liner at x 32 to 40 (was 0.55 to 1.6 mm); groove walls >= 1.297 on both fins |
| C8 front shell note | it said to tap the M12 seat | display build: leave the 10.5 mm bore for the push-in plug; v2: tap M12x1.5 | |
| C9 params pad_x | labelled as the slot mouth | kept (it feeds the fin study's interval checks), relabelled as the study's footprint envelope; the built mouth spans x 27.80 to 39.95 | |
| V1 rotor band (major) | stepped, acute spinner windows (up to 69 deg wide, flanks ~50 deg to the skin), hub steps and the recessed shank behind them; the tail read torn | spinner windows rebuilt (`rotor.window_outline0`, `spinner_window0_radial`): the swept root's footprint in the spinner material + RUN_CLEAR + 0.7, closed and smoothed, cut RADIALLY (flanks normal to the skin), open at the front edge for the slide-on over the stowed blades, rounded top. The hub raise and the shank fairing were not built (section 8) | window area in the spinner -45 % (theta -38 to +26 deg, aft end x 127.0, top bridge 12.2 mm, cone ring 1.0 mm); blade roots to the spinner 0.98 mm over the fold (was 0.69); open skin (`openarea.py`, packed display) x 108-118 46.1 % (was 57.1), x 118-128 36.7 % (was 50.5) |
| V2 nose plug (major) | the plug sat on a 5.1 mm plain cylinder boss: a two-step spigot 9 mm off the nose with a groove between seat chamfer and plug | the seat boss blends into the nose through an R6 root fillet (`shells._gland_root_fillet`; R7 failed in OCC), the seat lost its rim chamfer and the plug's base is flush with the seat edge (`GLAND_PLUG.base_r` 9.0): one smooth bump. Sinking the seat (moves the v2 interface) is a proposal for the owner | one bump: boss 5.3 mm, plug crown 9.1 mm above the skin (`plug_height.py`); rib press 0.192 mm3, seat face 0.000 mm |
| V3 exploded view (major) | diagonal composition, empty lower half, a blade across the aft shell rim | view 20 / 18, rotor turned 90 deg: one horizontal line in build order with the blades radially up and down; blades 60 mm out (was 34), plug 28 mm out of its seat (was 16); studio backdrop | 406.5 mm long, 2.28 x the egg |
| V4 pocket tips | the pocket's forward end was square (layout surface) while the CAD tip rounds off: an ~8 mm open triangle ahead of each stripe | over the tip zone the pocket takes the CAD blade's theta range + 2 mm, narrowing at the 45 deg rule: a chamfered corner that follows the rounded tip; the liner keeps the full bounds; every pocket check uses the CAD tip (`shells.blade_cloud`) | front row -theta bound -8.4 deg (was -18.0), 19 tip rows; the open triangle ahead of the tip (up to ~8 mm along x) is now a thin crescent under 4 mm; all pocket checks pass on the CAD tip |
| V5 hero framing | the near fin dominated the hero and cutaway foreground | cameras at azimuth -18 (`build.py VIEWS`) | renders |
| V6 shadow stipple | key-light self-shadowing on the blades | `key.shadow.normalBias = 0.5` | renders |
| V7 backdrop | ivory on ivory, contrast 1.10 | a `theme=studio` warm-grey backdrop (#d9d4cb to #b8b2a8) for the render set; light stays the viewer default | renders |
| P4 fin on plate 3 | 73.7 mm tower on a 3.6 mm first layer, packed with its thin axis along Y | the print plan packs the fins with their chord along Y (an A1's bed motion axis); notes: 5 mm brim all round, outer walls <= 50 mm/s, 6 to 8 s minimum layer time | plate 3: fin_0 (5, 5, turned 90), fin_1 (29.1, 5, turned 90) |
| P5 strip tip | `strip_tip_min` was declared but unused: a 0 to 0.47 mm knife edge | the strip's entry wedge is cut where thinner than 0.6 mm, with a 0.3 x 45 deg break; print the strips at 0.08 to 0.12 mm layers | strip front end 0.62 mm thick at its first sample (was 0.00 to 0.47) |
| P6 blade cradle | the rail top sat 0.47 to 0.6 mm under the shank ramp, nothing under the face next to the knuckle | stations every 0.25 mm (`BLADE_CRADLE station_step`, was 1.0) | 328 stations in 2 rails; top 0.20 to 0.33 mm under the blade along the span (`te_gap.py`, was up to 0.57); 71 % of the faces beyond 45 deg carried (was 57 %) |
| P7 coupons | the fin key printed like the strip, no 49 deg disc seat | a standing tongue key printed like the fin (`tolerance_coupon_key_tongue_up`, on a 10 mm stub, tongue up) and D4.0 / 4.1 / 4.2 seats on a sloped +X end of coupon B with their axes 49 deg from vertical; how_to_use rewritten | key 0.70 cm3, 10 x 5.8 x 15.2 mm; coupon B 18.40 cm3, all coupon features pass |
| P8 spinner seat ceiling | the D22.6 clamp ceiling is pierced by 4 holes (short cantilevers) | documented: face the clamp ring flat on sandpaper over glass before assembly (step 10, print notes) | |
| P9 thin webs | stator replica holes 0.24 mm from the chamfered rim; carrier web 0.56 mm in its top layer | stator: D1.6 x 0.6 dimples on PCD 39 instead of the D2.5 holes on PCD 40; carrier: D6.0 counterbores and no centre-bore chamfer on that face | stator dimples 1.19 mm from the chamfered rim in the first layer (was 0.24); carrier web 1.01 mm in every layer (was 0.56 at the top) |
| P10 plate 1 colours | the blue replica shared the ivory shell plate | the generator replica has its own plate (5, blue PLA) | plate 5: generator_bell (5, 5), generator_stator (59.5, 5) |

## 10. Params write-back (applied 2026-10-05; review 2 rows at the end of the table)

Every deviation that sections 9 and 9.1 and the fins round listed as a proposal is now in
`ovo/params.py`, so params is again the single source of truth for the CAD. The CAD modules read
the new keys instead of literals (geometry identical, except where noted). Both layout studies were
re-run on the new params (`python -m ovo.layout`, then `python -m ovo.layout_internal`):
`data/layout.json` passes 30/30 interface checks and both cross-study checks, no check limit was
loosened, and the claw B running check now reads `RUN_CLEAR` (still 0.30).

| params key | now | was | why |
|---|---|---|---|
| `CLEAR_ROTATE` | 0.15, RADIAL per side for rotating fits (knuckle bore 3.3 on the 3.0 pin) | 0.30, ambiguous | engineering call |
| `RUN_CLEAR` (new) | 0.30: swept / running clearance between moving bodies | (CLEAR_ROTATE) | engineering call; rotor.py `RUN`, drivetrain.py, assembly.py and layout.py read it |
| `ROTOR_HUB.ear_hole_d` | 2.9 = PIN_D + 2 PRESS_FIT (printed press bore) | 3.0 | T31 |
| `ROTOR_HUB.spinner_wall` | 1.25 | 1.2 | T01: the 1.2 shell measured 1.167 on the mesh; now 1.217 min (p01 1.236) |
| `ROTOR_HUB.pin_channel_tau` | [11, 16.876] | [11, 16.941] | layout re-run with the 1.25 wall |
| `ROTOR_HUB.pin_filler`, `pin_filler_l` (new) | ISO 8734 3 x 4, 4.0 | | pin retention |
| `DETENTS` deployed pair | `deployed_magnet`: N52 block 6 x 2 x 2, seat 6.1 x 2.1 x 2.2, r_e 8.2, recess 0.2, gap 0.4, qty 6; hub centre [112.3, 8.2, 0.0], blade [-5.8, 0, -1.3] | D4 x 2 disc at r_e 8.4388 / [-5.5612, 0, -1.3] | T03 |
| `DETENTS` stowed pair | leaned boss (owner call T15): seated blade centre [10.899, 7.1625, 7.0023], blade axis [0.7734, -0.2306, -0.5905], boss on the face normal (lean 49.03 deg, max 49.5), boss r 3.25, face clear 0.1, recess 0.05, contoured top, chin {6.5, 1.25, 0.25}, rest zone r 8.25, static centre [103.259, -6.5972, 19.539] and axis [0.7549, -0.2349, 0.6124] at theta 0, pair face gap 0.413 (0.420 measured) | radial pair, boss top r 21.196 | T15 owner call; `stowed_blade_magnet_blade_frame` [11.75, 7.1989, 7.2073] stays (blade.py anchors the seat on it, the seat stays where it is); `stowed_static_magnet_r` and `stowed_boss_top_r` kept, marked superseded |
| `MAGNET` | D4 x 2 N52, qty 6 (stowed pair only) | 16 per build | the deployed and fin magnets are their own blocks: 18 magnets per build |
| `FIN_MAGNET` (new) | N52 4 x 2 x 2, seat 4.1 x 2.1 x 2.05, recess 0.05, 1.75 ahead of the stop, qty 6 | | recessed rail (owner decision) |
| `FIN` | 30 x 70 Clark Y at 12.5 deg, theta 87.5 / 267.5, yaw pivot 0.05 c (x 27.6756), x_le 26.2111, root on the skin (r 32.05 / 37.6), recessed half-dovetail rail (neck 5.2, floor 6.5, slope -0.25, stop x 39.4), boss walls 1.6 / 3.45, magnet retention | 40 x 65 at 12 deg, theta 93 / 273, proud 10 mm pads, M3 grub | owner decision; fin study re-run with the root on the skin |
| `SHAFT.hub_flat_x` | [108.4, 128.0]; new `collar_flat_x` [55, 62], `collar_flat_depth` 0.5 | 111.6 | T31: the D-flat runs out through the aft end so the hub slides on |
| `COUPLER` | claw A x [46.5, 52.3], L 5.8, OD 32.4, flange 2.0, dogs to the OD; claw B x [49.5, 62.5], L 13.0, floor 1.2, cone x [60.5, 61.5] | claw A x [46.5, 51.5], claw B x [52.5, 63] | as built (T31, T20) |
| `CARRIER` | `pad_window` (the hardware.json cad_window merged with the notch, open to r 12.0); engraving chevron only | | T31, T23 |
| `INSERT_M3` | stator 'none', `stator_pilot_d` 2.5 | | T21 |
| `JOINT` | `seal_lead` [0.85, 18.0], `seal_lead_mouth_d` 65.79 | 0.4 x 20 | T05 |
| `BAYONET.click` (new) | 3 tines, lock 120 / 240 / 0, arc 28, bump 0.5, dimples | | T17 |
| `TOW_EYE` | teardrop outer profile, apex r 1.5, hole apex flat 2.0, root fillet r 2.0 | | T06 |
| `GLAND` | `tap_drill_d` 10.5; `display`: blanking plug | | owner decision |
| `GLAND_PLUG` (new) | base_r 8.5, band 0.4, dome_r 6, spigot D10 x 6, chamfer 0.5, 3 ribs r 0.5 to D10.6, slot 10 x 1.2 x 1 | | owner decision |
| `SPINNER_ACCESS` | x [131.0, 148.354]; `countersink_x` [148.354, 151.725] | x [134, 165] | T02 |
| `SPINNER_CAP` (new) | x_joint 151.725, spigot 5.7 x 12, rib tips D6.1, rib_w 0.6, groove chamfer 0.8, x_tip_built 160.99 | | T02 |
| `BLADE_CRADLE` (new) | gap 0.2, min_w 1.6, station step 1.0, min_float 0.4 | | T13 |
| `POCKETS.cad_floor_print_deg` (new) | 44.5 (floor_r is a minimum on the root-notch ramp) | | the ruled loft twists exact 45 deg ramps past 45 |
| `ROTOR_R_ENVELOPE` (new) | 97.31, measured on the deployed blade corners | | `ROTOR_R_TIP_TARGET` stays 97.0: bem.py and layout.BladeDef read it, 97.31 would grow the blade |
| `HULL_L` | stays 165 (built tip 160.99, L/D 2.15) | | the Myring point; the cap ends at its R3 round |
| `DETENTS.stowed_x` | stays 104.25 | | the owner call keeps the seat |
| `BAYONET.click.gap_w`, `tie_h` (review 2) | 1.15 (0.81 normal), tie 0.8 | 0.6 | C1 / P1: wider underside slot, print-only breakaway ties |
| `GLAND.boss_root_fillet_r` (review 2) | 6.0 (R7 failed in OCC) | | V2: the seat boss blends into the nose |
| `GLAND_PLUG.base_r` (review 2) | 9.0 = seat_d / 2 (no seat chamfer) | 8.5 | V2: plug flush with the seat edge |
| `BLADE_CRADLE.station_step` (review 2) | 0.25 | 1.0 | P6 |
| `FIN.pad_x` (review 2) | value kept, comment: the study's footprint envelope (built mouth x 27.80 to 39.95) | | C9 (the fin study's interval checks read it) |

Module edits that came with it: `common.py` (DEP_MAGNET, STOW_* from DETENTS), `shells.py` (FIN =
params.FIN; FIN_RAIL, STOW_DETENT, BAYONET_CLICK, SEAL_LEAD, FLOOR_PRINT_DEG, EYE_*, GLAND_DRILL_D
from params; `fin_block_vs_params` now reports drift and is empty), `rotor.py` (RUN = RUN_CLEAR,
EAR_HOLE_D, SPINNER_CAP, FILLER_L, WIN_SKIN_N = 1.25), `drivetrain.py` (claw values, hub flat, collar
flat, stator pilot, PLUG, GLAND_DRILL_D, RUN_CLEAR), `blade.py` (CRADLE_*), `assembly.py`
(RUN_CLEAR, BOSS_ZONE_R = `stowed_rest_zone_r` 8.25), `layout.py` (spinner_wall, run_clear,
knuckle_bore_d = pin_d + 2 CLEAR_ROTATE = 3.3, ear_hole_d = pin_d + 2 PRESS_FIT).

Still different on purpose: the studies' own `proposed_param_updates` in the layout JSONs are
pre-CAD study outputs (DETENTS deployed disc at r_e 8.4388, the COUPLER claw x ranges, SHAFT
hub_flat_x 111.6, SPINNER_ACCESS x [134, 165]), and `data/layout.json` `part_split` still carries the
study's print orientations (carrier front face down, fins flat). params.py and `print-plan.md` are
the as-built truth (params.py header).
