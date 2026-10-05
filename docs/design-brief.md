# Ovo: design brief

Product name: **Ovo** ("Ovo, a folding hydrokinetic charger"). Portfolio project: a generic portable
hydrokinetic charger, modeled as parametric code-CAD (build123d), 3D-printable on a
Bambu Lab printer (256 mm cube).

## Decisions taken with the owner (grill, 2026-10-03)

| Topic | Decision |
|---|---|
| Scope | Display model with real mechanics, no electronics. Rotor spins on real bearings, blades really fold. |
| Architecture | Axial (horizontal-axis) turbine. **No duct**: it does not fold and would double the stowed diameter. |
| Form | Elongated egg, L/D ~ 2.2 (Myring profile): about 75 mm diameter x 165 mm. Stowed it is a clean egg; deployed the rotor is ~194-200 mm. |
| Blades | 3, folding. |
| Generator bay | "v2-ready": real interfaces for a low-KV gimbal outrunner, shaft seal pocket, equator O-ring groove, cable gland. The display build carries a printed replica with the same envelope. |
| Extras | Real bearings. Cutaway variant of the shells that shows shaft, bearings and generator. |
| CAD | build123d (Python), everything parametric from `ovo/params.py`. |
| Printer | Bambu Lab A1 / P1S / X1, 256 mm cube. Single colour (no AMS assumed). |
| Deliverables | Public GitHub repo (README, renders, exploded GIF, STL/STEP from CI), technical report (PDF), project web page with an interactive 3D viewer. BOM and assembly guide always. |
| Language | English for repo, web and report. |

## Decisions after the first CAD review (2026-10-05)

| Topic | Decision |
|---|---|
| Name | Ovo. Python package `ovo`; the local folder keeps its old name. |
| Fins | Recessed dovetail rail cut into a thickened inner wall of the front shell, closed by a flush blanking strip when the fins are off. Stowed, the egg stays at its 75 mm diameter. |
| Colours | All-ivory egg (shells, hub, spinner, fins) with orange blades: stowed it reads as one egg with three orange stripes. |
| Cable gland | Display build carries a low M12 blanking plug; the gland shows in the cutaway and the v2 documentation. |

## Engineering decisions (taken in design, open to change)

1. **Downstream rotor, blades fold forward.** The hinge axis is tangential on the
   hub. A folded blade points upstream along the hull, so the flow drags its tip
   aft and opens it; thrust then holds it against the deployed stop. With an
   upstream rotor the flow would hold the blades shut.
2. **Pre-cone derived from the hull.** Folding a blade by ~90 deg about a tangential
   axis maps "aft" to "outward". So a blade built with an aft pre-cone that copies
   the hull curve will hug the convex tail when stowed. The deployed rotor gets mild
   downstream coning for free (passive yaw stability, lower root bending).
3. **Shallow pockets in the static tail** take the stowed blades near flush and
   index the rotor so it cannot spin in the bag.
4. **Egg splits at the equator** (Kinder-egg style) with a bayonet and an O-ring
   groove (v2 sealing). The cutaway variant swaps in shells with a quarter removed.
5. **Drivetrain**: 8 mm stainless shaft on two bearings in the static tail; the hub is
   keyed to the shaft; the spinner retains the hinge pins. A lip-seal pocket aft of
   the rear bearing is there for v2.
6. **Roll stability**: generator reaction torque wants to spin the hull. Two small
   slide-in side fins with opposite incidence provide counter-torque (sized in the
   layout study). Tow/anchor eye at the nose; the cable exits through a nose gland
   (v2) next to it.
7. **Detents**: small N52 magnets hold blades both stowed and deployed, so the
   display model stays put and the mechanism needs no springs.
8. **Blade design is constraint-aware BEM**: Schmitz/BEM optimum, then linearised
   chord and twist clipped to printability (min chord, min trailing edge) and to
   the folding constraint (root twist). The report quantifies the Cp cost of each
   constraint. Polars from NeuralFoil (an XFOIL-trained model).

## Out of scope for v1

Electronics, a real waterproof test, CFD. The design keeps the interfaces so v2 can add
them without redesign.

## Repo map

- `ovo/params.py`: every dimension.
- `ovo/profile.py`: hull profile r(x).
- `ovo/*.py`: one module per subassembly; `build.py` exports everything.
- `data/`: design-study outputs (JSON) consumed by the CAD.
- `docs/`: brief, design notes, report sources, figures.
- `out/`: generated STL / STEP / GLB / renders (git-ignored; CI publishes them).
