// Ovo technical report. Source of docs/report/ovo-report.pdf.
// Build with: .venv/bin/python scripts/report.py   (regenerates data.typ, refs.bib, figures, images)
// Every number below comes from data.typ (generated from the study JSONs, the CAD checks and params.py).
#import "data.typ": d, tb

// ---------------------------------------------------------------- style
#let ink = rgb("#1d1d1b")
#let ink2 = rgb("#5a5955")
#let rule-c = rgb("#cfccc4")
#let tint = rgb("#f3f1ec")
#let orange = rgb("#e8622c")
#let blue = rgb("#2a78d6")

#set document(title: "Ovo: a folding hydrokinetic charger. Technical report", author: d.author)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.2cm, bottom: 2.2cm),
  header: context {
    if counter(page).get().first() > 1 {
      set text(8pt, fill: ink2)
      [Ovo: a folding hydrokinetic charger #h(1fr) Technical report, v1 display build]
      v(-0.55em)
      line(length: 100%, stroke: 0.4pt + rule-c)
    }
  },
  footer: context {
    if counter(page).get().first() > 1 {
      set text(8pt, fill: ink2)
      h(1fr)
      counter(page).display()
      h(1fr)
    }
  },
)
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en", fill: ink)
#set par(justify: true, leading: 0.6em, spacing: 0.85em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => {
  set text(14pt, weight: "bold")
  v(0.5em)
  it
  v(0.2em)
}
#show heading.where(level: 2): it => {
  set text(11.5pt, weight: "bold")
  v(0.25em)
  it
}
#show heading.where(level: 3): it => {
  set text(10.5pt, weight: "bold", style: "italic")
  it
}
#set math.equation(numbering: "(1)")
#show figure.caption: it => {
  set text(9pt)
  set par(justify: true)
  it
}
#show figure.where(kind: table): set figure.caption(position: top)
#set figure(placement: auto)
#show figure.caption: set align(left)
#show raw: set text(font: "DejaVu Sans Mono", size: 0.82em)
#set list(indent: 0.6em, spacing: 0.6em)
#set enum(indent: 0.6em, spacing: 0.6em)

// "est": the number before it is a model estimate (BEM, generator model, layout mass model), not a measurement
#let est = box(text(size: 0.62em, fill: orange.darken(15%), baseline: -0.5em, [ est]))

// booktabs-style table; rows whose last element is `true`/`false` treat it as a highlight flag
#let btab(cols, header, rows, align: left, size: 8.3pt, flag: false, markup: false) = {
  set text(size: size)
  set par(justify: false)
  let cell(c) = if markup and type(c) == str { eval(c, mode: "markup") } else { c }
  let body = rows.map(r => {
    if flag {
      let hl = r.last()
      r.slice(0, r.len() - 1).map(c => if hl { strong(cell(c)) } else { cell(c) })
    } else { r.map(cell) }
  })
  table(
    columns: cols,
    align: align,
    inset: (x: 3.5pt, y: 2.6pt),
    stroke: none,
    table.hline(stroke: 0.8pt + ink),
    table.header(..header.map(h => text(weight: "bold", h))),
    table.hline(stroke: 0.45pt + ink),
    ..body.flatten(),
    table.hline(stroke: 0.8pt + ink),
  )
}

#let note(body) = block(fill: tint, inset: (x: 8pt, y: 6pt), radius: 2pt, width: 100%, text(size: 9.2pt, body))

// ---------------------------------------------------------------- title page
#page(header: none, footer: none)[
  #v(-0.4cm)
  #text(9pt, fill: ink2, tracking: 0.08em)[TECHNICAL REPORT #h(0.5em) · #h(0.5em) V1 DISPLAY BUILD]
  #v(0.25cm)
  #text(30pt, weight: "bold")[Ovo]
  #v(-0.55cm)
  #text(15pt)[A folding hydrokinetic charger, designed end to end as parametric code-CAD]
  #v(0.25cm)
  #text(10pt)[#d.author #h(0.4em) · #h(0.4em) #d.affiliation]
  #v(0.25cm)
  // the render's empty paper above and below the rotor is trimmed so the abstract and both tile rows fit page 1
  #box(width: 100%, height: 9.5cm, clip: true, radius: 0pt)[#move(dy: -0.35cm, image("img/hero_deployed.jpg", width: 100%))]
  #v(-0.3cm)
  #text(8.5pt, fill: ink2)[Rendered from the CAD assembly: blades deployed on their stops, roll fins fitted. Rotor diameter #d.rotor_d~mm, egg #d.hull_d~mm in diameter.]
  #v(0.3cm)
  #block(width: 100%)[
    #set par(justify: true)
    #text(weight: "bold")[Abstract.]
    Ovo is a portable hydrokinetic charger for rivers and towing: a #d.hull_d by #d.x_tip_r~mm egg that
    unfolds into a #d.rotor_d~mm three-blade rotor. The rotor sits downstream of the hull and its blades fold
    forward into shallow pockets in the tail, so that by design the flow opens them and the rotor thrust holds
    them on their stops (a water test must show that the flow overcomes the stowed detent); stowed, the
    product is a closed egg. The
    blades come from a constraint-aware blade element momentum (BEM) design with Prandtl losses, the
    Buhl high-induction correction and NeuralFoil polars: every printability and folding constraint is
    priced in power coefficient, and together they cost #d.chain_total_abs~% against the unconstrained
    optimum. The
    model predicts a peak #box[$C_P$ = #d.cp_peak15]#est and #d.p_shaft15~W#est of shaft power at
    #d.design_flow~m/s. A generator trade study selects the CubeMars GL40 for a v2 build and shows that
    winding losses move the electrical optimum far above the rotor's #box[$C_P$] peak, to a load line that
    must also respect a cavitation limit. The whole machine (#d.n_parts_total parts) is built by one
    Python pipeline in build123d from a single parameter file and is verified automatically:
    interference in #d.n_states assembly states, a rotor spin sweep, fold and closing sweeps and
    #d.kc_n key clearances on fine meshes of the CAD solids, and watertight print meshes. The v1 is a display build with real
    bearings, folding blades and magnetic detents but no electronics. It has not been printed or tested
    in water yet; every performance figure in this report is a model estimate and is marked #est.
  ]
  #v(0.15cm)
  #let tile(v, l) = box(width: 100%, height: 1.42cm, inset: (x: 5pt, y: 5pt), fill: tint, radius: 2pt)[
    #set par(justify: false)
    #text(13pt, weight: "bold", v) \ #text(7.8pt, fill: ink2, l)
  ]
  #grid(
    columns: (1fr, 1fr, 1fr, 1fr), gutter: 5pt,
    tile[#d.hull_d × #d.x_tip_r~mm][stowed egg, blades flush in pockets],
    tile[#d.rotor_d~mm][rotor diameter, 3 folding blades],
    tile[#d.p_shaft15~W#est][shaft power, #d.design_flow~m/s, $C_P$ #d.cp_peak15],
    tile[#d.gl_td15~W#est][DC through a diode bridge, GL40 generator, #d.design_flow~m/s],
    tile[#d.n_parts_total][parts in the assembly],
    tile[#d.kc_ok / #d.kc_n][key clearances pass, measured on the CAD],
    tile[#d.print_g~g][printed filament, display build],
    tile[#d.n_stl][STL print files, all watertight, shells support-free],
  )
]

// ---------------------------------------------------------------- contents
#{
  set text(9.5pt)
  show outline.entry.where(level: 1): it => {
    v(0.25em)
    strong(it)
  }
  outline(title: [Contents], depth: 2, indent: 1.2em)
}

#note[
  *How to read the numbers.* Ovo exists as CAD, design studies and print files. Nothing has been
  printed, bench-tested or put in water. Values marked #est are outputs of models (BEM with
  NeuralFoil polars, a generator loss model, a layout mass and buoyancy model) and carry their
  uncertainty, discussed where they appear. Dimensions and clearances are measured on the CAD solids by
  the build's own checks. The design and analysis numbers in this report are pulled by `scripts/report.py`
  from the study and check JSON files and the parameter file of the repository, so the report cannot drift
  from the model.
]

= Introduction

== Problem

Moving water is a dense energy source. The kinetic power through a rotor disc of area $A$ is

$ P_"avail" = 1/2 rho A U^3, $ <eq-pavail>

so a #d.rotor_d~mm rotor in a #d.design_flow~m/s river current (water at $rho$ = 1000 kg/m³) sees
#d.p_avail15~W. A good small rotor extracts a third of that, which is enough to charge a phone or a power
bank from a canoe camp, a moored boat or a yacht under tow. The difficulty is size: a rotor large enough
to matter does not fit in a bag, and the commercial products reviewed in @sec-prior either use small
rotors or rigid housings. Ovo's goal is the largest practical rotor when running and a closed, compact
shape when carried.

== Requirements and scope

@tab-req lists the requirements set with the project owner at the start (`docs/design-brief.md`) and
where each one is met. The design flow speed is #d.design_flow~m/s, a decent river current; towing
covers 2 to 2.5 m/s. The working range studied is 1.0 to 2.5 m/s.

#figure(
  btab((auto, 1fr, 1.15fr),
    ([ID], [Requirement], [How the v1 design meets it]),
    (
      ([R1], [Stows as a closed shape that fits a bag; nothing protrudes], [#d.hull_d × #d.x_tip_r~mm egg; blades in tail pockets, at most #d.stow_flush~mm above the hull; fins slide off, flush strips close their rails]),
      ([R2], [Rotor as large as practical when deployed], [#d.rotor_d~mm swept diameter, #box[#d.rotor_ratio ×] the egg diameter]),
      ([R3], [Real mechanics in the display build], [rotor on a 608 and a 688 ball bearing; blades fold on steel pins; magnetic detents hold both positions]),
      ([R4], [Ready for a v2 with a real generator and sealing, without redesign], [GL40 generator bay and mount, lip seal pocket, O-ring equator joint, M12 cable gland seat on a v2 front shell]),
      ([R5], [Printable on a Bambu Lab printer (#d.bed~mm cube), one material per plate], [all #d.n_stl print files fit the bed; shells print without supports; blades on a printed cradle]),
      ([R6], [One parametric source of truth], [every dimension in `ovo/params.py`; studies write JSON, the CAD reads both]),
      ([R7], [Every geometric claim checked automatically], [interference, clearances, spin, fold and closing sweeps, mesh audits (@sec-cad)]),
    ),
  ),
  caption: [Requirements and their realisation. Out of scope for v1: electronics, a water test and CFD.],
) <tab-req>

The scope is deliberately limited. The v1 is a *display build*: it shows real mechanics, but the generator
is a printed replica with the real envelope and hole patterns, and there is no electronics. All
performance numbers are therefore predictions. @sec-limits lists what they rest on and the measurements
that would confirm them.

= Concept and key decisions <sec-concept>

The concept follows from R1 and R2: a large rotor that folds into a small closed body. @tab-decisions
summarises the decisions, the alternatives that were considered and what each decision costs.

#figure(
  image("figures/concept_fold.png", width: 100%),
  caption: [Why the rotor sits downstream and its blades fold forward (schematic). (a) With the rotor at
  the nose, folded blades lie aft along the hull and the drag on their tips holds them shut; deploying would
  need the tip to move upstream. (b) With the rotor at the tail, folded blades lie forward over the tail:
  drag on the tip is meant to swing them open, and in operation the rotor thrust presses each blade onto its
  deployed stop.],
) <fig-concept>

*No duct.* A duct can raise the flow through a given rotor, but it is a rigid ring at least as large as
the rotor. It cannot fold, and it would fix the stowed diameter near #d.rotor_d~mm instead of #d.hull_d~mm.
An open axial rotor folds naturally about its hub.

*Downstream rotor, blades folding forward* (@fig-concept). Each blade hinges about a tangential pin on the
hub. Folded, it points upstream along the tail, so the flow is meant to drag its tip aft and open it;
running, the thrust holds it on the stop. Deployment is passive by design, which a water test still has to
confirm: the magnet detent that holds a stowed blade may resist the flow (@sec-limits). The load path is
benign: the dominant blade load, the flapwise moment about the hinge, pushes the blade onto its stop. The
price is mechanical: the rotor sits behind the bearings (an overhung load, @sec-bearings), and the blades
can only fold home at one rotor angle in every 120° (the stowing index), where each blade lies over a pocket.

*Pre-cone taken from the hull.* A 90° fold about a tangential axis maps "aft" in the deployed blade to
"outward" in the stowed one. A blade whose pitch axis leans aft by the local hull slope therefore hugs the
convex tail when folded. Deployed, the same lean is a downstream cone of #d.cone_mean° on average from hinge to
tip (#d.axis_cone_root° locally at the root, #d.axis_cone_tip° at the tip, where the curved blade turns back
upstream). Coning costs nothing measurable in power (@sec-blade), and downstream coning
is expected to help passive yaw alignment; that benefit has not been analysed here.

*An egg.* A Myring body of revolution (@sec-hull) gives a smooth, closed form whose tail curve doubles
as the bed for the folded blades. The egg splits at its equator into two printed shells joined by a bayonet
with an O-ring groove. Roll fins slide into rails recessed in the front shell and come off for packing, so
the stowed egg keeps its #d.hull_d~mm diameter.

#figure(
  btab((1.05fr, 1.1fr, 1.5fr, 1.35fr),
    ([Decision], [Alternatives], [Why], [Price paid]),
    (
      ([open axial rotor, 3 blades], [ducted rotor; cross-flow rotor], [folds about its hub; no rigid ring], [no flow augmentation; tip losses]),
      ([rotor downstream, blades fold forward], [rotor at the nose, blades fold aft], [flow should open folded blades (design intent); thrust holds them on the stop], [overhung rotor; fold only at the stowing index]),
      ([pre-cone from the hull curve], [straight blade, deeper pockets], [stowed blade flush with the hull], [curved pitch axis with an in-plane lead]),
      ([Myring egg split at the equator], [cylinder with end caps], [closed smooth form; one profile for shells, hub and spinner], [short body (L/D #d.ld), drag not optimised]),
      ([magnetic detents], [springs, ball detents], [no wear parts; hold deployed and stowed], [holding force only modelled]),
      ([slide-in recessed fins], [fixed fins], [stowed egg stays #d.hull_d~mm], [rail play: couple ratio spread (@sec-roll)]),
      ([constraint-aware BEM], [textbook optimum blade], [printable, foldable, cost of each constraint known], [#d.chain_total_abs~% of $C_P$ at the design point#est]),
    ),
  ),
  caption: [Key design decisions.],
) <tab-decisions>

= Hull <sec-hull>

The hull is a Myring body of revolution @myring1976: an elliptic-type nose, a short parallel band and a cubic
tail that closes to a point with a prescribed half-angle. With the nose length $a$, band $b$, tail length $c$,
diameter $D$, nose exponent $n$ and tail half-angle $theta$:

$ r(x) = D/2 [1 - ((x - a)/a)^2]^(1 slash n), quad 0 <= x < a, $ <eq-myring-nose>
$ r(x) = D/2 - ((3 D)/(2 c^2) - (tan theta)/c) xi^2 + (D/c^3 - (tan theta)/c^2) xi^3, quad xi = x - a - b, $ <eq-myring-tail>

with $r = D slash 2$ on the band. Ovo uses $a$ = #d.nose_a, $b$ = #d.mid_b, $c$ = #d.tail_c~mm, $D$ = #d.hull_d~mm,
$n$ = #d.nose_n (an elliptic nose: the wide end of an egg) and $theta$ = #d.tail_theta°, which keeps the
tail full enough to house the rear bearing and the shaft seal: the radius is still #d.r_static_tail~mm at the
static tail end (x = #d.static_tail~mm). The profile length is #d.hull_l~mm (L/D #d.ld). The built egg ends
at x = #d.x_tip~mm, where the spinner cap rounds the Myring point to R#d.spinner_tip_r (L/D #d.ld_built).

#figure(
  image("figures/hull_profile.png", width: 86%),
  caption: [Hull profile from `ovo/profile.py`. Forward of x = #d.static_tail~mm the hull is the two printed
  shells; aft of it the same curve is continued by the rotating hub and spinner, so the stowed product
  reads as one egg.],
) <fig-hull>

The body encloses #d.hull_vol~cm³ (#d.hull_vol_l~L), wets #d.hull_wet~cm² and has a frontal area of
#d.hull_frontal~cm². Its length Reynolds number at #d.design_flow~m/s is #d.hull_re × 10#super[5], in the
transitional range. The shell wall is #d.wall~mm (six perimeters of a #d.nozzle~mm nozzle). The split plane is the
equator, x = #d.equator~mm, in the middle of the parallel band, so the O-ring seals on a cylindrical
surface. The hull's drag was not computed: the tow-load estimate uses $C_D$ = #d.hull_cd on the frontal area
(#d.hull_drag~N#est at 2.5 m/s), a stated assumption. Myring's own study found the drag of such bodies fairly
insensitive to the fullness of nose and tail within a moderate range @myring1976, which is what allows the
shape to be chosen for packaging here; whether Ovo's short, full egg lies inside that range at its much lower
Reynolds number was not checked. The nose carries a printed tow eye sized for #d.tow_design~N, three times the
working load of #d.tow_working~N#est at 2.5 m/s (rotor thrust at $C_T$ = #d.ct_tow plus hull and fin drag).

= Blade design <sec-blade>

== Method

The rotor is designed with blade element momentum theory @burton2011 @hansen2015, extended for a coned
rotor. The pitch axis leans aft by the cone angle $beta$; axial momentum goes through the projected
annulus $2 pi r dif r$ and the spanwise component of the inflow is neglected. With axial and tangential
induction $a$, $a'$, tip speed ratio $lambda = Omega R slash U$, solidity $sigma = B c slash (2 pi r)$ and twist $theta_t$:

$ tan phi = (U (1 - a) cos beta)/(Omega r (1 + a')), quad alpha = phi - theta_t, $ <eq-phi>
$ a/(1 - a) = (sigma C_n cos^2 beta)/(4 F sin^2 phi), quad a'/(1 + a') = (sigma C_t)/(4 F sin phi cos phi), $ <eq-induction>
$ C_n = C_l cos phi + C_d sin phi, quad C_t = C_l sin phi - C_d cos phi, $ <eq-cncl>
$ F = 2/pi arccos(e^(-f_"tip")) dot 2/pi arccos(e^(-f_"hub")), quad f_"tip" = B/2 (R - r)/(r sin phi), quad f_"hub" = B/2 (r - R_h)/(R_h sin phi). $ <eq-prandtl>

$F$ is the product of the Prandtl tip and hub loss factors; the hub-loss radius $R_h$ is the aerodynamic root,
r = #d.r_aero_root~mm. Above $a$ = 0.4 the momentum relation is replaced by Buhl's empirical branch
@buhl2005 for the local thrust coefficient, used in the closed form of Ning @ning2014:

$ C_T = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2, quad a > 0.4. $ <eq-buhl>

Thrust, torque and the power and thrust coefficients follow by integration over #d.n_el annular elements:

$ C_P = P/(1/2 rho A U^3), quad C_T = T/(1/2 rho A U^2), quad P = Omega Q. $ <eq-coeffs>

The equations are solved by a relaxed fixed-point iteration on $(a, a')$; elements that do not converge fall
back to Ning's bracketed residual in $phi$ @ning2014. Lift and drag come from NeuralFoil @neuralfoil, a
neural network trained on XFOIL @drela1989 analyses, evaluated on the exact section the CAD builds (NACA
44xx with the blunt printed trailing edge) at each element's Reynolds number, with $n_"crit"$ = 5
for disturbed river inflow and a printed surface. At #d.design_flow~m/s the elements run at
Re #d.re_lo to #d.re_hi, a transitional regime where XFOIL-type predictions are least certain; the
sensitivity study below prices that uncertainty.

Two reference blades bracket the design. The *BEM optimum* lets every annulus maximise its own tangential
force (annuli are independent in BEM), including drag and both loss factors. The *Schmitz* blade uses the
classical closed form @gasch2012

$ phi_"opt" = 2/3 arctan(1/lambda_r), quad c = (16 pi r)/(B C_l) sin^2(phi_"opt"/2), $ <eq-schmitz>

at the lift coefficient of best lift-to-drag ratio. The *printable blade* is then a straight-tapered,
linearly twisted blade whose four parameters (root and tip chord and twist) are optimised with Nelder-Mead
for $C_P$ at the design point, inside the printability and folding bounds.

== Design tip speed ratio

For each candidate design TSR from 3 to 6 the study computes the BEM optimum, the Schmitz blade and the
printable blade (@tab-tsr). Unconstrained, the optimum is flat between TSR 4 and 4.5. The printability
bounds change the picture: the optimum chord at 0.9 R shrinks with TSR, so the #d.min_chord~mm tip-chord floor
makes high-TSR blades increasingly over-solid outboard, and below TSR 4 the root runs into the #d.max_chord~mm
chord cap (three stowed blades must fit side by side around the tail). The printable blade designed for
TSR #d.design_tsr has the highest $C_P$ at its own design point and is selected. Its $C_P$ peak sits at
TSR #d.tsr_peak at #d.design_flow~m/s.

#figure(
  btab((auto, auto, auto, auto, auto, auto, auto, auto, auto),
    ([Design TSR], [$C_P$ optimum], [$C_P$ Schmitz], [$C_P$ printable], [printable peak (TSR)], [chord root / tip (mm)], [rpm at 1.5 m/s], [torque (N~m)], [copper-loss index]),
    tb.tsr_study, align: (center, right, right, right, right, center, right, right, right), flag: true),
  caption: [Design TSR study at #d.design_flow~m/s (BEM estimates#est). The copper-loss index is
  $(Q slash Q_"TSR 4")^2$ at equal generator: lower torque at higher TSR means lower winding loss. The
  chosen design is in bold.],
) <tab-tsr>

== What the constraints cost

The printable blade differs from the BEM optimum in four ways: a blunt trailing edge of #d.te~mm (1.5 nozzle
lines, #d.te_c_root~% of the chord at the root and #d.te_c_tip~% at the tip), straight chord and twist laws, a
chord window of #d.min_chord to #d.max_chord~mm and a root-twist cap of #d.max_root_twist° from the folding
geometry, plus the #d.cone° aft cone. @tab-chain adds them one at a time, re-optimising the linear laws under
all constraints so far, so that each step's cost belongs to its constraint and not to the fitting method.

#figure(
  btab((1.6fr, auto, auto, auto, auto, auto),
    ([Step (cumulative)], [$C_P$ at TSR 4], [step], [total], [chord root / tip (mm)], [twist root / tip (°)]),
    tb.chain, align: (left, right, right, right, center, center)),
  caption: [Constraint-cost chain at the design point (#d.design_flow~m/s, TSR #(d.design_tsr)\; BEM estimates#est).],
) <tab-chain>

The blunt trailing edge is the largest single cost (#d.chain_te_abs~%): its base drag is the price of
printing the blade. The linear laws cost #d.chain_lin_abs~%, because a straight line cannot follow the
optimum where the loss factors unload it at the root and at the tip (@fig-chord). The chord window and the
twist cap cost almost nothing at TSR 4, and the cone is neutral within the model. In total the printable,
foldable blade gives up #d.chain_total_abs~% of $C_P$ (#d.cp_opt to #d.cp_final_dtsr). For scale, merely
drawing straight lines through the optimum at 0.7 R and 0.9 R, without optimising them, would have lost
#d.fit_0709_abs~% on top of the blunt edge.

#figure(
  image("/docs/figures/bem_chord_twist.png", width: 100%),
  caption: [Chord and twist along the span: BEM optimum, Schmitz, the best unbounded linear law and the final
  printable blade. The shaded band is the shank between the hinge and the aerodynamic root.],
) <fig-chord>

The final blade: #d.n_blades blades, tip radius #d.r_tip~mm, aerodynamic root at #d.r_aero_root~mm, chord
#d.chord_root to #d.chord_tip~mm and twist #d.twist_root° to #d.twist_tip° (both linear in radius), NACA 44xx
sections from #d.tc_root~% to #d.tc_tip~% thickness with the #d.te~mm trailing edge, stacked on a pitch axis
at #d.pitch_axis of the chord.

== Performance

#figure(
  image("/docs/figures/bem_cp_tsr.png", width: 100%),
  caption: [Power coefficient against tip speed ratio (BEM#est). Left: the final blade against the two
  references at #d.design_flow~m/s. Right: the final blade from 1.0 to 2.5 m/s; $C_P$ rises with Reynolds
  number.],
) <fig-cp>

#figure(
  btab((auto, auto, auto, auto, auto, auto, auto, auto),
    ([Flow (m/s)], [TSR], [rpm], [$C_P$], [shaft power (W)], [torque (N~m)], [thrust (N)], [flap moment per blade (N~m)]),
    tb.perf, align: right),
  caption: [Final blade at its $C_P$ peak (maximum power point), BEM estimates#est. The flap moment is about
  the hinge pin and pushes the blade onto its deployed stop.],
) <tab-perf>

At #d.design_flow~m/s the model gives #d.p_shaft15~W#est at #d.rpm15~rpm and #d.q15~N~m, against
#d.p_avail15~W in the swept area (@eq-pavail). The peak $C_P$ climbs from #d.cp_peak10 at 1.0 m/s to #d.cp_peak25
at 2.5 m/s as the Reynolds number rises. Unloaded, the rotor runs away to #d.run_rpm25~rpm (TSR #d.run_tsr25)
at 2.5 m/s; over all TSRs up to runaway the largest thrust at 2.5 m/s is #d.tmax25~N and the largest flap
moment #d.flapmax25~N~m per blade, the design loads for the hub, bearings and tow eye.

#figure(
  image("/docs/figures/bem_power.png", width: 82%),
  caption: [Shaft power, rotor speed, torque and thrust against flow speed (BEM#est).],
) <fig-power>

== Generator coupling and the load line <sec-loadline>

A generator does not run the rotor at its $C_P$ peak. With the hobby convention for the motor constant
$K_V$ (rpm per volt of line-to-line peak back-EMF), the torque constant and the losses are

$ K_(t 0) = 60/(2 pi K_V), quad P_"cu" = R_"ll" ((Q - Q_d)/K_(t 0))^2, quad P_"net" = (Q - Q_d) Omega - P_"cu", $ <eq-gen>

where $R_"ll"$ is the line-to-line winding resistance and $Q_d$ a constant iron and friction drag torque.
Copper loss grows with the square of the torque, so a generator with appreciable resistance prefers a higher
TSR, where the same power comes at lower torque. For the recommended CubeMars GL40 (@sec-generator) the
estimate at #d.design_flow~m/s is stark: at the $C_P$ peak the copper loss is #d.gen_peak_cu~W and the net
output #d.gen_peak_net~W#est. The best net point is at TSR #d.gen_best_tsr, #d.gen_best_rpm~rpm and
#d.gen_best_q~N~m, giving #d.gen_best_net~W#est (#d.gen_best_shaft~W shaft, #d.gen_best_cu~W copper,
#d.gen_best_drag~W drag). The operating point is therefore set by the generator, not by the blade.

#figure(
  image("/docs/figures/bem_generator.png", width: 100%),
  caption: [Left: shaft and net electrical power of the final blade at #d.design_flow~m/s against the operating
  TSR for the generator readings of the trade study (estimates#est). Right: best net power of the printable
  blade designed for each TSR.],
) <fig-gen>

This raised a design question: should the blade itself be designed for TSR 5 to 6? The right panel of
@fig-gen answers it. A blade designed for the electrical optimum gains at most a few percent of net power,
less than the polar uncertainty, while it gives up root chord (a weaker lug transition) and $C_P$ peak. TSR 4
stays the right blade, and it remains the best choice if a lower-resistance generator moves the electrical
optimum back toward the $C_P$ peak.

== Cavitation cap <sec-cavcap>

Near the free surface a fast blade tip can cavitate. The inception criterion compares the cavitation number
with the minimum pressure coefficient of the section @carlton2019:

$ sigma_v = (p_"atm" + rho g h - p_v)/(1/2 rho W^2), quad m = sigma_v + C_(p,"min") >= #d.cav_margin_cap, $ <eq-cav>

with the tip at its highest point, #d.cav_depth~m deep, $W$ the local relative speed and $C_(p,"min")$ from
NeuralFoil. The margin $m$ = #d.cav_margin_cap covers the difference between NeuralFoil and an inviscid panel
solution of the tip sections (0.02 to 0.04). At the $C_P$ peak there is no cavitation at any flow speed
(margin #d.cav_peak25 at 2.5 m/s). But the cavitation-free TSR falls quickly with flow, from #d.cav_free15 at
1.5 m/s to #d.cav_free25 at 2.5 m/s, while the generator's electrical optimum rises with flow. They cross at
#d.cap_first_flow~m/s (@fig-loadline, @tab-loadline).

#figure(
  image("/docs/figures/bem_load_line.png", width: 100%),
  caption: [GL40 load line against flow speed with the cavitation cap (estimates#est). Top left: operating TSR,
  cavitation limits and runaway. Top right: net electrical power, uncapped, capped and with a single rpm
  ceiling. Bottom: the cap in the rotor speed and shaft torque plane, where a controller can apply it.],
) <fig-loadline>

#figure(
  btab((auto, auto, auto, auto, auto, auto, auto, auto, auto),
    ([Flow (m/s)], [$C_P$-peak TSR], [optimum TSR (rpm)], [net (W)], [cav. margin], [cav.-free to TSR], [capped TSR (rpm)], [capped net (W)], [cap cost]),
    tb.load_line, align: right),
  caption: [GL40 load line on the final blade, tip #d.cav_depth~m deep (estimates#est).],
) <tab-loadline>

The recommended control law tracks the electrical optimum but never exceeds the cap. It changes nothing up
to the crossing and costs #d.cap_cost25~% of net power at 2.5 m/s (#d.ll_net25_cap~W instead of
#d.ll_net25_un~W#est, at TSR #d.ll_tsr25_cap and #d.ll_rpm25_cap~rpm). The controller does not know the flow
speed, but the cap can be written without it, as a minimum shaft torque at a given speed: above
#d.cap_rpm0~rpm the shaft torque must stay above a curve rising from #d.cap_q0~N~m to #d.cap_q1~N~m at
#d.cap_rpm1~rpm. At a fixed speed less torque means less flow, a higher TSR and a stronger suction peak on
the tip's pressure side. This rule classifies all #d.cap_checks sweep points from the $C_P$ peak to runaway
correctly (#d.cap_bad misclassified). A simple rpm ceiling would be a poor substitute: it would cost
#d.ceil_cost25~% at 2.5 m/s. The rotor must never be left to free-wheel when towed at #d.cav_freewheel_flow~m/s or more; it
must stay loaded or be braked electrically (@sec-generator).

== Loads, structure and the hinge geometry

The flap moment about the hinge is the dominant blade load and it acts in the direction the folding concept
needs. In the #d.lug~mm lug that carries it the bending stress is #d.lug_stress15~MPa#est at #d.design_flow~m/s
and #d.lug_stress25~MPa#est at 2.5 m/s (bending across the 8 mm depth, no stress concentration). The blade
prints with its span in the layer plane so this stress runs along the filament. The study first used a
straight #d.cone° pitch axis; the layout (@sec-fold) then solved a curved, hull-following axis. Flying the
same blade on that axis leaves the $C_P$ peak unchanged (#d.cp_layout25 against #d.cp_nominal25 at
2.5 m/s#est) and lowers the largest flap moment from #d.flap_nominal25 to #d.flap_layout25~N~m#est\; the hub
and stop are sized with the latter.

== Verification and sensitivity

The solver was checked against theory and against itself:
- *Glauert optimum rotor* (ideal polar, no losses, TSR 4): BEM $C_P$ #d.glauert_bem against the theoretical
  #d.glauert_th.
- *Betz limit*: the largest $C_P$ in any sweep is #d.cp_max_any, below 0.593.
- *Induction*: mid-span axial induction #d.a_mid_lo to #d.a_mid_hi at the design point, near the ideal 1/3;
  the Buhl branch acts on #d.n_buhl of #d.n_el elements, at the root and tip where $F$ is small.
- *Convergence audit*: #d.conv_unconv unconverged element solves out of #d.conv_solves (#d.conv_ops operating
  points), all in the two root elements, which carry at most #d.conv_root_share~% of the torque in any
  power-producing point. None is at a design or maximum-power point of the final blade.
- *Polar confidence*: NeuralFoil's own confidence is at least #d.nf_conf over the design-point elements; the
  angle of attack at mid-span is #d.alpha_lo° to #d.alpha_hi°, at least #d.stall_margin° from stall.

#figure(
  btab((1.8fr, auto, auto, auto),
    ([Case (peak $C_P$ at #d.design_flow~m/s)], [$C_P$ peak], [TSR], [change]),
    tb.sens, align: (left, right, right, right), markup: true),
  caption: [Sensitivity of the peak power coefficient (BEM estimates#est).],
) <tab-sens>

The sensitivity study (@tab-sens) shows where the uncertainty lies. Cone angle, the curved axis and the
element count move $C_P$ by 0.2 % or less. Transition modelling dominates: a clean inflow ($n_"crit"$ = 9)
and forced transition at 5 % chord (a rough print) cost #d.sens_ncrit9_abs~% and #d.sens_forced_abs~%. The nominal
$C_P$ should therefore be read as #d.cp_peak15 with an uncertainty of about 6 % from the polars alone, before
any effect of the hull, which the model ignores (no blockage, no hull wake at the rotor plane). The hub-loss
placement is conservative: without it $C_P$ would be #d.sens_hub~% higher.

= Fold kinematics, pockets and detents <sec-fold>

== Hinge and the stowed blade

The hinge position was chosen by a sweep over its axial station, its radius and the end of the static tail
(`ovo/layout.py`). The rule admits only geometries in which the pin keeps at least 2 mm of hub material
around it, the folded heel stays inside the spinner, the heel face holds a magnet, the hub wall over the
shaft bore keeps #d.min_wall~mm and the folding root clears the static tail by 0.8 mm; among those it maximises
the smallest margin. The result is a pin at x = #d.hinge_x~mm, r = #d.hinge_r~mm, with the static tail ending at
x = #d.static_tail~mm: pin cover #d.pin_cover~mm, heel #d.heel_inside~mm inside the spinner when stowed,
#d.heel_face~mm of flat heel face, hub wall #d.hub_wall~mm. The blade folds #d.fold_deg° and its tip lands at
x = #d.stow_tip_x~mm.

The stowed shape is solved rather than drawn. For every station the aft offset of the pitch axis is chosen so
that the folded section's top sits on the hull surface: the stowed blade stands at most #d.stow_flush~mm proud
of the hull. An in-plane lead of the stacking axis, #d.lead_root~mm at the root falling to #d.lead_tip~mm at
the tip, shifts the sections along the hull so their twist does not dig deep pockets; the BEM model ignores this
lead (@sec-limits). The deployed cone that results averages #d.cone_mean° from hinge to tip, close to the
#d.cone° the blade was designed with.

#figure(
  image("img/fold_sequence.jpg", width: 100%),
  caption: [Fold sequence rendered from the CAD assembly at the stowing index: stowed, 60°, 30° and deployed
  (fins fitted; camera fixed on the static shells).],
) <fig-foldseq>

#figure(
  image("figures/fold_side_main.png", width: 92%),
  caption: [Fold kinematics in the plane of one blade (layout study, `docs/figures/layout_fold_side.png`):
  deployed, outlines at 30° and 60°, and stowed over the pocket floor (dashed) and its keep-out (dotted).],
) <fig-fold>

#figure(
  image("/docs/figures/layout_fold_end.png", width: 100%),
  caption: [The rotor seen from aft: stowed, the three blades stay inside the egg outline; deployed, they sweep
  #d.rotor_d~mm. Pins are tangential, so each blade folds in its own meridional plane.],
) <fig-foldend>

== Pockets

Three pockets in the static hull, from x = #d.pocket_x0 to #d.pocket_x1~mm, take the folded blades. They are
#d.pocket_eq~mm deep at the equator, at most #d.pocket_aero~mm deep under the aerodynamic span and
#d.pocket_max~mm in the root notch, with at least #d.pocket_land~mm of hull land between neighbours. Because
both shells print with +x down, each pocket floor is cut so that it rises at most #d.max_overhang° going
forward, and the pocket's forward end and the root-notch ramp are #d.max_overhang° ramps: the shells print
without supports. On the CAD solids the stowed aero span clears its pocket by #d.aero_clear~mm and the root
clears the static tail by #d.root_clear~mm (@sec-cad). The pockets also index the rotor in the bag: stowed,
it cannot turn.

== Detents and hub loads

Magnets replace springs. Deployed, an N52 block (6 × 2 × 2 mm) in the blade heel faces one in the hub's stop
face across #d.det_gap~mm, about #d.det_force~N#est\; stowed, a 4 × 2 mm disc in the blade shank faces one in a
boss of the aft shell, #d.stow_gap~mm apart and #d.stow_ang° out of line, about #d.stow_force~N#est. Forces
come from a dipole model scaled to the catalogue pull of the 4 × 2 disc (factor #d.det_k); they must be
confirmed on a test print. At 2.5 m/s the deployed stop carries #d.stop_force25~N#est at
#d.stop_p25~MPa contact pressure, the 3 mm pin #d.pin_shear25~MPa in double shear and the printed ears
#d.ear_bearing25~MPa in bearing, all from the curved-axis flap moment.

The blades cannot be folded home with the shells being closed over them. A closing sweep on meshes of the CAD
solids (bayonet entry, push and #d.bay_twist° twist, both senses) shows that the front shell passes the blades only
if they are folded no further than #d.close_fins_off° (opened at least #(90 - int(d.close_fins_off))°) with the fins
off, or #d.close_fins_on° (opened at least #(90 - int(d.close_fins_on))°) with the fins on; the assembly sequence
encodes this.

= Layout and drivetrain <sec-layout>

== Axial stack

@fig-section shows the merged layout. The generator sits in the front shell on a twist-lock carrier; the
shaft runs from a claw coupling at the equator through two bearings in a tube of the aft shell to the hub.
The main stations are: generator bay D#d.bay_d × L#d.bay_l at x = #d.bay_x0 to #d.bay_x1, equator joint at
#d.equator, 608 bearing at #d.b608_x0 to #d.b608_x1, 688 bearing at #d.b688_x0 to #d.b688_x1, lip seal at
#d.seal_x0 to #d.seal_x1, static tail end at #d.static_tail, hub from #d.hub_front to the shaft end at
#d.shaft_end~mm. The two halves of the layout (rotor and internal) are checked against each other by
#d.interface_n interface checks; #d.interface_pass pass.

#figure(
  image("figures/layout_section_drawing.png", width: 100%),
  caption: [Merged layout, side section (from `docs/figures/layout_section.png`). Top half in the plane of blade 0
  with the deepest pocket floor and the stowed blade (salmon); bottom half between pockets. Blue: GL40 generator
  in its D#d.bay_d × L#d.bay_l bay (dashed); purple: twist-lock carrier and adapter; green: claw coupling with
  the thrust collar; charcoal: 608 and 688 bearings; brown: lip seal; tan: socket and spigot with the O-ring
  (black); light grey: bearing tube and ribs; orange: blade knuckle on its hinge pin. The section shows the v2
  nose with the M12 cable gland; the display build has a clean nose with the tow eye only.],
) <fig-section>

#figure(
  image("img/cutaway.jpg", width: 80%),
  caption: [Cutaway render of the CAD assembly: generator replica (blue) on its carrier, claw coupling, shaft,
  608 and 688 bearings, lip seal, and the v2 cable gland in the nose.],
) <fig-cutaway>

== Bearings, shaft and thrust path <sec-bearings>

The rotor is downstream, so its thrust pulls the hub and shaft aft. A set-screw collar (Ruland MSC-8-SS)
inside claw B bears on the inner race of the forward 608-2RS, which locates the shaft; the rear 688-2RS
floats. Their bores step from 22 to 18 to 16 mm so both press in from the front on crush ribs, and the 16 mm
seat of the 688 is also the lip-seal bore. The shaft is an 8 mm SUS440C rod (MISUMI PSSFJ class, about 56 HRC,
g6) @misumi, hard enough for the seal lip.

The bearing span is #d.span~mm against a #d.overhang~mm overhang to the rotor's thrust centroid, a ratio of
#d.span_ratio where a common rule of thumb asks for #d.crit_ratio. The static tail cannot be lengthened without failing
the hinge sweep, so the rule was replaced by the explicit checks of @tab-bearings for the tow case: thrust
#d.thrust_tow~N#est from $C_T$ = #d.ct_tow at 2.5 m/s, a radial load of #d.radial_pct~% of it and one blade
#d.asym_pct~% off its share, all assumptions.

#figure(
  btab((1.6fr, auto, auto),
    ([Criterion], [limit], [value]),
    (
      ([bearing span], [at least 4 shaft diameters, #d.crit_span~mm], [#d.span~mm]),
      ([blade-tip axial runout from #d.crit_slop~mm seat slop], [at most #d.crit_runout~mm], [#d.tip_runout~mm]),
      ([static safety $C_0 slash P_0$, 608 / 688], [at least #d.crit_s0], [#d.s0_front / #d.s0_rear]),
      ([shaft deflection at the load plane], [at most #d.crit_defl~mm], [#d.shaft_defl~mm]),
      ([span / overhang], [#d.crit_ratio (rule of thumb, reported)], [#d.span_ratio]),
    ),
    align: (left, left, right)),
  caption: [Bearing checks replacing the span-to-overhang rule (layout study; loads are estimates#est).],
) <tab-bearings>

The loads are small for these bearings: their static ratings exceed the largest thrust at 2.5 m/s
#d.c0_ratio608 times (608, $C_0$ = #d.c0_608~N @bearing608) and #d.c0_ratio688 times (688 @bearing688). The open question is
stiffness, not strength: rotor wobble should be measured on the first print (@sec-limits).

== Seals, joint and bayonet

The v2 shaft seal is a TC 8×16×5 NBR lip seal fitted 0.5 mm inside the static tail face. Its lip runs at
#d.seal_speed~m/s at 750 rpm, so wear is negligible; the installation guideline asks for a shaft of at least
#d.seal_hrc HRC @lipseal_install, which the hardened rod meets.

The equator joint is a radial static gland. The aft shell carries a spigot with the O-ring groove, the front
shell a socket bore. The ring is an ISO 3601-141 (AS568-141, #d.or_id × #d.or_cs~mm, EPDM 70 for water)
@applerubber, the largest 2.62 mm section whose socket keeps at least #d.min_wall~mm of wall under the blade
pockets with a 1 mm reserve on the pocket depth: the wall under the pockets is #d.or_wall~mm. Groove bottom
D#d.or_groove, spigot OD #d.or_spigot, bore D#d.or_bore, radial gland depth #d.or_gland~mm. A rectangular
groove would leave an unprintable flat overhang on the spigot, which prints tip-up, so the groove has a
radial aft flank and a 45° tip-side flank (bottom #d.or_bw, top #d.or_tw~mm). With squeeze $s = (d_2 - h) slash d_2$
($d_2$ the cross-section, $h$ the gland depth) and fill the ring area over the groove area, the worst case
combines the ring tolerances, printed diameters within ±#d.print_tol~mm and full eccentricity:

- squeeze #d.or_sq_nom~% nominal, #d.or_sq_lo~% to #d.or_sq_hi~% worst case (guidance #d.og_sq0 to #d.og_sq1~% @oring_guide);
- fill #d.or_fill_nom~% nominal, #d.or_fill_worst~% worst case (guidance #d.og_fill0 to #d.og_fill1~%);
- stretch #d.or_stretch_lo~% to #d.or_stretch_hi~% on the groove bottom.

The worst-case minimum squeeze is below the #d.tab_sq0~% minimum of the gland table @applerubber. It is the price of printed tolerances; at a few
centimetres of submersion the joint sees well under 0.1 bar, and the squeeze must be checked on a printed
test ring. The layout study recomputes the gland independently and agrees with the hardware study to
#d.oring_dev~mm. Three bayonet lugs lock the shells after a #d.bay_twist° twist against a hard stop that also aligns
the pockets across the joint; three flexure tines in the spigot wall click into dimples at the lock angle.

#figure(
  image("figures/joint_section.png", width: 92%),
  caption: [Equator joint section with the -141 O-ring, the bayonet lug and the socket (layout study).],
) <fig-joint>

#figure(
  image("figures/generator_bay.png", width: 86%),
  caption: [Generator bay frontier: the largest bay diameter for each length (set by the nose wall and the
  carrier ring), the chosen D#d.bay_d × L#d.bay_l bay and the catalogue motors; the shaded candidate box spans
  D 40 to 56 mm and L 20 to 35 mm (layout study,
  `data/layout_internal.json`).],
) <fig-bay>

== Generator trade study and the GL40 <sec-generator>

Although v1 carries no electronics, the bay is sized for a real generator, so a v2 needs no redesign. The
study evaluated #d.n_candidates off-the-shelf motors (@tab-gen): gimbal outrunners, low-KV drone motors, a frameless torque
motor and an industrial flat motor. Each was coupled to the BEM rotor and its operating point searched for
maximum DC power, with three rectifier models: an analytic passive six-pulse bridge (block currents), an ideal
active rectifier, and for motors with a published inductance a time-domain bridge model with a capacitor
link. Iron and friction losses enter as a constant drag torque from the published no-load input; where none
is published the GM4108H value is assumed.

#figure(
  btab((1.45fr, auto, auto, auto, auto, auto, auto, auto, auto),
    ([Candidate], [D × L (mm)], [$K_V$], [$R_"ll"$ (Ω)], [$K_m$ (N m/√W)], [bridge (W)], [active (W)], [fits bay], [price]),
    tb.gen_candidates, align: (left, center, right, right, right, right, right, center, right), flag: true, size: 7.6pt),
  caption: [Generator candidates coupled to the rotor at #d.design_flow~m/s (estimates#est\; bridge: analytic
  passive bridge; active: ideal active rectifier). Ranges cover the conflicting datasheet readings. The bay
  is the D#d.bay_d × L#d.bay_l envelope of the layout. Recommended choice in bold.],
) <tab-gen>

*Why the GL40 KV70.* It is the only in-envelope candidate whose datasheet is complete and self-consistent
@gl40_page @gl40_drawing: $K_V$ #d.gl_kv rpm/V and the stated torque constant agree ($K_(t 0)$ = #d.gl_kt0~N m/A),
the rated #d.gl_trated~N~m at #d.gl_irated~A matches, resistance (#d.gl_r~Ω) and inductance (#d.gl_l~mH) are
given line to line, and a manufacturer drawing gives both mounting faces. Coupled to the rotor at
#d.design_flow~m/s it delivers #d.gl_td15~W#est DC through a passive diode bridge (time-domain model,
#d.gl_td15_v~V, #d.gl_td15_i~A, TSR #d.gl_td15_tsr) and #d.gl_act15~W#est with active rectification, at
#d.gl_load_frac~% of its rated torque. Its #d.gl_d × #d.gl_len~mm body and 8 mm through bore fit the bay
and the shaft.

*An open conflict on the towing range.* At 2.5 m/s the electrical optimum would exceed the GL40's continuous
rating. Held at the rated current (#d.gl_irated~A), the rotor would run at TSR #d.gl_td25_cap_tsr and give
#d.gl_td25_cap~W#est, but at 2.5 m/s the tips are cavitation-free only up to TSR #d.cav_free25
(@sec-cavcap): the rated-current point cavitates. The cavitation-capped load line avoids that by loading the
rotor harder, and from about #d.conflict_flow~m/s it needs more current than the rating allows
(#d.conflict_i~A, #d.conflict_x × the rated copper loss); at 2.5 m/s it needs #d.ll_i25_cap~A, inside the
#d.gl_ipeak_r~A peak rating but #d.ll_x25_cap × the rated copper loss#est. Above about #d.conflict_flow~m/s
the GL40 therefore cannot stay inside its continuous rating and keep the tips cavitation-free at the same
time. This is an open v2 issue (@sec-limits).

*Honest uncertainty.* The ranking among the in-envelope motors is not sharp: with $K_m$ between 0.07 and 0.11
they land within a few watts of each other, and what separates them is voltage, current and data quality.
The cheaper iPower GM4108H-120T (USD #d.gm_price against USD #d.gl_price) could deliver #d.gm_p_lo to
#d.gm_p_hi~W#est, but its product page @gm4108_page contradicts itself: the no-load line implies $K_V$
#d.gm_kv_lo to #d.gm_kv_hi, which predicts #d.gm_ratio_lo to #d.gm_ratio_hi times the load torque listed on the same
page; if the load point is right, its output at #d.design_flow~m/s collapses to #d.gm_collapse_lo to
#d.gm_collapse_hi~W#est. The GL40 has its own open points: no no-load loss is published (its drag torque of
#d.gl_drag~N~m is assumed equal to the GM4108H value), and which part of it carries the 8 mm bore is inferred
from the drawing. A ten-minute bench test (line resistance, line inductance, $K_V$ from the open-circuit
voltage at a measured speed, spin-down) settles all of it before v2 is frozen.

*Towing state.* Electrically, an unloaded rotor is harmless: the GL40 reaches #d.gl_voc25~V#est line-to-line
peak at the 2.5 m/s runaway speed (the GM4108H #d.gm_voc_lo to #d.gm_voc_hi~V), so a charger input rated for
#d.input_v~V survives open phases. Hydrodynamically it is not: a free-wheeling rotor cavitates from
#d.cav_freewheel_flow~m/s (@sec-cavcap) and lets the fins roll the hull (@sec-roll). When towed, the charger
should therefore keep the rotor loaded. The hardware study sizes a #d.dump_r~Ω dump resistor across the DC link
for a full battery, absorbing #d.dump_p~W#est at 2.5 m/s outside the motor; it is sized for the rated-current
point (TSR #d.gl_td25_cap_tsr), which lies in the cavitating region, so it inherits the conflict above and is
not a settled design. A continuous phase short is not an acceptable
brake: at 2.5 m/s it would dissipate #d.gl_short25~W#est in the windings, #d.gl_short25_x × the rated copper
loss, inside a sealed hull.

= Roll stability and fins <sec-roll>

The generator's reaction torque acts on the hull, which hangs on an axial tow line and has almost nothing to
resist roll: its centre of gravity is within a fraction of a millimetre of the axis, so the pendulum
righting moment is at most #d.roll_pend~N~m against a generator torque of #d.roll_q15~N~m#est at
#d.design_flow~m/s. Two fins with opposite incidence supply a counter-torque:

$ Q_f = 2 dot 1/2 rho V_f^2 S C_L r_"cp" , quad Q_f slash Q_"gen" in [#d.fin_lim0, #d.fin_lim1] "at" #d.design_flow "m/s", $ <eq-fin>

with $V_f$ the local flow speed at the fin, $S$ the area of one fin and $r_"cp"$ the radius of its centre of
pressure. The band asks for 10 % more torque than the GL40 load line needs, a margin for a printed Clark Y
at low Reynolds number and for the estimated generator drag. The fin study (`ovo/layout_internal.py`) picks
the smallest tip-to-tip span that meets it and clears the pockets and the folding blades: two Clark Y fins of
#d.fin_chord × #d.fin_span~mm (#d.fin_area~cm² each, aspect ratio #d.fin_ar) at θ = #d.fin_theta0° and
#d.fin_theta1°, set at #d.fin_inc° incidence, #d.fin_tip2tip~mm tip to tip in the study (#d.fin_plan~mm
across the fins in plan as built). Their couple is #d.fin_torque~N~m#est,
#d.fin_ratio × the load-line torque, with $C_L$ #d.fin_cl15 against an estimated three-dimensional maximum of
#d.fin_clmax15 at Re #d.fin_re15.

#figure(
  image("figures/fins_plan.png", width: 74%),
  caption: [Roll fins from above with the blade fold torus (layout study, `docs/figures/layout_section_fins.png`,
  which also shows the view from aft). Folding at the stowing index, the blades clear the recessed rail slots by
  at least #d.fin_slot_fold_clear~mm and the fins by #d.fin_fin_fold_clear~mm in 3D.],
) <fig-fins>

A fixed fin set cannot balance every flow speed: the load line raises the TSR with flow, so the generator
torque does not grow with $V^2$ as the fin lift does. The couple ratio is #d.fin_ratio10 at 1.0 m/s and
#d.fin_ratio25 at 2.5 m/s#est, and since nothing restores roll, any mismatch becomes a slow hull rotation:
#d.roll_rate15~rpm#est at #d.design_flow~m/s, #d.roll_rate25~rpm#est at 2.5 m/s. With no electrical load the
fins alone would roll the hull at #d.freewheel_roll~rpm#est. Two consequences for v2: the generator must stay
loaded, and the tow line needs a roll reference, a two-leg bridle on two nose eyes (about #d.bridle_heel15°
of heel#est at #d.design_flow~m/s) or a swivel.

The fins slide aft into half-dovetail rails sunk into bosses on the inside of the front shell wall, so nothing
stands proud of the hull, and an N52 block pair holds each one. Flush blanking strips close the rails when
the fins are off. As built, the tongue has #d.clear_slide~mm side clearance, and the measured yaw play of a seated fin
(+#d.fin_yaw_p° / −#d.fin_yaw_m°) spreads the couple ratio over #d.fin_range0 to #d.fin_range1#est around the
nominal #d.fin_ratio_built: outside the band at both ends. That has no effect on the display build; for v2,
crush ribs on the tongue would locate it. Printed tip-down, the fin root sees #d.fin_root_stress~MPa#est in
the tow case, #d.fin_root_util of an assumed #d.fin_interlayer~MPa interlayer strength.

= CAD architecture and verification <sec-cad>

== Code-CAD from a single source of truth

Ovo is modelled entirely in Python with build123d @build123d, a boundary-representation library on the
OpenCascade kernel. There is no hand-edited CAD file. Every dimension lives in `ovo/params.py`; the design
studies read it and write their results as JSON; the CAD modules read both; one script builds, checks and
exports everything (@fig-pipeline). Deviations found during CAD reviews are written back into the
parameter file, so parameters and solids agree.

#let node(t, c: tint) = box(fill: c, stroke: 0.5pt + rule-c, inset: (x: 5pt, y: 5pt), radius: 2pt, text(size: 8pt, t))
#let arr = text(size: 11pt, fill: ink2, sym.arrow.r)
#figure(
  kind: image,
  block(width: 100%, inset: (y: 4pt))[
    #set par(justify: false)
    #grid(columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr), align: horizon + center, gutter: 4pt,
      node[*params.py*\ design brief], arr,
      node[*studies*\ `bem.py`, `layout.py`,\ `layout_internal.py`,\ hardware calc], arr,
      node[*data/\*.json*\ blade stations, pockets,\ axis, hardware with sources], arr,
      node[*CAD modules*\ shells, blade, rotor,\ drivetrain, fins, coupon],
    )
    #v(4pt)
    #grid(columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr), align: horizon + center, gutter: 4pt,
      node[*report*\ `report.py` → this PDF], text(size: 11pt, fill: ink2, sym.arrow.l),
      node[*outputs*\ STEP, STL, GLB,\ renders, viewer], text(size: 11pt, fill: ink2, sym.arrow.l),
      node[*checks JSON*\ interference, clearances,\ spin, meshes, print plan], text(size: 11pt, fill: ink2, sym.arrow.l),
      node[*assembly*\ #d.n_parts_total parts, states,\ poses, BOM],
    )
  ],
  caption: [Build pipeline. `scripts/build.py` runs the six part modules, the assembly with its checks and
  exports, and the renders; `scripts/report.py` turns the JSON outputs into this report.],
) <fig-pipeline>

The part modules are independent: each builds its solids from the parameters and the study tables, runs its
own checks and writes them to `out/cad/checks/<module>.json`. The assembly module collects the parts, poses
them in each state, adds the detent magnets that no builder models, and runs the integration checks. The
blades carry their hinge as data (origin, axis and deployed angle), so the same solid is posed stowed,
folded or deployed by one rotation; the integrator checks that this rotation reproduces the deployed pose to
10#super[#d.pivot_exp] in the transform matrices. A part cache keyed by a SHA-1 hash of the part modules, the parameter
and profile files, the layout code and every `data/*.json` file makes stale geometry impossible to reuse, and a
`--fresh` build deletes every output the current parts no longer produce. The studies are deterministic in the
same sense: re-running the rotor layout reproduced its JSON byte for byte.

== Parts and states

The assembly has #d.n_parts_total parts, #d.n_printed_total printed (from #d.n_stl STL files: the three blades
share one) and #d.n_bought_total hardware parts (bought items and their sub-parts): printed shells, hub, spinner and cap, three blades, fins and their
blanking strips, a printed claw coupling, motor carrier and generator replica, and the bought parts
(bearings, shaft, collar, pins, screws, magnets, O-ring, lip seal, cable gland). Seven states are modelled:
stowed, deployed, folded 30° and 60°, cutaway (a quarter of the shells removed), packed (fins off, strips in)
and v2 (front shell with the gland seat and the gland). The display and v2 builds differ only in the front
shell: the display shell has a clean nose with the tow eye, the v2 shell adds the M12 gland seat and tap
bore.

#if d.has_details [
  #figure(
    grid(columns: (1fr, 1fr), gutter: 6pt,
      image("img/details_hub.jpg", width: 100%),
      image("img/details_fin_rail.jpg", width: 100%),
    ),
    caption: [Details of the CAD model (renders, staged poses). Left: the hub with the spinner, its cap and
    screw hidden and the blades part-folded, showing the ears, the tangential hinge pins, the shaft bore and the
    stop faces that carry the deployed-detent magnets. Right: a fin slid partly out of its recessed rail to
    show the half-dovetail slot and the floor magnet; seated, nothing stands proud of the hull (the rail check
    finds #d.rail_proud~mm³ of shell more than 0.05 mm outside it).],
  ) <fig-details>
]

#figure(
  image("img/exploded.jpg", width: 86%),
  caption: [Exploded view in assembly order: front shell and fins, generator stack and claw A, claw B and
  bearings, aft shell with the shaft, hub with the blades and pins, spinner and cap.],
) <fig-exploded>

== Verification

Clearances and overlaps are measured on fine meshes of the CAD solids with manifold3d @manifold booleans and
minimum-distance queries (#d.mesh_chord~mm chord, #d.mesh_chord_shell~mm on the shells). In every state, every pair of parts whose
bounding boxes come within 3 mm is tested (@tab-interference). Intended overlaps (press fits on crush ribs,
pins in printed ears, screws in their threads) are listed in the code with their reason and reported
separately; anything else above #d.overlap_tol~mm³ is a defect.

#figure(
  btab((1.4fr, auto, auto, auto, auto, auto, auto),
    ([State], [fold (°)], [parts], [pairs within 3 mm], [intended overlaps], [zero-gap contacts], [defects]),
    tb.interference, align: (left, right, right, right, right, right, right)),
  caption: [Interference matrix by assembly state (`out/cad/checks/assembly.json`). Static pairs are computed
  once in the stowed state; the other states add the pairs that involve a part that moves or changes.],
) <tab-interference>

#figure(
  btab((1.55fr, 1.6fr, 1.25fr, auto),
    ([Check (CAD meshes)], [measured], [limit], [result]),
    tb.clearances, align: (left, left, left, center), size: 7.6pt),
  caption: [Selected key clearances (#d.kc_rows_note checks; all #d.kc_ok of #d.kc_n pass).],
) <tab-clear>

Further automatic checks:
- *Spin sweep.* The rotating parts are turned through 360° in #d.spin_step° steps against every static part
  nearby: #d.spin_pairs pairs, #d.spin_rub rubbing; the closest running gap is #d.spin_min~mm
  (#d.spin_min_pair).
- *Fold and closing sweeps.* Blades are folded from 0 to 90° in 3° steps against the hub, spinner, pins, shells
  and neighbours, and against the fins at the stowing index; the shells are closed over folded blades along
  the real bayonet path (@sec-fold).
- *Print meshes.* The STL audit covers #d.n_stl_entries printed parts in #d.n_stl files (the three
  blades share `blade.stl`): #d.stl_watertight of #d.n_stl_entries are watertight with consistent winding,
  have #d.stl_folded folded sliver pairs in all, match their B-rep volume within 1 %
  (#d.stl_volmatch / #d.n_stl_entries) and fit the bed (#d.stl_fit / #d.n_stl_entries).
- *Web model.* The viewer scene has #d.glb_meshes meshes and #d.glb_tris triangles (#d.glb_mb MB) with
  #d.glb_folded folded pairs; a headless browser check confirms that the viewer's folded and deployed blade
  poses match the CAD within #d.viewer_diff~mm.
- *Module checks.* #d.mod_total_ok of #d.mod_total module interface and feature checks pass
  (#d.mod_checks), on top of the #d.interface_pass layout interface checks.

These checks found real defects during development: blades classified inside-out by the B-rep classifier
because of a hooked surface parameterisation, a gland boss that reached into the generator bay, a fin pad
that overlapped a pocket once its yaw was modelled, and a coupon check that read the wrong part. Each was
fixed at its root and is now guarded by a check.

= Printability and print plan <sec-print>

The design rules are those of a #d.nozzle~mm nozzle at #d.layer~mm layers: walls of at least #d.min_wall~mm (three
perimeters), overhangs of at most #d.max_overhang° from vertical without support, and fits taken from printing
guides @prusa_fits @creative3dp: #d.press_fit~mm per side for press fits on three crush ribs, #d.clear_slide~mm per
side for sliding fits, #d.clear_rotate~mm for the blade knuckle on its pin and #d.run_clear~mm between moving
bodies. Two tolerance coupons are printed first, in each material and profile; they carry the bearing
seats, pin holes, magnet seats in the three print orientations, the hub D-bore, an M12 tap drill and keys for
the bayonet and the fin rail, and their result adjusts the fits before the real parts are printed.

Every part has a planned orientation. Both shells print with +x down without supports: the pocket floors and
ramps respect the 45° rule, the O-ring groove has its 45° flank, and the only large bridge is a #d.nose_bridge~mm
disc inside the nose (the three bayonet-click tines of the aft shell print with short breakaway ties across their
slits). The blades print chord-vertical with the trailing edge down on a sacrificial printed cradle whose
rails stand #d.cradle_gap~mm under the blade, which holds the #d.blade_te_min~mm trailing edge straight without
slicer supports. Display parts are PLA; parts that carry torque, slide or are clamped by a set screw, and the
fins, are PETG, the material recommended for wetted parts @prusa_petg. A v2 wet build moves the shells and rotor
to PETG.

#figure(
  btab((auto, 1.3fr, 2.4fr, auto, auto),
    ([Plate], [content], [parts], [max height (mm)], [filament (g)]),
    tb.plates, align: (center, left, left, right, right), size: 7.6pt),
  caption: [Print plates on the #d.bed~mm bed (`out/cad/checks/assembly.json`). Filament masses are estimates
  from the shell area, the wall count and the infill.],
) <tab-plates>

One display build needs about #d.print_g~g of filament (#d.print_pla~g PLA, #d.print_petg~g PETG#est), plus brims
and the coupons; the cutaway shells and the v2 front shell are alternatives outside that total.

#figure(
  grid(columns: (1fr, 1fr), gutter: 6pt,
    image("img/stowed.jpg", width: 100%),
    image("img/packed_kit.jpg", width: 100%),
  ),
  caption: [Packed state. Left: blades in their pockets, fins removed and blanking strips in place. Right: the
  packed egg with the two loose fins.],
) <fig-packed>

= Bill of materials and cost <sec-bom>

@tab-bom lists the bought parts of the as-built assembly, grouped, with unit prices where the hardware study
sourced one (USD converted at #d.fx USD per EUR, ECB reference rate of #d.fx_date @ecb). Printed parts and
filament are not priced.

#figure(
  btab((2.2fr, auto, auto, auto, auto, 1.3fr),
    ([Bought part], [qty], [unit (EUR)], [total (EUR)], [build], [price basis]),
    tb.bom, align: (left, right, right, right, left, left), size: 7.6pt),
  caption: [Bought parts of the as-built assembly (`out/cad/checks/assembly.json`), priced where a source
  exists (`data/hardware.json`).],
) <tab-bom>

The priced display-build lines of this as-built list total EUR #d.bom_asbuilt_known pro rata;
#d.bom_asbuilt_unpriced lines have no sourced price yet (pins, block magnets, shaft, collar, shim,
screws, O-ring, lip seal and the seal-fitting sleeve). The hardware study's own bill of materials, drawn up before the CAD, gives EUR
#d.bom_disp_pro pro rata and EUR #d.bom_disp_buy cash outlay for the display build (minimum orders: 50
magnets, 100 inserts) and EUR #d.bom_v2_pro / #d.bom_v2_buy for the v2 additions, with #d.bom_unpriced lines
unpriced. It predates three CAD changes: one 608 and one 688 instead of two 608s, #d.n_magnets magnets of
three kinds instead of sixteen discs, and two inserts instead of eight. The generator dominates the v2 cost.

= Prior art <sec-prior>

@tab-prior compares Ovo with the portable and towed water generators for which data could be found. Few
publish the flow speed at their rated power, so the last column inverts @eq-pavail: the flow each rotor
needs for its rated power at an electrical power coefficient of 0.20 to 0.30.

#figure(
  btab((0.95fr, auto, 1.15fr, auto, auto, 1.75fr, 0.85fr, auto),
    ([Product], [rotor D], [rated power], [mass], [stowed / body], [stowage], [price], [flow for rating (m/s)]),
    tb.prior_art + ((
      strong[Ovo (this work)], [#d.rotor_d~mm], [#d.gl_td15~W DC#est at #d.design_flow~m/s (GL40, bridge)], [n/a], [#d.hull_d × #d.x_tip_r~mm],
      [blades fold forward into tail pockets], [not a product (@sec-bom)], [#(d.ovo_vneed_dc_lo)-#d.ovo_vneed_dc_hi],
    ),),
    align: (left, right, left, right, left, left, left, right), size: 7.4pt),
  caption: [Prior art (sources: Estream @estream, WaterLily @waterlily @waterlily_review, Blue Freedom
  @bluefreedom @bluefreedom_review, Aquair UW @aquair). The last column assumes an electrical power coefficient
  of 0.30 to 0.20; for Ovo it is inverted, like the products, from its rated electrical power, the #d.gl_td15~W
  DC estimate (Ovo's own model reaches that power at #d.design_flow~m/s, an electrical power coefficient of
  #d.gl_td15_cpel).],
) <tab-prior>

Read with care: the data quality varies (some values come from reviews, and the WaterLily rotor size is
bounded by its housing). The comparison supports three statements. First, per swept area Ovo is in the class
of the Enomad Estream, the closest product: both have about 0.2 m rotors and folding blades; the Estream's
rating speed is unpublished. Second, Ovo's distinctive point is stowage: a closed #d.hull_vol_l~L egg (a
#d.hull_cyl_l~L bounding cylinder) against the Estream's #d.estream_len × #d.estream_d~mm tube
(#d.estream_cyl_l~L). Third, its predicted rotor performance is plausible but not exceptional: a full-scale
1.58 m river turbine measured $C_P$ #d.chica_cp @chica2025, against the #d.cp_peak15#est predicted here for a
far smaller, printed rotor at lower Reynolds number. Portable hydrokinetic prototypes are also studied in the
literature @riglin2016. What Ovo adds is an open, fully parametric and documented design.

= Limitations and v2 roadmap <sec-limits>

== Limitations

*Models, not measurements.* Every power, load and mass figure marked #est comes from a model:
- BEM is a quasi-steady annulus model: no hull blockage or hull wake, no yaw, shear or turbulence, rigid
  blades, coning by the independence principle, and no in-plane lead of the stacking axis.
- The polars are an XFOIL-trained network in a transitional Reynolds range (#d.re_all_lo to #d.re_all_hi over
  all maximum-power points); transition modelling alone moves $C_P$ by about 6 %, and printed roughness is
  the main unknown.
- The generator figures rest on datasheet values and an assumed iron and friction drag for the GL40; rectifier
  switching and cable losses are not modelled.
- Cavitation inception is taken at $sigma_v = -C_(p,"min")$ with no nucleation delay and no tip-vortex
  cavitation.
- Structure is checked by bending stresses at the lug and the fin root only: no stress concentration, fatigue
  or printed-strength data were sourced.

*Design margins that are thin or open.*
- Bearing span to overhang is #(d.span_ratio)\; the explicit checks pass, but wobble must be measured.
- Worst-case O-ring squeeze is #d.or_sq_lo~%, below the #d.tab_sq0~% table minimum.
- The v2 hull is nearly neutral (#d.buoy_v2~N#est net buoyancy for #d.mass_v2~g#est): reaching the +#d.v2_target~N
  target needs #d.v2_shed_g~g less mass or about #d.v2_foam_cc~cm³ of foam float (four-perimeter shells alone
  save #d.v2_4per_g~g#est).
- The fin couple ratio spreads over #d.fin_range0 to #d.fin_range1#est with the rail play, and residual roll
  needs a bridle or swivel.
- Detent forces are modelled; the stowed hold might resist the flow opening the blades, so flow-driven
  deployment is a design intent until a water test shows it.
- Towing above about #d.conflict_flow~m/s is unresolved: the cavitation cap needs more torque than the GL40's
  continuous rating there (#d.ll_x25_cap × the rated copper loss at 2.5 m/s#est), and the rated-current point
  with its #d.dump_r~Ω dump resistor cavitates. Options for v2: a larger generator, accepting a short-time
  overload inside the #d.gl_ipeak_r~A peak rating, a mechanical brake, or limiting towing to about
  #d.conflict_flow~m/s.
- The cost has #d.bom_asbuilt_unpriced unpriced lines; the hardware study's bill predates the CAD.

== Roadmap

+ *Print and measure the display build.* Coupons first, then the shells, rotor and drivetrain. Measure the
  fits, the detent forces and the bayonet click; check that the blades fold and deploy by hand, the rotor
  spins freely and wobble at the blade tips stays near the predicted #d.tip_runout~mm.
+ *Bench-test the generator.* Line resistance and inductance, $K_V$ from the open-circuit voltage, no-load
  loss from a spin-down; re-run the coupled study with the measured values.
+ *Electronics (v2).* A three-phase rectifier and a buck-boost stage to 5 V USB with an input rated for #d.input_v~V;
  maximum-power tracking along the load line, with the cavitation cap as a minimum shaft torque above
  #d.cap_rpm0~rpm, a current limit at the GL40 rating, and a dump load instead of a phase short. The cap and
  the current limit conflict above about #d.conflict_flow~m/s, so the towing strategy has to be settled first
  (@sec-limits).
+ *Water test.* Tow tests from 1.0 to 2.5 m/s logging flow speed, rotor speed, DC voltage and current to measure
  the electrical power coefficient; observe blade deployment, roll angle with and without a bridle, and
  cavitation near the surface.
+ *Design updates from the data.* Float or thinner shells for v2 buoyancy, crush ribs on the fin tongues,
  and a re-run of the blade design with measured polars or a corrected transition model.

= Conclusion

Ovo shows in CAD that a rotor #box[#d.rotor_ratio ×] the diameter of its body can fold flush into a printable egg
when the rotor is placed downstream and its blades fold forward, so that the flow is meant to deploy them and
the thrust to lock them; whether the flow does deploy the blades past the stowed detent is the first thing the
water test must show. The design is quantified rather than sketched: the cost of every printability constraint
is known in power coefficient, the generator choice and its operating law follow from a coupled loss model with
a cavitation limit (which leaves the towing range above about #d.conflict_flow~m/s as an open issue), and the
geometric claims about the #(d.n_parts_total)-part assembly are checked by the build itself on fine meshes of the
CAD solids. Its predictions are still predictions; the next step is to print it, bench-test the generator and
put it in the water.

#v(0.8em)
#block(width: 100%, inset: (top: 6pt), stroke: (top: 0.4pt + rule-c))[
  #set text(8.5pt, fill: ink2)
  #set par(justify: false)
  *Colophon.* © 2026 #d.author. This report, the CAD models and their exports are licensed under
  #d.licence (the code under MIT). Built #d.date by `scripts/report.py` from repository commit
  #raw(d.rev.split(" ").first())#if d.rev.contains("local") [ plus local changes] and CAD build #raw(d.cad_key).
]

#v(0.6em)
#set text(9pt)
#bibliography("refs.bib", title: [References], style: "ieee")
