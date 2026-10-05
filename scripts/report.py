"""Build the Ovo technical report: docs/report/ovo-report.pdf.

Every number in the report comes from the design-study outputs, the CAD checks and the
parameters, never from memory:

  data/blade_design.json, data/hardware.json, data/bom.json, data/layout.json
  out/cad/checks/*.json (assembly, blade, rotor, shells, drivetrain, fins, coupon)
  out/web/manifest.json, ovo/params.py, ovo/profile.py

The script
  1. extracts and formats those numbers into docs/report/data.typ (generated, do not edit),
  2. writes the bibliography docs/report/refs.bib (curated references plus the sourced
     catalogue and prior-art URLs of data/hardware.json),
  3. draws the report's own figures into docs/report/figures/ (hull profile, concept sketch),
  4. copies the renders as compressed JPEGs into docs/report/img/ (docs/media/* when present,
     else out/renders/*),
  5. compiles docs/report/ovo-report.typ with the typst Python package (embedded fonts only,
     so the PDF does not depend on the machine's fonts).

  .venv/bin/python scripts/report.py                  # data, figures, images, PDF
  .venv/bin/python scripts/report.py --png DIR        # also every page as PNG (DIR/page-NN.png)
  .venv/bin/python scripts/report.py --render-fold    # first re-render the fold-sequence strip
                                                      # through scripts/render.py (takes the heavy lock)

Requires the docs dependency group (typst). Run the CAD build first (scripts/build.py): the
checks and renders under out/ are inputs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from ovo import params as P
from ovo import profile

REPORT = REPO / "docs" / "report"
FIG = REPORT / "figures"
IMG = REPORT / "img"
SRC_TYP = REPORT / "ovo-report.typ"
PDF = REPORT / "ovo-report.pdf"
FOLD_FRAMES = REPO / "out" / "work" / "report" / "fold"

AUTHOR = "Luca Ramírez Olivo"
LICENCE = "CC BY 4.0"
AFFILIATION = "Computer Engineering (final year), Universidad Europea de Madrid"
EM_DASH = "\u2014"


# --------------------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------------------
def load(rel: str) -> dict:
    return json.loads((REPO / rel).read_text())


def f(v, nd: int = 1) -> str:
    """Fixed decimals; integers >= 10000 get a thousands comma."""
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    if nd == 0:
        iv = round(v)
        return f"{iv:,}" if abs(iv) >= 10000 else str(iv)
    return f"{v:.{nd}f}"


def pct(v, nd: int = 1, sign: bool = False) -> str:
    s = f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"
    return s


def g(v, sig: int = 3) -> str:
    """Compact general format (for mixed magnitudes)."""
    return f"{v:.{sig}g}"


def git_rev() -> str:
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO,
                               capture_output=True, text=True, check=False).stdout.strip()
        return rev + (" + local changes" if dirty else "")
    except Exception:  # noqa: BLE001
        return "unknown"


# --------------------------------------------------------------------------------------
# typst literal writer
# --------------------------------------------------------------------------------------
_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_-]*$")


def tlit(v, ind: int = 0) -> str:
    pad = "  " * ind
    if v is None:
        return "none"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return repr(v)
    if isinstance(v, str):
        s = v.replace(EM_DASH, ", ").replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
        return f'"{s}"'
    if isinstance(v, dict):
        if not v:
            return "(:)"
        items = []
        for k, x in v.items():
            key = k if _IDENT.match(k) else tlit(k)
            items.append(f"{pad}  {key}: {tlit(x, ind + 1)}")
        return "(\n" + ",\n".join(items) + f",\n{pad})"
    if isinstance(v, (list, tuple)):
        if not v:
            return "()"
        flat = all(not isinstance(x, (list, tuple, dict)) for x in v)
        if flat:
            return "(" + ", ".join(tlit(x) for x in v) + ",)"
        return "(\n" + ",\n".join(f"{pad}  {tlit(x, ind + 1)}" for x in v) + f",\n{pad})"
    raise TypeError(type(v))


# --------------------------------------------------------------------------------------
# data extraction
# --------------------------------------------------------------------------------------
def by_flow(rows: list[dict], flow: float, key: str = "flow_m_s") -> dict:
    return next(r for r in rows if abs(r[key] - flow) < 1e-6)


def hull_numbers() -> dict:
    x = np.linspace(0.0, P.HULL_L, 20001)
    r = profile.radius(x)
    vol = np.trapezoid(math.pi * r**2, x) / 1000.0                     # cm3
    drdx = np.gradient(r, x)
    wet = np.trapezoid(2 * math.pi * r * np.sqrt(1 + drdx**2), x) / 100.0  # cm2
    frontal = math.pi * (P.HULL_D / 2) ** 2 / 100.0                     # cm2
    x_tip = P.SPINNER_CAP["x_tip_built"]
    cyl = math.pi * (P.HULL_D / 20) ** 2 * (x_tip / 10) / 1000.0        # L, bounding cylinder
    re_l = P.DESIGN_FLOW_V * P.HULL_L * 1e-3 / P.WATER_NU
    return dict(vol_cm3=vol, wet_cm2=wet, frontal_cm2=frontal, cyl_L=cyl, re_l=re_l, x_tip=x_tip,
                ld=P.HULL_L / P.HULL_D, ld_built=x_tip / P.HULL_D,
                r_static_tail=float(profile.radius(P.STATIC_TAIL_X)),
                r_hinge=float(profile.radius(P.HINGE_X)))


def clean_check(v: str) -> str:
    """Typeset-friendly check strings: units with superscripts, comparison signs, no internal ticket ids."""
    v = re.sub(r"owner call T\d+:\s*", "", v)
    v = v.replace(">=", "≥").replace("<=", "≤").replace("+-", "±").replace("mm3", "mm³").replace("mm2", "mm²")
    return v.replace(" deg", "°")


def build_data() -> tuple[dict, dict]:
    bd = load("data/blade_design.json")
    hw = load("data/hardware.json")
    bom = load("data/bom.json")
    lay = load("data/layout.json")
    asm = load("out/cad/checks/assembly.json")
    fins = load("out/cad/checks/fins.json")
    blade_ck = load("out/cad/checks/blade.json")
    mods = {m: load(f"out/cad/checks/{m}.json") for m in ("shells", "rotor", "drivetrain", "fins", "coupon")}
    manifest = load("out/web/manifest.json")

    d: dict = {}
    t: dict = {}

    # ---- meta --------------------------------------------------------------------------
    d["date"] = dt.datetime.now().astimezone().date().isoformat()
    d["rev"] = git_rev()
    d["author"] = AUTHOR
    d["affiliation"] = AFFILIATION

    # ---- hull ---------------------------------------------------------------------------
    hn = hull_numbers()
    d.update(
        hull_d=f(P.HULL_D, 0), nose_a=f(P.HULL_NOSE_A, 0), nose_n=f(P.HULL_NOSE_N, 0),
        mid_b=f(P.HULL_MID_B, 0), tail_c=f(P.HULL_TAIL_C, 0), tail_theta=f(P.HULL_TAIL_THETA, 0),
        hull_l=f(P.HULL_L, 0), x_tip=f(hn["x_tip"], 2), ld=f(hn["ld"], 2), ld_built=f(hn["ld_built"], 2),
        wall=f(P.HULL_WALL, 1), equator=f(P.EQUATOR_X, 1), static_tail=f(P.STATIC_TAIL_X, 0),
        hub_gap=f(P.HUB_GAP, 1), spinner_tip_r=f(P.SPINNER_TIP_R, 0),
        hull_vol=f(hn["vol_cm3"], 0), hull_vol_l=f(hn["vol_cm3"] / 1000, 2), hull_wet=f(hn["wet_cm2"], 0),
        hull_frontal=f(hn["frontal_cm2"], 1), hull_cyl_l=f(hn["cyl_L"], 2),
        hull_re=f(hn["re_l"] / 1e5, 1), r_static_tail=f(hn["r_static_tail"], 1),
    )

    # ---- rotor and BEM ---------------------------------------------------------------------
    rot = bd["rotor"]
    laws = bd["laws"]
    sec = bd["section"]
    mppt = bd["performance"]["mppt"]
    runaway = bd["performance"]["runaway"]
    env = bd["performance"]["envelope_up_to_runaway"]
    p15 = by_flow(mppt, 1.5)
    p25 = by_flow(mppt, 2.5)
    p10 = by_flow(mppt, 1.0)
    run25 = by_flow(runaway, 2.5)
    run15 = by_flow(runaway, 1.5)
    env25 = by_flow(env, 2.5)
    R = rot["r_tip_mm"] * 1e-3
    area = math.pi * R**2
    p_avail15 = 0.5 * P.WATER_RHO * area * 1.5**3
    d.update(
        n_blades=str(rot["n_blades"]), r_tip=f(rot["r_tip_mm"], 0), rotor_d=f(rot["swept_diameter_mm"], 0),
        r_aero_root=f(rot["r_aero_root_mm"], 0), hinge_x=f(rot["hinge_x_mm"], 0), hinge_r=f(rot["hinge_r_mm"], 0),
        cone=f(rot["cone_deg"], 0), design_flow=f(rot["design_flow_m_s"], 1), design_tsr=f(rot["design_tsr"], 0),
        tsr_peak=f(rot["tsr_at_cp_peak_1p5_m_s"], 2),
        chord_root=f(laws["chord_mm"]["root"], 1), chord_tip=f(laws["chord_mm"]["tip"], 1),
        twist_root=f(laws["twist_deg"]["root"], 1), twist_tip=f(laws["twist_deg"]["tip"], 1),
        tc_root=f(sec["thickness_ratio_root"] * 100, 0), tc_tip=f(sec["thickness_ratio_tip"] * 100, 0),
        te=f(sec["te_thickness_mm"], 1), pitch_axis=f(sec["pitch_axis_chord_fraction"], 2),
        te_c_root=f(bd["stations"]["te_thickness_ratio"][0] * 100, 1),
        te_c_tip=f(bd["stations"]["te_thickness_ratio"][-1] * 100, 1),
        n_el=str(bd["model"]["elements"]), polars=bd["model"]["polars"],
        cp_peak15=f(p15["cp"], 3), p_shaft15=f(p15["power_W"], 1), rpm15=f(p15["rpm"], 0),
        q15=f(p15["torque_Nm"], 3), t15=f(p15["thrust_N"], 1), flap15=f(p15["flap_moment_hinge_per_blade_Nm"], 3),
        cp_peak10=f(p10["cp"], 3), p_shaft10=f(p10["power_W"], 1),
        cp_peak25=f(p25["cp"], 3), p_shaft25=f(p25["power_W"], 1), rpm25=f(p25["rpm"], 0),
        t25=f(p25["thrust_N"], 1),
        run_rpm25=f(run25["rpm"], 0), run_tsr25=f(run25["tsr"], 2), run_t25=f(run25["thrust_N"], 1),
        run_rpm15=f(run15["rpm"], 0), run_tsr15=f(run15["tsr"], 2),
        tmax25=f(env25["max_thrust_N"], 1), flapmax25=f(env25["max_flap_moment_hinge_per_blade_Nm"], 3),
        p_avail15=f(p_avail15, 1), area_m2=f(area, 4),
    )
    t["perf"] = [[f(r["flow_m_s"], 1), f(r["tsr"], 2), f(r["rpm"], 0), f(r["cp"], 3), f(r["power_W"], 1),
                  f(r["torque_Nm"], 3), f(r["thrust_N"], 1), f(r["flap_moment_hinge_per_blade_Nm"], 3)]
                 for r in mppt]
    t["limits"] = [[f(r["flow_m_s"], 1), f(r["tsr"], 2), f(r["rpm"], 0), f(r["thrust_N"], 1),
                    f(e["max_thrust_N"], 1), f(e["max_flap_moment_hinge_per_blade_Nm"], 3)]
                   for r, e in zip(runaway, env)]
    best_tsr = max(bd["tsr_study"], key=lambda r: r["cp_constrained"])
    t["tsr_study"] = [[f(r["design_tsr"], 1), f(r["cp_bem_optimum"], 4), f(r["cp_schmitz"], 4),
                       f(r["cp_constrained"], 4), f"{f(r['cp_constrained_peak'], 4)} ({f(r['tsr_constrained_peak'], 2)})",
                       f"{f(r['constrained_x'][0], 1)} / {f(r['constrained_x'][1], 1)}",
                       f(r["rpm_at_1p5_m_s"], 0), f(r["torque_at_1p5_m_s_Nm"], 3), f(r["copper_loss_index"], 2),
                       r is best_tsr] for r in bd["tsr_study"]]
    chain = {c["key"]: c for c in bd["constraint_chain"]}
    order = ["optimum", "te", "linear", "chord", "twist", "final"]
    labels = {"optimum": "BEM optimum (sharp TE, no cone)", "te": f"+ blunt TE {P.MIN_TE:g} mm",
              "linear": "+ linear chord and twist laws",
              "chord": f"+ chord within [{P.MIN_CHORD:g}, {P.BLADE_MAX_CHORD:g}] mm",
              "twist": f"+ root twist at most {P.MAX_ROOT_TWIST:g}°",
              "final": f"+ {rot['cone_deg']:g}° aft cone (final blade)"}
    cp0 = chain["optimum"]["cp_design_tsr"]
    rows, prev = [], None
    for k in order:
        c = chain[k]
        step = "" if prev is None else pct(100 * (c["cp_design_tsr"] - prev) / prev, 1, True) + " %"
        tot = "" if prev is None else pct(100 * (c["cp_design_tsr"] - cp0) / cp0, 1, True) + " %"
        rows.append([labels[k], f(c["cp_design_tsr"], 4), step, tot,
                     f"{f(c['chord_root_mm'], 1)} / {f(c['chord_tip_mm'], 1)}",
                     f"{f(c['twist_root_deg'], 1)} / {f(c['twist_tip_deg'], 1)}"])
        prev = c["cp_design_tsr"]
    t["chain"] = rows
    d["chain_total_pct"] = pct(100 * (chain["final"]["cp_design_tsr"] - cp0) / cp0, 1)
    d["cp_opt"] = f(cp0, 4)
    d["cp_final_dtsr"] = f(chain["final"]["cp_design_tsr"], 4)
    d["chain_te_pct"] = pct(100 * (chain["te"]["cp_design_tsr"] - cp0) / cp0, 1)
    d["chain_lin_pct"] = pct(100 * (chain["linear"]["cp_design_tsr"] - chain["te"]["cp_design_tsr"])
                             / chain["te"]["cp_design_tsr"], 1)
    d["cp_schmitz"] = f(chain["schmitz"]["cp_design_tsr"], 4)
    for k in ("chain_total_pct", "chain_te_pct", "chain_lin_pct"):
        d[k.replace("_pct", "_abs")] = d[k].lstrip("-")
    d["fit_0709_pct"] = pct(100 * (chain["fit_0p7_0p9"]["cp_design_tsr"] - chain["te"]["cp_design_tsr"])
                            / chain["te"]["cp_design_tsr"], 1)
    d["fit_0709_abs"] = d["fit_0709_pct"].lstrip("-")
    # first column as typst markup (subscripts), evaluated by btab(markup: true)
    t["sens"] = [[s["case"].replace("n_crit", '$n_"crit"$').replace(" deg", "°"), f(s["cp_peak"], 4),
                  f(s["tsr_peak"], 2), pct(s["delta_cp_pct"], 1, True) + " %"] for s in bd["sensitivity"]]
    sens = {s["case"]: s for s in bd["sensitivity"]}
    d["sens_ncrit9"] = pct(sens["n_crit 9 (clean inflow, smooth blade)"]["delta_cp_pct"], 1)
    d["sens_forced"] = pct(sens["forced transition x/c 0.05 (rough print)"]["delta_cp_pct"], 1)
    d["sens_hub"] = pct(sens["no hub loss"]["delta_cp_pct"], 1)
    d["sens_ncrit9_abs"], d["sens_forced_abs"] = d["sens_ncrit9"].lstrip("-"), d["sens_forced"].lstrip("-")
    ck = bd["checks"]
    ca = bd["convergence_audit"]
    d.update(
        glauert_bem=f(bd["glauert_verification"]["cp_bem"], 5), glauert_th=f(bd["glauert_verification"]["cp_theory"], 5),
        cp_max_any=f(ck["max_cp_any_sweep"], 4),
        a_mid_lo=f(ck["a_midspan_range"][0], 2), a_mid_hi=f(ck["a_midspan_range"][1], 2),
        n_buhl=str(ck["n_elements_buhl_design"]),
        re_lo=f(ck["re_design_range"][0], 0), re_hi=f(ck["re_design_range"][1], 0),
        re_all_lo=f(ck["re_range_all_mppt"][0], 0), re_all_hi=f(ck["re_range_all_mppt"][1], 0),
        alpha_lo=f(ck["alpha_design_range_deg"][0], 1), alpha_hi=f(ck["alpha_design_range_deg"][1], 1),
        stall_margin=f(ck["stall_margin_min_deg"], 1), nf_conf=f(ck["nf_confidence_min_design"], 3),
        conv_ops=f(ca["operating_points"], 0), conv_solves=f(ca["element_solves"], 0),
        conv_unconv=str(len(ca["unconverged"])), conv_fallback=str(ca["fallback_elements"]),
        conv_root_share=f(100 * ca["root_elements"]["max_torque_share_power_producing_sweeps"], 2),
    )
    # cavitation and load line
    cav = {r["flow_m_s"]: r for r in bd["cavitation_vs_tsr"]["rows"]}
    gce = bd["generator_coupling_estimate"]
    ll = gce["load_line_by_flow"]
    gl40_ll = next(gn for gn in ll["generators"] if gn["short"] == "GL40")
    pts = {round(pp["flow_m_s"], 1): pp for pp in gl40_ll["points"]}
    lims = {round(lm["flow_m_s"], 1): lm for lm in ll["limits"]}
    rows = []
    for v in (1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5):
        pp, lm = pts[v], lims[v]
        un, cp_ = pp["unconstrained"], pp["capped"]
        rows.append([f(v, 1), f(lm["tsr_cp_peak"], 2), f"{f(un['tsr'], 2)} ({f(un['rpm'], 0)})",
                     f(un["net_power_W"], 1), pct(un["cav_margin"], 2, True),
                     f(lm["tsr_cavitation_free_to"], 2) if lm["tsr_cavitation_free_to"] else "none",
                     f"{f(cp_['tsr'], 2)} ({f(cp_['rpm'], 0)})", f(cp_["net_power_W"], 1),
                     ("0" if abs(pp["cap_cost_pct"]) < 0.05 else pct(pp["cap_cost_pct"], 1) + " %")])
    t["load_line"] = rows
    cap = ll["cap_torque_curve"]
    d.update(
        cav_margin_cap=f(ll["cav_margin_cap"], 2), cav_depth=f(ll["depth_m"], 1),
        cap_rpm0=f(cap["rpm"][0], 0), cap_rpm1=f(cap["rpm"][-1], 0),
        cap_q0=f(cap["min_shaft_torque_Nm"][0], 3), cap_q1=f(cap["min_shaft_torque_Nm"][-1], 3),
        cap_checks=str(cap["check_points"]),
        cap_bad=str(cap["check_allowed_but_below_cap"] + cap["check_forbidden_but_above_cap"]),
        cav_free15=f(cav[1.5]["tsr_free_to"], 2), cav_free20=f(cav[2.0]["tsr_free_to"], 2),
        cav_free25=f(cav[2.5]["tsr_free_to"], 2), cav_peak25=f(cav[2.5]["margin_at_cp_peak"], 2),
        cap_cost25=f(abs(pts[2.5]["cap_cost_pct"]), 1), ll_net25_un=f(pts[2.5]["unconstrained"]["net_power_W"], 1),
        ll_net25_cap=f(pts[2.5]["capped"]["net_power_W"], 1), ll_rpm25_cap=f(pts[2.5]["capped"]["rpm"], 0),
        ll_tsr25_cap=f(pts[2.5]["capped"]["tsr"], 2), ceil_cost25=f(abs(pts[2.5]["ceiling_cost_pct"]), 1),
        ll_net10=f(pts[1.0]["capped"]["net_power_W"], 1), ll_net15=f(pts[1.5]["capped"]["net_power_W"], 1),
        ll_net20=f(pts[2.0]["capped"]["net_power_W"], 1),
        cap_first_flow=f(min(v for v, pp in pts.items() if pp["unconstrained"]["cav_margin"] < 0), 1),
    )
    # the cavitation cap against the GL40's continuous rating (A1): from which flow the capped load line needs
    # more than the rated current, and by how much at 2.5 m/s (copper loss over the rated copper loss)
    gl_rating = next(c for c in load("data/hardware.json")["generator"]["candidates"]
                     if c["name"] == "CubeMars GL40 KV70")["rating"]["value"]
    over = sorted(v for v, pp in pts.items() if pp["capped"]["current_A"] > gl_rating["I_rated_A"])
    d.update(
        conflict_flow=f(over[0], 1), conflict_i=f(pts[over[0]]["capped"]["current_A"], 2),
        conflict_x=f(pts[over[0]]["capped"]["copper_loss_W"] / gl_rating["P_cu_rated_W"], 2),
        ll_i25_cap=f(pts[2.5]["capped"]["current_A"], 1),
        ll_x25_cap=f(pts[2.5]["capped"]["copper_loss_W"] / gl_rating["P_cu_rated_W"], 2),
        gl_ipeak_r=f(gl_rating["I_peak_A"], 1),
    )
    gl40 = next(r for r in gce["rows"] if r["short"] == "GL40" and r["blade"] == "final")
    bn, atp = gl40["best_net"], gl40["at_cp_peak"]
    d.update(
        gen_best_tsr=f(bn["tsr"], 2), gen_best_rpm=f(bn["rpm"], 0), gen_best_q=f(bn["torque_Nm"], 3),
        gen_best_net=f(bn["net_power_W"], 1), gen_best_shaft=f(bn["shaft_power_W"], 1),
        gen_best_cu=f(bn["copper_loss_W"], 1), gen_best_drag=f(bn["drag_loss_W"], 1),
        gen_best_i=f(bn["current_A"], 2), gen_peak_net=f(atp["net_power_W"], 1), gen_peak_cu=f(atp["copper_loss_W"], 1),
    )
    se = bd["structure_estimates"]
    lac = {r["flow_m_s"]: r for r in bd["layout_axis_check"]["rows"]}
    d.update(
        lug_stress25=f(next(c for c in se["cases"] if c["case"] == "mppt_2p5_m_s")["flap_stress_lug_thin_depth_MPa"], 1),
        lug_stress15=f(next(c for c in se["cases"] if c["case"] == "mppt_1p5_m_s")["flap_stress_lug_thin_depth_MPa"], 1),
        cp_layout25=f(lac[2.5]["cp_peak"], 4), cp_nominal25=f(lac[2.5]["cp_peak_nominal_axis"], 4),
        flap_layout25=f(lac[2.5]["max_flap_moment_hinge_Nm"], 3), flap_nominal25=f(lac[2.5]["max_flap_moment_hinge_nominal_axis_Nm"], 3),
    )

    # ---- layout: rotor half --------------------------------------------------------------
    ls = lay["summary"]
    lr = lay["rotor"]
    psum = lr["pocket_summary"]
    pick = lr["selection"]["rule_pick"]
    hub_lc = lr["hub"]["load_cases"]["2.5"]
    det = lr["detents"]
    d.update(
        hub_front=f(ls["hub_front_x_mm"], 1), shaft_end=f(ls["shaft_end_x_mm"], 0), fold_deg=f(ls["fold_deg"], 0),
        pocket_x0=f(ls["pocket_x_mm"][0], 2), pocket_x1=f(ls["pocket_x_mm"][1], 0),
        pocket_eq=f(ls["pocket_depth_at_equator_mm"], 2), pocket_aero=f(ls["pocket_depth_aero_max_mm"], 2),
        pocket_max=f(ls["pocket_depth_max_mm"], 2), pocket_land=f(psum["land_min_mm"], 1),
        pocket_width=f(lr["pocket_width_mm"], 1),
        stow_flush=f(lr["stowed_max_over_hull"]["aero_mm"], 3),
        lead_root=f(P.BLADE_AXIS["lead_mm"][0], 1), lead_tip=f(P.BLADE_AXIS["lead_mm"][1], 1),
        axis_cone_root=f(bd["layout_axis_check"]["axis_cone_local_deg_root_tip"][0], 1),
        axis_cone_tip=f(bd["layout_axis_check"]["axis_cone_local_deg_root_tip"][1], 1),
        pin_cover=f(pick["pin_cover_in_hub_core_mm"], 2), heel_inside=f(pick["heel_stowed_inside_spinner_mm"], 2),
        heel_face=f(pick["heel_flat_face_mm"], 2), stow_tip_x=f(pick["stowed_tip_x_mm"], 2),
        hub_wall=f(pick["hub_wall_over_bore_mm"], 2),
        stop_force25=f(hub_lc["stop_force_N"], 0), stop_p25=f(hub_lc["stop_contact_pressure_MPa"], 2),
        pin_shear25=f(hub_lc["pin_double_shear_MPa"], 1), ear_bearing25=f(hub_lc["ear_bearing_MPa"], 1),
        det_force=f(det["deployed"]["force_N"], 2), det_gap=f(det["deployed"]["gap_mm"], 1),
        det_k=f(det["calibration"]["k_design"], 2),
    )

    # ---- layout: internal half -----------------------------------------------------------
    li = lay["internal"]
    br = li["drivetrain"]["bearings"]
    roll = li["roll"]
    fn = li["fins"]
    mb = li["mass_buoyancy"]
    loads = li["nose"]["loads"]
    cs = {c["check"].split(" ")[0]: c for c in lay["cross_study_checks"]}
    fin_cs = next(c for c in lay["cross_study_checks"] if c["check"].startswith("fin couple"))
    oring_cs = next(c for c in lay["cross_study_checks"] if c["check"].startswith("equator O-ring"))
    del cs
    d.update(
        bay_d=f(ls["bay_D_mm"], 1), bay_l=f(ls["bay_L_mm"], 0), bay_x0=f(ls["bay_x_mm"][0], 1), bay_x1=f(ls["bay_x_mm"][1], 1),
        b608_x0=f(ls["bearing_608_x_mm"][0], 0), b608_x1=f(ls["bearing_608_x_mm"][1], 0),
        b688_x0=f(ls["bearing_688_x_mm"][0], 0), b688_x1=f(ls["bearing_688_x_mm"][1], 0),
        seal_x0=f(ls["lip_seal_x_mm"][0], 1), seal_x1=f(ls["lip_seal_x_mm"][1], 1),
        span=f(ls["bearing_span_mm"], 1), overhang=f(ls["rotor_overhang_mm"], 1), span_ratio=f(ls["span_over_overhang"], 2),
        tip_runout=f(ls["tip_runout_mm"], 2), s0_front=f(br["static_safety_front_s0"], 1), s0_rear=f(br["static_safety_rear_s0"], 1),
        shaft_defl=f(br["shaft_deflection_at_load_plane_mm"], 3), thrust_tow=f(br["thrust_tow_N"], 1),
        interface_pass=str(ls["interface_checks_passed"]), interface_n=str(ls["interface_checks"]),
        tow_design=f(ls["tow_design_load_N"], 0), tow_working=f(loads["working_load_N"], 1),
        hull_drag=f(loads["hull_drag_est_N"], 2), hull_cd=f(loads["hull_cd_est"], 2),
        mass_display=f(mb["display_PLA"]["mass_g"], 0), buoy_display=f(mb["display_PLA"]["net_buoyancy_N"], 2),
        mass_v2=f(mb["v2_PETG"]["mass_g"], 0), buoy_v2=f(mb["v2_PETG"]["net_buoyancy_N"], 3),
        disp_cc=f(mb["v2_PETG"]["displacement_cc"], 0),
        oring_dev=f(oring_cs["max_dimension_deviation_mm"], 1),
    )
    q = roll["torque"]
    bal = roll["balance"]
    d.update(
        roll_q10=f(q["1_m_s"]["Q_Nm"], 4), roll_q15=f(q["1.5_m_s"]["Q_Nm"], 4), roll_q25=f(q["2.5_m_s"]["Q_Nm"], 4),
        roll_pend=f(roll["pendulum_righting_max_Nm"], 4), roll_rate15=f(bal["1.5_m_s"]["free_roll_rate_rpm"], 1),
        roll_ex15=f(bal["1.5_m_s"]["excess_pct"], 1), roll_ex10=f(bal["1_m_s"]["excess_pct"], 1),
        roll_ex25=f(bal["2.5_m_s"]["excess_pct"], 1), roll_rate10=f(abs(bal["1_m_s"]["free_roll_rate_rpm"]), 1),
        roll_rate25=f(bal["2.5_m_s"]["free_roll_rate_rpm"], 1), bridle_heel15=f(bal["1.5_m_s"]["bridle_heel_deg"], 1),
        freewheel_roll=f(roll["freewheel_roll_rate_1p5_m_s_rpm"], 0),
        fin_chord=f(fn["chord_mm"], 0), fin_span=f(fn["span_mm"], 0), fin_ar=f(fn["AR"], 2),
        fin_area=f(fn["area_each_cm2"], 1), fin_inc=f(fn["incidence_deg"], 1),
        fin_theta0=f(fn["theta_deg"][0], 1), fin_theta1=f(fn["theta_deg"][1], 1),
        fin_rcp=f(fn["r_cp_mm"], 1), fin_tip2tip=f(fn["tip_to_tip_mm"], 1),
        fin_x0=f(fn["x_le_mm"], 1), fin_x1=f(fn["x_te_mm"], 1),
        fin_torque=f(fin_cs["fin_torque_Nm"], 4), fin_ratio=f(fin_cs["ratio"], 3),
        fin_lim0=f(fin_cs["limit"][0], 2), fin_lim1=f(fin_cs["limit"][1], 3),
        fin_ratio10=f(fin_cs["off_design"]["1.0_m_s"]["ratio"], 2), fin_ratio25=f(fin_cs["off_design"]["2.5_m_s"]["ratio"], 2),
        fin_cl15=f(fn["margins"]["1.5_m_s"]["CL"], 2), fin_clmax15=f(fn["margins"]["1.5_m_s"]["CLmax_3D"], 2),
        fin_re15=f(fn["margins"]["1.5_m_s"]["Re"], 0), fin_drag=f(fn["drag_tow_N"], 2),
        fin_trim=f(100 * fn["trim"]["ratio_change_per_deg_at_design_flow"], 1),
        fin_clear_pockets=f(fn["clearance"]["stowed_pockets"]["min_clear_mm"], 2),
        fin_slot_fold_clear=f(fn["clearance"]["fold_at_index"]["pad_clear_min_mm"], 1),
        fin_fin_fold_clear=f(fn["clearance"]["fold_at_index"]["fin_clear_min_mm"], 1),
    )
    cr = fins["couple_ratio_as_built"]
    yp = fins["yaw_play"]
    st = fins["structure"]
    d.update(
        fin_ratio_built=f(cr["ratio"], 3), fin_range0=f(cr["range"][0], 2), fin_range1=f(cr["range"][1], 2),
        fin_yaw_p=f(yp["free_plus_deg"], 2), fin_yaw_m=f(yp["free_minus_deg"], 2),
        fin_root_stress=f(st["sigma_MPa"], 1), fin_root_util=f(st["utilisation"], 2), fin_interlayer=f(st["interlayer_MPa"], 0),
    )

    # ---- hardware: generator --------------------------------------------------------------
    cands = hw["generator"]["candidates"]
    fits = {c["name"]: c for c in li["generator_bay"]["candidate_fits"]}
    rows = []
    for c in cands:
        studies = c["study"]
        kv = sorted({s["KV_rpm_per_V"] for s in studies})
        rr = sorted({s["R_ll_ohm"] for s in studies})
        km = [s["Km_Nm_per_sqrtW"] for s in studies]
        pb = [s["coupled_best"]["bridge_1.5_m_s"]["P_el_dc_W"] for s in studies]
        pa = [s["coupled_best"]["active_1.5_m_s"]["P_el_dc_W"] for s in studies]
        D = c["D"]["value"] if isinstance(c.get("D"), dict) else None
        L = c["L"]["value"] if isinstance(c.get("L"), dict) else None

        def rng(vals, nd):
            lo, hi = min(vals), max(vals)
            return f(lo, nd) if abs(hi - lo) < 1e-9 else f"{f(lo, nd)}-{f(hi, nd)}"

        fit = fits.get(c["name"])
        price = c.get("price", {}).get("value", {})
        pstr = (f"USD {f(price['USD'], 2)}" if "USD" in price else f"EUR {f(price['EUR'], 2)}" if "EUR" in price else "n/a")
        name = re.sub(r"\s*\((outside envelope|frameless)\)", "", c["name"])
        rows.append([name, f"{f(D, 1) if D else 'n/a'} x {f(L, 1) if L else 'n/a'}",
                     "/".join(f(k, 1).rstrip("0").rstrip(".") for k in kv), "/".join(f(x, 2).rstrip("0").rstrip(".") for x in rr),
                     rng(km, 3), rng(pb, 2), rng(pa, 2),
                     ("yes" if fit and fit["fits_bay"] else "no") if fit else "n/a", pstr,
                     c["name"] == hw["generator"]["recommendation"]["choice"]])
    t["gen_candidates"] = rows
    gl = next(c for c in cands if c["name"] == "CubeMars GL40 KV70")
    gs = gl["study"][0]
    td = gs["coupled_best_time_domain_bridge"]
    act = gs["coupled_best"]["active_1.5_m_s"]
    rating = gl["rating"]["value"]
    d.update(
        gl_kv=f(gl["KV"]["value"], 1), gl_r=f(gl["R"]["value"], 1), gl_l=f(td["td_L_ll_H"] * 1e3, 1),
        gl_kt0=f(60 / (2 * math.pi * gl["KV"]["value"]), 4), gl_km=f(gs["Km_Nm_per_sqrtW"], 4),
        gl_mass=f(gl["mass"]["value"], 0), gl_price=f(gl["price"]["value"]["USD"], 2),
        gl_d=f(gl["D"]["value"], 1), gl_len=f(gl["L"]["value"], 1), gl_pp=gl["config"]["value"].split()[0],
        gl_trated=f(rating["T_rated_Nm"], 2), gl_irated=f(rating["I_rated_A"], 2), gl_ipeak=f(rating["I_peak_A"], 1),
        gl_pcu_rated=f(rating["P_cu_rated_W"], 2), gl_drag=f(gl["iron_friction_drag_Nm"]["value"], 4),
        gl_td15=f(td["td_bridge_1.5_m_s"]["P_el_dc_W"], 2), gl_td15_v=f(td["td_bridge_1.5_m_s"]["V_dc_V"], 2),
        gl_td15_i=f(td["td_bridge_1.5_m_s"]["I_A"], 3), gl_td15_rpm=f(td["td_bridge_1.5_m_s"]["rpm"], 0),
        gl_td15_tsr=f(td["td_bridge_1.5_m_s"]["tsr"], 2), gl_td15_cpel=f(td["td_bridge_1.5_m_s"]["Cp_el"], 3),
        gl_block15=f(gs["coupled_best"]["bridge_1.5_m_s"]["P_el_dc_W"], 2),
        gl_act15=f(act["P_el_dc_W"], 2), gl_act15_tsr=f(act["tsr"], 2),
        gl_td10=f(td["td_bridge_1.0_m_s"]["P_el_dc_W"], 2), gl_td10_v=f(td["td_bridge_1.0_m_s"]["V_dc_V"], 2),
        gl_td20=f(td["td_bridge_2.0_m_s"]["P_el_dc_W"], 2), gl_td25=f(td["td_bridge_2.5_m_s"]["P_el_dc_W"], 2),
        gl_td25_cap=f(td["td_bridge_capped_2.5_m_s"]["P_el_dc_W"], 2), gl_td25_cap_tsr=f(td["td_bridge_capped_2.5_m_s"]["tsr"], 2),
        gl_td25_tau=f(td["td_bridge_2.5_m_s"]["tau_em_Nm"], 4), gl_td25_pcu=f(td["td_bridge_2.5_m_s"]["P_cu_W"], 2),
        gl_voc25=f(gs["Voc_freewheel_runaway"]["2.5_m_s"]["Voc_ll_peak_V"], 1),
        gl_short25=f(gs["phase_short_brake_by_flow"]["2.5_m_s"]["P_cu_W"], 1),
        gl_short25_x=f(gs["phase_short_brake_by_flow"]["2.5_m_s"]["P_cu_over_rated"], 2),
        gl_short15=f(gs["phase_short_brake_by_flow"]["1.5_m_s"]["P_cu_W"], 1),
        gl_load_frac=f(100 * td["td_bridge_1.5_m_s"]["tau_em_Nm"] / rating["T_rated_Nm"], 0),
    )
    gm = next(c for c in cands if c["name"].startswith("iPower GM4108H"))
    gmv = [s["Voc_freewheel_runaway"]["2.5_m_s"]["Voc_ll_peak_V"] for s in gm["study"]]
    gmp = [s["coupled_best"]["bridge_1.5_m_s"]["P_el_dc_W"] for s in gm["study"]]
    gchk = hw["generator"]["gm4108_datasheet_check"]["value"]
    rat = [r_ for k_ in gchk["kv_derived"] for r_ in k_["ratio_to_datasheet"]]
    collapse = [c_["bridge_1p5"]["P_el_dc_W"] for c_ in gchk["coupled_if_Kt_from_load_point"]]
    d.update(gm_ratio_lo=f(min(rat), 1), gm_ratio_hi=f(max(rat), 1), gm_kv_lo=f(min(k_["KV"] for k_ in gchk["kv_derived"]), 0),
             gm_kv_hi=f(max(k_["KV"] for k_ in gchk["kv_derived"]), 0), gm_collapse_lo=f(min(collapse), 1),
             gm_collapse_hi=f(max(collapse), 1))
    d.update(gm_voc_lo=f(min(gmv), 1), gm_voc_hi=f(max(gmv), 1), gm_p_lo=f(min(gmp), 2), gm_p_hi=f(max(gmp), 2),
             gm_price=f(gm["price"]["value"]["USD"], 2))
    tp = hw["generator"]["recommendation"]["tow_protection"]["value"]["optional_dump_load"]
    d.update(dump_r=f(tp["R_dump_ohm"], 2), dump_p=f(tp["P_dump_W"], 2))

    # ---- hardware: bearings, seal, O-ring ------------------------------------------------
    b608 = hw["bearings"]["b608_2rs"]["stainless_S608_2RS"]
    b688 = hw["bearings"]["b688_2rs"]
    lc = hw["bearings"]["load_check"]["value"]
    orr = hw["oring_equator"]["recommended"]["value"]
    ogd = hw["oring_equator"]["gland_design"]["value"]
    operf = hw["oring_equator"]["performance"]["value"]
    d.update(
        c0_608=f(b608["C0"]["value"], 0), c0_688=f(b688["C0"]["value"], 0), c0_ratio608=f(lc["C0_over_T_2p5_S608"], 1),
        c0_ratio688=f(lc["C0_over_T_2p5_688"], 1), seal_speed=f(hw["lip_seal"]["lip_surface_speed_at_750rpm"]["value"], 3),
        or_id=f(orr["ID"], 2), or_cs=f(orr["cs"], 2), or_groove=f(ogd["spigot_groove_bottom_D"], 2),
        or_bore=f(ogd["bore_D"], 2), or_spigot=f(ogd["spigot_OD"], 2), or_gland=f(ogd["gland_depth"], 3),
        or_bw=f(ogd["groove_bottom_width"], 1), or_tw=f(ogd["groove_top_width"], 3),
        or_sq_nom=f(100 * operf["squeeze_nominal"], 1), or_sq_lo=f(100 * operf["squeeze_worst"][0], 1),
        or_sq_hi=f(100 * operf["squeeze_worst"][1], 1), or_fill_nom=f(100 * operf["fill_nominal"], 1),
        or_fill_worst=f(100 * operf["fill_worst"], 1), or_stretch_lo=f(100 * operf["stretch_range"][0], 1),
        or_stretch_hi=f(100 * operf["stretch_range"][1], 1), or_wall=f(ogd["female_wall_min_under_pockets"], 2),
        or_reserve=f(ogd["pocket_reserve"], 2),
    )

    # ---- prior art ------------------------------------------------------------------------
    pac = {r["product"]: r for r in hw["prior_art_power_check"]["value"]}
    rows = []
    for pa_ in hw["prior_art"]:
        name = pa_["product"]
        rd = pa_["rotor_D"]["value"]
        pw = pa_["power"]["value"]
        if isinstance(pw, list):
            pws = f"{f(pw[0], 1)}-{f(pw[1], 1)} W"
        elif isinstance(pw, dict):
            pws = ", ".join(f"{f(v, 0)} W at {k.split('_')[0]} mph" for k, v in pw.items())
        else:
            pws = f"{f(pw, 0)} W"
        mass = pa_.get("mass", {}).get("value")
        size = pa_.get("stowed", pa_.get("body", {})).get("value") if ("stowed" in pa_ or "body" in pa_) else None
        price = pa_.get("price", {}).get("value", {})
        prs = ", ".join(f"{k} {f(v, 0)}" for k, v in price.items()) if price else "n/a"
        short = name.split(" (")[0]
        pc_key = {"Enomad": "Enomad", "WaterLily": "WaterLily", "Blue": "Blue Freedom", "Ampair": "Aquair UW @6"}[short.split()[0]]
        vneed = next(v for k, v in pac.items() if k.startswith(pc_key))
        rows.append([short, f"{f(rd, 0)} mm" if rd else "n/a (unit 178 mm)", pws,
                     f"{f(mass, 0)} g" if mass else "n/a",
                     f"{f(size[0], 0)} x {f(size[1], 0)} mm" if size else "n/a",
                     pa_.get("stowage", ""), prs,
                     f"{f(vneed['V_needed_at_Cp_el_0.30_m_s'], 2)}-{f(vneed['V_needed_at_Cp_el_0.20_m_s'], 2)}"])
    t["prior_art"] = rows
    ovo_pc = next(v for k, v in pac.items() if k.startswith("Ovo"))
    d.update(ovo_vneed_lo=f(ovo_pc["V_needed_at_Cp_el_0.30_m_s"], 2), ovo_vneed_hi=f(ovo_pc["V_needed_at_Cp_el_0.20_m_s"], 2),
             ovo_pavail=f(ovo_pc["P_avail_W"], 1))
    est = {p_["product"]: p_ for p_ in hw["prior_art"]}["Enomad Estream"]
    est_cyl = math.pi * (est["stowed"]["value"][1] / 20) ** 2 * (est["stowed"]["value"][0] / 10) / 1000.0
    d.update(estream_cyl_l=f(est_cyl, 2), estream_len=f(est["stowed"]["value"][0], 0), estream_d=f(est["stowed"]["value"][1], 0))
    refs = hw["references"]
    d["chica_cp"] = re.search(r"Cp ([0-9.]+)", refs[0]["value"]).group(1)

    # ---- BOM and cost -----------------------------------------------------------------------
    tot = bom["totals"]
    fx = bom["meta"]["fx"]["usd_per_eur"]
    d.update(bom_disp_pro=f(tot["display_build_known_eur"], 2), bom_disp_buy=f(tot["display_build_purchase_eur"], 2),
             bom_v2_pro=f(tot["v2_additional_known_eur"], 2), bom_v2_buy=f(tot["v2_additional_purchase_eur"], 2),
             bom_unpriced=str(tot["unpriced_items"]), fx=f(fx, 4),
             fx_date=re.search(r"\((\d{4}-\d\d-\d\d)\)", hw["meta"]["fx"]["unit"]).group(1))
    # as-built bought parts (assembly.json BOM), grouped by part type
    groups: dict[str, dict] = {}

    def add(key, label, qty, unit_eur=None, v2=False, basis=""):
        gk = groups.setdefault(key, dict(label=label, qty=0, unit=unit_eur, v2=v2, basis=basis))
        gk["qty"] += qty

    def usd(x: float) -> float:
        return x / fx

    p608 = usd(b608["price"]["value"]["USD"])
    p688 = usd(b688["price"]["value"]["stainless_SS688_2RS_USD"])
    pmag = usd(hw["magnets"]["n52_4x2"]["price"]["value"]["USD_each"])
    pins = hw["heat_set_insert_m3"]["price"]["value"]["EUR_per_100_3djake"] / 100
    pgland = usd(hw["cable_gland"]["price"]["value"]["USD_per_10"]) / 10
    pgl40 = usd(gl["price"]["value"]["USD"])
    for b in asm["bom"]:
        n, qty = b["name"], b["qty"]
        if b.get("part_of"):
            continue
        if n.startswith("hinge_pin"):
            add("pin22", "Dowel pin ISO 8734 3 x 22", qty)
        elif n.startswith("pin_filler"):
            add("pin4", "Dowel pin ISO 8734 3 x 4 (pin retention)", qty)
        elif n.startswith(("hub_magnet", "blade_magnet_heel")):
            add("mag6", "N52 block magnet 6 x 2 x 2 (deployed detent)", qty)
        elif n.startswith(("blade_magnet_stow", "tail_magnet")):
            add("mag4", "N52 disc magnet 4 x 2 (stowed detent)", qty, pmag, basis="USD 0.30 each, MOQ 50")
        elif n.startswith(("rail_magnet", "fin_magnet", "fin_blank_magnet")):
            add("mag442", "N52 block magnet 4 x 2 x 2 (fin rails)", qty)
        elif n == "bearing_608":
            add("b608", "Ball bearing S608-2RS, 8 x 22 x 7, stainless", qty, p608, basis="USD 5.63")
        elif n == "bearing_688":
            add("b688", "Ball bearing SS688-2RS, 8 x 16 x 5, stainless", qty, p688,
                basis=f"USD {b688['price']['value']['stainless_SS688_2RS_USD']:.2f}")
        elif n == "shaft":
            add("shaft", "Shaft 8 x 75.5, SUS440C g6, machined flats", qty)
        elif n == "shaft_collar":
            add("collar", "Set-screw collar Ruland MSC-8-SS", qty)
        elif n.startswith("shim"):
            add("shim", "Shim 8 x 12 x 0.5, stainless", qty)
        elif n in ("hub_grub", "spinner_screw") or n.startswith("screws"):
            add("screws", "M3 screws (grub, cap, button head)", qty)
        elif n.startswith("inserts"):
            add("ins", "Heat-set insert M3 x 5.7", qty, pins, basis="EUR 12.09 per 100")
        elif n == "lip_seal_tc8x16x5":
            add("seal", "Lip seal TC 8x16x5 NBR (v2)", qty, v2=True)
        elif n == "oring_141":
            add("oring", "O-ring ISO 3601-141, 58.42 x 2.62", qty)
        elif n == "cable_gland_m12":
            add("gland", "Cable gland Hummel HSK-K M12x1.5 (v2)", qty, pgland, v2=True, basis="USD 10.48 per 10")
        elif n == "tool_seal_sleeve":
            add("sleeve", "Brass seal-fitting sleeve (tool, v2)", qty, v2=True)
    add("gl40", "Generator CubeMars GL40 KV70 (v2; display: printed replica)", 1, pgl40, v2=True, basis="USD 72.99")
    rows, disp_known, n_unpriced = [], 0.0, 0
    for k, gk in groups.items():
        unit = gk["unit"]
        tot_ = unit * gk["qty"] if unit is not None else None
        if unit is None:
            n_unpriced += 1
        elif not gk["v2"]:
            disp_known += tot_
        rows.append([gk["label"], str(gk["qty"]), f(unit, 2) if unit is not None else "not priced",
                     f(tot_, 2) if tot_ is not None else "", "v2" if gk["v2"] else "display + v2", gk["basis"]])
    t["bom"] = rows
    d.update(bom_asbuilt_known=f(disp_known, 2), bom_asbuilt_unpriced=str(n_unpriced),
             n_magnets=str(groups["mag6"]["qty"] + groups["mag4"]["qty"] + groups["mag442"]["qty"]))

    # ---- CAD: parts, checks, outputs ----------------------------------------------------------
    parts = asm["parts"]
    stl = asm["stl"]
    pp_ = asm["print_plan"]
    inter = asm["interference"]
    kc = asm["key_clearances"]
    ss = asm["spin_sweep"]
    web = asm["web"]
    state_names = {"stowed": "stowed (fold 90)", "deployed": "deployed (fold 0)", "fold_30": "fold 30",
                   "fold_60": "fold 60", "cutaway": "cutaway", "packed": "packed (fins off)", "v2": "v2 nose"}
    t["interference"] = [[state_names.get(k, k), f(v["fold_deg"], 0), str(v["n_parts"]), str(v["n_pairs_near"]),
                          str(len(v["expected_overlaps"])), str(len(v["contacts"])), str(len(v["defects"]))]
                         for k, v in inter.items()]
    want = [
        ("hub gap", "Hub gap, static tail to hub front"),
        ("stowed aero span", "Stowed aero span to pocket floors and walls"),
        ("stowed root", "Stowed root to static tail"),
        ("stowed detent magnet pair", "Stowed detent magnet pair: gap, angle, force"),
        ("deployed detent magnet pair", "Deployed detent pair: gap, force"),
        ("deployed heel", "Deployed heel on hub stop face"),
        ("blades folding 0..90 deg at the stowing index to the fins", "Folding blades to fins"),
        ("blades folding 6..90 deg to the hub", "Folding blades to hub (run clearance)"),
        ("blades folding 0..90 deg to the hinge pins", "Knuckle bore on hinge pin"),
        ("neighbour blades", "Neighbouring blades while folding"),
        ("claw A to claw B", "Claw coupling dog flanks"),
        ("front / aft shell locked", "Shells locked at the bayonet stop"),
        ("closing path", "Closing the shells over folded blades"),
        ("front shell material proud", "Front shell proud of the hull at the fin rails"),
    ]
    rows = []
    for key, label in want:
        row = next((k for k in kc if k["check"].startswith(key)), None)
        if row is None:
            raise KeyError(key)
        rows.append([label, clean_check(row["value"]), clean_check(re.sub(r"\s*\((?![^)]*tessellation)[^)]*\)", "", str(row.get("limit") or ""))),  # keep the mesh allowance visible
                     "pass" if row["ok"] else "FAIL"])
    t["clearances"] = rows
    n_kc_ok = sum(1 for k in kc if k["ok"])
    mod_n = {"shells": len(mods["shells"]["interface_checks"]), "rotor": len(mods["rotor"]["interface_checks"]),
             "drivetrain": len(mods["drivetrain"]["interfaces"]), "fins": len(mods["fins"]["interfaces"]),
             "coupon": len(mods["coupon"]["interfaces"]) + len(mods["coupon"]["features"])}
    mod_ok = {"shells": sum(c["ok"] for c in mods["shells"]["interface_checks"]),
              "rotor": sum(c["ok"] for c in mods["rotor"]["interface_checks"]),
              "drivetrain": sum(c["ok"] for c in mods["drivetrain"]["interfaces"]),
              "fins": sum(c["ok"] for c in mods["fins"]["interfaces"]),
              "coupon": sum(c["ok"] for c in mods["coupon"]["interfaces"]) + sum(c["ok"] for c in mods["coupon"]["features"])}
    spin_min = min(ss["pairs"], key=lambda p_: p_["min_gap_mm"] if not p_["turns_with_shaft"] and p_["min_gap_mm"] > 0 else 1e9)
    printed = [p_ for p_ in pp_["parts"]]
    d.update(
        n_parts_audit=str(len(parts)), n_parts_total=str(len(parts) + 1), n_stl=str(len({s["stl"] for s in stl})),
        n_stl_entries=str(len(stl)), n_printed_total=str(len(stl)), n_bought_total=str(len(parts) + 1 - len(stl)),
        stl_watertight=str(sum(1 for s in stl if s["watertight"])), stl_folded=str(sum(s["folded_pairs_gt170"] for s in stl)),
        stl_volmatch=str(sum(1 for s in stl if s["volume_match"])), stl_fit=str(sum(1 for s in stl if s["fits_bed"])),
        n_states=str(len(inter)), n_defects=str(sum(len(v["defects"]) for v in inter.values())),
        n_pairs_total=str(sum(v["n_pairs_near"] for v in inter.values())),
        kc_ok=str(n_kc_ok), kc_n=str(len(kc)),
        spin_pairs=str(len(ss["pairs"])), spin_rub=str(len(ss["rubbing"])), spin_step=f(ss["step_deg"], 0),
        spin_min=f(spin_min["min_gap_mm"], 3), spin_min_pair=f"{spin_min['rotor']} / {spin_min['static']}",
        glb_mb=f(web["glb_mb"], 2), glb_tris=f(web["triangles"], 0), glb_meshes=str(web["n_parts"]),
        glb_folded=str(web["mesh_quality"]["glb"]["folded_pairs_gt170_total"]),
        viewer_diff=f(asm["viewer_fold_check"]["max_abs_diff_mm"], 4),
        pivot_exp=str(math.ceil(math.log10(max(p_["max_matrix_diff"] for p_ in asm["pivot_proof"])))).replace("-", "−"),
        explode_ratio=f(asm["exploded_extent"]["ratio"], 2),
        n_step_parts=str(inter["stowed"]["n_parts"]), n_manifest=str(len(manifest["parts"])),
        print_g=f(pp_["total_grams_est"], 0), print_pla=f(pp_["grams_by_material"]["PLA"], 0),
        print_petg=f(pp_["grams_by_material"]["PETG"], 0), n_plates=str(len(pp_["plates"])),
        n_printed=str(len(printed)),
        assembly_s=f(asm["elapsed_s"], 0),
        blade_vol=f(blade_ck["mass"]["volume_mm3"] / 1000, 2), blade_mass=f(blade_ck["mass"]["mass_g_PLA_solid"], 1),
        blade_te_min=f(blade_ck["trailing_edge"]["min_mm"], 3),
        mod_checks=", ".join(f"{k} {mod_ok[k]}/{mod_n[k]}" for k in mod_n),
        mod_total_ok=str(sum(mod_ok.values())), mod_total=str(sum(mod_n.values())),
        n_steps=str(len(asm["assembly_steps"])),
    )
    t["plates"] = [[p_["plate"].split(" ", 1)[0], p_["plate"].split(" ", 1)[1],
                    ", ".join(x["name"].replace("tolerance_coupon", "coupon") for x in p_["parts"]),
                    f(p_["max_height_mm"], 1), f(p_["grams_est"], 0)]
                   for p_ in pp_["plates"]]
    d["kc_rows_note"] = f"{len(want)} of {len(kc)}"

    # numbers quoted in the prose from the key-clearance rows (values are strings in assembly.json)
    def kcv(prefix: str) -> str:
        return next(k["value"] for k in kc if k["check"].startswith(prefix))

    num = r"([0-9]+(?:\.[0-9]+)?)"
    m = re.search(num + r" / " + num + r" mm, " + num + r" deg.*?about " + num + " N", kcv("stowed detent magnet pair"))
    d.update(stow_gap=m[1], stow_min=m[2], stow_ang=m[3], stow_force=m[4])
    d["fold_fin_clear"] = re.search(num + " mm", kcv("blades folding 0..90 deg at the stowing index to the fins"))[1]
    d["aero_clear"] = re.search(num + " mm", kcv("stowed aero span"))[1]
    d["root_clear"] = re.search(num + " mm", kcv("stowed root"))[1]
    d["boss_clear"] = re.search(num + " mm", kcv("stowed blade to its detent boss"))[1]
    d["knuckle_clear"] = re.search(num + " mm", kcv("blades folding 0..90 deg to the hinge pins"))[1]
    d["claw_flank"] = re.search(num + " mm", kcv("claw A to claw B"))[1]
    m = re.search(r"fold <= " + num + r" deg.*fold <= " + num + " deg", kcv("closing path"))
    d.update(close_fins_off=m[1], close_fins_on=m[2])
    d["hub_d_clear"] = re.search(num + " mm", kcv("shaft to hub D-bore"))[1]

    # assorted values quoted in the prose
    bd_ = load("data/blade_design.json")
    lay_ = load("data/layout.json")
    hw_ = load("data/hardware.json")
    v2n = lay_["internal"]["mass_buoyancy"]["v2_needs"]
    crit = lay_["internal"]["drivetrain"]["bearings"]["criteria"]
    capc = bd_["generator_coupling_estimate"]["load_line_by_flow"]["cap_torque_curve"]
    tow_rec = hw_["generator"]["recommendation"]["tow_protection"]["value"]["recommended"]
    meth = asm["conventions"]["method"]
    og = hw_["oring_equator"]["design_guidance"]["value"]
    otab = hw_["oring_equator"]["gland_table_2p62"]["value"]
    d.update(
        v2_target=f(v2n["target_net_buoyancy_N"], 1), v2_shed_g=f(v2n["mass_to_shed_or_float_g"], 1),
        v2_foam_cc=f(v2n["option_external_foam_cc"], 0), v2_4per_g=f(v2n["option_4_perimeter_shells_saves_g"], 1),
        rotor_ratio=f(P.ROTOR_R_TIP_TARGET * 2 / P.HULL_D, 1), x_tip_r=f(P.SPINNER_CAP["x_tip_built"], 0),
        cone_mean=f(lay_["rotor"]["axis"]["cone_mean_hinge_to_tip_deg"], 1),
        lug=f"{f(bd_['shank']['lug_chord_mm'], 0)} × {f(bd_['shank']['lug_thickness_mm'], 0)}",
        bay_twist=f(P.BAYONET["twist_deg"], 0), nose_bridge=f(P.NOSE_PLUG["bridge_d"], 1),
        cradle_gap=f(P.BLADE_CRADLE["gap"], 1), n_candidates=str(len(hw_["generator"]["candidates"])),
        cav_freewheel_flow=f(capc["flow_m_s"][0], 1),
        input_v=re.search(r">= ?([0-9]+) V", tow_rec)[1],
        crit_span=f(crit["span_min_mm"], 0), crit_runout=f(crit["tip_runout_max_mm"], 1),
        crit_s0=f(crit["s0_min"], 0), crit_defl=f(crit["shaft_deflection_max_mm"], 2),
        crit_ratio=f(crit["span_over_overhang_reference"], 2), crit_slop=f(crit["seat_slop_mm"], 2),
        mesh_chord=re.search(r"\(([0-9.]+) mm chord", meth)[1], mesh_chord_shell=re.search(r"shells ([0-9.]+)\)", meth)[1],
        overlap_tol=re.search(r"overlaps <= ([0-9.]+) mm3", meth)[1],
        og_sq0=str(og["squeeze_pct"][0]), og_sq1=str(og["squeeze_pct"][1]), og_fill0=str(og["fill_pct"][0]),
        og_fill1=str(og["fill_pct"][1]), tab_sq0=str(otab["squeeze_pct"][0]),
        seal_hrc=str(hw_["lip_seal"]["installation"]["value"]["shaft_hardness_HRC_min"]),
        print_tol=f(hw_["meta"]["assumptions"]["print_tol_diameter_mm"], 2),
        radial_pct=re.search(r"radial load ([0-9]+) %", lay_["internal"]["drivetrain"]["bearings"]["load_basis"])[1],
        asym_pct=re.search(r"one blade ([0-9]+) % off", lay_["internal"]["drivetrain"]["bearings"]["load_basis"])[1],
        ct_tow=re.search(r"C_T ([0-9.]+)", lay_["internal"]["drivetrain"]["bearings"]["load_basis"])[1],
    )

    # as-built fin span in plan (the layout study's tip-to-tip value is quoted beside it, A10)
    fin_y = max(max(abs(p_["bbox_mm"][1]), abs(p_["bbox_mm"][4])) for p_ in parts if p_["name"] in ("fin_0", "fin_1"))
    d["fin_plan"] = f(2 * fin_y, 0)
    # prior art: Ovo's flow for rating inverted from the same quantity as the products' ratings, its DC
    # output (A12), not from the shaft power
    p_dc = td["td_bridge_1.5_m_s"]["P_el_dc_W"]
    vneed = [(p_dc / (cp_el * 0.5 * P.WATER_RHO * area)) ** (1 / 3) for cp_el in (0.30, 0.20)]
    d.update(ovo_vneed_dc_lo=f(vneed[0], 2), ovo_vneed_dc_hi=f(vneed[1], 2))
    # which CAD build the numbers and renders come from (H2): the assembly's cache key
    d["cad_key"] = asm.get("cache_key", "unknown")
    d["licence"] = LICENCE
    d["rail_proud"] = re.search(r"([0-9]+(?:\.[0-9]+)?) mm3", kcv("front shell material proud"))[1]

    # printing and fit parameters (params.py)
    d.update(clear_slide=f(P.CLEAR_SLIDE, 2), clear_rotate=f(P.CLEAR_ROTATE, 2), run_clear=f(P.RUN_CLEAR, 2),
             press_fit=f(P.PRESS_FIT, 2), min_wall=f(P.MIN_WALL, 1), max_overhang=f(P.MAX_OVERHANG, 0),
             bed=f(P.BED[0], 0), nozzle=f(P.NOZZLE, 1), layer=f(P.LAYER, 1), min_chord=f(P.MIN_CHORD, 0),
             max_chord=f(P.BLADE_MAX_CHORD, 0), max_root_twist=f(P.MAX_ROOT_TWIST, 0),
             pin=f"{f(P.PIN_D, 0)} x {f(P.PIN_L, 0)}", shaft_l=f(P.SHAFT["L"], 1))
    return d, t


# --------------------------------------------------------------------------------------
# bibliography
# --------------------------------------------------------------------------------------
CURATED = r"""
@article{myring1976,
  author = {Myring, D. F.},
  title = {A theoretical study of body drag in subcritical axisymmetric flow},
  journal = {Aeronautical Quarterly},
  volume = {27},
  number = {3},
  pages = {186--194},
  year = {1976},
}
@book{burton2011,
  author = {Burton, Tony and Jenkins, Nick and Sharpe, David and Bossanyi, Ervin},
  title = {Wind Energy Handbook},
  edition = {2},
  publisher = {John Wiley \& Sons},
  year = {2011},
}
@book{hansen2015,
  author = {Hansen, Martin O. L.},
  title = {Aerodynamics of Wind Turbines},
  edition = {3},
  publisher = {Routledge},
  year = {2015},
}
@techreport{buhl2005,
  author = {Buhl, Marshall L.},
  title = {A New Empirical Relationship between Thrust Coefficient and Induction Factor for the Turbulent Windmill State},
  institution = {National Renewable Energy Laboratory},
  number = {NREL/TP-500-36834},
  year = {2005},
}
@article{ning2014,
  author = {Ning, S. Andrew},
  title = {A simple solution method for the blade element momentum equations with guaranteed convergence},
  journal = {Wind Energy},
  volume = {17},
  number = {9},
  pages = {1327--1345},
  year = {2014},
}
@article{neuralfoil,
  author = {Sharpe, Peter D. and Hansman, R. John},
  title = {{NeuralFoil}: An Airfoil Aerodynamics Analysis Tool Using Physics-Informed Machine Learning},
  journal = {arXiv preprint arXiv:2503.16323},
  year = {2025},
  url = {https://github.com/peterdsharpe/NeuralFoil},
}
@inproceedings{drela1989,
  author = {Drela, Mark},
  title = {{XFOIL}: An Analysis and Design System for Low Reynolds Number Airfoils},
  booktitle = {Low Reynolds Number Aerodynamics},
  series = {Lecture Notes in Engineering},
  volume = {54},
  pages = {1--12},
  publisher = {Springer},
  year = {1989},
}
@book{gasch2012,
  author = {Gasch, Robert and Twele, Jochen},
  title = {Wind Power Plants: Fundamentals, Design, Construction and Operation},
  edition = {2},
  publisher = {Springer},
  year = {2012},
}
@book{carlton2019,
  author = {Carlton, John},
  title = {Marine Propellers and Propulsion},
  edition = {4},
  publisher = {Butterworth-Heinemann},
  year = {2019},
}
@misc{build123d,
  title = {build123d: a Python parametric boundary-representation CAD library on the OpenCascade kernel},
  url = {https://github.com/gumyr/build123d},
  note = {Version 0.13},
}
@misc{manifold,
  title = {Manifold: geometry library for topological robustness},
  url = {https://github.com/elalish/manifold},
  note = {Python bindings manifold3d},
}
@misc{typst,
  title = {Typst: a markup-based typesetting system},
  url = {https://typst.app},
}
"""


def bib_escape(s: str) -> str:
    return s.replace("&", r"\&").replace("%", r"\%").replace("#", r"\#").replace("_", r"\_")


def write_bib() -> dict:
    """Curated references plus the sourced URLs of data/hardware.json, keyed for the report."""
    hw = load("data/hardware.json")
    acc = hw["meta"]["date"]
    gl = next(c for c in hw["generator"]["candidates"] if c["name"] == "CubeMars GL40 KV70")
    gm = next(c for c in hw["generator"]["candidates"] if c["name"].startswith("iPower GM4108H"))
    entries = {
        "gl40_drawing": ("CubeMars GL40 gimbal motor, 2D drawing", gl["mount"]["source_url"]),
        "gl40_page": ("CubeMars GL40 KV70 product page (KV, resistance, inductance, ratings)", gl["KV"]["source_url"]),
        "gm4108_page": ("iFlight iPower GM4108H-120T product page", gm["KV"]["source_url"] if isinstance(gm.get("KV"), dict) else hw["generator"]["gm4108_datasheet_check"]["source_url"]),
        "applerubber": ("Apple Rubber, ISO 3601 metric size O-rings and static radial gland design", hw["oring_equator"]["recommended"]["source_url"]),
        "oring_guide": ("Global O-Ring and Seal, O-ring groove design guide", hw["oring_equator"]["design_guidance"]["source_url"]),
        "bearing608": ("BearingsDirect, S608-2RS stainless steel bearing 8x22x7 (load ratings, price)", hw["bearings"]["b608_2rs"]["d"]["source_url"]),
        "bearing688": ("BearingsDirect, 688-2RS miniature bearing 8x16x5 (load ratings, price)", hw["bearings"]["b688_2rs"]["price"]["source_url"]),
        "misumi": ("MISUMI, linear shafts catalogue (PSSFJ, SUS440C, g6)", hw["shaft"]["recommended"]["source_url"]),
        "lipseal_install": ("Oil Seal Stocks, TC rotary shaft seal installation guideline", hw["lip_seal"]["installation"]["source_url"]),
        "hummel": ("Hummel HSK-K M12x1.5 cable gland, product data", hw["cable_gland"]["hummel_hsk_k_m12"]["source_url"]),
        "iso8734": ("ISO 8734:1997, Parallel pins of hardened steel and martensitic stainless steel (sample)", hw["dowel_pins"]["standard"]["source_url"]),
        "insert": ("CNC Kitchen M3 x 5.7 heat-set inserts (3DJake listing)", hw["heat_set_insert_m3"]["cnc_kitchen_standard"]["source_url"]),
        "magnets": ("JC Magnetics, N52 4 x 2 mm disc magnets", hw["magnets"]["n52_4x2"]["price"]["source_url"]),
        "prusa_fits": ("Prusa Knowledge Base, Modeling with 3D printing in mind", hw["print_fits"]["prusa_moving_parts_gap"]["source_url"]),
        "creative3dp": ("Creative3DP, Press-fit tolerances for 3D printing", hw["print_fits"]["creative3dp_ladder_diametral"]["source_url"]),
        "prusa_petg": ("Prusa Knowledge Base, PETG", hw["materials"]["petg"]["source_url"]),
        "bambu_x1c": ("Bambu Lab X1 Carbon technical specifications", hw["printers"]["x1c"]["source_url"]),
        "ecb": ("European Central Bank, euro foreign exchange reference rates (2026-10-02)", hw["meta"]["fx"]["source_url"]),
    }
    pa = {p["product"]: p for p in hw["prior_art"]}
    entries["estream"] = ("New Atlas, Estream portable water turbine charger (Enomad)", pa["Enomad Estream"]["rotor_D"]["source_url"])
    entries["waterlily"] = ("WaterLily Turbine, paddling charger product page", pa["WaterLily (discontinued)"]["power"]["source_url"])
    entries["waterlily_review"] = ("World Water Reserve, WaterLily turbine review", pa["WaterLily (discontinued)"]["mass"]["source_url"])
    entries["bluefreedom"] = ("Blue Freedom, Portable hydro power plant product page", pa["Blue Freedom Portable"]["rotor_D"]["source_url"])
    entries["bluefreedom_review"] = ("Gigazine, Blue Freedom review (machine translated)", pa["Blue Freedom Portable"]["speed"]["source_url"])
    entries["aquair"] = ("No Outage, Ampair Aquair UW water generator", pa["Ampair Aquair UW (towed/stream, not for sale)"]["rotor_D"]["source_url"])
    lines = [CURATED.strip(), ""]
    for k, (title, url) in entries.items():
        lines.append(f"@misc{{{k},\n  title = {{{bib_escape(title)}}},\n  url = {{{url}}},\n  urldate = {{{acc}}},\n}}")
    r0, r1 = hw["references"]
    # titles and authors checked against the publisher / RePEc / ScienceDirect records (2026-10-05)
    lines.append("@article{chica2025,\n  author = {Chica, Edwin and Vel{\\'a}squez, Laura and Rubio-Clemente, Ainhoa},\n"
                 "  title = {Full-Scale Experimental Assessment of a Horizontal-Axis Hydrokinetic Turbine for River "
                 "Applications: A Challenge for Developing Countries},\n"
                 "  journal = {Energies},\n  volume = {18},\n  number = {7},\n  pages = {1657},\n  year = {2025},\n"
                 f"  doi = {{{r0['source_url'].split('doi.org/')[1]}}},\n}}")
    lines.append("@article{riglin2016,\n  author = {Riglin, Jacob and Carter, Fred and Oblas, Nicholas and Schleicher, W. Chris "
                 "and Daskiran, Cosan and Oztekin, Alparslan},\n"
                 "  title = {Experimental and numerical characterization of a full-scale portable hydrokinetic turbine "
                 "prototype for river applications},\n"
                 "  journal = {Renewable Energy},\n  volume = {99},\n  pages = {772--783},\n  year = {2016},\n"
                 f"  doi = {{{r1['source_url'].split('doi.org/')[1]}}},\n}}")
    (REPORT / "refs.bib").write_text("\n".join(lines) + "\n")
    return entries


# --------------------------------------------------------------------------------------
# figures (style shared with ovo/bem.py)
# --------------------------------------------------------------------------------------
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e2e1dc", "#fcfcfb"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
HULL_FILL, HUB_FILL = "#e6e5e1", "#f2f1ee"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelsize": 9, "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.8,
        "lines.solid_capstyle": "round", "legend.frameon": False, "legend.fontsize": 8.5,
        "font.size": 9, "text.color": INK,
    })
    return plt


def fig_hull(path: Path) -> None:
    plt = _style()
    x = np.linspace(0, P.HULL_L, 1200)
    r = profile.radius(x)
    ri = profile.inner_radius(x)
    a, b, c = P.HULL_NOSE_A, P.HULL_MID_B, P.HULL_TAIL_C
    R = P.HULL_D / 2
    fig, ax = plt.subplots(figsize=(7.4, 3.25))
    static = x <= P.STATIC_TAIL_X
    ax.fill_between(x[static], -r[static], r[static], color=HULL_FILL, lw=0, zorder=1)
    ax.fill_between(x[~static], -r[~static], r[~static], color=HUB_FILL, lw=0, zorder=1, hatch="////",
                    edgecolor="#d8d7d2")
    ax.plot(x, r, color=INK, lw=1.4, zorder=3)
    ax.plot(x, -r, color=INK, lw=1.4, zorder=3)
    ax.plot(x, ri, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=3)
    ax.plot(x, -ri, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=3)
    ax.annotate(f"inner wall\n({P.HULL_WALL:g} mm shell)", xy=(14, -float(profile.inner_radius(14.0))),
                xytext=(-17, -R - 8), fontsize=7.5, color=INK2, va="top",
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    for xs in (a, a + b):
        ax.plot([xs, xs], [R - 4, R + 6], color=INK2, lw=0.6, zorder=2)
    yl = R + 6
    ax.annotate("", xy=(0, yl), xytext=(a, yl), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
    ax.text(0, yl + 2.2, f"nose: super-ellipse\na {a:g}, n {P.HULL_NOSE_N:g}", ha="left", va="bottom",
            fontsize=7.5, color=INK2)
    ax.annotate("", xy=(a + b, yl), xytext=(P.HULL_L, yl), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
    ax.text(P.HULL_L, yl + 2.2, f"tail: Myring cubic\nc {c:g}, tip half-angle {P.HULL_TAIL_THETA:g} deg",
            ha="right", va="bottom", fontsize=7.5, color=INK2)
    ax.text(a + b / 2, R - 5, f"b {b:g}", ha="center", va="top", fontsize=7, color=INK2)
    for xs, lab, col in ((P.EQUATOR_X, f"equator x {P.EQUATOR_X:g}", S1),
                         (P.STATIC_TAIL_X, f"static tail x {P.STATIC_TAIL_X:g}", S2)):
        ax.plot([xs, xs], [-R - 2, R - 6], color=col, lw=1.0, ls=(0, (4, 2)), zorder=2)
        ax.text(xs - 1.3, -R + 1, lab, color=col, fontsize=7.5, rotation=90, ha="right", va="bottom", zorder=4)
    xt = P.SPINNER_CAP["x_tip_built"]
    ax.plot([xt], [0], marker="o", ms=3.5, color=INK, zorder=4)
    ax.annotate(f"built tip x {xt:.2f} (R{P.SPINNER_TIP_R:g})", xy=(xt, 0), xytext=(128, -30), fontsize=7.5,
                color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.plot([P.HINGE_X], [P.HINGE_R], marker="o", ms=4.5, color=S2, zorder=4)
    ax.annotate(f"blade hinge ({P.HINGE_X:g}, {P.HINGE_R:g})", xy=(P.HINGE_X, P.HINGE_R), xytext=(124, 30),
                fontsize=7.5, color=S2, arrowprops=dict(arrowstyle="-", color=S2, lw=0.6))
    ax.annotate("", xy=(30, R), xytext=(30, -R), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
    ax.text(31.5, 2, f"D {P.HULL_D:g}", fontsize=7.5, color=INK2)
    ax.text(82, 6, "static hull\n(two printed\nshells)", fontsize=7.5, color=INK2, ha="center")
    ax.text(124, -1, "rotating\nhub, spinner", fontsize=7, color=INK2, ha="center", va="top")
    ax.annotate("", xy=(-4, 14), xytext=(-17, 14), arrowprops=dict(arrowstyle="-|>", color=S1, lw=1.2))
    ax.text(-17, 17, "flow", color=S1, fontsize=7.5)
    ax.set_xlim(-19, P.HULL_L + 5)
    ax.set_ylim(-R - 20, R + 21)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm), nose at 0")
    ax.set_ylabel("r (mm)")
    ax.set_title(f"Hull profile r(x): L {P.HULL_L:g} mm, D {P.HULL_D:g} mm, L/D {P.HULL_L / P.HULL_D:.2f}")
    fig.savefig(path, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def fig_concept(path: Path) -> None:
    """Schematic: why the rotor sits downstream and its blades fold forward."""
    plt = _style()
    from matplotlib.patches import FancyArrowPatch
    x = np.linspace(0, P.HULL_L, 600)
    r = profile.radius(x)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.7))

    def hull(ax):
        ax.fill_between(x, -r, r, color=HULL_FILL, lw=0)
        ax.plot(x, r, color=INK, lw=1.0)
        ax.plot(x, -r, color=INK, lw=1.0)
        ax.annotate("", xy=(-4, -20), xytext=(-26, -20), arrowprops=dict(arrowstyle="-|>", color=S1, lw=1.3))
        ax.text(-26, -16, "flow", color=S1, fontsize=8)
        ax.set_xlim(-30, 190)
        ax.set_ylim(-42, 112)
        ax.set_aspect("equal")
        ax.axis("off")

    def lying(ax, x0, x1, hinge):
        xs = np.linspace(x0, x1, 60)
        ys = profile.radius(xs) + 1.5
        if hinge[0] < x0:
            px, py = np.r_[hinge[0], xs], np.r_[hinge[1], ys]
        else:
            px, py = np.r_[xs, hinge[0]], np.r_[ys, hinge[1]]
        ax.plot(px, py, color=S2, lw=4.2, solid_capstyle="round", zorder=3)
        ax.plot(*hinge, "o", color=INK, ms=3.5, zorder=4)

    def radial(ax, hinge, ang, length=80):
        a_ = math.radians(ang)
        p1 = (hinge[0] + length * math.cos(a_), hinge[1] + length * math.sin(a_))
        ax.plot([hinge[0], p1[0]], [hinge[1], p1[1]], color=S2, lw=2.6, ls=(0, (2, 1.5)), alpha=0.85, zorder=3)
        return p1

    # (a) rotor at the nose: folded blades lie aft along the hull, drag holds them there
    ax = axs[0]
    hull(ax)
    hinge = (10.0, 14.0)
    lying(ax, 14, 88, hinge)
    tip = (88, float(profile.radius(88.0)) + 1.5)
    ax.annotate("", xy=(tip[0] + 20, tip[1]), xytext=tip, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    ax.text(tip[0] + 4, tip[1] + 5, "drag pushes\nthe tip aft:\nblade held shut", fontsize=7.5, color=INK)
    dep = radial(ax, hinge, 95)
    ax.text(dep[0] + 4, dep[1] - 10, "to deploy, the tip\nmust swing upstream,\nagainst the drag", fontsize=7.5, color=INK2)
    ax.set_title("(a) rotor upstream, blades fold aft", fontsize=9)
    # (b) Ovo: rotor at the tail, blades fold forward into pockets; drag opens, thrust holds
    ax = axs[1]
    hull(ax)
    hinge = (P.HINGE_X, P.HINGE_R)
    lying(ax, 32, 104, hinge)
    tip = (32, float(profile.radius(32.0)) + 1.5)
    ax.annotate("", xy=(tip[0] + 18, tip[1] + 2), xytext=tip, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    ax.text(tip[0] - 26, tip[1] + 8, "drag on the\nfolded tip", fontsize=7.5, color=INK)
    dep = radial(ax, hinge, 80)
    arc = FancyArrowPatch(posA=(58, 52), posB=(dep[0] - 6, dep[1] - 12), connectionstyle="arc3,rad=-0.3",
                          arrowstyle="-|>", mutation_scale=9, color=S3, lw=1.3)
    ax.add_patch(arc)
    ax.text(2, 60, "opens the blade", fontsize=7.5, color=S3)
    ax.annotate("", xy=(dep[0] + 18, dep[1] - 16), xytext=(dep[0], dep[1] - 16),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    ax.text(dep[0] + 4, dep[1] - 12, "thrust holds it\non the stop", fontsize=7.5, color=INK)
    ax.set_title("(b) Ovo: rotor downstream, blades fold forward", fontsize=9)
    fig.tight_layout(w_pad=0.5)
    fig.savefig(path, dpi=200, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def fig_bay(path: Path) -> None:
    """Generator bay frontier, redrawn for an A4 column from data/layout_internal.json (the study figure's motor
    labels collide at print size): the largest bay diameter for each length, the chosen bay, the candidates."""
    plt = _style()
    gb = load("data/layout_internal.json")["generator_bay"]
    pick, fr = gb["pick"], gb["frontier"]
    fig, ax = plt.subplots(figsize=(6.6, 3.7))
    ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
        (20, 40), 15, 16, color="#eeede9", lw=0, zorder=0))
    ax.text(20.25, 40.4, "candidate box", ha="left", va="bottom", fontsize=8, color=INK2)
    Ls = [f_["L_mm"] for f_ in fr]
    Ds = [f_["D_max_mm"] for f_ in fr]
    ax.step(Ls + [Ls[-1] + 0.5], Ds + [Ds[-1]], where="post", color=S1, lw=2.0, zorder=2)
    ax.text(Ls[0], Ds[0] + 0.6, "largest D for each L", color=S1, fontsize=8.5, va="bottom")
    lmax = gb["rotor_face_x_max_mm"] - pick["x_fwd_mm"]
    ax.plot([lmax, lmax], [pick["bay_D_mm"], 62], color=S1, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.text(lmax + 0.3, 57.6, f"longer bays fail the\nbearing span (L <= {lmax:g})", color=S1, fontsize=8, va="top")
    bd, bl = pick["bay_D_mm"], pick["bay_L_mm"]
    ax.plot([17, bl, bl], [bd, bd, 38], color=S2, lw=1.8, zorder=2)
    ax.text(17.2, bd + 0.35, f"chosen bay D{bd:g} x L{bl:g}", color=S2, fontsize=8.5, va="bottom")
    # label offsets (points) chosen so no two labels touch at this size
    off = {"GL40 KV70": (6, -9), "MN4006 KV380": (6, -6), "GB4106 KV53": (-6, 5), "GM4108H-120T": (6, 3),
           "GBM4114HS 100T": (-6, -10), "GL40 II KV82.5": (6, 4), "RI50 KV100": (6, -3), "MN5008 KV170": (-6, -11),
           "GM5208-24": (-6, -11), "GB54-2 KV26": (6, -3)}
    for c in gb["candidate_fits"]:
        nm = re.sub(r"^(iPower|Rctimer|CubeMars|T-motor) ", "", c["name"].split(" (")[0])
        x, y = c["L_mm"], c["D_mm"]
        if c["fits_bay"]:
            ax.plot(x, y, "o", ms=6, color=INK, zorder=3)
        else:
            ax.plot(x, y, "x", ms=6.5, mew=1.8, color="#c8352b", zorder=3)
        dx, dy = off.get(nm, (6, 3))
        ax.annotate(nm, (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=8, color=INK2,
                    ha="left" if dx >= 0 else "right", va="center")
    ax.plot([], [], "o", ms=6, color=INK, label="fits the chosen bay")
    ax.plot([], [], "x", ms=6.5, mew=1.8, color="#c8352b", label="does not fit")
    ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax.set_xlim(17, 36.6)
    ax.set_ylim(38, 62)
    ax.set_xlabel("generator length L (mm)")
    ax.set_ylabel("generator diameter D (mm)")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def crop_figures() -> None:
    """Report-sized crops of two layout figures (the originals in docs/figures stay untouched)."""
    from PIL import Image
    src = REPO / "docs" / "figures"
    Image.open(src / "layout_fold_side.png").crop((90, 0, 1500, 900)).save(FIG / "fold_side_main.png", optimize=True)
    Image.open(src / "layout_section.png").crop((246, 51, 1725, 815)).save(FIG / "layout_section_drawing.png", optimize=True)
    # the fins figure's plan view only, so its legend prints at a legible size (the aft view's numbers go in the
    # caption); the joint figure split into its two panels, each at full width
    fins = Image.open(src / "layout_section_fins.png")
    fins.crop((0, 30, round(fins.width * 0.505), fins.height)).save(FIG / "fins_plan.png", optimize=True)
    jb = Image.open(src / "layout_section_joint_bay.png")
    cut = round(jb.width * 0.555)
    jb.crop((0, 0, cut, jb.height)).save(FIG / "joint_section.png", optimize=True)
    (FIG / "bay_frontier.png").unlink(missing_ok=True)       # redrawn by fig_bay instead


# --------------------------------------------------------------------------------------
# renders
# --------------------------------------------------------------------------------------
RENDER_KEYS = {
    "hero_deployed": ("hero_deployed", "hero_light", "hero", "deployed"),
    "stowed": ("stowed", "packed"),
    "cutaway": ("cutaway",),
    "exploded": ("exploded",),
    "packed_kit": ("packed_kit", "kit"),
}


def find_render(key: str) -> Path | None:
    media = REPO / "docs" / "media"
    if media.is_dir():
        for stem in RENDER_KEYS[key]:
            for ext in (".png", ".jpg", ".jpeg", ".webp"):
                hits = sorted(media.rglob(stem + ext))
                if hits:
                    return hits[0]
    p = REPO / "out" / "renders" / f"{key}.png"
    return p if p.exists() else None


DETAIL_KEYS = ("details_hub", "details_fin_rail")
RENDER_GROUPS = (("stowed", "packed_kit"),)    # shown side by side: take both from the same render set


def prepare_images() -> dict:
    from PIL import Image
    IMG.mkdir(parents=True, exist_ok=True)
    used = {}
    srcs = {key: find_render(key) for key in RENDER_KEYS}
    for group in RENDER_GROUPS:
        if len({srcs[k] is not None and "docs/media" in str(srcs[k]) for k in group}) > 1:
            for k in group:
                fallback = REPO / "out" / "renders" / f"{k}.png"
                srcs[k] = fallback if fallback.exists() else srcs[k]
    for key in RENDER_KEYS:
        src = srcs[key]
        if src is None:
            raise FileNotFoundError(f"no render for {key}: run scripts/build.py (out/renders) or the media script")
        im = Image.open(src).convert("RGB")
        if im.width > 1600:
            im = im.resize((1600, round(im.height * 1600 / im.width)), Image.LANCZOS)
        dst = IMG / f"{key}.jpg"
        im.save(dst, quality=86, optimize=True, progressive=True)
        used[key] = str(src.relative_to(REPO))
    # optional close-ups from docs/media (used only when every one of them exists)
    details = {k: REPO / "docs" / "media" / f"{k}.png" for k in DETAIL_KEYS}
    if all(p.exists() for p in details.values()):
        for k, src in details.items():
            im = Image.open(src).convert("RGB")
            if im.width > 1600:
                im = im.resize((1600, round(im.height * 1600 / im.width)), Image.LANCZOS)
            im.save(IMG / f"{k}.jpg", quality=86, optimize=True, progressive=True)
            used[k] = str(src.relative_to(REPO))
    return used


FOLD_QUERY = ("view=-40,22&theme=paper&rig=media&fit=tight&fill=0.40&focus=front_shell,aft_shell"
              "&fold={fold:.4f}")
FOLD_CROP = (190, 90, 950, 830)          # px in the 1000 x 1000 frames (camera fitted to the shells)


def render_fold() -> None:
    """Fold sequence 90 / 60 / 30 / 0 deg (stowed to deployed) with the camera fitted to the static shells."""
    FOLD_FRAMES.mkdir(parents=True, exist_ok=True)
    jobs = [{"out": str(FOLD_FRAMES / f"fold_{i}.png"), "q": FOLD_QUERY.format(fold=v)}
            for i, v in enumerate((0.0, 1 / 3, 2 / 3, 1.0))]       # stowed -> deployed, the site's order
    jf = FOLD_FRAMES / "jobs.json"
    jf.write_text(json.dumps(jobs, indent=1))
    subprocess.run([sys.executable, str(REPO / "scripts" / "render.py"), "--batch", str(jf), "-W", "1000", "-H", "1000"],
                   cwd=REPO, check=True)
    compose_fold([Path(j["out"]) for j in jobs])


def compose_fold(frames: list[Path]) -> None:
    from PIL import Image
    ims = [Image.open(p).convert("RGB").crop(FOLD_CROP) for p in frames]
    w, h = ims[0].size
    gap = 14
    strip = Image.new("RGB", (w * len(ims) + gap * (len(ims) - 1), h), "white")
    for i, im in enumerate(ims):
        strip.paste(im, (i * (w + gap), 0))
    strip = strip.resize((2400, round(h * 2400 / strip.width)), Image.LANCZOS)
    IMG.mkdir(parents=True, exist_ok=True)
    strip.save(IMG / "fold_sequence.jpg", quality=86, optimize=True, progressive=True)


# --------------------------------------------------------------------------------------
# compile
# --------------------------------------------------------------------------------------
def check_em_dash(paths: list[Path]) -> None:
    bad = [p for p in paths if EM_DASH in p.read_text()]
    if bad:
        raise SystemExit(f"em-dash (U+2014) found in {', '.join(str(b) for b in bad)}")
    src = SRC_TYP.read_text()
    if "---" in re.sub(r"//.*", "", src):
        raise SystemExit("'---' in the typst source renders as an em-dash")


def compile_pdf(png_dir: Path | None) -> None:
    import typst
    check_em_dash([SRC_TYP, REPORT / "data.typ", REPORT / "refs.bib"])
    typst.compile(str(SRC_TYP), output=str(PDF), root=str(REPO), ignore_system_fonts=True)
    print(f"wrote {PDF.relative_to(REPO)} ({PDF.stat().st_size / 1e6:.2f} MB)")
    if png_dir is not None:
        png_dir.mkdir(parents=True, exist_ok=True)
        typst.compile(str(SRC_TYP), output=str(png_dir / "page-{0p}.png"), root=str(REPO), format="png", ppi=110,
                      ignore_system_fonts=True)
        print(f"pages as PNG in {png_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--png", type=Path, help="also export every page as PNG into this directory")
    ap.add_argument("--render-fold", action="store_true", help="re-render the fold-sequence strip (heavy)")
    ap.add_argument("--no-compile", action="store_true", help="only regenerate data, bibliography and figures")
    a = ap.parse_args()

    REPORT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    d, t = build_data()
    write_bib()
    fig_hull(FIG / "hull_profile.png")
    fig_concept(FIG / "concept_fold.png")
    fig_bay(FIG / "generator_bay.png")
    crop_figures()
    if a.render_fold:
        render_fold()
    if not (IMG / "fold_sequence.jpg").exists():
        frames = sorted(FOLD_FRAMES.glob("fold_?.png"))
        if len(frames) == 4:
            compose_fold(frames)
        else:
            raise SystemExit("docs/report/img/fold_sequence.jpg missing: run with --render-fold")
    d["render_sources"] = prepare_images()
    d["has_details"] = all(k in d["render_sources"] for k in DETAIL_KEYS)
    mj = REPO / "docs" / "media" / "media.json"
    if mj.exists() and any("docs/media" in v for v in d["render_sources"].values()):
        mk = json.loads(mj.read_text()).get("cad_key")
        if mk != d["cad_key"]:
            print(f"WARNING: docs/media comes from CAD build {mk}, out/cad/checks from {d['cad_key']}: "
                  "run `make media` before the report")
    header = ("// GENERATED by scripts/report.py from data/*.json, out/cad/checks/*.json, out/web/manifest.json\n"
              "// and ovo/params.py. Do not edit: re-run the script.\n")
    (REPORT / "data.typ").write_text(header + f"#let d = {tlit(d)}\n\n#let tb = {tlit(t)}\n")
    print(f"wrote docs/report/data.typ ({len(d)} values, {len(t)} tables)")
    if not a.no_compile:
        compile_pdf(a.png)


if __name__ == "__main__":
    main()
