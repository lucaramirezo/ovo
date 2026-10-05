"""Regenerate every HydroEgg CAD output from scratch with one command.

  .venv/bin/python scripts/build.py                 # part modules, assembly, renders, summary table
  .venv/bin/python scripts/build.py --fresh         # also drop the shells and assembly BREP caches
  .venv/bin/python scripts/build.py --only-assembly # reuse the part modules' outputs (out/cad/*)
  .venv/bin/python scripts/build.py --summary       # only print the table from the existing checks
  .venv/bin/python scripts/build.py --no-render     # skip out/renders
  .venv/bin/python scripts/build.py --previews      # also each module's own review renders

Steps (each geometry job goes through scripts/cadrun.sh, at most 2 OCC jobs machine-wide):
  1. part modules, in dependency order: shells, blade, rotor, drivetrain, fins, coupon. Each writes
     out/cad/step/<part>.step, out/cad/stl/<part>.stl (print orientation) and out/cad/checks/<module>.json
  2. hydroegg.cad.assembly --force: rebuilds the part cache from the modules, interference matrix
     (stowed / deployed / fold 30 / fold 60 / cutaway / packed), key clearances, rotor spin sweep, STL audit,
     out/cad/step/hydroegg_assembly_{stowed,deployed}.step, out/web/hydroegg.glb + manifest.json,
     viewer fold check, out/cad/checks/assembly.json
  3. renders through the viewer (1600 x 1000): out/renders/{hero_deployed,stowed,packed_kit,cutaway,exploded,side,front}.png
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = str(REPO / ".venv" / "bin" / "python")
CADRUN = str(REPO / "scripts" / "cadrun.sh")

# (module, flags that skip its review renders)
MODULES = [("shells", ["--no-render"]), ("blade", ["--no-render"]), ("rotor", ["--no-preview"]),
           ("drivetrain", ["--no-render"]), ("fins", ["--no-render"]), ("coupon", ["--no-render"])]

VIEWS = [  # out/renders/<name>.png, viewer query (manifest out/web/manifest.json unless m= is given)
    ("hero_deployed", "view=hero&fold=1&zoom=1.45"),
    ("stowed", "view=hero&fold=0&packed=1&zoom=1.55"),          # PACKED: fins off, blanking strips in (T08)
    ("packed_kit", "m=../out/web/packed_kit.json&view=hero&zoom=1.4"),   # the packed egg with its fins beside it
    ("cutaway", "view=hero&fold=1&cutaway=1&zoom=1.45"),
    ("exploded", "view=55,24&fold=1&explode=1&zoom=1.15"),       # all three blades in view (T29)
    ("side", "view=side&fold=1&zoom=1.5"),
    ("front", "view=front&fold=1&zoom=1.4"),
]

DENSITY = {"PLA": 1.24, "PETG": 1.27}


def run(cmd: list[str], log: Path) -> tuple[int, float]:
    t0 = time.time()
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
    return rc, time.time() - t0


def step(label: str, cmd: list[str], log: Path, timing: list) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {label} ...", flush=True)
    rc, dt = run(cmd, log)
    timing.append((label, rc, dt, log))
    print(f"    {'ok' if rc == 0 else 'FAILED (rc %d)' % rc} in {dt:.0f} s, log {log.relative_to(REPO)}", flush=True)
    if rc != 0:
        print("".join(log.read_text().splitlines(True)[-25:]))
        sys.exit(f"build stopped: {label} failed")


def table(rows: list[list], head: list[str]) -> str:
    w = [max(len(str(x)) for x in col) for col in zip(head, *rows)]
    line = lambda r: "  ".join(str(x).ljust(n) for x, n in zip(r, w))
    return "\n".join([line(head), line(["-" * n for n in w])] + [line(r) for r in rows])


def summary(timing: list | None = None) -> None:
    f = REPO / "out" / "cad" / "checks" / "assembly.json"
    if not f.exists():
        sys.exit("no out/cad/checks/assembly.json: run the build first")
    A = json.loads(f.read_text())
    stl = {r["name"]: r for r in A["stl"]}
    plan = {r["name"]: r for r in A.get("print_plan", {}).get("parts", [])}
    rows = []
    for p in A["parts"]:
        s = stl.get(p["name"], {})
        pl = plan.get(p["name"], {})
        rows.append([p["module"], p["name"], "yes" if p["valid"] else "NO", p["n_solids"], f"{p['volume_cm3']:.2f}",
                     f"{s['stl_volume_cm3']:.2f}" + ("" if s.get("volume_match", True) else " (!)") if s else "-",
                     "bought" if p.get("bought") else ({True: "yes", False: "NO"}.get(s.get("watertight"), "-")),
                     "x".join(f"{v:.0f}" for v in s["print_size_mm"]) if s else "-",
                     pl.get("material", "-"), f"{pl['grams_est']:.1f}" if pl else "-"])
    print()
    print(table(rows, ["module", "part", "valid", "solids", "B-rep cm3", "STL cm3", "STL watertight", "print mm",
                       "mat", "g est"]))
    if any(not r.get("volume_match", True) for r in A["stl"]):
        print("(!) STL and B-rep volumes differ by more than 1 %: see out/cad/checks/assembly.json stl[]")
    print()
    irows = [[st, r["fold_deg"], r["n_parts"], r["n_pairs_near"], len(r["expected_overlaps"]), len(r["contacts"]),
              len(r["defects"])] for st, r in A["interference"].items()]
    print(table(irows, ["state", "fold", "parts", "pairs <3 mm", "intended overlaps", "contacts", "DEFECTS"]))
    print()
    krows = [["PASS" if k["ok"] else "FAIL", k["check"][:78], str(k["value"])[:60],
              "" if k.get("layout") is None else k["layout"]] for k in A["key_clearances"]]
    print(table(krows, ["", "key clearance (real solids)", "measured", "layout"]))
    sp = A["spin_sweep"]
    print(f"\nrotor spin sweep: {len(sp['pairs'])} rotor/static pairs within 3 mm, {len(sp['rubbing'])} rubbing")
    vf = A.get("viewer_fold_check")
    if vf:
        print(f"viewer fold check (browser vs CAD): max box difference {vf['max_abs_diff_mm']:.4f} mm "
              f"({'ok' if vf['ok'] else 'MISMATCH'})")
    if A.get("print_plan"):
        pp = A["print_plan"]
        print(f"print plan: {len(pp['plates'])} plates, {pp['total_grams_est']:.0f} g estimated "
              f"({', '.join(f'{k} {v:.0f} g' for k, v in pp['grams_by_material'].items())})")
    outs = [REPO / "out/cad/step/hydroegg_assembly_stowed.step", REPO / "out/cad/step/hydroegg_assembly_deployed.step",
            REPO / "out/web/hydroegg.glb", REPO / "out/web/manifest.json"] + \
           [REPO / "out/renders" / f"{n}.png" for n, _ in VIEWS]
    print("\noutputs:")
    for o in outs:
        print(f"  {'%7.2f MB' % (o.stat().st_size / 1e6) if o.exists() else 'MISSING   '}  {o.relative_to(REPO)}")
    if timing:
        print("\nsteps:")
        for label, rc, dt, _ in timing:
            print(f"  {dt:7.0f} s  {label}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--only-assembly", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--no-render", action="store_true")
    ap.add_argument("--previews", action="store_true")
    a = ap.parse_args()
    if a.summary:
        summary()
        return
    logs = REPO / "out" / "work" / "build"
    timing: list = []
    t0 = time.time()
    if a.fresh:
        for d in (REPO / "out/work/cad/shells/cache", REPO / "out/work/cad/integrate/cache"):
            shutil.rmtree(d, ignore_errors=True)
    if not a.only_assembly:
        for mod, quiet in MODULES:
            flags = [] if a.previews else list(quiet)
            if mod == "shells" and a.fresh:
                flags.append("--force")
            step(f"module {mod}", [CADRUN, PY, "-m", f"hydroegg.cad.{mod}", *flags], logs / f"{mod}.log", timing)
    step("assembly", [CADRUN, PY, "-m", "hydroegg.cad.assembly", "--force"], logs / "assembly.log", timing)
    if not a.no_render:
        jobs = [{"out": f"out/renders/{n}.png", "q": q} for n, q in VIEWS]
        jf = REPO / "out" / "work" / "cad" / "integrate" / "views.json"
        jf.parent.mkdir(parents=True, exist_ok=True)
        jf.write_text(json.dumps(jobs, indent=1))
        step("renders", [PY, str(REPO / "scripts" / "render.py"), "--batch", str(jf), "-W", "1600", "-H", "1000"],
             logs / "renders.log", timing)
    summary(timing)
    print(f"\nbuild done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
