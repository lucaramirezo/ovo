"""Regenerate every Ovo CAD output from scratch with one command.

  .venv/bin/python scripts/build.py                 # part modules, assembly, renders, summary table
  .venv/bin/python scripts/build.py --fresh         # also drop the shells and assembly BREP caches and every
                                                    # output file that the current parts no longer produce
  .venv/bin/python scripts/build.py --only-assembly # reuse the part modules' outputs (out/cad/*)
  .venv/bin/python scripts/build.py --summary       # only print the table from the existing checks
  .venv/bin/python scripts/build.py --no-render     # skip out/renders
  .venv/bin/python scripts/build.py --previews      # also each module's own review renders

Steps (each geometry job goes through scripts/cadrun.sh, at most 2 OCC jobs machine-wide):
  1. part modules, in dependency order: shells, blade, rotor, drivetrain, fins, coupon. Each writes
     out/cad/step/<part>.step, out/cad/stl/<part>.stl (print orientation) and out/cad/checks/<module>.json
  2. ovo.cad.assembly --force: rebuilds the part cache from the modules, interference matrix
     (stowed / deployed / fold 30 / fold 60 / cutaway / packed), key clearances, rotor spin sweep, STL audit,
     out/cad/step/ovo_assembly_{stowed,deployed}.step, out/web/ovo.glb + manifest.json,
     viewer fold check, out/cad/checks/assembly.json
  3. renders through the viewer (1600 x 1000): out/renders/{hero_deployed,stowed,packed_kit,cutaway,exploded,side,front}.png
  4. --fresh only: prune out/cad/step, out/cad/stl, out/cad/checks, out/web and out/renders of files the
     build did not produce (parts removed since, test artefacts, files from before the rename)
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
    # fit=tight frames the projected parts (fill = share of the frame they span); shadow=1 frames the
    # ground shadow too (the level side and front views, whose shadow otherwise runs off the bottom);
    # ground=0 drops the ground shadow (the exploded cloud would cast a scattered one).
    # theme=studio (review V7): a mid warm-grey backdrop for this portfolio set, so the all-ivory egg, its
    # strips and the fins separate from it (the light theme stays the viewer's default for the docs).
    # Review V5: the hero and cutaway cameras sit at azimuth -40, elevation 22: the two fins read as a level
    # pair across the egg instead of the near fin filling the foreground.
    ("hero_deployed", "view=-40,22&fold=1&fit=tight&fill=0.86&theme=studio"),
    ("stowed", "view=hero&fold=0&packed=1&fit=tight&fill=0.8&theme=studio"),   # PACKED: fins off, strips in (T08)
    ("packed_kit", "m=../out/web/packed_kit.json&view=hero&fit=tight&fill=0.86&theme=studio"),   # loose fins
    ("cutaway", "view=-40,22&fold=1&cutaway=1&fit=tight&fill=0.86&theme=studio"),
    # review V3: one horizontal line in build order, blades radially up and down, the bottom half used
    ("exploded", "view=20,18&spin=30&fold=1&explode=1&fit=tight&fill=0.92&ground=0&theme=studio"),
    ("side", "view=side&fold=1&fit=tight&fill=0.88&shadow=1&theme=studio"),
    ("front", "view=front&fold=1&fit=tight&fill=0.9&shadow=1&theme=studio"),
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
    outs = [REPO / "out/cad/step/ovo_assembly_stowed.step", REPO / "out/cad/step/ovo_assembly_deployed.step",
            REPO / "out/web/ovo.glb", REPO / "out/web/manifest.json"] + \
           [REPO / "out/renders" / f"{n}.png" for n, _ in VIEWS]
    print("\noutputs:")
    for o in outs:
        print(f"  {'%7.2f MB' % (o.stat().st_size / 1e6) if o.exists() else 'MISSING   '}  {o.relative_to(REPO)}")
    if timing:
        print("\nsteps:")
        for label, rc, dt, _ in timing:
            print(f"  {dt:7.0f} s  {label}")


def produced_outputs() -> set[Path]:
    """Every output file the current build writes: one STEP per module part (print pose for the
    blade and coupon), one STL per printed part (the blades share blade.stl), the module checks,
    the assembly STEPs, the web scenes and the renders."""
    meta = sorted((REPO / "out/work/cad/integrate/cache").glob("*/meta.json"), key=lambda f: f.stat().st_mtime)
    keep: set[Path] = set()
    if meta:
        parts = json.loads(meta[-1].read_text())["parts"]
        for d in parts:
            if d.get("module") == "assembly":
                continue
            n = d["name"]
            stem = "blade" if n.startswith("blade_") and n != "blade_cradle" else n
            keep.add(REPO / "out/cad/step" / f"{stem}.step")
            if d.get("has_print"):
                keep.add(REPO / "out/cad/stl" / f"{stem}.stl")
    keep |= {REPO / "out/cad/checks" / f"{m}.json" for m, _ in MODULES} | {REPO / "out/cad/checks/assembly.json"}
    keep |= {REPO / "out/cad/step" / f"ovo_assembly_{s}.step" for s in ("stowed", "deployed")}
    keep |= {REPO / "out/web" / f for f in ("ovo.glb", "manifest.json", "packed_kit.glb", "packed_kit.json")}
    keep |= {REPO / "out/renders" / f"{n}.png" for n, _ in VIEWS}
    return keep


def prune_outputs() -> list[str]:
    """--fresh: delete files in the output folders that the build did not produce (parts removed since,
    test artefacts such as out/web/test.glb and manifest_test.json, names from before the rename)."""
    keep = produced_outputs()
    if not any(p.parent.name == "step" for p in keep):
        return []                                         # no part cache: never prune blindly
    removed = []
    for d in ("out/cad/step", "out/cad/stl", "out/cad/checks", "out/web", "out/renders"):
        for f in sorted((REPO / d).glob("*")):
            if f.is_file() and f not in keep:
                f.unlink()
                removed.append(str(f.relative_to(REPO)))
    return removed


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
            step(f"module {mod}", [CADRUN, PY, "-m", f"ovo.cad.{mod}", *flags], logs / f"{mod}.log", timing)
    step("assembly", [CADRUN, PY, "-m", "ovo.cad.assembly", "--force"], logs / "assembly.log", timing)
    if not a.no_render:
        jobs = [{"out": f"out/renders/{n}.png", "q": q} for n, q in VIEWS]
        jf = REPO / "out" / "work" / "cad" / "integrate" / "views.json"
        jf.parent.mkdir(parents=True, exist_ok=True)
        jf.write_text(json.dumps(jobs, indent=1))
        step("renders", [PY, str(REPO / "scripts" / "render.py"), "--batch", str(jf), "-W", "1600", "-H", "1000"],
             logs / "renders.log", timing)
    if a.fresh:
        removed = prune_outputs()
        print(f"\npruned {len(removed)} stale output file(s)" + "".join(f"\n  {r}" for r in removed))
    summary(timing)
    print(f"\nbuild done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
