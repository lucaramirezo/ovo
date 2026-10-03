"""Quick visual check of any part(s) while modeling.

  .venv/bin/python scripts/preview.py hydroegg.blade:blade_stowed hydroegg.hull:front_shell \
      --out out/work/preview_blade.png --q "view=iso"

Each target is module:function; the function is called with no arguments and must
return a build123d Shape, or a list of dicts in export_web.write_scene format.
Writes a throwaway GLB + manifest under out/work/preview/ and renders it through the
viewer (takes the shared heavy lock). Then open the PNG with your image reader.
"""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from hydroegg.export_web import write_scene  # noqa: E402

PALETTE = ["shell", "blade", "hub", "generator", "steel", "fin", "stand", "rubber"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("targets", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--q", default="view=iso")
    ap.add_argument("-W", type=int, default=1400)
    ap.add_argument("-H", type=int, default=900)
    a = ap.parse_args()

    parts = []
    for i, t in enumerate(a.targets):
        mod, _, fn = t.partition(":")
        res = getattr(importlib.import_module(mod), fn)()
        if isinstance(res, list):
            parts.extend(res)
        else:
            parts.append({"name": f"{fn}_{i}", "shape": res, "material": PALETTE[i % len(PALETTE)]})
    stem = Path(a.out).stem
    write_scene(parts, REPO / "out/work/preview", glb_name=f"{stem}.glb")
    q = f"m=../out/work/preview/{stem}.json&{a.q}"
    subprocess.run([sys.executable, str(REPO / "scripts/render.py"), a.out, q, "-W", str(a.W), "-H", str(a.H)],
                   check=True, cwd=REPO)


if __name__ == "__main__":
    main()
