"""Build the project web page's assets under site/ (the page itself, site/index.html, is hand-written).

  .venv/bin/python scripts/site.py                  # GLB, manifest, media, report, data blocks, number check
  .venv/bin/python scripts/site.py --stills         # also re-render the hero posters and the four fold stills
                                                    # (headless Chrome via agent-browser, takes out/work/heavy.lock)
  .venv/bin/python scripts/site.py --check          # only check the numbers quoted in site/index.html
  .venv/bin/python scripts/site.py --links          # check that every outward link and the og:image answer 200
  .venv/bin/python scripts/site.py --artifact DIR   # a copy of site/ for a Claude Artifact publish (see below)

What it writes (every path relative to site/, so the page works from GitHub Pages or any static host):
  assets/ovo.glb          out/web/ovo.glb without the per-mesh 'extras', meshopt-compressed and quantized with
                          @gltf-transform/cli (node names kept: the page finds parts by name)
  assets/manifest.json    out/web/manifest.json (parts, roles, pivots, explode vectors, materials), minified
  media/fold.mp4          docs/media/fold.mp4 when the media step has made it  } local preview copies, git-ignored:
  media/og_card.png       docs/media/og_card.png when present (social card)     } the pages workflow copies them
  report/ovo-report.pdf   docs/report/ovo-report.pdf when the report step has made it }  from docs/ at deploy time
  assets/poster-*.webp, assets/fold-*.webp   with --stills only: transparent stills drawn by the page's own
                          renderer (window.__ovo.capture), so the posters match the live 3D view exactly; re-run
                          --stills after a CAD change
  docs/media/fold_kinematics_{light,dark}.png   with --stills only: the page's fold chart (traced from the GLB)
                          as README images, light and dark
and it rewrites three generated blocks inside index.html: the chart data (<!-- ovo:data -->, including the fold
chart's blade outlines traced from out/web/ovo.glb), the print plates table (<!-- ovo:plates -->) and the fold media
(<!-- ovo:fold-media -->: the video if there is one). It ends with a freshness check: the media, the report and the
site assets must come from the CAD build that out/cad/checks/assembly.json describes, in that order.

--artifact DIR: the Artifact host wraps a page in its own <html>/<head>/<body>, and its CSP may refuse
WebAssembly (the meshopt decoder). The copy in DIR has the page without those wrapper tags and carries
assets/ovo-q.glb.gz, a quantize-only GLB (no decoder) that the page falls back to through DecompressionStream.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import functools
import gzip
import http.server
import json
import os
import re
import shutil
import socketserver
import struct
import subprocess
import sys
import threading
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SITE = REPO / "site"
WORK = REPO / "out/work/deliverables/site"
GLTF_CLI = "@gltf-transform/cli@4.5.1"
AB = shutil.which("agent-browser")
if AB is None:
    raise SystemExit("agent-browser not on PATH: npm install -g agent-browser && agent-browser install")


def load(rel: str):
    return json.loads((REPO / rel).read_text())


# --------------------------------------------------------------------------- GLB + manifest

def strip_extras(src: Path, dst: Path) -> None:
    """Copy a GLB without the per-mesh / per-node 'extras' (trimesh component lists, 1.7 MB of JSON)."""
    d = src.read_bytes()
    jlen = struct.unpack("<I", d[12:16])[0]
    j = json.loads(d[20:20 + jlen])
    rest = d[20 + jlen:]
    for key in ("meshes", "nodes"):
        for item in j.get(key, []):
            item.pop("extras", None)
    js = json.dumps(j, separators=(",", ":")).encode()
    js += b" " * ((4 - len(js) % 4) % 4)
    head = struct.pack("<III", 0x46546C67, 2, 12 + 8 + len(js) + len(rest)) + struct.pack("<II", len(js), 0x4E4F534A)
    dst.write_bytes(head + js + rest)


def gltf(*args: str) -> None:
    subprocess.run(["npx", "--yes", GLTF_CLI, *args], check=True, cwd=WORK, capture_output=True, text=True)


def build_glb() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (SITE / "assets").mkdir(parents=True, exist_ok=True)
    src = REPO / "out/web/ovo.glb"
    clean = WORK / "ovo_clean.glb"
    strip_extras(src, clean)
    out = SITE / "assets/ovo.glb"
    gltf("meshopt", str(clean), str(out), "--level", "high")
    man = load("out/web/manifest.json")
    man["glb"] = "ovo.glb"
    man.pop("cache", None)
    man["cad_key"] = load("out/cad/checks/assembly.json").get("cache_key")
    (SITE / "assets/manifest.json").write_text(json.dumps(man, separators=(",", ":")))
    print(f"assets/ovo.glb {out.stat().st_size / 1e6:.2f} MB (from {src.stat().st_size / 1e6:.2f} MB)")


# --------------------------------------------------------------------------- media + report

def copy_media() -> dict:
    found = {}
    media = SITE / "media"
    media.mkdir(parents=True, exist_ok=True)
    (media / "fold-kinematics.png").unlink(missing_ok=True)    # replaced by the page's fold chart
    for name in ("fold.mp4", "fold_dark.mp4", "og_card.png"):
        src = REPO / "docs/media" / name
        dst = media / name
        if src.exists():
            shutil.copy2(src, dst)
            found[name] = dst.stat().st_size
        elif dst.exists():
            dst.unlink()
    rep = REPO / "docs/report/ovo-report.pdf"
    if rep.exists():
        (SITE / "report").mkdir(exist_ok=True)
        shutil.copy2(rep, SITE / "report/ovo-report.pdf")
        found["ovo-report.pdf"] = rep.stat().st_size
    for k, v in found.items():
        print(f"copied {k} ({v / 1e6:.2f} MB)")
    for k in ("fold.mp4", "ovo-report.pdf"):
        if k not in found:
            print(f"note: {k} not found yet (run make media / make report, then this script again)")
    return found


# --------------------------------------------------------------------------- data blocks

def _rdp(pts, eps: float):
    """Ramer-Douglas-Peucker simplification of a polyline (N x 2)."""
    import numpy as np

    if len(pts) < 3:
        return pts
    a, b = pts[0], pts[-1]
    ab = b - a
    L = float(np.hypot(*ab))
    d = (np.hypot(*(pts - a).T) if L < 1e-12
         else np.abs(ab[0] * (pts[:, 1] - a[1]) - ab[1] * (pts[:, 0] - a[0])) / L)
    i = int(np.argmax(d))
    if d[i] > eps:
        return np.vstack([_rdp(pts[: i + 1], eps)[:-1], _rdp(pts[i:], eps)])
    return np.vstack([a, b])


def fold_chart() -> dict:
    """Blade 0 in its own meridional plane (x, z: the side view, blade 0 lies at theta 0) at fold 0 / 30 / 60 / 90,
    traced from the CAD mesh in out/web/ovo.glb: the stowed blade turned about its pivot (manifest), every triangle
    projected and rasterised at 0.2 mm, the outline (and the hinge bore) traced and simplified to 0.12 mm."""
    import math

    import contourpy
    import numpy as np
    import trimesh
    from PIL import Image, ImageDraw

    sys.path.insert(0, str(REPO))
    from ovo import params as P
    from ovo import profile

    man = load("out/web/manifest.json")
    piv = next(p_ for p_ in man["parts"] if p_["name"] == "blade_0")["pivot"]
    o, ax = np.array(piv["origin"], float), np.array(piv["axis"], float)
    assert np.allclose(ax, [0, 1, 0]), "blade 0 should fold about +Y (tangential at theta 0)"
    sc = trimesh.load(REPO / "out/web/ovo.glb", force="scene")
    T, gname = sc.graph["blade_0"]
    mesh = sc.geometry[gname].copy()
    mesh.apply_transform(T)                                   # the stowed pose
    px = 5.0
    blades, tips = {}, {}
    for fold in (0, 30, 60, 90):
        th = math.radians(piv["deployed_deg"] * (90 - fold) / 90)
        R = np.array([[math.cos(th), 0, math.sin(th)], [0, 1, 0], [-math.sin(th), 0, math.cos(th)]])
        V = ((mesh.vertices - o) @ R.T + o)[:, [0, 2]]
        lo = V.min(0) - 2
        size = np.ceil((V.max(0) + 2 - lo) * px).astype(int)
        im = Image.new("L", (int(size[0]), int(size[1])), 0)
        dr = ImageDraw.Draw(im)
        Q = (V - lo) * px
        for f_ in mesh.faces:
            dr.polygon([tuple(Q[i]) for i in f_], fill=255)
        lines = contourpy.contour_generator(z=np.asarray(im, float) / 255.0).lines(0.5)
        polys = []
        for ln in lines:
            if len(ln) < 8:
                continue
            xy = ln / px + lo
            xy = np.vstack([xy, xy[:1]])
            sp = _rdp(xy, 0.12)
            area = 0.5 * abs(np.dot(sp[:-1, 0], sp[1:, 1]) - np.dot(sp[1:, 0], sp[:-1, 1]))
            if area > 2.0:
                polys.append((area, sp[:-1]))
        polys.sort(key=lambda t: -t[0])
        blades[str(fold)] = [np.round(pp, 1).tolist() for _, pp in polys]
        outer = polys[0][1]
        k = int(np.argmax(np.hypot(*(outer - o[[0, 2]]).T)))
        tips[str(fold)] = np.round(outer[k], 1).tolist()
    xt = P.SPINNER_CAP["x_tip_built"]
    # cosine spacing on the elliptic nose (vertical tangent at x = 0), then even steps to the static tail end
    xn = P.HULL_NOSE_A * (1 - np.cos(np.linspace(0, math.pi / 2, 36)))
    xs = np.r_[xn, np.linspace(P.HULL_NOSE_A, P.STATIC_TAIL_X, 70)[1:]]
    xh = np.r_[np.linspace(P.STATIC_TAIL_X, xt, 40)]
    return {
        "source": "out/web/ovo.glb (blade_0 and its pivot), ovo/profile.py (hull), ovo/params.py",
        "hull": np.round(np.c_[xs, profile.radius(xs)], 2).tolist(),
        "hub": np.round(np.c_[xh, profile.radius(np.minimum(xh, P.HULL_L))], 2).tolist(),
        "blades": blades, "tips": tips, "hinge": [round(float(o[0]), 2), round(float(o[2]), 2)],
        "static_tail": P.STATIC_TAIL_X, "x_tip": round(xt, 2),
    }


def chart_data() -> dict:
    bd = load("data/blade_design.json")
    hw = load("data/hardware.json")
    ct = bd["cp_tsr"]
    pk = next(r for r in bd["performance"]["mppt"] if abs(r["flow_m_s"] - 1.5) < 1e-9)
    gens = []
    for c in hw["generator"]["candidates"]:
        if c.get("in_envelope_D56_L35") is False:
            continue                                   # the two out-of-envelope references
        rows = []
        for s in c["study"]:
            b = s["coupled_best"]["bridge_1.5_m_s"]
            rows.append({"kv": s["KV_rpm_per_V"], "r": s["R_ll_ohm"], "p": b["P_el_dc_W"], "v": b["V_dc_V"],
                         "i": b["I_A"], "tsr": b["tsr"], "rpm": b["rpm"]})
        gens.append({"name": c["name"].split(" (")[0], "cls": c["cls"],
                     "envelope": c.get("in_envelope_D56_L35"), "readings": rows})
    gl = next(c for c in hw["generator"]["candidates"] if c["id"] == "cubemars_gl40_kv70")["study"][0]
    td = gl["coupled_best_time_domain_bridge"]["td_bridge_1.5_m_s"]
    act = gl["coupled_best"]["active_1.5_m_s"]
    chain = [{"key": r["key"], "label": r["label"], "cp": r["cp_design_tsr"]} for r in bd["constraint_chain"]
             if r["key"] in ("optimum", "te", "linear", "chord", "twist", "final")]
    return {
        "source": "data/blade_design.json (cp_tsr, performance.mppt), data/hardware.json (generator.candidates)",
        "cp": {"flow": ct["flow_m_s"], "tsr": ct["tsr"], "final": ct["cp_final"], "optimum": ct["cp_bem_optimum"],
               "peak": {"tsr": round(pk["tsr"], 3), "cp": round(pk["cp"], 4)},
               "gen_tsr": td["tsr"], "design_tsr": bd["rotor"]["design_tsr"]},
        "gen": gens,
        "gl40": {"td_W": td["P_el_dc_W"], "td_V": td["V_dc_V"], "td_A": td["I_A"], "td_tsr": td["tsr"],
                 "active_W": act["P_el_dc_W"], "active_tsr": act["tsr"]},
        "chain": chain,
        "fold": fold_chart(),
    }


OPTIONAL_PLATES = {"4", "6"}   # alternatives outside the display build totals (cutaway shells, v2 front shell)


def plates_rows() -> str:
    pp = load("out/cad/checks/assembly.json")["print_plan"]
    rows = []
    for pl in pp["plates"]:
        label = pl["plate"]
        num, _, rest = label.partition(" ")
        title, _, mat = rest.partition(" (")
        mat = mat.rstrip(")")
        names = [p["name"] for p in pl["parts"] if not p["name"].startswith("tolerance_coupon_key")]
        keys = len(pl["parts"]) - len(names)
        parts = ", ".join(dict.fromkeys(names)) + (f" and {keys} fit keys" if keys else "")
        cls = ' class="opt"' if num in OPTIONAL_PLATES else ""
        rows.append(
            f'<tr{cls}><td class="pl">{num}</td><td class="job">{title}</td><td class="mat">{mat}</td>'
            f'<td class="parts">{parts}</td><td class="num h">{pl["max_height_mm"]:.1f}</td>'
            f'<td class="num g">{pl["grams_est"]:.1f}</td></tr>')
    return "\n".join(rows)


FOLD_VIDEO = """<figure class="fold-video">
  <video id="fold-video" src="media/fold.mp4" data-light="media/fold.mp4" data-dark="media/fold_dark.mp4" muted loop playsinline preload="metadata" aria-describedby="fold-video-cap"></video>
  <figcaption id="fold-video-cap">Animation rendered from the CAD assembly: the stowed egg unfolds its three blades and the rotor spins.</figcaption>
</figure>"""


def replace_block(html: str, name: str, body: str) -> str:
    pat = re.compile(rf"(<!-- ovo:{name} -->)(.*?)(<!-- /ovo:{name} -->)", re.DOTALL)
    if not pat.search(html):
        raise SystemExit(f"index.html has no <!-- ovo:{name} --> block")
    return pat.sub(lambda m: f"{m.group(1)}\n{body}\n{m.group(3)}", html)


def write_blocks(found: dict) -> None:
    page = SITE / "index.html"
    html = page.read_text()
    data = json.dumps(chart_data(), separators=(",", ":"))
    html = replace_block(html, "data", f'<script type="application/json" id="ovo-data">{data}</script>')
    html = replace_block(html, "plates", plates_rows())
    html = replace_block(html, "fold-media", FOLD_VIDEO if (SITE / "media/fold.mp4").exists() else "")
    page.write_text(html)


# --------------------------------------------------------------------------- number check

def fmt(x: float, nd: int) -> str:
    return f"{x:.{nd}f}"


def checks() -> list[tuple[str, str]]:
    """(claim, the exact text that must appear in index.html), computed from the sources."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("params", REPO / "ovo/params.py")
    P = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(P)
    bd = load("data/blade_design.json")
    a = load("out/cad/checks/assembly.json")
    bom = load("data/bom.json")["totals"]
    hw = load("data/hardware.json")
    pk = next(r for r in bd["performance"]["mppt"] if abs(r["flow_m_s"] - 1.5) < 1e-9)
    gl = next(c for c in hw["generator"]["candidates"] if c["id"] == "cubemars_gl40_kv70")["study"][0]
    td = gl["coupled_best_time_domain_bridge"]["td_bridge_1.5_m_s"]
    act = gl["coupled_best"]["active_1.5_m_s"]
    chain = {r["key"]: r["cp_design_tsr"] for r in bd["constraint_chain"]}
    kc = {k["check"]: k for k in a["key_clearances"]}
    pp = a["print_plan"]
    L = P.SPINNER_CAP["x_tip_built"]
    out = [
        ("hull diameter", f"{P.HULL_D:.0f}"),
        ("built egg length", f"{L:.0f} mm"),
        ("L/D", f"L/D {L / P.HULL_D:.2f}"),
        ("rotor diameter", f"{2 * P.ROTOR_R_TIP_TARGET:.0f} mm"),
        ("Cp peak", f"Cp {pk['cp']:.3f}"),
        ("TSR at Cp peak", f"TSR {pk['tsr']:.1f}"),
        ("shaft power at the peak, 1.5 m/s", f"{pk['power_W']:.1f} W"),
        ("rpm at the peak", f"{pk['rpm']:.0f} rpm"),
        ("Cp optimum (design TSR)", f"{chain['optimum']:.3f}"),
        ("constraint cost %", f"{100 * (1 - chain['final'] / chain['optimum']):.1f} %"),
        ("Cp final at the design TSR", f"to {chain['final']:.3f}"),
        ("chord root", f"{bd['laws']['chord_mm']['root']} mm"),
        ("twist root", f"{bd['laws']['twist_deg']['root']}"),
        ("GL40 TD bridge W", f"{td['P_el_dc_W']:.2f} W"),
        ("GL40 TD bridge V", f"{td['V_dc_V']:.2f} V"),
        ("GL40 TD TSR", f"TSR {td['tsr']:.2f}"),
        ("GL40 active W", f"{act['P_el_dc_W']:.2f} W"),
        ("GL40 freewheel at 2.5 m/s", f"{gl['Voc_freewheel_runaway']['2.5_m_s']['Voc_ll_peak_V']:.1f} V"),
        ("candidates screened", f"{len(hw['generator']['candidates'])} "),
        ("stowed span to pocket", kc["stowed aero span (x < 100) to the pocket floors and walls, 3 blades"]["value"].split()[0] + " mm"),
        ("stowed root to tail", kc["stowed root (x > 100, more than 8.25 mm from the detent boss face centre) to the static tail"]["value"].split()[0] + " mm"),
        ("hub gap", kc["hub gap: static tail end (aft shell) to hub front"]["value"]),
        ("key clearance checks", f"{len(a['key_clearances'])} key clearances"),
        ("interference states", f"{len(a['interference'])} states"),
        ("STL files (distinct; the blades share one)", f"{len({x['stl'] for x in a['stl']})} print STLs"),
        ("spin sweep pairs", f"{len(a['spin_sweep']['pairs'])} rotor"),
        ("viewer fold check", f"{a['viewer_fold_check']['max_abs_diff_mm']:.3f} mm"),
        ("filament total", f"{pp['total_grams_est']:.0f} g"),
        ("PLA grams", f"PLA {pp['grams_by_material']['PLA']:.0f} g"),
        ("PETG grams", f"PETG {pp['grams_by_material']['PETG']:.0f} g"),
        ("BOM display purchase", f"&euro;{bom['display_build_purchase_eur']:.2f}"),
        ("BOM display pro-rata", f"&euro;{bom['display_build_known_eur']:.2f}"),
        ("BOM v2 purchase", f"&euro;{bom['v2_additional_purchase_eur']:.2f}"),
        ("unpriced BOM lines", f"{bom['unpriced_items']} lines"),
        ("fin chord x span", f"{P.FIN['chord']:.0f} &times; {P.FIN['span']:.0f} mm"),
        ("fin incidence", f"{P.FIN['incidence_deg']}&deg;"),
        ("design flow", f"{P.DESIGN_FLOW_V} m/s"),
        ("hinge radius", f"r {P.HINGE_R:.0f} mm"),
        ("Glauert check", f"{bd['glauert_verification']['cp_bem']:.4f}"),
        ("assembly STEP parts", f"{a['interference']['stowed']['n_parts']} named parts"),
        ("as-built priced lines, pro rata (report)", f"&euro;{asbuilt_cost()} pro rata"),
        ("GL40 active TSR", f"{act['tsr']:.2f} with active"),
        ("stowed blade to its detent boss",
         kc["stowed blade to its detent boss (blade within 8.25 mm of the boss face centre)"]["value"].split()[0] + " mm"),
        ("printed parts in the display build", f"{display_printed(pp)} printed parts"),
        ("assembly steps", f"{len(a['assembly_steps'])} steps"),
    ]
    return out


def asbuilt_cost() -> str:
    """The report's as-built priced total (scripts/report.py writes it to docs/report/data.typ)."""
    m = re.search(r'bom_asbuilt_known: "([0-9.]+)"', (REPO / "docs/report/data.typ").read_text())
    return m[1] if m else "?"


def display_printed(pp: dict) -> int:
    """Printed parts of one display build: the plates without the coupons and the alternatives (cradles apart)."""
    return sum(len(pl["parts"]) for pl in pp["plates"]
               if pl["plate"].split(" ", 1)[0] not in OPTIONAL_PLATES | {"0"})


def run_check() -> bool:
    html = (SITE / "index.html").read_text()
    bad = [(c, s) for c, s in checks() if s not in html]
    for c, s in bad:
        print(f"CHECK FAIL: {c}: '{s}' not in site/index.html")
    if "\u2014" in html:   # the em dash, written as an escape
        print("CHECK FAIL: em dash in site/index.html")
        bad.append(("em dash", ""))
    print(f"number check: {len(checks()) - len([b for b in bad if b[0] != 'em dash'])} / {len(checks())} ok")
    return not bad


# --------------------------------------------------------------------------- stills (posters, fold strip)

@contextlib.contextmanager
def heavy_lock():
    lock = REPO / "out/work/heavy.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        yield


@contextlib.contextmanager
def serve(root: Path):
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    handler.log_message = lambda *a, **k: None
    with socketserver.TCPServer(("127.0.0.1", 0), handler) as httpd:
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{httpd.server_address[1]}"
        finally:
            httpd.shutdown()


STILLS = [
    # name, w, h, state, fit state (all four fold stills share the deployed fit so they keep one scale), view,
    # extra capture options (the strip is framed tighter and drawn without the ground shadow)
    ("poster-wide.png", 1600, 1200, {"fold": 1}, None, "hero", {}),
    ("poster-square.png", 1000, 1000, {"fold": 1}, None, "hero", {}),
    ("fold-90.png", 720, 540, {"fold": 0}, {"fold": 1}, "strip", {"fill": 0.94, "ground": False}),
    ("fold-60.png", 720, 540, {"fold": 1 / 3}, {"fold": 1}, "strip", {"fill": 0.94, "ground": False}),
    ("fold-30.png", 720, 540, {"fold": 2 / 3}, {"fold": 1}, "strip", {"fill": 0.94, "ground": False}),
    ("fold-0.png", 720, 540, {"fold": 1}, {"fold": 1}, "strip", {"fill": 0.94, "ground": False}),
]


def render_stills() -> None:
    import base64

    session = f"ovo-site-{os.getpid()}"

    def ab(*args: str) -> str:
        r = subprocess.run([AB, "--session", session, *args], capture_output=True, text=True, timeout=240, check=False,
                           env={**os.environ, "AGENT_BROWSER_SESSION": session})
        if r.returncode != 0:
            raise RuntimeError(f"agent-browser {args[0]} failed: {r.stderr or r.stdout}")
        return r.stdout.strip()

    with heavy_lock(), serve(SITE) as base:
        try:
            ab("set", "viewport", "1200", "900")
            ab("open", f"{base}/index.html")
            ab("eval", "window.__ovoBoot && window.__ovoBoot()")
            ab("wait", "--fn", "window.__ovo && window.__ovo.ready === true")
            for name, w, h, st, fit, view, opts in STILLS:
                args = {"w": w, "h": h, "state": st, "fit": fit, "view": view, **opts}
                js = f"window.__ovo.capture({json.dumps(args)})"
                out = ab("eval", js).strip().strip('"')
                b64 = out.split(",", 1)[1]
                (SITE / "assets" / name).write_bytes(base64.b64decode(b64))
                print(f"assets/{name}")
            fold_chart_pngs(ab)
        finally:
            subprocess.run([AB, "--session", session, "close"], capture_output=True, check=False)
    to_webp([SITE / "assets" / n for n, *_ in STILLS])


def fold_chart_pngs(ab) -> None:
    """The page's fold chart as README images (docs/media/fold_kinematics_{light,dark}.png): drawn at 900 CSS px
    wide with device pixel ratio 2, so its 14 px labels stay legible across GitHub's README column."""
    import time

    from PIL import Image

    ab("set", "viewport", "1240", "1000", "2")
    ab("eval", "(() => { const h = document.querySelector('#chart-fold'); h.style.width = '900px';"
               " h.style.padding = '18px 18px 8px'; h.style.background = 'var(--bg)'; h.style.boxSizing = 'content-box';"
               " h.scrollIntoView({block: 'center'}); return 1; })()")
    time.sleep(0.8)
    for theme in ("light", "dark"):
        ab("eval", f"document.documentElement.dataset.theme = '{theme}'; 1")
        time.sleep(0.6)
        out = REPO / "docs/media" / f"fold_kinematics_{theme}.png"
        # element screenshots clip the wrong region here: shoot the viewport and crop the box (CSS px x dpr 2)
        box = json.loads(json.loads(ab("eval", "JSON.stringify((() => { const r = document.querySelector('#chart-fold')"
                                              ".getBoundingClientRect(); return [r.left, r.top, r.width, r.height]; })())")))
        shot = WORK / f"fold_chart_{theme}.png"
        ab("screenshot", str(shot))
        full = Image.open(shot)
        k = full.width / 1240
        im = full.crop(tuple(round(v * k) for v in (box[0], box[1], box[0] + box[2], box[1] + box[3])))
        im.save(out, optimize=True)
        print(f"docs/media/{out.name}: {im.size[0]} x {im.size[1]}, {out.stat().st_size / 1e3:.0f} KB")


def to_webp(paths: list[Path], quality: int = 90) -> None:
    """Convert the transparent PNG stills to WebP (alpha kept, about a fifth of the PNG size) and drop the PNGs."""
    from PIL import Image

    for p in paths:
        out = p.with_suffix(".webp")
        Image.open(p).convert("RGBA").save(out, "WEBP", quality=quality, alpha_quality=100, method=6)
        p.unlink()
        print(f"{out.name}: {out.stat().st_size / 1e3:.0f} KB")


# --------------------------------------------------------------------------- freshness (H1, H2, P4)

def freshness() -> bool:
    """The media, the report and the site assets must all come from the CAD build in out/cad/checks/assembly.json,
    and be made in that order (media -> report -> site): a stale output is reported here, not shipped silently."""
    cad = load("out/cad/checks/assembly.json").get("cache_key")
    warn = []
    mj = REPO / "docs/media/media.json"
    if mj.exists() and json.loads(mj.read_text()).get("cad_key") != cad:
        warn.append(f"docs/media was rendered from CAD build {json.loads(mj.read_text()).get('cad_key')}, not {cad}: make media")
    dt = REPO / "docs/report/data.typ"
    m = re.search(r'cad_key: "([^"]*)"', dt.read_text()) if dt.exists() else None
    if not m or m[1] != cad:
        warn.append(f"the report's numbers come from CAD build {m[1] if m else '?'}, not {cad}: make report")
    pdf = REPO / "docs/report/ovo-report.pdf"
    media = [x for x in (REPO / "docs/media").glob("*") if x.suffix in (".png", ".mp4")]
    newer = [x.name for x in media if pdf.exists() and x.stat().st_mtime > pdf.stat().st_mtime
             and not x.name.startswith("fold_kinematics")]
    if newer:
        warn.append(f"docs/media is newer than the report PDF ({', '.join(sorted(newer)[:4])}...): make report")
    sm = SITE / "assets/manifest.json"
    if sm.exists() and json.loads(sm.read_text()).get("cad_key") != cad:
        warn.append("site/assets/manifest.json comes from another CAD build: make site")
    glb = REPO / "out/web/ovo.glb"
    poster = SITE / "assets/poster-wide.webp"
    if poster.exists() and glb.exists() and poster.stat().st_mtime < glb.stat().st_mtime:
        warn.append("the site posters are older than out/web/ovo.glb: scripts/site.py --stills")
    for a, b in (("docs/media/fold.mp4", "site/media/fold.mp4"), ("docs/media/fold_dark.mp4", "site/media/fold_dark.mp4"),
                 ("docs/media/og_card.png", "site/media/og_card.png"),
                 ("docs/report/ovo-report.pdf", "site/report/ovo-report.pdf")):
        if (REPO / a).exists() and (REPO / b).exists() and (REPO / a).read_bytes() != (REPO / b).read_bytes():
            warn.append(f"{b} differs from {a}: make site")
    for w in warn:
        print(f"FRESHNESS: {w}")
    print(f"freshness: {'ok' if not warn else f'{len(warn)} stale output(s)'} (CAD build {cad})")
    return not warn


# --------------------------------------------------------------------------- outward links (P1, H10)

PROJECT_URL = re.compile(r"https://(?:github\.com/lucaramirezo/ovo|lucaramirezo\.github\.io/ovo)[^\s\"'<>)]*")


def project_urls() -> list[str]:
    urls = []
    for f in ("site/index.html", "README.md", "pyproject.toml"):
        urls += PROJECT_URL.findall((REPO / f).read_text())
    return list(dict.fromkeys(u.rstrip(".,") for u in urls))


def url_status(url: str) -> int:
    import urllib.error
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "ovo-site-check"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            return r.status
    except urllib.error.HTTPError as e:
        return e.code
    except Exception:  # noqa: BLE001  (no network, DNS, timeout)
        return 0


def check_links() -> bool:
    """Every project link the page, the README and pyproject.toml point at must answer 200 before the page or the CV
    link is shared: the repository, its release, the Pages site and its og:image."""
    bad = []
    for u in project_urls():
        st = url_status(u)
        print(f"{st or 'ERR':>4}  {u}")
        if st != 200:
            bad.append(u)
    print(f"links: {'all answer 200' if not bad else f'{len(bad)} do not answer 200 yet (push the repo, tag a v* release, deploy Pages)'}")
    return not bad


# --------------------------------------------------------------------------- artifact copy

def artifact_copy(dest: Path) -> None:
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(SITE, dest)
    html = (dest / "index.html").read_text()
    html = re.sub(r"^\s*<!doctype html>\s*", "", html, flags=re.IGNORECASE)
    for tag in (r"<html[^>]*>", r"</html>", r"<head>", r"</head>", r"<body[^>]*>", r"</body>"):
        html = re.sub(tag, "", html, count=1)
    # project links that do not answer yet (repository not pushed, no release) are dropped from the copy
    for u in project_urls():
        if u.endswith(".png") or url_status(u) == 200:
            continue
        # buttons and file cards go; an inline link keeps its text
        pat = re.compile(rf'<a\b([^>]*)href="{re.escape(u)}"([^>]*)>(.*?)</a>(\s*)', re.DOTALL)
        html, k = pat.subn(lambda m: "" if re.search(r'class="(btn|file)', m[1] + m[2]) else m[3] + m[4], html)
        if k:
            print(f"artifact copy: dropped or unlinked {k} link(s) to {u} (not live yet)")
    (dest / "index.html").write_text(html.strip() + "\n")
    WORK.mkdir(parents=True, exist_ok=True)
    clean = WORK / "ovo_clean.glb"
    strip_extras(REPO / "out/web/ovo.glb", clean)
    q = WORK / "ovo_quant.glb"
    gltf("quantize", str(clean), str(q))
    with open(q, "rb") as fi, gzip.open(dest / "assets/ovo-q.glb.gz", "wb", compresslevel=9) as fo:
        shutil.copyfileobj(fi, fo)
    print(f"artifact copy in {dest} (fallback GLB {(dest / 'assets/ovo-q.glb.gz').stat().st_size / 1e6:.2f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stills", action="store_true", help="re-render the posters and fold stills")
    ap.add_argument("--check", action="store_true", help="only check the numbers in index.html")
    ap.add_argument("--artifact", metavar="DIR", help="write an Artifact-ready copy of site/ to DIR")
    ap.add_argument("--links", action="store_true", help="check the outward project links (needs network)")
    a = ap.parse_args()
    if a.check:
        sys.exit(0 if run_check() else 1)
    if a.links:
        sys.exit(0 if check_links() else 1)
    if a.artifact:
        artifact_copy(Path(a.artifact).resolve())
        return
    build_glb()
    found = copy_media()
    write_blocks(found)
    if a.stills:
        render_stills()
    ok = run_check()
    ok = freshness() and ok
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
