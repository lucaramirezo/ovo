"""Docs media for Ovo: stills, close-ups, social card, technical view and animations, all rendered
from out/web/ovo.glb through viewer/index.html, into docs/media/.

  .venv/bin/python scripts/media.py                      # everything (about 15 min with a GPU)
  .venv/bin/python scripts/media.py --only hero,fold     # by name or group (see --list)
  .venv/bin/python scripts/media.py --list
  .venv/bin/python scripts/media.py --sheet out.png "view=-40,22&fold=1" "view=40,20&fold=1" ...
                                                          # contact sheet of viewer queries (look dev)
  MEDIA_PREVIEW=0,45,90 .venv/bin/python scripts/media.py --only fold   # a few animation frames, no encoding

Outputs (docs/media/, listed with their viewer queries in media.json): hero_light / hero_dark (deployed,
tail-side 3/4), stowed (packed: fins off, blanking strips in), cutaway, exploded, details_{hub, fin_rail,
nose, blade}, og_card (1200 x 630), technical (first-angle 3-view, dimensioned from ovo/params.py and the
model's part boxes), packed_kit (the packed egg with its two loose fins, for the report), fold.{mp4,gif}
(stowed -> unfold -> spin-up -> fade back to the start, motion-blurred, a seamless loop), fold_dark.mp4 (the same
on the `ink` backdrop, for the project page's dark theme) and exploded.mp4 (assemble / disassemble loop).
media.json records the CAD build (assembly cache key) the files come from and the label font used.

Fonts: the labels use Inter when it is installed (~/.local/share/fonts or /usr/share/fonts), else matplotlib's
bundled DejaVu Sans, which always exists; media.json says which one a run used.

One visual language: every render uses the `paper` backdrop (hero_dark: `ink`), the `media` light rig,
the viewer's 28 deg lens and its fixed world lights (see viewer/index.html, presentation-only params).

Rendering: Chrome through agent-browser (scripts/render.py helpers, its own session `ovo-media`),
under out/work/heavy.lock. Chrome is launched with ANGLE on the host GL (`--use-angle=gl-egl`; on WSL2
that is Mesa's D3D12 driver, about 0.1 s a frame); if that fails, or with --cpu, it falls back to
Chrome's default software GL (SwiftShader, about 5 s a frame at 1200 x 750). Frames go to
out/work/deliverables/media/ and are deleted after encoding (--keep keeps them).

Size budget for the committed files: PNG <= 600 KB (losslessly optimised, else 256-colour quantised
with dithering), GIF <= 8 MB, MP4 <= 6 MB. Every performance figure in this project is a model
estimate; the media show geometry only and carry no performance claims.
"""

from __future__ import annotations

import argparse
import functools
import io
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "docs/media"
WORK = REPO / "out/work/deliverables/media"
os.environ.setdefault("AGENT_BROWSER_SESSION", "ovo-media")
sys.path.insert(0, str(REPO / "scripts"))
import render as R  # noqa: E402  (ab, heavy_lock, serve)

GPU_ARGS = "--use-angle=gl-egl,--ignore-gpu-blocklist"
PNG_BUDGET = 600_000
GIF_BUDGET = 8_000_000
MP4_BUDGET = 6_000_000

LIGHT = "theme=paper&rig=media"
DARK = "theme=ink&rig=media"
PAPER_BG = (0xE3, 0xDF, 0xD8)       # viewer body.paper --bg0: the panels' and the drawing's paper
INK = (0x1D, 0x1D, 0x1B)
ORANGE = (0xE8, 0x62, 0x2C)         # manifest blade colour


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------------------------
# Browser
# --------------------------------------------------------------------------------------------
class Viewer:
    """One agent-browser session on viewer/index.html (served by render.serve)."""

    def __init__(self, base: str, gpu: bool = True):
        self.base, self.gpu, self.page = base, gpu, None
        self.renderer = "?"

    def _launch_env(self) -> None:
        if self.gpu:
            os.environ["AGENT_BROWSER_ARGS"] = GPU_ARGS
        else:
            os.environ.pop("AGENT_BROWSER_ARGS", None)

    def open(self, q: str, w: int, h: int) -> None:
        page = (q, w, h)
        if page == self.page:
            return
        self._launch_env()
        url = f"{self.base}/viewer/index.html?w={w}&h={h}&dpr=1&{q}"
        try:
            R.ab("set", "viewport", str(w), str(h))
            R.ab("open", url)
            R.ab("wait", "#ready")
        except RuntimeError as e:
            if not self.gpu:
                raise
            log(f"GPU launch failed ({str(e).strip()[:120]}); falling back to software GL")
            R.ab("close", check=False)
            self.gpu = False
            self._launch_env()
            R.ab("set", "viewport", str(w), str(h))
            R.ab("open", url)
            R.ab("wait", "#ready")
        if self.page is None:
            self.renderer = self.js("""(() => { const gl = document.querySelector('canvas').getContext('webgl2');
                const x = gl.getExtension('WEBGL_debug_renderer_info');
                return x ? gl.getParameter(x.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER); })()""")
            log(f"renderer: {self.renderer}")
        self.page = page

    @property
    def soft(self) -> bool:
        """Software GL (SwiftShader): about 50 x slower, so the animations drop to fewer, smaller sub-frames."""
        return not self.gpu or "swiftshader" in str(self.renderer).lower()

    def anim_open(self, q: str) -> tuple[int, int]:
        """Open an animation page: ANIM_RENDER on a GPU, the output size on software GL."""
        self.open(q, *ANIM_RENDER)
        if self.soft:
            log("software GL: animation frames at output size, at most 3 sub-frames each")
            self.open(q, *ANIM_OUT)
            return ANIM_OUT
        return ANIM_RENDER

    def js(self, expr: str):
        out = R.ab("eval", f"JSON.stringify({expr})")
        v = json.loads(out) if out else None
        return json.loads(v) if isinstance(v, str) else v

    def run(self, expr: str) -> None:
        R.ab("eval", f"{expr}; 0")

    def shot(self, path: Path) -> Image.Image:
        path.parent.mkdir(parents=True, exist_ok=True)
        R.ab("screenshot", str(path))
        return Image.open(path).convert("RGB")

    def close(self) -> None:
        R.ab("close", check=False)


# --------------------------------------------------------------------------------------------
# Image helpers
# --------------------------------------------------------------------------------------------
def save_png(im: Image.Image, path: Path, budget: int = PNG_BUDGET) -> int:
    """Optimised lossless PNG; over budget, a 256-colour PNG from ffmpeg's palettegen / paletteuse (sierra
    dithering: on these renders it is hard to tell from the lossless file at about half the size)."""
    import imageio_ffmpeg
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    data = buf.getvalue()
    if len(data) > budget:
        src, dst = WORK / "q_src.png", WORK / "q_dst.png"
        im.convert("RGB").save(src, compress_level=1)
        subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error", "-i", str(src), "-vf",
                        "split[a][b];[a]palettegen=max_colors=256:stats_mode=full:reserve_transparent=0[p];"
                        "[b][p]paletteuse=dither=sierra2_4a", "-pix_fmt", "pal8", str(dst)], check=True)
        buf = io.BytesIO()
        with Image.open(dst) as q:
            q.save(buf, "PNG", optimize=True)
        data = buf.getvalue()
        src.unlink()
        dst.unlink()
    path.write_bytes(data)
    return len(data)


def down(im: Image.Image, size: tuple[int, int]) -> Image.Image:
    return im if im.size == tuple(size) else im.resize(size, Image.Resampling.LANCZOS)


def contact_sheet(viewer: Viewer, queries: list[str], out: Path, w: int = 960, h: int = 600, cols: int = 3) -> None:
    tiles = []
    for i, q in enumerate(queries):
        viewer.open(q, w, h)
        im = viewer.shot(WORK / "sheet" / f"{i:02d}.png")
        d = ImageDraw.Draw(im)
        d.text((10, 8), q[:150], fill=(0, 0, 0))
        tiles.append(im)
    rows = math.ceil(len(tiles) / cols)
    sheet = Image.new("RGB", (cols * w, rows * h), "white")
    for i, t in enumerate(tiles):
        sheet.paste(t, ((i % cols) * w, (i // cols) * h))
    out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out)
    log(f"sheet {out} ({len(tiles)} tiles)")


# --------------------------------------------------------------------------------------------
# Fonts and labels
# --------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def font(size: int, weight: str = "Regular") -> ImageFont.FreeTypeFont:
    """Inter (variable font) when installed, else matplotlib's bundled DejaVu Sans (always there)."""
    for c in (Path.home() / ".local/share/fonts/Inter.ttf", Path("/usr/share/fonts/truetype/inter/Inter.ttf"),
              Path("/usr/share/fonts/opentype/inter/Inter.ttf")):
        if c.exists():
            f = ImageFont.truetype(str(c), size)
            try:
                f.set_variation_by_name(weight)
            except (OSError, ValueError):
                pass
            return f
    import matplotlib
    name = "DejaVuSans-Bold" if weight in ("Bold", "SemiBold", "Black") else "DejaVuSans"
    return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / f"fonts/ttf/{name}.ttf"), size)


def label(im: Image.Image, text: str, xy: tuple[int, int], size: int, color=INK, weight: str = "Medium") -> None:
    ImageDraw.Draw(im).text(xy, text, font=font(size, weight), fill=color)


# --------------------------------------------------------------------------------------------
# Stills
# --------------------------------------------------------------------------------------------
HERO_VIEW = "view=50,18"      # tail-side 3/4: all three blades, spinner, egg and a fin read at once
HERO = f"{HERO_VIEW}&fold=1&fit=tight&fill=0.86&shadow=1"
BLADE_ONLY = "only=blade_0,blade_magnet_heel_0,blade_magnet_stow_0"
NOSE_BOX = "fbox=-20,-40,-40,25,40,40"

# name -> (viewer query, render size, output size)
STILLS = {
    "hero_light": (f"{HERO}&{LIGHT}", (2400, 1500), (1600, 1000)),
    "hero_dark": (f"{HERO}&{DARK}", (2400, 1500), (1600, 1000)),
    # packed: blades in their pockets, fins off, blanking strips in (the egg as carried)
    "stowed": (f"view=-35,22&fold=0&packed=1&fit=tight&fill=0.8&shadow=1&{LIGHT}", (2400, 1500), (1600, 1000)),
    "cutaway": (f"view=-50,28&fold=1&cutaway=1&fit=tight&fill=0.86&shadow=1&{LIGHT}", (2400, 1500), (1600, 1000)),
    # one line in build order (assembly.EXPLODE_STATIONS), rotor turned 30 deg: one blade up, one down
    "exploded": (f"view=20,18&spin=30&fold=1&explode=1&fit=tight&fill=0.92&ground=0&{LIGHT}", (2400, 1500), (1600, 1000)),
    # hub with the spinner off: knuckles on their pins, heel stops, magnets; blades part-folded
    "details_hub": (f"view=40,32&fold=0.55&hide=spinner,spinner_cap,spinner_screw&internals=1&focus=hub&fit=tight"
                    f"&fill=0.62&{LIGHT}", (2400, 1500), (1600, 1000)),
    # the packed egg with its two loose fins (out/web/packed_kit.glb), same camera as `stowed`: the report pairs them
    "packed_kit": (f"m=../out/web/packed_kit.json&view=-35,22&fit=tight&fill=0.8&shadow=1&{LIGHT}", (2400, 1500), (1600, 1000)),
    # a fin part-way out of its recessed half-dovetail rail (the slot mouth and the floor magnet show)
    "details_fin_rail": (f"view=-62,16&fold=1&pexp=fin_1:0.35,fin_magnet_1:0.35&fbox=0,-115,-30,70,-15,20&fit=tight"
                         f"&fill=0.9&{LIGHT}", (2400, 1500), (1600, 1000)),
}
# two portrait panels side by side, with a caption each
PANELS = {
    "details_nose": [
        (f"view=-70,15&fold=1&{NOSE_BOX}&fit=tight&fill=0.8&{LIGHT}", "Display build: clean nose, tow eye"),
        (f"view=-70,15&fold=1&v2=1&{NOSE_BOX}&fit=tight&fill=0.8&{LIGHT}", "v2 front shell: M12 cable gland"),
    ],
    "details_blade": [
        (f"view=100,8&fold=1&{BLADE_ONLY}&fit=tight&fill=0.84&ground=0&{LIGHT}", "Blade planform, seen from downstream"),
        (f"view=150,60&fold=1&{BLADE_ONLY}&fit=tight&fill=0.84&ground=0&{LIGHT}", "From above and downstream: the twist"),
    ],
}
PANEL_RENDER, PANEL_OUT = (1200, 1500), (800, 1000)


def still(v: Viewer, name: str, manifest: dict) -> None:
    q, size, out = STILLS[name]
    v.open(q, *size)
    im = down(v.shot(WORK / "raw" / f"{name}.png"), out)
    n = save_png(im, OUT / f"{name}.png")
    manifest[f"{name}.png"] = {"viewer": q, "render_px": list(size), "px": list(out), "bytes": n}
    log(f"{name}.png {out[0]}x{out[1]} {n / 1e3:.0f} KB")


def panels(v: Viewer, name: str, manifest: dict) -> None:
    tiles = []
    for i, (q, cap) in enumerate(PANELS[name]):
        v.open(q, *PANEL_RENDER)
        tiles.append(down(v.shot(WORK / "raw" / f"{name}_{i}.png"), PANEL_OUT))
    pw, ph = PANEL_OUT
    im = Image.new("RGB", (pw * len(tiles), ph), PAPER_BG)
    for i, t in enumerate(tiles):
        im.paste(t, (i * pw, 0))
    d = ImageDraw.Draw(im)
    for i, (_, cap) in enumerate(PANELS[name]):
        label(im, cap, (i * pw + 36, ph - 64), 26)
        if i:
            d.line([(i * pw, 40), (i * pw, ph - 40)], fill=(0xC9, 0xC3, 0xB9), width=2)
    n = save_png(im, OUT / f"{name}.png")
    manifest[f"{name}.png"] = {"viewer": [q for q, _ in PANELS[name]], "render_px": list(PANEL_RENDER),
                               "px": list(im.size), "bytes": n}
    log(f"{name}.png {im.size[0]}x{im.size[1]} {n / 1e3:.0f} KB")


# --------------------------------------------------------------------------------------------
# Social card (1200 x 630)
# --------------------------------------------------------------------------------------------
CARD_Q = f"{HERO}&{LIGHT}&offset=0.2,0".replace("fill=0.86", "fill=0.84")
CARD_LINE = "A folding hydrokinetic charger,\ndesigned end to end in code-CAD."


def og_card(v: Viewer, manifest: dict) -> None:
    v.open(CARD_Q, 2400, 1260)
    im = down(v.shot(WORK / "raw" / "og_card.png"), (1200, 630))
    d = ImageDraw.Draw(im)
    x0 = 72
    d.rectangle([x0, 200, x0 + 56, 207], fill=ORANGE)
    d.text((x0 - 6, 222), "Ovo", font=font(132, "Bold"), fill=INK)
    d.multiline_text((x0, 392), CARD_LINE, font=font(30, "Regular"), fill=(0x3A, 0x39, 0x36), spacing=10)
    n = save_png(im, OUT / "og_card.png")
    manifest["og_card.png"] = {"viewer": CARD_Q, "render_px": [2400, 1260], "px": [1200, 630], "bytes": n,
                               "text": ["Ovo", CARD_LINE.replace("\n", " ")]}
    log(f"og_card.png 1200x630 {n / 1e3:.0f} KB")


# --------------------------------------------------------------------------------------------
# Technical view: first-angle 3-view of the deployed machine, dimensioned
# --------------------------------------------------------------------------------------------
KR = 6.0          # render px per mm (drawn at this scale, then halved: 3 px per mm in the PNG)
TECH_M = 6.0      # mm of paper around each view's geometry
TXT = 9.6         # label height in sheet mm: about 29 px in the 1700 px wide PNG, legible at README width


def _params():
    sys.path.insert(0, str(REPO))
    import ovo.params as P
    return P


class Sheet:
    """Drawing at KR px/mm; mm coordinates in, pixels out."""

    def __init__(self, w_mm: float, h_mm: float):
        self.im = Image.new("RGB", (round(w_mm * KR), round(h_mm * KR)), PAPER_BG)
        self.d = ImageDraw.Draw(self.im)
        self.dim_c, self.lw = (0x4A, 0x47, 0x42), 4

    def px(self, p):
        return (p[0] * KR, p[1] * KR)

    def line(self, a, b, w=None, c=None):
        self.d.line([self.px(a), self.px(b)], fill=c or self.dim_c, width=w or self.lw)

    def dashed(self, pts, dash=3.0, gap=2.0, w=None):
        """Polyline (mm) drawn dashed."""
        acc, on = 0.0, True
        for a, b in zip(pts, pts[1:]):
            seg = math.dist(a, b)
            t = 0.0
            while t < seg:
                run = (dash if on else gap) - acc
                t1 = min(seg, t + run)
                if on:
                    f0, f1 = t / seg, t1 / seg
                    self.line((a[0] + (b[0] - a[0]) * f0, a[1] + (b[1] - a[1]) * f0),
                              (a[0] + (b[0] - a[0]) * f1, a[1] + (b[1] - a[1]) * f1), w=w)
                acc = 0.0 if t1 - t >= run else acc + (t1 - t)
                if t1 - t >= run:
                    on = not on
                t = t1

    def arrow(self, tip, frm, size=4.2):
        ang = math.atan2(tip[1] - frm[1], tip[0] - frm[0])
        a1 = (tip[0] - size * math.cos(ang - 0.3), tip[1] - size * math.sin(ang - 0.3))
        a2 = (tip[0] - size * math.cos(ang + 0.3), tip[1] - size * math.sin(ang + 0.3))
        self.d.polygon([self.px(tip), self.px(a1), self.px(a2)], fill=self.dim_c)

    def text(self, p, s, size=TXT, weight="Medium", c=None, anchor="mm", bg=True):
        f = font(round(size * KR), weight)
        if bg:
            box = self.d.textbbox(self.px(p), s, font=f, anchor=anchor)
            pad = 0.8 * KR
            self.d.rectangle([box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad], fill=PAPER_BG)
        self.d.text(self.px(p), s, font=f, fill=c or self.dim_c, anchor=anchor)

    def vtext(self, p, s, size=TXT, weight="Medium", c=None):
        """Text along a vertical line, reading bottom to top, centred on p, on a paper patch."""
        f = font(round(size * KR), weight)
        box = self.d.textbbox((0, 0), s, font=f, anchor="lt")
        pad = round(0.8 * KR)
        tile = Image.new("RGBA", (box[2] + 2 * pad, box[3] + 2 * pad), PAPER_BG + (255,))
        ImageDraw.Draw(tile).text((pad, pad), s, font=f, fill=(c or self.dim_c) + (255,), anchor="lt")
        tile = tile.rotate(90, expand=True)
        x, y = self.px(p)
        self.im.paste(tile, (round(x - tile.width / 2), round(y - tile.height / 2)), tile)

    def dim(self, a, b, off, label, ext_from=(None, None), text_size=None, beside=False):
        """Linear dimension between sheet points a and b, offset `off` mm along the normal (left of a->b);
        extension lines start at ext_from points (default a and b)."""
        text_size = text_size or TXT
        dx, dy = b[0] - a[0], b[1] - a[1]
        L = math.hypot(dx, dy)
        nx, ny = -dy / L, dx / L
        a2 = (a[0] + nx * off, a[1] + ny * off)
        b2 = (b[0] + nx * off, b[1] + ny * off)
        for src, dst in ((ext_from[0] or a, a2), (ext_from[1] or b, b2)):
            gap, over = 1.5, 2.0
            vx, vy = dst[0] - src[0], dst[1] - src[1]
            vl = math.hypot(vx, vy) or 1.0
            self.line((src[0] + vx / vl * gap, src[1] + vy / vl * gap), (dst[0] + vx / vl * over, dst[1] + vy / vl * over), w=3)
        self.line(a2, b2)
        self.arrow(a2, b2)
        self.arrow(b2, a2)
        mid = ((a2[0] + b2[0]) / 2, (a2[1] + b2[1]) / 2)
        if abs(dx) < 1e-9 and not beside:   # vertical dimension: text along the line
            self.vtext(mid, label, size=text_size)
        elif beside:   # text outside the line instead of on it (side = the offset direction)
            sgn = 1 if off > 0 else -1
            self.text((mid[0] + nx * sgn * 3, mid[1] + ny * sgn * 3), label, size=text_size, bg=False,
                      anchor="lm" if nx * sgn > 0.5 else "rm" if nx * sgn < -0.5 else "mm")
        else:
            self.text(mid, label, size=text_size)


def technical(v: Viewer, manifest: dict) -> None:
    P = _params()
    base = "fold=1&ground=0&ortho=1&fit=tight&theme=paper&rig=media"
    v.open(base + "&view=0,0", 800, 500)
    boxes = [b for b in v.js("window.partBoxes()") if b["visible"]]
    lo = np.min([b["min"] for b in boxes], axis=0)
    hi = np.max([b["max"] for b in boxes], axis=0)
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    m = TECH_M
    # each view: (key, az, el, u range, v range, CAD point at the tile centre as f(uc, vc))
    views = {
        "side": (0, 0, (x0, x1), (z0, z1), lambda uc, vc: (uc, 0.0, vc)),
        "plan": (0, 89.9, (x0, x1), (y0, y1), lambda uc, vc: (uc, vc, 0.0)),
        "upstream": (-90, 0, (-y1, -y0), (z0, z1), lambda uc, vc: (0.0, -uc, vc)),
    }
    tiles = {}
    for key, (az, el, ur, vr, centre) in views.items():
        tw, th = (ur[1] - ur[0]) + 2 * m, (vr[1] - vr[0]) + 2 * m
        W, H = round(tw * KR), round(th * KR)
        uc, vc = (ur[0] + ur[1]) / 2, (vr[0] + vr[1]) / 2
        cx, cy, cz = centre(uc, vc)
        v.open(f"{base}&view={az},{el}", W, H)
        v.run("document.body.classList.add('flat')")
        v.run(f"window.setCamera({{c: [{cx}, {cz}, {-cy}], half: {th / 2}, az: {az}, el: {el}, r: 800}})")
        # check the mapping: CAD centre -> tile centre, +10 mm along u -> +10 KR px
        du = {"side": (10, 0, 0), "plan": (10, 0, 0), "upstream": (0, -10, 0)}[key]
        pa, pb = v.js(f"window.projectCad([[{cx}, {cy}, {cz}], [{cx + du[0]}, {cy + du[1]}, {cz + du[2]}]])")
        assert abs(pa[0] - W / 2) < 1.5 and abs(pa[1] - H / 2) < 1.5 and abs(pb[0] - pa[0] - 10 * KR) < 1.5, (key, pa, pb)
        tiles[key] = (v.shot(WORK / "raw" / f"technical_{key}.png"), ur, vr, tw, th)

    # sheet layout (mm)
    ML, MT, GX, GY, MR, MB = 50.0, 40.0, 56.0, 62.0, 52.0, 34.0
    sw, sh = tiles["side"][3], tiles["side"][4]
    lw, ph = tiles["upstream"][3], tiles["plan"][4]
    S = Sheet(ML + sw + GX + lw + MR, MT + sh + GY + ph + MB)
    org = {"side": (ML, MT), "upstream": (ML + sw + GX, MT), "plan": (ML, MT + sh + GY)}
    for key, (im, *_r) in tiles.items():
        S.im.paste(im, (round(org[key][0] * KR), round(org[key][1] * KR)))

    def to_sheet(key, u, vv):
        _, ur, vr, tw, th = tiles[key]
        ox, oy = org[key]
        return (ox + m + (u - ur[0]), oy + m + (vr[1] - vv))

    side = lambda x, z: to_sheet("side", x, z)          # noqa: E731
    plan = lambda x, y: to_sheet("plan", x, y)          # noqa: E731
    ups = lambda y, z: to_sheet("upstream", -y, z)      # noqa: E731

    tip = P.SPINNER_CAP["x_tip_built"]
    D = P.HULL_D
    swept = 2 * P.ROTOR_R_TIP_TARGET
    eye = x0                                    # tow eye: the front-most point of the model
    span = y1 - y0                              # over the fins
    assert abs(x1 - tip) < 0.2, (x1, tip)       # the model ends where params says the egg ends

    # side elevation: hull and overall length below, diameter at the left
    S.dim(side(0, z0), side(tip, z0), 14, f"{tip:.0f}", ext_from=(side(0, 0), side(tip, 0)))
    S.dim(side(eye, z0), side(tip, z0), 33, f"{tip - eye:.0f} with tow eye", ext_from=(side(eye, 0), side(tip, 0)))
    S.dim(side(eye, -D / 2), side(eye, D / 2), -18, f"Ø{D:.0f}",
          ext_from=(side(P.EQUATOR_X, -D / 2), side(P.EQUATOR_X, D / 2)))
    # flow arrow along the axis, ahead of the nose
    fa, fb = side(eye - 2, z1 - 6), side(eye + 30, z1 - 6)
    S.line(fa, fb, w=5)
    S.arrow(fb, fa, size=5)
    S.text((fb[0] + 2.5, fb[1]), "flow", anchor="lm", bg=False)

    # view from upstream: swept circle, rotation sense
    cxy = ups(0, 0)
    R_ = swept / 2
    circ = [(cxy[0] + R_ * math.cos(t), cxy[1] + R_ * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 241)]
    S.dashed(circ, dash=3.0, gap=2.2, w=3)
    a = math.radians(-38)
    ca, sa = math.cos(a), math.sin(a)
    pa_ = (cxy[0] - R_ * ca, cxy[1] - R_ * sa)
    pb_ = (cxy[0] + R_ * ca, cxy[1] + R_ * sa)
    pc_ = (cxy[0] + (R_ + 9) * ca, cxy[1] + (R_ + 9) * sa)
    S.line(pa_, pc_, w=3)
    S.arrow(pa_, pb_)
    S.arrow(pb_, pa_)
    S.line(pc_, (pc_[0] + 6, pc_[1]), w=3)
    S.text((pc_[0] + 7.5, pc_[1]), f"Ø{swept:.0f} swept", anchor="lm", bg=False)
    # clockwise seen from upstream (omega along +X): an arc outside the disc, screen-clockwise
    rr = R_ + 7
    arc = [(cxy[0] + rr * math.cos(t), cxy[1] + rr * math.sin(t)) for t in np.linspace(math.radians(-150), math.radians(-105), 30)]
    for p0, p1 in zip(arc, arc[1:]):
        S.line(p0, p1, w=3)
    S.arrow(arc[-1], arc[-4], size=4.6)
    S.text((arc[0][0] - 2, arc[0][1] - 3), "clockwise", size=TXT * 0.9, anchor="rm", bg=False)

    # plan: fin span at the right
    fin = next(b for b in boxes if b["name"] == "fin_0")
    S.dim(plan(tip, y1), plan(tip, y0), -18, f"{span:.0f} over fins",
          ext_from=(plan(fin["max"][0], y1), plan(fin["max"][0], y0)))

    # view labels
    caps = {"side": (side((x0 + x1) / 2, z0)[0], side(0, z0)[1] + 57), "upstream": (cxy[0], cxy[1] + R_ + 14),
            "plan": (plan((x0 + x1) / 2, y0)[0], plan(0, y0)[1] + 18)}
    for key, cap in (("side", "Side"), ("upstream", "From upstream"), ("plan", "Plan")):
        S.text(caps[key], cap.upper(), size=TXT * 0.85, weight="SemiBold", c=INK, bg=False)

    # title block
    tx, ty = org["upstream"][0], org["plan"][1] + 18
    S.d.rectangle([S.px((tx + 1, ty - 8)), S.px((tx + 17, ty - 6))], fill=ORANGE)
    S.d.text(S.px((tx, ty)), "Ovo", font=font(round(24 * KR), "Bold"), fill=INK)
    lines = [
        "Folding hydrokinetic charger, deployed,",
        "fins fitted. First-angle projection, mm,",
        "measured on the CAD model.",
        "",
        f"Stowed: Ø{D:.0f} × {tip:.0f} ({tip - eye:.0f} with tow eye)",
        f"Rotor: {P.N_BLADES} blades, Ø{swept:.0f} swept,",
        "folding forward into the tail",
        "",
        "Display build: not yet printed",
        "or tested in water.",
    ]
    S.d.multiline_text(S.px((tx + 0.5, ty + 36)), "\n".join(lines), font=font(round(TXT * 0.9 * KR), "Regular"),
                       fill=(0x3A, 0x39, 0x36), spacing=round(3.6 * KR))

    out = S.im.resize((round(S.im.width / 2), round(S.im.height / 2)), Image.Resampling.LANCZOS)
    n = save_png(out, OUT / "technical.png")
    manifest["technical.png"] = {"viewer": base, "px": list(out.size), "bytes": n, "scale_px_per_mm": KR / 2,
                                 "dims_mm": {"hull_length": tip, "overall_length": round(tip - eye, 1), "hull_d": D,
                                             "swept_d": swept, "fin_span": round(span, 1)}}
    log(f"technical.png {out.size[0]}x{out.size[1]} {n / 1e3:.0f} KB")


# --------------------------------------------------------------------------------------------
# Animations
# --------------------------------------------------------------------------------------------
FPS = 30
ANIM_RENDER, ANIM_OUT, GIF_W = (1920, 1200), (1280, 800), 800
PAPER_SHADOW = 0.24     # viewer: ground shadow opacity of theme=paper
INK_SHADOW = 0.42       # viewer: ground shadow opacity of theme=ink


def smooth(t: float) -> float:
    t = min(1.0, max(0.0, t))
    return t * t * t * (t * (6 * t - 15) + 10)


def ramp(t: float, a: float, b: float) -> float:
    return smooth((t - a) / (b - a))


def fit(v: Viewer, patch: dict) -> dict:
    """Refit the viewer camera for a state; returns {c, r, az, el} (the lights and ground follow each fit)."""
    v.run(f"window.setState({json.dumps(patch)})")
    return v.js("window.camInfo()")


def cam_at(c0: dict, c1: dict, p: float, pa: float | None = None, daz: float = 0.0) -> dict:
    """Camera between two fits: centre and distance at progress p, azimuth / elevation at pa (default p)."""
    pa = p if pa is None else pa
    lerp = lambda a, b, t: a + (b - a) * t       # noqa: E731
    return {"c": [lerp(a, b, p) for a, b in zip(c0["c"], c1["c"])], "r": c0["r"] * (c1["r"] / c0["r"]) ** p,
            "az": lerp(c0["az"], c1["az"], pa) + daz, "el": lerp(c0["el"], c1["el"], pa)}


def render_frames(v: Viewer, frames: list[dict], out_dir: Path, tag: str) -> None:
    """frames: [{cam, states: [state patch per sub-frame], look?}] -> out_dir/0000.png. The sub-frames are
    averaged in the page (window.accumFrame), so a frame costs one screenshot however many there are."""
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    t0 = time.time()
    for i, f in enumerate(frames):
        if f.get("look"):
            v.run(f"window.setLook({json.dumps(f['look'])})")
        v.run(f"window.accumFrame({json.dumps(f['states'])}, {json.dumps(f['cam'])})")
        im = v.shot(WORK / "frame.png")
        down(im, ANIM_OUT).save(out_dir / f"{i:04d}.png", compress_level=1)
        if i % 15 == 0:
            log(f"{tag}: frame {i + 1}/{len(frames)} ({time.time() - t0:.0f} s)")


def preview_only(v: Viewer, frames: list[dict], tag: str) -> bool:
    """MEDIA_PREVIEW=0,45,90 (frame indices): render just those frames to out/work/.../preview, no encoding."""
    sel = os.environ.get("MEDIA_PREVIEW")
    if not sel:
        return False
    idx = [int(x) for x in sel.split(",") if int(x) < len(frames)]
    render_frames(v, [frames[i] for i in idx], WORK / f"preview/{tag}", tag)
    log(f"{tag}: preview frames {idx} of {len(frames)} in {WORK / 'preview' / tag}")
    return True


def encode(frame_dir: Path, name: str, manifest: dict, gif: bool = True) -> None:
    import imageio_ffmpeg
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    mp4 = OUT / f"{name}.mp4"
    for crf in (16, 19, 22, 25, 28):
        subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", str(frame_dir / "%04d.png"),
                        "-c:v", "libx264", "-preset", "slow", "-crf", str(crf), "-pix_fmt", "yuv420p",
                        "-movflags", "+faststart", "-an", str(mp4)], check=True)
        if mp4.stat().st_size <= MP4_BUDGET:
            break
    nfr = len(list(frame_dir.glob("*.png")))
    manifest[mp4.name] = {"px": list(ANIM_OUT), "fps": FPS, "frames": nfr, "seconds": round(nfr / FPS, 2),
                          "codec": "H.264 yuv420p", "crf": crf, "bytes": mp4.stat().st_size}
    log(f"{mp4.name} {mp4.stat().st_size / 1e6:.2f} MB (crf {crf})")
    if not gif:
        return
    out = OUT / f"{name}.gif"
    for w, colors in ((GIF_W, 256), (GIF_W, 192), (720, 160), (640, 128)):
        vf = (f"scale={w}:-2:flags=lanczos,split[a][b];[a]palettegen=max_colors={colors}:stats_mode=full[p];"
              f"[b][p]paletteuse=dither=sierra2_4a:diff_mode=rectangle")
        subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", str(frame_dir / "%04d.png"),
                        "-vf", vf, "-loop", "0", str(out)], check=True)
        if out.stat().st_size <= GIF_BUDGET:
            break
    manifest[out.name] = {"px": [w, round(w * ANIM_OUT[1] / ANIM_OUT[0])], "fps": FPS, "frames": nfr,
                          "colors": colors, "bytes": out.stat().st_size}
    log(f"{out.name} {w} px, {colors} colours, {out.stat().st_size / 1e6:.2f} MB")


def spin_angle(t: float, t_on: float, t_full: float, w_max: float) -> float:
    """Rotor angle (deg) at t: angular speed eases from 0 at t_on to w_max (deg/s) at t_full, then holds."""
    if t <= t_on:
        return 0.0
    n = 400
    ts = np.linspace(t_on, t, n)
    w = np.array([w_max * ramp(x, t_on, t_full) for x in ts])
    return float(np.trapezoid(w, ts))


def crossfade_tail(frame_dir: Path, n_fade: int) -> None:
    """Blend the last n_fade frames toward frame 0 so the loop closes without a hard cut (P10): frame N-n_fade+j
    becomes (1 - a) * itself + a * frame 0 with a = (j + 1) / (n_fade + 1); the loop then restarts on frame 0."""
    frames = sorted(frame_dir.glob("*.png"))
    first = np.asarray(Image.open(frames[0]).convert("RGB"), dtype=np.float32)
    for j, f in enumerate(frames[-n_fade:]):
        a = (j + 1) / (n_fade + 1)
        a = a * a * (3 - 2 * a)                   # smoothstep
        im = np.asarray(Image.open(f).convert("RGB"), dtype=np.float32)
        Image.fromarray(np.clip((1 - a) * im + a * first + 0.5, 0, 255).astype(np.uint8)).save(f, compress_level=1)


def anim_fold(v: Viewer, manifest: dict, keep: bool, dark: bool = False) -> None:
    """Stowed egg -> blades unfold -> rotor spins up; the camera swings from the nose side to the tail side, holds,
    and the last half second cross-fades back to the first frame, so the GIF and the page's looping video have no
    seam. The spin speed is a presentation speed (90 rpm), not the rotor's running speed. dark=True renders the
    same animation on the `ink` backdrop as fold_dark.mp4 (no GIF), for the project page's dark theme."""
    T, step = 5.3, 0.3                            # seconds; motion blur: one sub-frame per 0.3 deg of travel
    t_fold, t_spin, w_max = (0.6, 2.0), (1.85, 3.0), 540.0
    n_fade = 16                                   # 0.53 s cross-fade to frame 0 at the end (after a 0.2 s hold)
    shutter = 0.5 / FPS                           # 180 deg shutter
    look = DARK if dark else LIGHT
    name = "fold_dark" if dark else "fold"
    P = _params()
    v.anim_open(f"view=-38,20&fold=0&fit=tight&fill=0.62&{look}")
    max_sub = 3 if v.soft else 36
    c0 = fit(v, {"fold": 0, "spin": 0, "view": "-38,20", "fill": 0.54, "fitspin": 0})
    c1 = fit(v, {"fold": 1, "spin": 0, "view": "50,18", "fill": 0.86, "fitspin": 12})
    fit(v, {"fold": 1, "spin": 0, "view": "6,19", "fitspin": 12})   # lights for the middle of the swing
    # ground under the swept disc: a spinning blade passes through the bottom of it
    v.run(f"window.setLook({{groundz: {-P.ROTOR_R_ENVELOPE - 3}, shadowop: {(INK_SHADOW if dark else PAPER_SHADOW) * 0.65}}})")
    frames = []
    for i in range(round(T * FPS)):
        t = i / FPS
        # the dolly-out is done by the time the blades are open; the swing to the tail side runs on
        cam = cam_at(c0, c1, ramp(t, 0.25, 1.75), pa=ramp(t, 0.15, 4.0))
        a, b = t - shutter / 2, t + shutter / 2
        travel = max(abs(spin_angle(b, *t_spin, w_max) - spin_angle(a, *t_spin, w_max)),
                     90.0 * abs(ramp(b, *t_fold) - ramp(a, *t_fold)))      # deg the blades turn in the shutter
        n = min(max_sub, max(1, math.ceil(travel / step)))
        ts = [t] if n == 1 else [a + (b - a) * k / (n - 1) for k in range(n)]
        states = [{"fold": round(ramp(x, *t_fold), 5), "spin": round(spin_angle(x, *t_spin, w_max), 3)} for x in ts]
        frames.append({"cam": cam, "states": states})
    d = WORK / f"frames/{name}"
    if preview_only(v, frames, name):
        return
    render_frames(v, frames, d, name)
    crossfade_tail(d, n_fade)
    encode(d, name, manifest, gif=not dark)
    note = ("stowed -> unfold -> spin-up -> cross-fade back to the first frame (seamless loop); spin shown at a "
            "presentation speed (90 rpm), not the running speed")
    manifest[f"{name}.mp4"]["note"] = note
    if not dark:
        manifest["fold.gif"]["note"] = note
    if not keep:
        shutil.rmtree(d)


def anim_exploded(v: Viewer, manifest: dict, keep: bool) -> None:
    """Assembled -> exploded (build order) -> assembled, a seamless 5 s loop with a slow camera sway."""
    T = 5.0
    v.anim_open(f"view=20,18&spin=30&fold=1&fit=tight&fill=0.86&{LIGHT}")
    c1 = fit(v, {"explode": 1, "view": "20,18", "fill": 0.8, "shadow": False})
    c0 = fit(v, {"explode": 0, "view": "20,18", "fill": 0.7, "shadow": True})   # last: ground and lights
    frames = []
    for i in range(round(T * FPS)):
        t = i / FPS
        e = ramp(t, 0.5, 2.0) if t < 2.5 else 1.0 - ramp(t, 3.0, 4.5)
        cam = cam_at(c0, c1, e, daz=5.0 * math.sin(2 * math.pi * t / T))
        op = PAPER_SHADOW * (1.0 - ramp(e, 0.02, 0.3))
        frames.append({"cam": cam, "states": [{"explode": round(e, 5)}], "look": {"shadowop": op, "ground": op > 0.002}})
    d = WORK / "frames/exploded"
    if preview_only(v, frames, "exploded"):
        return
    render_frames(v, frames, d, "exploded")
    encode(d, "exploded", manifest, gif=False)
    (OUT / "exploded.gif").unlink(missing_ok=True)
    manifest.pop("exploded.gif", None)
    if not keep:
        shutil.rmtree(d)


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
JOBS = {name: (lambda n: (lambda v, m, k: still(v, n, m)))(name) for name in STILLS}
JOBS.update({name: (lambda n: (lambda v, m, k: panels(v, n, m)))(name) for name in PANELS})
JOBS.update({"og_card": lambda v, m, k: og_card(v, m), "technical": lambda v, m, k: technical(v, m),
             "fold": anim_fold, "fold_dark": lambda v, m, k: anim_fold(v, m, k, dark=True),
             "exploded_anim": anim_exploded})
GROUPS = {"hero": ["hero_light", "hero_dark"], "details": ["details_hub", "details_fin_rail", "details_nose", "details_blade"],
          "stills": list(STILLS) + list(PANELS) + ["og_card", "technical"], "anim": ["fold", "fold_dark", "exploded_anim"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", help="comma list of jobs or groups (see --list)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--cpu", action="store_true", help="software GL (no GPU flags)")
    ap.add_argument("--keep", action="store_true", help="keep the animation frames")
    ap.add_argument("--sheet", nargs="+", metavar=("OUT", "QUERY"), help="contact sheet: OUT.png then viewer queries")
    a = ap.parse_args()
    if a.list:
        for k in JOBS:
            print(k)
        print("groups:", ", ".join(f"{g} = {'+'.join(v)}" for g, v in GROUPS.items()))
        return
    names = list(JOBS)
    if a.only:
        names = []
        for n in a.only.split(","):
            names += GROUPS.get(n, [n])
        bad = [n for n in names if n not in JOBS]
        if bad:
            sys.exit(f"unknown job(s): {bad}; see --list")
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    mpath = OUT / "media.json"
    manifest = json.loads(mpath.read_text()).get("files", {}) if mpath.exists() else {}
    t0 = time.time()
    with R.heavy_lock(), R.serve() as base:
        v = Viewer(base, gpu=not a.cpu)
        try:
            if a.sheet:
                contact_sheet(v, a.sheet[1:], Path(a.sheet[0]))
                return
            for n in names:
                log(f"-- {n}")
                JOBS[n](v, manifest, a.keep)
        finally:
            v.close()
    if not a.sheet:
        cad_key = json.loads((REPO / "out/cad/checks/assembly.json").read_text()).get("cache_key")
        mpath.write_text(json.dumps({
            "generated_by": "scripts/media.py", "source": "out/web/ovo.glb + manifest.json via viewer/index.html",
            "cad_key": cad_key, "label_font": Path(font(12).path).name,
            "look": {"light": LIGHT, "dark": DARK, "lens_deg": 28},
            "note": "geometry renders of the CAD model; the display build has not been printed or tested in water",
            "files": dict(sorted(manifest.items()))}, indent=1, ensure_ascii=False) + "\n")
        log(f"done in {time.time() - t0:.0f} s -> {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
