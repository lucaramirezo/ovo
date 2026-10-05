"""Headless renders of the Ovo assembly through viewer/index.html.

Serves the repo over HTTP, drives Chrome with agent-browser (own named session so
it never collides with other browser users), and screenshots each requested view.

  .venv/bin/python scripts/render.py out/renders/hero.png "fold=1&view=hero"
  .venv/bin/python scripts/render.py --batch views.json      # [{"out":..., "q":...}, ...]
  .venv/bin/python scripts/render.py --frames out/frames/fold "view=hero" \
        --vary fold 0 1 48                                    # animation frames

Heavy-ish (Chrome + SwiftShader): takes the out/work/heavy.lock so only one heavy
job runs at a time on this machine.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import functools
import http.server
import json
import os
import shutil
import socketserver
import subprocess
import threading
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SESSION = os.environ.get("AGENT_BROWSER_SESSION", "ovo-render")
AB = shutil.which("agent-browser")
if AB is None:
    raise SystemExit("agent-browser not on PATH: npm install -g agent-browser && agent-browser install")


def ab(*args: str, check: bool = True) -> str:
    env = {**os.environ, "AGENT_BROWSER_SESSION": SESSION}
    r = subprocess.run([AB, "--session", SESSION, *args], capture_output=True, text=True, env=env, timeout=180)
    if check and r.returncode != 0:
        raise RuntimeError(f"agent-browser {' '.join(args)} failed: {r.stderr or r.stdout}")
    return r.stdout.strip()


@contextlib.contextmanager
def heavy_lock():
    lock = REPO / "out/work/heavy.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        yield


@contextlib.contextmanager
def serve():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(REPO))
    handler.log_message = lambda *a, **k: None
    with socketserver.TCPServer(("127.0.0.1", 0), handler) as httpd:
        t = threading.Thread(target=httpd.serve_forever, daemon=True)
        t.start()
        try:
            yield f"http://127.0.0.1:{httpd.server_address[1]}"
        finally:
            httpd.shutdown()


def shoot(base: str, jobs: list[dict], w: int, h: int) -> None:
    ab("set", "viewport", str(w), str(h))
    loaded = None
    for job in jobs:
        out = Path(job["out"])
        out.parent.mkdir(parents=True, exist_ok=True)
        page = job.get("page_q", "")
        url = f"{base}/viewer/index.html?w={w}&h={h}&dpr=1&{page}"
        if url != loaded:
            ab("open", url)
            ab("wait", "#ready")
            loaded = url
        if job.get("state") is not None:
            ab("eval", f"window.setState({json.dumps(job['state'])})")
        ab("screenshot", str(out.resolve()))
        print(out)


def parse_q(qs: str) -> dict:
    out = {}
    for kv in filter(None, qs.split("&")):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?")
    ap.add_argument("q", nargs="?", default="")
    ap.add_argument("--batch")
    ap.add_argument("--frames")
    ap.add_argument("--vary", nargs=4, metavar=("KEY", "FROM", "TO", "N"))
    ap.add_argument("-W", type=int, default=1600)
    ap.add_argument("-H", type=int, default=1000)
    a = ap.parse_args()

    jobs: list[dict] = []
    if a.batch:
        for j in json.loads(Path(a.batch).read_text()):
            jobs.append({"out": j["out"], "page_q": j.get("q", "")})
    elif a.frames:
        key, lo, hi, n = a.vary[0], float(a.vary[1]), float(a.vary[2]), int(a.vary[3])
        for i in range(n):
            v = lo + (hi - lo) * i / max(n - 1, 1)
            jobs.append({"out": f"{a.frames}_{i:03d}.png", "page_q": a.q, "state": {key: v}})
    else:
        jobs.append({"out": a.out, "page_q": a.q})

    with heavy_lock(), serve() as base:
        try:
            shoot(base, jobs, a.W, a.H)
        finally:
            ab("close", check=False)


if __name__ == "__main__":
    main()
