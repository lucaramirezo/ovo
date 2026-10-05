"""Ovo assembly (integrator).

Collects every part module's ``parts()`` (see common.py for the part dict contract), adds the
detent magnets no builder models, poses the rotor in its states and checks the whole machine:

  STOWED    blades folded FOLD_DEG (90) into the pockets, rotor at the stowing index (spin 0:
            blade i at psi 120 i sits over the pocket at theta = -psi)
  DEPLOYED  blades rotated by their pivot (blade module) to fold 0
  FOLD_30 / FOLD_60  intermediate folds (fold angle as in blade.blade_world: 0 deployed)
  CUTAWAY   cutaway shells + internals, blades deployed (display variant; the cutaway front shell is cut
            from the v2 front shell, so it shows the gland seat with the gland in it)
  PACKED    STOWED with the fins taken off and the flush blanking strips in the rail grooves
            (finding T08 / owner decision: the clean 75 mm egg; used for the stowed STEP and render)
  V2        STOWED with the v2 hardware: the v2 front shell (gland seat, tapped bore) and the M12 gland in
            place of the display front shell (the viewer's v2 mode; only pairs with a v2 part are new)

Variants (owner decision 2026-10-05, nose round): a part flagged v2=True (the M12 cable gland, the v2 front
shell front_shell_v2, role "v2") exists only in CUTAWAY / V2 (role "v2" never in CUTAWAY: the cutaway shells
replace it there); v2=False (the display front shell, clean nose) everywhere else. The display build has no
gland and no plug. packed=False parts (fins) leave the PACKED state, packed=True parts (blanking strips)
exist only there.

Outputs (scripts/build.py runs this after the part modules):
  out/cad/checks/assembly.json              interference matrix, key clearances, STL audit, BOM
  out/cad/step/ovo_assembly_{stowed,deployed}.step
  out/web/ovo.glb + out/web/manifest.json   (viewer: roles, rotor flags, explode, pivots)
  out/work/cad/integrate/                   part cache (BREP, keyed on the sources), meshes

Interference is computed on manifold3d meshes (0.02 mm chord, shells 0.03): exact OCC booleans
against the full shells take minutes each. Overlaps below ~0.5 mm3 on curved contacts are
facets of the tessellation; flat contacts tessellate exactly.

  scripts/cadrun.sh .venv/bin/python -m ovo.cad.assembly [--force] [--no-web] [--no-verify]
"""

from __future__ import annotations

import fnmatch
import hashlib
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from build123d import Compound, Location, Plane, Shape, Solid, export_brep, export_step, import_brep

from .. import params as P
from . import common as C

REPO = Path(__file__).resolve().parents[2]
WORK = REPO / "out" / "work" / "cad" / "integrate"
CACHE = WORK / "cache"
OUT_STEP = REPO / "out" / "cad" / "step"
OUT_STL = REPO / "out" / "cad" / "stl"
OUT_CHK = REPO / "out" / "cad" / "checks"
OUT_WEB = REPO / "out" / "web"

MODULES = ("shells", "rotor", "blade", "drivetrain", "fins", "coupon")
MESH_TOL = {"default": (0.02, 0.1), "shell": (0.03, 0.1)}
# web meshes: ABSOLUTE chord deflection (export_web.occ_mesh), sub-pixel at 1600 px; blades finer so
# the trailing edge reads straight (recheck R1)
WEB_TOL = {"default": (0.05, 0.25), "shell": (0.05, 0.25), "blade": (0.03, 0.2)}
SEARCH = 3.0                       # min_gap search length (mm): gaps above it read as ">= 3"
DEFECT_MM3 = 0.5                   # unexpected overlap above this is a defect
FOLD_STATES = {"stowed": 90.0, "fold_60": 60.0, "fold_30": 30.0, "deployed": 0.0, "packed": 90.0, "v2": 90.0}
# Engineering call 2026-10-05: CLEAR_ROTATE is radial per side for rotating fits (knuckle on pin 0.15);
# swept / running clearances between moving bodies are RUN_CLEAR 0.3 (params).
RUN_CLEAR = float(P.RUN_CLEAR)


# --------------------------------------------------------------------------
# Part collection (cached as BREP, keyed on every source the parts depend on)
# --------------------------------------------------------------------------

def _sources() -> list[Path]:
    src = sorted(p for p in (REPO / "ovo" / "cad").glob("*.py") if p.name != "assembly.py")
    src += [REPO / "ovo" / n for n in ("params.py", "profile.py", "layout.py", "layout_internal.py")]
    src += sorted((REPO / "data").glob("*.json"))
    return [p for p in src if p.exists()]


def cache_key() -> str:
    h = hashlib.sha1()
    for p in _sources():
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()[:12]


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def _build_fresh(cdir: Path) -> None:
    """Call every module's parts(), write world and print BREPs and the metadata."""
    cdir.mkdir(parents=True, exist_ok=True)
    meta, timing = [], {}
    for mod in MODULES:
        t0 = time.time()
        m = importlib.import_module(f"ovo.cad.{mod}")
        for p in m.parts():
            n = p["name"]
            export_brep(p["shape"], str(cdir / f"{n}.brep"))
            if p.get("print") is not None:
                export_brep(p["print"], str(cdir / f"{n}.print.brep"))
            d = {k: _jsonable(v) for k, v in p.items() if k not in ("shape", "print")}
            d.update(module=mod, has_print=p.get("print") is not None)
            meta.append(d)
        timing[mod] = round(time.time() - t0, 1)
        print(f"  parts {mod}: {timing[mod]} s", flush=True)
    # geometry the integrator needs beyond the part dicts
    from . import blade as BL
    from . import shells as SH
    info = BL._build()[1]
    _, lay, _ = BL.context()
    blades = []
    for i in range(lay.n_blades):
        blades.append(dict(T_stowed=BL.blade_transform(i, P.FOLD_DEG).tolist(),
                           T_deployed=BL.blade_transform(i, 0.0).tolist(), pivot=BL.pivot(i)))
        export_brep(BL.blade_world(i, 0.0), str(cdir / f"blade_{i}.deployed.brep"))
    extras = dict(blade_seats=_jsonable({k: info[k] for k in ("deployed_seat", "stowed_seat")}),
                  blades=blades, stow_detent=_jsonable({k: v for k, v in SH.stow_detent_geometry().items() if k != "top_grid"}),
                  blade_theta=list(SH.BLADE_THETA), timing_s=timing)
    (cdir / "meta.json").write_text(json.dumps(dict(parts=meta, extras=extras), indent=1))


_COLLECTED: dict = {}
_ACTIVE_CDIR: Path | None = None


def collect(force: bool = False) -> tuple[list[dict], dict]:
    """(parts, extras): every module part plus the integrator's magnets, shapes in the WORLD
    stowed pose, 'print' loaded when the part has one."""
    key = cache_key()
    if key in _COLLECTED and not force:
        return _COLLECTED[key]
    cdir = CACHE / key
    if force or not (cdir / "meta.json").exists():
        import shutil
        for old in CACHE.glob("*"):
            if old.is_dir():
                shutil.rmtree(old, ignore_errors=True)            # BREPs and meshes, every key
        print(f"building parts (cache {key})", flush=True)
        _build_fresh(cdir)
    data = json.loads((cdir / "meta.json").read_text())
    parts = []
    for d in data["parts"]:
        p = dict(d)
        p["shape"] = import_brep(str(cdir / f"{d['name']}.brep"))
        p["print"] = import_brep(str(cdir / f"{d['name']}.print.brep")) if d["has_print"] else None
        parts.append(p)
    extras = data["extras"]
    extras["cache_dir"] = str(cdir.relative_to(REPO))
    global _ACTIVE_CDIR
    _ACTIVE_CDIR = cdir           # meshes go next to these BREPs even if a source changes mid-run
    extras["deployed_ref"] = [import_brep(str(cdir / f"blade_{i}.deployed.brep"))
                              for i in range(len(extras["blades"]))]
    mo = extras["stow_detent"]["mouth"]
    extras["stow_detent_theta"] = math.degrees(math.atan2(mo[1], mo[2]))      # blade 0's boss azimuth
    parts += magnets(extras)
    parts += oring_cutaway(parts)
    _COLLECTED[key] = (parts, extras)
    return parts, extras


# --------------------------------------------------------------------------
# Detent magnets (no builder models them): 3 in the blade heels, 3 in the blade shanks
# (stowed detent), 3 in the aft shell bosses. The hub's 3 come from rotor.py.
# --------------------------------------------------------------------------

MAG_D, MAG_H = P.MAGNET["D"], P.MAGNET["H"]


def _loc(T) -> Location:
    from OCP.gp import gp_Trsf
    T = np.asarray(T, float)
    tr = gp_Trsf()
    tr.SetValues(*(float(T[i, j]) for i in range(3) for j in range(4)))
    return Location(tr)


def block_solid(centre, size_s_c_n, T=None) -> Solid:
    """N52 block (finding T03) in the blade frame: centre, sizes along (s, c, n), 0.1 mm edge chamfers."""
    a, b, c = (float(v) for v in size_s_c_n)
    s = Solid.make_box(a, b, c).moved(Location((float(centre[0]) - a / 2, float(centre[1]) - b / 2,
                                                float(centre[2]) - c / 2)))
    try:
        s = s.chamfer(0.1, None, s.edges())
    except Exception:  # noqa: BLE001
        pass
    return s


def magnet_solid(back, axis) -> Solid:
    """N52 disc D4 x 2 with 0.1 mm edge rounds, back face at ``back``, along unit ``axis``."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    s = Solid.make_cylinder(MAG_D / 2, MAG_H, Plane(origin=tuple(map(float, back)), z_dir=tuple(map(float, a))))
    try:
        s = s.chamfer(0.1, None, s.edges())
    except Exception:  # noqa: BLE001
        pass
    return s


def _own(shape) -> Solid:
    """Independent copy (own TShape): STEP writes shapes that share a TShape as instances of ONE
    product, so the nine magnets would lose their names."""
    import copy
    return copy.deepcopy(shape)


def magnets(extras: dict) -> list[dict]:
    seats = extras["blade_seats"]
    dep, sto = seats["deployed_seat"], seats["stowed_seat"]
    a_d = np.asarray(dep["axis"], float)
    c_d = np.asarray(dep["centre"], float)
    heel = block_solid(c_d, dep["size_s_c_n"])                                   # blade frame (T03)
    a_s = np.asarray(sto["axis_out_blade_frame"], float)
    m_s = np.asarray(sto["mouth_blade_frame_mm"], float)
    shank = magnet_solid(m_s + (sto["t_face_mm"] - MAG_H) * a_s, a_s)
    g = extras["stow_detent"]
    n = np.asarray(g["n"], float)
    tail0 = magnet_solid(np.asarray(g["magnet_centre"], float) - MAG_H / 2 * n, n)  # theta 0 (blade 0)
    out = []
    for i, b in enumerate(extras["blades"]):
        T = np.asarray(b["T_stowed"], float)
        piv = b["pivot"]
        out.append(dict(
            name=f"blade_magnet_heel_{i}", shape=_own(heel).moved(_loc(T)), print=None, material="magnet",
            role="body", rotor=True, qty=1, bought=True, explode=[0, 0, 0], pivot=piv, module="assembly",
            axis_blade_frame=a_d.tolist(), has_print=False,
            notes="N52 block 6 x 2 x 2 (6 along the pin, magnetised through 2 along the blade normal), "
                  "CA-glued in the heel stop-face slot (6.1 x 2.1, face 0.2 below the stop face); pairs with "
                  "hub_magnet_%d when deployed (attracting, 0.4 mm gap)." % i))
        out.append(dict(
            name=f"blade_magnet_stow_{i}", shape=_own(shank).moved(_loc(T)), print=None, material="magnet",
            role="body", rotor=True, qty=1, bought=True, explode=[0, 0, 0], pivot=piv, module="assembly",
            axis_blade_frame=a_s.tolist(), has_print=False,
            notes="N52 disc 4 x 2 axial, CA-glued in the shank inner-face seat (D4.1, along the face "
                  "normal, face 0.1 below the seat mouth); pairs with tail_magnet_%d when stowed (attracting)." % i))
        th = extras["blade_theta"][i]
        lo, hi = _cutaway_theta()
        boss = ((th + extras["stow_detent_theta"] + 180.0) % 360.0) - 180.0
        out.append(dict(
            name=f"tail_magnet_{i}", shape=_own(C.at_azimuth(tail0, th) if th else tail0), print=None,
            material="magnet", role="shell" if lo <= boss <= hi else "body",   # boss cut away in CUTAWAY
            rotor=False, qty=1, bought=True, explode=[0, 0, 0],
            module="assembly", has_print=False, theta_deg=th,
            notes="N52 disc 4 x 2 axial, CA-glued in the aft shell's stowed-detent boss (D4.1 x 2.05, axis on the "
                  "blade face normal, 49 deg from radial; face 0.05 below the lowest blade-driven boss top over the "
                  "hole, so it may stand up to ~0.2 mm proud of the boss top on one side); holds blade %d stowed "
                  "(pair ~0.42 mm apart, ~1.7 deg, ~2.3 N)." % i))
    return out


def oring_cutaway(parts: list[dict]) -> list[dict]:
    """Finding T09: the O-ring is a full ring, so in the CUTAWAY state it arced across the removed
    sector like a stray cable. It now takes the shell role (hidden with the full shells) and a copy
    cut by the same sector solid as the aft shell (shells.cutaway_sector) takes shell_cutaway."""
    from . import shells as SH
    B = _by(parts)
    o = B.get("oring_141")
    if o is None:
        return []
    o["role"] = "shell"
    cut = o["shape"].cut(SH.cutaway_sector(False))
    sols = cut.solids()
    cut = sols[0] if len(sols) == 1 else Compound(list(sols))
    q = {k: v for k, v in o.items() if k not in ("shape", "print")}
    q.update(name="oring_141_cutaway", shape=cut, print=None, role="shell_cutaway", module="assembly",
             has_print=False, bom_parent="oring_141",
             notes="oring_141 with the cutaway sector theta %g..%g removed (display only)." % SH.CUTAWAY_THETA)
    return [q]


def cutaway_crossers(parts: list[dict], r_min: float = 25.0) -> list[dict]:
    """Internal parts with material inside the cutaway sector beyond r_min (the T09 check: what else
    would float across the removed sector in the CUTAWAY render): mesh vertices in the deployed pose."""
    from . import shells as SH
    th0, th1 = SH.CUTAWAY_THETA
    out = []
    for p in parts:
        if p.get("role") != "internal":
            continue
        T = fold_T(p, 0.0)
        V = np.asarray(part_mesh(p).vertices)
        if T is not None:
            V = V @ T[:3, :3].T + T[:3, 3]
        r = np.hypot(V[:, 1], V[:, 2])
        th = np.degrees(np.arctan2(V[:, 1], V[:, 2]))
        m = (r > r_min) & (th >= th0) & (th <= th1)
        if m.any():
            out.append(dict(name=p["name"], n_vertices=int(m.sum()), r_max_mm=round(float(r[m].max()), 2),
                            x_range_mm=[round(float(V[m, 0].min()), 2), round(float(V[m, 0].max()), 2)]))
    return out


def _cutaway_theta() -> tuple[float, float]:
    from . import shells as SH
    return tuple(float(v) for v in SH.CUTAWAY_THETA)


# --------------------------------------------------------------------------
# Poses
# --------------------------------------------------------------------------

def rot_T(origin, axis, deg: float) -> np.ndarray:
    """4x4 right-hand rotation by deg about the axis through origin (three.js setFromAxisAngle)."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    t = math.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    R = np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * K @ K
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = np.asarray(origin, float) - R @ np.asarray(origin, float)
    return T


def fold_T(part: dict, fold_deg: float) -> np.ndarray | None:
    """Transform from the stowed shape to blade fold ``fold_deg`` (None: part does not fold)."""
    piv = part.get("pivot")
    if not piv:
        return None
    return rot_T(piv["origin"], piv["axis"], piv["deployed_deg"] * (P.FOLD_DEG - fold_deg) / P.FOLD_DEG)


def posed_shape(part: dict, fold_deg: float) -> Shape:
    T = fold_T(part, fold_deg)
    return part["shape"] if T is None else part["shape"].moved(_loc(T))


def state_parts(parts: list[dict], state: str) -> list[dict]:
    """Parts present in a state (no hidden parts; cutaway shells only in CUTAWAY; full shells and the role "v2"
    shell never in CUTAWAY; parts flagged packed=True only in PACKED, packed=False never in PACKED; v2=True only
    in CUTAWAY / V2, v2=False never there)."""
    cut = state == "cutaway"
    packed = state == "packed"
    v2 = cut or state == "v2"
    out = []
    for p in parts:
        r = p.get("role", "body")
        if r == "hidden" or (r in ("shell", "v2") and cut) or (r == "shell_cutaway" and not cut):
            continue
        if (p.get("packed") is True and not packed) or (p.get("packed") is False and packed):
            continue
        if (p.get("v2") is True and not v2) or (p.get("v2") is False and v2):
            continue
        out.append(p)
    return out


def state_fold(state: str) -> float:
    return 0.0 if state == "cutaway" else FOLD_STATES[state]


def _is_packed_part(p: dict) -> bool:
    return p.get("packed") is True


# --------------------------------------------------------------------------
# Meshes (cached per part, world stowed pose)
# --------------------------------------------------------------------------

_MESH: dict = {}


def part_mesh(p: dict):
    """trimesh of the part's world stowed shape (cached next to the BREPs)."""
    import trimesh
    n = p["name"]
    if n in _MESH:
        return _MESH[n]
    cdir = _ACTIVE_CDIR or CACHE / cache_key()
    tol, ang = MESH_TOL.get(p.get("material") if p.get("material") == "shell" else "default")
    f = cdir / f"{n}.mesh_{tol:g}_{ang:g}.npz"
    if f.exists():
        z = np.load(f)
        m = trimesh.Trimesh(z["v"], z["f"], process=False)
    else:
        _clean(p["shape"])
        v, fc = p["shape"].tessellate(tol, ang)
        m = trimesh.Trimesh(np.array([(q.X, q.Y, q.Z) for q in v]), np.array(fc), process=True)
        m.merge_vertices()
        C.closed_mesh(m)                      # sliver gaps the mesher leaves at some fillet junctions
        if p.get("module") != "assembly":
            np.savez_compressed(f, v=m.vertices, f=m.faces)
    _MESH[n] = m
    return m


def _clean(shape) -> None:
    """Drop any triangulation the shape carries (a BREP keeps the mesh of whoever tessellated it
    first, and OCC reuses it whatever tolerance is asked next)."""
    from OCP.BRepTools import BRepTools
    BRepTools.Clean_s(shape.wrapped)


def manifold(mesh, T: np.ndarray | None = None):
    import manifold3d as m3
    mm = m3.Mesh(vert_properties=np.asarray(mesh.vertices, np.float32), tri_verts=np.asarray(mesh.faces, np.uint32))
    mm.merge()
    m = m3.Manifold(mm)
    if m.is_empty() or m.status() != m3.Error.NoError:      # an open mesh would read as "no overlap"
        raise RuntimeError(f"mesh is not a closed manifold ({m.status()}, {len(mesh.faces)} faces)")
    if T is not None:
        m = m.transform(np.asarray(T, float)[:3, :4])
    return m


_MF: dict = {}


def part_manifold(p: dict, fold_deg: float = P.FOLD_DEG):
    key = (p["name"], round(fold_deg, 6) if p.get("pivot") else None)
    if key not in _MF:
        _MF[key] = manifold(part_mesh(p), fold_T(p, fold_deg))
    return _MF[key]


# --------------------------------------------------------------------------
# Interference matrix
# --------------------------------------------------------------------------

# Intended overlaps and contacts: (pattern a, pattern b, kind, note). Kinds: press (designed
# interference), thread (thread modelled at the major / hole diameter), seal (elastomer squeeze),
# part (sub-part of the same bought item), contact (designed zero-gap face contact).
EXPECTED = [
    ("hinge_pin_*", "hub", "press", "3 mm pin pressed into the D2.90 printed ears"),
    ("hub_magnet_*", "hub", "press", "6 x 2 x 2 block glued in the 6.1 x 2.1 stop-face seat (no overlap expected)"),
    ("spinner_cap", "spinner", "press", "tip cap: 3 crush ribs (D6.1 over the tips) in the D6 access bore"),
    ("bearing_608", "aft_shell*", "press", "608 outer race on 3 crush ribs (tips D21.9)"),
    ("bearing_688", "aft_shell*", "press", "688 outer race on 3 crush ribs (tips D15.9)"),
    ("shaft_collar", "claw_b", "press", "Ruland collar pressed into claw B's 3 crush ribs (D15.9)"),
    ("lip_seal_tc8x16x5", "aft_shell*", "press", "seal OD pressed into the D16 bore"),
    ("lip_seal_tc8x16x5", "shaft", "seal", "seal lip on the shaft"),
    ("oring_141*", "*_shell*", "seal", "O-ring squeezed in the equator gland"),
    ("inserts_m3_bell", "generator_bell", "press", "heat-set inserts melted into the bell"),
    ("screws_carrier_m3x8", "generator_stator", "thread", "M3 x 8 self-tapped into the D2.5 stator pilots (T21)"),
    ("screws_claw_a_m3x4", "inserts_m3_bell", "thread", "M3 screws in the inserts"),
    ("spinner_screw", "shaft", "thread", "M3 x 8 in the shaft's tapped end"),
    ("hub_grub", "hub", "thread", "M3 grub tapped into the PLA hub"),
    ("cable_gland_m12", "front_shell*", "thread", "M12x1.5 gland thread in the tapped seat"),
    ("bearing_608", "shaft", "contact", "inner race on the shaft (zero clearance)"),
    ("bearing_688", "shaft", "contact", "inner race on the shaft (slip fit)"),
    ("bearing_*_seals", "bearing_*", "part", "rubber seals of the same bearing"),
    ("lip_seal_spring", "lip_seal_tc8x16x5", "part", "garter spring of the same seal"),
]


def expected(a: str, b: str):
    for pa, pb, kind, note in EXPECTED:
        if (fnmatch.fnmatch(a, pa) and fnmatch.fnmatch(b, pb)) or (fnmatch.fnmatch(b, pa) and fnmatch.fnmatch(a, pb)):
            return kind, note
    return None


def _bbox(m):
    b = m.bounding_box()
    return np.array(b[:3]), np.array(b[3:])


_PAIR: dict = {}


def pair_metrics(pa: dict, pb: dict, fold: float) -> dict | None:
    """(overlap mm3, min gap mm) of two parts at a blade fold; None when the boxes are SEARCH apart."""
    ka = (pa["name"], fold if pa.get("pivot") else None)
    kb = (pb["name"], fold if pb.get("pivot") else None)
    key = (ka, kb) if ka <= kb else (kb, ka)
    if key in _PAIR:
        return _PAIR[key]
    ma, mb = part_manifold(pa, fold), part_manifold(pb, fold)
    (a0, a1), (b0, b1) = _bbox(ma), _bbox(mb)
    if np.any(a0 > b1 + SEARCH) or np.any(b0 > a1 + SEARCH):
        _PAIR[key] = None
        return None
    ov = 0.0
    if not (np.any(a0 > b1) or np.any(b0 > a1)):
        ov = float((ma ^ mb).volume())
    gap = 0.0 if ov > 1e-3 else float(ma.min_gap(mb, SEARCH))
    _PAIR[key] = dict(overlap_mm3=ov, min_gap_mm=gap)
    return _PAIR[key]


def interference(parts: list[dict], states=("stowed", "deployed", "fold_30", "fold_60", "cutaway", "packed", "v2")) -> dict:
    """Pairwise overlap volume and min gap for every pair of parts present in each state."""
    res = {}
    for st in states:
        t0 = time.time()
        fold = state_fold(st)
        S = state_parts(parts, st)
        rows = []
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                a, b = S[i], S[j]
                if st == "packed" and not (_is_packed_part(a) or _is_packed_part(b)):
                    continue                     # identical to the stowed row
                if st == "v2" and not (a.get("v2") is True or b.get("v2") is True):
                    continue                     # identical to the stowed row
                if st not in ("stowed", "packed", "v2") and not (a.get("pivot") or b.get("pivot")) and \
                        not (st == "cutaway" and ("cutaway" in a["name"] or "cutaway" in b["name"])):
                    continue                     # static pair: identical to the stowed row
                m = pair_metrics(a, b, fold)
                if m is None:
                    continue
                e = expected(a["name"], b["name"])
                rows.append(dict(a=a["name"], b=b["name"], overlap_mm3=round(m["overlap_mm3"], 4),
                                 min_gap_mm=round(m["min_gap_mm"], 4),
                                 expected=e[0] if e else None, note=e[1] if e else None))
        defects = [r for r in rows if r["overlap_mm3"] > DEFECT_MM3 and not r["expected"]]
        res[st] = dict(fold_deg=fold, n_parts=len(S), n_pairs_near=len(rows), rows=rows,
                       defects=defects,
                       expected_overlaps=[r for r in rows if r["expected"] and r["overlap_mm3"] > 0],
                       contacts=[r for r in rows if not r["expected"] and r["overlap_mm3"] <= DEFECT_MM3
                                 and (r["overlap_mm3"] > 0 or r["min_gap_mm"] < 0.02)],
                       seconds=round(time.time() - t0, 1))
        print(f"  matrix {st}: {len(S)} parts, {len(rows)} near pairs, {len(defects)} defects "
              f"({res[st]['seconds']} s)", flush=True)
    return res


# --------------------------------------------------------------------------
# Key clearances on the real solids (vs data/layout.json interface checks)
# --------------------------------------------------------------------------

def _by(parts):
    return {p["name"]: p for p in parts}


def gap_location(ma, mb, n: int = 60000) -> list | None:
    """Approximate midpoint of the closest approach of two trimeshes (surface samples)."""
    import trimesh
    from scipy.spatial import cKDTree
    lo = np.maximum(ma.bounds[0], mb.bounds[0]) - SEARCH
    hi = np.minimum(ma.bounds[1], mb.bounds[1]) + SEARCH
    pa = np.vstack([trimesh.sample.sample_surface(ma, n, seed=1)[0], ma.vertices])
    pb = np.vstack([trimesh.sample.sample_surface(mb, 2 * n, seed=2)[0], mb.vertices])
    pa = pa[np.all((pa >= lo) & (pa <= hi), axis=1)]
    pb = pb[np.all((pb >= lo) & (pb <= hi), axis=1)]
    if not len(pa) or not len(pb):
        return None
    d, k = cKDTree(pb).query(pa)
    j = int(np.argmin(d))
    return [round(float(v), 2) for v in (pa[j] + pb[k[j]]) / 2]


def _gap(pa, pb, fold=P.FOLD_DEG) -> float:
    return float(part_manifold(pa, fold).min_gap(part_manifold(pb, fold), SEARCH))


def _sphere(c, r: float):
    import manifold3d as m3
    return m3.Manifold.sphere(r, 96).translate(tuple(map(float, c)))


def hull_outset(dr: float = 0.05):
    """Manifold of the hull r(x) + dr (x 0..HULL_L), to isolate what stands proud of it."""
    import manifold3d as m3
    from .. import profile
    xs = (1 - np.cos(np.linspace(0, math.pi, 320))) / 2 * P.HULL_L
    rs = np.clip(profile.radius(xs), 0, None) + dr
    poly = [(0.0, -1.0), (float(rs[0]), -1.0)] + [(float(r), float(x)) for x, r in zip(xs, rs)] + \
           [(float(rs[-1]), P.HULL_L + 1.0), (0.0, P.HULL_L + 1.0)]
    cs = m3.CrossSection([poly])
    return m3.Manifold.revolve(cs, 360).rotate((0.0, 90.0, 0.0))     # chord error 0.0015 mm at r 37.5


def _force(gap: float) -> float:
    d = json.loads((REPO / "data" / "layout_rotor.json").read_text())["detents"]["force_vs_gap"]
    return float(np.interp(gap, d["gap_mm"], d["force_N"]))


BOSS_ZONE_R = float(P.DETENTS["stowed_rest_zone_r"])   # 8.25: the stowed rest (leaned boss + chin) around
                                                       # the boss face centre (shells.stow_detent_geometry)


def key_clearances(parts: list[dict], extras: dict) -> list[dict]:
    B = _by(parts)
    lay_checks = {c["check"]: c for c in json.loads((REPO / "data" / "layout.json").read_text())["interface_checks"]}
    out = []

    def add(check, value, ok, layout=None, limit=None, where=None, note=None):
        out.append(dict(check=check, value=value, ok=bool(ok), layout=layout, limit=limit, where=where, note=note))
        print(("  PASS " if ok else "  FAIL ") + check + ": " + str(value), flush=True)

    blades = [B[f"blade_{i}"] for i in range(P.N_BLADES)]
    aft, front, hub, spin = B["aft_shell"], B["front_shell"], B["hub"], B["spinner"]

    # 1. hub gap, lip seal to hub, rotor to static tail
    g = _gap(aft, hub)
    add("hub gap: static tail end (aft shell) to hub front", f"{g:.3f} mm", g >= P.HUB_GAP - 0.01,
        lay_checks["hub gap = HUB_GAP"]["value_mm"], f">= {P.HUB_GAP}",
        gap_location(part_mesh(aft), part_mesh(hub)))
    g2 = _gap(aft, spin)
    add("spinner to aft shell", f"{g2:.3f} mm", g2 >= P.HUB_GAP - 0.01, None, f">= {P.HUB_GAP}")
    g = _gap(B["lip_seal_tc8x16x5"], hub)
    add("lip seal outer face to hub front (axial)", f"{g:.3f} mm", g >= 0.6,
        lay_checks["lip seal outer face ahead of the hub front (axial gap)"]["value_mm"], ">= 0.6")
    g = _gap(B["shaft"], hub)
    add("shaft to hub D-bore (radial)", f"{g:.3f} mm, 0 mm3", g > 0, 0.075, "> 0")

    # 2. stowed blade vs shells: aero span over the pocket floors, root in the tail notch,
    #    the stowed-detent boss
    g_det = extras["stow_detent"]
    face0 = np.asarray(g_det["face"], float)
    aero, root, boss = [], [], []
    shells_m = {"front_shell": part_manifold(front), "aft_shell": part_manifold(aft)}
    for i, b in enumerate(blades):
        th = extras["blade_theta"][i]
        R = rot_T((0, 0, 0), (1, 0, 0), -th)                       # at_azimuth: -theta about +X
        fc = R[:3, :3] @ face0
        mb = part_manifold(b)
        span = mb.trim_by_plane((-1.0, 0.0, 0.0), -100.0)          # x < 100
        rt = mb.trim_by_plane((1.0, 0.0, 0.0), 100.0) - _sphere(fc, BOSS_ZONE_R)
        bs = mb ^ _sphere(fc, BOSS_ZONE_R)
        aero.append(min(float(shells_m[k].min_gap(span, SEARCH)) for k in shells_m))
        root.append(float(shells_m["aft_shell"].min_gap(rt, SEARCH)))
        boss.append(float(shells_m["aft_shell"].min_gap(bs, SEARCH)))
    add("stowed aero span (x < 100) to the pocket floors and walls, 3 blades",
        f"{min(aero):.3f} mm (blades {', '.join(f'{v:.3f}' for v in aero)})", min(aero) >= P.POCKET_CLEAR - 0.05,
        0.300, f">= {P.POCKET_CLEAR} (-0.05 tessellation)",
        note="real CAD blade vs real shells, 0.02 / 0.03 mm meshes (chord error lowers the reading up to 0.05)")
    add(f"stowed root (x > 100, more than {BOSS_ZONE_R:g} mm from the detent boss face centre) to the static tail",
        f"{min(root):.3f} mm (blades {', '.join(f'{v:.3f}' for v in root)})", min(root) >= P.ROOT_CLEAR - 0.05,
        P.ROOT_CLEAR, f">= {P.ROOT_CLEAR}",
        note="the leaned (49 deg) boss top follows the blade face about 0.09 mm off it out to about 7 mm from its "
             "centre, so that zone is the boss check below, not root clearance")
    add(f"stowed blade to its detent boss (blade within {BOSS_ZONE_R:g} mm of the boss face centre)",
        f"{min(boss):.3f} mm (blades {', '.join(f'{v:.3f}' for v in boss)})", min(boss) > 0.0,
        None, "> 0 (shells design 0.1 to the layout surface)")

    # 3. detent magnet pairs
    sto = extras["blade_seats"]["stowed_seat"]
    a_s = np.asarray(sto["axis_out_blade_frame"], float)
    m_s = np.asarray(sto["mouth_blade_frame_mm"], float)
    T0 = np.asarray(extras["blades"][0]["T_stowed"], float)
    bf = T0[:3, :3] @ (m_s + sto["t_face_mm"] * a_s) + T0[:3, 3]      # blade magnet face centre (world)
    ba = T0[:3, :3] @ a_s
    n = np.asarray(g_det["n"], float)
    tf = np.asarray(g_det["magnet_centre"], float) + MAG_H / 2 * n       # tail magnet face centre
    face_gap = float(np.linalg.norm(bf - tf))
    ang = math.degrees(math.acos(min(1.0, abs(float(ba @ n)))))
    lat = float(np.linalg.norm(np.cross(bf - tf, n)))
    mg = [_gap(B[f"blade_magnet_stow_{i}"], B[f"tail_magnet_{i}"]) for i in range(P.N_BLADES)]
    lr = json.loads((REPO / "data" / "layout_rotor.json").read_text())["detents"]
    f_st = _force(face_gap)
    add("stowed detent magnet pair (blade shank / tail boss): face-centre gap, min gap, axis angle, lateral",
        f"{face_gap:.3f} / {min(mg):.3f} mm, {ang:.2f} deg, {lat:.2f} mm; about {f_st:.2f} N "
        f"(layout {lr['stowed']['force_N']:.2f} N at {lr['stowed']['gap_mm']} mm)",
        face_gap <= 0.45 and ang <= 2.0 and f_st >= 2.0, lr["stowed"]["gap_mm"],
        "owner call T15: face gap <= 0.45 mm, axis angle <= 2.0 deg, force >= 2.0 N",
        note="owner call T15 built: the aft-shell boss leans %.2f deg from radial on the blade face normal (the blade "
             "seat stays where it is); holding torque about %.1f N mm vs the %.2f N mm gravity moment"
             % (float(g_det.get("lean_deg", float("nan"))), f_st * lr["stowed"]["lever_mm"], lr["gravity_moment_max_Nmm"]))
    dg = [_gap(B[f"blade_magnet_heel_{i}"], B[f"hub_magnet_{i}"], 0.0) for i in range(P.N_BLADES)]
    add("deployed detent magnet pair (blade heel / hub stop face): min gap",
        f"{min(dg):.3f} mm, about {_force(min(dg)):.2f} N", abs(min(dg) - lr["deployed"]["gap_mm"]) <= 0.15,
        lr["deployed"]["gap_mm"], "0.4 +- 0.15 (seat depth H + 0.1: 0.3 is slightly stronger)")
    ov = sum(float((part_manifold(b, 0.0) ^ part_manifold(hub)).volume()) for b in blades)
    gh = min(_gap(b, hub, 0.0) for b in blades)
    add("deployed heel on the hub stop face: gap / overlap (3 blades)", f"{gh:.3f} mm / {ov:.3f} mm3",
        gh < 0.02 and ov < DEFECT_MM3, 0.0, "contact, no overlap")

    # 4. folding sweep at the stowing index. Hub, spinner, pins, shells and the pockets are 3-fold
    #    symmetric, so blade 0 stands for all three there (its root, inside a 26 mm ball around
    #    the hinge, against hub / spinner / pins); the fins and pads are not, so all 3 blades go there.
    import manifold3d as m3
    side = m3.Manifold.cube((50.0, 40.0, 40.0)).translate((10.0, 22.0, -20.0))     # fins at theta ~90 / 270
    pads = (part_manifold(front) - hull_outset(0.05)) ^ (side + side.mirror((0.0, 1.0, 0.0)))
    fin_parts = [q for q in parts if q.get("packed") is False and q.get("print") is not None]
    fins = m3.Manifold.batch_boolean([part_manifold(q) for q in fin_parts], m3.OpType.Add)
    pins = m3.Manifold.batch_boolean([part_manifold(B[f"hinge_pin_{i}"]) for i in range(3)], m3.OpType.Add)
    hinge0 = np.asarray(blades[0]["pivot"]["origin"], float)
    root0 = part_manifold(blades[0]) ^ _sphere(hinge0, 26.0)
    sym = {"hub": (part_manifold(hub), True), "spinner": (part_manifold(spin), True),
           "hinge pins": (pins, True), "front_shell": (shells_m["front_shell"], False),
           "aft_shell": (shells_m["aft_shell"], False)}
    far = {"fins": 10.0, "fin pads (front shell proud of the hull)": 6.0}
    sweep = {k: (9e9, None) for k in list(sym) + list(far) + ["neighbour blade"]}

    def upd(k, g, f):
        if g < sweep[k][0]:
            sweep[k] = (g, float(f))

    folds = list(np.arange(0.0, P.FOLD_DEG + 1e-9, 3.0))
    t0 = time.time()
    for f in folds:
        T = [fold_T(b, f) for b in blades]
        m0 = part_manifold(blades[0]).transform(T[0][:3, :4] @ np.linalg.inv(fold_T(blades[0], P.FOLD_DEG)))
        r0 = root0.transform(T[0][:3, :4] @ np.linalg.inv(fold_T(blades[0], P.FOLD_DEG)))
        for k, (t, use_root) in sym.items():
            if k == "hub" and f < 6.0:
                continue                 # heel lifting off its stop face (gap grows from 0): rotor.py rule
            upd(k, float(t.min_gap(r0 if use_root else m0, SEARCH)), f)
        m1 = part_manifold(blades[1]).transform(T[1][:3, :4] @ np.linalg.inv(fold_T(blades[1], P.FOLD_DEG)))
        upd("neighbour blade", float(m0.min_gap(m1, SEARCH)), f)
        for i, b in enumerate(blades):
            mi = m0 if i == 0 else (m1 if i == 1 else
                                    part_manifold(b).transform(T[i][:3, :4] @ np.linalg.inv(fold_T(b, P.FOLD_DEG))))
            upd("fins", float(fins.min_gap(mi, far["fins"])), f)
            k = "fin pads (front shell proud of the hull)"
            if not pads.is_empty():
                upd(k, float(pads.min_gap(mi, far[k])), f)
    print(f"  fold sweep {len(folds)} folds: {time.time() - t0:.0f} s", flush=True)

    def fmt(k):
        lim = far.get(k, SEARCH)
        return (f"{sweep[k][0]:.3f} mm at fold {sweep[k][1]:g}" if sweep[k][0] < lim - 1e-6
                else f">= {lim:g} mm at every fold")

    fin_l = lay_checks["yawed fins and pads to blades folding at the stowing index (3D)"]
    add("blades folding 0..90 deg at the stowing index to the fins", fmt("fins"),
        sweep["fins"][0] >= 2.0, 8.44, ">= 2.0", note="layout: " + fin_l["note"])
    k = "fin pads (front shell proud of the hull)"
    pad_mm3 = float(pads.volume())
    if pad_mm3 > 1.0:
        add("blades folding 0..90 deg at the stowing index to the fin pads", fmt(k),
            sweep[k][0] >= 2.0, fin_l["value_mm"], ">= 2.0")
    else:
        add("front shell material proud of the hull at the fin rails (recessed rail, owner decision)",
            f"{pad_mm3:.3f} mm3 more than 0.05 mm outside the hull", True, None, "flush egg (75 mm stowed)")
    add("blades folding 6..90 deg to the hub (run clearance; below 6 deg the heel lifts off its stop)", fmt("hub"),
        sweep["hub"][0] >= RUN_CLEAR - 0.01, None, f">= RUN_CLEAR {RUN_CLEAR} (-0.01 tessellation)")
    add("blades folding 0..90 deg to the spinner", fmt("spinner"),
        sweep["spinner"][0] >= RUN_CLEAR - 0.01, None, f">= RUN_CLEAR {RUN_CLEAR}")
    add("blades folding 0..90 deg to the hinge pins (knuckle bore)", fmt("hinge pins"),
        sweep["hinge pins"][0] > 0.1, 0.15, "> 0.1 (bore 3.3 on a 3.0 pin)")
    add("blades folding 0..90 deg to the shells (minimum is the stowed boss)",
        f"front {fmt('front_shell')}, aft {fmt('aft_shell')}",
        min(sweep["front_shell"][0], sweep["aft_shell"][0]) > 0.0, None, "> 0")
    add("neighbour blades folding together", fmt("neighbour blade"),
        sweep["neighbour blade"][0] >= 1.0, None, ">= 1.0 (layout rule)")

    # 5. drivetrain
    g = _gap(B["claw_b"], aft)
    add("claw B (rotating) in the static tube mouth", f"{g:.3f} mm", g >= 0.3,
        lay_checks["claw B (rotating) to the static tube mouth, radial"]["value_mm"], ">= 0.3",
        gap_location(part_mesh(B["claw_b"]), part_mesh(aft)))
    g = _gap(B["claw_a"], B["claw_b"])
    add("claw A to claw B dog flanks", f"{g:.3f} mm", g >= 0.25, 0.3, ">= 0.25")
    # owner decision 2026-10-05 (nose round): the display nose is clean (no plug); the gland sits in the v2 shell
    v2s, gl_ = B.get("front_shell_v2"), B["cable_gland_m12"]
    if v2s is not None:
        ov_g = float((part_manifold(gl_) ^ part_manifold(v2s)).volume())
        add("v2 front shell / M12 cable gland: thread core overlap, gap (seated on the spot face)",
            f"{ov_g:.3f} mm3, {_gap(gl_, v2s) if ov_g < 1e-3 else 0.0:.3f} mm", ov_g < DEFECT_MM3, None,
            "no overlap, contact on the spot face")
        ov_d = float((part_manifold(gl_) ^ part_manifold(front)).volume())
        add("display front shell has no gland seat: the gland placed on it would overlap the clean nose by", f"{ov_d:.1f} mm3",
            ov_d > 100.0, None, "> 100 mm3 (the gland is v2 only, never with front_shell)")
    g = float(part_manifold(B["cable_gland_m12"]).min_gap(part_manifold(B["motor_carrier"]), 10.0))
    gl = lay_checks["gland seat block and bore (cable exit) ahead of the carrier plate"]
    add("cable gland to motor carrier", f"{g:.3f} mm" if g < 10.0 else ">= 10 mm", g >= gl["limit_mm"],
        gl["value_mm"], f">= {gl['limit_mm']}", note="layout value is the gland seat block / bore to the carrier "
        "plate; the bought gland body sits further forward")
    gb = min(float(part_manifold(B[a]).min_gap(part_manifold(B[b]), 10.0))
             for a in ("generator_stator", "generator_bell") for b in ("front_shell", "aft_shell"))
    add("generator replica to the shells", f"{gb:.3f} mm", gb >= 0.5, 0.637, ">= 0.5 (bay to shell)")
    ov = float((part_manifold(front) ^ part_manifold(aft)).volume())
    add("front / aft shell locked: overlap, gap", f"{ov:.4f} mm3, {_gap(front, aft):.3f} mm (lugs on the hard stop)",
        ov < DEFECT_MM3, 0.0, "no overlap")

    # 6. closing the shells (recheck R3): the front shell enters at the bayonet entry (rotated by the
    #    twist about X, both senses), is pushed on from 40 mm and twisted to the lock, against the three
    #    blades at a range of folds, without and with the fins fitted
    cs = closing_sweep(parts)
    out.append(dict(check="closing path (bayonet entry, push, twist) vs the blades", value=cs["summary"], ok=True,
                    layout=None, limit="0 mm3 overlap", where=None, note=cs["note"], detail=cs))
    print("  INFO closing path: " + cs["summary"], flush=True)
    return out


CLOSE_FOLDS = (90.0, 85.0, 80.0, 78.0, 75.0, 72.0, 70.0, 65.0, 60.0, 50.0, 45.0, 30.0, 0.0)


def closing_sweep(parts: list[dict], twist: float | None = None) -> dict:
    """Peak overlap (mm3) between the front shell (alone, then with the fins fitted) and the three blades
    while the shells close: entry at +-twist about X, pushed on along X from 40 mm, then twisted to
    the lock in 2.5 deg steps; at each fold of CLOSE_FOLDS (rotor at the stowing index). Reports the
    most folded blade angle that closes with 0 mm3 in both twist senses."""
    import manifold3d as m3
    B = _by(parts)
    if twist is None:
        twist = float(getattr(P, "BAYONET", {}).get("twist_deg", 35.0))
    blades = [B[f"blade_{i}"] for i in range(P.N_BLADES)]
    front = part_manifold(B["front_shell"])
    fins = [q for q in parts if q.get("packed") is False and q.get("rotor") is not True]
    with_fins = m3.Manifold.batch_boolean([front] + [part_manifold(q) for q in fins], m3.OpType.Add)
    res = {"no_fins": {}, "with_fins": {}}
    for label, F in (("no_fins", front), ("with_fins", with_fins)):
        for fold in CLOSE_FOLDS:
            Bu = m3.Manifold.batch_boolean([part_manifold(b, fold) for b in blades], m3.OpType.Add)
            worst = (0.0, None)
            for sg in (1.0, -1.0):
                for dx in (40.0, 25.0, 15.0, 10.0, 6.0, 3.0, 1.0, 0.0):
                    T = rot_T((0, 0, 0), (1, 0, 0), -sg * twist)
                    T[0, 3] -= dx
                    v = float((F.transform(T[:3, :4]) ^ Bu).volume())
                    if v > worst[0]:
                        worst = (v, f"push dx {dx:g}, sense {sg:+g}")
                for a in np.arange(twist, -0.01, -2.5):
                    T = rot_T((0, 0, 0), (1, 0, 0), -sg * a)
                    v = float((F.transform(T[:3, :4]) ^ Bu).volume())
                    if v > worst[0]:
                        worst = (v, f"twist {a:g} deg to lock, sense {sg:+g}")
            res[label][f"{fold:g}"] = dict(peak_overlap_mm3=round(worst[0], 3), at=worst[1])
    ok_n = [f for f in CLOSE_FOLDS if res["no_fins"][f"{f:g}"]["peak_overlap_mm3"] < 0.01]
    ok_f = [f for f in CLOSE_FOLDS if res["with_fins"][f"{f:g}"]["peak_overlap_mm3"] < 0.01]
    # the most folded angle below which every tested fold closes cleanly
    def limit(oks):
        lim = None
        for f in sorted(CLOSE_FOLDS):
            if f in oks:
                lim = f
            else:
                break
        return lim
    ln, lf = limit(ok_n), limit(ok_f)
    res.update(twist_deg=twist, close_fold_max_no_fins=ln, close_fold_max_with_fins=lf,
               fins_in_check=[q["name"] for q in fins])
    res["summary"] = (f"fins off: closes at fold <= {ln:g} deg (unfold >= {P.FOLD_DEG - ln:g} deg); fins on: fold "
                      f"<= {lf:g} deg (unfold >= {P.FOLD_DEG - lf:g} deg)" if ln is not None and lf is not None
                      else f"fins off limit {ln}, fins on limit {lf}")
    res["note"] = ("manifold meshes; overlap > 0.01 mm3 at any push or twist step fails that fold; the fold "
                   "sweep is 0 / 30 / 45 / 50 / 60 / 65 / 70 / 72 / 75 / 78 / 80 / 85 / 90 deg")
    return res


# --------------------------------------------------------------------------
# Rotor spin sweep (deployed): every rotating part vs every static part
# --------------------------------------------------------------------------

SPIN_CONTACT = {("bearing_608", "shaft"), ("bearing_688", "shaft"), ("bearing_608", "shim_8x12x0p5"),
                ("lip_seal_tc8x16x5", "shaft")}     # inner races and the seal lip turn with the shaft


def spin_sweep(parts: list[dict], step: float = 10.0) -> dict:
    S = [p for p in state_parts(parts, "deployed")]
    rot = [p for p in S if p.get("rotor")]
    sta = [p for p in S if not p.get("rotor")]
    pairs = []
    for r in rot:
        for s_ in sta:
            m = pair_metrics(r, s_, 0.0)
            if m is not None and m["min_gap_mm"] < SEARCH:
                pairs.append((r, s_))
    res = []
    for r, s_ in pairs:
        ms = part_manifold(s_, 0.0)
        worst = (9e9, None, 0.0)
        for a in np.arange(0.0, 360.0, step):
            mr = part_manifold(r, 0.0).rotate((float(a), 0.0, 0.0))
            ov = 0.0
            g = float(mr.min_gap(ms, SEARCH))
            if g < 1e-3:
                ov = float((mr ^ ms).volume())
            if g < worst[0] or ov > worst[2]:
                worst = (min(g, worst[0]), float(a), max(ov, worst[2]))
        key = tuple(sorted((r["name"], s_["name"])))
        res.append(dict(rotor=r["name"], static=s_["name"], min_gap_mm=round(worst[0], 4), at_spin_deg=worst[1],
                        max_overlap_mm3=round(worst[2], 4), turns_with_shaft=key in SPIN_CONTACT))
    bad = [x for x in res if not x["turns_with_shaft"] and (x["min_gap_mm"] < 0.1 or x["max_overlap_mm3"] > 0)]
    print(f"  spin sweep: {len(res)} rotor/static pairs within {SEARCH} mm, {len(bad)} rubbing", flush=True)
    return dict(step_deg=step, pairs=res, rubbing=bad)


# --------------------------------------------------------------------------
# Exploded view (finding T29 / polish): one compact axial line in assembly order, applied by the
# viewer as explode * vector in the DEPLOYED pose. Each station is placed behind the previous one
# (its front face GAP mm behind the previous station's aft face, extents measured on the real parts),
# so a geometry change never makes two stations collide. Nested on purpose: the shaft stays in the aft
# shell (its aft end shows behind the tail, with the lip seal on it), screws beside their holes on
# their own axes, pins out of their ears, magnets out of their seats, blades radially out.
# --------------------------------------------------------------------------

# (station, member patterns, gap in front of it mm; None = anchor). Front to aft.
EXPLODE_STATIONS = [
    ("front", ["front_shell*", "fin_*", "rail_magnet_*", "cable_gland_m12"], None),
    ("carrier", ["motor_carrier"], 15.0),           # room for the carrier screws in front of it
    ("spacer", ["adapter_spacer"], 2.5),
    ("stator", ["generator_stator"], 2.5),
    ("bell", ["generator_bell", "inserts_m3_bell"], 2.5),
    ("claw_a", ["claw_a"], 2.5),
    ("claw_b", ["claw_b"], 10.0),                   # room for claw A's screws behind it
    ("collar", ["shaft_collar"], 2.5),
    ("shim", ["shim_*"], 2.5),
    ("b608", ["bearing_608*"], 2.5),
    ("b688", ["bearing_688*"], 2.5),
    ("oring", ["oring_141*"], 3.0),
    ("aft", ["aft_shell*", "tail_magnet_*", "shaft"], 3.0),
    ("hub", ["hub", "hub_magnet_*", "hinge_pin_*", "pin_filler_*", "hub_grub", "blade_*"], -6.0),
    ("spinner", ["spinner"], 3.0),
    ("spinner_screw", ["spinner_screw"], 2.5),
    ("cap", ["spinner_cap"], 2.5),
]
# station whose extent is measured without some of its members (they stick out on purpose)
EXTENT_ONLY = {"aft": ["aft_shell"], "hub": ["hub"], "front": ["front_shell"]}
FRONT_DX = -60.0          # the anchor: front shell forward
BLADE_OUT = 60.0          # blades radially out of the hub slots (review V3: at 34 the deployed blades
                          # were drawn across the aft shell rim)
PIN_OUT = 13.0            # hinge pins along the pin axis out of the ears
FILLER_OUT = 9.0          # pin fillers beyond their pin
HUB_MAGNET_OUT = 16.0     # hub magnets radially out of the stop-face seats (through the open slot)
BLADE_MAGNET_OUT = 6.0    # blade magnets out of their seats along the seat axis
TAIL_MAGNET_OUT = 8.0     # tail magnets out of the bosses along the boss axis
SCREW_OUT = {"screws_carrier_m3x8": -11.0, "screws_claw_a_m3x4": 8.0}   # along X from their station
PLUG_OUT = 28.0           # v2 gland out of the nose seat along the gland axis (V3: 16 barely separated)
LIP_GAP = 3.0             # lip seal this far behind the static tail end, on the protruding shaft


def _x_extent(ps: list[dict], fold: float = 0.0) -> tuple[float, float]:
    lo, hi = np.inf, -np.inf
    for p in ps:
        V = np.asarray(part_mesh(p).vertices)
        T = fold_T(p, fold)
        if T is not None:
            V = V @ T[:3, :3].T + T[:3, 3]
        lo, hi = min(lo, float(V[:, 0].min())), max(hi, float(V[:, 0].max()))
    return lo, hi


def explode_map(parts: list[dict], extras: dict) -> dict:
    """{part name: [dx, dy, dz]} for the viewer's explode slider (deployed pose)."""
    from .. import layout as LY
    member: dict = {}
    for st, pats, _ in EXPLODE_STATIONS:
        for p in parts:
            if p["name"] not in member and p.get("role") != "hidden" and any(fnmatch.fnmatch(p["name"], pt)
                                                                           for pt in pats):
                member[p["name"]] = st
    B = _by(parts)
    dx, prev_hi = {}, None
    for st, pats, gap in EXPLODE_STATIONS:
        ps = [B[n] for n, s_ in member.items() if s_ == st]
        if not ps:
            continue
        keep = EXTENT_ONLY.get(st)
        ext = [p for p in ps if keep is None or p["name"] in keep] or ps
        lo, hi = _x_extent(ext)
        if gap is None:
            d = FRONT_DX
        else:
            d = prev_hi + gap - lo
        if st == "hub":
            # the hub's bore takes the last -gap mm of the shaft end (it sticks out of the tail with the lip
            # seal on it), and stays 2 mm behind the lip seal
            sh_hi = _x_extent([B["shaft"]])[1] if "shaft" in B else lo
            d = sh_hi + dx.get("aft", 0.0) + gap - lo
            if "lip_seal_tc8x16x5" in B and "aft_shell" in B:
                l_lo, l_hi = _x_extent([B["lip_seal_tc8x16x5"]])
                lip_hi = dx["aft"] + _x_extent([B["aft_shell"]])[1] + LIP_GAP + (l_hi - l_lo)
                d = max(d, lip_hi + 2.0 - lo)
        dx[st] = d
        prev_hi = hi + d
        if st == "hub":
            prev_hi = max(prev_hi, _x_extent([p for p in ps if p["name"] in ("hub",)] or ps)[1] + d)
    out = {}
    aft_hi = _x_extent([B["aft_shell"]])[1] if "aft_shell" in B else 0.0
    for p in parts:
        n = p["name"]
        st = member.get(n)
        if st is None and not n.startswith(("screws_", "lip_seal")):
            continue
        v = np.array([dx.get(st, 0.0), 0.0, 0.0])
        if n in SCREW_OUT:
            host = "carrier" if "carrier" in n else "claw_a"
            v = np.array([dx[host] + SCREW_OUT[n], 0.0, 0.0])
        elif n.startswith("lip_seal"):
            lo, _ = _x_extent([B["lip_seal_tc8x16x5"]]) if "lip_seal_tc8x16x5" in B else _x_extent([p])
            v = np.array([dx["aft"] + aft_hi + LIP_GAP - lo, 0.0, 0.0])
        elif n.startswith(("fin_", "rail_magnet_")):
            v += np.asarray(p.get("explode") or [0.0, 0.0, 0.0], float)      # the fins module's slide-out path
        elif n == "cable_gland_m12":
            g = P.GLAND
            t = math.radians(g["theta_deg"])
            a = math.radians(g["axis_angle_from_x_deg"])
            v += PLUG_OUT * np.array([math.cos(a), math.sin(a) * math.sin(t), math.sin(a) * math.cos(t)])
        elif n.startswith("tail_magnet_"):
            th = float(p.get("theta_deg", 0.0))
            v += TAIL_MAGNET_OUT * (rot_T((0, 0, 0), (1, 0, 0), -th)[:3, :3] @ np.asarray(extras["stow_detent"]["n"], float))
        elif n.startswith(("blade_", "hinge_pin_", "pin_filler_", "hub_magnet_")):
            i = int(n.split("_")[-1][0])
            _, e_r, e_t = LY.rotor_axes(120.0 * i)
            if n.startswith(("blade_magnet_heel", "blade_magnet_stow")):
                ax = rot_T((0, 0, 0), *_pivot_axis(p))[:3, :3] @ (np.asarray(extras["blades"][i]["T_stowed"])[:3, :3]
                                                                   @ np.asarray(p["axis_blade_frame"]))
                v += BLADE_OUT * e_r + BLADE_MAGNET_OUT * ax
            elif n.startswith("blade_"):
                v += BLADE_OUT * e_r
            elif n.startswith("hinge_pin_"):
                v += PIN_OUT * e_t
            elif n.startswith("pin_filler_"):
                side = 1.0 if n.endswith("a") else -1.0
                v += PIN_OUT * e_t + side * FILLER_OUT * e_t
            elif n.startswith("hub_magnet_"):
                v += HUB_MAGNET_OUT * e_r
        elif n == "hub_grub":
            psi = math.radians(float(P.ROTOR_HUB["grub_psi_deg"]))
            v += 12.0 * np.array([0.0, -math.sin(psi), math.cos(psi)])
        out[n] = [round(float(c), 3) for c in v]
    for n, st in member.items():                      # any remaining member keeps its station offset
        out.setdefault(n, [round(float(dx.get(st, 0.0)), 3), 0.0, 0.0])
    return out


def exploded_extent(parts: list[dict], extras: dict) -> dict:
    """Axial extent of the exploded deployed scene vs the egg's own length (target <= ~2.2 x)."""
    ex = explode_map(parts, extras)
    S = [p for p in state_parts(parts, "deployed") if p["name"] in ex]
    lo, hi = np.inf, -np.inf
    for p in S:
        V = np.asarray(part_mesh(p).vertices)
        T = fold_T(p, 0.0)
        if T is not None:
            V = V @ T[:3, :3].T + T[:3, 3]
        V = V + np.asarray(ex[p["name"]])
        lo, hi = min(lo, float(V[:, 0].min())), max(hi, float(V[:, 0].max()))
    e_lo, e_hi = _x_extent([p for p in S if p["name"] in ("front_shell", "spinner_cap")])
    return dict(exploded_x_mm=[round(lo, 1), round(hi, 1)], exploded_length_mm=round(hi - lo, 1),
                egg_length_mm=round(e_hi - e_lo, 1), ratio=round((hi - lo) / (e_hi - e_lo), 3))


def _pivot_axis(p: dict):
    """(axis, deployed_deg) of a pivot part: rotating vectors from the stowed to the deployed pose."""
    piv = p["pivot"]
    return piv["axis"], piv["deployed_deg"]


# --------------------------------------------------------------------------
# Web scene (GLB + manifest) and its browser check
# --------------------------------------------------------------------------

_WEBQ: dict = {}
_WEBM: dict = {}


def web_mesh(p: dict):
    """Viewer mesh: absolute chord deflection, sliver-free (export_web.clean_mesh), area-weighted
    normals; its clean-up report goes to _WEBQ[name]."""
    from ..export_web import to_trimesh
    if p["name"] in _WEBM:
        return _WEBM[p["name"]]
    mat = p.get("material", "default")
    tol, ang = WEB_TOL.get(mat if mat in WEB_TOL else "default")
    if p["name"].startswith("blade"):
        tol, ang = WEB_TOL["blade"]
    _clean(p["shape"])
    rep: dict = {}
    m = to_trimesh(p["shape"], tol, ang, relative=False, report=rep)
    _WEBQ[p["name"]] = dict(rep.get("output", {}), input_folded_pairs=rep.get("input", {}).get("folded_pairs_gt170"),
                            flips=rep.get("flips"), collapses=rep.get("collapses"))
    _WEBM[p["name"]] = m
    return m


def glb_quality(path: Path) -> dict:
    """Folded pairs (> 170 deg) per mesh of a written GLB (float32 positions, welded)."""
    import trimesh
    from ..export_web import mesh_quality
    sc = trimesh.load(str(path))
    per = {n: mesh_quality(g) for n, g in sc.geometry.items()}
    bad = {n: q["folded_pairs_gt170"] for n, q in per.items() if q["folded_pairs_gt170"]}
    return dict(meshes=len(per), triangles=int(sum(q["triangles"] for q in per.values())),
                folded_pairs_gt170_total=int(sum(bad.values())), meshes_with_folds=bad,
                min_corner_angle_deg=min(q["min_corner_angle_deg"] for q in per.values()))


SECTION_MATERIAL = {"color": "#C98B5E", "roughness": 0.7, "metalness": 0.0}   # cut faces of the cutaway shells


def section_meshes(parts: list[dict], lift: float = 0.03) -> list[dict]:
    """Finding T28: the planar faces of the cutaway shells that lie in the cut half-planes (theta
    CUTAWAY_THETA) and the front shell's conical cut end, as their own 'section' meshes (lifted
    `lift` mm off the face so they win the depth test), shown with the cutaway shells."""
    import trimesh
    from build123d import GeomType
    from ..export_web import to_trimesh
    from . import shells as SH
    th = [math.radians(t) for t in SH.CUTAWAY_THETA]
    out = []
    for p in parts:
        if p.get("role") != "shell_cutaway" or not p["name"].endswith("shell_cutaway"):
            continue
        meshes = []
        for f in p["shape"].faces():
            c = f.center()
            n = f.normal_at(c)
            keep = False
            if f.geom_type == GeomType.PLANE:
                for t in th:
                    nn = np.array([0.0, math.cos(t), -math.sin(t)])           # normal of the half-plane at t
                    er = np.array([0.0, math.sin(t), math.cos(t)])
                    cc = np.array([c.X, c.Y, c.Z])
                    if abs(abs(float(np.dot([n.X, n.Y, n.Z], nn))) - 1.0) < 1e-6 and abs(float(cc @ nn)) < 1e-3 \
                            and float(cc @ er) > 0:
                        keep = True
            elif f.geom_type == GeomType.CONE and p["name"].startswith("front"):
                keep = 10.0 < c.X < 25.0 and abs(f.area) > 1.0                 # the cut end at x ~16
            if keep:
                m = to_trimesh(f, 0.05, 0.25)
                m.vertices = m.vertices + lift * np.asarray(m.vertex_normals)
                meshes.append(m)
        if meshes:
            out.append(dict(name=p["name"] + "_section", material="section", role="shell_cutaway", rotor=False,
                            module="assembly", mesh=trimesh.util.concatenate(meshes)))
    return out


def write_web(parts: list[dict], extras: dict) -> dict:
    from ..export_web import MATERIALS, write_scene
    ex = explode_map(parts, extras)
    scene = []
    for p in parts:
        if p.get("role") == "hidden":
            continue
        e = {k: p[k] for k in ("name", "material", "role", "rotor", "qty", "bought", "module") if k in p}
        e["explode"] = ex.get(p["name"], [0.0, 0.0, 0.0])
        if p.get("pivot"):
            e["pivot"] = p["pivot"]
        if p.get("packed") is not None:
            e["packed"] = p["packed"]
        if p.get("v2") is not None:
            e["v2"] = p["v2"]
        e["mesh"] = web_mesh(p)
        scene.append(e)
    for e in section_meshes(parts):
        e["explode"] = ex.get(e["name"].replace("_section", ""), [0.0, 0.0, 0.0])
        scene.append(e)
    path = write_scene(scene, OUT_WEB, glb_name="ovo.glb", materials=dict(MATERIALS, section=SECTION_MATERIAL),
                       extra={"generated_by": "ovo.cad.assembly", "cache": extras["cache_dir"],
                              "states": {"stowed": "fold=0", "packed": "fold=0&packed=1", "deployed": "fold=1",
                                         "cutaway": "cutaway=1", "exploded": "explode=1", "v2": "v2=1"},
                              "stowing_index": "spin 0: blade i (psi 120 i) over the pocket at theta -psi"})
    glb = OUT_WEB / "ovo.glb"
    gq = glb_quality(glb)
    return dict(glb=str(glb.relative_to(REPO)), manifest=str(path.relative_to(REPO)),
                glb_mb=round(glb.stat().st_size / 1e6, 2), n_parts=len(scene),
                triangles=int(sum(len(e["mesh"].faces) for e in scene)), explode=ex,
                mesh_quality=dict(glb=gq, per_part={k: v for k, v in _WEBQ.items()}))


def _flat_pose(m) -> np.ndarray:
    """4x4 that lays a loose part on the largest flat facet of its convex hull (the face it rests on),
    long side of the footprint along X."""
    h = m.convex_hull
    n, a = np.asarray(h.face_normals), np.asarray(h.area_faces)
    best, best_a = None, -1.0
    for i in np.argsort(-a):
        tot = float(a[n @ n[i] > 0.99995].sum())         # coplanar hull faces = one facet
        if tot > best_a:
            best, best_a = n[i], tot
    z = np.array([0.0, 0.0, -1.0])
    v = np.cross(best, z)
    s_, c_ = np.linalg.norm(v), float(best @ z)
    if s_ < 1e-12:
        T0 = np.eye(4) if c_ > 0 else rot_T((0, 0, 0), (1, 0, 0), 180.0)
    else:
        T0 = rot_T((0, 0, 0), v / s_, math.degrees(math.atan2(s_, c_)))
    V = np.asarray(m.vertices) @ T0[:3, :3].T
    xy = V[:, :2] - V[:, :2].mean(0)
    w, U = np.linalg.eigh(xy.T @ xy)
    ang = math.atan2(U[1, -1], U[0, -1])                 # main axis of the footprint -> X
    return rot_T((0, 0, 0), (0, 0, 1), -math.degrees(ang)) @ T0


def write_kit_scene(parts: list[dict]) -> dict:
    """Finding T08 / R5: the packed egg (fins off, blanking strips in) with the loose fins (and the magnets
    glued in them) lying flat on the ground in front of the nose, for out/renders/packed_kit.png (a
    separate light scene). Each fin lies on its most stable face, long side along the egg axis, its
    lowest point on the egg's ground plane (the egg's lowest point)."""
    import trimesh
    from ..export_web import MATERIALS, write_scene
    keep = [p for p in state_parts(parts, "packed") if p.get("role") in ("body", "shell")]
    scene = []
    for p in keep:
        scene.append(dict(name=p["name"], material=p.get("material", "default"), role="body", mesh=web_mesh(p)))
    egg = trimesh.util.concatenate([e["mesh"] for e in scene])
    z_ground, x_nose = float(egg.bounds[0][2]), float(egg.bounds[0][0])
    loose = [p for p in parts if p.get("packed") is False and p.get("role") != "hidden"]
    hosts = [p for p in loose if p.get("print") is not None]
    placed = []
    for k, h in enumerate(sorted(hosts, key=lambda q: q["name"])):
        mh = web_mesh(h)
        T = _flat_pose(mh)
        V = np.asarray(mh.vertices) @ T[:3, :3].T + T[:3, 3]
        lo, hi = V.min(0), V.max(0)
        # side by side in front of the nose, 14 mm clear of it and 10 mm apart
        dy = (hi[1] - lo[1] + 10.0) * (k - (len(hosts) - 1) / 2.0)
        sh = np.array([x_nose - 14.0 - hi[0], dy - 0.5 * (lo[1] + hi[1]), z_ground - lo[2]])
        T[:3, 3] += sh
        members = [h] + [q for q in loose if q is not h and q.get("print") is None
                         and abs(float(q.get("theta_deg", 1e9)) - float(h.get("theta_deg", -1e9))) < 1e-6]
        for q in members:
            m = (mh if q is h else web_mesh(q)).copy()
            m.apply_transform(T)
            scene.append(dict(name=q["name"] + "_loose", material=q.get("material", "fin"), role="body", mesh=m))
        placed.append(dict(name=h["name"], lowest_z=round(float(z_ground), 3),
                           footprint_mm=[round(float(hi[0] - lo[0]), 1), round(float(hi[1] - lo[1]), 1)],
                           height_mm=round(float(hi[2] - lo[2]), 1), with_parts=[q["name"] for q in members]))
    path = write_scene(scene, OUT_WEB, glb_name="packed_kit.glb", materials=MATERIALS,
                       extra={"generated_by": "ovo.cad.assembly.write_kit_scene"})
    return dict(manifest=str(path.relative_to(REPO)), n_parts=len(scene), loose=placed,
                ground_z=round(z_ground, 3))


def verify_viewer_fold(parts: list[dict], extras: dict) -> dict:
    """Open the viewer headless at fold=1 and compare each blade's posed vertex box with the
    blade module's deployed solid (blade_world(i, 0)); fold=0 against the stowed shapes."""
    sys.path.insert(0, str(REPO / "scripts"))
    import render as R
    names = [p["name"] for p in parts if p.get("pivot")]
    B = _by(parts)
    ref_dep = {}
    for i, sh in enumerate(extras["deployed_ref"]):
        bb = sh.bounding_box()
        ref_dep[f"blade_{i}"] = np.array([bb.min.X, bb.min.Y, bb.min.Z, bb.max.X, bb.max.Y, bb.max.Z])
    res = {}
    with R.heavy_lock(), R.serve() as base:
        try:
            R.ab("set", "viewport", "800", "500")
            for fold in (1, 0):
                R.ab("open", f"{base}/viewer/index.html?w=800&h=500&dpr=1&fold={fold}")
                R.ab("wait", "#ready")
                raw = R.ab("eval", f"JSON.stringify(window.partBoxes({json.dumps(names)}))")
                boxes = json.loads(json.loads(raw) if raw.startswith('"') else raw)
                rows = []
                for b in boxes:
                    got = np.array(b["min"] + b["max"])
                    n = b["name"]
                    if fold == 1 and n in ref_dep:
                        ref, src = ref_dep[n], "blade.blade_world(i, 0) exact box"
                    else:
                        sh = posed_shape(B[n], 0.0 if fold == 1 else P.FOLD_DEG)
                        bb = sh.bounding_box()
                        ref = np.array([bb.min.X, bb.min.Y, bb.min.Z, bb.max.X, bb.max.Y, bb.max.Z])
                        src = "assembly pose of the CAD solid"
                    rows.append(dict(name=n, max_abs_diff_mm=round(float(np.abs(got - ref).max()), 4), ref=src))
                res[f"fold_{fold}"] = rows
        finally:
            R.ab("close", check=False)
    worst = max(r["max_abs_diff_mm"] for v in res.values() for r in v)
    res["max_abs_diff_mm"] = worst
    res["ok"] = worst < 0.1
    print(f"  viewer fold check: max box difference {worst:.4f} mm", flush=True)
    return res


# --------------------------------------------------------------------------
# STEP assemblies, STL audit, BOM
# --------------------------------------------------------------------------

def export_assembly_steps(parts: list[dict]) -> list[str]:
    out = []
    for st, fold in (("stowed", P.FOLD_DEG), ("deployed", 0.0)):
        kids = []
        for p in state_parts(parts, "packed" if st == "stowed" else st):   # T08: stowed = packed (fins off)
            sh = posed_shape(p, fold)
            sh = sh.moved(Location()) if sh is p["shape"] else sh
            sh.label = p["name"]
            kids.append(sh)
        comp = Compound(children=kids)
        comp.label = f"ovo_{st}"
        f = OUT_STEP / f"ovo_assembly_{st}.step"
        export_step(comp, str(f))
        out.append(str(f.relative_to(REPO)))
        print(f"  wrote {f.relative_to(REPO)} ({f.stat().st_size / 1e6:.1f} MB)", flush=True)
    return out


def stl_path(name: str) -> Path:
    f = OUT_STL / f"{name}.stl"
    if not f.exists() and name.startswith("blade_"):
        f = OUT_STL / "blade.stl"                    # 3 identical blades share one print file
    return f


def stl_audit(parts: list[dict]) -> list[dict]:
    """Every printed part's STL: present, watertight, free of folded sliver pairs (recheck R1: adjacent
    normals more than 170 deg apart); a missing STL, or one with folded pairs, is (re)written from the
    print shape through export_web.write_clean_stl (absolute 0.03 mm deflection). Then volume vs the
    BREP print solid, print size, fits the bed."""
    import trimesh
    from ..export_web import mesh_quality, write_clean_stl
    rows = []
    done: dict = {}
    for p in parts:
        if p.get("print") is None:
            continue
        f = stl_path(p["name"])
        regenerated, q0 = False, None
        if f in done:
            q0, regenerated = done[f]
        elif not f.exists():
            if not write_clean_stl(p["print"], f, tol=0.03, angular=0.2)["watertight"]:
                write_clean_stl(p["print"], f, tol=0.05, angular=0.25)
            regenerated = "missing"
        else:
            q0 = mesh_quality(trimesh.load(f))
            if q0["folded_pairs_gt170"] or not q0["watertight"]:
                old_bytes = f.read_bytes()
                q1 = write_clean_stl(p["print"], f, tol=0.03, angular=0.2)
                if not q1["watertight"]:                        # OCC gap at 0.03: the 0.05 mesh closes
                    q1 = write_clean_stl(p["print"], f, tol=0.05, angular=0.25)
                if (q0["watertight"] and not q1["watertight"]) or q1["folded_pairs_gt170"] > q0["folded_pairs_gt170"]:
                    f.write_bytes(old_bytes)            # the rewrite is not better: keep the module's file
                    regenerated = "kept the module STL (rewrite: folded %d, watertight %s)" % (
                        q1["folded_pairs_gt170"], q1["watertight"])
                else:
                    regenerated = "rewritten (module STL: folded pairs %d, watertight %s)" % (
                        q0["folded_pairs_gt170"], q0["watertight"])
        done[f] = (q0, regenerated)
        m = trimesh.load(f)
        mq = mesh_quality(m)
        bb = m.bounds
        v_brep = float(p["print"].volume) / 1000.0
        v_stl = float(m.volume) / 1000.0
        sh = p["shape"]
        rows.append(dict(
            name=p["name"], module=p["module"], stl=str(f.relative_to(REPO)), regenerated=regenerated,
            module_stl_folded_pairs=None if q0 is None else q0["folded_pairs_gt170"],
            folded_pairs_gt170=mq["folded_pairs_gt170"], min_corner_angle_deg=round(mq["min_corner_angle_deg"], 4),
            valid=bool(sh.is_valid), n_solids=len(sh.solids()), watertight=bool(m.is_watertight),
            winding_consistent=bool(m.is_winding_consistent), stl_volume_cm3=round(float(m.volume) / 1000, 3),
            brep_volume_cm3=round(v_brep, 3), area_cm2=round(float(m.area) / 100, 2),
            volume_match=bool(abs(v_stl - v_brep) <= 0.01 * max(v_brep, 1e-9)),
            print_size_mm=[round(float(v), 2) for v in (bb[1] - bb[0])],
            on_bed=bool(abs(bb[0][2]) < 0.05), fits_bed=bool(C.fits_bed(p["print"])),
            role=p.get("role"), material=p.get("material"), print_with=p.get("print_with"), qty=p.get("qty", 1),
            notes=p.get("notes", "")))
    return rows


def bom(parts: list[dict]) -> list[dict]:
    rows = []
    for p in parts:
        if not p.get("bought") or p.get("role") == "hidden" and p["name"] != "tool_seal_sleeve":
            continue
        rows.append(dict(name=p["name"], qty=p.get("qty", 1), module=p["module"], part_of=p.get("bom_parent"),
                         notes=p.get("notes", "")))
    return rows


# --------------------------------------------------------------------------
# Print plan: material, walls / infill, filament estimate, plates on the 256 bed
# --------------------------------------------------------------------------

DENSITY = {"PLA": 1.24, "PETG": 1.27}          # g/cm3
LINE_W = 0.42                                    # mm, 0.4 nozzle
# display build (v1). Walls and infill follow the module notes where they give them.
PLATE_COUPON = "0 tolerance coupon (PLA; print a PETG copy too)"
PLATE_SHELLS = "1 shells and strips (ivory PLA)"
PLATE_GEN = "5 generator replica (blue PLA)"
PLATE_ROTOR = "2 rotor (PLA, each blade on its cradle)"
PLATE_PETG = "3 drivetrain + fins (PETG)"
PLATE_CUT = "4 cutaway shells (PLA, optional display variant)"
PLATE_V2 = "6 v2 front shell (PETG, v2 wet build only: replaces front_shell)"
PRINT_SETTINGS = {
    "front_shell": dict(mat="PLA", walls=6, infill=0.15, supports="none (38.6 mm nose-plug bridge)", brim=0, plate=PLATE_SHELLS),
    "aft_shell": dict(mat="PLA", walls=6, infill=0.15, supports="none (0.6 mm seal stop ring; detent boss + chin tops "
                      "lean <= 49 deg, no island at 0.05 / 0.1 / 0.15 layer phases; click tines bridge to breakaway "
                      "ties: cut the 3 ties through the slits after printing)", brim=0, plate=PLATE_SHELLS),
    "front_shell_v2": dict(mat="PETG", walls=6, infill=0.15, supports="none (38.6 mm nose-plug bridge; gland seat boss "
                           "with a concave R6 root fillet)", brim=0, optional=True, plate=PLATE_V2),
    "front_shell_cutaway": dict(mat="PLA", walls=6, infill=0.15, supports="none", brim=0, optional=True, plate=PLATE_CUT),
    "aft_shell_cutaway": dict(mat="PLA", walls=6, infill=0.15, supports="none", brim=0, optional=True, plate=PLATE_CUT),
    "hub": dict(mat="PLA", walls=4, infill=0.40, supports="none (bridges <= 14.6 mm)", brim=0, plate=PLATE_ROTOR),
    "spinner": dict(mat="PLA", walls=3, infill=0.30, supports="none (D22.6 seat bridge on a full cone ring, pierced "
                    "by 4 holes: face the clamp ring flat; window tops: rounded, a ~12 mm bridge)", brim=5,
                    plate=PLATE_ROTOR),
    "spinner_cap": dict(mat="PLA", walls=3, infill=1.0, supports="none (45 deg shoulder)", brim=3, plate=PLATE_ROTOR),
    "fin_blank": dict(mat="PLA", walls=3, infill=1.0, supports="none (floor face down, magnet seat on the bed; "
                      "0.08 to 0.12 mm layers for the sloped flush face)", brim=0, plate=PLATE_SHELLS),
    "blade": dict(mat="PLA", walls=4, infill=0.40, supports="none: print on blade_cradle.stl (TE, shank ramp, heel; "
                  "top 0.2 under the blade); fallback tree supports, top Z 0.2", brim=5, plate=PLATE_ROTOR),
    "blade_cradle": dict(mat="PLA", walls=2, infill=0.15, supports="none (it is the support)", brim=0,
                         plate=PLATE_ROTOR),
    "claw_a": dict(mat="PETG", walls=4, infill=0.40, supports="none", brim=0, plate=PLATE_PETG),
    "claw_b": dict(mat="PETG", walls=4, infill=0.40, supports="none (D16.2 bridge on a 0.2 mm membrane, drill it out)",
                   brim=0, plate=PLATE_PETG),
    "motor_carrier": dict(mat="PETG", walls=4, infill=0.40, supports="none", brim=0, plate=PLATE_PETG),
    "adapter_spacer": dict(mat="PETG", walls=3, infill=0.40, supports="none", brim=0, plate=PLATE_PETG),
    "generator_stator": dict(mat="PLA", walls=3, infill=0.15, supports="none", brim=0, plate=PLATE_GEN),
    "generator_bell": dict(mat="PLA", walls=3, infill=0.15, supports="none", brim=0, plate=PLATE_GEN),
    "fin": dict(mat="PETG", walls=3, infill=0.25, supports="none (tip down, tongue up: 73.7 mm on a 3.6 mm wide first "
                "layer: 5 mm brim all round, outer walls <= 50 mm/s, 6-8 s minimum layer time; on an A1 (bed moves in "
                "Y) the chord runs along Y)", brim=5, plate=PLATE_PETG, chord_along_y=True),
    "tolerance_coupon": dict(mat="PLA", walls=3, infill=0.15, supports="none", brim=0, first=True, plate=PLATE_COUPON),
    "tolerance_coupon_key_lug": dict(mat="PLA", walls=3, infill=1.0, supports="none", brim=0, first=True,
                                     plate=PLATE_COUPON),
    "tolerance_coupon_key_tongue": dict(mat="PLA", walls=3, infill=1.0, supports="none", brim=0, first=True,
                                        plate=PLATE_COUPON),
    "tolerance_coupon_key_tongue_up": dict(mat="PLA", walls=3, infill=1.0, supports="none (stands on its 10 mm stub, "
                                           "tongue up, as the fin)", brim=3, first=True, plate=PLATE_COUPON),
    "tolerance_coupon_b": dict(mat="PLA", walls=3, infill=0.15, supports="none", brim=0, first=True, plate=PLATE_COUPON),
}
BED_MARGIN, BED_SPACING = 5.0, 8.0


def _settings(name: str) -> dict:
    for k in (name, name.rstrip("0123456789").rstrip("_")):
        if k in PRINT_SETTINGS:
            return PRINT_SETTINGS[k]
    return dict(mat="PLA", walls=3, infill=0.15, supports="see notes", brim=0)


def _shelf_pack(items: list[dict], bed: float) -> list[list[dict]]:
    """Greedy shelf packing of footprint rectangles (w, d, plus brim) on square plates."""
    plates, cur, x, y, row_d = [], [], BED_MARGIN, BED_MARGIN, 0.0
    usable = bed - BED_MARGIN
    for it in sorted(items, key=lambda r: -r["d"]):
        w, d = it["w"], it["d"]
        if x + w > usable:                                  # next shelf
            x, y, row_d = BED_MARGIN, y + row_d + BED_SPACING, 0.0
        if y + d > usable:                                  # next plate
            plates.append(cur)
            cur, x, y, row_d = [], BED_MARGIN, BED_MARGIN, 0.0
        cur.append(dict(it, x=round(x, 1), y=round(y, 1)))
        x += w + BED_SPACING
        row_d = max(row_d, d)
    if cur:
        plates.append(cur)
    return plates


def print_plan(stl_rows: list[dict]) -> dict:
    """Per printed part: material, orientation note, supports, filament estimate (solid shell of
    walls x line width over the STL area, infill inside), and plates per material."""
    rows = []
    for r in stl_rows:
        st = _settings(r["name"])
        if r.get("print_with"):
            st = dict(st, plate=_settings(r["print_with"]).get("plate"))
        V = r["stl_volume_cm3"] * 1000.0
        A = r["area_cm2"] * 100.0
        w = st["walls"] * LINE_W
        solid = min(V, A * w)
        g_est = DENSITY[st["mat"]] * (solid + st["infill"] * max(0.0, V - solid)) / 1000.0
        n_copies = int(r.get("qty", 1)) if r.get("print_with") else 1       # one cradle per blade
        rows.append(dict(name=r["name"], material=st["mat"], plate=st.get("plate"), walls=st["walls"], infill=st["infill"],
                         print_with=r.get("print_with"),
                         supports=st["supports"], brim_mm=st["brim"], optional=bool(st.get("optional")),
                         first=bool(st.get("first")), print_size_mm=r["print_size_mm"],
                         copies=n_copies, grams_solid=round(n_copies * DENSITY[st["mat"]] * V / 1000.0, 1),
                         grams_est=round(n_copies * g_est, 1),
                         stl=r["stl"], notes=r["notes"]))
    groups = {}
    for r in rows:
        if r.get("print_with"):
            continue                       # printed under its host part (same footprint), not packed apart
        key = _settings(r["name"]).get("plate", r["material"])
        b = 2 * r["brim_mm"]
        sx, sy = r["print_size_mm"][0] + b, r["print_size_mm"][1] + b
        if _settings(r["name"]).get("chord_along_y"):
            # review P4: tall slender part, its long footprint side along Y (an A1's bed moves in Y): packed with
            # its short side along the shelf (x), turned 90 deg when the STL's long side is along X
            groups.setdefault(key, []).append(dict(name=r["name"], w=min(sx, sy), d=max(sx, sy), h=r["print_size_mm"][2],
                                                    rotated=sx > sy, chord_along_y=True))
            continue
        groups.setdefault(key, []).append(dict(name=r["name"], w=max(sx, sy), d=min(sx, sy), h=r["print_size_mm"][2],
                                                rotated=sy > sx))
    plates = []
    for key in sorted(groups):
        packed = _shelf_pack(groups[key], P.BED[0])
        for i, pl in enumerate(packed):
            plates.append(dict(plate=key + (f" #{i + 1}" if len(packed) > 1 else ""), parts=pl,
                               max_height_mm=round(max(it["h"] for it in pl), 1),
                               grams_est=round(sum(next(r["grams_est"] for r in rows if r["name"] == it["name"])
                                                   for it in pl), 1)))
    core = [r for r in rows if not r["optional"] and not r["first"]]
    by_mat = {}
    for r in core:
        by_mat[r["material"]] = round(by_mat.get(r["material"], 0.0) + r["grams_est"], 1)
    return dict(parts=rows, plates=plates, grams_by_material=by_mat,
                total_grams_est=round(sum(by_mat.values()), 1),
                method="ESTIMATE: grams = density x (min(V, A x walls x %.2f mm) + infill x rest); PLA %.2f, "
                       "PETG %.2f g/cm3; excludes supports, brims and purge" % (LINE_W, DENSITY["PLA"], DENSITY["PETG"]))


def assembly_steps(kc: list[dict]) -> list[dict]:
    """drivetrain.assembly_sequence() with step 12 restated with the measured closing limits."""
    from . import drivetrain as DT
    steps = DT.assembly_sequence()
    cs = next((k.get("detail") for k in kc if k["check"].startswith("closing path")), None)
    if cs and cs.get("close_fold_max_no_fins") is not None:
        lf0 = 0.0 if cs["close_fold_max_with_fins"] is None else cs["close_fold_max_with_fins"]
        # instruct one sweep step below each measured limit (the sweep only proves its own grid)
        below = lambda f: max([c for c in CLOSE_FOLDS if c < f] or [0.0])
        ln, lf = below(cs["close_fold_max_no_fins"]), below(lf0)
        for st in steps:
            if st["step"] == 12:
                st["action"] = (f"Take the fins off and unfold every blade by at least {P.FOLD_DEG - ln:g} deg (fold "
                                f"{ln:g} deg or less; with the fins left on, at least {P.FOLD_DEG - lf:g} deg, "
                                f"fold {lf:g} deg or less; the closing sweep is clean from fold "
                                f"{cs['close_fold_max_no_fins']:g} / {lf0:g} down, one sweep step of margin) before "
                                "joining the shells: the stowed blades lie in pockets "
                                "of both shells, 35 deg apart at the bayonet entry. Align the bayonet entry marks and "
                                "push the shells together; the claws find each other on the V lead-ins (turn the rotor "
                                "a little if two dog tips meet), twist 35 deg to the hard stop; only then fold the "
                                "blades home and refit the fins or the blanking strips. Nothing stops the shaft toward "
                                "-X but the hub on the static tail: after stowing, pull the spinner aft so the hub "
                                "runs free.")
                st["check"] = (st.get("check", "") + f" Closing path measured (assembly.closing_sweep): {cs['summary']}.")
    return steps


# --------------------------------------------------------------------------
# Previews for scripts/preview.py
# --------------------------------------------------------------------------

def _preview(state: str) -> list[dict]:
    ps, ex = collect()
    fold = state_fold(state)
    return [dict(name=p["name"], shape=posed_shape(p, fold), material=p.get("material", "default"),
                 role="body") for p in state_parts(ps, state)]


def preview_stowed() -> list[dict]:
    return _preview("stowed")


def preview_deployed() -> list[dict]:
    return _preview("deployed")


def preview_cutaway() -> list[dict]:
    return _preview("cutaway")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def run(force: bool = False, web: bool = True, verify: bool = True, steps: bool = True) -> dict:
    t0 = time.time()
    OUT_CHK.mkdir(parents=True, exist_ok=True)
    ps, ex = collect(force=force)
    print(f"collected {len(ps)} parts ({time.time() - t0:.1f} s)", flush=True)
    rep = dict(module="ovo.cad.assembly", cache=ex["cache_dir"], cache_key=cache_key(),
               build_timing_s=ex.get("timing_s"),
               conventions=dict(stowing_index="rotor spin 0: blade i at psi 120 i lies over the pocket at theta -psi",
                                fold="0 deployed .. 90 stowed (blade.blade_world); viewer fold = 1 - fold/90",
                                states={k: v for k, v in FOLD_STATES.items()},
                                method="manifold3d meshes (0.02 mm chord, shells 0.03); overlaps <= %.1f mm3 on "
                                       "curved contacts are facets" % DEFECT_MM3))
    pivots = []
    for i, b in enumerate(ex["blades"]):
        Tp = rot_T(b["pivot"]["origin"], b["pivot"]["axis"], b["pivot"]["deployed_deg"])
        d = np.abs(Tp @ np.asarray(b["T_stowed"]) - np.asarray(b["T_deployed"])).max()
        pivots.append(dict(blade=i, pivot=b["pivot"], max_matrix_diff=float(d)))
    rep["pivot_proof"] = pivots
    print(f"  pivot proof: max |R(pivot) T_stowed - T_deployed| = {max(x['max_matrix_diff'] for x in pivots):.2e}")
    rep["parts"] = []
    for p in ps:
        if p.get("role") == "hidden" and p.get("bought"):
            continue
        rep["parts"].append(dict(C.check(p["shape"], p["name"]), module=p["module"], role=p.get("role"),
                                 material=p.get("material"), rotor=p.get("rotor"), qty=p.get("qty", 1),
                                 bought=p.get("bought")))
    rep["stl"] = stl_audit(ps)
    rep["print_plan"] = print_plan(rep["stl"])
    rep["exploded_extent"] = exploded_extent(ps, ex)
    print(f"  exploded: {rep['exploded_extent']}", flush=True)
    rep["interference"] = interference(ps)
    rep["cutaway_crossers_r_gt_25"] = cutaway_crossers(ps)
    rep["key_clearances"] = key_clearances(ps, ex)
    rep["spin_sweep"] = spin_sweep(ps)
    rep["bom"] = bom(ps)
    rep["assembly_steps"] = assembly_steps(rep["key_clearances"])
    if steps:
        rep["step_assemblies"] = export_assembly_steps(ps)
    if web:
        w = write_web(ps, ex)
        print(f"  web: {w['glb']} {w['glb_mb']} MB, {w['n_parts']} parts, {w['triangles']} triangles", flush=True)
        rep["web"] = w
        rep["kit_scene"] = write_kit_scene(ps)
        if verify:
            rep["viewer_fold_check"] = verify_viewer_fold(ps, ex)
    rep["elapsed_s"] = round(time.time() - t0, 1)
    (OUT_CHK / "assembly.json").write_text(json.dumps(_jsonable(rep), indent=1))
    n_def = sum(len(v["defects"]) for v in rep["interference"].values())
    n_fail = sum(not k["ok"] for k in rep["key_clearances"])
    print(f"assembly: {len(ps)} parts, {n_def} interference defects, {n_fail} key clearance(s) off design, "
          f"{len(rep['spin_sweep']['rubbing'])} rubbing pairs; wrote out/cad/checks/assembly.json "
          f"({rep['elapsed_s']} s)", flush=True)
    return rep


if __name__ == "__main__":
    a = sys.argv[1:]
    run(force="--force" in a, web="--no-web" not in a, verify="--no-verify" not in a, steps="--no-step" not in a)
