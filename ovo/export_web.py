"""Export build123d shapes to the web viewer, plus the mesh tools every exporter shares.

The viewer gets one GLB with a named mesh per part, and manifest.json with materials, roles,
explode vectors and hinge pivots.

Manifest contract (read by viewer/index.html):
  parts[i] = {
    name:      unique mesh name
    material:  key into materials
    role:      body | shell | shell_cutaway | internal | stand | hidden
               (shell is hidden in cutaway mode, shell_cutaway only shown in it,
                internal shown in cutaway or exploded mode)
    rotor:     true if it spins with the rotor about the X axis
    explode:   [dx, dy, dz] mm at explode = 1 (CAD frame)
    pivot:     optional {origin: [x,y,z], axis: [ax,ay,az], deployed_deg}; the mesh is
               exported in its STOWED pose and the viewer rotates it by
               deployed_deg * fold about the axis through origin (right-hand rule).
    packed:    optional; false = taken off for packing (fins), true = only in the packed egg
               (blanking strips)
    v2:        optional; true = v2 / cutaway hardware (the M12 cable gland), shown only with
               cutaway=1 or v2=1; false = its display-build counterpart (front_shell, clean nose), hidden then
  }

Mesh quality (recheck finding R1). OCC's mesher leaves thousands of zero-area "cap" slivers on the
blade's B-spline skin (one corner near 180 deg, normal arbitrary). trimesh weights vertex normals by
corner ANGLE, so each sliver's 180 deg corner dragged its vertex normal around and the blades
showed streaks. Every exporter now goes through ``clean_mesh`` (weld in float64, collapse needles,
flip caps, until no adjacent pair folds more than 170 deg; the mesh stays watertight) and the
viewer normals are AREA-weighted (``smooth_normals``). ``mesh_quality`` reports the folded pairs.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import trimesh

_SHELL = {"color": "#ECE7DD", "roughness": 0.52, "metalness": 0.0, "clearcoat": 0.15}
MATERIALS = {
    # owner decision 2026-10-05: an all-ivory egg. Shells, hub, spinner, spinner cap, fins and
    # blanking strips share the shell material; the blades stay orange.
    "shell": dict(_SHELL),
    "hub": dict(_SHELL),
    "fin": dict(_SHELL),
    "blade": {"color": "#E8622C", "roughness": 0.45, "metalness": 0.0, "clearcoat": 0.25},
    "drive": {"color": "#2F3337", "roughness": 0.5, "metalness": 0.0},      # printed drivetrain internals
    "steel": {"color": "#C9CDD2", "roughness": 0.28, "metalness": 1.0},
    "rubber": {"color": "#1A1A1A", "roughness": 0.8, "metalness": 0.0},
    "magnet": {"color": "#A7ADB4", "roughness": 0.3, "metalness": 1.0},
    "generator": {"color": "#3A6EA5", "roughness": 0.45, "metalness": 0.2},
    "stand": {"color": "#8A8378", "roughness": 0.7, "metalness": 0.0},
    "default": {"color": "#9A9A9A", "roughness": 0.6, "metalness": 0.0},
}


# --------------------------------------------------------------------------
# Triangulation
# --------------------------------------------------------------------------

def occ_mesh(shape, tol: float = 0.05, angular: float = 0.3, relative: bool = False) -> trimesh.Trimesh:
    """Triangulate a build123d shape with OCC's BRepMesh and return a welded float64 trimesh.

    relative=False (default) makes ``tol`` an ABSOLUTE chord deflection in mm. build123d's own
    ``tessellate`` and ``export_stl`` pass isRelative=True, where the deflection scales with each
    edge's length, so on long faces only the angular limit applies. Any triangulation the shape
    already carries is dropped first, so the requested tolerance is the one used."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.BRepTools import BRepTools
    from OCP.IMeshTools import IMeshTools_Parameters
    from OCP.TopAbs import TopAbs_REVERSED
    from OCP.TopLoc import TopLoc_Location
    BRepTools.Clean_s(shape.wrapped)
    prm = IMeshTools_Parameters()
    prm.Deflection, prm.Angle, prm.Relative, prm.InParallel = float(tol), float(angular), bool(relative), True
    BRepMesh_IncrementalMesh(shape.wrapped, prm)
    V, F, off = [], [], 0
    for f in shape.faces():
        loc = TopLoc_Location()
        t = BRep_Tool.Triangulation_s(f.wrapped, loc)
        if t is None:
            continue
        tr = loc.Transformation()
        nodes = np.array([(q.X(), q.Y(), q.Z()) for q in (t.Node(i).Transformed(tr)
                                                          for i in range(1, t.NbNodes() + 1))], float)
        tris = np.array([t.Triangle(i).Get() for i in range(1, t.NbTriangles() + 1)], np.int64) - 1
        if f.wrapped.Orientation() == TopAbs_REVERSED:
            tris = tris[:, [0, 2, 1]]
        V.append(nodes)
        F.append(tris + off)
        off += len(nodes)
    BRepTools.Clean_s(shape.wrapped)
    return trimesh.Trimesh(np.vstack(V), np.vstack(F), process=True)


# --------------------------------------------------------------------------
# Mesh quality and sliver removal
# --------------------------------------------------------------------------

def mesh_quality(m: trimesh.Trimesh, fold_deg: float = 170.0) -> dict:
    """Folded adjacent pairs (normals more than fold_deg apart) on the welded mesh, plus the
    smallest face area and corner angle, watertightness and winding."""
    w = trimesh.Trimesh(np.asarray(m.vertices, float), np.asarray(m.faces), process=True)
    a = np.degrees(w.face_adjacency_angles) if len(w.face_adjacency) else np.zeros(0)
    ang = np.degrees(trimesh.triangles.angles(w.triangles)) if len(w.faces) else np.zeros((0, 3))
    return dict(triangles=int(len(w.faces)), folded_pairs_gt170=int((a > fold_deg).sum()),
                pairs_gt120=int((a > 120.0).sum()), min_face_area_mm2=float(w.area_faces.min()) if len(w.faces) else 0.0,
                min_corner_angle_deg=float(ang.min()) if len(ang) else 0.0,
                max_corner_angle_deg=float(ang.max()) if len(ang) else 0.0,
                watertight=bool(w.is_watertight), winding_consistent=bool(w.is_winding_consistent))


def _vertex_faces(F: np.ndarray, nv: int):
    """CSR vertex -> faces lookup: (order, start) so faces of v are order[start[v]:start[v+1]]."""
    flat = F.ravel()
    order = np.argsort(flat, kind="stable")
    start = np.searchsorted(flat[order], np.arange(nv + 1))
    return order // 3, start


def _normal(p0, p1, p2):
    n = np.cross(p1 - p0, p2 - p0)
    ln = np.linalg.norm(n)
    return n / ln if ln > 0 else n


def _min_angle(P3: np.ndarray) -> float:
    return float(np.degrees(trimesh.triangles.angles(P3[None]))[0].min())


def clean_mesh(m: trimesh.Trimesh, simplify: float | None = None, **kw) -> tuple[trimesh.Trimesh, dict]:
    """Sliver removal (``_clean_mesh_once``); when folded pairs survive, retry after a manifold3d
    simplify at 0.002 / 0.005 / 0.01 / 0.02 mm (every surface moves less than that) and keep the
    first fold-free watertight result, else the best one. Returns (mesh, report)."""
    out, rep = _clean_mesh_once(m, simplify=simplify, **kw)
    if rep["output"]["folded_pairs_gt170"] == 0:
        return out, rep
    best = (out, rep)
    for eps in (0.002, 0.005, 0.01, 0.02):
        if simplify and eps <= simplify:
            continue
        o2, r2 = _clean_mesh_once(m, simplify=eps, **kw)
        q2, qb = r2["output"], best[1]["output"]
        if (q2["watertight"] or not rep["input"]["watertight"]) and q2["folded_pairs_gt170"] < qb["folded_pairs_gt170"]:
            best = (o2, dict(r2, retried_with_simplify_mm=eps))
        if best[1]["output"]["folded_pairs_gt170"] == 0:
            break
    return best


def _clean_mesh_once(m: trimesh.Trimesh, fold_deg: float = 170.0, cap_deg: float = 170.0, min_edge: float = 1e-3,
                     max_pass: int = 40, simplify: float | None = None) -> tuple[trimesh.Trimesh, dict]:
    """Remove sliver triangles from a closed mesh without opening it.

    Welds in float64, drops faces with a repeated vertex, then repeats passes of independent local
    operations on the faces that take part in a folded pair (adjacent normals more than fold_deg
    apart), that have a corner above cap_deg or an edge shorter than min_edge:
      * cap (largest corner > 120 deg, the vertex almost on the opposite edge): flip that long
        edge, only if both new faces face the same way as the pair's larger face and the smallest
        corner of the pair grows (no ping-pong);
      * needle (otherwise): collapse its shortest edge to its midpoint, only when the one-ring link
        condition holds and no neighbouring face flips.
    simplify (mm): first run manifold3d's simplify at that tolerance (keeps a subset of the original
    vertices, every surface moves less than it), which removes most slivers and many triangles.
    The result is kept only if it has no more folded pairs than the input and stays watertight when
    the input was. Returns (mesh, report)."""
    V = np.asarray(m.vertices, float).copy()
    F = np.asarray(m.faces, np.int64).copy()
    w = trimesh.Trimesh(V, F, process=True)
    V, F = np.asarray(w.vertices, float).copy(), np.asarray(w.faces, np.int64).copy()
    F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])]
    base = trimesh.Trimesh(V.copy(), F.copy(), process=False)
    rep = dict(input=mesh_quality(base), flips=0, collapses=0, passes=0, simplified=False, kept_input=False)
    if simplify:
        try:
            import manifold3d as m3
            mm = m3.Mesh64(vert_properties=np.ascontiguousarray(V), tri_verts=np.ascontiguousarray(F, np.uint64))
            mf = m3.Manifold(mm)
            if mf.status() == m3.Error.NoError and not mf.is_empty():
                o = mf.simplify(float(simplify)).to_mesh64()
                V2 = np.asarray(o.vert_properties)[:, :3].astype(float)
                F2 = np.asarray(o.tri_verts).astype(np.int64)
                if len(F2):
                    V, F = V2.copy(), F2.copy()
                    rep["simplified"] = True
        except Exception as e:  # noqa: BLE001 (keep the unsimplified mesh)
            rep["simplify_error"] = f"{type(e).__name__}: {e}"
    hist = []
    for it in range(max_pass):
        tm = trimesh.Trimesh(V, F, process=False)
        adj = tm.face_adjacency
        fa = np.degrees(tm.face_adjacency_angles) if len(adj) else np.zeros(0)
        tri = V[F]
        A = np.degrees(trimesh.triangles.angles(tri))
        E = np.stack([np.linalg.norm(tri[:, 1] - tri[:, 2], axis=1), np.linalg.norm(tri[:, 2] - tri[:, 0], axis=1),
                      np.linalg.norm(tri[:, 0] - tri[:, 1], axis=1)], axis=1)        # edge opposite corner k
        area = tm.area_faces
        fnorm = tm.face_normals
        folds = adj[fa > fold_deg]
        hist.append(len(folds))
        bad = set(np.nonzero((A.max(1) > cap_deg) | (E.min(1) < min_edge))[0].tolist())
        if len(folds):
            f_, g_ = folds[:, 0], folds[:, 1]
            q = area / np.maximum(E.max(1) ** 2, 1e-30)
            bad.update(np.where(q[f_] <= q[g_], f_, g_).tolist())
        if not bad or (len(hist) > 4 and hist[-1] >= hist[-4] and hist[-1] > 0 and it > 8):
            break
        rep["passes"] = it + 1
        order, start = _vertex_faces(F, len(V))
        order, start = order.tolist(), start.tolist()
        Fl = F.tolist()
        Al, El, arl = A.tolist(), E.tolist(), area.tolist()
        alive = [True] * len(Fl)
        locked = bytearray(len(V))
        vfaces = lambda v: [h for h in order[start[v]:start[v + 1]] if alive[h]]
        changed = False
        for f in sorted(bad, key=lambda k: arl[k]):
            fv = Fl[f]
            if not alive[f] or locked[fv[0]] or locked[fv[1]] or locked[fv[2]]:
                continue
            Af, Ef = Al[f], El[f]
            k_cap = int(np.argmax(Af))
            if Af[k_cap] > 120.0 and min(Ef) >= min_edge:
                # flip the long edge (a, b) opposite the cap corner c
                c, a, b = fv[k_cap], fv[(k_cap + 1) % 3], fv[(k_cap + 2) % 3]
                g = [h for h in vfaces(a) if h != f and b in Fl[h]]
                if len(g) != 1:
                    continue
                g = g[0]
                gv = Fl[g]
                d = [v for v in gv if v not in (a, b)]
                if len(d) != 1 or any(locked[v] for v in gv):
                    continue
                d = d[0]
                if any(d in Fl[h] for h in vfaces(c)):
                    continue                                   # edge c-d exists already
                nref = fnorm[g] if arl[g] >= arl[f] else fnorm[f]
                f1, f2 = [a, d, c], [d, b, c]
                n1, n2 = _normal(*V[f1]), _normal(*V[f2])
                if nref @ n1 <= 0.0 or nref @ n2 <= 0.0:
                    continue
                if min(_min_angle(V[f1]), _min_angle(V[f2])) <= min(min(Af), min(Al[g])) + 1e-9:
                    continue
                Fl[f], Fl[g] = f1, f2
                for v in (a, b, c, d):
                    locked[v] = 1
                rep["flips"] += 1
                changed = True
            else:
                # collapse the shortest edge (a, b) into its midpoint
                k = int(np.argmin(Ef))
                a, b = fv[(k + 1) % 3], fv[(k + 2) % 3]
                fa_, fb_ = vfaces(a), vfaces(b)
                shared = [h for h in fa_ if h in fb_]
                if len(shared) != 2:
                    continue
                ring = set(fa_) | set(fb_)
                ring_v = {v for h in ring for v in Fl[h]}
                if any(locked[v] for v in ring_v):
                    continue
                na = {v for h in fa_ for v in Fl[h]} - {a}
                nb = {v for h in fb_ for v in Fl[h]} - {b}
                opp = {v for h in shared for v in Fl[h]} - {a, b}
                if (na & nb) != opp:
                    continue                                   # link condition
                pm = 0.5 * (V[a] + V[b])
                ok = True
                for h in ring:
                    if h in shared:
                        continue
                    P1 = np.array([pm if v in (a, b) else V[v] for v in Fl[h]])
                    n1 = _normal(*P1)
                    if not n1.any() or (arl[h] > 1e-6 and fnorm[h] @ n1 < 0.3):
                        ok = False
                        break
                if not ok:
                    continue
                V[a] = pm
                for h in fb_:
                    Fl[h] = [a if v == b else v for v in Fl[h]]
                for h in shared:
                    alive[h] = False
                for v in ring_v:
                    locked[v] = 1
                rep["collapses"] += 1
                changed = True
        F = np.asarray([fl for fl, al in zip(Fl, alive) if al], np.int64)
        if not changed:
            break
    rep["folded_pairs_per_pass"] = hist
    out = trimesh.Trimesh(V, F, process=False)
    out.remove_unreferenced_vertices()
    q = mesh_quality(out)
    if q["folded_pairs_gt170"] > rep["input"]["folded_pairs_gt170"] or \
            (rep["input"]["watertight"] and not q["watertight"]):
        out, q, rep["kept_input"] = base, rep["input"], True
    rep["output"] = q
    return out, rep


def smooth_normals(m: trimesh.Trimesh, smooth_deg: float = 28.0) -> trimesh.Trimesh:
    """Split the mesh where adjacent faces meet at more than smooth_deg (crisp trailing edges and
    chamfers), then give every vertex the AREA-weighted mean of its faces' normals."""
    s = trimesh.graph.smooth_shade(m, angle=np.radians(smooth_deg))
    fn = np.asarray(s.face_normals) * np.asarray(s.area_faces)[:, None]
    vn = np.zeros_like(np.asarray(s.vertices, float))
    for k in range(3):
        np.add.at(vn, s.faces[:, k], fn)
    ln = np.linalg.norm(vn, axis=1)
    bad = ln < 1e-20
    vn[~bad] /= ln[~bad, None]
    if bad.any():
        vn[bad] = np.asarray(s.vertex_normals)[bad]
    s.vertex_normals = vn
    return s


def to_trimesh(shape, tol: float = 0.02, angular: float = 0.12, smooth_deg: float = 28.0,
               relative: bool = True, simplify: float | None = None, report: dict | None = None) -> trimesh.Trimesh:
    """Tessellate a build123d shape for the viewer: sliver-free (``clean_mesh``), normals smooth
    across edges flatter than smooth_deg and split at sharper ones, area-weighted (so curved hulls
    shade smoothly while trailing edges and chamfers stay crisp). relative=True keeps the
    build123d-style deflection (scaled by edge length, the angular limit governs); when folded
    pairs survive the clean-up it re-meshes with the absolute deflection and keeps the better one.
    ``report`` (a dict) receives the clean_mesh report."""
    m = occ_mesh(shape, tol, angular, relative=relative)
    c, rep = clean_mesh(m, simplify=simplify)
    if relative and rep["output"]["folded_pairs_gt170"]:
        c2, rep2 = clean_mesh(occ_mesh(shape, tol, angular, relative=False), simplify=simplify)
        o, o2 = rep["output"], rep2["output"]
        if o2["folded_pairs_gt170"] < o["folded_pairs_gt170"] and (o2["watertight"] or not o["watertight"]):
            c, rep = c2, dict(rep2, fallback="absolute deflection")
    if report is not None:
        report.update(rep)
    return smooth_normals(c, smooth_deg)


def write_clean_stl(shape, path, tol: float = 0.05, angular: float = 0.3, simplify: float | None = 1e-3) -> dict:
    """Print STL from an absolute-deflection mesh, sliver-free and watertight; returns the
    mesh_quality report (plus the clean_mesh counts)."""
    m = occ_mesh(shape, tol, angular, relative=False)
    m, rep = clean_mesh(m, simplify=simplify)
    if not m.is_watertight:                 # OCC sometimes leaves sliver gaps where blends meet
        from .cad.common import close_mesh_slivers
        close_mesh_slivers(m, 12.0)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    m.export(str(path))
    q = mesh_quality(trimesh.load(str(path)))
    q.update(sliver_flips=rep["flips"], sliver_collapses=rep["collapses"], input_folded_pairs=rep["input"]["folded_pairs_gt170"],
             simplified_mm=simplify if rep["simplified"] else None, deflection_mm=tol, angular_rad=angular)
    return q


# --------------------------------------------------------------------------
# Scene
# --------------------------------------------------------------------------

def write_scene(parts: list[dict], out_dir: str | Path, glb_name: str = "ovo.glb",
                materials: dict | None = None, extra: dict | None = None) -> Path:
    """parts: [{name, shape | mesh, material, role?, rotor?, explode?, pivot?, packed?, v2?}, ...]"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scene = trimesh.Scene()
    manifest_parts = []
    for p in parts:
        mesh = p["mesh"] if "mesh" in p else to_trimesh(p["shape"])
        scene.add_geometry(mesh, node_name=p["name"], geom_name=p["name"])
        entry = {k: v for k, v in p.items() if k not in ("shape", "mesh")}
        entry.setdefault("role", "body")
        entry.setdefault("material", "default")
        manifest_parts.append(entry)
    scene.export(out_dir / glb_name, file_type="glb", include_normals=True)
    manifest = {"units": "mm", "up": "Z", "glb": glb_name,
                "materials": materials or MATERIALS, "parts": manifest_parts, **(extra or {})}
    path = out_dir / ("manifest.json" if glb_name == "ovo.glb" else glb_name.replace(".glb", ".json"))
    path.write_text(json.dumps(manifest, indent=2))
    return path
