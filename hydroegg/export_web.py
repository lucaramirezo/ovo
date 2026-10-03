"""Export build123d shapes to the web viewer: one GLB with a named mesh per part
plus manifest.json with materials, roles, explode vectors and hinge pivots.

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
  }
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import trimesh

MATERIALS = {
    "shell": {"color": "#ECE7DD", "roughness": 0.52, "metalness": 0.0, "clearcoat": 0.15},
    "blade": {"color": "#E8622C", "roughness": 0.45, "metalness": 0.0, "clearcoat": 0.25},
    "hub": {"color": "#2F3337", "roughness": 0.5, "metalness": 0.0},
    "fin": {"color": "#2F3337", "roughness": 0.5, "metalness": 0.0},
    "steel": {"color": "#C9CDD2", "roughness": 0.28, "metalness": 1.0},
    "rubber": {"color": "#1A1A1A", "roughness": 0.8, "metalness": 0.0},
    "magnet": {"color": "#A7ADB4", "roughness": 0.3, "metalness": 1.0},
    "generator": {"color": "#3A6EA5", "roughness": 0.45, "metalness": 0.2},
    "stand": {"color": "#8A8378", "roughness": 0.7, "metalness": 0.0},
    "default": {"color": "#9A9A9A", "roughness": 0.6, "metalness": 0.0},
}


def to_trimesh(shape, tol: float = 0.02, angular: float = 0.12, smooth_deg: float = 28.0) -> trimesh.Trimesh:
    """Tessellate a build123d shape. Normals are smooth across edges flatter than
    smooth_deg and split at sharper ones, so curved hulls shade smoothly while
    trailing edges and chamfers stay crisp."""
    verts, faces = shape.tessellate(tol, angular)
    mesh = trimesh.Trimesh(
        np.array([(v.X, v.Y, v.Z) for v in verts], dtype=np.float64),
        np.array(faces, dtype=np.int64),
        process=True,
    )
    mesh = trimesh.graph.smooth_shade(mesh, angle=np.radians(smooth_deg))
    _ = mesh.vertex_normals  # populate the cache so the GLB carries normals
    return mesh


def write_scene(parts: list[dict], out_dir: str | Path, glb_name: str = "hydroegg.glb",
                materials: dict | None = None, extra: dict | None = None) -> Path:
    """parts: [{name, shape, material, role?, rotor?, explode?, pivot?}, ...]"""
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
    path = out_dir / ("manifest.json" if glb_name == "hydroegg.glb" else glb_name.replace(".glb", ".json"))
    path.write_text(json.dumps(manifest, indent=2))
    return path
