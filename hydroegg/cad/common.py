"""Shared CAD helpers for every part module.

Conventions (see hydroegg/params.py): mm, X = body axis with the flow toward +X, Z up.
Every part module exposes ``parts() -> list[dict]`` with one dict per printed or bought
part, in the format consumed by hydroegg/cad/assembly.py and export_web.write_scene:

    name       unique snake_case id (blade_0, front_shell, bearing_608, ...)
    shape      build123d Shape in WORLD coordinates, assembled pose (blades: STOWED pose)
    print      Shape already oriented for printing, resting on z = 0 (None for bought parts)
    material   key of export_web.MATERIALS (shell, blade, hub, fin, steel, rubber, magnet, ...)
    role       body | shell | shell_cutaway | internal | stand | hidden
    rotor      True if it spins with the rotor
    qty        how many of this exact part the build needs (default 1)
    bought     True for purchased parts (BOM, never printed)
    explode    [dx, dy, dz] mm offset for the exploded view (assembly.py may override)
    pivot      blades only: {origin, axis, deployed_deg} so that rotating the STOWED shape by
               deployed_deg about axis through origin gives the deployed pose
    notes      short print notes (orientation, supports, material)
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np
from build123d import (
    Axis,
    BuildLine,
    BuildPart,
    BuildSketch,
    Cylinder,
    Location,
    Plane,
    Polyline,
    Rotation,
    Shape,
    Solid,
    Spline,
    Vector,
    Line,
    make_face,
    revolve,
)

from .. import params as P
from .. import profile


# --------------------------------------------------------------------------
# Bodies of revolution about the X axis
# --------------------------------------------------------------------------

def revolve_xr(pts: Sequence[Sequence[float]], spline: bool = True) -> Solid:
    """Revolve a closed (x, r) outline (r >= 0, drawn in the XZ half-plane) about X.

    pts is the outer edge from one end to the other; the outline is closed back along
    the axis (r = 0) when the first and last points sit on it, otherwise by straight
    segments down to the axis."""
    pts = [(float(x), float(r)) for x, r in pts]
    with BuildPart() as bp:
        with BuildSketch(Plane.XZ):
            with BuildLine():
                if spline:
                    Spline(*pts)
                else:
                    Polyline(*pts)
                tail = []
                if pts[-1][1] > 1e-9:
                    tail.append((pts[-1][0], 0.0))
                if pts[0][1] > 1e-9:
                    tail.append((pts[0][0], 0.0))
                Polyline(pts[-1], *tail, pts[0])
            make_face()
        revolve(axis=Axis.X)
    return bp.part


def hull_xr(x0: float, x1: float, n: int = 160, inset: float = 0.0) -> np.ndarray:
    """(x, r) points of the hull between x0 and x1; inset > 0 moves the curve inward
    along the surface normal (inner surface of a shell of that wall)."""
    s = (1 - np.cos(np.linspace(0, math.pi, n))) / 2
    xs = x0 + (x1 - x0) * s
    rs = profile.radius(xs) if inset == 0 else profile.inner_radius(xs, wall=inset)
    pts = np.column_stack([xs, np.clip(rs, 0, None)])
    # An inset curve hits the axis before the hull ends (r = 0 runs near the nose and
    # the tail). Keep a single on-axis point at each end, or the spline doubles back
    # along the axis and the outline no longer closes.
    pos = np.nonzero(pts[:, 1] > 1e-6)[0]
    if len(pos) and (pos[0] > 1 or pos[-1] < len(pts) - 2):
        pts = pts[max(pos[0] - 1, 0): min(pos[-1] + 2, len(pts))]
    return pts


def hull_solid(x0: float = 0.0, x1: float | None = None, inset: float = 0.0, n: int = 160) -> Solid:
    """Solid of revolution bounded by the hull (or an inset copy of it) from x0 to x1,
    capped flat at both ends."""
    x1 = P.HULL_L if x1 is None else x1
    pts = hull_xr(x0, x1, n, inset)
    pts[0, 1] = max(pts[0, 1], 0.0)
    return revolve_xr(pts)


def ring_xr(x: Sequence[float], r_out: Sequence[float], r_in: Sequence[float]) -> Solid:
    """Axisymmetric ring from piecewise-linear outer and inner radius tables over x
    (straight segments, good for tubes, seats, spigots and grooves)."""
    outer = list(zip(x, r_out))
    inner = list(zip(x, r_in))[::-1]
    with BuildPart() as bp:
        with BuildSketch(Plane.XZ):
            with BuildLine():
                Polyline(*outer, *inner, outer[0])
            make_face()
        revolve(axis=Axis.X)
    return bp.part


# --------------------------------------------------------------------------
# Primitives placed between points
# --------------------------------------------------------------------------

def cyl(p0: Sequence[float], p1: Sequence[float], r: float) -> Solid:
    """Cylinder of radius r from p0 to p1 (any direction)."""
    a, b = Vector(*p0), Vector(*p1)
    d = b - a
    h = d.length
    pl = Plane(origin=a, z_dir=d.normalized())
    return pl * Cylinder(r, h, align=None)


def axis_cyl(x0: float, x1: float, r: float) -> Solid:
    """Cylinder on the body axis from x0 to x1."""
    return cyl((x0, 0, 0), (x1, 0, 0), r)


def at_azimuth(shape: Shape, theta_deg: float) -> Shape:
    """Carry a shape built at theta = 0 (on the +Z side) to azimuth theta, where theta is
    measured about the body axis FROM +Z TOWARD +Y (a point at theta sits at
    y = r sin(theta), z = r cos(theta)). That is a rotation of -theta about +X.
    Layout tables state their own azimuth convention (e.g. pockets: theta = -psi):
    convert explicitly at the call site and assert a sample point."""
    return shape.rotate(Axis.X, -theta_deg)


# --------------------------------------------------------------------------
# Airfoils
# --------------------------------------------------------------------------

def naca4(m: float, p: float, t: float, n: int = 80, te: float = 0.0) -> np.ndarray:
    """Closed NACA 4-digit loop for a UNIT chord: m camber (fraction), p camber position
    (fraction), t thickness (fraction). te adds a linear thickness ramp so the trailing
    edge is te thick (as a fraction of chord), keeping the section printable.

    Returns (2n-1, 2) points: from the trailing edge along the upper surface to the
    leading edge and back along the lower surface (first and last points differ by the
    blunt trailing edge, so close the loop with a straight segment)."""
    beta = np.linspace(0, math.pi, n)
    x = (1 - np.cos(beta)) / 2
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    yt = yt + 0.5 * te * x                         # blunt TE ramp, total te at x = 1
    if m > 0 and p > 0:
        yc = np.where(x < p, m / p**2 * (2 * p * x - x**2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))
        dyc = np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    upper = np.column_stack([xu, yu])[::-1]        # TE -> LE
    lower = np.column_stack([xl, yl])[1:]          # LE -> TE (skip duplicate LE)
    return np.vstack([upper, lower])


def naca_code(code: str) -> tuple[float, float, float]:
    """'4412' -> (0.04, 0.4, 0.12)."""
    code = code.strip().lower().replace("naca", "")
    return int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100


# --------------------------------------------------------------------------
# Checks and print orientation
# --------------------------------------------------------------------------

def check(shape: Shape, name: str = "") -> dict:
    """Validity summary: valid flag, solid count, volume (cm3), bounding box (mm)."""
    bb = shape.bounding_box()
    solids = shape.solids() if hasattr(shape, "solids") else []
    return {
        "name": name,
        "valid": bool(shape.is_valid),
        "n_solids": len(solids),
        "volume_cm3": round(shape.volume / 1000.0, 3),
        "bbox_mm": [round(v, 2) for v in (bb.min.X, bb.min.Y, bb.min.Z, bb.max.X, bb.max.Y, bb.max.Z)],
        "size_mm": [round(v, 2) for v in (bb.size.X, bb.size.Y, bb.size.Z)],
    }


def on_bed(shape: Shape, rot: Iterable[float] = (0, 0, 0)) -> Shape:
    """Rotate by (rx, ry, rz) degrees about the world axes (in that order), then move so
    the part rests on z = 0 centred on x = y = 0. Use for the 'print' entry of a part."""
    rx, ry, rz = rot
    s = shape
    if rx:
        s = s.rotate(Axis.X, rx)
    if ry:
        s = s.rotate(Axis.Y, ry)
    if rz:
        s = s.rotate(Axis.Z, rz)
    bb = s.bounding_box()
    return s.moved(Location((-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z)))


def fits_bed(shape: Shape, margin: float = 3.0) -> bool:
    bb = shape.bounding_box()
    return bb.size.X <= P.BED[0] - 2 * margin and bb.size.Y <= P.BED[1] - 2 * margin and bb.size.Z <= P.BED[2]


def interference(a: Shape, b: Shape) -> float:
    """Overlap volume (mm3) between two shapes; 0 when they do not intersect."""
    try:
        common = a & b
        return float(common.volume) if common is not None else 0.0
    except Exception:  # noqa: BLE001 (OCC boolean failure: report as NaN, never as 0)
        return float("nan")


def min_distance(a: Shape, b: Shape) -> float:
    """Minimum distance between two shapes (mm), 0 if they touch or overlap."""
    return float(a.distance_to(b))
