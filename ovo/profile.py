"""Hull profile: Myring body of revolution tuned to read as an egg (L/D ~ 2.2).

Myring, P.B. (1976) "A theoretical study of body drag in subcritical axisymmetric
flow". Nose is a super-ellipse, then a short parallel section, then a cubic tail
that closes to a point with half-angle theta. The static nacelle is the hull from
x = 0 to STATIC_TAIL_X; aft of that the same curve is continued by the rotating
hub and spinner, so stowed the whole thing reads as one egg.

Frame: X is the body axis, flow travels +X, nose at x = 0.
"""

from __future__ import annotations

import math

import numpy as np

from . import params as P


def radius(x: float | np.ndarray) -> float | np.ndarray:
    """Outer hull radius r(x) in mm. Zero outside [0, HULL_L]."""
    a, b, c, d, n = P.HULL_NOSE_A, P.HULL_MID_B, P.HULL_TAIL_C, P.HULL_D, P.HULL_NOSE_N
    tan_t = math.tan(math.radians(P.HULL_TAIL_THETA))
    x = np.asarray(x, dtype=float)
    r = np.zeros_like(x)

    nose = (x >= 0) & (x < a)
    r[nose] = (d / 2) * np.power(np.clip(1 - ((x[nose] - a) / a) ** 2, 0, None), 1 / n)

    mid = (x >= a) & (x <= a + b)
    r[mid] = d / 2

    tail = (x > a + b) & (x <= a + b + c)
    t = x[tail] - a - b
    r[tail] = (
        d / 2
        - (3 * d / (2 * c**2) - tan_t / c) * t**2
        + (d / c**3 - tan_t / c**2) * t**3
    )
    r = np.clip(r, 0, None)
    return float(r) if r.ndim == 0 else r


def slope(x: float, h: float = 1e-3) -> float:
    """dr/dx by central difference."""
    return float((radius(x + h) - radius(x - h)) / (2 * h))


def sample(x0: float = 0.0, x1: float | None = None, n: int = 400) -> np.ndarray:
    """(n, 2) array of (x, r) points along the outer profile, cosine-spaced so the
    blunt nose and the tail tip get more points."""
    x1 = P.HULL_L if x1 is None else x1
    s = (1 - np.cos(np.linspace(0, math.pi, n))) / 2
    xs = x0 + (x1 - x0) * s
    return np.column_stack([xs, radius(xs)])


def inner_radius(x: float | np.ndarray, wall: float | None = None) -> float | np.ndarray:
    """Approximate inner surface radius of a constant-thickness shell (wall measured
    along the surface normal). Good to ~0.1 mm where the slope is moderate."""
    wall = P.HULL_WALL if wall is None else wall
    x = np.asarray(x, dtype=float)
    h = 1e-3
    drdx = (radius(x + h) - radius(x - h)) / (2 * h)
    r = radius(x) - wall * np.sqrt(1 + drdx**2)
    r = np.clip(r, 0, None)
    return float(r) if r.ndim == 0 else r


if __name__ == "__main__":
    for xi in (0, 10, 25, 40, 55, 60, 75, 90, 100, 110, 118, 124, 130, 145, 160, 165):
        print(f"x={xi:6.1f}  r={radius(xi):6.2f}  r_in={inner_radius(xi):6.2f}  dr/dx={slope(xi):+.3f}")
