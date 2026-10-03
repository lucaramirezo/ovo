# Blade design: constraint-aware BEM

Status: DRAFT, results being regenerated (`python -m hydroegg.bem`).

Source: `hydroegg/bem.py`. Outputs: `data/blade_design.json`, `docs/figures/bem_*.png`.

## Method

Blade Element Momentum on 24 annular elements (cosine-clustered toward root and tip)
from the aerodynamic root at r = 30 mm to the tip at r = 97 mm (perpendicular radius).

- Prandtl tip loss and Prandtl hub loss (hub-loss radius = aerodynamic root).
- High induction: Buhl (2005) empirical branch above a = 0.4, in the closed form of
  Ning (2014), continuous with momentum theory at a = 0.4.
- Coning: the pitch axis leans aft by the cone angle beta. Momentum is taken through
  the projected annulus 2 pi r dr; the blade sees U (1 - a) cos(beta) normal to its span;
  the spanwise inflow component is neglected (independence principle).
- Per-element Reynolds number from the local relative speed W and chord.
- Polars: NeuralFoil (xlarge), an XFOIL-trained model, evaluated on the exact section
  that the CAD builds (NACA 44xx plus the blunt trailing edge), n_crit 5, free
  transition. Sweeps and optimisation use a per-station table (linear in alpha and
  TE thickness, cubic in log Re, log Cd); every reported number is re-solved with
  direct NeuralFoil calls.
- Solver: relaxed fixed point on (a, a'), iterated per element until converged; the
  elements that do not converge fall back to Ning's bracketed residual in phi.
  Elements that still fail are flagged and counted in every result.
