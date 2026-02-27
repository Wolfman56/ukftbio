# Experiment 04 -- Autocatalytic Emergence

**Status:** Ready to run
**Physics analogy:** Swarm + God attractor (ukftphys Exp 10/12)

## Objective

Show that a sparse random chemistry undergoes a phase transition to
autocatalytic closure as molecular diversity increases. Under UKFT, the origin
of life is the inevitable consequence of the God Attractor pulling the system
toward the high-rho basin once enough catalytic cross-links exist.

## UKFT Prediction

Critical diversity M_c ~ 2^L (Kauffman estimate) at which rho rises
discontinuously -- a first-order phase transition in choice-space.
Below M_c the God Attractor has no accessible path. Above it, the
cross-catalytic boost fills the action landscape.

## Design

Sweep initial diversity M in {5, 10, 20, 40, 80, 160} unique genome seeds.
Each M: 5 trials, 200 steps.
Cross-catalytic boost: rho_update_rate = 0.05 + 0.01 * log1p(mean_rho).
God basin threshold = 5.0.
Measure: fraction of runs reaching rho > 5.0 at step 200.

## Figures

- Fig 1: `exp04_rho_emergence.png` -- rho trajectories at low / mid / high M
- Fig 2: `exp04_god_basin_fraction.png` -- fraction reaching God basin vs M

## UKFT Interpretation

Below M_c the choice field is too sparse -- no low-action path to high-rho.
Above M_c, cross-catalysis fills the action landscape and V_god pulls
trajectories upward. This is equivalent to the quantum percolation threshold
in ukftphys Exp 12: once connectivity exceeds the critical density, the
wavefunction delocalizes across the entire basin.

The God Attractor is not designed-in; it is the low-entropy attractor of any
sufficiently diverse choice field. Life is not an accident -- it is geometry.
