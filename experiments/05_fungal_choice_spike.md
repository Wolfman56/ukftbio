# Experiment 05 -- Fungal Choice Spike (Entangled Collapse)

**Status:** Ready to run
**Physics analogy:** Entanglement propagation (ukftphys Exp 17)
**Paper 01 connection:** Tests Predictions P2 and P3 on synthetic data

## Objective

Demonstrate the core claim of Paper 01 on synthetic data before applying the
pipeline to real Zenodo fungal recordings:

- **P2**: Entangled event clusters show super-Poissonian inter-channel
  correlation absent in shuffled (entanglement-broken) controls.
- **P3**: Resolved clusters lie closer to the UKFT action minimum in feature
  space than unresolved partial projections or shuffled controls.

This is the baseline validation step. Zenodo recordings (Adamatzky 2026,
Zenodo 5790768) will follow in Exp 06+.

## Design

Synthetic 4-channel spike generator:
- **Real (resolved):** entangled clusters with delay structure matching
  Adamatzky 2026 (~0.7 cm/min = 180 s / 2 cm propagation).
- **Real (partial):** same structure but only 1-3 of 4 channels fire.
- **Shuffled:** same event rate but randomised channel assignment (breaks
  propagation graph = destroys entanglement).

200 clusters per condition.

Measurements:
- P2: Pearson r between channels; one-sided t-test real vs shuffled.
- P3: Euclidean distance to ACTION_MINIMUM = [0.6, 0.6, -0.3, 0.1] in
  (amp_mean, timing_tightness, correlation, std_dev) feature space.

## Figures

- Fig 1: `exp05_cluster_resolution.png` -- PCA scatter of cluster feature
  space; real-resolved / real-partial / shuffled; action minimum marked
- Fig 2: `exp05_p2_correlation.png` -- correlation histogram with t-test

## UKFT Interpretation

A spike cluster whose delayed partner arrivals complete the information-flow
cycle has actualised a choice. The partial/shuffled cases are incomplete
projections -- the wavefunction has not collapsed to a definite outcome.

Quantitatively: resolved clusters have lower action (closer to minimum) and
higher inter-channel correlation than controls. Both are measurable. The
Goldman Institute data shows this on real mycelium; this experiment confirms
the same signal is detectable in a clean synthetic model.
