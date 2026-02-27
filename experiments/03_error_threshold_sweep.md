# Experiment 03 -- Error Threshold Sweep

**Status:** Ready to run
**Physics analogy:** Entropic alpha sweep (ukftphys Exp 03)

## Objective

Map the phase boundary between replication survival and error catastrophe as a
function of fidelity. Show that the sharp transition is a fundamental property
of S_bio, not a model artifact.

## UKFT Prediction

Critical fidelity f_c ~ 1 - 1/(L * sigma_selection).
The transition is sharp (first-order-like): rho_final vs fidelity produces a
step function, not a smooth decay.

## Design

Sweep fidelity in [0.70, 0.995] at 22 evenly-spaced values.
Each level: 5 independent populations of 50 agents for 200 steps.
Measure: mean final rho_bio across surviving agents.

## Figures

- Fig 1: `exp03_phase_diagram.png` -- mean final rho vs fidelity, f_c marked
- Fig 2: `exp03_h_rep_landscape.png` -- H_rep(fidelity, L) nonlinear entropy cost

## UKFT Interpretation

The phase diagram is the bio analogue of the entropic alpha sweep in ukftphys
Exp 03. Below f_c, the replication entropy H_rep overwhelms the selection
gradient and the population loses all heritable structure.

The sharpness of the transition is a UKFT fingerprint: entropy gaps do not
degrade gracefully -- they collapse. This is the analogue of the quantum
decoherence threshold in the physical UKFT.
