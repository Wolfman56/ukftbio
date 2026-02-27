# Experiments — UKFT Biology

Numbered experiments following the ukftphys convention.
Each experiment has a `.py` script and a `.md` companion explainer.

> **Design philosophy (mirrors ukftphys):**  
> Establish the canonical **baseline set (Exp 01–05)** first — confirming the  
> framework works at the most fundamental level before entering the  
> paper / hypothesis / proof cycle.  
> Exp 02 (Levinthal) is the bio double-slit.

---

## Baseline Set — confirm before advancing

| Exp | Title | Physics Analogy | Key UKFT Concept | Status |
|-----|-------|----------------|-----------------|--------|
| 01 | Minimal Entropic Replicator | Free particle | Error threshold as S_bio floor | ✅ Ready |
| 02 | Levinthal Paradox | **Double slit** | Folding = choice collapse, not random search | ✅ Ready |
| 03 | Error Threshold Sweep | Entropic alpha sweep | fidelity ↔ ρ phase boundary | ✅ Ready |
| 04 | Autocatalytic Emergence | Swarm / God attractor | Life as rho phase transition | Ready |
| 05 | Fungal Choice Spike | Entanglement propagation | Entangled cluster = EPR-bounded collapse | Ready |

## Extended Series

| Exp | Title | Status |
|-----|-------|--------|
| 06 | Entropic Genetic Code | Next |
| 07 | Morphogen-Guided Development | Future |
| 08 | Dissipation-Driven Adaptation | Future |
| 09 | Molecular Prophet | Future |
| 10 | Collective Choice (CLKOS) | Future |

## Running the baseline

    conda activate ukftbio
    python experiments/01_minimal_entropic_replicator.py
    python experiments/02_levinthal_folding_paradox.py
    python experiments/03_error_threshold_sweep.py
    python experiments/04_autocatalytic_emergence.py
    python experiments/05_fungal_choice_spike.py

---

## Canonical Figures

| Exp | Primary Figure | What it proves |
|-----|---------------|----------------|
| 01 | `exp01_rho_trajectory.png` | ρ survives (f=0.99) vs collapses (f=0.80) — error threshold |
| 02 | `exp02_folding_landscape.png` | UKFT collapses to native basin; random walk lost in Levinthal space |
| 03 | `exp03_phase_diagram.png` | Sharp f_critical boundary in (fidelity, ρ_final) space |
| 04 | `exp04_rho_emergence.png` | ρ jumps discontinuously at autocatalytic threshold → God basin |
| 05 | `exp05_cluster_resolution.png` | Real clusters resolve at action minimum; shuffled controls do not |
