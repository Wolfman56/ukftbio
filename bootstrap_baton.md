# UKFT Biological Simulation: Bootstrap Baton

**⚡ START HERE (AI Agents Only)**

If you have been tasked to initialize the `ukftbio` codebase, this is your primary directive. Our goal is to replicate the successful "Choice-Guided" simulation architecture of `ukftphys` but applied to biological systems.

## Mission Objective

Create a functional, modular, and scientifically rigorous simulation environment for **UKFT Biology**. We believe biological complexity (folding, synthesis, evolution) is an emergent property of *collective choice* minimizing entropy in a shared knowledge field.

## Setup Instructions

### 1. Directory Structure
Initialize the following folder hierarchy (matching `ukftphys`):

```
ukftbio/
├── ukft_bio/                 # Core Python package
│   ├── __init__.py
│   ├── folding.py            # Protein folding logic (entropic choice)
│   ├── synthesis.py          # Molecular synthesis logic
│   └── evolution.py          # Evolutionary dynamics
├── experiments/              # Reproducible experiment scripts
│   ├── README.md
│   ├── 01_folding_basics.py  # Initial Proof of Concept
│   └── ...
├── data/                     # Data storage (git-ignored large files)
├── feedback/                 # Agent collaboration logs
│   └── README.md
├── references/               # Academic papers & theory
├── results/                  # Simulation outputs (git-ignored)
├── tools/                    # Utility scripts
├── environment.yaml          # Conda environment
├── requirements.txt          # Pip requirements
└── README.md                 # Project overview
```

### 2. Core Dependencies
Create `requirements.txt` and `environment.yaml` with essential libraries:
*   `numpy`, `scipy`, `pandas` (Numerics)
*   `biopython` (Structural biology tools)
*   `torch` or `jax` (If neural potentials are needed)
*   `plotly` (Visualization)
*   `networkx` (Graph theory for interaction networks)

### 3. Verification Protocol
1.  **Unit Tests**: Create a `tests/` folder. Ensure basic physics/chemistry sanity checks pass.
2.  **Experiment 01**: Create a simple script (`experiments/01_folding_basics.py`) that demonstrates a "peptide chain" folding into a lower-entropy state using a *toy* choice-minimization algorithm (not full MD, but UKFT-style discrete choices).

### 4. Feedback Integration
Establish the `feedback/` directory immediately. Agents should log their "Thought Process" and "Architecture Decisions" there, just like in `ukftphys`.

## Theoretical Alignment
Review `grok_x_bytedance_mol-syn.md` in the root. This is our "North Star". We are not just simulating biology; we are simulating **Biology as a Language of Choice**.

**Proceed with initialization.**
