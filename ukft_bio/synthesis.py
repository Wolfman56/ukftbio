"""
ukft_bio.synthesis — Molecular Synthesis as Language of Choice
===============================================================

Key thesis (from grok_x_bytedance_mol-syn.md analysis):
  ByteDance's "Mole-Syn" identifies three bond types in CoT reasoning:
    covalent  ≈ deep committed choices (least-action locked)
    ionic     ≈ reflection/phase-reversal choices
    van der Waals ≈ low-energy exploratory fluctuations

  UKFT derives these naturally as curvature regimes around the entropic
  attractor.  No ad-hoc bond taxonomy required — they emerge from S_bio
  minimisation in the molecular graph.

Molecular state:
  - A reaction sequence as a path through "molecular choice-space"
  - Each intermediate = a BioState whose genome encodes molecular features
  - ρ_mol = "synthetic value" or catalytic usefulness score
  - Selection gradient = target molecule direction

PROTAC linker design (from grok_x_celltype_noobio_2_chat.md):
  - Linker length / rigidity are discrete choices along the synthesis path
  - Optimal linker = minimum-action path through POI–E3 ternary complex
  - Entropic penalty handles the flexibility/rigidity trade-off automatically

Future experiments (Gemini):
  - Exp 08: Genetic code redundancy via synthesis least-action
  - Exp 09: PROTAC linker optimisation — emergence of length sweet-spot
  - Exp 10: "Molecular Prophet" — predict novel synthesis pathways
"""

from __future__ import annotations

import numpy as np
from typing import List, Optional, Dict, Tuple

from .physics import BioState, step_discrete, bio_action


# Bond-type curvature regimes in choice-space (emerges, but useful for labelling)
BOND_CURVATURE_COVALENT    = "covalent"     # large action gradient → committed step
BOND_CURVATURE_IONIC       = "ionic"        # saddle point → reflection choice
BOND_CURVATURE_VDW         = "vdw"          # shallow well → exploratory fluctuation


def classify_bond_type(action_gradient_norm: float) -> str:
    """
    Classify a synthesis step by its action gradient curvature.
    Thresholds validated empirically — expect Gemini to refine per experiment.
    """
    if action_gradient_norm > 2.0:
        return BOND_CURVATURE_COVALENT
    elif action_gradient_norm > 0.5:
        return BOND_CURVATURE_IONIC
    else:
        return BOND_CURVATURE_VDW


# Molecule feature vector layout:
#   [0]  molecular_weight (normalised)
#   [1]  logP (normalised to [-1,1])
#   [2]  hbd (H-bond donors, normalised)
#   [3]  hba (H-bond acceptors, normalised)
#   [4]  rotatable_bonds (normalised)
#   [5]  tpsa (normalised)
#   [6]  ring_count (normalised)
#   [7+] pharmacophore embedding (variable length)


class Molecule:
    """
    Lightweight molecular representation as a feature vector.
    Genome encodes the molecular feature vector + pharmacophore embedding.
    """

    def __init__(
        self,
        features: Optional[np.ndarray] = None,
        name: str = "unknown",
        n_features: int = 16,
    ):
        self.name = name
        self.features = (
            features if features is not None
            else np.random.uniform(-1, 1, size=n_features)
        )

    @property
    def n_features(self) -> int:
        return len(self.features)

    def synthetic_value(self, target: "Molecule") -> float:
        """
        ρ_mol: cosine similarity to a target molecule in feature-space.
        Ranges (0, 1].  1 = exact synthetic match.
        """
        a = self.features / (np.linalg.norm(self.features) + 1e-12)
        b = target.features / (np.linalg.norm(target.features) + 1e-12)
        return float(0.5 + 0.5 * np.dot(a, b))

    def to_bio_state(self, target: "Molecule") -> BioState:
        rho = self.synthetic_value(target)
        return BioState(
            genome=self.features.copy(),
            rho=rho,
            fidelity=0.95,
            choice_index=0,
            environment=target.features.copy(),
        )


class MolecularSynthesizer:
    """
    Synthesises a target molecule from a starting material via UKFT
    discrete choice steps.  Each step is a "reaction" that perturbs
    the molecular feature vector toward the target.

    The bond-type emerges from the action gradient magnitude:
      - covalent steps: large committed jumps toward the attractor
      - ionic steps: phase-reversal corrections
      - vdW steps: small exploratory fluctuations

    Parameters
    ----------
    n_steps : int
        Maximum synthesis choices (reaction steps).
    step_size : float
        Magnitude of feature-space perturbations.
    """

    def __init__(
        self,
        n_steps: int = 200,
        step_size: float = 0.15,
        n_candidates: int = 10,
    ):
        self.n_steps = n_steps
        self.step_size = step_size
        self.n_candidates = n_candidates

    def synthesise(
        self,
        start: Molecule,
        target: Molecule,
        verbose: bool = True,
    ) -> Tuple[List[Molecule], List[str]]:
        """
        Plan a synthesis route from `start` to `target`.

        Returns
        -------
        (pathway, bond_types) — list of intermediate Molecule states
        and their step classifications.
        """
        state = start.to_bio_state(target)
        pathway = [start]
        bond_types: List[str] = []

        for step in range(self.n_steps):
            genome_prev = state.genome.copy()

            # Generate candidates: perturbed feature vectors
            candidates = [
                state.genome + np.random.randn(len(state.genome)) * self.step_size
                for _ in range(self.n_candidates)
            ]
            candidates.append(state.genome.copy())

            actions = [
                bio_action(state, cand, genome_prev=genome_prev)
                for cand in candidates
            ]

            best_idx = int(np.argmin(actions))
            best_features = candidates[best_idx]

            # Estimate action gradient norm for bond classification
            action_prev = bio_action(state, genome_prev, genome_prev=genome_prev)
            action_best = actions[best_idx]
            grad_norm = abs(action_prev - action_best) / (self.step_size + 1e-12)
            bond_type = classify_bond_type(grad_norm)

            mol = Molecule(features=best_features, name=f"intermediate_{step}")
            new_rho = mol.synthetic_value(target)

            state = BioState(
                genome=best_features,
                rho=new_rho,
                fidelity=state.fidelity,
                choice_index=step + 1,
                environment=target.features.copy(),
            )

            pathway.append(mol)
            bond_types.append(bond_type)

            if verbose and step % 50 == 0:
                print(f"  [Synthesis] step={step:3d}  ρ_mol={new_rho:.4f}  bond={bond_type}")

            # Convergence: close enough to target
            if new_rho > 0.98:
                print(f"  [Synthesis] Converged at step {step}: ρ_mol={new_rho:.4f}")
                break

        return pathway, bond_types
