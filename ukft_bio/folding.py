"""
ukft_bio.folding — Protein Folding as Choice-Guided Least-Action
================================================================

UKFT resolution of the Levinthal Paradox:
  Protein folding is NOT random search.  It is guided navigation through
  choice-space attracted to the minimum-action (native) conformation.

The native fold is an entropic attractor — the high-ρ basin corresponding to
maximum catalytic context gain (the protein performs its function = high ρ).

Folding state:
  - Backbone dihedral angles (φ, ψ) per residue → the "genome" in angle-space
  - ρ_fold = inverse RMSD to native (or proxy: contact order score)
  - Each discrete step selects the dihedral perturbation minimising S_bio

Future work (Gemini experimentalist):
  - Exp 06: "Molecular Prophet" — predict novel fold from sequence alone
  - Exp 07: Entropic regularisation vs Brownian MD comparison
"""

from __future__ import annotations

import numpy as np
from typing import Optional, List, Tuple

from .physics import BioState, step_discrete, bio_action


class ResidueChain:
    """
    Simplified backbone representation.
    Each residue i has (phi_i, psi_i) dihedral angles in radians.
    Genome = flattened [phi_0, psi_0, phi_1, psi_1, ...].
    """

    def __init__(self, n_residues: int, sequence: Optional[np.ndarray] = None):
        self.n_residues = n_residues
        self.sequence = (
            sequence
            if sequence is not None
            else np.random.randint(0, 20, size=n_residues)  # 20 amino acids
        )
        # Random initialisation (Levinthal starting point)
        self.dihedrals = np.random.uniform(-np.pi, np.pi, size=2 * n_residues)

    @property
    def phi(self) -> np.ndarray:
        return self.dihedrals[0::2]

    @property
    def psi(self) -> np.ndarray:
        return self.dihedrals[1::2]

    def contact_order_score(self) -> float:
        """
        Proxy for fold quality: reward for long-range backbone contacts.
        Rough measure: variance of cos(phi+psi) — well-packed folds have
        low variance (regular secondary structure).
        Higher score = better fold (higher ρ_fold).
        """
        angles = self.phi + self.psi
        return 1.0 / (np.var(np.cos(angles)) + 0.1)

    def to_bio_state(self, environment: Optional[np.ndarray] = None) -> BioState:
        """Convert chain to BioState with fold quality as ρ."""
        rho = self.contact_order_score()
        env = environment if environment is not None else np.zeros(4)
        return BioState(
            genome=self.dihedrals.copy(),
            rho=rho,
            fidelity=0.999,  # folding is deterministic (not mutating)
            choice_index=0,
            environment=env,
        )


class ProteinFolder:
    """
    Folds a ResidueChain by UKFT discrete choice minimisation.

    Each choice step perturbs a random subset of dihedral angles and
    selects the perturbation that minimises S_bio (dominated here by the
    catalytic context gain term, which rewards high ρ_fold = better packed).

    Parameters
    ----------
    n_steps : int
        Maximum folding choices (analogous to Monte Carlo steps).
    perturb_fraction : float
        Fraction of dihedrals to perturb per candidate (0 < f ≤ 1).
    perturb_scale : float
        Std-dev of dihedral perturbation in radians.
    """

    def __init__(
        self,
        n_steps: int = 1000,
        perturb_fraction: float = 0.2,
        perturb_scale: float = 0.3,
        n_candidates: int = 12,
    ):
        self.n_steps = n_steps
        self.perturb_fraction = perturb_fraction
        self.perturb_scale = perturb_scale
        self.n_candidates = n_candidates

    def _perturb(self, dihedrals: np.ndarray) -> np.ndarray:
        d = dihedrals.copy()
        n_perturb = max(1, int(self.perturb_fraction * len(d)))
        idx = np.random.choice(len(d), size=n_perturb, replace=False)
        d[idx] += np.random.randn(n_perturb) * self.perturb_scale
        # Wrap to [-π, π]
        d = ((d + np.pi) % (2 * np.pi)) - np.pi
        return d

    def fold(self, chain: ResidueChain, verbose: bool = True) -> Tuple[ResidueChain, List[float]]:
        """
        Run the folding trajectory.

        Returns
        -------
        (folded_chain, rho_history)
        """
        state = chain.to_bio_state()
        rho_history = [state.rho]

        for step in range(self.n_steps):
            # Generate candidates
            candidates = [self._perturb(state.genome) for _ in range(self.n_candidates)]
            candidates.append(state.genome.copy())

            # Score by contact order (ρ_fold proxy) — no full physics needed here
            scores = []
            for cand in candidates:
                tmp_chain = ResidueChain(chain.n_residues, chain.sequence)
                tmp_chain.dihedrals = cand
                scores.append(-tmp_chain.contact_order_score())  # minimise → maximise score

            best_idx = int(np.argmin(scores))
            best_dihedrals = candidates[best_idx]

            tmp_chain = ResidueChain(chain.n_residues, chain.sequence)
            tmp_chain.dihedrals = best_dihedrals
            new_rho = tmp_chain.contact_order_score()

            state = BioState(
                genome=best_dihedrals,
                rho=new_rho,
                fidelity=state.fidelity,
                choice_index=step + 1,
                environment=state.environment,
            )
            rho_history.append(new_rho)

            if verbose and step % 100 == 0:
                print(f"  [Folding] step={step:4d}  ρ_fold={new_rho:.4f}")

        # Write result back to chain
        chain.dihedrals = state.genome
        return chain, rho_history
