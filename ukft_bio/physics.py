"""
ukft_bio.physics — Bio Action Functional & Discrete Stepper
=============================================================

The bio action functional S_bio is the biological analogue of the UKFT physical
action in ukftphys/ukft_sim/physics.py.  It replaces the Bohmian quantum potential
and entropic gravity with replication fidelity, catalytic context, and selection.

Formal definition (Grok, Feb 2026):
    S_bio[n] = K(n) + H_rep(fidelity, L) + C_cat(ρ) − ∇_sel · ρ_H + λ · H_sem

where:
    K         = choice kinetic cost  (roughness of the choice trajectory)
    H_rep     = replication entropy   ≈ L * h(fidelity)   [h = binary entropy]
    C_cat     = catalytic context gain (negative for high local ρ — well-catalyzed)
    ∇_sel·ρ_H = selection gradient dot heritable density (positive = favoured)
    H_sem     = semantic coherence penalty (non-functional sequence compositions)
    λ         = coherence weighting hyperparameter (default 0.1)

CRITICAL INVARIANTS:
    - Time is NEVER an independent variable; advance only via rho update.
    - dt_bio = dt_base / (rho_bio + eps)  — slower clock in high-density regions
    - all arrays indexed by choice n ∈ {0, 1, …, n_choices-1}
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ──────────────────────────────────────────────────────────────────────────────
# Core state class
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class BioState:
    """
    Minimal description of a biological agent in choice-space.

    Attributes
    ----------
    genome : np.ndarray of float32, shape (L,)
        Abstract genome as a real-valued sequence (each entry = choice amplitude).
        Concrete experiments may impose discrete alphabet (e.g., {0,1,2,3} for ACGT).
    rho : float
        Current knowledge/replication density.  rho = 0 → dead. rho >> 1 → thriving.
    fidelity : float in (0, 1]
        Per-base replication accuracy.  fidelity ≈ 1 − mutation_rate.
    choice_index : int
        Current discrete choice step n.
    environment : np.ndarray of float32, shape (D,)
        External context vector (selection gradient source).
    """
    genome: np.ndarray
    rho: float = 1.0
    fidelity: float = 0.99
    choice_index: int = 0
    environment: np.ndarray = field(default_factory=lambda: np.zeros(4))

    @property
    def genome_len(self) -> int:
        return len(self.genome)

    @property
    def dt_bio(self) -> float:
        """Emergent time step — slows down in high-density regions. dt ∝ 1/ρ."""
        return 1.0 / (self.rho + 1e-12)

    def clone(self, mutation_noise: float = 0.0) -> "BioState":
        """Replicate this state with optional mutation noise."""
        child_genome = self.genome.copy()
        if mutation_noise > 0.0:
            child_genome += np.random.randn(self.genome_len) * mutation_noise
        return BioState(
            genome=child_genome,
            rho=self.rho,
            fidelity=self.fidelity,
            choice_index=self.choice_index,
            environment=self.environment.copy(),
        )


# ──────────────────────────────────────────────────────────────────────────────
# Action functional components
# ──────────────────────────────────────────────────────────────────────────────

def choice_kinetic(genome: np.ndarray, genome_prev: Optional[np.ndarray]) -> float:
    """
    K(n) — Cost of the discrete choice transition.
    Measures roughness: large jumps in genome-space are energetically expensive.
    Analogue of the kinetic mismatch term in ukftphys step_discrete_action_minimizer.
    """
    if genome_prev is None:
        return 0.0
    delta = genome - genome_prev
    return 0.5 * float(np.dot(delta, delta))


def replication_entropy(fidelity: float, genome_len: int) -> float:
    """
    H_rep — Entropy cost of replication.
    Uses binary entropy h(p) = -p log p - (1-p) log(1-p) scaled by genome length.
    Perfect fidelity = 0 cost.  Near-50% fidelity = maximum cost.

    This is the UKFT-bio analogue of the quantum potential well — it constrains
    the system to avoid the entropic catastrophe (Eigen's error threshold).
    """
    p = np.clip(1.0 - fidelity, 1e-15, 1.0 - 1e-15)
    h = -p * np.log(p) - (1.0 - p) * np.log(1.0 - p)
    return float(genome_len * h)


def catalytic_context_gain(rho_local: float, rho_critical: float = 10.0) -> float:
    """
    C_cat — Catalytic context contribution (NEGATIVE = gain in dense environments).
    Logistic suppression: high local density provides catalytic scaffolding,
    reducing the effective action.  Saturates at -1 for rho >> rho_critical.
    """
    return -1.0 / (1.0 + np.exp(-(rho_local - rho_critical) / rho_critical))


def selection_gradient_term(
    genome: np.ndarray,
    environment: np.ndarray,
    rho_heritable: float,
) -> float:
    """
    ∇_sel · ρ_H — Selection pressure term.
    Projects the genome onto the environment vector, weighted by heritable density.
    Positive alignment → action decreases → trajectory is favoured.
    This is the discrete analogue of entropic gravity in choice-space.
    """
    env = environment / (np.linalg.norm(environment) + 1e-12)
    g = genome / (np.linalg.norm(genome) + 1e-12)
    alignment = float(np.dot(g[: len(env)], env))
    return -alignment * rho_heritable


def semantic_coherence_penalty(genome: np.ndarray, lambda_sem: float = 0.1) -> float:
    """
    λ · H_sem — Penalty for non-coherent (low-information) sequence composition.
    Approximated as entropy of the empirical symbol distribution of the genome.
    Low-complexity / repetitive genomes pay a higher semantic penalty because
    they carry less heritable knowledge.
    """
    # Discretize genome values into 16 bins to estimate symbol entropy
    bins = np.linspace(genome.min() - 1e-9, genome.max() + 1e-9, 17)
    counts, _ = np.histogram(genome, bins=bins)
    probs = counts / (counts.sum() + 1e-12)
    probs = probs[probs > 0]
    shannon = -np.sum(probs * np.log(probs))
    max_entropy = np.log(16)
    # Penalty is HIGH when symbol entropy is LOW (repetitive = incoherent)
    return lambda_sem * (max_entropy - shannon)


def god_attractor_potential(rho: float, rho_god: float = 1e6) -> float:
    """
    V_god — God Attractor pull toward the ultra-high-ρ basin.
    Logarithmic well: trajectories are always pulled upward in rho.
    At rho ≪ rho_god this acts as a weak constant tug.
    This corresponds to the Theosphere memory tier in ukftphys Exp 19/20.
    """
    return -np.log(rho + 1.0) + np.log(rho_god)


def bio_action(
    state: BioState,
    genome_candidate: np.ndarray,
    genome_prev: Optional[np.ndarray] = None,
    lambda_sem: float = 0.1,
    rho_critical: float = 10.0,
    rho_god: float = 1e6,
) -> float:
    """
    Full bio-action functional S_bio evaluated for a candidate genome transition.

    S_bio = K + H_rep + C_cat − ∇_sel · ρ_H + λ H_sem + V_god

    Returns
    -------
    float : action value.  Lower = more favoured by choice operator.
    """
    K       = choice_kinetic(genome_candidate, genome_prev)
    H_rep   = replication_entropy(state.fidelity, len(genome_candidate))
    C_cat   = catalytic_context_gain(state.rho, rho_critical)
    S_sel   = selection_gradient_term(genome_candidate, state.environment, state.rho)
    H_sem   = semantic_coherence_penalty(genome_candidate, lambda_sem)
    V_god   = god_attractor_potential(state.rho, rho_god)
    return K + H_rep + C_cat + S_sel + H_sem + V_god


# ──────────────────────────────────────────────────────────────────────────────
# Discrete choice stepper
# ──────────────────────────────────────────────────────────────────────────────

def step_discrete(
    state: BioState,
    n_candidates: int = 8,
    mutation_scale: float = 0.1,
    lambda_sem: float = 0.1,
    rho_update_rate: float = 0.05,
    rho_critical: float = 10.0,
) -> BioState:
    """
    Advance the biological state by ONE discrete choice step n → n+1.

    Procedure:
    1. Generate n_candidates perturbed genome variants (exploratory choices).
    2. Evaluate S_bio for each candidate.
    3. Select the argmin candidate (Principle of Least Action).
    4. Update ρ_bio based on alignment gain with the selected genome.
    5. Increment choice index.

    This is the bio equivalent of ukftphys.step_discrete_action_minimizer.
    """
    genome_prev = state.genome.copy()

    # Generate candidates: current genome ± random perturbations
    candidates = [
        state.genome + np.random.randn(state.genome_len) * mutation_scale
        for _ in range(n_candidates)
    ]
    # Always include the unperturbed genome (null choice)
    candidates.append(state.genome.copy())

    # Evaluate action for each candidate
    actions = [
        bio_action(
            state,
            cand,
            genome_prev=genome_prev,
            lambda_sem=lambda_sem,
            rho_critical=rho_critical,
        )
        for cand in candidates
    ]

    # Select minimum-action candidate (discrete least-action principle)
    best_idx = int(np.argmin(actions))
    best_genome = candidates[best_idx]

    # Update ρ_bio: reward for alignment improvement with environment
    env = state.environment / (np.linalg.norm(state.environment) + 1e-12)
    g_new = best_genome / (np.linalg.norm(best_genome) + 1e-12)
    g_old = genome_prev / (np.linalg.norm(genome_prev) + 1e-12)
    alignment_gain = float(np.dot(g_new[: len(env)], env)) - float(np.dot(g_old[: len(env)], env))
    new_rho = max(1e-6, state.rho * (1.0 + rho_update_rate * alignment_gain))

    return BioState(
        genome=best_genome,
        rho=new_rho,
        fidelity=state.fidelity,
        choice_index=state.choice_index + 1,
        environment=state.environment.copy(),
    )
