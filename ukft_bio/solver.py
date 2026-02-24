"""
ukft_bio.solver — Simulation Loop
===================================

Drives a biological system through N discrete choice steps, collecting observables.
Bio analogue of ukftphys.solver.  Core invariant: advance by choice index n, not time t.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Callable

from .physics import BioState, step_discrete


@dataclass
class TrajectoryRecord:
    """Lightweight record of one choice step for analysis / plotting."""
    n: int
    rho: float
    action: float
    fidelity: float
    genome_norm: float
    dt_bio: float


class BioSolver:
    """
    Run a BioState trajectory for `n_steps` discrete choice steps.

    Parameters
    ----------
    n_steps : int
        Total discrete choice steps to simulate.
    n_candidates : int
        Branching factor per choice step (more = more exploratory).
    mutation_scale : float
        Std-dev of genome perturbation per candidate generation.
    lambda_sem : float
        Semantic coherence penalty weight λ.
    rho_critical : float
        Catalytic half-saturation density.
    graduation_sigma : float
        If σ_rho < graduation_sigma → the agent has "graduated" (mastered the state).
        Mirroring Noogine Graduation concept (σ < 5.0).
    callback : optional callable(state, record)
        Called after every step for live reporting or custom logging.
    """

    def __init__(
        self,
        n_steps: int = 500,
        n_candidates: int = 8,
        mutation_scale: float = 0.1,
        lambda_sem: float = 0.1,
        rho_critical: float = 10.0,
        graduation_sigma: float = 0.5,
        callback: Optional[Callable[[BioState, TrajectoryRecord], None]] = None,
    ):
        self.n_steps = n_steps
        self.n_candidates = n_candidates
        self.mutation_scale = mutation_scale
        self.lambda_sem = lambda_sem
        self.rho_critical = rho_critical
        self.graduation_sigma = graduation_sigma
        self.callback = callback

    def run(self, initial_state: BioState) -> tuple[BioState, List[TrajectoryRecord]]:
        """
        Simulate from `initial_state` for `self.n_steps` choice steps.

        Returns
        -------
        (final_state, trajectory)
        """
        state = initial_state
        trajectory: List[TrajectoryRecord] = []
        rho_window: List[float] = []

        for n in range(self.n_steps):
            state = step_discrete(
                state,
                n_candidates=self.n_candidates,
                mutation_scale=self.mutation_scale,
                lambda_sem=self.lambda_sem,
                rho_critical=self.rho_critical,
            )

            rho_window.append(state.rho)
            if len(rho_window) > 50:
                rho_window.pop(0)

            sigma_rho = float(np.std(rho_window)) if len(rho_window) > 1 else 99.0

            record = TrajectoryRecord(
                n=n,
                rho=state.rho,
                action=0.0,  # recompute if needed by experiment
                fidelity=state.fidelity,
                genome_norm=float(np.linalg.norm(state.genome)),
                dt_bio=state.dt_bio,
            )
            trajectory.append(record)

            if self.callback is not None:
                self.callback(state, record)

            # Graduation check
            if sigma_rho < self.graduation_sigma and n > 50:
                print(
                    f"  [Graduation] Step n={n}: σ_ρ={sigma_rho:.4f} < "
                    f"{self.graduation_sigma} → agent mastered this level."
                )
                break

        return state, trajectory
