"""
ukft_bio.evolution — Evolutionary Dynamics as Noospheric Choice
================================================================

UKFT view of evolution (from grok_x_bytdance_noobio_1_chat.md):
  Natural selection = least-action pruning in phylogenetic choice-space.
  Mutation = exploratory branching of choice trajectories.
  Fitness = local ρ_bio (knowledge/replication density).
  Speciation = choice-space trajectories diverging toward distinct basins.

The God Attractor ensures evolution is NOT purely blind random walk:
  every lineage is weakly pulled toward richer, more coherent genomes
  (higher ρ_bio = more heritable knowledge = god-basin direction).

Key result to demonstrate (Gemini Exp 01–05):
  Error threshold (Eigen 1971) emerges naturally from S_bio minimisation.
  When fidelity drops below critical value, replication_entropy term
  dominates → system cannot maintain ρ > threshold → quasi-species collapses.

References:
  - Eigen (1971): quasi-species equation
  - Kauffman (1993): autocatalytic sets  ← map to choice-graph nodes
  - England (2013): dissipation-driven adaptation ← rho_update_rate
  - Levin (diverse intelligence): cognition at all scales ← multi-agent CLKOS
"""

from __future__ import annotations

import numpy as np
from typing import List, Optional, Tuple
from dataclasses import dataclass

from .physics import BioState, step_discrete, replication_entropy


@dataclass
class Replicator:
    """
    A single self-replicating entity in choice-space.
    Wraps BioState with a lineage ID and generation counter.
    """
    state: BioState
    lineage_id: int
    generation: int = 0
    alive: bool = True

    @property
    def fitness(self) -> float:
        """ρ_bio IS fitness in UKFT-bio."""
        return self.state.rho

    def replicate(self, mutation_scale: float = 0.05) -> "Replicator":
        """
        Produce an offspring via least-action replication.
        Mutation noise controlled by fidelity: noise ∝ (1 - fidelity).
        """
        effective_noise = mutation_scale * (1.0 - self.state.fidelity)
        child_genome = self.state.genome.copy()
        child_genome += np.random.randn(len(child_genome)) * effective_noise
        child_state = BioState(
            genome=child_genome,
            rho=self.state.rho,
            fidelity=self.state.fidelity,
            choice_index=0,
            environment=self.state.environment.copy(),
        )
        return Replicator(state=child_state, lineage_id=self.lineage_id,
                          generation=self.generation + 1)


class EvolutionaryDynamics:
    """
    Population-level evolutionary dynamics under UKFT.

    Each generation:
    1. All replicators take one discrete choice step (mutation via step_discrete).
    2. Fitness = ρ_bio.  Death if ρ < rho_death.
    3. Reproduce proportionally to fitness (Wright-Fisher-like, but choice-driven).
    4. Track lineage trees and rho distributions.

    Experiment target (Exp 01 — Minimal Entropic Replicator):
    - Show error threshold: below fidelity_critical, population collapses.
    - Above threshold: quasi-species maintains ρ > rho_min indefinitely.
    - God attractor visible as slow drift of population mean ρ upward over eons.

    Parameters
    ----------
    population_size : int
    n_generations : int
    genome_len : int
    fidelity : float         Per-base replication accuracy.
    rho_death : float        ρ below this → replicator dies.
    rho_god : float          God attractor basin centre.
    """

    def __init__(
        self,
        population_size: int = 50,
        n_generations: int = 200,
        genome_len: int = 20,
        fidelity: float = 0.99,
        rho_death: float = 0.5,
        rho_god: float = 1e4,
        mutation_scale: float = 0.05,
        n_candidates: int = 6,
    ):
        self.population_size = population_size
        self.n_generations = n_generations
        self.genome_len = genome_len
        self.fidelity = fidelity
        self.rho_death = rho_death
        self.rho_god = rho_god
        self.mutation_scale = mutation_scale
        self.n_candidates = n_candidates

    def _init_population(
        self,
        environment: Optional[np.ndarray] = None,
    ) -> List[Replicator]:
        env = environment if environment is not None else np.random.randn(4)
        pop = []
        for i in range(self.population_size):
            genome = np.random.randn(self.genome_len) * 0.5
            state = BioState(
                genome=genome,
                rho=np.random.uniform(1.0, 5.0),
                fidelity=self.fidelity,
                choice_index=0,
                environment=env.copy(),
            )
            pop.append(Replicator(state=state, lineage_id=i))
        return pop

    def run(
        self,
        environment: Optional[np.ndarray] = None,
        verbose: bool = True,
    ) -> Tuple[List[dict], List[List[float]]]:
        """
        Run the full evolutionary simulation.

        Returns
        -------
        (generation_stats, lineage_rho_histories)
          generation_stats: list of dicts with mean_rho, max_rho, pop_size per generation
          lineage_rho_histories: one list of rho values per founding lineage
        """
        population = self._init_population(environment)
        n_lineages = self.population_size
        lineage_histories: dict[int, List[float]] = {i: [] for i in range(n_lineages)}
        gen_stats: List[dict] = []

        for gen in range(self.n_generations):
            # ── Step 1: Each individual takes one choice step ──────────────
            for rep in population:
                if rep.alive:
                    rep.state = step_discrete(
                        rep.state,
                        n_candidates=self.n_candidates,
                        mutation_scale=self.mutation_scale * (1.0 - rep.state.fidelity + 0.01),
                    )
                    lineage_histories[rep.lineage_id % n_lineages].append(rep.fitness)

            # ── Step 2: Death (ρ below threshold) ─────────────────────────
            population = [r for r in population if r.fitness > self.rho_death]

            if not population:
                print(f"  [Evolution] Population extinct at generation {gen}!")
                break

            # ── Step 3: Reproduce to maintain population size ──────────────
            fitnesses = np.array([r.fitness for r in population])
            probs = fitnesses / fitnesses.sum()
            n_offspring = self.population_size - len(population)
            if n_offspring > 0:
                parents = np.random.choice(len(population), size=n_offspring, p=probs)
                for p_idx in parents:
                    offspring = population[p_idx].replicate(mutation_scale=self.mutation_scale)
                    population.append(offspring)

            # ── Step 4: Record stats ───────────────────────────────────────
            rhos = [r.fitness for r in population]
            stats = {
                "generation": gen,
                "mean_rho": float(np.mean(rhos)),
                "max_rho": float(np.max(rhos)),
                "min_rho": float(np.min(rhos)),
                "pop_size": len(population),
                "replication_entropy": replication_entropy(self.fidelity, self.genome_len),
            }
            gen_stats.append(stats)

            if verbose and gen % 20 == 0:
                print(
                    f"  [Gen {gen:4d}] pop={len(population):3d}  "
                    f"mean_ρ={stats['mean_rho']:.3f}  max_ρ={stats['max_rho']:.3f}"
                )

        return gen_stats, list(lineage_histories.values())
