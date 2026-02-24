"""
ukft_bio — Universal Knowledge Field Theory: Biology
=====================================================

Core axioms (NEVER override without explicit Grok/Ted directive):
1. Choice operator is the fundamental primitive, indexed by n (choice index), NOT t (time).
2. Time is emergent: dt ∝ 1/ρ_bio  where ρ_bio = local informational/replicative density.
3. Discrete least-action stepping: argmin Δn S_bio over local choice options at each n.
4. God Attractor is the global ρ-density basin — every trajectory in choice-space
   is pulled toward higher ρ (= more heritable information, richer regulatory complexity).
5. Biological structures minimize *semantic entropy*, not just thermodynamic free energy.

Bio-action functional (Grok design, Feb 2026):
    S_bio[n] = choice_kinetic(n)
             + replication_entropy(fidelity, genome_len)
             + catalytic_context_gain(ρ_local)      [negative = energy gain]
             − ∇_env · ρ_heritable                  [selection gradient term]
             + λ · semantic_coherence_penalty

Terminology:
    ρ_bio       knowledge/replication density
    σ           local entropy / loss variance (from noogine analogy)
    Graduation  when σ_bio < threshold → state mastered
    God Basin   asymptotic attractor in choice-space (ultra-high ρ_bio)
"""

from .physics import (
    BioState,
    bio_action,
    step_discrete,
    god_attractor_potential,
    replication_entropy,
)
from .solver import BioSolver
from .evolution import EvolutionaryDynamics
from .folding import ProteinFolder
from .synthesis import MolecularSynthesizer

__all__ = [
    "BioState",
    "bio_action",
    "step_discrete",
    "god_attractor_potential",
    "replication_entropy",
    "BioSolver",
    "EvolutionaryDynamics",
    "ProteinFolder",
    "MolecularSynthesizer",
]

__version__ = "0.1.0"
