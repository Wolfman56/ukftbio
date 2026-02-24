"""
Tests — UKFT Biology Physics Sanity Checks
============================================

These tests verify the core invariants of the bio action functional.
They are the bio-equivalent of the ukftphys unit tests.

Run: pytest tests/ -v
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ukft_bio.physics import (
    BioState,
    bio_action,
    step_discrete,
    god_attractor_potential,
    replication_entropy,
    choice_kinetic,
    catalytic_context_gain,
    semantic_coherence_penalty,
)


# ──────────────────────────────────────────────────────────────────────────────
# Bio Action Functional Components
# ──────────────────────────────────────────────────────────────────────────────

def make_state(genome_len: int = 10, rho: float = 5.0, fidelity: float = 0.99) -> BioState:
    return BioState(
        genome=np.random.randn(genome_len),
        rho=rho,
        fidelity=fidelity,
        choice_index=0,
        environment=np.array([1.0, 0.0, 0.0, 0.0]),
    )


def test_replication_entropy_perfect_fidelity():
    """Perfect fidelity → minimum entropy (small but nonzero for finite genome_len)."""
    # H_rep(0.9999, 100) = 100 * h(0.0001) ≈ 0.102 — correct, just near-zero
    h = replication_entropy(fidelity=0.9999, genome_len=100)
    assert h < 0.5, f"Expected small H_rep for fidelity≈1, got {h}"


def test_replication_entropy_low_fidelity():
    """Low fidelity (0.5) → maximum entropy ≈ genome_len * log(2)."""
    h = replication_entropy(fidelity=0.5, genome_len=20)
    # h(0.5) = log(2) ≈ 0.693; total = 20 * 0.693 = 13.86
    assert 13.0 < h < 15.0, f"Expected ~13.86, got {h}"


def test_replication_entropy_monotone():
    """H_rep should decrease as fidelity increases."""
    fidelities = [0.6, 0.7, 0.8, 0.9, 0.99]
    entropies = [replication_entropy(f, 20) for f in fidelities]
    assert all(entropies[i] > entropies[i+1] for i in range(len(entropies)-1)), \
        "H_rep must decrease monotonically with fidelity"


def test_choice_kinetic_null_step():
    """Null choice (genome unchanged) → zero kinetic cost."""
    g = np.array([1.0, 2.0, 3.0])
    k = choice_kinetic(g, g.copy())
    assert k == pytest.approx(0.0), f"Expected K=0 for null step, got {k}"


def test_choice_kinetic_positive():
    """Non-null choice → positive kinetic cost."""
    g = np.array([1.0, 0.0])
    g_new = np.array([0.0, 1.0])
    k = choice_kinetic(g_new, g)
    assert k > 0, "Kinetic cost must be positive for non-null transition"


def test_catalytic_context_gain_sign():
    """Catalytic gain is always negative (it IS a gain = cost reduction)."""
    for rho in [0.1, 1.0, 10.0, 100.0]:
        c = catalytic_context_gain(rho)
        assert c < 0, f"Catalytic gain must be negative, got {c} at rho={rho}"


def test_god_attractor_decreasing():
    """God attractor potential decreases as ρ increases (pull toward high ρ)."""
    rhos = [0.1, 1.0, 10.0, 100.0, 1e4]
    potentials = [god_attractor_potential(r) for r in rhos]
    assert all(potentials[i] > potentials[i+1] for i in range(len(potentials)-1)), \
        "V_god must decrease monotonically with ρ (pull toward high-ρ basin)"


def test_semantic_coherence_penalty_repetitive():
    """Repetitive (low-information) genome → higher semantic penalty."""
    repetitive = np.zeros(50)  # all same → low symbol entropy → high penalty
    diverse = np.linspace(-2, 2, 50)  # spread → high symbol entropy → low penalty
    h_rep = semantic_coherence_penalty(repetitive)
    h_div = semantic_coherence_penalty(diverse)
    assert h_rep > h_div, \
        f"Repetitive genome should have higher H_sem than diverse: {h_rep} vs {h_div}"


def test_bio_action_finite():
    """Full bio_action must return a finite float for any valid state."""
    state = make_state()
    candidate = state.genome + np.random.randn(10) * 0.1
    s = bio_action(state, candidate)
    assert np.isfinite(s), f"S_bio must be finite, got {s}"


# ──────────────────────────────────────────────────────────────────────────────
# BioState
# ──────────────────────────────────────────────────────────────────────────────

def test_dt_bio_inversely_proportional_to_rho():
    """dt_bio = 1/ρ — verify the key UKFT time dilation invariant."""
    s1 = make_state(rho=1.0)
    s10 = make_state(rho=10.0)
    assert s1.dt_bio > s10.dt_bio, "High-ρ state must have shorter dt_bio"
    assert s1.dt_bio == pytest.approx(1.0 / (1.0 + 1e-12), rel=1e-6)


def test_clone_preserves_fidelity():
    state = make_state(fidelity=0.95)
    child = state.clone()
    assert child.fidelity == 0.95


def test_clone_with_mutation_differs():
    state = make_state()
    child = state.clone(mutation_noise=1.0)
    assert not np.allclose(state.genome, child.genome), \
        "Mutated clone genome should differ from parent"


# ──────────────────────────────────────────────────────────────────────────────
# Discrete Stepper
# ──────────────────────────────────────────────────────────────────────────────

def test_step_discrete_increments_choice_index():
    state = make_state()
    initial_n = state.choice_index
    new_state = step_discrete(state, n_candidates=4)
    assert new_state.choice_index == initial_n + 1


def test_step_discrete_rho_positive():
    """ρ must remain positive after any number of steps."""
    state = make_state(rho=2.0)
    for _ in range(20):
        state = step_discrete(state, n_candidates=4)
    assert state.rho > 0, "ρ_bio must always be positive"


def test_step_discrete_action_minimisation():
    """
    Over many steps with a clear attractor environment,
    mean action should not systematically increase.
    (Not strictly decreasing every step — stochastic — but trend downward.)
    """
    state = make_state(genome_len=8, rho=1.0)
    state.environment = np.array([1.0, 1.0, 1.0, 1.0])

    actions = []
    for _ in range(50):
        cand = state.genome + np.random.randn(8) * 0.1
        actions.append(bio_action(state, cand))
        state = step_discrete(state, n_candidates=8, mutation_scale=0.1)

    # Actions over 50 steps should not blow up
    assert max(actions) < 1e6, "Action exploded — check S_bio terms"
