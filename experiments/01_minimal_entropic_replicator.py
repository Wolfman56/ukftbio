"""
Experiment 01 — Minimal Entropic Replicator
=============================================

Goal
----
Demonstrate that Eigen's error threshold (1971) emerges naturally from
S_bio minimisation — without imposing it from outside.

UKFT prediction:
  When fidelity > fidelity_critical:
    Replicators converge to quasi-species attractor (stable high-ρ basin).
  When fidelity < fidelity_critical:
    Replication entropy term H_rep dominates S_bio →
    choice operator cannot maintain ρ > rho_death →
    population collapses (error catastrophe).

The critical fidelity is:
    f_c ≈ 1 − 1 / (genome_len * selection_coeff)

Experimental design (mirrors ukftphys Exp 01 free particle → bio hello world):
  - Population of 50 replicators, genome length 20
  - Two conditions: fidelity = 0.99 (above threshold), 0.80 (below threshold)
  - Track mean ρ_bio over 200 generations
  - Visualise: lineage tree in ρ-space + phase portrait

Expected output (for Gemini to validate/falsify):
  - Condition A (f=0.99): mean_ρ grows → plateau ~God-attractor pull visible
  - Condition B (f=0.80): mean_ρ collapses → population extinction by gen ~100

Run:
    conda activate ukftbio
    python experiments/01_minimal_entropic_replicator.py
"""

import sys
import os
import numpy as np

# Allow running from repo root OR experiments/ directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ukft_bio.evolution import EvolutionaryDynamics
from ukft_bio.vis import plot_rho_trajectory, plot_phase_portrait, plot_lineage_tree
from ukft_bio.solver import TrajectoryRecord


def run_condition(label: str, fidelity: float, n_gen: int = 200) -> None:
    print(f"\n{'='*60}")
    print(f"Exp 01 — Condition {label}: fidelity={fidelity}")
    print(f"{'='*60}")

    # Fixed environment: simple 4D "selection gradient"
    environment = np.array([1.0, 0.5, -0.3, 0.2])

    evo = EvolutionaryDynamics(
        population_size=50,
        n_generations=n_gen,
        genome_len=20,
        fidelity=fidelity,
        rho_death=0.3,
        rho_god=1e4,
        mutation_scale=0.05,
        n_candidates=6,
    )

    gen_stats, lineage_histories = evo.run(environment=environment, verbose=True)

    if not gen_stats:
        print(f"  Condition {label}: population went extinct immediately.")
        return

    # ── Summary ──────────────────────────────────────────────────────────────
    final = gen_stats[-1]
    print(f"\n  Final generation {final['generation']}:")
    print(f"    Population size : {final['pop_size']}")
    print(f"    Mean ρ_bio      : {final['mean_rho']:.4f}")
    print(f"    Max ρ_bio       : {final['max_rho']:.4f}")
    print(f"    H_rep (S_bio cost) : {final['replication_entropy']:.4f}")

    # ── UKFT analysis ─────────────────────────────────────────────────────────
    mean_rhos = [s["mean_rho"] for s in gen_stats]
    if len(mean_rhos) > 10:
        early_mean = np.mean(mean_rhos[:10])
        late_mean  = np.mean(mean_rhos[-10:])
        trend = "GROWING (→ God Basin)" if late_mean > early_mean else "COLLAPSING (← Error Catastrophe)"
        print(f"    ρ trend         : {trend}  (early={early_mean:.3f}, late={late_mean:.3f})")

    # ── Visualise ─────────────────────────────────────────────────────────────
    os.makedirs("results", exist_ok=True)

    # Convert gen_stats to TrajectoryRecord-like objects for plot_rho_trajectory
    records = [
        type("R", (), {"n": s["generation"], "rho": s["mean_rho"]})()
        for s in gen_stats
    ]
    plot_rho_trajectory(
        records,
        title=f"Exp 01 {label}: Mean ρ_bio  (fidelity={fidelity})",
        output_html=f"results/exp01_{label.lower()}_rho.html",
    )

    rho_vals     = [s["mean_rho"]  for s in gen_stats]
    fidelity_vals = [fidelity] * len(gen_stats)
    plot_phase_portrait(
        rho_vals, fidelity_vals,
        title=f"Exp 01 {label}: Phase Portrait",
        output_html=f"results/exp01_{label.lower()}_phase.html",
    )

    # Trim lineages for readability
    trimmed = [hist for hist in lineage_histories if len(hist) > 5][:8]
    if trimmed:
        plot_lineage_tree(
            trimmed,
            title=f"Exp 01 {label}: Lineage Trees  (fidelity={fidelity})",
            output_html=f"results/exp01_{label.lower()}_lineages.html",
        )


def main():
    print("UKFT Biology — Experiment 01: Minimal Entropic Replicator")
    print("Biological hello-world: error threshold from S_bio minimisation")

    # Condition A: above error threshold → quasi-species attractor
    run_condition("A_HighFidelity", fidelity=0.99, n_gen=200)

    # Condition B: below error threshold → error catastrophe
    run_condition("B_LowFidelity", fidelity=0.80, n_gen=200)

    print("\n\nExperiment 01 complete.  Results in results/exp01_*.html")
    print("Hand-off to Gemini for Exp 02 (entropic genetic code).")


if __name__ == "__main__":
    main()
