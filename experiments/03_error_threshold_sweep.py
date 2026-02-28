"""
Experiment 03 -- Error Threshold Sweep
Phase boundary in (fidelity, rho_final) space.

Run:
    python experiments/03_error_threshold_sweep.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ukft_bio.physics import replication_entropy
from ukft_bio.evolution import EvolutionaryDynamics

RESULTS = Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

FIDELITIES = np.linspace(0.70, 0.995, 22)
N_POP      = 50
N_GEN      = 200
N_TRIALS   = 5
GENOME_LEN = 20
ENV        = np.array([1.0, 0.5, -0.3, 0.2])
np.random.seed(7)

print("Exp 03 -- Error Threshold Sweep")
print(f"  {len(FIDELITIES)} fidelity levels x {N_TRIALS} trials x {N_POP} agents x {N_GEN} gens")

mean_rho_final, std_rho_final = [], []

for fid in FIDELITIES:
    trial_rhos = []
    for _ in range(N_TRIALS):
        dyn = EvolutionaryDynamics(
            population_size=N_POP,
            n_generations=N_GEN,
            genome_len=GENOME_LEN,
            fidelity=float(fid),
            rho_death=0.3,
            rho_god=1e4,
            mutation_scale=0.05,
            n_candidates=6,
        )
        gen_stats, _ = dyn.run(environment=ENV.copy(), verbose=False)
        if gen_stats:
            trial_rhos.append(gen_stats[-1]["mean_rho"])
        else:
            trial_rhos.append(0.0)
    mean_rho_final.append(np.mean(trial_rhos))
    std_rho_final.append(np.std(trial_rhos))
    print(f"  fidelity={fid:.3f}  mean_rho_final={mean_rho_final[-1]:.4f}")

mean_rho_final = np.array(mean_rho_final)
std_rho_final  = np.array(std_rho_final)

rho_half  = mean_rho_final.max() / 2.0
cross_idx = np.where(mean_rho_final > rho_half)[0]
f_c_est   = FIDELITIES[cross_idx[0]] if len(cross_idx) > 0 else None

# -- Figure 1: Phase diagram --
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(FIDELITIES, mean_rho_final, "o-", color="steelblue", lw=2, ms=6)
ax.fill_between(FIDELITIES,
                mean_rho_final - std_rho_final,
                mean_rho_final + std_rho_final,
                alpha=0.25, color="steelblue", label="+/-1 SD")
if f_c_est is not None:
    ax.axvline(f_c_est, color="tomato", ls="--", lw=1.8,
               label=f"f_c estimate ~ {f_c_est:.3f}")
ax.set_xlabel("Replication fidelity", fontsize=12)
ax.set_ylabel("Mean final rho_bio", fontsize=12)
ax.set_title("Exp 03 -- Error Threshold Phase Diagram", fontsize=13)
ax.legend()
fig.tight_layout()
fig.savefig(RESULTS / "exp03_phase_diagram.png", dpi=150)
plt.close()
print(f"  Fig 1 saved.  f_c estimate = {f_c_est}")

# -- Figure 2: H_rep landscape --
fid_grid   = np.linspace(0.50, 0.9999, 200)
h_rep_vals = np.array([replication_entropy(f, GENOME_LEN) for f in fid_grid])
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.plot(fid_grid, h_rep_vals, color="darkorange", lw=2)
if f_c_est is not None:
    ax2.axvline(f_c_est, color="tomato", ls="--", lw=1.5,
                label=f"f_c ~ {f_c_est:.3f}")
ax2.set_xlabel("Replication fidelity", fontsize=12)
ax2.set_ylabel(f"H_rep (L={GENOME_LEN})", fontsize=12)
ax2.set_title("Exp 03 -- Replication Entropy Cost", fontsize=12)
ax2.legend()
fig2.tight_layout()
fig2.savefig(RESULTS / "exp03_h_rep_landscape.png", dpi=150)
plt.close()
print("  Fig 2 saved.")
print("Exp 03 DONE.")

