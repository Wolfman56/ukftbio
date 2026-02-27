"""
Experiment 02 -- Levinthal Paradox (Bio Double-Slit)
UKFT choice-guided folding vs random dihedral walk.

Run:
    python experiments/02_levinthal_folding_paradox.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ukft_bio.physics import BioState, step_discrete
from ukft_bio.folding import ResidueChain

RESULTS = Path(__file__).resolve().parent.parent / "results"
RESULTS.mkdir(exist_ok=True)

N_RESIDUES   = 20
N_STEPS_UKFT = 500
N_STEPS_RAND = 5000
N_TRIALS     = 5
np.random.seed(42)


def random_fold_walk(chain, n_steps):
    """Random dihedral walk - returns rho_fold history."""
    rho_history = []
    dihedrals = chain.dihedrals.copy()
    for _ in range(n_steps):
        dihedrals += np.random.randn(len(dihedrals)) * 0.1
        chain.dihedrals = dihedrals
        rho_history.append(chain.contact_order_score())
    return rho_history


def ukft_fold(chain, n_steps):
    """UKFT least-action guided folding - returns rho_fold trajectory."""
    env = np.array([1.0, 0.5, -0.3, 0.2])
    state = chain.to_bio_state(environment=env)
    rho_history = [state.rho]
    for _ in range(n_steps):
        state = step_discrete(state, n_candidates=12, mutation_scale=0.15,
                              rho_update_rate=0.1)
        chain.dihedrals = state.genome
        rho_history.append(chain.contact_order_score())
    return rho_history


print("Exp 02 -- Levinthal Paradox")
ukft_all, rand_all = [], []
for t in range(N_TRIALS):
    chain_u = ResidueChain(N_RESIDUES)
    chain_r = ResidueChain(N_RESIDUES)
    chain_r.dihedrals = chain_u.dihedrals.copy()
    ukft_all.append(ukft_fold(chain_u, N_STEPS_UKFT))
    rand_all.append(random_fold_walk(chain_r, N_STEPS_RAND))

ukft_arr = np.array(ukft_all)
rand_arr = np.array(rand_all)

# -- Figure 1: rho trajectories --
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Exp 02 -- Levinthal Paradox: UKFT vs Random Folding", fontsize=13)

ax = axes[0]
for row in ukft_arr:
    ax.plot(row, alpha=0.4, color="steelblue", lw=1)
ax.plot(ukft_arr.mean(axis=0), color="steelblue", lw=2.5, label="UKFT mean")
ax.set_xlabel("Choice step n")
ax.set_ylabel("rho_fold (contact score)")
ax.set_title("Condition A: UKFT least-action")
ax.legend()

ax = axes[1]
for row in rand_arr:
    ax.plot(row, alpha=0.3, color="tomato", lw=1)
ax.plot(rand_arr.mean(axis=0), color="tomato", lw=2.5, label="Random mean")
ax.set_xlabel("Random walk step")
ax.set_title(f"Condition B: Random walk ({N_STEPS_RAND} steps)")
ax.legend()

ukft_max = ukft_arr.max(axis=1).mean()
rand_max = rand_arr.max(axis=1).mean()
fig.text(0.5, 0.01,
         f"Mean peak rho -- UKFT: {ukft_max:.3f}  |  Random: {rand_max:.3f}  "
         f"({'UKFT wins' if ukft_max > rand_max else 'Random wins'})",
         ha="center", fontsize=10)
plt.tight_layout()
fig.savefig(RESULTS / "exp02_folding_landscape.png", dpi=150)
plt.close()
print(f"  Fig 1 saved.  UKFT peak rho {ukft_max:.3f}  vs  Random {rand_max:.3f}")

# -- Figure 2: S_bio landscape (2D PCA) --
print("  Building action landscape (PCA slice)...")
N_SAMPLES = 2000
chain_probe = ResidueChain(N_RESIDUES)
genomes = [np.random.uniform(-np.pi, np.pi, 2 * N_RESIDUES) for _ in range(N_SAMPLES)]
rhos = []
for g in genomes:
    chain_probe.dihedrals = g
    rhos.append(chain_probe.contact_order_score())
rhos = np.array(rhos)

G = np.array(genomes)
G_c = G - G.mean(axis=0)
_, _, Vt = np.linalg.svd(G_c, full_matrices=False)
pca = G_c @ Vt[:2].T

fig2, ax2 = plt.subplots(figsize=(7, 6))
sc = ax2.scatter(pca[:, 0], pca[:, 1], c=rhos, cmap="viridis", s=8, alpha=0.6)
plt.colorbar(sc, ax=ax2, label="rho_fold (contact score)")
ax2.set_xlabel("PCA component 1 (dihedral space)")
ax2.set_ylabel("PCA component 2")
ax2.set_title("Exp 02 -- S_bio landscape (rho proxy, 2D PCA slice)")
fig2.savefig(RESULTS / "exp02_action_landscape.png", dpi=150)
plt.close()
print("  Fig 2 saved.")
print("Exp 02 DONE.")
