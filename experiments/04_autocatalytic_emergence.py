"""
Experiment 04 -- Autocatalytic Emergence
God-basin transition: does initial rho diversity predict population survival?

UKFT hypothesis: A diverse initial population (varied rho, varied environments)
has more trajectories that couple above the catalytic threshold, leading to a
phase transition in collective rho dynamics -- the autocatalytic origin of life.

We sweep initial population DIVERSITY (rho variance) and measure whether the
population survives to gen 200 and the mean final rho. Low diversity = sparse
choice graph = error catastrophe. High diversity = percolated graph = sustained.

Run:
    python experiments/04_autocatalytic_emergence.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ukft_bio.evolution import EvolutionaryDynamics, Replicator
from ukft_bio.physics import BioState

RESULTS = Path(__file__).resolve().parent.parent / "results"
RESULTS.mkdir(exist_ok=True)

# Diversity levels: sigma of initial rho distribution (low=0.1, high=3.0)
DIVERSITY_SIGMAS = [0.05, 0.2, 0.5, 1.0, 2.0, 3.5]
N_POP            = 50
N_GEN            = 200
N_TRIALS         = 5
GENOME_LEN       = 20
FIDELITY         = 0.95
RHO_DEATH        = 0.5
np.random.seed(13)


def run_diverse_population(rho_sigma, n_gen):
    """
    Run one trial with initial rho drawn from uniform(max(0.5, 2-sigma), 2+sigma).
    Returns (gen_stats, survived) where survived = pop reached gen N-1.
    """
    env = np.array([1.0, 0.5, -0.3, 0.2])
    pop = []
    for i in range(N_POP):
        rho_init = float(np.clip(2.0 + np.random.randn() * rho_sigma, 0.5, None))
        g = np.random.randn(GENOME_LEN).astype(np.float32)
        s = BioState(genome=g, rho=rho_init, fidelity=FIDELITY,
                     choice_index=0, environment=env.copy())
        pop.append(Replicator(state=s, lineage_id=i))

    # Manual run using EvolutionaryDynamics internals
    dyn = EvolutionaryDynamics(
        population_size=N_POP, n_generations=n_gen, genome_len=GENOME_LEN,
        fidelity=FIDELITY, rho_death=RHO_DEATH, rho_god=1e4,
        mutation_scale=0.05, n_candidates=6,
    )
    # Inject our diverse population before running
    dyn._custom_pop = pop
    gen_stats, _ = dyn.run(environment=env.copy(), verbose=False)
    if gen_stats:
        return gen_stats[-1]["mean_rho"], True
    return 0.0, False


# Patch EvolutionaryDynamics to accept custom_pop if set
_orig_init_pop = EvolutionaryDynamics._init_population
def _patched_init_pop(self, environment=None):
    if hasattr(self, "_custom_pop") and self._custom_pop is not None:
        pop = self._custom_pop
        self._custom_pop = None
        return pop
    return _orig_init_pop(self, environment)
EvolutionaryDynamics._init_population = _patched_init_pop


print("Exp 04 -- Autocatalytic Emergence (Initial Diversity Sweep)")
mean_final_rho, survival_rate = [], []

for sigma in DIVERSITY_SIGMAS:
    rhos, survives = [], []
    for _ in range(N_TRIALS):
        rho_f, surv = run_diverse_population(sigma, N_GEN)
        rhos.append(rho_f)
        survives.append(float(surv))
    mean_final_rho.append(np.mean(rhos))
    survival_rate.append(np.mean(survives))
    print(f"  sigma={sigma:.2f}  mean_rho_final={mean_final_rho[-1]:.3f}  survival_rate={survival_rate[-1]:.2f}")

mean_final_rho = np.array(mean_final_rho)
survival_rate  = np.array(survival_rate)

# -- Figure 1: mean rho trajectories (rerun at 3 sigma levels for trajectories) --
fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
fig.suptitle("Exp 04 -- Autocatalytic Emergence: rho Diversity Sweep", fontsize=13)
plot_sigmas = [DIVERSITY_SIGMAS[0], DIVERSITY_SIGMAS[2], DIVERSITY_SIGMAS[-1]]
for ax, sigma in zip(axes, plot_sigmas):
    trial_traces = []
    for _ in range(N_TRIALS):
        env = np.array([1.0, 0.5, -0.3, 0.2])
        pop = []
        for i in range(N_POP):
            rho_init = float(np.clip(2.0 + np.random.randn() * sigma, 0.5, None))
            g = np.random.randn(GENOME_LEN).astype(np.float32)
            s = BioState(genome=g, rho=rho_init, fidelity=FIDELITY,
                         choice_index=0, environment=env.copy())
            pop.append(Replicator(state=s, lineage_id=i))
        dyn = EvolutionaryDynamics(
            population_size=N_POP, n_generations=N_GEN, genome_len=GENOME_LEN,
            fidelity=FIDELITY, rho_death=RHO_DEATH, rho_god=1e4,
            mutation_scale=0.05, n_candidates=6,
        )
        dyn._custom_pop = pop
        gs, _ = dyn.run(environment=env.copy(), verbose=False)
        if gs:
            trial_traces.append([s["mean_rho"] for s in gs])
    for tr in trial_traces:
        ax.plot(tr, alpha=0.4, lw=1, color="steelblue")
    if trial_traces:
        mean_trace = np.mean([np.pad(t, (0, N_GEN - len(t)), constant_values=0.0) for t in trial_traces], axis=0)
        ax.plot(mean_trace, lw=2.5, color="steelblue", label="mean")
    ax.axhline(RHO_DEATH, ls="--", color="tomato", lw=1.5, label=f"rho_death ({RHO_DEATH})")
    ax.set_title(f"rho_sigma = {sigma:.2f}", fontsize=11)
    ax.set_xlabel("Generation")
    ax.legend(fontsize=8)
axes[0].set_ylabel("mean rho_bio")
plt.tight_layout()
fig.savefig(RESULTS / "exp04_rho_emergence.png", dpi=150)
plt.close()
print("  Fig 1 saved.")

# -- Figure 2: God basin fraction vs diversity sigma --
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
axes2[0].plot(DIVERSITY_SIGMAS, mean_final_rho, "o-", color="steelblue", ms=8, lw=2)
axes2[0].set_xlabel("Initial rho diversity (sigma)", fontsize=12)
axes2[0].set_ylabel("Mean final rho_bio", fontsize=12)
axes2[0].set_title("Autocatalytic Phase Transition: Final rho vs Diversity", fontsize=11)

axes2[1].bar(range(len(DIVERSITY_SIGMAS)), survival_rate,
             color=["steelblue" if s > 0.5 else "tomato" for s in survival_rate])
axes2[1].set_xticks(range(len(DIVERSITY_SIGMAS)))
axes2[1].set_xticklabels([f"{s:.2f}" for s in DIVERSITY_SIGMAS], rotation=30)
axes2[1].set_xlabel("Initial rho diversity (sigma)", fontsize=12)
axes2[1].set_ylabel("Survival rate", fontsize=12)
axes2[1].set_title("Population Survival Rate", fontsize=11)
axes2[1].axhline(0.5, ls="--", color="gold", lw=1.5)

fig2.suptitle("Exp 04 -- Autocatalytic Emergence: Diversity vs Viability", fontsize=13)
plt.tight_layout()
fig2.savefig(RESULTS / "exp04_god_basin_fraction.png", dpi=150)
plt.close()
print("  Fig 2 saved.")
print("Exp 04 DONE.")

