"""
Experiment 05 -- Fungal Choice Spike (Entangled Collapse)
Synthetic 4-channel spike model testing Paper 01 Predictions P2 & P3.

Run:
    python experiments/05_fungal_choice_spike.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RESULTS = Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

N_CHANNELS = 4
N_CLUSTERS = 200
PROP_SD    = 0.15   # propagation jitter (fraction of root amplitude)
np.random.seed(99)

# UKFT action minimum in 4D feature space:
# (amp_mean, timing_tightness, inter_channel_corr, amp_std)
ACTION_MINIMUM = np.array([0.6, 0.6, -0.3, 0.1])


def _feature(amps, timing_spread):
    """Compute 4D feature vector for one cluster."""
    all_amps = np.array(amps, dtype=float)
    amp_mean = float(all_amps.mean())
    ref = np.ones(len(all_amps)) / len(all_amps)
    r, _ = stats.pearsonr(all_amps, ref + np.random.randn(len(all_amps)) * 0.01)
    r = float(np.clip(r, -1, 1))
    return np.array([amp_mean, 1.0 - timing_spread, r, float(all_amps.std())])


def make_entangled(resolved=True):
    """Generate one synthetic spike cluster feature vector."""
    root_amp = np.random.exponential(1.0)
    if resolved:
        delay_noise  = np.random.randn(N_CHANNELS - 1) * PROP_SD
        partner_amps = root_amp * (1.0 - np.abs(delay_noise)) + np.random.randn(N_CHANNELS - 1) * 0.05
        timing_spread = float(np.abs(np.random.randn()) * 0.3)
    else:
        n_fired      = np.random.randint(1, N_CHANNELS)
        partner_amps = np.zeros(N_CHANNELS - 1)
        partner_amps[:n_fired] = np.random.exponential(0.5, size=n_fired)
        timing_spread = float(0.8 + np.random.rand() * 0.2)

    all_amps = [root_amp] + list(partner_amps)
    feat = _feature(all_amps, timing_spread)
    am_n = ACTION_MINIMUM / np.linalg.norm(ACTION_MINIMUM)
    fn   = feat / (np.linalg.norm(feat) + 1e-9)
    return {"corr": feat[2], "feature": feat,
            "action_bias": float(np.dot(fn, am_n)), "resolved": resolved}


def make_shuffled():
    """Shuffled control: random amps with no propagation structure."""
    amps = np.random.exponential(1.0, size=N_CHANNELS)
    timing_spread = float(np.random.rand())
    feat = _feature(amps, timing_spread)
    am_n = ACTION_MINIMUM / np.linalg.norm(ACTION_MINIMUM)
    fn   = feat / (np.linalg.norm(feat) + 1e-9)
    return {"corr": feat[2], "feature": feat,
            "action_bias": float(np.dot(fn, am_n)), "resolved": True}


print("Exp 05 -- Fungal Choice Spike")
real_res = [make_entangled(resolved=True)  for _ in range(N_CLUSTERS)]
real_par = [make_entangled(resolved=False) for _ in range(N_CLUSTERS)]
shuffled = [make_shuffled()                for _ in range(N_CLUSTERS)]

# -- P2: inter-channel correlation t-test --
real_corr = np.array([c["corr"] for c in real_res])
shuf_corr = np.array([c["corr"] for c in shuffled])
t2, p2 = stats.ttest_ind(real_corr, shuf_corr)
print(f"  P2  t={t2:.3f}  p={p2:.4f}  real_mean_r={real_corr.mean():.3f}  shuffled_mean_r={shuf_corr.mean():.3f}")

# -- P3: distance to action minimum --
am_n = ACTION_MINIMUM / np.linalg.norm(ACTION_MINIMUM)

def dist_to_min(clusters):
    return np.array([
        np.linalg.norm(c["feature"] / (np.linalg.norm(c["feature"]) + 1e-9) - am_n)
        for c in clusters
    ])

dist_res = dist_to_min(real_res)
dist_par = dist_to_min(real_par)
dist_shu = dist_to_min(shuffled)
t3, p3 = stats.ttest_ind(dist_res, dist_shu)
print(f"  P3  dist: resolved={dist_res.mean():.3f}  partial={dist_par.mean():.3f}  shuffled={dist_shu.mean():.3f}")
print(f"      t={t3:.3f}  p={p3:.4f}")

# -- Figure 1: PCA scatter in feature space --
all_feats = np.vstack([c["feature"] for c in real_res + real_par + shuffled])
labels    = ["real-resolved"] * N_CLUSTERS + ["real-partial"] * N_CLUSTERS + ["shuffled"] * N_CLUSTERS
G_c = all_feats - all_feats.mean(axis=0)
_, _, Vt = np.linalg.svd(G_c, full_matrices=False)
pca = G_c @ Vt[:2].T

fig, ax = plt.subplots(figsize=(8, 7))
for label, colour, marker in [
    ("real-resolved", "steelblue", "o"),
    ("real-partial",  "skyblue",   "^"),
    ("shuffled",      "tomato",    "s"),
]:
    idx = [i for i, l in enumerate(labels) if l == label]
    ax.scatter(pca[idx, 0], pca[idx, 1], c=colour, marker=marker,
               s=18, alpha=0.6, label=label)

am_proj = (ACTION_MINIMUM - all_feats.mean(axis=0)) @ Vt[:2].T
ax.scatter(*am_proj, marker="*", s=300, color="gold", zorder=10, label="Action minimum")
ax.set_xlabel("PCA 1")
ax.set_ylabel("PCA 2")
ax.set_title(
    f"Exp 05 -- Cluster Feature Space\n"
    f"P2: corr t={t2:.2f} p={p2:.3f}  |  P3: dist t={t3:.2f} p={p3:.3f}"
)
ax.legend(markerscale=1.5, fontsize=9)
fig.tight_layout()
fig.savefig(RESULTS / "exp05_cluster_resolution.png", dpi=150)
plt.close()
print("  Fig 1 saved.")

# -- Figure 2: correlation histogram (P2) --
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.hist(real_corr, bins=25, alpha=0.65, color="steelblue", label="Real entangled")
ax2.hist(shuf_corr, bins=25, alpha=0.65, color="tomato",    label="Shuffled control")
ax2.axvline(real_corr.mean(), color="steelblue", ls="--", lw=2,
            label=f"Real mean r={real_corr.mean():.3f}")
ax2.axvline(shuf_corr.mean(), color="tomato",    ls="--", lw=2,
            label=f"Shuffled mean r={shuf_corr.mean():.3f}")
ax2.set_xlabel("Inter-channel correlation r")
ax2.set_ylabel("Count")
ax2.set_title(f"Exp 05 -- P2 Test: t={t2:.2f}  p={p2:.4f}")
ax2.legend()
fig2.tight_layout()
fig2.savefig(RESULTS / "exp05_p2_correlation.png", dpi=150)
plt.close()
print("  Fig 2 saved.")
print("Exp 05 DONE.")
