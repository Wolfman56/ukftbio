"""
Experiment 19 — CV-Stratified JEPA Transfer
=============================================

Exp 16 established that Pleurotus → multispecies JEPA transfer quality (η)
tracks AR1 similarity rather than ecology or phylogeny (all η < 0.70).
Exp 16's temporal structure analysis also revealed a stark CV split:

    Bursty tier  (CV > 0.5) : Pleurotus (2.03), Schizophyllum (1.53)
    Smooth tier  (CV ≤ 0.5) : Cordyceps (0.21), Enoki (0.21), Omphalotus (0.06)

Exp 16 only tested one source (Pleurotus).  Exp 19 runs a full N×N transfer
matrix: every species trains a TinyJEPA, applies it frozen to every other
species, and computes η = native_S / transfer_S.

Central question: is CV the true grammar boundary for JEPA transferability?

HYPOTHESES
   H19a:  Intra-tier η (bursty→bursty, smooth→smooth) is significantly higher
          than cross-tier η (bursty→smooth, smooth→bursty).
   H19b:  Across all directed pairs, η anti-correlates more strongly with
          |ΔCV| than with |ΔAR1| — CV is the dominant grammar axis, not AR1.
   H19c:  Pleurotus → Schizophyllum η is the single highest transfer value
          across the 5×4 = 20 directed pairs (both high-CV bursty).

Figures
-------
  19_cv_tier_scatter.png      AR1 vs CV scatter, species labelled, tiers coloured
  19_transfer_matrix.png      N×N heatmap of η (source rows, target cols)
  19_intra_vs_cross_tier.png  Boxplot: intra-tier η vs cross-tier η + Mann-Whitney
  19_cv_ar1_boundary.png      η vs |ΔCV| and η vs |ΔAR1| scatter (2-panel)
  19_report.json              Full η matrix + hypothesis test results
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.ndimage import uniform_filter1d

# ── Paths ─────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
PLEU_EMB    = MYCO / "tools/data/pleurotus_spike_emb_40d.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants ──────────────────────────────────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600
DT_S             = 1.0

# ── JEPA parameters (consistent with Exps 15–18) ─────────────────────────────
SEQ_LEN      = 16
PRED_HORIZON = 4
HIDDEN_DIM   = 32
LR           = 3e-3
N_EPOCHS     = 60
MIN_WINDOWS  = 30
DEGENERATE   = 1e-4

# ── Species metadata (from Exps 14–16) ───────────────────────────────────────
SPECIES_META = {
    "Pleurotus":     {"ar1": 0.769, "cv": 2.03, "stress": 0.35, "evo_mya": None},
    "Schizophyllum": {"ar1": 0.944, "cv": 1.53, "stress": 0.90, "evo_mya": 250.0},
    "Cordyceps":     {"ar1": 0.911, "cv": 0.21, "stress": 0.70, "evo_mya": 800.0},
    "Enoki":         {"ar1": 0.595, "cv": 0.21, "stress": 0.50, "evo_mya": 100.0},
    "Omphalotus":    {"ar1": 0.617, "cv": 0.06, "stress": 0.20, "evo_mya": 150.0},
}

CV_THRESHOLD  = 0.5
TIER_BURSTY   = {sp for sp, m in SPECIES_META.items() if m["cv"] > CV_THRESHOLD}
TIER_SMOOTH   = {sp for sp, m in SPECIES_META.items() if m["cv"] <= CV_THRESHOLD}

MULTISPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
    "Enoki":         "Enoki fungi Flammulina velutipes.txt",
    "Omphalotus":    "Ghost Fungi Omphalotus nidiformis.txt",
}

SPECIES_COLOURS = {
    "Pleurotus":     "#2ecc71",
    "Schizophyllum": "#9b59b6",
    "Cordyceps":     "#e67e22",
    "Enoki":         "#27ae60",
    "Omphalotus":    "#2980b9",
}

TIER_COLOURS = {"bursty": "#e74c3c", "smooth": "#3498db"}


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADERS
# ─────────────────────────────────────────────────────────────────────────────

def load_pleurotus_rho() -> np.ndarray:
    """Combined density from pleurotus_spike_emb_40d.ndjson (rho field)."""
    records = [json.loads(l) for l in open(PLEU_EMB)]
    records.sort(key=lambda r: r["t_start_s"])
    rho = np.array([r["rho"] for r in records], dtype=np.float32)
    print(f"  Pleurotus: {len(rho)} windows from NDJSON  rho mean={rho.mean():.4f}")
    return rho


def load_multispecies_rho(species: str) -> np.ndarray:
    """Combined (rho_local + rho_global) density series from raw .txt."""
    path = DATA_DIR / MULTISPECIES_FILES[species]
    df = pd.read_csv(str(path), sep="\t", header=0, dtype=np.float32,
                     na_values=["NaN", "", " "], engine="c", on_bad_lines="skip")
    df.dropna(how="all", inplace=True)
    df = df.select_dtypes(include=[np.floating, np.integer])
    arr = df.values.astype(np.float32)
    n, nchan = arr.shape

    det = np.empty_like(arr)
    for ch in range(nchan):
        col = arr[:, ch].astype(np.float64)
        nm  = np.isnan(col)
        if nm.any():
            ix = np.arange(n); col[nm] = np.interp(ix[nm], ix[~nm], col[~nm])
        det[:, ch] = (col - uniform_filter1d(col, size=600, mode="mirror")).astype(np.float32)

    abs_v  = np.abs(det)
    sigma  = np.nanmedian(abs_v, axis=0) / 0.6745
    floor  = (NOISE_FLOOR_K * sigma).astype(np.float32)
    noise  = abs_v < floor[np.newaxis, :]
    glob_m = abs_v >= GLOBAL_THRESH_MV
    loc_m  = (~noise) & (~glob_m)
    n_win  = n // WINDOW_S
    tr_l   = loc_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    tr_g   = glob_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    rho    = (tr_l.mean(axis=(1, 2)) + tr_g.mean(axis=(1, 2))).astype(np.float32)
    print(f"  {species}: {n:,} rows → {n_win} windows  rho mean={rho.mean():.4f}")
    return rho


# ─────────────────────────────────────────────────────────────────────────────
# TINY JEPA (consistent with Exps 15–18)
# ─────────────────────────────────────────────────────────────────────────────

class TinyJEPA:
    def __init__(self):
        np.random.seed(42)
        s1 = np.sqrt(2.0 / SEQ_LEN); s2 = np.sqrt(2.0 / HIDDEN_DIM)
        self.We = np.random.randn(HIDDEN_DIM, SEQ_LEN).astype(np.float32) * s1
        self.be = np.zeros(HIDDEN_DIM, dtype=np.float32)
        self.Wp = np.random.randn(HIDDEN_DIM, HIDDEN_DIM).astype(np.float32) * s2
        self.bp = np.zeros(HIDDEN_DIM, dtype=np.float32)
        self.Wo = np.random.randn(1, HIDDEN_DIM).astype(np.float32) * s2
        self.bo = np.zeros(1, dtype=np.float32)

    # Deepcopy weights for frozen transfer
    def clone(self) -> "TinyJEPA":
        other = TinyJEPA.__new__(TinyJEPA)
        for attr in ("We", "be", "Wp", "bp", "Wo", "bo"):
            setattr(other, attr, getattr(self, attr).copy())
        return other

    def _enc(self, x): return np.maximum(0, self.We @ x + self.be)
    def _pred(self, h): return np.maximum(0, self.Wp @ h + self.bp)
    def _dec(self, h): return float((self.Wo @ h + self.bo)[0])

    def _step(self, ctx, tgt):
        h_e = self._enc(ctx); h_p = self._pred(h_e)
        p = self._dec(h_p); e = p - tgt
        dWo = e * h_p[np.newaxis, :]; dbo = np.array([e])
        dp = (self.Wo.T * e).flatten() * (h_p > 0)
        dWp = np.outer(dp, h_e); dbp = dp
        de = (self.Wp.T @ dp) * (h_e > 0)
        dWe = np.outer(de, ctx); dbe = de
        r = 1e-4
        self.Wo -= LR * (dWo + r * self.Wo); self.bo -= LR * dbo
        self.Wp -= LR * (dWp + r * self.Wp); self.bp -= LR * dbp
        self.We -= LR * (dWe + r * self.We); self.be -= LR * dbe

    def _norm(self, s): return (s - s.mean()) / (s.std() + 1e-8)

    def train(self, series: np.ndarray) -> None:
        s = self._norm(series)
        idx = list(range(len(s) - SEQ_LEN - PRED_HORIZON))
        for _ in range(N_EPOCHS):
            np.random.shuffle(idx)
            for i in idx:
                self._step(s[i:i+SEQ_LEN], float(s[i+SEQ_LEN+PRED_HORIZON-1]))

    def mean_surprise(self, series: np.ndarray) -> float:
        """Mean squared prediction error over the series."""
        s = self._norm(series)
        n = len(s)
        if n < SEQ_LEN + PRED_HORIZON + 1:
            return float("nan")
        errs = [(self._dec(self._pred(self._enc(s[i:i+SEQ_LEN])))
                 - float(s[i+SEQ_LEN+PRED_HORIZON-1])) ** 2
                for i in range(n - SEQ_LEN - PRED_HORIZON)]
        return float(np.mean(errs))


def train_jepa(rho: np.ndarray, label: str) -> TinyJEPA | None:
    """Train TinyJEPA on a density series.  Returns None if degenerate."""
    if len(rho) < MIN_WINDOWS or rho.std() < DEGENERATE:
        print(f"    {label}: degenerate (n={len(rho)}, std={rho.std():.4f}) → skip")
        return None
    model = TinyJEPA()
    model.train(rho)
    native_s = model.mean_surprise(rho)
    print(f"    {label}: n={len(rho)}  native_S={native_s:.4f}")
    return model


def transfer_surprise(model: TinyJEPA, target_rho: np.ndarray) -> float:
    """Apply a frozen model to target series; return mean surprise."""
    frozen = model.clone()
    return frozen.mean_surprise(target_rho)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_cv_tier_scatter() -> None:
    fig, ax = plt.subplots(figsize=(8, 6))

    for sp, meta in SPECIES_META.items():
        tier = "bursty" if meta["cv"] > CV_THRESHOLD else "smooth"
        c    = SPECIES_COLOURS[sp]
        ax.scatter(meta["ar1"], meta["cv"], s=200, color=c,
                   edgecolor="black", linewidth=1.2, zorder=5)
        ax.annotate(sp, (meta["ar1"], meta["cv"]),
                    xytext=(6, 4), textcoords="offset points", fontsize=10)

    ax.axhline(CV_THRESHOLD, color="gray", linestyle="--",
               linewidth=1, label=f"CV threshold = {CV_THRESHOLD}")
    ax.set_xlabel("AR1 (lag-1 autocorrelation)", fontsize=12)
    ax.set_ylabel("CV (coefficient of variation)", fontsize=12)
    ax.set_title("CV-Tier Classification\nBursty tier: CV > 0.5  |  Smooth tier: CV ≤ 0.5",
                 fontsize=11)

    bursty_patch = mpatches.Patch(color=TIER_COLOURS["bursty"], alpha=0.3,
                                   label=f"Bursty: {', '.join(sorted(TIER_BURSTY))}")
    smooth_patch = mpatches.Patch(color=TIER_COLOURS["smooth"], alpha=0.3,
                                   label=f"Smooth: {', '.join(sorted(TIER_SMOOTH))}")
    ax.legend(handles=[bursty_patch, smooth_patch,
                        plt.Line2D([0], [0], color="gray", linestyle="--",
                                   label=f"CV threshold = {CV_THRESHOLD}")],
              fontsize=9, loc="upper left")

    # Shade tier regions
    ax.axhspan(CV_THRESHOLD, ax.get_ylim()[1] if ax.get_ylim()[1] > CV_THRESHOLD else 3.0,
               alpha=0.07, color=TIER_COLOURS["bursty"])
    ax.axhspan(0, CV_THRESHOLD, alpha=0.07, color=TIER_COLOURS["smooth"])
    ax.set_ylim(-0.1, max(m["cv"] for m in SPECIES_META.values()) * 1.15)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "19_cv_tier_scatter.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 19_cv_tier_scatter.png")


def fig_transfer_matrix(species: list[str], eta_matrix: np.ndarray) -> None:
    n = len(species)
    fig, ax = plt.subplots(figsize=(8, 7))
    masked = np.ma.masked_where(np.eye(n, dtype=bool), eta_matrix)
    im = ax.imshow(masked, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto")
    plt.colorbar(im, ax=ax, label="η = native_S / transfer_S", fraction=0.04)

    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(eta_matrix[i, j]):
                val = eta_matrix[i, j]
                colour = "black" if 0.25 < val < 0.85 else "white"
                text = f"{val:.3f}"
                ax.text(j, i, text, ha="center", va="center",
                        fontsize=9, color=colour, fontweight="bold")
            elif i == j:
                ax.text(j, i, "—", ha="center", va="center", fontsize=10, color="gray")

    # Draw tier boundary boxes
    tiers = ["bursty" if sp in TIER_BURSTY else "smooth" for sp in species]
    for i in range(n):
        for j in range(n):
            if i != j:
                same_tier = tiers[i] == tiers[j]
                if same_tier:
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                               fill=False, edgecolor="white",
                                               linewidth=2.5, zorder=10))

    labels = [f"{sp}\n({'B' if sp in TIER_BURSTY else 'S'})" for sp in species]
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Target species", fontsize=11)
    ax.set_ylabel("Source species (model trained here)", fontsize=11)
    ax.set_title(
        "JEPA Transfer Efficiency Matrix\n"
        "η = native_S / transfer_S  (white border = same CV tier)",
        fontsize=11,
    )
    plt.tight_layout()
    out = str(RESULTS_DIR / "19_transfer_matrix.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 19_transfer_matrix.png")


def fig_intra_vs_cross_tier(eta_matrix: np.ndarray, species: list[str],
                             mw_stat: float, mw_p: float) -> None:
    n = len(species)
    intra_eta: list[float] = []
    cross_eta:  list[float] = []

    for i in range(n):
        for j in range(n):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            si = species[i]; sj = species[j]
            same = (si in TIER_BURSTY and sj in TIER_BURSTY) or \
                   (si in TIER_SMOOTH and sj in TIER_SMOOTH)
            if same:
                intra_eta.append(eta_matrix[i, j])
            else:
                cross_eta.append(eta_matrix[i, j])

    fig, ax = plt.subplots(figsize=(7, 5))
    parts = ax.violinplot([intra_eta, cross_eta], positions=[0, 1],
                          showmedians=True, showextrema=True)
    colours = [TIER_COLOURS["bursty"], "#888"]
    for body, c in zip(parts["bodies"], colours):
        body.set_facecolor(c); body.set_alpha(0.7)
    parts["cmedians"].set_color("black")
    parts["cmaxes"].set_color("black")
    parts["cmins"].set_color("black")
    parts["cbars"].set_color("black")

    # Overlay individual points
    np.random.seed(0)
    for xi, vals in [(0, intra_eta), (1, cross_eta)]:
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(np.full(len(vals), xi) + jitter, vals,
                   color="black", s=30, alpha=0.6, zorder=5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"Intra-tier\n(n={len(intra_eta)} pairs)", f"Cross-tier\n(n={len(cross_eta)} pairs)"],
        fontsize=11,
    )
    ax.set_ylabel("η = native_S / transfer_S", fontsize=11)
    ax.set_title(
        f"JEPA Transfer Efficiency: Intra-tier vs Cross-tier\n"
        f"Mann-Whitney U: p = {mw_p:.4f}  "
        f"{'H19a ✓ confirmed' if mw_p < 0.05 else 'H19a ✗ not confirmed'}",
        fontsize=11,
    )
    for xi, vals, c in [(0, intra_eta, TIER_COLOURS["bursty"]),
                         (1, cross_eta, "#888")]:
        ax.text(xi, max(vals) + 0.02, f"med={np.median(vals):.3f}",
                ha="center", fontsize=9, color=c)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "19_intra_vs_cross_tier.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 19_intra_vs_cross_tier.png")


def fig_cv_ar1_boundary(eta_matrix: np.ndarray, species: list[str]) -> None:
    n = len(species)
    delta_cv:  list[float] = []
    delta_ar1: list[float] = []
    eta_vals:  list[float] = []
    pair_labels:list[str]  = []

    for i in range(n):
        for j in range(n):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            si = species[i]; sj = species[j]
            delta_cv.append(abs(SPECIES_META[si]["cv"] - SPECIES_META[sj]["cv"]))
            delta_ar1.append(abs(SPECIES_META[si]["ar1"] - SPECIES_META[sj]["ar1"]))
            eta_vals.append(eta_matrix[i, j])
            pair_labels.append(f"{si[:4]}→{sj[:4]}")

    delta_cv  = np.array(delta_cv)
    delta_ar1 = np.array(delta_ar1)
    eta_vals  = np.array(eta_vals)

    r_cv,  p_cv  = stats.spearmanr(delta_cv,  eta_vals)
    r_ar1, p_ar1 = stats.spearmanr(delta_ar1, eta_vals)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, dx, xlabel, r, p, label in [
        (axes[0], delta_cv,  "|ΔCV|",  r_cv,  p_cv,  "ΔCV"),
        (axes[1], delta_ar1, "|ΔAR1|", r_ar1, p_ar1, "ΔAR1"),
    ]:
        colours_pts = [SPECIES_COLOURS.get(species[i], "#888")
                       for i in range(n) for j in range(n)
                       if i != j and not np.isnan(eta_matrix[i, j])]
        ax.scatter(dx, eta_vals, c=colours_pts, s=60, alpha=0.8, edgecolor="black",
                   linewidth=0.5, zorder=5)

        # Regression line
        slope, intercept, *_ = stats.linregress(dx, eta_vals)
        xs = np.linspace(dx.min(), dx.max(), 100)
        ax.plot(xs, slope * xs + intercept, "k--", linewidth=1.5, alpha=0.7)

        for i, lbl in enumerate(pair_labels):
            ax.annotate(lbl, (dx[i], eta_vals[i]),
                        xytext=(3, 2), textcoords="offset points", fontsize=6, alpha=0.7)

        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("η" if ax is axes[0] else "", fontsize=12)
        ax.set_title(f"η vs {label}\nSpearman rs = {r:.3f}  p = {p:.4f}",
                     fontsize=11)
        ax.grid(True, alpha=0.3)

    plt.suptitle(
        "Which axis governs JEPA transfer quality?\n"
        f"H19b: |ΔCV| stronger than |ΔAR1|  →  "
        f"{'✓ confirmed' if abs(r_cv) > abs(r_ar1) else '✗ not confirmed'}",
        fontsize=11, y=1.02,
    )
    plt.tight_layout()
    out = str(RESULTS_DIR / "19_cv_ar1_boundary.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 19_cv_ar1_boundary.png")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 19 — CV-Stratified JEPA Transfer")
    print("=" * 70)
    print(f"\n  Bursty tier (CV > {CV_THRESHOLD}): {sorted(TIER_BURSTY)}")
    print(f"  Smooth tier (CV ≤ {CV_THRESHOLD}): {sorted(TIER_SMOOTH)}")

    # ── 1. Load density series ────────────────────────────────────────────
    print("\n--- Step 1: Load density series ---")
    rho: dict[str, np.ndarray] = {}
    rho["Pleurotus"] = load_pleurotus_rho()
    for sp in MULTISPECIES_FILES:
        print(f"  {sp} …", end=" ", flush=True)
        rho[sp] = load_multispecies_rho(sp)
    species = list(rho.keys())   # consistent order
    print(f"\n  Loaded {len(species)} species")

    # ── 2. Train one TinyJEPA per species (native) ────────────────────────
    print("\n--- Step 2: Train native TinyJEPA per species ---")
    models: dict[str, TinyJEPA | None] = {}
    native_S: dict[str, float] = {}
    for sp in species:
        print(f"  [{sp}]")
        m = train_jepa(rho[sp], sp)
        models[sp] = m
        if m is not None:
            native_S[sp] = m.mean_surprise(rho[sp])
        else:
            native_S[sp] = float("nan")

    # ── 3. Transfer matrix ────────────────────────────────────────────────
    print("\n--- Step 3: Transfer matrix (N×N) ---")
    n = len(species)
    transfer_S = np.full((n, n), np.nan)
    eta_matrix  = np.full((n, n), np.nan)

    for i, src in enumerate(species):
        if models[src] is None:
            continue
        for j, tgt in enumerate(species):
            if i == j:
                continue
            ts = transfer_surprise(models[src], rho[tgt])
            transfer_S[i, j] = ts
            ns = native_S[tgt]
            if not np.isnan(ns) and ts > 0:
                eta_matrix[i, j] = ns / ts
            print(f"  {src:15s} → {tgt:15s}  transfer_S={ts:.4f}  "
                  f"η={eta_matrix[i, j]:.4f}" if not np.isnan(eta_matrix[i, j])
                  else f"  {src:15s} → {tgt:15s}  transfer_S={ts:.4f}  η=N/A",
                  flush=True)

    # ── 4. Hypothesis tests ───────────────────────────────────────────────
    print("\n--- Step 4: Hypothesis tests ---")

    # H19a: intra-tier vs cross-tier η (Mann-Whitney)
    intra_eta: list[float] = []
    cross_eta:  list[float] = []
    for i in range(n):
        for j in range(n):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            si = species[i]; sj = species[j]
            same = (si in TIER_BURSTY and sj in TIER_BURSTY) or \
                   (si in TIER_SMOOTH and sj in TIER_SMOOTH)
            (intra_eta if same else cross_eta).append(eta_matrix[i, j])

    if intra_eta and cross_eta:
        mw_stat, mw_p = stats.mannwhitneyu(intra_eta, cross_eta,
                                            alternative="greater")
    else:
        mw_stat, mw_p = float("nan"), float("nan")

    print(f"  H19a — intra-tier η > cross-tier η:")
    print(f"    intra  n={len(intra_eta)}  median={np.median(intra_eta):.4f}")
    print(f"    cross  n={len(cross_eta)}  median={np.median(cross_eta):.4f}")
    print(f"    Mann-Whitney U={mw_stat:.1f}  p={mw_p:.4f}  "
          f"→ {'✓ CONFIRMED' if mw_p < 0.05 else '✗ NOT CONFIRMED'}")

    # H19b: |ΔCV| vs |ΔAR1| correlation with η
    delta_cv_all:  list[float] = []
    delta_ar1_all: list[float] = []
    eta_all:       list[float] = []
    for i in range(n):
        for j in range(n):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            si = species[i]; sj = species[j]
            delta_cv_all.append(abs(SPECIES_META[si]["cv"] - SPECIES_META[sj]["cv"]))
            delta_ar1_all.append(abs(SPECIES_META[si]["ar1"] - SPECIES_META[sj]["ar1"]))
            eta_all.append(eta_matrix[i, j])

    r_cv,  p_cv  = stats.spearmanr(delta_cv_all,  eta_all)
    r_ar1, p_ar1 = stats.spearmanr(delta_ar1_all, eta_all)
    h19b = abs(r_cv) > abs(r_ar1)
    print(f"\n  H19b — |ΔCV| stronger predictor of η than |ΔAR1|:")
    print(f"    |ΔCV|  rs={r_cv:.4f}  p={p_cv:.4f}")
    print(f"    |ΔAR1| rs={r_ar1:.4f}  p={p_ar1:.4f}")
    print(f"    → {'✓ CONFIRMED (|ΔCV| stronger)' if h19b else '✗ NOT CONFIRMED (|ΔAR1| stronger)'}")

    # H19c: Pleurotus → Schizophyllum is highest single-pair η
    pair_etas = [(species[i], species[j], float(eta_matrix[i, j]))
                 for i in range(n) for j in range(n)
                 if i != j and not np.isnan(eta_matrix[i, j])]
    pair_etas.sort(key=lambda x: -x[2])
    print(f"\n  H19c — Pleurotus → Schizophyllum is highest η:")
    print(f"    Top-5 pairs:")
    for src, tgt, e in pair_etas[:5]:
        marker = "← ★" if src == "Pleurotus" and tgt == "Schizophyllum" else ""
        print(f"    {src:15s} → {tgt:15s}  η={e:.4f}  {marker}")
    top_src, top_tgt, top_e = pair_etas[0]
    h19c = (top_src == "Pleurotus" and top_tgt == "Schizophyllum")
    print(f"    Highest pair: {top_src} → {top_tgt}  η={top_e:.4f}")
    print(f"    → {'✓ CONFIRMED' if h19c else '✗ NOT CONFIRMED'}")

    # ── 5. Figures ────────────────────────────────────────────────────────
    print("\n--- Step 5: Figures ---")
    fig_cv_tier_scatter()
    fig_transfer_matrix(species, eta_matrix)
    fig_intra_vs_cross_tier(eta_matrix, species, mw_stat, mw_p)
    fig_cv_ar1_boundary(eta_matrix, species)

    # ── 6. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 19 SUMMARY")
    print("=" * 70)
    print(f"\n  Species: {species}")
    print(f"\n  Native JEPA surprise:")
    for sp in species:
        tier = "bursty" if sp in TIER_BURSTY else "smooth"
        print(f"    {sp:15s}  native_S={native_S[sp]:.4f}  "
              f"AR1={SPECIES_META[sp]['ar1']:.3f}  CV={SPECIES_META[sp]['cv']:.2f}  "
              f"[{tier}]")

    print(f"\n  η matrix (source → target):")
    header = "   " + " ".join(f"{s[:8]:>8}" for s in species)
    print(header)
    for i, src in enumerate(species):
        row = f"  {src[:8]:8s}"
        for j in range(n):
            if i == j:
                row += "       — "
            elif np.isnan(eta_matrix[i, j]):
                row += "      N/A"
            else:
                row += f"  {eta_matrix[i, j]:7.4f}"
        print(row)

    print(f"\n  H19a: intra-tier η > cross-tier η  p={mw_p:.4f}  "
          f"{'✓' if mw_p < 0.05 else '✗'}")
    print(f"  H19b: |ΔCV| stronger than |ΔAR1|   rs_cv={r_cv:.3f}, rs_ar1={r_ar1:.3f}  "
          f"{'✓' if h19b else '✗'}")
    print(f"  H19c: Pleurotus→Schizophyllum highest  {'✓' if h19c else f'✗ (actual: {top_src}→{top_tgt})'}")

    # ── 7. Save report ────────────────────────────────────────────────────
    report = {
        "species":      species,
        "cv_tiers":     {"bursty": sorted(TIER_BURSTY), "smooth": sorted(TIER_SMOOTH)},
        "native_S":     {sp: round(native_S[sp], 6) for sp in species},
        "eta_matrix":   {species[i]: {species[j]: (round(float(eta_matrix[i, j]), 6)
                                                    if not np.isnan(eta_matrix[i, j]) else None)
                                      for j in range(n)}
                         for i in range(n)},
        "top5_pairs":   [{"source": s, "target": t, "eta": round(e, 6)}
                         for s, t, e in pair_etas[:5]],
        "hypotheses": {
            "H19a": {"result": "confirmed" if mw_p < 0.05 else "not_confirmed",
                     "mw_p": round(float(mw_p), 6) if not np.isnan(mw_p) else None,
                     "intra_median": round(float(np.median(intra_eta)), 4),
                     "cross_median": round(float(np.median(cross_eta)), 4)},
            "H19b": {"result": "confirmed" if h19b else "not_confirmed",
                     "rs_cv": round(float(r_cv), 4),
                     "rs_ar1": round(float(r_ar1), 4),
                     "p_cv": round(float(p_cv), 4),
                     "p_ar1": round(float(p_ar1), 4)},
            "H19c": {"result": "confirmed" if h19c else "not_confirmed",
                     "top_pair": f"{top_src}→{top_tgt}",
                     "top_eta": round(top_e, 4)},
        },
    }
    out = str(RESULTS_DIR / "19_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report → 19_report.json")
    print("Done.  Figures → experiments/results/19_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
