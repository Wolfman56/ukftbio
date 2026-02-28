"""
Experiment 21 — Normalisation Ablation: Killing the Kurtosis Advantage
=======================================================================

Exp 20 found that |Δkurtosis| is the only statistically significant predictor
of JEPA transfer efficiency η (rs = +0.480, p = 0.032).  The proposed mechanism:

  Omphalotus density (κ = 73.703) is massively sparse-spike — nearly all windows
  near zero with rare extremes.  TinyJEPA trains a "quiet predictor" biased toward
  zero.  Cordyceps (κ = −1.149, near-Gaussian) happens to be predicted well by this
  zero-attractor, because its normalised values also cluster near zero.

TinyJEPA already z-scores each series globally before training (in _norm()).
Z-scoring preserves kurtosis: the distribution SHAPE is scale-invariant.
To truly kill the kurtosis advantage, we need a transformation that equalises
distribution *shape*, not just mean and variance.

Exp 21 applies **rank normalisation** before JEPA training:
  1. Convert each rho window series to uniform quantile ranks ∈ [0, 1]
  2. Apply inverse-normal transform (probit) → Gaussian marginal for all species
  3. Re-run the full N×N JEPA transfer matrix
  4. Compare η matrix to Exp 19 (raw z-score) to measure the kurtosis contribution

If the kurtosis mechanism is the driver, rank normalisation should:
  - Collapse the Omphalotus→Cordyceps η advantage (ideally below 0.70)
  - Reduce the overall η range (max − min shrinks)
  - Narrow the gap between Omphalotus-as-source (best) and Enoki-as-source (worst)
  - Potentially improve Schizophyllum as source (de-trending its H > 1 structure)

HYPOTHESES
   H21a:  Omphalotus→Cordyceps η drops after rank-normalisation (≤ 0.70),
          confirming that kurtosis contrast mechanically created the η=0.790 result.
   H21b:  The η range (max − min over 20 pairs) shrinks by ≥ 30% vs Exp 19.
   H21c:  The Spearman correlation between Exp 19 η and Exp 21 η across all 20
          pairs is ≥ 0.75 — the rank *ordering* of pairs is preserved even if
          the absolute values shift.
   H21d:  Schizophyllum becomes a better source after normalisation (row mean η
          increases relative to Exp 19), consistent with its superdiffusive H=1.216
          trend being the transfer barrier.

Figures
-------
  21_eta_comparison.png     Side-by-side heatmap: Exp 19 (raw) vs Exp 21 (rank-norm)
  21_delta_eta.png          Δη = Exp21 − Exp19 per pair (who gains, who loses)
  21_source_row_means.png   Per-species row-mean η: Exp 19 vs Exp 21 (grouped bar)
  21_range_collapse.png     Distribution of η: Exp 19 vs Exp 21 (overlapping KDE)
  21_report.json            Full result + hypothesis test results
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
from scipy.special import ndtri  # probit / inverse-normal

# ── Paths ─────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
PLEU_EMB    = MYCO / "tools/data/pleurotus_spike_emb_40d.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (consistent with Exps 14–20) ────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600

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

# ── JEPA parameters (unchanged from Exps 15–20) ──────────────────────────────
SEQ_LEN      = 16
PRED_HORIZON = 4
HIDDEN_DIM   = 32
LR           = 3e-3
N_EPOCHS     = 60
MIN_WINDOWS  = 30
DEGENERATE   = 1e-4

# ── Exp 19 baseline η matrix (raw z-score normalisation) ─────────────────────
EXP19_SPECIES = ["Pleurotus", "Schizophyllum", "Cordyceps", "Enoki", "Omphalotus"]
EXP19_ETA = np.array([
    [np.nan, 0.2866, 0.5495, 0.3368, 0.5076],
    [0.5920, np.nan, 0.3500, 0.3555, 0.4895],
    [0.6397, 0.4992, np.nan, 0.3737, 0.7263],
    [0.5405, 0.2379, 0.3700, np.nan, 0.4320],
    [0.6320, 0.5286, 0.7904, 0.3912, np.nan],
])


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADERS
# ─────────────────────────────────────────────────────────────────────────────

def load_pleurotus_rho() -> np.ndarray:
    records = [json.loads(l) for l in open(PLEU_EMB)]
    records.sort(key=lambda r: r["t_start_s"])
    rho = np.array([r["rho"] for r in records], dtype=np.float32)
    print(f"  Pleurotus: {len(rho)} windows  rho mean={rho.mean():.4f}")
    return rho


def load_multispecies_rho(species: str) -> np.ndarray:
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
# RANK NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────

def rank_normalise(x: np.ndarray) -> np.ndarray:
    """
    Rank-normalise a 1D array:
      1. Convert to fractional ranks ∈ (0, 1) with van der Waerden correction
         (avoid exact 0 and 1 at boundaries).
      2. Apply probit (inverse-normal CDF) to map to Gaussian marginal.

    After transformation every species has:
      mean ≈ 0, std ≈ 1, skewness ≈ 0, kurtosis ≈ 0
    regardless of the original distribution shape.

    This is the minimal transformation that equalises distribution *shape*
    while preserving monotone rank ordering of the original values.
    """
    n = len(x)
    ranks = stats.rankdata(x, method="average")          # 1..n (average for ties)
    quantiles = ranks / (n + 1.0)                         # van der Waerden: ∈ (0,1)
    quantiles = np.clip(quantiles, 1e-9, 1 - 1e-9)        # numerical safety
    return ndtri(quantiles).astype(np.float32)             # probit → N(0,1)


# ─────────────────────────────────────────────────────────────────────────────
# TINY JEPA (identical weights + structure to Exps 15–20)
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

    def clone(self) -> "TinyJEPA":
        other = TinyJEPA.__new__(TinyJEPA)
        for attr in ("We", "be", "Wp", "bp", "Wo", "bo"):
            setattr(other, attr, getattr(self, attr).copy())
        return other

    # NOTE: _norm is a pass-through here — rank_normalise is applied at load time
    def _norm(self, s): return s.astype(np.float64)

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

    def train(self, series: np.ndarray) -> None:
        s = self._norm(series)
        idx = list(range(len(s) - SEQ_LEN - PRED_HORIZON))
        for _ in range(N_EPOCHS):
            np.random.shuffle(idx)
            for i in idx:
                self._step(s[i:i+SEQ_LEN], float(s[i+SEQ_LEN+PRED_HORIZON-1]))

    def mean_surprise(self, series: np.ndarray) -> float:
        s = self._norm(series)
        n = len(s)
        if n < SEQ_LEN + PRED_HORIZON + 1:
            return float("nan")
        errs = [(self._dec(self._pred(self._enc(s[i:i+SEQ_LEN])))
                 - float(s[i+SEQ_LEN+PRED_HORIZON-1])) ** 2
                for i in range(n - SEQ_LEN - PRED_HORIZON)]
        return float(np.mean(errs))


def train_jepa(rho: np.ndarray, label: str) -> tuple[TinyJEPA | None, float]:
    if len(rho) < MIN_WINDOWS or rho.std() < DEGENERATE:
        print(f"    {label}: degenerate → skip")
        return None, float("nan")
    model = TinyJEPA()
    model.train(rho)
    native_s = model.mean_surprise(rho)
    print(f"    {label}: n={len(rho)}  native_S={native_s:.4f}")
    return model, native_s


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_eta_comparison(species: list[str], eta19: np.ndarray, eta21: np.ndarray) -> None:
    n = len(species)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, eta, title in [
        (axes[0], eta19, "Exp 19 — raw z-score"),
        (axes[1], eta21, "Exp 21 — rank-normalised"),
    ]:
        masked = np.ma.masked_where(np.eye(n, dtype=bool), eta)
        im = ax.imshow(masked, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto")
        plt.colorbar(im, ax=ax, label="η", fraction=0.04)
        for i in range(n):
            for j in range(n):
                if i != j and not np.isnan(eta[i, j]):
                    val = eta[i, j]
                    c = "black" if 0.25 < val < 0.85 else "white"
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=9, color=c, fontweight="bold")
                elif i == j:
                    ax.text(j, i, "—", ha="center", va="center", fontsize=10, color="gray")

        labels = [f"{sp}\n({'B' if sp in {'Pleurotus','Schizophyllum'} else 'S'})"
                  for sp in species]
        ax.set_xticks(range(n)); ax.set_xticklabels(labels, fontsize=9)
        ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Target"); ax.set_ylabel("Source")
        ax.set_title(title, fontsize=11)

    plt.suptitle("JEPA Transfer Matrix: Raw vs Rank-Normalised\n"
                 "H21a: Does Omphalotus→Cordyceps drop below 0.70 after rank-normalisation?",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    out = str(RESULTS_DIR / "21_eta_comparison.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 21_eta_comparison.png")


def fig_delta_eta(species: list[str], eta19: np.ndarray, eta21: np.ndarray) -> None:
    n = len(species)
    delta = eta21 - eta19   # positive = rank-norm helped; negative = hurt

    fig, ax = plt.subplots(figsize=(8, 7))
    cmap = plt.cm.RdBu
    masked = np.ma.masked_where(np.eye(n, dtype=bool), delta)
    im = ax.imshow(masked, cmap=cmap, vmin=-0.35, vmax=0.35, aspect="auto")
    plt.colorbar(im, ax=ax, label="Δη = Exp21 − Exp19", fraction=0.04)

    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(delta[i, j]):
                val = delta[i, j]
                c = "black" if abs(val) < 0.25 else "white"
                sign = "+" if val >= 0 else ""
                ax.text(j, i, f"{sign}{val:.3f}", ha="center", va="center",
                        fontsize=9, color=c, fontweight="bold")
            elif i == j:
                ax.text(j, i, "—", ha="center", va="center", fontsize=10, color="gray")

    labels = [f"{sp[:8]}\n({'B' if sp in {'Pleurotus','Schizophyllum'} else 'S'})"
              for sp in species]
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Target", fontsize=11)
    ax.set_ylabel("Source", fontsize=11)
    ax.set_title("Δη = Rank-Normalised − Raw\nRed = rank-norm hurt  |  Blue = rank-norm helped",
                 fontsize=11)
    plt.tight_layout()
    out = str(RESULTS_DIR / "21_delta_eta.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 21_delta_eta.png")


def fig_source_row_means(species: list[str], eta19: np.ndarray, eta21: np.ndarray) -> None:
    n = len(species)
    means19 = [np.nanmean(eta19[i, :]) for i in range(n)]
    means21 = [np.nanmean(eta21[i, :]) for i in range(n)]

    x = np.arange(n); w = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    bars19 = ax.bar(x - w / 2, means19, w, color=[SPECIES_COLOURS[sp] for sp in species],
                    alpha=0.6, edgecolor="black", linewidth=0.7, label="Exp 19 (raw)")
    bars21 = ax.bar(x + w / 2, means21, w, color=[SPECIES_COLOURS[sp] for sp in species],
                    alpha=1.0, edgecolor="black", linewidth=0.7, label="Exp 21 (rank-norm)",
                    hatch="///")

    # Annotate delta
    for xi, (m19, m21) in enumerate(zip(means19, means21)):
        delta = m21 - m19
        sign  = "+" if delta >= 0 else ""
        ax.text(xi, max(m19, m21) + 0.015, f"{sign}{delta:.3f}",
                ha="center", fontsize=8, color="black")

    ax.set_xticks(x); ax.set_xticklabels(species, fontsize=9, rotation=15)
    ax.set_ylabel("Row-mean η (as source)", fontsize=11)
    ax.set_title("Per-Species η Row Mean: Raw vs Rank-Normalised\n"
                 "H21d: Does Schizophyllum row mean increase?  "
                 "H21a: Does Omphalotus lose its top rank?",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(0, 0.9)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "21_source_row_means.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 21_source_row_means.png")


SPECIES_COLOURS = {
    "Pleurotus":     "#2ecc71",
    "Schizophyllum": "#9b59b6",
    "Cordyceps":     "#e67e22",
    "Enoki":         "#27ae60",
    "Omphalotus":    "#2980b9",
}


def fig_range_collapse(eta19: np.ndarray, eta21: np.ndarray,
                       species: list[str]) -> None:
    n = len(species)
    vals19 = [float(eta19[i, j]) for i in range(n) for j in range(n)
              if i != j and not np.isnan(eta19[i, j])]
    vals21 = [float(eta21[i, j]) for i in range(n) for j in range(n)
              if i != j and not np.isnan(eta21[i, j])]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: overlapping histogram
    ax = axes[0]
    ax.hist(vals19, bins=10, alpha=0.5, color="#e74c3c", label=f"Exp 19 raw  range={max(vals19)-min(vals19):.3f}")
    ax.hist(vals21, bins=10, alpha=0.5, color="#3498db", label=f"Exp 21 rank  range={max(vals21)-min(vals21):.3f}")
    ax.axvline(np.median(vals19), color="#e74c3c", linestyle="--", linewidth=1.5)
    ax.axvline(np.median(vals21), color="#3498db", linestyle="--", linewidth=1.5)
    ax.set_xlabel("η", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("η Distribution: Raw vs Rank-Normalised", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: Exp 19 vs Exp 21 scatter (H21c: rank ordering preserved?)
    ax2 = axes[1]
    pair_colours = [SPECIES_COLOURS[species[i]]
                    for i in range(n) for j in range(n)
                    if i != j and not np.isnan(eta19[i, j])]
    ax2.scatter(vals19, vals21, c=pair_colours, s=60, alpha=0.85,
                edgecolor="black", linewidth=0.5, zorder=5)
    ax2.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.4, label="y=x (no change)")
    rs, p = stats.spearmanr(vals19, vals21)
    ax2.set_xlabel("η — Exp 19 (raw)", fontsize=11)
    ax2.set_ylabel("η — Exp 21 (rank-norm)", fontsize=11)
    ax2.set_title(f"Rank Ordering Preservation\n"
                  f"Spearman rs = {rs:.3f}  p = {p:.4f}  "
                  f"H21c: {'✓ rs≥0.75' if rs >= 0.75 else '✗ rs<0.75'}",
                  fontsize=11)
    patches = [mpatches.Patch(color=SPECIES_COLOURS[sp], label=sp)
               for sp in species]
    ax2.legend(handles=patches, fontsize=8, loc="upper left")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = str(RESULTS_DIR / "21_range_collapse.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 21_range_collapse.png")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 21 — Normalisation Ablation: Killing the Kurtosis Advantage")
    print("=" * 70)
    print("\nExp 19 η range: [0.238, 0.790]  Omphalotus κ=73.703  Cordyceps κ=-1.149")
    print("Test: does rank-normalisation (probit transform) destroy the advantage?")

    # ── 1. Load density series ────────────────────────────────────────────
    print("\n--- Step 1: Load density series ---")
    rho_raw: dict[str, np.ndarray] = {}
    rho_raw["Pleurotus"] = load_pleurotus_rho()
    for sp in MULTISPECIES_FILES:
        rho_raw[sp] = load_multispecies_rho(sp)
    species = list(rho_raw.keys())

    # ── 2. Rank-normalise each series ────────────────────────────────────
    print("\n--- Step 2: Rank-normalise (probit transform) ---")
    rho_ranked: dict[str, np.ndarray] = {}
    for sp, raw in rho_raw.items():
        ranked = rank_normalise(raw)
        rho_ranked[sp] = ranked
        orig_kurt = float(stats.kurtosis(raw))
        new_kurt  = float(stats.kurtosis(ranked))
        print(f"  {sp:15s}  κ_before={orig_kurt:8.3f} → κ_after={new_kurt:6.3f}  "
              f"n={len(raw)}")

    # ── 3. Train TinyJEPA on ranked series ───────────────────────────────
    print("\n--- Step 3: Train TinyJEPA on rank-normalised series ---")
    models: dict[str, TinyJEPA | None] = {}
    native_S: dict[str, float] = {}
    for sp in species:
        print(f"  [{sp}]")
        m, ns = train_jepa(rho_ranked[sp], sp)
        models[sp] = m; native_S[sp] = ns

    # ── 4. Transfer matrix (rank-norm) ───────────────────────────────────
    print("\n--- Step 4: Transfer matrix (rank-normalised) ---")
    n = len(species)
    transfer_S = np.full((n, n), np.nan)
    eta21      = np.full((n, n), np.nan)

    for i, src in enumerate(species):
        if models[src] is None:
            continue
        for j, tgt in enumerate(species):
            if i == j:
                continue
            frozen = models[src].clone()
            ts = frozen.mean_surprise(rho_ranked[tgt])
            transfer_S[i, j] = ts
            ns = native_S[tgt]
            if not np.isnan(ns) and ts > 0:
                eta21[i, j] = ns / ts
            tag = f"η={eta21[i,j]:.4f}" if not np.isnan(eta21[i,j]) else "η=N/A"
            print(f"  {src:15s} → {tgt:15s}  transfer_S={ts:.4f}  {tag}")

    # ── 5. Hypothesis tests ──────────────────────────────────────────────
    print("\n--- Step 5: Hypothesis tests ---")

    eta19 = EXP19_ETA

    # H21a: Omphalotus→Cordyceps drops to ≤ 0.70
    i_om = species.index("Omphalotus"); j_co = species.index("Cordyceps")
    i_co = species.index("Cordyceps");  j_om = species.index("Omphalotus")
    om_co_19 = eta19[i_om, j_co]; om_co_21 = eta21[i_om, j_co]
    co_om_19 = eta19[i_co, j_om]; co_om_21 = eta21[i_co, j_om]
    h21a = om_co_21 <= 0.70
    print(f"\n  H21a — Omphalotus→Cordyceps drops to ≤ 0.70:")
    print(f"    Exp 19 (raw):        η = {om_co_19:.4f}")
    print(f"    Exp 21 (rank-norm):  η = {om_co_21:.4f}")
    print(f"    Cordyceps→Omphalotus:  {co_om_19:.4f} → {co_om_21:.4f}")
    print(f"    → {'✓ CONFIRMED' if h21a else '✗ NOT CONFIRMED'}")

    # H21b: η range shrinks ≥ 30%
    vals19 = [float(eta19[i,j]) for i in range(n) for j in range(n)
              if i!=j and not np.isnan(eta19[i,j])]
    vals21 = [float(eta21[i,j]) for i in range(n) for j in range(n)
              if i!=j and not np.isnan(eta21[i,j])]
    range19 = max(vals19) - min(vals19)
    range21 = max(vals21) - min(vals21)
    shrink_pct = (range19 - range21) / range19 * 100
    h21b = shrink_pct >= 30
    print(f"\n  H21b — η range shrinks ≥ 30%:")
    print(f"    Exp 19 range: [{min(vals19):.4f}, {max(vals19):.4f}] = {range19:.4f}")
    print(f"    Exp 21 range: [{min(vals21):.4f}, {max(vals21):.4f}] = {range21:.4f}")
    print(f"    Shrinkage: {shrink_pct:.1f}%")
    print(f"    → {'✓ CONFIRMED (≥ 30%)' if h21b else '✗ NOT CONFIRMED'}")

    # H21c: rank ordering preserved (Spearman rs ≥ 0.75 between Exp19 η and Exp21 η)
    common_19 = []; common_21 = []
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if not np.isnan(eta19[i,j]) and not np.isnan(eta21[i,j]):
                common_19.append(float(eta19[i,j]))
                common_21.append(float(eta21[i,j]))
    rs_order, p_order = stats.spearmanr(common_19, common_21)
    h21c = rs_order >= 0.75
    print(f"\n  H21c — Rank ordering preserved (rs ≥ 0.75):")
    print(f"    Spearman rs = {rs_order:.4f}  p = {p_order:.4f}")
    print(f"    → {'✓ CONFIRMED' if h21c else '✗ NOT CONFIRMED'}")

    # H21d: Schizophyllum row mean increases
    sz_idx = species.index("Schizophyllum")
    sz_mean19 = np.nanmean(eta19[sz_idx, :])
    sz_mean21 = np.nanmean(eta21[sz_idx, :])
    h21d = sz_mean21 > sz_mean19
    print(f"\n  H21d — Schizophyllum row mean η increases:")
    print(f"    Exp 19: {sz_mean19:.4f}  Exp 21: {sz_mean21:.4f}"
          f"  Δ = {sz_mean21 - sz_mean19:+.4f}")
    print(f"    → {'✓ CONFIRMED' if h21d else '✗ NOT CONFIRMED'}")

    # ── 6. Full comparison table ──────────────────────────────────────────
    print("\n--- Full η comparison (Exp 19 raw → Exp 21 rank-norm) ---")
    header = f"  {'Pair':30s}  {'η_19':>7}  {'η_21':>7}  {'Δη':>7}"
    print(header)
    pair_data = []
    for i, si in enumerate(species):
        for j, sj in enumerate(species):
            if i == j: continue
            e19 = float(eta19[i, j]) if not np.isnan(eta19[i, j]) else None
            e21 = float(eta21[i, j]) if not np.isnan(eta21[i, j]) else None
            if e19 is not None and e21 is not None:
                d = e21 - e19
                flag = "★" if (si == "Omphalotus" and sj == "Cordyceps") else \
                       ("★" if (si == "Cordyceps" and sj == "Omphalotus") else "")
                print(f"  {si[:14]:14s} → {sj[:14]:14s}  "
                      f"{e19:7.4f}  {e21:7.4f}  {d:+7.4f}  {flag}")
                pair_data.append({"src": si, "tgt": sj, "eta19": round(e19, 5),
                                   "eta21": round(e21, 5), "delta": round(d, 5)})

    # ── 7. Figures ────────────────────────────────────────────────────────
    print("\n--- Step 6: Figures ---")
    fig_eta_comparison(species, eta19, eta21)
    fig_delta_eta(species, eta19, eta21)
    fig_source_row_means(species, eta19, eta21)
    fig_range_collapse(eta19, eta21, species)

    # ── 8. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 21 SUMMARY")
    print("=" * 70)
    print(f"\n  Omphalotus→Cordyceps:  η_19={om_co_19:.4f} → η_21={om_co_21:.4f}"
          f"  ({'kurtosis advantage REMOVED' if h21a else 'kurtosis advantage PERSISTS'})")
    print(f"  η range:               {range19:.4f} → {range21:.4f}"
          f"  ({shrink_pct:+.1f}% change)")
    print(f"  Rank ordering rs:      {rs_order:.4f}  p={p_order:.4f}")
    print(f"  Schizophyllum row Δη:  {sz_mean21 - sz_mean19:+.4f}")

    print(f"\n  H21a: Omphalotus→Cordyceps ≤ 0.70  {'✓' if h21a else '✗'}")
    print(f"  H21b: Range shrinks ≥ 30%         {'✓' if h21b else '✗'}")
    print(f"  H21c: Ranking preserved rs ≥ 0.75 {'✓' if h21c else '✗'}")
    print(f"  H21d: Schizophyllum improves       {'✓' if h21d else '✗'}")

    # ── 9. Save report ────────────────────────────────────────────────────
    report = {
        "exp19_eta_range": round(range19, 5),
        "exp21_eta_range": round(range21, 5),
        "range_shrinkage_pct": round(float(shrink_pct), 2),
        "ordering_rs": round(float(rs_order), 5),
        "ordering_p":  round(float(p_order), 5),
        "key_pairs": {
            "Omphalotus→Cordyceps": {"eta19": round(om_co_19, 5), "eta21": round(float(om_co_21), 5)},
            "Cordyceps→Omphalotus": {"eta19": round(co_om_19, 5), "eta21": round(float(co_om_21), 5)},
        },
        "schizophyllum_row_mean": {"eta19": round(float(sz_mean19), 5),
                                    "eta21": round(float(sz_mean21), 5)},
        "all_pairs": pair_data,
        "hypotheses": {
            "H21a": {"result": "confirmed" if h21a else "not_confirmed",
                     "om_co_eta21": round(float(om_co_21), 4)},
            "H21b": {"result": "confirmed" if h21b else "not_confirmed",
                     "range_shrinkage_pct": round(float(shrink_pct), 2)},
            "H21c": {"result": "confirmed" if h21c else "not_confirmed",
                     "ordering_rs": round(float(rs_order), 4)},
            "H21d": {"result": "confirmed" if h21d else "not_confirmed",
                     "sz_delta": round(float(sz_mean21 - sz_mean19), 4)},
        },
    }
    out = str(RESULTS_DIR / "21_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report → 21_report.json")
    print("Done.  Figures → experiments/results/21_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
