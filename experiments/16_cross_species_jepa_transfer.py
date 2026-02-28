"""
Experiment 16 — Cross-Species JEPA Transfer
============================================

A TinyJEPA predictor is trained once on Pleurotus ostreatus rho series
(the reference / best-characterised species, Exp 12) and then applied
WITHOUT any retraining to all four multispecies targets from Exp 14/15.

DESIGN
------
1. TRAIN:     TinyJEPA fitted on Pleurotus rho series from
              pleurotus_spike_aligned.ndjson (6507 × 120 s windows, 9 days).
              This encodes the Pleurotus "temporal grammar".

2. TRANSFER:  Trained model is frozen and scored on each target species's
              windowed combined density (ρ_local + ρ_global) WITHOUT
              retraining.  High transfer surprise = target dynamics are
              alien to the Pleurotus grammar.

3. NATIVE:    For comparison, the same TinyJEPA architecture is freshly
              trained on each target species's own data (identical hyper-
              parameters).  Native surprise is the achievable floor.

4. TRANSFER EFFICIENCY  η = native_S / transfer_S  ∈ [0, 1]
   η → 1.0  : target dynamics are fully captured by the Pleurotus model.
   η → 0.0  : target dynamics are completely alien (model surprise >> floor).

HYPOTHESES
   H16a:  Transfer surprise increases monotonically with ecological stress
          distance |Δ_stress| from Pleurotus (stress_Pleurotus = 0.35).
          Ecology shapes the temporal grammar more than taxonomy.
   H16b:  Transfer efficiency η correlates negatively with evolutionary
          distance from Pleurotus (Mya).  Close relatives → similar grammar.
   H16c:  At least one species achieves η > 0.70 (meaningful shared
          temporal structure exists across the fungal kingdom).

UKFT FRAMING
   Pleurotus represents a "reference attractor" in temporal density space.
   The JEPA transfer score measures how far each species's knowledge-density
   trajectory departs from this attractor — i.e., the entropic distance
   between species-specific choice-collapse patterns.

Data
   • Pleurotus:    tools/data/pleurotus_spike_aligned.ndjson (Exp 12 source)
   • Multispecies: tools/data/multispecies/*.txt (Exp 14/15 source)

Outputs (experiments/results/)
   16_pleurotus_training.png          Training curve + Pleurotus rho series
   16_transfer_vs_native.png          Native vs transfer surprise per species
   16_transfer_efficiency_vs_eco.png  η vs Δ_stress + evo. distance
   16_temporal_structure_matrix.png   Pairwise temporal stats heatmap
   16_surprise_ratio_vs_stress.png    Transfer/native ratio vs stress rank
   16_report.json                     All metrics
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import uniform_filter1d
from scipy.stats import spearmanr

# ── Paths ────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
ALIGNED_NDJ = MYCO / "tools/data/pleurotus_spike_aligned.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (same as Exp 14/15) ────────────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600      # 10-min windows for multispecies

# ── JEPA hyper-parameters (same as Exp 15) ───────────────────────────────────
SEQ_LEN         = 16
PRED_HORIZON    = 4
HIDDEN_DIM      = 32
LEARNING_RATE   = 3e-3
N_EPOCHS        = 60
MIN_WINDOWS     = 64
DEGENERATE_STD  = 1e-4

# ── Colour palette ────────────────────────────────────────────────────────────
SPECIES_COLOURS = {
    "Pleurotus":     "#2ecc71",
    "Schizophyllum": "#9b59b6",
    "Cordyceps":     "#e67e22",
    "Enoki":         "#27ae60",
    "Omphalotus":    "#2980b9",
}

# ─────────────────────────────────────────────────────────────────────────────
# SEMANTIC METADATA  (carried from Exp 15; Pleurotus added as reference)
# Evolutionary distances are FROM Pleurotus ostreatus.
# ─────────────────────────────────────────────────────────────────────────────
PLEUROTUS_META = {
    "full_name":          "Pleurotus ostreatus",
    "phylum":             "Basidiomycota",
    "order":              "Agaricales",
    "ecological_stress":  0.35,   # mild temperate saprophyte; opportunistic fruiter
    "colony_lifespan_yr": 5.0,
    "mating_types":       100,    # tetrapolar; moderate diversity
    "evolutionary_age_Ma":0.0,    # reference — all distances measured from here
    "broadcast_cost":     2.0,    # moderate metabolic rate; dense mycelium
}

SPECIES_META = {
    "Schizophyllum": {
        "full_name":          "Schizophyllum commune",
        "phylum":             "Basidiomycota",
        "order":              "Polyporales",
        "ecological_stress":  0.90,
        "colony_lifespan_yr": 30.0,
        "mating_types":       28000,
        "evolutionary_age_Ma":250.0,
        "broadcast_cost":     4.5,
        "local_global_ratio": 69.578,
    },
    "Cordyceps": {
        "full_name":          "Cordyceps militaris",
        "phylum":             "Ascomycota",
        "order":              "Hypocreales",
        "ecological_stress":  0.70,
        "colony_lifespan_yr": 0.2,
        "mating_types":       2,
        "evolutionary_age_Ma":800.0,
        "broadcast_cost":     3.0,
        "local_global_ratio": 16.618,
    },
    "Enoki": {
        "full_name":          "Flammulina velutipes",
        "phylum":             "Basidiomycota",
        "order":              "Agaricales",
        "ecological_stress":  0.50,
        "colony_lifespan_yr": 7.0,
        "mating_types":       150,
        "evolutionary_age_Ma":100.0,
        "broadcast_cost":     2.5,
        "local_global_ratio": 0.247,
    },
    "Omphalotus": {
        "full_name":          "Omphalotus nidiformis",
        "phylum":             "Basidiomycota",
        "order":              "Agaricales",
        "ecological_stress":  0.20,
        "colony_lifespan_yr": 2.0,
        "mating_types":       4,
        "evolutionary_age_Ma":150.0,
        "broadcast_cost":     1.5,
        "local_global_ratio": 0.003,
    },
}

SPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
    "Enoki":         "Enoki fungi Flammulina velutipes.txt",
    "Omphalotus":    "Ghost Fungi Omphalotus nidiformis.txt",
}


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_pleurotus_rho() -> np.ndarray:
    """
    Load Pleurotus ostreatus rho time series from aligned NDJSON.
    Records are already sorted by t_start_s (120 s apart, 6507 windows).
    Returns: shape (6507,) float32 rho array.
    """
    if not ALIGNED_NDJ.exists():
        raise FileNotFoundError(f"Pleurotus data not found: {ALIGNED_NDJ}")
    records = [json.loads(l) for l in open(ALIGNED_NDJ)]
    records.sort(key=lambda r: r["t_start_s"])
    rhos = np.array([r["rho"] for r in records], dtype=np.float32)
    print(f"  Pleurotus rho: {len(rhos)} windows, "
          f"min={rhos.min():.3f} max={rhos.max():.3f} "
          f"mean={rhos.mean():.3f} std={rhos.std():.3f}", flush=True)
    return rhos


def load_combined_density(name: str) -> Optional[np.ndarray]:
    """
    Load multispecies .txt, compute windowed combined density
    (ρ_local + ρ_global per window).  Returns 1D float32 array or None.
    """
    path = DATA_DIR / SPECIES_FILES[name]
    if not path.exists():
        print(f"  [SKIP] {name}: file not found", flush=True)
        return None
    print(f"  Loading {name} …", flush=True)
    df = pd.read_csv(str(path), sep="\t", header=0, dtype=np.float32,
                     na_values=["NaN", "", " "], engine="c",
                     on_bad_lines="skip")
    df.dropna(how="all", inplace=True)
    df = df.select_dtypes(include=[np.floating, np.integer])
    arr = df.values.astype(np.float32)
    n, nchan = arr.shape

    # Detrend
    det = np.empty_like(arr)
    for ch in range(nchan):
        col = arr[:, ch].astype(np.float64)
        nan_mask = np.isnan(col)
        if nan_mask.any():
            idxs = np.arange(n)
            col[nan_mask] = np.interp(idxs[nan_mask], idxs[~nan_mask], col[~nan_mask])
        trend = uniform_filter1d(col, size=600, mode="mirror")
        det[:, ch] = (col - trend).astype(np.float32)

    # Tier masks
    abs_v   = np.abs(det)
    sigma   = np.nanmedian(abs_v, axis=0) / 0.6745
    floor   = (NOISE_FLOOR_K * sigma).astype(np.float32)
    noise_m = abs_v < floor[np.newaxis, :]
    global_m = abs_v >= GLOBAL_THRESH_MV
    local_m  = (~noise_m) & (~global_m)

    # Windowed combined density
    n_win = n // WINDOW_S
    trim_l = local_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    trim_g = global_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    combined = (trim_l + trim_g).astype(np.float32).mean(axis=(1, 2))
    print(f"    → {n:,} rows × {nchan} ch  |  {n_win} windows, "
          f"combined density mean={combined.mean():.4f} std={combined.std():.4f}",
          flush=True)
    return combined


# ─────────────────────────────────────────────────────────────────────────────
# TinyJEPA  (same architecture as Exp 15)
# ─────────────────────────────────────────────────────────────────────────────

def relu(x: np.ndarray) -> np.ndarray:
    return np.maximum(0.0, x)


class TinyJEPA:
    """Minimal JEPA: linear encoder + predictor + output head, SGD training."""

    def __init__(self, seq_len: int = SEQ_LEN, hidden: int = HIDDEN_DIM,
                 horizon: int = PRED_HORIZON, lr: float = LEARNING_RATE):
        self.seq_len = seq_len
        self.hidden  = hidden
        self.horizon = horizon
        self.lr      = lr
        # Encoder  (seq_len,) → (hidden,)
        s1 = np.sqrt(2.0 / seq_len)
        self.We = np.random.randn(hidden, seq_len).astype(np.float32) * s1
        self.be = np.zeros(hidden, dtype=np.float32)
        # Predictor  (hidden,) → (hidden,)
        s2 = np.sqrt(2.0 / hidden)
        self.Wp = np.random.randn(hidden, hidden).astype(np.float32) * s2
        self.bp = np.zeros(hidden, dtype=np.float32)
        # Output  (hidden,) → scalar
        self.Wo = np.random.randn(1, hidden).astype(np.float32) * s2
        self.bo = np.zeros(1, dtype=np.float32)

    def encode(self, x: np.ndarray) -> np.ndarray:
        return relu(self.We @ x + self.be)

    def predict_latent(self, h: np.ndarray) -> np.ndarray:
        return relu(self.Wp @ h + self.bp)

    def decode(self, h: np.ndarray) -> float:
        return float((self.Wo @ h + self.bo)[0])

    def _forward(self, ctx: np.ndarray, target: float) -> tuple[float, float]:
        """Returns (prediction, squared error) — no in-place mutation."""
        h = self.predict_latent(self.encode(ctx))
        pred = self.decode(h)
        return pred, (pred - target) ** 2

    def _backward_step(self, ctx: np.ndarray, target: float) -> float:
        h_enc  = self.encode(ctx)
        h_pred = self.predict_latent(h_enc)
        pred   = self.decode(h_pred)
        err    = pred - target

        dWo = err * h_pred[np.newaxis, :]
        dbo = np.array([err])

        delta_pred = (self.Wo.T * err).flatten() * (h_pred > 0)
        dWp = np.outer(delta_pred, h_enc)
        dbp = delta_pred

        delta_enc = (self.Wp.T @ delta_pred) * (h_enc > 0)
        dWe = np.outer(delta_enc, ctx)
        dbe = delta_enc

        reg = 1e-4
        self.Wo -= self.lr * (dWo + reg * self.Wo)
        self.bo -= self.lr * dbo
        self.Wp -= self.lr * (dWp + reg * self.Wp)
        self.bp -= self.lr * dbp
        self.We -= self.lr * (dWe + reg * self.We)
        self.be -= self.lr * dbe
        return err ** 2

    def _normalise(self, series: np.ndarray) -> np.ndarray:
        return (series - series.mean()) / (series.std() + 1e-8)

    def train(self, series: np.ndarray, epochs: int = N_EPOCHS) -> list[float]:
        s = self._normalise(series).astype(np.float32)
        n = len(s)
        losses = []
        for _ in range(epochs):
            idx = list(range(n - self.seq_len - self.horizon))
            np.random.shuffle(idx)
            ep_loss = sum(
                self._backward_step(s[i:i + self.seq_len],
                                    float(s[i + self.seq_len + self.horizon - 1]))
                for i in idx
            )
            losses.append(ep_loss / max(len(idx), 1))
        return losses

    def surprise_scores(self, series: np.ndarray) -> np.ndarray:
        """Score a series using CURRENT (possibly pretrained) weights."""
        s = self._normalise(series).astype(np.float32)
        n = len(s)
        return np.array(
            [self._forward(s[i:i + self.seq_len],
                           float(s[i + self.seq_len + self.horizon - 1]))[1]
             for i in range(n - self.seq_len - self.horizon)],
            dtype=np.float32,
        )


# ─────────────────────────────────────────────────────────────────────────────
# JEPA UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def _degenerate_result(label: str, reason: str) -> dict:
    return {"label": label, "mean_surprise": float("nan"),
            "median_surprise": float("nan"), "p90_surprise": float("nan"),
            "degenerate": True, "degenerate_reason": reason}


def run_jepa_native(series: np.ndarray, label: str) -> tuple[dict, Optional["TinyJEPA"]]:
    """Train a fresh TinyJEPA on series.  Returns (result_dict, trained_model)."""
    if len(series) < MIN_WINDOWS:
        return _degenerate_result(label, "too_few_windows"), None
    if np.std(series) < DEGENERATE_STD:
        return _degenerate_result(label, f"constant_series_std={np.std(series):.2e}"), None
    np.random.seed(42)
    model = TinyJEPA()
    losses = model.train(series)
    if np.isnan(losses[-1]):
        return _degenerate_result(label, "nan_training_loss"), None
    scores = model.surprise_scores(series)
    print(f"    {label}: train_loss={losses[-1]:.4f}  "
          f"S_mean={scores.mean():.4f}  S_p90={np.percentile(scores, 90):.4f}",
          flush=True)
    return {
        "label": label,
        "mean_surprise":   float(scores.mean()),
        "median_surprise": float(np.median(scores)),
        "p90_surprise":    float(np.percentile(scores, 90)),
        "final_train_loss": float(losses[-1]),
        "n_windows":       len(series),
        "degenerate":      False,
    }, model


def run_jepa_transfer(pretrained: "TinyJEPA", series: np.ndarray, label: str) -> dict:
    """Score series with a PRETRAINED (frozen) model — no weight updates."""
    if len(series) < MIN_WINDOWS:
        return _degenerate_result(label, "too_few_windows")
    if np.std(series) < DEGENERATE_STD:
        return _degenerate_result(label, f"constant_series_std={np.std(series):.2e}")
    scores = pretrained.surprise_scores(series)
    if np.isnan(scores).all():
        return _degenerate_result(label, "nan_transfer_scores")
    print(f"    {label} (transfer): "
          f"S_mean={scores.mean():.4f}  S_p90={np.percentile(scores, 90):.4f}",
          flush=True)
    return {
        "label": label,
        "mean_surprise":   float(scores.mean()),
        "median_surprise": float(np.median(scores)),
        "p90_surprise":    float(np.percentile(scores, 90)),
        "n_windows":       len(series),
        "degenerate":      False,
    }


# ─────────────────────────────────────────────────────────────────────────────
# TEMPORAL STRUCTURE STATISTICS
# used to build the pairwise similarity matrix
# ─────────────────────────────────────────────────────────────────────────────

def temporal_stats(series: np.ndarray) -> dict:
    """
    Compute a small feature vector describing the temporal structure of a
    1D density series, after z-score normalisation.
    Features: AR1 (lag-1 autocorrelation), CV (coeff. of variation),
              skewness, kurtosis, spectral entropy.
    """
    s = series.astype(np.float64)
    mu, sig = s.mean(), s.std() + 1e-12
    z = (s - mu) / sig

    # AR1
    ar1 = float(np.corrcoef(z[:-1], z[1:])[0, 1]) if len(z) > 2 else 0.0

    # CV on raw series
    cv = float(sig / (abs(mu) + 1e-12))

    # Skewness
    skew = float(np.mean(z ** 3))

    # Excess kurtosis
    kurt = float(np.mean(z ** 4) - 3.0)

    # Spectral entropy (normalised power spectrum)
    ps = np.abs(np.fft.rfft(z)) ** 2
    ps /= (ps.sum() + 1e-12)
    spec_ent = float(-np.sum(ps * np.log(ps + 1e-12)) / np.log(len(ps) + 1))

    return {"ar1": ar1, "cv": cv, "skew": skew, "kurt": kurt, "spec_ent": spec_ent}


def temporal_distance(a: dict, b: dict) -> float:
    """Euclidean distance in normalised temporal-stat feature space."""
    keys = ["ar1", "cv", "skew", "kurt", "spec_ent"]
    va = np.array([a[k] for k in keys])
    vb = np.array([b[k] for k in keys])
    return float(np.linalg.norm(va - vb))


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_pleurotus_training(rho: np.ndarray, losses: list[float],
                           scores: np.ndarray) -> None:
    """Training curve + Pleurotus rho time series + self-surprise histogram."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    # Panel 1: rho time series (first 500 windows for readability)
    ax = axes[0]
    disp = rho[:500]
    ax.plot(disp, color=SPECIES_COLOURS["Pleurotus"], linewidth=0.7, alpha=0.8)
    ax.set_xlabel("Window index (120 s each)", fontsize=10)
    ax.set_ylabel("Pleurotus rho (spikes / window)", fontsize=10)
    ax.set_title("Pleurotus ostreatus — rho time series\n(first 500 of 6507 windows)",
                 fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 2: training loss curve
    ax = axes[1]
    ax.semilogy(losses, color=SPECIES_COLOURS["Pleurotus"], linewidth=1.5)
    ax.set_xlabel("Epoch", fontsize=10)
    ax.set_ylabel("Mean MSE loss (log scale)", fontsize=10)
    ax.set_title(f"JEPA training loss\nfinal = {losses[-1]:.5f}", fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 3: self-surprise distribution
    ax = axes[2]
    ax.hist(scores, bins=50, color=SPECIES_COLOURS["Pleurotus"],
            alpha=0.8, edgecolor="white")
    ax.axvline(scores.mean(), color="black", linestyle="--",
               linewidth=1.5, label=f"mean = {scores.mean():.4f}")
    ax.set_xlabel("JEPA surprise (MSE)", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title("Pleurotus self-surprise distribution\n(training set score)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Exp 16 — Pleurotus Reference JEPA (training baseline)",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    out = str(RESULTS_DIR / "16_pleurotus_training.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 16_pleurotus_training.png", flush=True)


def fig_transfer_vs_native(results: dict) -> None:
    """
    Side-by-side grouped bars: native vs transfer surprise per species.
    A gap between bars signals how much the Pleurotus model over-predicts
    the target species's surprise.
    """
    names    = list(results.keys())
    native_s = [results[n]["native"]["mean_surprise"]   for n in names]
    xfer_s   = [results[n]["transfer"]["mean_surprise"] for n in names]
    native_d = [results[n]["native"].get("degenerate", False)   for n in names]
    xfer_d   = [results[n]["transfer"].get("degenerate", False) for n in names]

    def safe(v): return 0.0 if (v is None or (isinstance(v, float) and np.isnan(v))) else v

    x     = np.arange(len(names))
    w     = 0.35
    fig, ax = plt.subplots(figsize=(11, 5))

    bars_n = ax.bar(x - w / 2, [safe(v) for v in native_s], w,
                    label="Native (species-specific training)",
                    color=[SPECIES_COLOURS.get(n, "#888") for n in names],
                    edgecolor="white", alpha=0.90)
    bars_t = ax.bar(x + w / 2, [safe(v) for v in xfer_s], w,
                    label="Transfer (Pleurotus model, no retraining)",
                    color=[SPECIES_COLOURS.get(n, "#888") for n in names],
                    edgecolor="black", linestyle="--", alpha=0.55,
                    hatch=["" if not d else "//" for d in xfer_d])

    # Annotate bars
    for bar, val, deg in zip(list(bars_n) + list(bars_t),
                              native_s + xfer_s,
                              native_d + xfer_d):
        txt = "degen." if deg else (f"{val:.3f}" if not np.isnan(val) else "nan")
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.001, txt,
                ha="center", va="bottom", fontsize=7,
                color="gray" if deg else "black")

    # Efficiency label above pairs
    for i, n in enumerate(names):
        nv, tv = native_s[i], xfer_s[i]
        if not (np.isnan(nv) or np.isnan(tv) or tv == 0):
            eta = nv / tv
            ax.text(i, max(safe(nv), safe(tv)) + 0.015,
                    f"η={eta:.2f}", ha="center", fontsize=8,
                    color=SPECIES_COLOURS.get(n, "black"), fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel("Mean JEPA surprise S", fontsize=11)
    ax.set_title(
        "Native vs Transfer JEPA Surprise per Species\n"
        "η = native / transfer  (higher η = better transfer from Pleurotus grammar)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "16_transfer_vs_native.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 16_transfer_vs_native.png", flush=True)


def fig_transfer_efficiency_vs_eco(results: dict, meta: dict) -> None:
    """
    Two-panel scatter:
    Left:  η vs |Δ_stress| from Pleurotus (H16a)
    Right: η vs evolutionary distance Mya   (H16b)
    """
    names = [n for n in results if not results[n]["transfer"].get("degenerate")
             and not results[n]["native"].get("degenerate")
             and not np.isnan(results[n]["native"]["mean_surprise"])
             and not np.isnan(results[n]["transfer"]["mean_surprise"])
             and results[n]["transfer"]["mean_surprise"] > 0]

    if not names:
        print("  [SKIP] fig_transfer_efficiency_vs_eco: no valid species", flush=True)
        return

    etas   = [results[n]["native"]["mean_surprise"] / results[n]["transfer"]["mean_surprise"]
              for n in names]
    dstress = [abs(meta[n]["ecological_stress"] - PLEUROTUS_META["ecological_stress"])
               for n in names]
    evod    = [meta[n]["evolutionary_age_Ma"] for n in names]
    colours = [SPECIES_COLOURS.get(n, "#888") for n in names]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, xvals, xlabel, htitle in [
        (axes[0], dstress,
         "|Δ ecological stress| from Pleurotus",
         "H16a: larger stress distance → lower η"),
        (axes[1], evod,
         "Evolutionary distance from Pleurotus (Mya)",
         "H16b: older divergence → lower η"),
    ]:
        for x, y, c, n in zip(xvals, etas, colours, names):
            ax.scatter(x, y, s=140, color=c, zorder=3)
            ax.text(x + (max(xvals) - min(xvals)) * 0.02, y + 0.005,
                    n, fontsize=9, color=c)
        if len(names) >= 2:
            rs, pval = spearmanr(xvals, etas)
            m, b = np.polyfit(xvals, etas, 1)
            xfit = np.linspace(min(xvals), max(xvals), 50)
            ax.plot(xfit, m * xfit + b, "k--", linewidth=1, alpha=0.6)
        else:
            rs, pval = float("nan"), float("nan")
        ax.axhline(0.70, color="steelblue", linestyle=":", linewidth=1.2,
                   alpha=0.8, label="η = 0.70 (H16c threshold)")
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Transfer efficiency η = S_native / S_transfer", fontsize=10)
        ax.set_title(f"{htitle}\nSpearman rs = {rs:.3f}  p = {pval:.3f}", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Exp 16 — Transfer Efficiency vs Ecological & Phylogenetic Distance",
                 fontsize=12)
    plt.tight_layout()
    out = str(RESULTS_DIR / "16_transfer_efficiency_vs_eco.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 16_transfer_efficiency_vs_eco.png", flush=True)


def fig_temporal_structure_matrix(ts_stats: dict) -> None:
    """
    Heatmap of pairwise temporal-stat distances.
    All 5 species including Pleurotus; lighter = more similar temporal grammar.
    """
    species_order = ["Pleurotus"] + list(SPECIES_META.keys())
    present = [s for s in species_order if s in ts_stats]
    n = len(present)
    if n < 2:
        print("  [SKIP] temporal structure matrix: too few species", flush=True)
        return

    dist_mat = np.zeros((n, n), dtype=np.float32)
    for i, a in enumerate(present):
        for j, b in enumerate(present):
            dist_mat[i, j] = temporal_distance(ts_stats[a], ts_stats[b])

    # Normalise to [0, 1]
    maxd = dist_mat.max() + 1e-9
    norm_mat = dist_mat / maxd

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(norm_mat, cmap="RdYlGn_r", vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Normalised temporal distance (0=identical)")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(present, rotation=30, ha="right", fontsize=10)
    ax.set_yticklabels(present, fontsize=10)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{norm_mat[i, j]:.2f}", ha="center", va="center",
                    fontsize=9, color="white" if norm_mat[i, j] > 0.6 else "black")
    ax.set_title(
        "Pairwise Temporal Structure Distance\n"
        "(AR1, CV, skewness, kurtosis, spectral entropy)",
        fontsize=11,
    )
    plt.tight_layout()
    out = str(RESULTS_DIR / "16_temporal_structure_matrix.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 16_temporal_structure_matrix.png", flush=True)


def fig_surprise_ratio_vs_stress(results: dict, meta: dict) -> None:
    """
    Transfer surprise / native surprise ratio (over-surprise factor)
    vs ecological stress rank.  Species ranked by stress; Pleurotus as anchor.
    """
    all_entries = {"Pleurotus": {"stress": PLEUROTUS_META["ecological_stress"]}}
    all_entries.update({n: {"stress": meta[n]["ecological_stress"]} for n in meta})

    plot_data = []
    for n in results:
        nv = results[n]["native"]["mean_surprise"]
        tv = results[n]["transfer"]["mean_surprise"]
        if not (np.isnan(nv) or np.isnan(tv) or nv == 0 or tv == 0
                or results[n]["native"].get("degenerate")
                or results[n]["transfer"].get("degenerate")):
            ratio = tv / nv   # >1 means transfer is worse; =1 means perfect
            stress = meta[n]["ecological_stress"]
            plot_data.append((stress, ratio, n))

    if not plot_data:
        print("  [SKIP] fig_surprise_ratio_vs_stress: no valid data", flush=True)
        return

    plot_data.sort(key=lambda t: t[0])
    stress_vals = [t[0] for t in plot_data]
    ratios      = [t[1] for t in plot_data]
    names_      = [t[2] for t in plot_data]
    colours_    = [SPECIES_COLOURS.get(n, "#888") for n in names_]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(stress_vals, ratios, s=150,
               color=colours_, zorder=3)
    for x, y, n in zip(stress_vals, ratios, names_):
        ax.text(x + 0.01, y + 0.02, n, fontsize=9,
                color=SPECIES_COLOURS.get(n, "black"))

    if len(plot_data) >= 2:
        rs, pval = spearmanr(stress_vals, ratios)
        m, b = np.polyfit(stress_vals, ratios, 1)
        xfit = np.linspace(min(stress_vals) - 0.05, max(stress_vals) + 0.05, 50)
        ax.plot(xfit, m * xfit + b, "k--", linewidth=1, alpha=0.6)
        rs_str = f"Spearman rs = {rs:.3f}  p = {pval:.3f}"
    else:
        rs_str = "n < 3 — Spearman not computed"

    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.9,
               label="ratio = 1 (perfect transfer)")
    ax.axvline(PLEUROTUS_META["ecological_stress"], color=SPECIES_COLOURS["Pleurotus"],
               linestyle="--", linewidth=1.0, alpha=0.7,
               label=f"Pleurotus stress = {PLEUROTUS_META['ecological_stress']}")
    ax.set_xlabel("Ecological stress (0 = benign, 1 = harsh)", fontsize=11)
    ax.set_ylabel("Transfer / native surprise ratio\n(1.0 = perfect transfer)", fontsize=11)
    ax.set_title(
        "Transfer Over-Surprise vs Ecological Stress\n" + rs_str,
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "16_surprise_ratio_vs_stress.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 16_surprise_ratio_vs_stress.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 16 — Cross-Species JEPA Transfer")
    print("Reference species: Pleurotus ostreatus")
    print("=" * 70)

    # ── 1. Load and train on Pleurotus ────────────────────────────────────
    print("\n--- Step 1: Train JEPA on Pleurotus rho series ---", flush=True)
    pleu_rho = load_pleurotus_rho()

    np.random.seed(42)
    pleu_model = TinyJEPA()
    print(f"  Training for {N_EPOCHS} epochs …", flush=True)
    pleu_losses = pleu_model.train(pleu_rho, epochs=N_EPOCHS)
    if np.isnan(pleu_losses[-1]):
        raise RuntimeError("Pleurotus JEPA training diverged — check data")
    print(f"  Final loss: {pleu_losses[-1]:.5f}", flush=True)

    pleu_scores = pleu_model.surprise_scores(pleu_rho)
    print(f"  Pleurotus self-surprise: "
          f"mean={pleu_scores.mean():.4f}  p90={np.percentile(pleu_scores, 90):.4f}",
          flush=True)

    # Temporal stats for Pleurotus
    ts_stats = {"Pleurotus": temporal_stats(pleu_rho)}

    # Figure 1: training summary
    fig_pleurotus_training(pleu_rho, pleu_losses, pleu_scores)

    # ── 2. Transfer + native scoring for each target species ─────────────
    print("\n--- Step 2: Transfer + native scoring ---", flush=True)
    results = {}

    for name in SPECIES_META:
        print(f"\n[{name}]", flush=True)
        combined = load_combined_density(name)
        if combined is None:
            results[name] = {
                "native":   _degenerate_result(name, "file_not_found"),
                "transfer": _degenerate_result(name, "file_not_found"),
            }
            continue

        ts_stats[name] = temporal_stats(combined)

        # Transfer: apply Pleurotus model (frozen)
        res_xfer = run_jepa_transfer(pleu_model, combined, name)

        # Native: train fresh model on this species
        res_native, _ = run_jepa_native(combined, name)

        # Transfer efficiency
        nv = res_native["mean_surprise"]
        tv = res_xfer["mean_surprise"]
        if not (np.isnan(nv) or np.isnan(tv) or tv == 0):
            eta = nv / tv
            h16c = "H16c ✓" if eta > 0.70 else "H16c ✗"
        else:
            eta, h16c = float("nan"), "—"
        print(f"  η = {eta:.3f}  {h16c}" if not np.isnan(eta) else "  η = nan",
              flush=True)

        results[name] = {
            "native":   res_native,
            "transfer": res_xfer,
            "efficiency_eta": float(eta) if not np.isnan(eta) else None,
        }

    # ── 3. Figures ────────────────────────────────────────────────────────
    print("\n--- Step 3: Figures ---", flush=True)
    fig_transfer_vs_native(results)
    fig_transfer_efficiency_vs_eco(results, SPECIES_META)
    fig_temporal_structure_matrix(ts_stats)
    fig_surprise_ratio_vs_stress(results, SPECIES_META)

    # ── 4. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 16 SUMMARY")
    print("=" * 70)
    print(f"\n  Pleurotus self-surprise (training baseline):")
    print(f"    mean = {pleu_scores.mean():.4f}   p90 = {np.percentile(pleu_scores, 90):.4f}")

    print(f"\n  {'Species':20s}  {'Native S':>9}  {'Transfer S':>11}  {'η':>7}  H16c")
    print("  " + "-" * 58)
    h16c_count = 0
    valid_count = 0
    for n in SPECIES_META:
        nv = results[n]["native"]["mean_surprise"]
        tv = results[n]["transfer"]["mean_surprise"]
        eta = results[n].get("efficiency_eta")
        if eta is not None:
            valid_count += 1
            h16c_check = "✓" if eta > 0.70 else "✗"
            if eta > 0.70:
                h16c_count += 1
            print(f"  {n:20s}  {nv:9.4f}  {tv:11.4f}  {eta:7.3f}  {h16c_check}")
        else:
            nv_str = f"{nv:.4f}" if not np.isnan(nv) else "  NaN  "
            tv_str = f"{tv:.4f}" if not np.isnan(tv) else "  NaN  "
            print(f"  {n:20s}  {nv_str:>9}  {tv_str:>11}  {'—':>7}  degenerate")

    print(f"\n  H16c confirmed in {h16c_count}/{valid_count} valid species "
          f"(η > 0.70)")

    # Ecological correlation
    valid_eco = [
        (abs(SPECIES_META[n]["ecological_stress"] - PLEUROTUS_META["ecological_stress"]),
         results[n]["transfer"]["mean_surprise"])
        for n in SPECIES_META
        if results[n].get("efficiency_eta") is not None
    ]
    if len(valid_eco) >= 3:
        ds_vals = [v[0] for v in valid_eco]
        ts_vals = [v[1] for v in valid_eco]
        rs_eco, p_eco = spearmanr(ds_vals, ts_vals)
        print(f"\n  H16a (transfer S ↑ with |Δ stress|): rs = {rs_eco:.3f}  p = {p_eco:.3f}")
    else:
        print(f"\n  H16a: n = {len(valid_eco)} — insufficient for Spearman")

    valid_evo = [
        (SPECIES_META[n]["evolutionary_age_Ma"],
         results[n].get("efficiency_eta"))
        for n in SPECIES_META
        if results[n].get("efficiency_eta") is not None
    ]
    if len(valid_evo) >= 3:
        ev_vals  = [v[0] for v in valid_evo]
        eta_vals = [v[1] for v in valid_evo]
        rs_evo, p_evo = spearmanr(ev_vals, eta_vals)
        print(f"  H16b (η ↓ with evo. distance):       rs = {rs_evo:.3f}  p = {p_evo:.3f}")
    else:
        print(f"  H16b: n = {len(valid_evo)} — insufficient for Spearman")

    print("=" * 70)

    # ── 5. Save report ────────────────────────────────────────────────────
    report = {
        "pleurotus_self_surprise": {
            "mean": float(pleu_scores.mean()),
            "p90":  float(np.percentile(pleu_scores, 90)),
        },
        "results": {
            n: {
                "native_S_mean":   results[n]["native"]["mean_surprise"],
                "transfer_S_mean": results[n]["transfer"]["mean_surprise"],
                "efficiency_eta":  results[n].get("efficiency_eta"),
                "native_degenerate":   results[n]["native"].get("degenerate"),
                "transfer_degenerate": results[n]["transfer"].get("degenerate"),
            }
            for n in SPECIES_META
        },
        "temporal_stats": {
            n: ts_stats[n] for n in ts_stats
        },
    }
    out = str(RESULTS_DIR / "16_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report saved → 16_report.json")
    print("Done.  Figures → experiments/results/16_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
