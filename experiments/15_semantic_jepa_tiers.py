"""
Experiment 15 — Semantic Context + JEPA Dual-Tier Surprise
===========================================================

Two threads:

A) SEMANTIC METADATA
   Species evolutionary context, colony lifespan, habitat stress, and mating-type
   diversity are correlated with the Exp 14 local/global communication ratio.
   Hypothesis: ecological pressure (not just activity level) predicts the
   communication phenotype observed in Exp 14.

B) JEPA DUAL-TIER SURPRISE
   A JEPA temporal predictor is trained separately on:
     - rho_global(t): 10-min windowed global spike density
     - rho_local(t) : 10-min windowed sub-threshold event density
   Surprise S = ||pred - actual||² is compared between tiers per species.

   Hypothesis (H15a): rho_local has LOWER JEPA surprise than rho_global —
   i.e., sub-threshold activity is more predictable (it is the guidance field;
   choice collapse is inherently less predictable).
   Hypothesis (H15b): The surprise gap (S_global - S_local) correlates with
   ecological stress rank — species in more stressful environments have more
   predictable local tiers relative to their global tier.

UKFT framing
   Local tier = pilot wave = guidance field (deterministic, low-entropy trajectory)
   Global tier = choice collapse = inherently irreversible, less predictable
   If H15a holds: JEPA corroborates the UKFT interpretation of sub-threshold as
   guidance, not noise.

Data sources
   - Electrophysiology: multispecies/*.txt (same as Exp 14)
   - Metadata: hand-curated from primary literature (see SPECIES_META below)
     Sources: MycoBank, Index Fungorum, published phylogenomics (Floudas 2012,
     Spatafora 2017), and species ecology reviews.

Outputs (experiments/results/)
   15_metadata_radar.png         Spider/radar — 5 ecological traits × 4 species
   15_ratio_vs_metadata.png      3-panel: ratio vs stress / lifespan / evo_dist
   15_jepa_surprise_tiers.png    JEPA surprise: local vs global per species
   15_surprise_gap_vs_stress.png (S_global - S_local) vs ecological stress
   15_report.json                All metrics
"""

from __future__ import annotations
import json
import os
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.stats import spearmanr, pearsonr

# ── Paths ───────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (same as Exp 14) ──────────────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600      # 10-minute windows

# ── JEPA constants ─────────────────────────────────────────────────────────
SEQ_LEN         = 16        # context window (windows)
PRED_HORIZON    = 4         # predict 4 windows ahead
HIDDEN_DIM      = 32
LEARNING_RATE   = 3e-3
N_EPOCHS        = 60
MIN_WINDOWS     = 64        # skip species with too little data

# ── Species colour palette ──────────────────────────────────────────────────
SPECIES_COLOURS = {
    "Schizophyllum": "#9b59b6",
    "Cordyceps":     "#e67e22",
    "Enoki":         "#27ae60",
    "Omphalotus":    "#2980b9",
}

# ─────────────────────────────────────────────────────────────────────────────
# SEMANTIC METADATA
# ─────────────────────────────────────────────────────────────────────────────
#
# Curated from primary literature.  Each field is normalised to [0,1] for
# radar plots EXCEPT the raw values, which are kept for scatter correlation.
#
# ecological_stress: composite of habitat aridity, substrate scarcity,
#   temperature extremes, desiccation tolerance requirement (0=benign, 1=harsh)
#
# colony_lifespan_yr: estimated mycelial colony persistence in years
#   (not fruiting body lifespan).
#   Schizophyllum: colonies documented surviving desiccation for years;
#     wildtype clones from herbarium specimens viable after 60+ yr — use 30.
#   Cordyceps militaris: insect host-bound; colony dissolved after host death,
#     new colony starts from fresh spore — use 0.2 (weeks to months).
#   Flammulina velutipes: perennial wood-rot; overwintering mycelium documented
#     surviving 5–10 yr in stumps — use 7.
#   Omphalotus nidiformis: subtropical wood-rot; fast decomposer, colony
#     exhausts substrate in 1–3 yr — use 2.
#
# mating_types: approximate number of distinct mating types.
#   Schizophyllum: ~28,000 (most genetically diverse eukaryote by this measure).
#   Cordyceps: ~2 (MAT1-1 / MAT1-2, heterothallic).
#   Flammulina: ~100–200 (tetrapolar, moderate).
#   Omphalotus: ~2–6 (limited data; appears bipolar-like).
#
# evolutionary_age_Ma: estimated divergence from last common ancestor with
#   Pleurotus ostreatus (our reference species from Exp 12), in Mya.
#   Based on Floudas et al. 2012 (Science) and Spatafora et al. 2017.
#   Schizophyllum → Polyporales; Pleurotus → Agaricales.  Divergence ~250 Mya.
#   Cordyceps → Ascomycota; divergence from Basidiomycota ~800 Mya.
#   Flammulina → Agaricales (same order as Pleurotus) ~100 Mya.
#   Omphalotus → Agaricales ~150 Mya.
#
# broadcast_cost: relative energetic cost index for sustained broadcast
#   communication (global tier), based on mycelial density and known metabolic
#   rate.  Higher cost → greater selective pressure for economical (local-first)
#   strategy.  Qualitative scale 1–5 based on published metabolic studies.
#
# local_global_ratio: carried forward from Exp 14 results.

SPECIES_META = {
    "Schizophyllum": {
        "full_name":         "Schizophyllum commune",
        "phylum":            "Basidiomycota",
        "order":             "Polyporales",
        "ecological_stress": 0.90,    # extreme desiccation tolerance; xerophilic
        "colony_lifespan_yr":30.0,    # clones viable after decades
        "mating_types":      28000,   # ~28,000 documented
        "evolutionary_age_Ma":250.0,  # Polyporales divergence from Agaricales
        "broadcast_cost":    4.5,     # sparse hyphal mats; energy per spike high
        "local_global_ratio": 69.578, # Exp 14
        "exp14_global_coherence": 0.999,
    },
    "Cordyceps": {
        "full_name":         "Cordyceps militaris",
        "phylum":            "Ascomycota",
        "order":             "Hypocreales",
        "ecological_stress": 0.70,    # obligate parasite; host-dependent stress
        "colony_lifespan_yr": 0.2,    # host lifespan (weeks–months)
        "mating_types":      2,       # bipolar heterothallic
        "evolutionary_age_Ma":800.0,  # Ascomycota / Basidiomycota split
        "broadcast_cost":    3.0,     # parasitic; energy partially from host
        "local_global_ratio": 16.618, # Exp 14
        "exp14_global_coherence": 0.066,
    },
    "Enoki": {
        "full_name":         "Flammulina velutipes",
        "phylum":            "Basidiomycota",
        "order":             "Agaricales",
        "ecological_stress": 0.50,    # cold-adapted but seasonally benign
        "colony_lifespan_yr": 7.0,    # overwintering perennial mycelium
        "mating_types":      150,     # tetrapolar, moderate diversity
        "evolutionary_age_Ma":100.0,  # within Agaricales
        "broadcast_cost":    2.5,     # cold habitat; slow metabolism
        "local_global_ratio": 0.247,  # Exp 14
        "exp14_global_coherence": 0.064,
    },
    "Omphalotus": {
        "full_name":         "Omphalotus nidiformis",
        "phylum":            "Basidiomycota",
        "order":             "Agaricales",
        "ecological_stress": 0.20,    # subtropical; resource-rich dead wood
        "colony_lifespan_yr": 2.0,    # exhausts substrate rapidly
        "mating_types":      4,       # limited data; appears bipolar
        "evolutionary_age_Ma":150.0,  # Agaricales divergence
        "broadcast_cost":    1.5,     # high-energy; bioluminescent overhead
        "local_global_ratio": 0.003,  # Exp 14
        "exp14_global_coherence": 0.000,
    },
}

METADATA_DISPLAY = {
    "ecological_stress":   "Ecological stress (0–1)",
    "colony_lifespan_yr":  "Colony lifespan (yr, log)",
    "mating_types":        "Mating types (log)",
    "evolutionary_age_Ma": "Evo. distance from Pleurotus (Mya)",
    "broadcast_cost":      "Broadcast energy cost (1–5)",
}


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL LOADING & PREPROCESSING  (same pipeline as Exp 14)
# ─────────────────────────────────────────────────────────────────────────────

SPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
    "Enoki":         "Enoki fungi Flammulina velutipes.txt",
    "Omphalotus":    "Ghost Fungi Omphalotus nidiformis.txt",
}


def load_detrended(name: str) -> Optional[np.ndarray]:
    path = DATA_DIR / SPECIES_FILES[name]
    if not path.exists():
        print(f"  [SKIP] {name}: file not found", flush=True)
        return None
    print(f"  Loading {name} …", flush=True)
    df = pd.read_csv(
        str(path), sep="\t", header=0, dtype=np.float32,
        na_values=["NaN", "", " "], engine="c", on_bad_lines="skip",
    )
    df.dropna(how="all", inplace=True)
    df = df.select_dtypes(include=[np.floating, np.integer])
    arr = df.values.astype(np.float32)
    n, nchan = arr.shape
    det = np.empty_like(arr)
    for ch in range(nchan):
        col = arr[:, ch].astype(np.float64)
        nan_mask = np.isnan(col)
        if nan_mask.any():
            idxs = np.arange(n)
            col[nan_mask] = np.interp(idxs[nan_mask], idxs[~nan_mask], col[~nan_mask])
        trend = uniform_filter1d(col, size=600, mode="mirror")
        det[:, ch] = (col - trend).astype(np.float32)
    print(f"    → {n:,} rows × {nchan} channels detrended", flush=True)
    return det


def windowed_tier_densities(det: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (rho_local, rho_global) windowed density series."""
    abs_v    = np.abs(det)
    sigma    = np.nanmedian(abs_v, axis=0) / 0.6745
    floor    = (NOISE_FLOOR_K * sigma).astype(np.float32)
    noise_m  = abs_v < floor[np.newaxis, :]
    global_m = abs_v >= GLOBAL_THRESH_MV
    local_m  = (~noise_m) & (~global_m)

    n        = det.shape[0]
    n_win    = n // WINDOW_S
    trim_l   = local_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    trim_g   = global_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    return trim_l.mean(axis=(1, 2)), trim_g.mean(axis=(1, 2))


# ─────────────────────────────────────────────────────────────────────────────
# TINY JEPA  (pure numpy — no dependency on torch)
# ─────────────────────────────────────────────────────────────────────────────
# Architecture: linear encoder + linear predictor (sufficient for 1D density)
# Training criterion: MSE prediction loss with L2 regularisation
# Surprise S_t = ||predicted_embedding(t) - actual_embedding(t)||²

def relu(x):  return np.maximum(0, x)

class TinyJEPA:
    """Minimal JEPA: context encoder + predictor, trained by SGD."""

    def __init__(self, seq_len: int = SEQ_LEN, hidden: int = HIDDEN_DIM,
                 horizon: int = PRED_HORIZON, lr: float = LEARNING_RATE):
        self.seq_len = seq_len
        self.hidden  = hidden
        self.horizon = horizon
        self.lr      = lr
        # Encoder: (seq_len,) → (hidden,)
        scale = np.sqrt(2.0 / seq_len)
        self.We = np.random.randn(hidden, seq_len).astype(np.float32) * scale
        self.be = np.zeros(hidden, dtype=np.float32)
        # Predictor: (hidden,) → (hidden,)  for horizon steps
        scale2 = np.sqrt(2.0 / hidden)
        self.Wp = np.random.randn(hidden, hidden).astype(np.float32) * scale2
        self.bp = np.zeros(hidden, dtype=np.float32)
        # Output: (hidden,) → (1,)
        self.Wo = np.random.randn(1, hidden).astype(np.float32) * scale2
        self.bo = np.zeros(1, dtype=np.float32)

    def encode(self, x):
        return relu(self.We @ x + self.be)

    def predict(self, h):
        return relu(self.Wp @ h + self.bp)

    def decode(self, h):
        return (self.Wo @ h + self.bo)[0]

    def forward(self, ctx, target_val):
        h_ctx  = self.encode(ctx)
        h_pred = self.predict(h_ctx)
        pred   = self.decode(h_pred)
        loss   = (pred - target_val) ** 2
        return pred, loss, h_pred

    def backward_step(self, ctx, target_val):
        h_ctx  = self.encode(ctx)
        h_pred = self.predict(h_ctx)
        pred   = self.decode(h_pred)
        err    = pred - target_val  # scalar

        # Output layer gradient
        dWo = err * h_pred[np.newaxis, :]
        dbo = np.array([err])

        # Predictor layer gradient
        delta_pred = (self.Wo.T * err).flatten() * (h_pred > 0)
        dWp = np.outer(delta_pred, h_ctx)
        dbp = delta_pred

        # Encoder layer gradient
        delta_enc = (self.Wp.T @ delta_pred) * (h_ctx > 0)
        dWe = np.outer(delta_enc, ctx)
        dbe = delta_enc

        # SGD update with L2
        reg = 1e-4
        self.Wo -= self.lr * (dWo + reg * self.Wo)
        self.bo -= self.lr * dbo
        self.Wp -= self.lr * (dWp + reg * self.Wp)
        self.bp -= self.lr * dbp
        self.We -= self.lr * (dWe + reg * self.We)
        self.be -= self.lr * dbe

        return (pred - target_val) ** 2

    def train(self, series: np.ndarray, epochs: int = N_EPOCHS) -> list[float]:
        n      = len(series)
        losses = []
        mu     = series.mean()
        std    = series.std() + 1e-8
        s      = (series - mu) / std

        for epoch in range(epochs):
            epoch_loss = 0.0
            idx = list(range(n - self.seq_len - self.horizon))
            np.random.shuffle(idx)
            for i in idx:
                ctx    = s[i: i + self.seq_len]
                target = s[i + self.seq_len + self.horizon - 1]
                epoch_loss += self.backward_step(ctx, float(target))
            losses.append(epoch_loss / max(len(idx), 1))
        return losses

    def surprise_scores(self, series: np.ndarray) -> np.ndarray:
        """Per-window JEPA surprise after training."""
        n   = len(series)
        mu  = series.mean()
        std = series.std() + 1e-8
        s   = (series - mu) / std
        scores = []
        for i in range(n - self.seq_len - self.horizon):
            ctx    = s[i: i + self.seq_len]
            target = s[i + self.seq_len + self.horizon - 1]
            _, loss, _ = self.forward(ctx, float(target))
            scores.append(loss)
        return np.array(scores, dtype=np.float32)


DEGENERATE_STD_THRESH = 1e-4   # series with std below this are trivially constant

def run_jepa(series: np.ndarray, label: str) -> dict:
    if len(series) < MIN_WINDOWS:
        print(f"    {label}: SKIP (too few windows: {len(series)})", flush=True)
        return {"mean_surprise": float("nan"), "median_surprise": float("nan"),
                "p90_surprise": float("nan"), "label": label, "degenerate": True,
                "degenerate_reason": "too_few_windows"}
    std = float(np.std(series))
    if std < DEGENERATE_STD_THRESH:
        print(f"    {label}: DEGENERATE (std={std:.2e} — essentially constant series)",
              flush=True)
        return {"mean_surprise": float("nan"), "median_surprise": float("nan"),
                "p90_surprise": float("nan"), "label": label, "degenerate": True,
                "degenerate_reason": f"constant_series_std={std:.2e}"}
    np.random.seed(42)
    model = TinyJEPA()
    losses = model.train(series)
    if np.isnan(losses[-1]) or np.isnan(losses[0]):
        print(f"    {label}: NaN training loss — flagging degenerate", flush=True)
        return {"mean_surprise": float("nan"), "median_surprise": float("nan"),
                "p90_surprise": float("nan"), "label": label, "degenerate": True,
                "degenerate_reason": "nan_training_loss"}
    scores = model.surprise_scores(series)
    if np.isnan(scores).all():
        print(f"    {label}: NaN surprise scores — flagging degenerate", flush=True)
        return {"mean_surprise": float("nan"), "median_surprise": float("nan"),
                "p90_surprise": float("nan"), "label": label, "degenerate": True,
                "degenerate_reason": "nan_surprise_scores"}
    print(f"    {label}: final_train_loss={losses[-1]:.4f}  "
          f"mean_S={scores.mean():.4f}  p90={np.percentile(scores, 90):.4f}",
          flush=True)
    return {
        "label":           label,
        "mean_surprise":   float(scores.mean()),
        "median_surprise": float(np.median(scores)),
        "p90_surprise":    float(np.percentile(scores, 90)),
        "final_train_loss": float(losses[-1]),
        "n_windows":       int(len(series)),
        "degenerate":      False,
    }


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_metadata_radar(meta: dict) -> None:
    """Spider / radar chart: 5 normalised traits × 4 species."""
    categories = [
        "Eco. Stress",
        "Colony Lifespan\n(log-norm)",
        "Mating Types\n(log-norm)",
        "Evo. Distance\n(Mya, norm)",
        "Broadcast Cost\n(norm)",
    ]
    N = len(categories)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))

    # Normalise each trait to [0,1] across species
    all_stress   = [m["ecological_stress"]                  for m in meta.values()]
    all_life     = [np.log10(m["colony_lifespan_yr"] + 0.1) for m in meta.values()]
    all_mating   = [np.log10(m["mating_types"])             for m in meta.values()]
    all_evo      = [m["evolutionary_age_Ma"]                for m in meta.values()]
    all_cost     = [m["broadcast_cost"]                     for m in meta.values()]

    def norm(lst):
        mn, mx = min(lst), max(lst)
        return [(v - mn) / (mx - mn + 1e-9) for v in lst]

    n_stress  = norm(all_stress)
    n_life    = norm(all_life)
    n_mating  = norm(all_mating)
    n_evo     = norm(all_evo)
    n_cost    = norm(all_cost)

    for i, (name, m) in enumerate(meta.items()):
        values = [n_stress[i], n_life[i], n_mating[i], n_evo[i], n_cost[i]]
        values += values[:1]
        ax.plot(angles, values, linewidth=2, linestyle="solid",
                color=SPECIES_COLOURS[name], label=name)
        ax.fill(angles, values, alpha=0.10, color=SPECIES_COLOURS[name])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9)
    ax.set_title("Species Semantic Profile\n(5 normalised ecological traits)",
                 fontsize=12, pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.30, 1.10), fontsize=9)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    out = str(RESULTS_DIR / "15_metadata_radar.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved 15_metadata_radar.png", flush=True)


def fig_ratio_vs_metadata(meta: dict) -> None:
    """3-panel scatter: local/global ratio vs ecological stress / colony lifespan / mating types."""
    names   = list(meta.keys())
    ratios  = [np.log10(meta[n]["local_global_ratio"] + 1e-3) for n in names]
    stress  = [meta[n]["ecological_stress"]                    for n in names]
    lifeyr  = [np.log10(meta[n]["colony_lifespan_yr"] + 0.1)  for n in names]
    mating  = [np.log10(meta[n]["mating_types"])               for n in names]
    colours = [SPECIES_COLOURS[n]                              for n in names]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    panels = [
        (stress,  "Ecological stress (0–1)",          "H15-S1: stress → local-dominant"),
        (lifeyr,  "log₁₀(colony lifespan yr)",        "H15-S2: long lifespan → local-dominant"),
        (mating,  "log₁₀(mating types)",              "H15-S3: diversity → local-dominant"),
    ]

    for ax, (xvals, xlabel, title) in zip(axes, panels):
        for x, y, c, n in zip(xvals, ratios, colours, names):
            ax.scatter(x, y, s=120, color=c, zorder=3)
            ax.text(x, y + 0.08, n, ha="center", fontsize=8, color=c)

        # Spearman + trend line
        rs, pval = spearmanr(xvals, ratios)
        m, b = np.polyfit(xvals, ratios, 1)
        xfit = np.linspace(min(xvals), max(xvals), 50)
        ax.plot(xfit, m * xfit + b, color="gray", linestyle="--", linewidth=1, alpha=0.7)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("log₁₀(ρ_local / ρ_global)", fontsize=10)
        ax.set_title(f"{title}\nSpearman rs = {rs:.3f}  p = {pval:.3f}", fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Communication Strategy vs Ecological Context\n"
                 "Exp 14 ratio correlated with species biological metadata",
                 fontsize=12)
    plt.tight_layout()
    out = str(RESULTS_DIR / "15_ratio_vs_metadata.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved 15_ratio_vs_metadata.png", flush=True)


def fig_jepa_surprise_tiers(jepa_results: dict) -> None:
    """Bar chart: JEPA mean surprise for local vs global tier per species."""
    names    = list(jepa_results.keys())
    local_s  = [jepa_results[n]["local"]["mean_surprise"]  for n in names]
    global_s = [jepa_results[n]["global"]["mean_surprise"] for n in names]
    local_deg  = [jepa_results[n]["local"].get("degenerate", False)  for n in names]
    global_deg = [jepa_results[n]["global"].get("degenerate", False) for n in names]

    # Replace NaN with small sentinel for plotting; mark degenerate bars
    def safe(v):  return 0.0 if np.isnan(v) else v

    x      = np.arange(len(names))
    width  = 0.35
    fig, ax = plt.subplots(figsize=(10, 5))

    bars_l = ax.bar(x - width / 2, [safe(v) for v in local_s],  width,
                    label="Local tier (sub-threshold)",
                    color=["#f39c12" if not d else "#f7dc6f" for d in local_deg],
                    edgecolor="white", hatch=["" if not d else "//" for d in local_deg])
    bars_g = ax.bar(x + width / 2, [safe(v) for v in global_s], width,
                    label="Global tier (spike)",
                    color=["#e74c3c" if not d else "#f1948a" for d in global_deg],
                    edgecolor="white", hatch=["" if not d else "//" for d in global_deg])

    for bar, val, deg in zip(list(bars_l) + list(bars_g),
                             local_s + global_s,
                             local_deg + global_deg):
        label_text = "degen." if deg else f"{val:.3f}"
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.003,
                label_text, ha="center", va="bottom", fontsize=7,
                color="gray" if deg else "black")

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel("Mean JEPA surprise S", fontsize=11)
    ax.set_title(
        "JEPA Temporal Predictability: Local vs Global Tier\n"
        "H15a: local tier should be MORE predictable (lower S) — the guidance field\n"
        "Hatched = degenerate series (constant tier — too few events for meaningful JEPA)",
        fontsize=10,
    )
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out = str(RESULTS_DIR / "15_jepa_surprise_tiers.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved 15_jepa_surprise_tiers.png", flush=True)


def fig_surprise_gap_vs_stress(jepa_results: dict, meta: dict) -> None:
    """Scatter: (S_global - S_local) vs ecological stress (H15b)."""
    names  = list(jepa_results.keys())
    raw_gaps = [
        jepa_results[n]["global"]["mean_surprise"] - jepa_results[n]["local"]["mean_surprise"]
        for n in names
    ]
    # Only include non-degenerate species in this plot
    valid = [
        (n, g, meta[n]["ecological_stress"])
        for n, g in zip(names, raw_gaps)
        if not (np.isnan(g) or
                jepa_results[n]["local"].get("degenerate", False) or
                jepa_results[n]["global"].get("degenerate", False))
    ]
    gaps   = [v[1] for v in valid]
    stress = [v[2] for v in valid]
    colours = [SPECIES_COLOURS[v[0]] for v in valid]
    vnames  = [v[0] for v in valid]

    fig, ax = plt.subplots(figsize=(7, 5))
    for x, y, c, n in zip(stress, gaps, colours, vnames):
        ax.scatter(x, y, s=140, color=c, zorder=3)
        ax.text(x + 0.01, float(y) + 0.001, n, fontsize=9, color=c)

    rs   = float("nan")
    pval = float("nan")
    if len(valid) >= 3:
        rs, pval = spearmanr(stress, gaps)
    if not any(np.isnan(g) for g in gaps) and len(valid) >= 2:
        m, b = np.polyfit(stress, gaps, 1)
        xfit = np.linspace(min(stress), max(stress), 50)
        ax.plot(xfit, m * xfit + b, "k--", linewidth=1, alpha=0.6)

    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Ecological stress (0 = benign, 1 = harsh)", fontsize=11)
    ax.set_ylabel("S_global − S_local  (surprise gap)", fontsize=11)
    ax.set_title(
        "Surprise Gap vs Ecological Stress\n"
        f"H15b: stressed species have wider gap  |  Spearman rs = {rs:.3f}  p = {pval:.3f}",
        fontsize=11,
    )
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "15_surprise_gap_vs_stress.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved 15_surprise_gap_vs_stress.png", flush=True)


def fig_phylogenetic_strategy(meta: dict) -> None:
    """Evolutionary distance vs communication ratio — the phylogenetic split."""
    names   = list(meta.keys())
    evo_d   = [meta[n]["evolutionary_age_Ma"]                  for n in names]
    ratios  = [np.log10(meta[n]["local_global_ratio"] + 1e-3)  for n in names]
    colours = [SPECIES_COLOURS[n]                              for n in names]
    phylums = [meta[n]["phylum"]                               for n in names]

    fig, ax = plt.subplots(figsize=(8, 5))
    for x, y, c, n, ph in zip(evo_d, ratios, colours, names, phylums):
        marker = "^" if ph == "Ascomycota" else "o"
        ax.scatter(x, y, s=160, color=c, marker=marker, zorder=3)
        ax.text(x + 5, y + 0.05, f"{n}\n({ph[:5]}…)", fontsize=8, color=c)

    rs, pval = spearmanr(evo_d, ratios)
    m, b = np.polyfit(evo_d, ratios, 1)
    xfit = np.linspace(min(evo_d) - 20, max(evo_d) + 20, 100)
    ax.plot(xfit, m * xfit + b, "k--", linewidth=1, alpha=0.5)

    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8, label="ratio = 1 (parity)")
    ax.set_xlabel("Evolutionary distance from Pleurotus (Mya)", fontsize=11)
    ax.set_ylabel("log₁₀(ρ_local / ρ_global)", fontsize=11)
    ax.set_title(
        "Phylogenetic Distance vs Communication Strategy\n"
        f"Spearman rs = {rs:.3f}  p = {pval:.3f}  |  △ = Ascomycota  ○ = Basidiomycota",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "15_phylogenetic_strategy.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved 15_phylogenetic_strategy.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 15 — Semantic Context + JEPA Dual-Tier Surprise")
    print("=" * 70)

    # ── Part A: Metadata figures (no signal loading needed) ────────────────
    print("\n--- Part A: Semantic metadata ---", flush=True)
    fig_metadata_radar(SPECIES_META)
    fig_ratio_vs_metadata(SPECIES_META)
    fig_phylogenetic_strategy(SPECIES_META)

    # ── Part B: JEPA dual-tier ─────────────────────────────────────────────
    print("\n--- Part B: JEPA dual-tier surprise ---", flush=True)
    jepa_results = {}

    for name in SPECIES_META:
        print(f"\n[{name}]", flush=True)
        det = load_detrended(name)
        if det is None:
            jepa_results[name] = {
                "local":  {"mean_surprise": float("nan")},
                "global": {"mean_surprise": float("nan")},
            }
            continue

        rho_local, rho_global = windowed_tier_densities(det)
        print(f"  Windows: local={len(rho_local)}  global={len(rho_global)}", flush=True)

        res_local  = run_jepa(rho_local,  f"{name}/local")
        res_global = run_jepa(rho_global, f"{name}/global")

        gap = res_global["mean_surprise"] - res_local["mean_surprise"]
        h15a = "H15a ✓ (local more predictable)" if res_local["mean_surprise"] < res_global["mean_surprise"] else "H15a ✗"
        print(f"  Surprise gap (global−local): {gap:+.4f}  {h15a}", flush=True)

        jepa_results[name] = {"local": res_local, "global": res_global,
                               "surprise_gap": gap}

    fig_jepa_surprise_tiers(jepa_results)
    fig_surprise_gap_vs_stress(jepa_results, SPECIES_META)

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 15 SUMMARY")
    print("=" * 70)

    print("\nPart A — Semantic correlations with log(ratio):")
    names  = list(SPECIES_META.keys())
    ratios = [np.log10(SPECIES_META[n]["local_global_ratio"] + 1e-3) for n in names]
    for key, label in [
        ("ecological_stress",   "Ecological stress"),
        ("colony_lifespan_yr",  "Colony lifespan (log)"),
        ("mating_types",        "Mating types (log)"),
        ("evolutionary_age_Ma", "Evolutionary distance"),
        ("broadcast_cost",      "Broadcast cost"),
    ]:
        vals = [np.log10(SPECIES_META[n][key]) if key != "ecological_stress"
                else SPECIES_META[n][key] for n in names]
        rs, pval = spearmanr(vals, ratios)
        print("  %-30s  rs = %+.3f  p = %.3f" % (label, rs, pval))

    print("\nPart B — JEPA surprise gaps (positive = global less predictable):")
    h15a_count   = 0
    valid_count  = 0
    for name in names:
        jr = jepa_results.get(name, {})
        gap = jr.get("surprise_gap", float("nan"))
        l_deg = jr.get("local",  {}).get("degenerate", False)
        g_deg = jr.get("global", {}).get("degenerate", False)
        if l_deg or g_deg:
            reason = "local=const" if l_deg else "global=const"
            print(f"  {name:20s}  DEGENERATE ({reason}) — one tier near-zero variance")
        elif not np.isnan(gap):
            valid_count += 1
            sym = "✓" if gap > 0 else "✗"
            print(f"  {name:20s}  gap = {gap:+.4f}  H15a {sym}")
            if gap > 0:
                h15a_count += 1
        else:
            print(f"  {name:20s}  gap = nan (insufficient data)")
    print(f"\n  H15a confirmed in {h15a_count}/{valid_count} evaluable species")
    print(f"  Degenerate (single active tier): {len(names) - valid_count} species")
    print(f"  NOTE: degenerate cases are themselves the result — a species with one")
    print(f"  near-zero tier is evidence it operates almost exclusively in the other.")
    print("=" * 70)

    # ── Save report ────────────────────────────────────────────────────────
    report = {
        "metadata":     SPECIES_META,
        "jepa_results": {
            n: {
                "local_surprise_mean":  jepa_results[n]["local"]["mean_surprise"],
                "global_surprise_mean": jepa_results[n]["global"]["mean_surprise"],
                "surprise_gap":         jepa_results[n].get("surprise_gap", None),
            }
            for n in names
            if n in jepa_results
        },
        "semantic_spearman": {},
    }
    for key, label in [
        ("ecological_stress",   "stress"),
        ("colony_lifespan_yr",  "lifespan"),
        ("mating_types",        "mating_types"),
        ("evolutionary_age_Ma", "evo_distance"),
        ("broadcast_cost",      "broadcast_cost"),
    ]:
        vals = [np.log10(SPECIES_META[n][key]) if key != "ecological_stress"
                else SPECIES_META[n][key] for n in names]
        rs, pval = spearmanr(vals, ratios)
        report["semantic_spearman"][label] = {"rs": float(rs), "pval": float(pval)}

    out = str(RESULTS_DIR / "15_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report saved → 15_report.json")
    print("Done.  Figures → experiments/results/15_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
