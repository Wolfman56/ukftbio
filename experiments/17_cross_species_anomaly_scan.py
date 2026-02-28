"""
Experiment 17 — Cross-Species Anomaly Scan (myco-explorer Borda fusion)
========================================================================

Runs the full myco-explorer UKFT anomaly scoring pipeline across all five
fungal species simultaneously (Pleurotus + 4 multispecies targets), then
performs a blind cross-species Borda scan to surface the kingdom-level
anomaly candidates.

PIPELINE (mirrors myco-explorer tools/ but self-contained)
-----------
1. EMBED   — Build 40D spike-histogram embeddings per species.
             Pleurotus: loaded from pleurotus_spike_emb_40d.ndjson (already run).
             Multispecies: inline spike detection + make_histogram from
                           tools/spike_embed.py + tools/export_spike_events.py.

2. S3      — Cosine distance of each window's embedding to the pooled
             CROSS-SPECIES centroid.  This is the kingdom-level anomaly flag:
             windows far from the shared background are candidates.

3. S1 (FMM) — Intra-species Fast Marching Method deviation score.
              kNN graph of 40D embeddings per species; deviation from
              rho-predicted edge weights (same algorithm as fmm_score.py).

4. S2 (JEPA) — TinyJEPA temporal surprise on species's windowed combined
               density series (rho_local + rho_global per 10-min window).
               Uses the same pure-numpy JEPA from Exp 15/16.

5. BORDA   — Three signals fused by Borda rank:
               borda = w1*B(S1) + w2*B(S2) + w3*B(S3),  w1=0.16, w2=0.08, w3=0.76
             Computed per-species (intra-species ranking) AND pooled cross-species.

HYPOTHESES
   H17a:  Schizophyllum top-K anomaly candidates cluster in the global tier
          (large global-spike windows) — its rare global firings are the
          kingdom-level choice collapse events.
   H17b:  Cordyceps contributes the most anomaly candidates to the cross-species
          top-20 relative to dataset size (anomaly density ↑ with ecological stress).
   H17c:  The cross-species centroid S3 signal is dominated by w3 Borda weight,
          and the top-K anomalies have significantly higher rho than background.

Outputs (experiments/results/)
   17_embedding_pca.png        PCA of all 5 species pooled in 40D, top anomalies marked
   17_borda_per_species.png    Per-species Borda score distributions (violin)
   17_top20_crossspecies.png   Ranked table: top-20 cross-species anomalies
   17_score_components.png     S1/S2/S3 contribution heatmap for top-20
   17_anomaly_temporal.png     Temporal position of top-10 per species + their tier
   17_report.json              All metrics and ranked lists
"""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import uniform_filter1d

# ── Paths ────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
PLEU_EMB    = MYCO / "tools/data/pleurotus_spike_emb_40d.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants ──────────────────────────────────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600       # 10-min non-overlapping windows (consistent with Exps 14-16)
DT_S             = 1.0       # sample interval for multispecies ADC-24 data (1 Hz)

# ── Borda weights (mirrors borda_rank.py) ────────────────────────────────────
W1, W2, W3 = 0.16, 0.08, 0.76   # FMM, JEPA, centroid-cosine

# ── FMM / kNN parameters (mirrors fmm_score.py) ──────────────────────────────
K_NEIGHBORS = 10
CHUNK_SIZE  = 512

# ── JEPA parameters (consistent with Exps 15/16) ────────────────────────────
SEQ_LEN        = 16
PRED_HORIZON   = 4
HIDDEN_DIM     = 32
LEARNING_RATE  = 3e-3
N_EPOCHS       = 60
MIN_WINDOWS    = 64
DEGENERATE_STD = 1e-4

# ── 40D histogram parameters (mirrors spike_embed.py) ───────────────────────
N_BINS    = 10
AMP_EDGES   = np.logspace(np.log10(0.03), np.log10(2.1), N_BINS + 1)
ISI_EDGES   = np.logspace(0, np.log10(3600), N_BINS + 1)
BURST_EDGES = np.array([1, 2, 3, 4, 5, 6, 8, 10, 14, 18, 100])

# ── Species ───────────────────────────────────────────────────────────────────
MULTISPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
    "Enoki":         "Enoki fungi Flammulina velutipes.txt",
    "Omphalotus":    "Ghost Fungi Omphalotus nidiformis.txt",
}

SPECIES_META_STRESS = {
    "Pleurotus":     0.35,
    "Schizophyllum": 0.90,
    "Cordyceps":     0.70,
    "Enoki":         0.50,
    "Omphalotus":    0.20,
}

SPECIES_COLOURS = {
    "Pleurotus":     "#2ecc71",
    "Schizophyllum": "#9b59b6",
    "Cordyceps":     "#e67e22",
    "Enoki":         "#27ae60",
    "Omphalotus":    "#2980b9",
}

TOP_K = 20   # cross-species anomaly candidates to report


# ─────────────────────────────────────────────────────────────────────────────
# 40D HISTOGRAM EMBEDDING  (inlined from spike_embed.py + export_spike_events.py)
# ─────────────────────────────────────────────────────────────────────────────

def make_histogram(spikes_in_window: list[dict]) -> np.ndarray:
    """
    Build L1-normalised 40D histogram embedding from a list of spike dicts.
    Layout:  dims  0- 9 amplitude deciles,
             dims 10-19 ISI deciles (log-spaced),
             dims 20-29 burst-size deciles,
             dims 30-39 channel activity (8 ch + total + burst-active).
    Mirrors tools/spike_embed.py::make_histogram exactly.
    """
    if not spikes_in_window:
        return np.zeros(40, dtype=np.float32)

    amps  = np.array([s["amplitude_mv"]  for s in spikes_in_window], dtype=np.float64)
    isis  = np.array([s["preceding_isi_s"] for s in spikes_in_window
                      if s.get("preceding_isi_s") is not None], dtype=np.float64)
    burst = np.array([3 if s["in_burst"] else 1 for s in spikes_in_window], dtype=np.float64)
    chans = np.array([s["channel"] for s in spikes_in_window], dtype=np.int32)

    h_amp  = np.histogram(amps, bins=AMP_EDGES)[0].astype(float)
    h_isi  = np.histogram(isis if len(isis) > 0 else np.array([1.0]),
                          bins=ISI_EDGES)[0].astype(float)
    h_bst  = np.histogram(burst, bins=BURST_EDGES)[0].astype(float)

    h_chan = np.zeros(10, dtype=float)
    for ch in range(8):
        h_chan[ch] = np.sum(chans == ch)
    h_chan[8] = float(len(spikes_in_window))
    h_chan[9] = float(sum(s["in_burst"] for s in spikes_in_window))

    h = np.concatenate([h_amp, h_isi, h_bst, h_chan])
    total = h.sum()
    if total > 0:
        h /= total
    return h.astype(np.float32)


def detect_spikes_from_voltage(det: np.ndarray, dt_s: float = 1.0) -> list[dict]:
    """
    Threshold-crossing spike detector on detrended voltage matrix.
    Returns list of spike dicts compatible with make_histogram.
    Mirrors the core of export_spike_events.py::detect_spikes.

    Parameters
    ----------
    det  : (n_samples, n_channels) float32 detrended voltage
    dt_s : sample interval in seconds
    """
    n_rows, n_chan = det.shape
    all_spikes: list[dict] = []

    for ch in range(n_chan):
        sig = det[:, ch].astype(np.float64)
        above = sig >= GLOBAL_THRESH_MV
        # Leading edges only
        padded = np.concatenate([[False], above])
        starts = np.where(np.diff(padded.astype(int)) == 1)[0]

        prev_t_s: Optional[float] = None
        for idx in starts:
            # Find peak within ±25 samples around the leading edge
            lo = max(0, idx - 25)
            hi = min(n_rows, idx + 100)
            local = sig[lo:hi]
            peak_off = int(np.argmax(local))
            amp = float(local[peak_off])
            t_s = (lo + peak_off) * dt_s
            isi = (t_s - prev_t_s) if prev_t_s is not None else None
            in_burst = (isi is not None and isi < 60.0)
            all_spikes.append({
                "channel":       ch,
                "t_s":           t_s,
                "amplitude_mv":  amp,
                "preceding_isi_s": isi,
                "in_burst":      in_burst,
            })
            prev_t_s = t_s

    all_spikes.sort(key=lambda s: s["t_s"])
    return all_spikes


def build_spike_embeddings(det: np.ndarray, dt_s: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert detrended voltage matrix to 40D histogram embeddings + density stats.

    Returns
    -------
    embeddings : (n_windows, 40) float32 — L1-normalised histogram per window
    rho_local  : (n_windows,) float32 — sub-threshold event density
    rho_global : (n_windows,) float32 — global spike density
    """
    n_rows, n_chan = det.shape
    n_win = n_rows // WINDOW_S
    print(f"    Detecting spikes …", end=" ", flush=True)
    spikes = detect_spikes_from_voltage(det, dt_s)
    print(f"{len(spikes)} spikes detected", flush=True)

    # Bucketise spikes into windows
    windows: list[list[dict]] = [[] for _ in range(n_win)]
    win_size_s = WINDOW_S * dt_s
    for sp in spikes:
        w_idx = int(sp["t_s"] / win_size_s)
        if w_idx < n_win:
            windows[w_idx].append(sp)

    # Build histograms
    embeddings = np.array([make_histogram(windows[i]) for i in range(n_win)],
                          dtype=np.float32)

    # Windowed density (same as Exps 14-16)
    abs_v    = np.abs(det)
    sigma    = np.nanmedian(abs_v, axis=0) / 0.6745
    floor    = (NOISE_FLOOR_K * sigma).astype(np.float32)
    noise_m  = abs_v < floor[np.newaxis, :]
    global_m = abs_v >= GLOBAL_THRESH_MV
    local_m  = (~noise_m) & (~global_m)

    trim_l = local_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    trim_g = global_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    rho_local  = trim_l.mean(axis=(1, 2)).astype(np.float32)
    rho_global = trim_g.mean(axis=(1, 2)).astype(np.float32)

    return embeddings, rho_local, rho_global


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADERS
# ─────────────────────────────────────────────────────────────────────────────

def load_pleurotus_embeddings() -> tuple[np.ndarray, np.ndarray]:
    """
    Load Pleurotus 40D histogram embeddings from pleurotus_spike_emb_40d.ndjson.
    Returns (embeddings (N,40), rho (N,)) sorted by t_start_s.
    """
    if not PLEU_EMB.exists():
        raise FileNotFoundError(f"Pleurotus 40D embeddings not found: {PLEU_EMB}")
    records = [json.loads(l) for l in open(PLEU_EMB)]
    records.sort(key=lambda r: r["t_start_s"])
    embeddings = np.array([r["embedding"] for r in records], dtype=np.float32)
    rho        = np.array([r["rho"]       for r in records], dtype=np.float32)
    print(f"    Loaded {len(records)} windows  embedding shape={embeddings.shape}", flush=True)
    return embeddings, rho


def load_multispecies_embeddings(name: str) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Load a multispecies .txt file, detrend, detect spikes, build 40D histograms.
    Returns (embeddings (N,40), rho_local (N,), rho_global (N,)) or None.
    """
    path = DATA_DIR / MULTISPECIES_FILES[name]
    if not path.exists():
        print(f"  [SKIP] {name}: file not found", flush=True)
        return None
    print(f"  Loading {name} …", flush=True)
    df = pd.read_csv(str(path), sep="\t", header=0, dtype=np.float32,
                     na_values=["NaN", "", " "], engine="c", on_bad_lines="skip")
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

    print(f"    {n:,} rows × {nchan} ch  |  {n // WINDOW_S} windows", flush=True)
    emb, rho_l, rho_g = build_spike_embeddings(det, dt_s=DT_S)
    return emb, rho_l, rho_g


# ─────────────────────────────────────────────────────────────────────────────
# FMM INTRA-SPECIES SCORE  (inlined from fmm_score.py)
# ─────────────────────────────────────────────────────────────────────────────

def build_knn_graph(embeddings: np.ndarray, k: int = K_NEIGHBORS):
    """
    kNN graph over N × 40D embeddings.
    Returns (knn_indices (N,k), knn_distances (N,k)).
    Mirrors fmm_score.py::build_knn_graph exactly.
    """
    n, d = embeddings.shape
    knn_idx  = np.zeros((n, k), dtype=np.int32)
    knn_dist = np.zeros((n, k), dtype=np.float32)

    for start in range(0, n, CHUNK_SIZE):
        end   = min(start + CHUNK_SIZE, n)
        chunk = embeddings[start:end]
        sq_a = (chunk ** 2).sum(axis=1, keepdims=True)
        sq_b = (embeddings ** 2).sum(axis=1, keepdims=True)
        dist2 = sq_a + sq_b.T - 2.0 * (chunk @ embeddings.T)
        np.clip(dist2, 0.0, None, out=dist2)
        for li in range(end - start):
            gi  = start + li
            row = dist2[li].copy()
            row[gi] = np.inf
            top = np.argpartition(row, k)[:k]
            ts  = top[np.argsort(row[top])]
            knn_idx[gi]  = ts
            knn_dist[gi] = np.sqrt(np.maximum(row[ts], 0.0))

    return knn_idx, knn_dist


def compute_fmm_scores(embeddings: np.ndarray, k: int = K_NEIGHBORS) -> np.ndarray:
    """
    Compute FMM deviation scores for a set of embeddings.
    Returns (N,) float32 array of per-window FMM anomaly scores.
    Mirrors fmm_score.py::main logic.
    """
    n = len(embeddings)
    if n <= k:
        return np.zeros(n, dtype=np.float32)

    knn_idx, knn_dist = build_knn_graph(embeddings, k)

    # rho_i = 1 / mean_kNN_distance
    mean_dist  = knn_dist.mean(axis=1)
    mean_dist  = np.where(mean_dist < 1e-12, 1e-12, mean_dist)
    rho        = 1.0 / mean_dist
    rho_max    = rho.max() + 1e-12

    fmm_scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        deviations = []
        for j in knn_idx[i]:
            w_exp = ((rho[i] + rho[j]) / 2.0) / rho_max
            dist  = knn_dist[i, np.where(knn_idx[i] == j)[0][0]] if j in knn_idx[i] else 0.0
            w_act = rho[i] / rho_max / (dist + 1e-12) * (knn_dist[i].mean() + 1e-12)
            if w_exp > 1e-12:
                deviations.append(abs(w_act - w_exp) / w_exp)
        fmm_scores[i] = float(np.mean(deviations)) if deviations else 0.0

    # Normalise to [0, 1] (same as fmm_score.py normalisation step)
    fmm_max = fmm_scores.max() + 1e-9
    return (fmm_scores / fmm_max).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# S3 — COSINE DISTANCE TO CROSS-SPECIES CENTROID  (from borda_rank.py)
# ─────────────────────────────────────────────────────────────────────────────

def cosine_scores_from_centroid(embeddings: np.ndarray,
                                centroid: np.ndarray) -> np.ndarray:
    """
    Cosine distance of each row of `embeddings` to the given unit-norm `centroid`.
    Zero-norm rows → distance = 0 (background, not anomalous).
    """
    norms = np.linalg.norm(embeddings, axis=1)
    zero  = norms < 1e-12
    safe  = norms.copy()
    safe[zero] = 1.0
    e_n = embeddings / safe[:, np.newaxis]
    cos = e_n @ centroid
    dist = 1.0 - cos
    dist[zero] = 0.0
    return dist.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# JEPA S2  (TinyJEPA from Exps 15/16)
# ─────────────────────────────────────────────────────────────────────────────

def relu(x): return np.maximum(0.0, x)


class TinyJEPA:
    def __init__(self):
        s1 = np.sqrt(2.0 / SEQ_LEN)
        self.We = np.random.randn(HIDDEN_DIM, SEQ_LEN).astype(np.float32) * s1
        self.be = np.zeros(HIDDEN_DIM, dtype=np.float32)
        s2 = np.sqrt(2.0 / HIDDEN_DIM)
        self.Wp = np.random.randn(HIDDEN_DIM, HIDDEN_DIM).astype(np.float32) * s2
        self.bp = np.zeros(HIDDEN_DIM, dtype=np.float32)
        self.Wo = np.random.randn(1, HIDDEN_DIM).astype(np.float32) * s2
        self.bo = np.zeros(1, dtype=np.float32)

    def _enc(self, x): return relu(self.We @ x + self.be)
    def _pred(self, h): return relu(self.Wp @ h + self.bp)
    def _dec(self, h): return float((self.Wo @ h + self.bo)[0])

    def _step(self, ctx, target):
        h_e = self._enc(ctx); h_p = self._pred(h_e)
        pred = self._dec(h_p); err = pred - target
        dWo = err * h_p[np.newaxis, :]; dbo = np.array([err])
        dp = (self.Wo.T * err).flatten() * (h_p > 0)
        dWp = np.outer(dp, h_e); dbp = dp
        de = (self.Wp.T @ dp) * (h_e > 0)
        dWe = np.outer(de, ctx); dbe = de
        reg = 1e-4
        self.Wo -= LEARNING_RATE * (dWo + reg * self.Wo)
        self.bo -= LEARNING_RATE * dbo
        self.Wp -= LEARNING_RATE * (dWp + reg * self.Wp)
        self.bp -= LEARNING_RATE * dbp
        self.We -= LEARNING_RATE * (dWe + reg * self.We)
        self.be -= LEARNING_RATE * dbe
        return err ** 2

    def _norm(self, s): return (s - s.mean()) / (s.std() + 1e-8)

    def train(self, series):
        s = self._norm(series).astype(np.float32)
        n = len(s)
        for _ in range(N_EPOCHS):
            idx = list(range(n - SEQ_LEN - PRED_HORIZON))
            np.random.shuffle(idx)
            for i in idx:
                self._step(s[i:i+SEQ_LEN], float(s[i+SEQ_LEN+PRED_HORIZON-1]))

    def surprise(self, series):
        s = self._norm(series).astype(np.float32)
        n = len(s)
        return np.array(
            [(self._dec(self._pred(self._enc(s[i:i+SEQ_LEN]))) -
              float(s[i+SEQ_LEN+PRED_HORIZON-1])) ** 2
             for i in range(n - SEQ_LEN - PRED_HORIZON)],
            dtype=np.float32,
        )


def jepa_surprise_scores(combined_density: np.ndarray) -> Optional[np.ndarray]:
    """
    Train TinyJEPA on combined rho and return per-window surprise.
    Returns None if series is degenerate.
    Pads output with zeros to match original series length.
    """
    n = len(combined_density)
    if n < MIN_WINDOWS or np.std(combined_density) < DEGENERATE_STD:
        return None
    np.random.seed(42)
    model = TinyJEPA()
    model.train(combined_density)
    scores = model.surprise(combined_density)
    n_scored = len(scores)
    # Pad front (no context for first SEQ_LEN + PRED_HORIZON windows) with mean
    pad = np.full(n - n_scored, scores.mean(), dtype=np.float32)
    return np.concatenate([pad, scores])


# ─────────────────────────────────────────────────────────────────────────────
# BORDA FUSION  (inlined from borda_rank.py)
# ─────────────────────────────────────────────────────────────────────────────

def borda_fuse(s1: np.ndarray, s2: np.ndarray, s3: np.ndarray) -> np.ndarray:
    """Borda rank fusion of three scores.  Returns fused score in [0, 1]."""
    n = len(s1)
    def to_borda(v):
        r = np.argsort(np.argsort(v)).astype(float)
        return r / (n - 1) if n > 1 else np.zeros(n)
    return W1 * to_borda(s1) + W2 * to_borda(s2) + W3 * to_borda(s3)


# ─────────────────────────────────────────────────────────────────────────────
# PCA (pure numpy — no sklearn dependency)
# ─────────────────────────────────────────────────────────────────────────────

def pca_2d(X: np.ndarray) -> np.ndarray:
    """Project (N, D) matrix to (N, 2) via top-2 PCA components."""
    X = X.astype(np.float64)
    X -= X.mean(axis=0)
    cov = X.T @ X / (len(X) - 1)
    vals, vecs = np.linalg.eigh(cov)
    top2 = vecs[:, np.argsort(vals)[::-1][:2]]
    return (X @ top2).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_embedding_pca(all_embeddings: dict[str, np.ndarray],
                      top_anomalies: list[dict]) -> None:
    pool = np.vstack([v for v in all_embeddings.values()])
    projected = pca_2d(pool)

    offset = 0
    species_order = list(all_embeddings.keys())
    fig, ax = plt.subplots(figsize=(11, 8))

    for sp in species_order:
        n  = len(all_embeddings[sp])
        xy = projected[offset:offset + n]
        ax.scatter(xy[:, 0], xy[:, 1], s=6, alpha=0.35,
                   color=SPECIES_COLOURS.get(sp, "#888"), label=sp)
        offset += n

    # Mark top anomalies
    for entry in top_anomalies[:TOP_K]:
        sp   = entry["species"]
        widx = entry["window_idx"]
        # Recompute offset for this species
        off = sum(len(all_embeddings[s]) for s in species_order
                  if s != sp and species_order.index(s) < species_order.index(sp))
        pt = projected[off + widx]
        ax.scatter(pt[0], pt[1], s=120, color="red",
                   edgecolor="black", linewidth=0.8, zorder=5)

    ax.set_xlabel("PC1", fontsize=11)
    ax.set_ylabel("PC2", fontsize=11)
    ax.set_title(
        "PCA of Pooled 40D Spike-Histogram Embeddings — All 5 Species\n"
        f"Red dots = top-{TOP_K} cross-species anomaly candidates",
        fontsize=11,
    )
    ax.legend(fontsize=9, markerscale=3)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    out = str(RESULTS_DIR / "17_embedding_pca.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 17_embedding_pca.png", flush=True)


def fig_borda_per_species(species_scores: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(11, 5))
    positions = list(range(len(species_scores)))
    names = list(species_scores.keys())
    data  = [species_scores[n] for n in names]
    colours = [SPECIES_COLOURS.get(n, "#888") for n in names]

    vp = ax.violinplot(data, positions=positions,
                       showmedians=True, showextrema=True)
    for i, (body, col) in enumerate(zip(vp["bodies"], colours)):
        body.set_facecolor(col)
        body.set_alpha(0.7)
    vp["cmedians"].set_color("black")
    vp["cmaxes"].set_color("black")
    vp["cmins"].set_color("black")
    vp["cbars"].set_color("black")

    # Annotate top-1 per species
    for i, (n, s) in enumerate(species_scores.items()):
        ax.text(i, s.max() + 0.01, f"max={s.max():.3f}",
                ha="center", fontsize=7, color=SPECIES_COLOURS.get(n, "black"))

    ax.set_xticks(positions)
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel("Intra-species Borda score", fontsize=11)
    ax.set_title(
        "Borda Anomaly Score Distributions per Species\n"
        "(w1=0.16 FMM + w2=0.08 JEPA + w3=0.76 centroid cosine)",
        fontsize=11,
    )
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "17_borda_per_species.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 17_borda_per_species.png", flush=True)


def fig_top20_crossspecies(top_anomalies: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(13, 8))
    ranks = list(range(1, len(top_anomalies) + 1))
    borda_scores = [e["borda_cross"] for e in top_anomalies]
    colours      = [SPECIES_COLOURS.get(e["species"], "#888") for e in top_anomalies]

    bars = ax.barh(ranks[::-1], borda_scores[::-1],
                   color=colours[::-1], edgecolor="white", height=0.7)

    for i, (entry, rank) in enumerate(zip(top_anomalies, ranks)):
        ax.text(entry["borda_cross"] + 0.002,
                len(top_anomalies) + 1 - rank,
                f"{entry['species']}  w{entry['window_idx']}  "
                f"S1={entry['s1']:.3f} S2={entry['s2']:.3f} S3={entry['s3']:.3f}",
                va="center", fontsize=7)

    ax.set_yticks(ranks[::-1])
    ax.set_yticklabels([f"#{r}" for r in ranks], fontsize=8)
    ax.set_xlabel("Cross-species Borda score", fontsize=11)
    ax.set_title(
        f"Top-{TOP_K} Cross-Species Anomaly Candidates\n"
        "Pooled Borda ranking across all 5 species",
        fontsize=11,
    )
    patches = [mpatches.Patch(color=SPECIES_COLOURS.get(n, "#888"), label=n)
               for n in SPECIES_COLOURS]
    ax.legend(handles=patches, fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "17_top20_crossspecies.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 17_top20_crossspecies.png", flush=True)


def fig_score_components(top_anomalies: list[dict]) -> None:
    """Stacked horizontal bar: S1, S2, S3 contribution per top-20."""
    n = len(top_anomalies)
    s1 = np.array([e["s1"] for e in top_anomalies])
    s2 = np.array([e["s2"] for e in top_anomalies])
    s3 = np.array([e["s3"] for e in top_anomalies])
    total = W1 * s1 + W2 * s2 + W3 * s3
    total = np.where(total < 1e-9, 1e-9, total)
    frac1 = W1 * s1 / total
    frac2 = W2 * s2 / total
    frac3 = W3 * s3 / total

    ys = list(range(1, n + 1))[::-1]
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.barh(ys, frac1, height=0.7, color="#3498db", label=f"S1 × {W1} (FMM deviation)")
    ax.barh(ys, frac2, height=0.7, color="#2ecc71", left=frac1,
            label=f"S2 × {W2} (JEPA surprise)")
    ax.barh(ys, frac3, height=0.7, color="#e74c3c", left=frac1 + frac2,
            label=f"S3 × {W3} (centroid cosine)")

    ax.set_yticks(ys)
    ax.set_yticklabels(
        [f"#{r+1} {top_anomalies[r]['species'][:5]}" for r in range(n)][::-1],
        fontsize=8,
    )
    ax.set_xlabel("Fractional score contribution", fontsize=11)
    ax.set_title(
        "Score Component Breakdown for Top-20 Cross-Species Anomalies\n"
        "Relative contribution of each signal to the Borda score",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(0, 1.05)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "17_score_components.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 17_score_components.png", flush=True)


def fig_anomaly_temporal(species_data: dict,
                         top_anomalies: list[dict]) -> None:
    """
    For each species: rho_combined timeline with top-5 anomaly windows marked.
    """
    names_with_data = [sp for sp in SPECIES_COLOURS if sp in species_data]
    n_sp = len(names_with_data)
    fig, axes = plt.subplots(n_sp, 1, figsize=(14, 3 * n_sp), sharex=False)
    if n_sp == 1:
        axes = [axes]

    for ax, sp in zip(axes, names_with_data):
        rho_l = species_data[sp]["rho_local"]
        rho_g = species_data[sp]["rho_global"]
        combined = (rho_l + rho_g).astype(np.float32)
        ts = np.arange(len(combined))
        ax.fill_between(ts, 0, rho_l, alpha=0.5,
                        color="#f39c12", label="ρ local")
        ax.fill_between(ts, rho_l, combined, alpha=0.5,
                        color="#e74c3c", label="ρ global")

        # Mark top-5 anomalies for this species
        sp_anomalies = [e for e in top_anomalies if e["species"] == sp][:5]
        for i, entry in enumerate(sp_anomalies):
            widx = entry["window_idx"]
            if widx < len(combined):
                ax.axvline(widx, color="black", linewidth=1.5, linestyle="--",
                           alpha=0.8)
                ax.text(widx + 0.5, combined.max() * 0.85,
                        f"#{entry['cross_rank']}", fontsize=7, color="black")

        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(f"{sp}  —  ρ_local + ρ_global  "
                     f"(vertical lines = top cross-species anomalies)", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.2)

    plt.suptitle("Anomaly Temporal Placement per Species",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    out = str(RESULTS_DIR / "17_anomaly_temporal.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 17_anomaly_temporal.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 17 — Cross-Species Anomaly Scan (myco-explorer Borda)")
    print("=" * 70)

    # ── 1. Build embeddings for all species ───────────────────────────────
    print("\n--- Step 1: Build 40D spike-histogram embeddings ---", flush=True)
    all_embeddings:  dict[str, np.ndarray] = {}
    all_rho_local:   dict[str, np.ndarray] = {}
    all_rho_global:  dict[str, np.ndarray] = {}

    # Pleurotus from precomputed NDJSON
    print(f"\n[Pleurotus]", flush=True)
    pleu_emb, pleu_rho = load_pleurotus_embeddings()
    all_embeddings["Pleurotus"] = pleu_emb
    # Pleurotus doesn't have tier split in the 40D file; use rho as combined proxy
    all_rho_local["Pleurotus"]  = pleu_rho * 0.5   # proxy: assume equal split
    all_rho_global["Pleurotus"] = pleu_rho * 0.5

    # Multispecies from raw .txt
    for name in MULTISPECIES_FILES:
        print(f"\n[{name}]", flush=True)
        result = load_multispecies_embeddings(name)
        if result is None:
            continue
        emb, rho_l, rho_g = result
        all_embeddings[name]  = emb
        all_rho_local[name]   = rho_l
        all_rho_global[name]  = rho_g
        print(f"    Embeddings: {emb.shape}  rho_local mean={rho_l.mean():.4f}  "
              f"rho_global mean={rho_g.mean():.4f}", flush=True)

    present = list(all_embeddings.keys())
    total_windows = sum(len(all_embeddings[sp]) for sp in present)
    print(f"\n  Total windows across {len(present)} species: {total_windows:,}", flush=True)

    # ── 2. Cross-species centroid ──────────────────────────────────────────
    print("\n--- Step 2: Cross-species centroid → S3 ---", flush=True)
    pool_emb = np.vstack([all_embeddings[sp] for sp in present])
    c = pool_emb.mean(axis=0)
    c_norm = np.linalg.norm(c)
    centroid_unit = (c / c_norm) if c_norm > 1e-12 else c

    all_s3: dict[str, np.ndarray] = {}
    for sp in present:
        all_s3[sp] = cosine_scores_from_centroid(all_embeddings[sp], centroid_unit)
        print(f"  {sp:20s}  S3 mean={all_s3[sp].mean():.4f}  "
              f"max={all_s3[sp].max():.4f}", flush=True)

    # ── 3. FMM S1 per species ─────────────────────────────────────────────
    print("\n--- Step 3: Intra-species FMM (S1) ---", flush=True)
    all_s1: dict[str, np.ndarray] = {}
    for sp in present:
        n = len(all_embeddings[sp])
        print(f"  {sp}: {n} windows …", end=" ", flush=True)
        all_s1[sp] = compute_fmm_scores(all_embeddings[sp])
        print(f"mean={all_s1[sp].mean():.4f}", flush=True)

    # ── 4. JEPA S2 per species ────────────────────────────────────────────
    print("\n--- Step 4: JEPA surprise (S2) ---", flush=True)
    all_s2: dict[str, np.ndarray] = {}
    for sp in present:
        combined = all_rho_local[sp] + all_rho_global[sp]
        sc = jepa_surprise_scores(combined)
        if sc is None:
            n = len(combined)
            print(f"  {sp}: degenerate — using zeros", flush=True)
            sc = np.zeros(n, dtype=np.float32)
        else:
            print(f"  {sp}: mean S2={sc.mean():.4f}", flush=True)
        all_s2[sp] = sc

    # ── 5. Intra-species Borda ────────────────────────────────────────────
    print("\n--- Step 5: Borda fusion (intra- and cross-species) ---", flush=True)
    intra_borda: dict[str, np.ndarray] = {}
    for sp in present:
        # Align lengths (JEPA padding may create length mismatches)
        n = len(all_s1[sp])
        s1 = all_s1[sp][:n]
        s2 = all_s2[sp][:n] if len(all_s2[sp]) >= n else np.pad(all_s2[sp], (0, n - len(all_s2[sp])))
        s3 = all_s3[sp][:n]
        intra_borda[sp] = borda_fuse(s1, s2, s3)
        print(f"  {sp:20s}  intra-borda max={intra_borda[sp].max():.4f}  "
              f"mean={intra_borda[sp].mean():.4f}", flush=True)

    # ── 6. Cross-species pooled Borda ─────────────────────────────────────
    # Pool all windows and re-rank in the combined space
    all_s1_pool = np.concatenate([all_s1[sp] for sp in present])
    all_s2_pool = np.concatenate([all_s2[sp][:len(all_s1[sp])] for sp in present])
    all_s3_pool = np.concatenate([all_s3[sp] for sp in present])
    cross_borda  = borda_fuse(all_s1_pool, all_s2_pool, all_s3_pool)

    # Label each window with species + index
    labels = []
    for sp in present:
        labels.extend([(sp, i) for i in range(len(all_s1[sp]))])

    ranked_cross = sorted(enumerate(cross_borda.tolist()),
                          key=lambda x: -x[1])

    top_anomalies: list[dict] = []
    for cross_rank_0, (idx, borda_val) in enumerate(ranked_cross[:TOP_K]):
        sp, widx = labels[idx]
        n = len(all_s1[sp])
        entry = {
            "cross_rank":  cross_rank_0 + 1,
            "species":     sp,
            "window_idx":  widx,
            "borda_cross": round(borda_val, 6),
            "s1": round(float(all_s1[sp][widx]), 6),
            "s2": round(float(all_s2[sp][widx] if widx < len(all_s2[sp]) else 0), 6),
            "s3": round(float(all_s3[sp][widx]), 6),
            "rho_combined": round(float(all_rho_local[sp][widx] +
                                         all_rho_global[sp][widx]), 6),
            "intra_borda": round(float(intra_borda[sp][widx]), 6),
        }
        top_anomalies.append(entry)

    # ── 7. Figures ────────────────────────────────────────────────────────
    print("\n--- Step 6: Figures ---", flush=True)
    species_data = {sp: {"rho_local": all_rho_local[sp],
                         "rho_global": all_rho_global[sp]} for sp in present}
    fig_embedding_pca(all_embeddings, top_anomalies)
    fig_borda_per_species(intra_borda)
    fig_top20_crossspecies(top_anomalies)
    fig_score_components(top_anomalies)
    fig_anomaly_temporal(species_data, top_anomalies)

    # ── 8. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 17 SUMMARY")
    print("=" * 70)
    print(f"\n  Species in scan:  {len(present)}")
    print(f"  Total windows:    {total_windows:,}")
    print(f"\n  Top-{TOP_K} cross-species anomaly candidates:")
    print(f"  {'#':>3}  {'Species':20s}  {'Win':>6}  {'Borda':>7}  "
          f"{'S1':>7}  {'S2':>7}  {'S3':>7}  {'ρ_comb':>8}")
    print("  " + "-" * 73)

    sp_counts: dict[str, int] = {}
    for e in top_anomalies:
        sp_counts[e["species"]] = sp_counts.get(e["species"], 0) + 1
        print(f"  #{e['cross_rank']:>2}  {e['species']:20s}  {e['window_idx']:>6}  "
              f"{e['borda_cross']:>7.4f}  {e['s1']:>7.4f}  {e['s2']:>7.4f}  "
              f"{e['s3']:>7.4f}  {e['rho_combined']:>8.4f}")

    print(f"\n  Species contribution to top-{TOP_K}:")
    for sp in present:
        n_sp  = len(all_s1[sp])
        cnt   = sp_counts.get(sp, 0)
        dens  = cnt / n_sp * 1000  # per 1000 windows
        print(f"    {sp:20s}  {cnt:>2}/{TOP_K}  (anomaly density = {dens:.2f} per 1000 windows)")

    # Hypothesis checks
    print(f"\n  H17a — Schizophyllum top anomalies in global tier:")
    schiz_top = [e for e in top_anomalies if e["species"] == "Schizophyllum"]
    for e in schiz_top[:3]:
        w = e["window_idx"]
        g = all_rho_global["Schizophyllum"][w]
        l = all_rho_local["Schizophyllum"][w]
        tier = "global" if g > l else "local"
        print(f"    #{e['cross_rank']:>2}  w{w:>4}  ρ_global={g:.4f}  ρ_local={l:.4f}  → {tier}")

    most_dense = max(present, key=lambda s: sp_counts.get(s, 0) / len(all_s1[s]))
    print(f"\n  H17b — Species with highest anomaly density: {most_dense}")
    for sp in sorted(present, key=lambda s: -sp_counts.get(s, 0) / len(all_s1[s])):
        cnt  = sp_counts.get(sp, 0)
        dens = cnt / len(all_s1[sp]) * 1000
        stress = SPECIES_META_STRESS[sp]
        print(f"    {sp:20s}  density={dens:.3f}/1000  stress={stress}")

    print("\n" + "=" * 70)

    # ── 9. Save report ────────────────────────────────────────────────────
    report = {
        "species_window_counts": {sp: int(len(all_s1[sp])) for sp in present},
        "total_windows": int(total_windows),
        "top_anomalies": top_anomalies,
        "species_anomaly_counts": sp_counts,
        "intra_borda_stats": {
            sp: {"mean": float(intra_borda[sp].mean()),
                 "max":  float(intra_borda[sp].max()),
                 "p90":  float(np.percentile(intra_borda[sp], 90))}
            for sp in present
        },
        "s3_stats": {
            sp: {"mean": float(all_s3[sp].mean()),
                 "max":  float(all_s3[sp].max())}
            for sp in present
        },
    }
    out = str(RESULTS_DIR / "17_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report saved → 17_report.json")
    print("Done.  Figures → experiments/results/17_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
