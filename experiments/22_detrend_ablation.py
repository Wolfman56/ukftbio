"""
Experiment 22 — Linear De-Trending Ablation
============================================

Exp 21 confirmed that the η matrix rank-ordering is intrinsic (rs = 0.790 vs
Exp 19), surviving rank/probit normalisation.  But rank normalisation preserves
the *temporal ordering* of values — including any long-range superdiffusive trend.
Schizophyllum's H = 1.216 (Exp 20) means its rho series has a monotonic
territorial-expansion trend baked into the rank sequence; any model trained on a
stationary species cannot predict it.

Exp 21's largest effect was Schizophyllum-as-TARGET collapsing (Omphalotus→SZ:
0.529 → 0.222, Δ = −0.306).  But we cannot distinguish two explanations:

  (A) Schizophyllum's H > 1 trend is hard to predict → high transfer_S → low η
  (B) Schizophyllum has deep intrinsic grammar incompatibility regardless of trend

Exp 22 applies **linear de-trending** first, then rank normalisation:
  1. scipy.signal.detrend (type="linear") on each rho series
  2. rank_normalise (van der Waerden → probit) — same as Exp 21
  3. Re-run full N×N JEPA transfer matrix
  4. Three-way comparison: Exp 19 (raw) / Exp 21 (rank-norm) / Exp 22 (detrend+rank-norm)

If explanation (A) is correct:
  - Schizophyllum-as-target η should recover substantially toward Exp 19 levels
  - Schizophyllum Hurst H should drop toward 0.5 after detrend
  - Omphalotus→Cordyceps should stay high (their compatibility is not trend-driven)

If explanation (B) is correct:
  - Schizophyllum-as-target η stays low even after detrend
  - The grammar incompatibility is temporal-sequence-intrinsic, not trend-derived

HYPOTHESES
   H22a: Schizophyllum column mean η (mean over 4 incoming pairs) increases
         vs Exp 21 after detrend+rank-norm (trend was causing the collapse)
   H22b: Omphalotus→Cordyceps η remains ≥ 0.70 (intrinsic, not trend-mediated)
   H22c: Rank ordering of 20 pairs vs Exp 19 is preserved (rs ≥ 0.70)
   H22d: Schizophyllum H drops toward 0.5 after linear detrend (confirming
         that detrending actually removes the superdiffusion)

Figures
-------
  22_three_way_comparison.png  Three heatmaps: Exp 19 / Exp 21 / Exp 22
  22_sz_target_recovery.png    Schizophyllum column η across all three experiments
  22_hurst_before_after.png    Hurst exponent per species before and after detrend
  22_rank_ordering.png         Scatter: η_19 vs η_22 for all 20 pairs + rs annotation
  22_report.json               Full results + hypothesis tests
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
from scipy import stats
from scipy.ndimage import uniform_filter1d
from scipy.signal import detrend as scipy_detrend
from scipy.special import ndtri

# ── Paths ─────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
PLEU_EMB    = MYCO / "tools/data/pleurotus_spike_emb_40d.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (consistent with Exps 14–21) ────────────────────────────
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

# ── JEPA parameters (unchanged from Exps 15–21) ──────────────────────────────
SEQ_LEN      = 16
PRED_HORIZON = 4
HIDDEN_DIM   = 32
LR           = 3e-3
N_EPOCHS     = 60
MIN_WINDOWS  = 30
DEGENERATE   = 1e-4

# ── Baseline η matrices from previous experiments ─────────────────────────────
SPECIES = ["Pleurotus", "Schizophyllum", "Cordyceps", "Enoki", "Omphalotus"]

EXP19_ETA = np.array([
    [np.nan, 0.2866, 0.5495, 0.3368, 0.5076],
    [0.5920, np.nan, 0.3500, 0.3555, 0.4895],
    [0.6397, 0.4992, np.nan, 0.3737, 0.7263],
    [0.5405, 0.2379, 0.3700, np.nan, 0.4320],
    [0.6320, 0.5286, 0.7904, 0.3912, np.nan],
])

EXP21_ETA = np.array([
    [np.nan, 0.1731, 0.6232, 0.3513, 0.5340],
    [0.5090, np.nan, 0.4358, 0.2840, 0.4255],
    [0.6150, 0.2313, np.nan, 0.3898, 0.6910],
    [0.5540, 0.1700, 0.3730, np.nan, 0.3280],
    [0.6450, 0.2223, 0.7401, 0.3680, np.nan],
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
# TRANSFORMS
# ─────────────────────────────────────────────────────────────────────────────

def rank_normalise(x: np.ndarray) -> np.ndarray:
    """Van der Waerden rank transform → probit → Gaussian marginal."""
    n = len(x)
    ranks = stats.rankdata(x, method="average")
    quantiles = np.clip(ranks / (n + 1.0), 1e-9, 1 - 1e-9)
    return ndtri(quantiles).astype(np.float32)


def hurst_exponent(x: np.ndarray) -> float:
    """Hurst exponent via R/S analysis (log-log regression over multiple lags)."""
    n = len(x)
    lags = []
    rs_vals = []
    for lag in range(10, n // 2, max(1, n // 40)):
        chunks = [x[i:i+lag] for i in range(0, n - lag, lag)]
        if len(chunks) < 2:
            continue
        rs_chunk = []
        for c in chunks:
            m = np.mean(c)
            s = np.std(c, ddof=1)
            if s < 1e-10:
                continue
            dev = np.cumsum(c - m)
            R = dev.max() - dev.min()
            rs_chunk.append(R / s)
        if rs_chunk:
            lags.append(lag)
            rs_vals.append(np.mean(rs_chunk))
    if len(lags) < 4:
        return float("nan")
    log_lags = np.log(lags)
    log_rs   = np.log(rs_vals)
    H, _, _, _, _ = stats.linregress(log_lags, log_rs)
    return float(H)


def linear_detrend_then_rank(rho: np.ndarray) -> np.ndarray:
    """
    Step 1: Remove linear trend from rho series (scipy.signal.detrend type='linear').
    Step 2: Rank-normalise the residuals to Gaussian marginal.
    The detrending eliminates the H > 1 superdiffusive trend while preserving
    the residual spike structure.  Rank-norm then equalises marginal shape.
    """
    detrended = scipy_detrend(rho.astype(np.float64), type="linear").astype(np.float32)
    return rank_normalise(detrended)


# ─────────────────────────────────────────────────────────────────────────────
# TINY JEPA (identical to Exps 15–21; _norm is pass-through)
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

    # All normalisation done at load time — this is a pass-through
    def _norm(self, s): return s.astype(np.float64)

    def _enc(self, x): return np.maximum(0, self.We @ x + self.be)
    def _pred(self, h): return np.maximum(0, self.Wp @ h + self.bp)
    def _dec(self, h): return float((self.Wo @ h + self.bo)[0])

    def _step(self, ctx, tgt):
        h_e = self._enc(ctx); h_p = self._pred(h_e)
        p = self._dec(h_p); e = p - tgt
        dWo = e * h_p[np.newaxis, :]; dbo = np.array([e])
        dp  = (self.Wo.T * e).flatten() * (h_p > 0)
        dWp = np.outer(dp, h_e); dbp = dp
        de  = (self.Wp.T @ dp) * (h_e > 0)
        dWe = np.outer(de, ctx); dbe = de
        r   = 1e-4
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
# HURST DIAGNOSTIC: measure before and after detrend
# ─────────────────────────────────────────────────────────────────────────────

def measure_hurst_table(raw_series: dict[str, np.ndarray]) -> dict:
    table = {}
    for sp, rho in raw_series.items():
        H_before = hurst_exponent(rho)
        detrended = scipy_detrend(rho.astype(np.float64), type="linear").astype(np.float32)
        H_after   = hurst_exponent(detrended)
        table[sp]  = {"H_before": H_before, "H_after": H_after,
                      "delta_H": H_after - H_before}
        print(f"  {sp:14s}: H_before={H_before:.3f}  H_after={H_after:.3f}  "
              f"ΔH={H_after - H_before:+.3f}")
    return table


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def _heatmap(ax, eta: np.ndarray, species: list[str], title: str,
             vmin: float = 0.0, vmax: float = 1.0) -> None:
    n = len(species)
    masked = np.ma.masked_where(np.eye(n, dtype=bool), eta)
    im = ax.imshow(masked, cmap="RdYlGn", vmin=vmin, vmax=vmax, aspect="auto")
    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(eta[i, j]):
                val = eta[i, j]
                c   = "black" if 0.25 < val < 0.85 else "white"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=8, color=c, fontweight="bold")
            elif i == j:
                ax.text(j, i, "—", ha="center", va="center", fontsize=9, color="gray")
    labels = [f"{sp}\n({'B' if sp in {'Pleurotus','Schizophyllum'} else 'S'})"
              for sp in species]
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Target", fontsize=8); ax.set_ylabel("Source", fontsize=8)
    ax.set_title(title, fontsize=10)
    return im


def fig_three_way(eta22: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, eta, title in [
        (axes[0], EXP19_ETA, "Exp 19\n(raw z-score)"),
        (axes[1], EXP21_ETA, "Exp 21\n(rank-norm only)"),
        (axes[2], eta22,     "Exp 22\n(detrend + rank-norm)"),
    ]:
        im = _heatmap(ax, eta, SPECIES, title)

    # Shared colorbar on the right
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="η (JEPA transfer efficiency)")

    plt.suptitle("Three-way η Ablation: Raw → Rank-norm → Detrend+Rank-norm\n"
                 "H22a: Does detrending rescue Schizophyllum as target?", fontsize=12, y=1.02)
    out = str(RESULTS_DIR / "22_three_way_comparison.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 22_three_way_comparison.png")


def fig_sz_target_recovery(eta22: np.ndarray) -> None:
    """Focus on Schizophyllum column (target) across three experiments."""
    sz_idx = SPECIES.index("Schizophyllum")
    sources = [sp for sp in SPECIES if sp != "Schizophyllum"]
    src_idx = [SPECIES.index(s) for s in sources]

    eta19_col = EXP19_ETA[src_idx, sz_idx]
    eta21_col = EXP21_ETA[src_idx, sz_idx]
    eta22_col = eta22[src_idx, sz_idx]

    x  = np.arange(len(sources))
    w  = 0.25
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - w, eta19_col, w, label="Exp 19 (raw)", color="#3498db", alpha=0.85)
    ax.bar(x,     eta21_col, w, label="Exp 21 (rank-norm)", color="#e74c3c", alpha=0.85)
    ax.bar(x + w, eta22_col, w, label="Exp 22 (detrend+rank-norm)", color="#2ecc71", alpha=0.85)

    # Mean lines
    ax.axhline(eta19_col.mean(), color="#3498db", lw=1.5, ls="--", alpha=0.7)
    ax.axhline(eta21_col.mean(), color="#e74c3c", lw=1.5, ls="--", alpha=0.7)
    ax.axhline(eta22_col.mean(), color="#2ecc71", lw=1.5, ls="--", alpha=0.7)

    ax.set_xticks(x); ax.set_xticklabels(sources, fontsize=10)
    ax.set_ylabel("η (transfer efficiency into Schizophyllum)")
    ax.set_title("Schizophyllum as TARGET — η across Three Experiments\n"
                 "H22a: Does detrend+rank-norm recover vs Exp 21?", fontsize=11)
    ax.set_ylim(0, 1.0); ax.axhline(0.70, color="gray", ls=":", lw=1.2, label="η = 0.70")
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = str(RESULTS_DIR / "22_sz_target_recovery.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 22_sz_target_recovery.png")


def fig_hurst_before_after(hurst_table: dict) -> None:
    sps   = list(hurst_table.keys())
    H_bef = [hurst_table[s]["H_before"] for s in sps]
    H_aft = [hurst_table[s]["H_after"]  for s in sps]

    x   = np.arange(len(sps))
    w   = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - w/2, H_bef, w, label="H before detrend (raw rho)", color="#e74c3c", alpha=0.85)
    ax.bar(x + w/2, H_aft, w, label="H after linear detrend",      color="#2ecc71", alpha=0.85)

    ax.axhline(0.5, color="black", ls="--", lw=1.2, label="H = 0.5 (random walk)")
    ax.axhline(1.0, color="gray",  ls=":",  lw=1.0, label="H = 1.0 (superdiffusive boundary)")

    for i, sp in enumerate(sps):
        dH = hurst_table[sp]["delta_H"]
        y  = max(H_bef[i], H_aft[i]) + 0.03
        ax.text(i, y, f"ΔH={dH:+.3f}", ha="center", fontsize=8,
                color="navy" if dH < 0 else "darkgreen")

    ax.set_xticks(x); ax.set_xticklabels(sps, fontsize=10)
    ax.set_ylabel("Hurst Exponent H")
    ax.set_title("Hurst Exponent Before and After Linear Detrend\n"
                 "H22d: Does Schizophyllum H drop toward 0.5?", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = str(RESULTS_DIR / "22_hurst_before_after.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 22_hurst_before_after.png")


def fig_rank_ordering(eta22: np.ndarray) -> None:
    """Scatter of η_19 vs η_22 for all 20 non-diagonal pairs."""
    n = len(SPECIES)
    pairs_19, pairs_22, labels = [], [], []
    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(EXP19_ETA[i, j]):
                pairs_19.append(EXP19_ETA[i, j])
                pairs_22.append(eta22[i, j])
                labels.append(f"{SPECIES[i][:3]}→{SPECIES[j][:3]}")

    pairs_19 = np.array(pairs_19); pairs_22 = np.array(pairs_22)
    rs, pval = stats.spearmanr(pairs_19, pairs_22)

    fig, ax = plt.subplots(figsize=(8, 7))
    for i, (x, y, lab) in enumerate(zip(pairs_19, pairs_22, labels)):
        sp_src = SPECIES[[SPECIES[k][:3] == lab[:3] for k in range(len(SPECIES))].index(True)]
        c = SPECIES_COLOURS.get(sp_src, "#888888")
        ax.scatter(x, y, color=c, s=55, zorder=3)
        # label the Schizophyllum-target pairs (most interesting)
        if "Sch" in lab and "→" in lab and lab.index("Sch") > 0:
            ax.annotate(lab, (x, y), textcoords="offset points", xytext=(6, 3), fontsize=7)

    # highlight Omphalotus→Cordyceps
    om_co_19 = EXP19_ETA[SPECIES.index("Omphalotus"), SPECIES.index("Cordyceps")]
    om_co_22 = eta22[SPECIES.index("Omphalotus"), SPECIES.index("Cordyceps")]
    ax.scatter(om_co_19, om_co_22, color="gold", s=120, zorder=5,
               edgecolors="black", linewidths=1.5, label="Omphalotus→Cordyceps ★")

    lims = [min(pairs_19.min(), pairs_22.min()) - 0.05,
            max(pairs_19.max(), pairs_22.max()) + 0.05]
    ax.plot(lims, lims, "k--", lw=1, alpha=0.4, label="η_19 = η_22 (no change)")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("η (Exp 19 — raw z-score)")
    ax.set_ylabel("η (Exp 22 — detrend + rank-norm)")

    ax.axhline(0.70, color="gray", ls=":", lw=1, alpha=0.5)
    ax.axvline(0.70, color="gray", ls=":", lw=1, alpha=0.5)

    ax.set_title(f"Rank-Ordering Stability: Exp 19 vs Exp 22\n"
                 f"Spearman rs = {rs:.4f}  p = {pval:.2e}  (H22c threshold: rs ≥ 0.70)",
                 fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = str(RESULTS_DIR / "22_rank_ordering.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 22_rank_ordering.png")
    return rs, pval


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Experiment 22 — Linear De-Trending Ablation")
    print("=" * 70)

    # ── 1. Load raw rho series ─────────────────────────────────────────────
    print("\n[1/5] Loading raw rho series...")
    raw_rho: dict[str, np.ndarray] = {}
    raw_rho["Pleurotus"] = load_pleurotus_rho()
    for sp in ["Schizophyllum", "Cordyceps", "Enoki", "Omphalotus"]:
        raw_rho[sp] = load_multispecies_rho(sp)

    # ── 2. Hurst diagnostic before/after detrend ──────────────────────────
    print("\n[2/5] Measuring Hurst exponents before and after linear detrend...")
    hurst_table = measure_hurst_table(raw_rho)

    # ── 3. Apply detrend + rank-norm to all species ────────────────────────
    print("\n[3/5] Applying detrend + rank-norm to all species...")
    norm_rho: dict[str, np.ndarray] = {}
    for sp, rho in raw_rho.items():
        norm_rho[sp] = linear_detrend_then_rank(rho)
        k_before = float(stats.kurtosis(rho, fisher=True))
        k_after  = float(stats.kurtosis(norm_rho[sp], fisher=True))
        print(f"  {sp:14s}: κ_raw={k_before:+.3f}  κ_after={k_after:+.3f}")

    # ── 4. Train JEPA models on detrended+rank-normalised series ──────────
    print("\n[4/5] Training JEPA models (detrend+rank-norm) and computing N×N matrix...")
    n = len(SPECIES)
    eta22     = np.full((n, n), np.nan)
    native_s  = {}
    transfer_s = {}
    models = {}

    # Train each source model
    for sp in SPECIES:
        print(f"\n  Training source: {sp}")
        model, ns = train_jepa(norm_rho[sp], label=f"{sp} (detrend+rank)")
        models[sp] = model
        native_s[sp] = ns

    # Evaluate cross-species transfer
    print("\n  Evaluating transfer matrix...")
    for i, src in enumerate(SPECIES):
        if models[src] is None:
            continue
        for j, tgt in enumerate(SPECIES):
            if i == j:
                continue
            ts = models[src].mean_surprise(norm_rho[tgt])
            ns = native_s[tgt]
            if np.isnan(ts) or np.isnan(ns) or ns < DEGENERATE:
                continue
            eta22[i, j] = float(ns / ts)
            transfer_s[(src, tgt)] = ts
            print(f"    {src:14s} → {tgt:14s}  transfer_S={ts:.4f}  "
                  f"native_S={ns:.4f}  η={eta22[i, j]:.4f}")

    # ── 5. Hypothesis tests ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("HYPOTHESIS TESTS")
    print("=" * 70)

    sz_idx = SPECIES.index("Schizophyllum")
    src_idx = [i for i, s in enumerate(SPECIES) if s != "Schizophyllum"]

    # H22a: Schizophyllum-as-target column mean recovers vs Exp 21
    sz_col_21 = np.nanmean(EXP21_ETA[src_idx, sz_idx])
    sz_col_22 = np.nanmean(eta22[src_idx, sz_idx])
    h22a_pass = sz_col_22 > sz_col_21
    print(f"\nH22a: SZ column mean η: Exp21={sz_col_21:.4f}  Exp22={sz_col_22:.4f}  "
          f"Δ={sz_col_22-sz_col_21:+.4f}  → {'✓ Confirmed' if h22a_pass else '✗ Not confirmed'}")
    print(f"      (Detrend {'DID' if h22a_pass else 'did NOT'} rescue Schizophyllum as target)")

    # H22b: Omphalotus→Cordyceps remains ≥ 0.70
    om_idx = SPECIES.index("Omphalotus"); co_idx = SPECIES.index("Cordyceps")
    om_co_22 = eta22[om_idx, co_idx]
    h22b_pass = (not np.isnan(om_co_22)) and om_co_22 >= 0.70
    print(f"\nH22b: Omphalotus→Cordyceps: Exp19=0.7904  Exp21=0.7401  Exp22={om_co_22:.4f}  "
          f"≥0.70? → {'✓ Confirmed' if h22b_pass else '✗ Not confirmed'}")

    # H22c: Rank ordering vs Exp 19 preserved (rs ≥ 0.70)
    pairs_19, pairs_22 = [], []
    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(EXP19_ETA[i, j]) and not np.isnan(eta22[i, j]):
                pairs_19.append(EXP19_ETA[i, j])
                pairs_22.append(eta22[i, j])
    rs_22_19, p_rs = stats.spearmanr(pairs_19, pairs_22)
    h22c_pass = rs_22_19 >= 0.70
    print(f"\nH22c: Rank ordering rs(Exp19, Exp22) = {rs_22_19:.4f}  p = {p_rs:.2e}  "
          f"≥0.70? → {'✓ Confirmed' if h22c_pass else '✗ Not confirmed'}")
    print(f"      (For reference: Exp 21 rs = 0.7895)")

    # H22d: Schizophyllum Hurst drops toward 0.5
    sz_H_before = hurst_table["Schizophyllum"]["H_before"]
    sz_H_after  = hurst_table["Schizophyllum"]["H_after"]
    h22d_pass   = sz_H_after < sz_H_before
    closer_05   = abs(sz_H_after - 0.5) < abs(sz_H_before - 0.5)
    print(f"\nH22d: Schizophyllum H: before={sz_H_before:.3f}  after={sz_H_after:.3f}  "
          f"closer to 0.5? {'yes' if closer_05 else 'no'}  "
          f"→ {'✓ Confirmed' if h22d_pass else '✗ Not confirmed'}")

    # ── 6. Summary table ──────────────────────────────────────────────────
    print("\n" + "─" * 70)
    print("EXP 22 η MATRIX (detrend + rank-norm):")
    print(f"{'':16s}" + "".join(f" {sp:>13s}" for sp in SPECIES))
    for i, src in enumerate(SPECIES):
        row = f"  {src:14s}"
        for j in range(n):
            v = eta22[i, j]
            if np.isnan(v):
                row += f"        —    "
            elif i == j:
                row += f"        —    "
            else:
                star = " ★" if (i == om_idx and j == co_idx) else "  "
                row += f"    {v:.4f}{star}"
        print(row)

    print("\nΔη = Exp22 − Exp19 (top changes):")
    deltas = []
    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(EXP19_ETA[i, j]) and not np.isnan(eta22[i, j]):
                d = float(eta22[i, j] - EXP19_ETA[i, j])
                deltas.append((d, SPECIES[i], SPECIES[j]))
    deltas.sort(key=lambda x: x[0])
    for d, src, tgt in deltas[:5]:
        print(f"  {src:14s} → {tgt:14s}  Δη = {d:+.4f}")
    print("  ...")
    for d, src, tgt in deltas[-3:]:
        print(f"  {src:14s} → {tgt:14s}  Δη = {d:+.4f}")

    print("\nη range:")
    valid = eta22[~np.isnan(eta22)]
    rng19 = float(np.nanmax(EXP19_ETA) - np.nanmin(EXP19_ETA[~np.eye(n, dtype=bool)]))
    rng21 = float(np.nanmax(EXP21_ETA) - np.nanmin(EXP21_ETA[~np.eye(n, dtype=bool)]))
    rng22 = float(valid.max() - valid.min())
    print(f"  Exp 19: [{float(np.nanmin(EXP19_ETA[~np.eye(n, dtype=bool)])):.3f}, "
          f"{float(np.nanmax(EXP19_ETA)):.3f}]  range={rng19:.4f}")
    print(f"  Exp 21: [{float(np.nanmin(EXP21_ETA[~np.eye(n, dtype=bool)])):.3f}, "
          f"{float(np.nanmax(EXP21_ETA)):.3f}]  range={rng21:.4f}")
    print(f"  Exp 22: [{float(valid.min()):.3f}, {float(valid.max()):.3f}]"
          f"  range={rng22:.4f}")

    # ── 7. Figures ────────────────────────────────────────────────────────
    print("\n[5/5] Generating figures...")
    fig_three_way(eta22)
    fig_sz_target_recovery(eta22)
    fig_hurst_before_after(hurst_table)
    rs_fig, p_fig = fig_rank_ordering(eta22)

    # ── 8. Save JSON report ───────────────────────────────────────────────
    report = {
        "experiment": 22,
        "title":      "Linear De-Trending Ablation",
        "pipeline":   "detrend(linear) → rank_normalise → TinyJEPA",
        "species":    SPECIES,
        "eta_matrix": {
            SPECIES[i]: {SPECIES[j]: float(eta22[i, j]) if not np.isnan(eta22[i, j]) else None
                         for j in range(n)}
            for i in range(n)
        },
        "hurst_table": hurst_table,
        "key_pairs": {
            "Omphalotus_Cordyceps": {
                "exp19": 0.7904, "exp21": 0.7401, "exp22": float(om_co_22)
            },
            "Schizophyllum_column_mean": {
                "exp21": float(sz_col_21), "exp22": float(sz_col_22),
                "delta": float(sz_col_22 - sz_col_21)
            },
        },
        "hypotheses": {
            "H22a": {"result": "✓" if h22a_pass else "✗",
                     "description": "Schizophyllum-as-target column mean recovers vs Exp 21",
                     "sz_col_mean_exp21": float(sz_col_21),
                     "sz_col_mean_exp22": float(sz_col_22),
                     "delta": float(sz_col_22 - sz_col_21)},
            "H22b": {"result": "✓" if h22b_pass else "✗",
                     "description": "Omphalotus→Cordyceps η remains ≥ 0.70",
                     "eta_exp22": float(om_co_22)},
            "H22c": {"result": "✓" if h22c_pass else "✗",
                     "description": "Rank ordering rs(Exp19, Exp22) ≥ 0.70",
                     "rs": float(rs_22_19), "p": float(p_rs)},
            "H22d": {"result": "✓" if h22d_pass else "✗",
                     "description": "Schizophyllum H drops after linear detrend",
                     "H_before": float(sz_H_before), "H_after": float(sz_H_after),
                     "delta_H": float(sz_H_after - sz_H_before)},
        },
        "eta_range": {"exp19": rng19, "exp21": rng21, "exp22": rng22},
        "rank_ordering_rs_vs_exp19": float(rs_22_19),
    }
    out_json = str(RESULTS_DIR / "22_report.json")
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"  Saved 22_report.json")

    print("\n" + "=" * 70)
    print("DONE")
    print(f"  H22a (SZ target recovers):           {'✓' if h22a_pass else '✗'}")
    print(f"  H22b (OM→CO stays ≥ 0.70):          {'✓' if h22b_pass else '✗'}")
    print(f"  H22c (rank ordering rs ≥ 0.70):     {'✓' if h22c_pass else '✗'}")
    print(f"  H22d (SZ Hurst drops after detrend): {'✓' if h22d_pass else '✗'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
