"""
Experiment 20 — Temporal Grammar Fingerprint
=============================================

Exp 19 found 20 directed JEPA transfer pairs across 5 fungal species.
Only two pairs cleared η > 0.70 (Omphalotus↔Cordyceps), and neither
CV nor AR1 distance explained this at any significant level (both |rs| < 0.19).

Exp 20 asks: **why does that one pair work?**  We compute a 6-feature temporal
grammar fingerprint for each species — AR1, CV, Hurst exponent H (DFA), spectral
entropy S_ent, skewness γ, and kurtosis κ — then test which pairwise feature
distance best predicts η across all 20 directed pairs.

HYPOTHESES
   H20a:  Hurst exponent distance |ΔH| anti-correlates with η (rs < -0.4)
          more strongly than any other single feature — persistence structure
          is the true grammar fingerprint, not burstiness (CV) or inertia (AR1).
   H20b:  Omphalotus and Cordyceps share the most similar Hurst exponent
          across all 10 within-tier pairs, explaining the η = 0.790 breakthrough.
   H20c:  native_S(target) alone explains > 50% of η variance (Spearman)
          — target flatness is the dominant mechanical driver of low-η pairs
          (Schizophyllum native_S = 0.064 mechanically suppresses its column).

Secondary:
   H20d:  A composite Euclidean feature distance (all 6 normalised) predicts
          η better than any individual feature — no single axis dominates.

Figures
-------
  20_temporal_fingerprint_radar.png   Spider chart: 6 features per species
  20_feature_correlation_matrix.png   Spearman rs(|Δfeature|, η) bar chart
  20_hurst_eta_scatter.png            |ΔHurst| vs η scatter (all 20 pairs)
  20_omphalotus_cordyceps_spotlight.png  Feature profile of top pair vs worst pair
  20_report.json                      Full fingerprint + hypothesis test results
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
from scipy.signal import welch

# ── Paths ─────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
PLEU_EMB    = MYCO / "tools/data/pleurotus_spike_emb_40d.ndjson"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (consistent with Exps 14–19) ────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600
DT_S             = 1.0

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

# ── η matrix from Exp 19 ──────────────────────────────────────────────────────
EXP19_SPECIES = ["Pleurotus", "Schizophyllum", "Cordyceps", "Enoki", "Omphalotus"]
EXP19_ETA = np.array([
    # PL      SZ      CO      EN      OM
    [np.nan, 0.2866, 0.5495, 0.3368, 0.5076],   # Pleurotus
    [0.5920, np.nan, 0.3500, 0.3555, 0.4895],   # Schizophyllum
    [0.6397, 0.4992, np.nan, 0.3737, 0.7263],   # Cordyceps
    [0.5405, 0.2379, 0.3700, np.nan, 0.4320],   # Enoki
    [0.6320, 0.5286, 0.7904, 0.3912, np.nan],   # Omphalotus
])

EXP19_NATIVE_S = {
    "Pleurotus":     0.6541,
    "Schizophyllum": 0.0639,
    "Cordyceps":     0.1595,
    "Enoki":         0.4015,
    "Omphalotus":    0.3190,
}

CV_THRESHOLD  = 0.5
TIER_BURSTY   = {"Pleurotus", "Schizophyllum"}
TIER_SMOOTH   = {"Cordyceps", "Enoki", "Omphalotus"}


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
# TEMPORAL FEATURE EXTRACTION
# ─────────────────────────────────────────────────────────────────────────────

def dfa_hurst(x: np.ndarray, min_box: int = 4, max_box: int | None = None,
              n_scales: int = 20) -> float:
    """
    Detrended Fluctuation Analysis → Hurst exponent.
    H ≈ 0.5: uncorrelated; H > 0.5: persistent; H < 0.5: anti-persistent.
    """
    n = len(x)
    if max_box is None:
        max_box = n // 4
    max_box = max(max_box, min_box * 2)

    y = np.cumsum(x - x.mean())
    scales = np.unique(np.logspace(np.log10(min_box), np.log10(max_box),
                                    n_scales).astype(int))
    scales = scales[scales >= min_box]

    fluctuations = []
    valid_scales = []
    for s in scales:
        n_seg = n // s
        if n_seg < 2:
            continue
        segments = y[:n_seg * s].reshape(n_seg, s)
        t = np.arange(s)
        A = np.vstack([t, np.ones(s)]).T
        rms_list = []
        for seg in segments:
            c, _, _, _ = np.linalg.lstsq(A, seg, rcond=None)
            trend = A @ c
            rms_list.append(np.sqrt(np.mean((seg - trend) ** 2)))
        fluctuations.append(np.mean(rms_list))
        valid_scales.append(s)

    if len(valid_scales) < 3:
        return float("nan")

    log_s = np.log10(valid_scales)
    log_f = np.log10(fluctuations)
    slope, *_ = np.polyfit(log_s, log_f, 1)
    return float(slope)


def spectral_entropy(x: np.ndarray, fs: float = 1.0) -> float:
    """
    Normalised spectral entropy of the power spectral density.
    Low = dominated by a few frequencies; high = broad-band.
    """
    f, psd = welch(x, fs=fs, nperseg=min(256, len(x) // 4))
    psd = psd[psd > 0]
    psd /= psd.sum()
    return float(-np.sum(psd * np.log2(psd)) / np.log2(len(psd)))


def temporal_features(rho: np.ndarray, label: str) -> dict[str, float]:
    """
    Compute 6-feature temporal grammar fingerprint from a density window series.

    Features:
        ar1       lag-1 autocorrelation
        cv        coefficient of variation (std / mean)
        hurst     DFA Hurst exponent
        spec_ent  normalised spectral entropy
        skewness  Fisher skewness of rho distribution
        kurtosis  excess kurtosis (0 = normal)
    """
    x = rho.astype(np.float64)
    x_std = x.std()
    mean  = x.mean()

    ar1 = float(np.corrcoef(x[:-1], x[1:])[0, 1]) if len(x) > 2 else float("nan")
    cv  = float(x_std / mean) if mean > 0 else float("nan")
    H   = dfa_hurst(x)
    s_e = spectral_entropy(x)
    skw = float(stats.skew(x))
    krt = float(stats.kurtosis(x))

    feats = {"ar1": ar1, "cv": cv, "hurst": H, "spec_ent": s_e,
             "skewness": skw, "kurtosis": krt}
    print(f"  {label:15s}  AR1={ar1:.3f}  CV={cv:.3f}  H={H:.3f}  "
          f"S_ent={s_e:.3f}  skew={skw:.3f}  kurt={krt:.3f}")
    return feats


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

FEATURE_LABELS = {
    "ar1":      "AR1",
    "cv":       "CV",
    "hurst":    "Hurst H",
    "spec_ent": "Spectral\nEntropy",
    "skewness": "Skewness",
    "kurtosis": "Kurtosis",
}

def fig_radar(features: dict[str, dict]) -> None:
    """Spider/radar chart of normalised temporal features per species."""
    feat_keys = list(FEATURE_LABELS.keys())
    n_feat = len(feat_keys)
    species = list(features.keys())

    # Normalise features to [0, 1] across species
    raw = np.array([[features[sp][f] for f in feat_keys] for sp in species])
    mn  = raw.min(axis=0); mx = raw.max(axis=0)
    rng = np.where(mx - mn > 0, mx - mn, 1.0)
    norm = (raw - mn) / rng

    angles = np.linspace(0, 2 * np.pi, n_feat, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"polar": True})

    for sp, row in zip(species, norm):
        vals = row.tolist() + row[:1].tolist()
        c = SPECIES_COLOURS[sp]
        ax.plot(angles, vals, "-o", color=c, linewidth=2, markersize=5,
                label=sp, alpha=0.9)
        ax.fill(angles, vals, color=c, alpha=0.08)

    labels = [FEATURE_LABELS[f] for f in feat_keys]
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=7, alpha=0.6)
    ax.set_title("Temporal Grammar Fingerprint\n(6-feature radar, normalised per axis)",
                 fontsize=12, pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=9)
    plt.tight_layout()
    out = str(RESULTS_DIR / "20_temporal_fingerprint_radar.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 20_temporal_fingerprint_radar.png")


def fig_feature_correlation(species_list: list[str],
                             features: dict[str, dict],
                             eta_matrix: np.ndarray) -> dict[str, tuple]:
    """
    Spearman rs between |Δfeature| and η for all pairs.
    Also test native_S(target) vs η (H20c).
    Returns dict of feature → (rs, p) tuples.
    """
    n = len(species_list)
    feat_keys  = list(FEATURE_LABELS.keys())

    # Build pair lists
    pair_eta: list[float] = []
    pair_feat: dict[str, list[float]] = {f: [] for f in feat_keys}
    pair_native_s: list[float] = []

    for i, si in enumerate(species_list):
        for j, sj in enumerate(species_list):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            pair_eta.append(float(eta_matrix[i, j]))
            pair_native_s.append(EXP19_NATIVE_S[sj])
            for f in feat_keys:
                vi = features[si][f]; vj = features[sj][f]
                if np.isnan(vi) or np.isnan(vj):
                    pair_feat[f].append(np.nan)
                else:
                    pair_feat[f].append(abs(vi - vj))

    results: dict[str, tuple] = {}

    # Spearman for each |Δfeature|
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    rs_vals  = []; feat_names = []
    for f in feat_keys:
        vals = np.array(pair_feat[f])
        eta  = np.array(pair_eta)
        ok   = ~np.isnan(vals) & ~np.isnan(eta)
        if ok.sum() < 4:
            results[f] = (float("nan"), float("nan"))
            continue
        rs, p = stats.spearmanr(vals[ok], eta[ok])
        results[f] = (float(rs), float(p))
        rs_vals.append(rs); feat_names.append(FEATURE_LABELS[f].replace("\n", " "))

    # Panel 1: bar of rs values
    ax = axes[0]
    colours = ["#e74c3c" if r < 0 else "#27ae60" for r in rs_vals]
    bars = ax.barh(feat_names, rs_vals, color=colours, alpha=0.8, edgecolor="black",
                   linewidth=0.7)
    ax.axvline(0, color="black", linewidth=1)
    for bar, (f, (rs, p)) in zip(bars, [(fn, results[fn]) for fn in FEATURE_LABELS
                                         if fn in results and not np.isnan(results[fn][0])]):
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        ax.text(bar.get_width() + (0.01 if rs >= 0 else -0.01),
                bar.get_y() + bar.get_height() / 2,
                f"{rs:+.3f}{sig}", va="center", ha="left" if rs >= 0 else "right",
                fontsize=9)
    ax.set_xlabel("Spearman rs(|Δfeature|, η)", fontsize=11)
    ax.set_title("Which feature distance predicts η?\n(red=negative, green=positive)",
                 fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    ax.set_xlim(-0.7, 0.9)

    # Panel 2: native_S(target) vs η
    ax2 = axes[1]
    eta_arr = np.array(pair_eta)
    nat_arr = np.array(pair_native_s)
    rs_nat, p_nat = stats.spearmanr(nat_arr, eta_arr)
    results["native_S_target"] = (float(rs_nat), float(p_nat))
    sp_colours = []
    sp_idx_j = []
    idx = 0
    for i in range(n):
        for j in range(n):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            sp_colours.append(SPECIES_COLOURS[species_list[j]])
            idx += 1

    ax2.scatter(nat_arr, eta_arr, c=sp_colours, s=60, alpha=0.85,
                edgecolor="black", linewidth=0.5)
    slope, intercept, *_ = stats.linregress(nat_arr, eta_arr)
    xs = np.linspace(nat_arr.min(), nat_arr.max(), 100)
    ax2.plot(xs, slope * xs + intercept, "k--", linewidth=1.5, alpha=0.7)
    ax2.set_xlabel("native_S(target species)", fontsize=11)
    ax2.set_ylabel("η", fontsize=11)
    ax2.set_title(
        f"Target flatness vs η\n"
        f"Spearman rs = {rs_nat:+.3f}  p = {p_nat:.4f}  "
        f"H20c: {'✓ |rs|>0.5' if abs(rs_nat) > 0.5 else '✗ |rs|≤0.5'}",
        fontsize=11,
    )
    # Legend patches
    patches = [mpatches.Patch(color=SPECIES_COLOURS[sp], label=sp)
               for sp in species_list]
    ax2.legend(handles=patches, fontsize=8, loc="lower right")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = str(RESULTS_DIR / "20_feature_correlation_matrix.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 20_feature_correlation_matrix.png")
    return results


def fig_hurst_eta(species_list: list[str], features: dict[str, dict],
                  eta_matrix: np.ndarray, rs: float, p: float) -> None:
    n = len(species_list)
    delta_h: list[float] = []
    eta_vals: list[float] = []
    pair_labels: list[str] = []
    pair_colours: list[str] = []

    for i, si in enumerate(species_list):
        for j, sj in enumerate(species_list):
            if i == j or np.isnan(eta_matrix[i, j]):
                continue
            hi = features[si]["hurst"]; hj = features[sj]["hurst"]
            if not (np.isnan(hi) or np.isnan(hj)):
                delta_h.append(abs(hi - hj))
                eta_vals.append(float(eta_matrix[i, j]))
                pair_labels.append(f"{si[:4]}→{sj[:4]}")
                pair_colours.append(SPECIES_COLOURS[si])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(delta_h, eta_vals, c=pair_colours, s=80, alpha=0.85,
               edgecolor="black", linewidth=0.6, zorder=5)

    slope, intercept, *_ = stats.linregress(delta_h, eta_vals)
    xs = np.linspace(min(delta_h), max(delta_h), 100)
    ax.plot(xs, slope * xs + intercept, "k--", linewidth=1.5, alpha=0.7,
            label=f"rs = {rs:+.3f}  p = {p:.4f}")

    for i, lbl in enumerate(pair_labels):
        ax.annotate(lbl, (delta_h[i], eta_vals[i]),
                    xytext=(3, 2), textcoords="offset points", fontsize=7, alpha=0.8)

    ax.set_xlabel("|ΔHurst H|  (DFA)", fontsize=12)
    ax.set_ylabel("η = native_S / transfer_S", fontsize=12)
    ax.set_title(
        f"|ΔHurst| vs JEPA Transfer Efficiency\n"
        f"Spearman rs = {rs:+.3f}  p = {p:.4f}  "
        f"H20a: {'✓ strongest predictor + p<0.05' if p < 0.05 and rs < -0.3 else '✗'}",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "20_hurst_eta_scatter.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 20_hurst_eta_scatter.png")


def fig_spotlight(features: dict[str, dict]) -> None:
    """
    Side-by-side feature profile for the top pair (Omphalotus, Cordyceps)
    and worst pair (Enoki, Schizophyllum).
    """
    feat_keys = list(FEATURE_LABELS.keys())
    pairs = [
        ("Omphalotus", "Cordyceps", 0.7904, "#27ae60", "Top pair (η=0.790)"),
        ("Enoki",      "Schizophyllum", 0.2379, "#e74c3c", "Worst pair (η=0.238)"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)

    for ax, (sp1, sp2, eta, colour, title) in zip(axes, pairs):
        f1 = [features[sp1][f] for f in feat_keys]
        f2 = [features[sp2][f] for f in feat_keys]
        labels = [FEATURE_LABELS[f].replace("\n", " ") for f in feat_keys]

        x  = np.arange(len(feat_keys))
        w  = 0.35
        ax.bar(x - w / 2, f1, w, color=SPECIES_COLOURS[sp1], alpha=0.8,
               edgecolor="black", linewidth=0.6, label=sp1)
        ax.bar(x + w / 2, f2, w, color=SPECIES_COLOURS[sp2], alpha=0.8,
               edgecolor="black", linewidth=0.6, label=sp2)

        # Annotate |Δ|
        for xi, (v1, v2) in enumerate(zip(f1, f2)):
            if not (np.isnan(v1) or np.isnan(v2)):
                delta = abs(v1 - v2)
                ymax  = max(abs(v1), abs(v2))
                ax.text(xi, ymax + 0.05 * abs(ymax + 1e-9), f"|Δ|={delta:.3f}",
                        ha="center", fontsize=7, color="black", alpha=0.7)

        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
        ax.set_title(f"{title}\n{sp1} vs {sp2}", fontsize=11)
        ax.legend(fontsize=9)
        ax.axhline(0, color="black", linewidth=0.7, alpha=0.5)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylabel("Feature value (raw)")

    plt.suptitle(
        "Feature Profile: Top vs Worst Transfer Pair\nDo similar Hurst / spectral entropy → high η?",
        fontsize=12, y=1.02,
    )
    plt.tight_layout()
    out = str(RESULTS_DIR / "20_omphalotus_cordyceps_spotlight.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 20_omphalotus_cordyceps_spotlight.png")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EXPERIMENT 20 — Temporal Grammar Fingerprint")
    print("=" * 70)
    print("\nExp 19 η matrix loaded (20 pairs).")
    print(f"Golden pair: Omphalotus→Cordyceps η=0.7904, Cordyceps→Omphalotus η=0.7263")

    # ── 1. Load density series ────────────────────────────────────────────
    print("\n--- Step 1: Load density series ---")
    rho: dict[str, np.ndarray] = {}
    rho["Pleurotus"] = load_pleurotus_rho()
    for sp in MULTISPECIES_FILES:
        rho[sp] = load_multispecies_rho(sp)
    species = list(rho.keys())
    print(f"\n  Loaded {len(species)} species")

    # ── 2. Compute temporal fingerprints ────────────────────────────────
    print("\n--- Step 2: Temporal grammar fingerprints ---")
    features: dict[str, dict[str, float]] = {}
    for sp in species:
        features[sp] = temporal_features(rho[sp], sp)

    # ── 3. Hypothesis tests ──────────────────────────────────────────────
    print("\n--- Step 3: Hypothesis tests ---")
    n = len(species)

    # Collect pair lists for all features
    feat_keys = list(FEATURE_LABELS.keys())
    pair_data: dict[str, list[float]] = {f: [] for f in feat_keys}
    pair_native_s: list[float] = []
    pair_eta:  list[float] = []
    pair_names:list[str]  = []

    for i, si in enumerate(species):
        for j, sj in enumerate(species):
            if i == j or np.isnan(EXP19_ETA[i, j]):
                continue
            pair_eta.append(float(EXP19_ETA[i, j]))
            pair_native_s.append(EXP19_NATIVE_S[sj])
            pair_names.append(f"{si}→{sj}")
            for f in feat_keys:
                vi = features[si][f]; vj = features[sj][f]
                pair_data[f].append(float("nan") if (np.isnan(vi) or np.isnan(vj))
                                    else abs(vi - vj))

    eta_arr = np.array(pair_eta)

    feat_results: dict[str, tuple] = {}
    print("\n  Feature Spearman rs with η (20 directed pairs):")
    for f in feat_keys:
        vals = np.array(pair_data[f])
        ok   = ~np.isnan(vals)
        if ok.sum() < 4:
            feat_results[f] = (float("nan"), float("nan"))
            print(f"    {f:12s}  rs=N/A (insufficient data)")
            continue
        rs, p = stats.spearmanr(vals[ok], eta_arr[ok])
        feat_results[f] = (float(rs), float(p))
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        print(f"    {f:12s}  rs={rs:+.4f}  p={p:.4f}  {sig}")

    # H20a: Is Hurst the best predictor?
    valid_feats = {f: feat_results[f] for f in feat_keys
                   if not np.isnan(feat_results[f][0])}
    best_feat = min(valid_feats.items(), key=lambda kv: kv[1][0])
    h_rs, h_p = feat_results.get("hurst", (float("nan"), float("nan")))
    h20a = (not np.isnan(h_rs)) and (best_feat[0] == "hurst") and (h_p < 0.05)
    print(f"\n  H20a — Hurst distance is best predictor (rs most negative, p<0.05):")
    print(f"    Best predictor: {best_feat[0]}  rs={best_feat[1][0]:+.4f}  p={best_feat[1][1]:.4f}")
    print(f"    Hurst: rs={h_rs:+.4f}  p={h_p:.4f}")
    print(f"    → {'✓ CONFIRMED' if h20a else '✗ NOT CONFIRMED'}")

    # H20b: Omphalotus–Cordyceps have smallest |ΔHurst| within smooth tier
    smooth_pairs: list[tuple[str, str, float]] = []
    for i, si in enumerate(species):
        for j, sj in enumerate(species):
            if i == j or si not in TIER_SMOOTH or sj not in TIER_SMOOTH:
                continue
            hi = features[si]["hurst"]; hj = features[sj]["hurst"]
            if not (np.isnan(hi) or np.isnan(hj)):
                smooth_pairs.append((si, sj, abs(hi - hj)))

    smooth_pairs.sort(key=lambda x: x[2])
    print(f"\n  H20b — Omphalotus–Cordyceps have smallest |ΔH| within smooth tier:")
    for si, sj, dh in smooth_pairs[:5]:
        marker = "← ★" if {si, sj} == {"Omphalotus", "Cordyceps"} else ""
        print(f"    {si:15s} ↔ {sj:15s}  |ΔH|={dh:.4f}  {marker}")
    h20b_met = any({si, sj} == {"Omphalotus", "Cordyceps"}
                   for si, sj, _ in smooth_pairs[:2])
    print(f"    → {'✓ CONFIRMED (top-2)' if h20b_met else '✗ NOT CONFIRMED'}")

    # H20c: native_S(target) explains > 50% of η variance
    rs_nat, p_nat = stats.spearmanr(pair_native_s, pair_eta)
    h20c = abs(rs_nat) > 0.5
    print(f"\n  H20c — native_S(target) explains >50% η variance (|rs|>0.5):")
    print(f"    Spearman rs={rs_nat:+.4f}  p={p_nat:.4f}")
    print(f"    → {'✓ CONFIRMED' if h20c else '✗ NOT CONFIRMED'}")

    # H20d: composite distance beats individual features
    # Build composite: Euclidean over all normalised features
    feat_matrix = np.array([[features[sp][f] for f in feat_keys] for sp in species])
    mn = feat_matrix.min(axis=0); mx = feat_matrix.max(axis=0)
    rng = np.where(mx - mn > 0, mx - mn, 1.0)
    norm = (feat_matrix - mn) / rng

    composite_d: list[float] = []
    for i, si in enumerate(species):
        for j, sj in enumerate(species):
            if i == j or np.isnan(EXP19_ETA[i, j]):
                continue
            d = float(np.linalg.norm(norm[i] - norm[j]))
            composite_d.append(d)

    rs_comp, p_comp = stats.spearmanr(composite_d, pair_eta)
    best_individual_rs = min(v[0] for v in valid_feats.values())  # most negative
    h20d = abs(rs_comp) > abs(best_individual_rs)
    print(f"\n  H20d — Composite distance (L2 all features) beats best individual:")
    print(f"    Composite: rs={rs_comp:+.4f}  p={p_comp:.4f}")
    print(f"    Best individual: {best_feat[0]} rs={best_feat[1][0]:+.4f}")
    print(f"    → {'✓ CONFIRMED' if h20d else '✗ NOT CONFIRMED'}")

    # ── 4. Figures ────────────────────────────────────────────────────────
    print("\n--- Step 4: Figures ---")
    fig_radar(features)
    corr_results = fig_feature_correlation(species, features, EXP19_ETA)
    fig_hurst_eta(species, features, EXP19_ETA, h_rs, h_p)
    fig_spotlight(features)

    # ── 5. Print fingerprint table ────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 20 SUMMARY — Temporal Grammar Fingerprints")
    print("=" * 70)
    header = f"  {'Species':15s}  {'AR1':>6}  {'CV':>6}  {'Hurst':>6}  {'S_ent':>6}  {'Skew':>7}  {'Kurt':>7}  Tier"
    print(header)
    for sp in species:
        f = features[sp]
        tier = "bursty" if sp in TIER_BURSTY else "smooth"
        print(f"  {sp:15s}  {f['ar1']:6.3f}  {f['cv']:6.3f}  {f['hurst']:6.3f}  "
              f"{f['spec_ent']:6.3f}  {f['skewness']:7.3f}  {f['kurtosis']:7.3f}  {tier}")

    print(f"\n  H20a: Hurst best predictor  {'✓' if h20a else '✗'}")
    print(f"  H20b: Omphalotus–Cordyceps min |ΔH| in smooth tier  {'✓' if h20b_met else '✗'}")
    print(f"  H20c: native_S(target) rs={rs_nat:+.3f} |rs|>0.5  {'✓' if h20c else '✗'}")
    print(f"  H20d: Composite beats individual  {'✓' if h20d else '✗'}")

    # ── 6. Save report ────────────────────────────────────────────────────
    report = {
        "fingerprints": {sp: {k: round(v, 5) if not np.isnan(v) else None
                              for k, v in features[sp].items()}
                         for sp in species},
        "feature_spearman_vs_eta": {f: {"rs": round(rs, 5), "p": round(p, 5)}
                                    for f, (rs, p) in feat_results.items()
                                    if not np.isnan(rs)},
        "native_s_spearman_vs_eta": {"rs": round(float(rs_nat), 5),
                                     "p": round(float(p_nat), 5)},
        "composite_spearman_vs_eta": {"rs": round(float(rs_comp), 5),
                                      "p": round(float(p_comp), 5)},
        "hypotheses": {
            "H20a": {"result": "confirmed" if h20a else "not_confirmed",
                     "best_predictor": best_feat[0],
                     "hurst_rs": round(float(h_rs), 4) if not np.isnan(h_rs) else None,
                     "hurst_p":  round(float(h_p),  4) if not np.isnan(h_p)  else None},
            "H20b": {"result": "confirmed" if h20b_met else "not_confirmed",
                     "top_smooth_pairs": [{"src": si, "tgt": sj, "delta_hurst": round(dh, 5)}
                                          for si, sj, dh in smooth_pairs[:5]]},
            "H20c": {"result": "confirmed" if h20c else "not_confirmed",
                     "rs_native_s": round(float(rs_nat), 4),
                     "p": round(float(p_nat), 4)},
            "H20d": {"result": "confirmed" if h20d else "not_confirmed",
                     "composite_rs": round(float(rs_comp), 4),
                     "composite_p": round(float(p_comp), 4)},
        },
    }
    out = str(RESULTS_DIR / "20_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report → 20_report.json")
    print("Done.  Figures → experiments/results/20_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
