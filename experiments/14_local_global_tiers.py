"""
Experiment 14 — Two-Tier Communication: Local vs Global Signal Tiers
=====================================================================

Hypothesis: fungal electrical activity operates on two distinct communication tiers.

  Global tier  — spike/burst events (|v_detrended| >= 0.5 mV)
                 High-amplitude, network-propagating events. Previously counted
                 as the *only* meaningful signal.

  Local tier   — sub-threshold events (noise_floor < |v_detrended| < 0.5 mV)
                 Low-amplitude, discarded in prior experiments as "background".
                 UKFT hypothesis: this is the pilot-wave / pre-collapse activity,
                 carrying local coordination before a choice collapse (spike) fires.

UKFT framing
  - Global tier = choice collapse at high information density rho
  - Local tier  = guidance field (pilot wave), sub-threshold rho fluctuations
  - Schizophyllum (98.5% silent) predicted to operate mostly in local tier
  - Omphalotus  (0.03% silent)   predicted to broadcast globally most of the time

Four hypotheses
  H14a  Local event density has a Zipf / power-law distribution
        (different exponent from global, indicating a different generative process)
  H14b  rho_local peaks PRECEDE global burst onsets (lag tau > 0 in cross-correlation)
        = Granger-style evidence that local activity primes global broadcast
  H14c  Spatial correlation matrix for local tier differs from global tier
        (local activity has a different network topology)
  H14d  Schizophyllum rho_local/rho_global  >  Omphalotus rho_local/rho_global
        (species operating in silence still communicate locally)

Data
  tools/data/multispecies/*.txt   Raw continuous differential electrode recordings
  Format: tab-separated, 1 header row, 7 differential channels, NaN possible
  Sample rate: ~1 sample/second

Outputs (experiments/results/)
  14_tier_distribution.png        Stacked bar: noise / local / global fractions × 4 species
  14_ratio_per_species.png        rho_local / rho_global ratio per species
  14_crosscorr_lag.png            Cross-correlation of rho_local(t) vs rho_global(t+tau)
  14_spatial_topology.png         Local vs global spatial correlation matrices
  14_report.json                  All computed metrics
"""

from __future__ import annotations
import os
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import uniform_filter1d
from scipy.stats import pearsonr, linregress

# ── Paths ─────────────────────────────────────────────────────────────────────
# Convention (matches Exp 11): this script lives in ukftbio/experiments/
# Data lives in sibling repo noosphere/apps/myco-explorer/tools/data/
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = str(MYCO / "tools/data/multispecies")
RESULTS_DIR = str(Path(__file__).resolve().parent / "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
GLOBAL_THRESH_MV  = 0.5         # |v_detrended| >= this → global event
NOISE_FLOOR_K     = 2.0         # noise: |v_detrended| < K * sigma_noise
WINDOW_S          = 600         # 10-minute analysis window (samples @ 1 Hz)
DETREND_WINDOW    = 600         # rolling median window for baseline removal
MAX_LAG_WINDOWS   = 300         # cross-correlation lag range (in windows)
POWER_LAW_BINS    = 30          # bins for density histogram / Zipf fit

SPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
    "Enoki":         "Enoki fungi Flammulina velutipes.txt",
    "Omphalotus":    "Ghost Fungi Omphalotus nidiformis.txt",
}

# Colour palette per species (consistent with prior experiments)
SPECIES_COLOURS = {
    "Schizophyllum": "#9b59b6",   # purple  (silent phenotype)
    "Cordyceps":     "#e67e22",   # orange
    "Enoki":         "#27ae60",   # green
    "Omphalotus":    "#2980b9",   # blue    (active phenotype)
}

TIER_COLOURS = {
    "noise":  "#bdc3c7",
    "local":  "#f39c12",
    "global": "#e74c3c",
}


# ─────────────────────────────────────────────────────────────────────────────
# LOADING & PREPROCESSING
# ─────────────────────────────────────────────────────────────────────────────

def load_species(name: str) -> np.ndarray:
    """Load raw voltage array (n_samples, n_channels) for one species.

    Returns float32 array with NaN rows dropped.
    Columns: 7 differential channels.
    """
    path = os.path.join(DATA_DIR, SPECIES_FILES[name])
    print(f"  Loading {name} from {os.path.basename(path)} …", flush=True)

    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        dtype=np.float32,
        na_values=["NaN", "", " ", "nan"],
        engine="c",
        on_bad_lines="skip",
    )

    # Drop rows where ALL channels are NaN
    df.dropna(how="all", inplace=True)

    # Keep only numeric columns
    df = df.select_dtypes(include=[np.floating, np.integer])

    arr = df.values.astype(np.float32)
    print(f"    → {arr.shape[0]:,} rows × {arr.shape[1]} channels", flush=True)
    return arr


def rolling_median_detrend(data: np.ndarray, window: int = DETREND_WINDOW) -> np.ndarray:
    """Subtract rolling median baseline from each channel.

    Uses a fast approach: smooth with a uniform filter as an approximation
    for the low-frequency trend, then subtract.  For biological signals at
    1 Hz, a 600-sample (10-min) window captures slow electrode drift.

    For memory-efficiency we operate in float32 throughout.
    """
    n, nchan = data.shape
    detrended = np.empty_like(data)
    for ch in range(nchan):
        col = data[:, ch].astype(np.float64)
        # Replace NaN with local interpolation before filtering
        nan_mask = np.isnan(col)
        if nan_mask.any():
            idxs = np.arange(n)
            col[nan_mask] = np.interp(idxs[nan_mask], idxs[~nan_mask], col[~nan_mask])
        trend = uniform_filter1d(col, size=window, mode="mirror")
        detrended[:, ch] = (col - trend).astype(np.float32)
    return detrended


def estimate_noise_floor(detrended: np.ndarray) -> np.ndarray:
    """Estimate per-channel noise floor sigma using Median Absolute Deviation.

    Returns array of shape (n_channels,).
    """
    med = np.nanmedian(np.abs(detrended), axis=0)
    # MAD → sigma conversion (assuming Gaussian noise): sigma = MAD / 0.6745
    return med / 0.6745


def classify_tiers(
    detrended: np.ndarray,
    sigma_noise: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return boolean masks for each tier.

    noise  = |v| < NOISE_FLOOR_K * sigma_noise
    local  = NOISE_FLOOR_K * sigma_noise <= |v| < GLOBAL_THRESH_MV
    global = |v| >= GLOBAL_THRESH_MV

    Returns dict with keys 'noise', 'local', 'global', each (n_samples, n_channels) bool.
    """
    abs_v = np.abs(detrended)
    floor = (NOISE_FLOOR_K * sigma_noise).astype(np.float32)  # shape (n_channels,)

    noise_mask  = abs_v < floor[np.newaxis, :]
    global_mask = abs_v >= GLOBAL_THRESH_MV
    local_mask  = (~noise_mask) & (~global_mask)

    return {"noise": noise_mask, "local": local_mask, "global": global_mask}


# ─────────────────────────────────────────────────────────────────────────────
# METRICS
# ─────────────────────────────────────────────────────────────────────────────

def tier_fractions(masks: dict[str, np.ndarray]) -> dict[str, float]:
    """Fraction of all (sample × channel) observations in each tier."""
    total = masks["noise"].size
    return {
        tier: float(np.sum(m)) / total
        for tier, m in masks.items()
    }


def windowed_density(mask: np.ndarray, window: int = WINDOW_S) -> np.ndarray:
    """Average per-window event density (fraction of active samples).

    mask: (n_samples, n_channels) bool
    Returns: (n_windows,) float — mean across channels, fraction in window
    """
    n = mask.shape[0]
    n_windows = n // window
    # Reshape to (n_windows, window, n_channels), mean over window × channels
    trimmed = mask[: n_windows * window]
    shaped  = trimmed.reshape(n_windows, window, -1)
    return shaped.mean(axis=(1, 2))


def cross_correlation_lag(
    x: np.ndarray,
    y: np.ndarray,
    max_lag: int = MAX_LAG_WINDOWS,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalised cross-correlation of x (local) vs y (global) at lags.

    Positive lag: x leads y (local precedes global).
    Returns (lags, xcorr_values).
    """
    n = min(len(x), len(y))
    x = x[:n] - x[:n].mean()
    y = y[:n] - y[:n].mean()

    norm = (np.std(x) * np.std(y) * n)
    if norm < 1e-12:
        lags = np.arange(-max_lag, max_lag + 1)
        return lags, np.zeros_like(lags, dtype=float)

    lags = np.arange(-max_lag, max_lag + 1)
    xcorr = np.array([
        np.dot(x[:n - abs(lag)], y[abs(lag):] if lag >= 0 else y[:n - abs(lag)])
        if lag >= 0
        else np.dot(x[abs(lag):], y[:n - abs(lag)])
        for lag in lags
    ]) / norm

    return lags, xcorr


def spatial_correlation_matrix(masks: dict[str, np.ndarray], tier: str) -> np.ndarray:
    """Compute channel × channel Pearson correlation for a given tier.

    For each channel, the signal is the per-WINDOW_S density time series.
    """
    mask = masks[tier]
    n = mask.shape[0]
    n_windows = n // WINDOW_S
    trimmed = mask[: n_windows * WINDOW_S]
    shaped  = trimmed.reshape(n_windows, WINDOW_S, -1).mean(axis=1)  # (n_windows, n_chan)
    n_chan  = shaped.shape[1]

    corr = np.full((n_chan, n_chan), np.nan)
    for i in range(n_chan):
        for j in range(n_chan):
            if np.std(shaped[:, i]) < 1e-10 or np.std(shaped[:, j]) < 1e-10:
                corr[i, j] = 0.0
            else:
                corr[i, j], _ = pearsonr(shaped[:, i], shaped[:, j])
    return corr


def fit_power_law(density: np.ndarray) -> tuple[float, float]:
    """Fit log-log slope to the density distribution (H14a test).

    Returns (slope, r_squared).  Slope < -1 is consistent with Zipf.
    """
    pos = density[density > 0]
    if len(pos) < 10:
        return 0.0, 0.0
    counts, edges = np.histogram(pos, bins=POWER_LAW_BINS)
    centres = 0.5 * (edges[1:] + edges[:-1])
    mask = counts > 0
    log_x = np.log10(centres[mask])
    log_y = np.log10(counts[mask].astype(float))
    slope, intercept, r, _, _ = linregress(log_x, log_y)
    return float(slope), float(r ** 2)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS LOOP
# ─────────────────────────────────────────────────────────────────────────────

def analyse_species(name: str) -> dict:
    print(f"\n[{name}]", flush=True)
    raw = load_species(name)

    print("  Detrending …", flush=True)
    detrended = rolling_median_detrend(raw)

    print("  Classifying tiers …", flush=True)
    sigma = estimate_noise_floor(detrended)
    masks = classify_tiers(detrended, sigma)

    fractions = tier_fractions(masks)
    print(f"  Fractions → noise={fractions['noise']:.3f}  "
          f"local={fractions['local']:.3f}  global={fractions['global']:.3f}", flush=True)

    # Windowed densities (for cross-correlation and ratio)
    rho_local  = windowed_density(masks["local"])
    rho_global = windowed_density(masks["global"])

    mean_local  = float(rho_local.mean())
    mean_global = float(rho_global.mean())
    ratio = mean_local / mean_global if mean_global > 1e-10 else float("inf")
    print(f"  rho_local={mean_local:.4f}  rho_global={mean_global:.4f}  ratio={ratio:.3f}", flush=True)

    # Cross-correlation
    print("  Computing cross-correlation …", flush=True)
    lags, xcorr = cross_correlation_lag(rho_local, rho_global)
    peak_lag = int(lags[np.argmax(xcorr)])
    peak_val = float(np.max(xcorr))
    print(f"  XCorr peak: lag={peak_lag} windows  value={peak_val:.4f}", flush=True)

    # Spatial topology
    print("  Computing spatial correlations …", flush=True)
    corr_local  = spatial_correlation_matrix(masks, "local")
    corr_global = spatial_correlation_matrix(masks, "global")

    # Off-diagonal means (connectivity proxy)
    n_chan = corr_local.shape[0]
    off_diag = np.triu(np.ones((n_chan, n_chan), dtype=bool), k=1)
    mean_conn_local  = float(np.nanmean(np.abs(corr_local[off_diag])))
    mean_conn_global = float(np.nanmean(np.abs(corr_global[off_diag])))
    print(f"  Mean |r| local={mean_conn_local:.3f}  global={mean_conn_global:.3f}", flush=True)

    # Power-law fit (H14a)
    slope_local,  r2_local  = fit_power_law(rho_local)
    slope_global, r2_global = fit_power_law(rho_global)
    print(f"  Power-law slope: local={slope_local:.2f} (R²={r2_local:.2f})  "
          f"global={slope_global:.2f} (R²={r2_global:.2f})", flush=True)

    return {
        "name":          name,
        "n_samples":     int(raw.shape[0]),
        "n_channels":    int(raw.shape[1]),
        "sigma_noise_mean": float(sigma.mean()),
        "fractions":     fractions,
        "rho_local_mean":  mean_local,
        "rho_global_mean": mean_global,
        "ratio_local_global": ratio,
        "xcorr_peak_lag":    peak_lag,
        "xcorr_peak_val":    peak_val,
        "lags":          lags.tolist(),
        "xcorr":         xcorr.tolist(),
        "rho_local_series":  rho_local.tolist(),
        "rho_global_series": rho_global.tolist(),
        "mean_conn_local":   mean_conn_local,
        "mean_conn_global":  mean_conn_global,
        "corr_local":    corr_local.tolist(),
        "corr_global":   corr_global.tolist(),
        "powerlaw_local":  {"slope": slope_local,  "r2": r2_local},
        "powerlaw_global": {"slope": slope_global, "r2": r2_global},
    }


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def fig_tier_distribution(results: list[dict]) -> None:
    """Fig A: Stacked bar chart of noise/local/global fractions."""
    species_names = [r["name"] for r in results]
    noise_fracs  = [r["fractions"]["noise"]  for r in results]
    local_fracs  = [r["fractions"]["local"]  for r in results]
    global_fracs = [r["fractions"]["global"] for r in results]

    x = np.arange(len(species_names))
    fig, ax = plt.subplots(figsize=(8, 5))

    b1 = ax.bar(x, noise_fracs,  label="Noise (<2σ)",       color=TIER_COLOURS["noise"])
    b2 = ax.bar(x, local_fracs,  bottom=noise_fracs,         label="Local (sub-threshold)",
                color=TIER_COLOURS["local"])
    b3 = ax.bar(x, global_fracs, bottom=[n+l for n, l in zip(noise_fracs, local_fracs)],
                label="Global (spike ≥0.5 mV)", color=TIER_COLOURS["global"])

    ax.set_xticks(x)
    ax.set_xticklabels(species_names, fontsize=11)
    ax.set_ylabel("Fraction of observations", fontsize=11)
    ax.set_title("Two-Tier Signal Distribution Across Species\n"
                 "(after rolling-median detrend, per-channel adaptive noise floor)",
                 fontsize=12)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))

    for i, (nl, ll, gl) in enumerate(zip(noise_fracs, local_fracs, global_fracs)):
        if ll > 0.01:
            ax.text(i, nl + ll / 2, f"{ll:.1%}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        if gl > 0.01:
            ax.text(i, nl + ll + gl / 2, f"{gl:.1%}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")

    plt.tight_layout()
    out = os.path.join(RESULTS_DIR, "14_tier_distribution.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {os.path.basename(out)}", flush=True)


def fig_ratio_per_species(results: list[dict]) -> None:
    """Fig B: rho_local / rho_global ratio per species."""
    species_names = [r["name"] for r in results]
    ratios = [r["ratio_local_global"] for r in results]
    colours = [SPECIES_COLOURS[n] for n in species_names]

    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(species_names, ratios, color=colours, edgecolor="white", linewidth=0.5)

    for bar, val in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{val:.2f}x", ha="center", va="bottom", fontsize=10)

    ax.set_ylabel("ρ_local / ρ_global", fontsize=11)
    ax.set_title("Local-to-Global Communication Ratio per Species\n"
                 "H14d: Schizophyllum (silent) predicted high; Omphalotus (active) predicted low",
                 fontsize=11)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="ratio = 1 (parity)")
    ax.legend(fontsize=9)

    plt.tight_layout()
    out = os.path.join(RESULTS_DIR, "14_ratio_per_species.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {os.path.basename(out)}", flush=True)


def fig_crosscorr_lag(results: list[dict]) -> None:
    """Fig C: Cross-correlation of rho_local(t) vs rho_global(t+lag), all species."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ax, r in zip(axes_flat, results):
        lags   = np.array(r["lags"])
        xcorr  = np.array(r["xcorr"])
        colour = SPECIES_COLOURS[r["name"]]

        ax.plot(lags, xcorr, color=colour, linewidth=1.5, alpha=0.9)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
        peak_lag = r["xcorr_peak_lag"]
        peak_val = r["xcorr_peak_val"]
        ax.axvline(peak_lag, color=colour, linestyle=":", linewidth=1.5)
        ax.set_title(f"{r['name']}\npeak lag = {peak_lag} win  r = {peak_val:.3f}",
                     fontsize=10)
        ax.set_xlabel("Lag (10-min windows)", fontsize=9)
        ax.set_ylabel("Norm. cross-corr", fontsize=9)

        direction = "LOCAL LEADS" if peak_lag > 0 else ("GLOBAL LEADS" if peak_lag < 0 else "SIMULTANEOUS")
        colour_dir = "#27ae60" if peak_lag > 0 else "#e74c3c"
        ax.text(0.97, 0.95, direction, transform=ax.transAxes,
                ha="right", va="top", fontsize=8, color=colour_dir,
                fontweight="bold")

    fig.suptitle("Cross-Correlation: ρ_local(t) vs ρ_global(t+lag)\n"
                 "H14b: positive peak lag → local activity precedes global bursts",
                 fontsize=12)
    plt.tight_layout()
    out = os.path.join(RESULTS_DIR, "14_crosscorr_lag.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {os.path.basename(out)}", flush=True)


def fig_spatial_topology(results: list[dict]) -> None:
    """Fig D: Local vs Global spatial correlation matrices (for each species)."""
    fig = plt.figure(figsize=(14, 8))
    gs  = gridspec.GridSpec(2, len(results), figure=fig, hspace=0.4, wspace=0.3)

    for col, r in enumerate(results):
        corr_local  = np.array(r["corr_local"])
        corr_global = np.array(r["corr_global"])
        name = r["name"]

        for row, (mat, label) in enumerate([
            (corr_local,  "Local tier"),
            (corr_global, "Global tier"),
        ]):
            ax = fig.add_subplot(gs[row, col])
            im = ax.imshow(mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
            ax.set_title(f"{name}\n{label}\n"
                         f"mean|r|={np.nanmean(np.abs(mat[np.triu(np.ones_like(mat, dtype=bool), k=1)])):.2f}",
                         fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
            if col == len(results) - 1:
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Spatial Correlation Network: Local Tier vs Global Tier\n"
                 "H14c: different topologies predicted (local = pre-collapse coordination)",
                 fontsize=12)
    out = os.path.join(RESULTS_DIR, "14_spatial_topology.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {os.path.basename(out)}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY & REPORT
# ─────────────────────────────────────────────────────────────────────────────

def print_summary(results: list[dict]) -> None:
    print("\n" + "=" * 70)
    print("EXPERIMENT 14 — TWO-TIER COMMUNICATION SUMMARY")
    print("=" * 70)

    # H14d check
    print("\nH14d: rho_local/rho_global ratio per species")
    sorted_r = sorted(results, key=lambda r: r["ratio_local_global"], reverse=True)
    for r in sorted_r:
        bar = "#" * int(r["ratio_local_global"] * 10)
        print(f"  {r['name']:20s}  {r['ratio_local_global']:6.2f}x  {bar}")

    schiz = next((r for r in results if "Schizophyllum" in r["name"]), None)
    omph  = next((r for r in results if "Omphalotus"    in r["name"]), None)
    if schiz and omph:
        h14d_pass = schiz["ratio_local_global"] > omph["ratio_local_global"]
        print(f"\n  H14d {'CONFIRMED' if h14d_pass else 'NOT CONFIRMED'}: "
              f"Schizophyllum ratio ({schiz['ratio_local_global']:.2f}x) "
              f"{'>' if h14d_pass else '<='} "
              f"Omphalotus ({omph['ratio_local_global']:.2f}x)")

    # H14b check
    print("\nH14b: cross-correlation peak lag (positive = local leads global)")
    for r in results:
        direction = "LOCAL LEADS  ✓" if r["xcorr_peak_lag"] > 0 else (
                    "GLOBAL LEADS ✗" if r["xcorr_peak_lag"] < 0 else "SIMULTANEOUS")
        print(f"  {r['name']:20s}  lag={r['xcorr_peak_lag']:+5d} win  "
              f"r={r['xcorr_peak_val']:.4f}  {direction}")

    # H14a check
    print("\nH14a: power-law slope (local density distribution; slope < -1 → Zipf)")
    for r in results:
        sl = r["powerlaw_local"]["slope"]
        r2 = r["powerlaw_local"]["r2"]
        zipf = sl < -1.0 and r2 > 0.5
        print(f"  {r['name']:20s}  slope={sl:6.2f}  R²={r2:.2f}  "
              f"{'Zipf-like ✓' if zipf else 'not Zipf'}")

    # H14c check
    print("\nH14c: spatial topology (|mean_conn_local| vs |mean_conn_global|)")
    for r in results:
        diff = r["mean_conn_local"] - r["mean_conn_global"]
        print(f"  {r['name']:20s}  local={r['mean_conn_local']:.3f}  "
              f"global={r['mean_conn_global']:.3f}  diff={diff:+.3f}")

    print("\n" + "=" * 70)


def save_report(results: list[dict]) -> None:
    """Save full results to JSON (omitting long series for readability)."""
    report = []
    for r in results:
        rec = {k: v for k, v in r.items()
               if k not in ("lags", "xcorr", "rho_local_series", "rho_global_series",
                            "corr_local", "corr_global")}
        report.append(rec)

    out = os.path.join(RESULTS_DIR, "14_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Report saved → {os.path.basename(out)}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 70)
    print("EXPERIMENT 14 — Two-Tier Communication: Local vs Global Signal Tiers")
    print("=" * 70)
    print(f"\nGlobal threshold : {GLOBAL_THRESH_MV} mV")
    print(f"Noise floor      : {NOISE_FLOOR_K}× σ_noise (per-channel adaptive)")
    print(f"Analysis window  : {WINDOW_S} s (10 min)")
    print(f"Cross-corr range : ±{MAX_LAG_WINDOWS} windows")
    print()

    all_results = []
    for name in SPECIES_FILES:
        try:
            result = analyse_species(name)
            all_results.append(result)
        except FileNotFoundError:
            print(f"  [SKIP] {name} — data file not found", flush=True)
        except Exception as exc:
            print(f"  [ERROR] {name}: {exc}", flush=True)
            raise

    if not all_results:
        print("No species processed — check DATA_DIR path.", flush=True)
        return

    print("\n--- Generating figures ---", flush=True)
    fig_tier_distribution(all_results)
    fig_ratio_per_species(all_results)
    fig_crosscorr_lag(all_results)
    fig_spatial_topology(all_results)

    print_summary(all_results)
    save_report(all_results)

    print("\nDone.  Figures → experiments/results/14_*.png", flush=True)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
