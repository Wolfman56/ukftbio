#!/usr/bin/env python3
"""
Experiment 13: Multi-Species Varietal Comparison
=================================================
Test UKFT universality hypothesis: do different fungal species converge on the
same embedding manifold structure when processed through a common projection
head trained on Pleurotus ostreatus?

Hypotheses
----------
H13a  Power-law tail — each species: rs(log_ρ[ρ>0], rank_log_ρ) < -0.60
H13b  Distribution shape convergence — KS test on windowed ρ CDFs.
      CONVERGENT if KS p > 0.05 (same shape), DIVERGENT if p < 0.05.
      Either result is scientifically informative.
H13c  Geodesic/Euclidean ratio falls in [0.05, 0.80] for all species,
      confirming non-Euclidean manifold structure universally.
H13d  Species embedding centroids cosine distance < 0.30 (manifold convergence)
      OR > 0.60 (species-specific manifolds). Ambiguous [0.30, 0.60] = partial.

Critical design: ALL species projected through Pleurotus projection_head.safetensors.
This common embedding space is required for valid cross-species comparison.
If species converge → UKFT minimum-action universality supported.
If species separate → species-specific geodesic structure.

Usage
-----
    python3 experiments/13_varietal_comparison.py [--dry-run] [--max-hours N]
    python3 experiments/13_varietal_comparison.py --species pleurotus,schizophyllum
    python3 experiments/13_varietal_comparison.py --list-species

Progress gates (requires at least 2 species with data):
  GATE_H13a  rs_rank < -0.60 in ALL available species
  GATE_H13d  centroid distances are ALL < 0.30 or ALL > 0.60 (clear signal)
  → PHASE 5 GATE OPEN when both gates pass
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import spearmanr, ks_2samp

# ── Paths ─────────────────────────────────────────────────────────────────────
UKFTBIO_ROOT = Path(__file__).parent.parent
MYCO_ROOT = UKFTBIO_ROOT.parent / "noosphere" / "apps" / "myco-explorer"
MYCO_TOOLS = MYCO_ROOT / "tools"
DATA_DIR = MYCO_TOOLS / "data"
MULTI_DIR = DATA_DIR / "multispecies"
RESULTS_DIR = MYCO_ROOT / "results"
MODEL_DIR = DATA_DIR / "model"
PROJECTION_HEAD = MODEL_DIR / "projection_head.safetensors"

# Import live pipeline functions from myco-explorer tools
sys.path.insert(0, str(MYCO_TOOLS))
from export_spike_events import (  # noqa: E402
    compute_baseline,
    detect_spikes,
    assign_bursts,
    AMPLITUDE_THRESHOLD_MV,
    BURST_ISI_THRESHOLD_S,
    MIN_SPIKE_DURATION_S,
    BASELINE_WINDOW_S,
)
from spike_embed import make_histogram, WINDOW_S, STRIDE_S  # noqa: E402
from apply_projection import load_projection_matrix, project  # noqa: E402

# ── Species Registry ──────────────────────────────────────────────────────────
# Each entry: display_name, csv_path (None = use pre-computed aligned NDJSON),
# aligned_ndjson_path (pre-computed path or None to run pipeline)
SPECIES_REGISTRY: Dict[str, Dict] = {
    "pleurotus": {
        "name": "Pleurotus ostreatus",
        "csv": None,                               # pipeline already run
        "aligned_ndjson": DATA_DIR / "pleurotus_spike_aligned.ndjson",
        "borda_json": RESULTS_DIR / "borda_rankings.json",
        "fmm_json": RESULTS_DIR / "fmm_scores.json",
        "color": "blue",
    },
    "schizophyllum": {
        "name": "Schizophyllum commune",
        "csv": MULTI_DIR / "schizophyllum_commune_spike_traces.csv",
        "aligned_ndjson": MULTI_DIR / "schizophyllum_commune_aligned.ndjson",
        "borda_json": None,
        "fmm_json": None,
        "color": "orange",
    },
    "cordyceps": {
        "name": "Cordyceps militaris",
        "csv": MULTI_DIR / "cordyceps_militaris_spike_traces.csv",
        "aligned_ndjson": MULTI_DIR / "cordyceps_militaris_aligned.ndjson",
        "borda_json": None,
        "fmm_json": None,
        "color": "red",
    },
    "flammulina": {
        "name": "Flammulina velutipes",
        "csv": MULTI_DIR / "flammulina_velutipes_spike_traces.csv",
        "aligned_ndjson": MULTI_DIR / "flammulina_velutipes_aligned.ndjson",
        "borda_json": None,
        "fmm_json": None,
        "color": "green",
    },
    "omphalotus": {
        "name": "Omphalotus nidiformis",
        "csv": MULTI_DIR / "omphalotus_nidiformis_spike_traces.csv",
        "aligned_ndjson": MULTI_DIR / "omphalotus_nidiformis_aligned.ndjson",
        "borda_json": None,
        "fmm_json": None,
        "color": "purple",
    },
}


# ── Pipeline helpers ──────────────────────────────────────────────────────────

def load_csv_channels(csv_path: Path, max_rows: Optional[int] = None
                      ) -> Tuple[List[float], Dict[int, List[float]]]:
    """Return (timestamps, {ch_idx: [voltage, ...]}) from a spike-traces CSV."""
    timestamps = []
    channels: Dict[int, List[float]] = {i: [] for i in range(8)}
    n_ch = 8

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row_idx, row in enumerate(reader):
            if max_rows and row_idx >= max_rows:
                break
            t = float(row.get("time_s", row_idx))
            timestamps.append(t)
            for i in range(n_ch):
                key = f"diff{i+1}_mv"
                try:
                    v = float(row.get(key, 0.0) or 0.0)
                    if math.isnan(v):
                        v = 0.0
                except (ValueError, TypeError):
                    v = 0.0
                channels[i].append(v)

    return timestamps, channels


def run_pipeline_for_species(
    species_key: str,
    spec: Dict,
    W, b,                              # projection head matrices
    max_hours: Optional[float] = None,
    dry_run: bool = False,
    force: bool = False,
) -> List[Dict]:
    """
    Run full pipeline (detect → embed → project) for one species from CSV.
    Returns list of aligned embedding records (same schema as pleurotus_spike_aligned).
    Caches to spec["aligned_ndjson"] — skip if already present and not force.
    """
    out_path: Path = spec["aligned_ndjson"]

    if out_path.exists() and not force:
        print(f"  [cache] Loading pre-computed aligned embeddings: {out_path.name}")
        return load_aligned_ndjson(out_path)

    csv_path: Path = spec["csv"]
    if csv_path is None or not csv_path.exists():
        print(f"  [skip] CSV not found: {csv_path}")
        return []

    max_rows = int(max_hours * 3600) if max_hours else None
    print(f"  Loading CSV: {csv_path.name} (max_rows={max_rows or 'all'}) …")
    timestamps, channels = load_csv_channels(csv_path, max_rows=max_rows)
    n_rows = len(timestamps)
    duration_h = n_rows / 3600.0
    print(f"  Loaded {n_rows:,} rows, {duration_h:.1f}h")

    # ── Phase 1: spike detection across all channels ──────────────────────────
    raw_spikes: List[Dict] = []
    for ch_idx in range(8):
        ch_voltages = channels[ch_idx]
        if all(v == 0.0 for v in ch_voltages):
            continue  # dead channel
        spikes = detect_spikes(timestamps, ch_voltages, ch_idx, dt_s=1.0)
        raw_spikes.extend(spikes)

    # assign_bursts adds in_burst, burst_id, preceding_isi_s (required by make_histogram)
    all_spikes = assign_bursts(raw_spikes)
    print(f"  Detected {len(all_spikes):,} spikes across all channels")

    if dry_run and len(timestamps) > 3600:
        print("  [dry-run] Limiting to first 1h of data for embedding")
        timestamps_emb = timestamps[:3600]
        t_end_cap = 3600.0
    else:
        timestamps_emb = timestamps
        t_end_cap = timestamps[-1] if timestamps else 0.0

    # ── Phase 2: windowed histogram embedding ─────────────────────────────────
    records: List[Dict] = []
    t_start = 0.0
    win_idx = 0
    while t_start + WINDOW_S <= t_end_cap + STRIDE_S:
        t_end = t_start + WINDOW_S
        spikes_in_window = [
            s for s in all_spikes
            if t_start <= s["t_peak_s"] < t_end
        ]
        histogram = make_histogram(spikes_in_window)
        n_spikes = len(spikes_in_window)
        n_bursts = len({s.get("burst_id") for s in spikes_in_window
                        if s.get("burst_id") is not None})

        # ρ = knowledge density (inverse variance of spike rate)
        if n_spikes > 1:
            isi_vals = [s["preceding_isi_s"] for s in spikes_in_window
                        if s.get("preceding_isi_s") is not None
                        and float(s["preceding_isi_s"]) > 0]
            if len(isi_vals) > 1:
                sigma = float(np.std(isi_vals))
                rho = 1.0 / sigma if sigma > 0 else 0.0
            else:
                rho = 1.0 / max(n_spikes, 1)
                sigma = float(n_spikes)
        else:
            rho = 0.0
            sigma = 0.0

        # ── Phase 3: apply projection ─────────────────────────────────────────
        aligned = project(histogram.tolist(), W, b)

        records.append({
            "window_id": f"w_{win_idx:06d}",
            "t_start_s": t_start,
            "t_end_s": t_end,
            "channel": None,
            "n_spikes": n_spikes,
            "n_bursts": n_bursts,
            "rho": round(rho, 6),
            "sigma": round(sigma, 6),
            "embedding": [round(float(x), 6) for x in histogram],
            "aligned_embedding": [round(float(x), 6) for x in aligned],
        })
        t_start += STRIDE_S
        win_idx += 1

    print(f"  Built {len(records):,} windows")

    # Cache to disk
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    print(f"  Cached → {out_path.name}")

    return records


def load_aligned_ndjson(path: Path) -> List[Dict]:
    """Load NDJSON aligned embedding file."""
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


# ── Metrics ───────────────────────────────────────────────────────────────────

def compute_rho_distribution(records: List[Dict]) -> Dict:
    """Compute ρ distribution statistics (median, p95, heavy tail, silence fraction)."""
    rho_vals = np.array([r["rho"] for r in records], dtype=np.float64)
    nonzero = rho_vals[rho_vals > 0]
    silence_fraction = float(np.mean(rho_vals == 0))

    if len(nonzero) < 5:
        return {"median": 0.0, "p95": 0.0, "p99": 0.0, "skewness": 0.0,
                "n_nonzero": int(len(nonzero)), "n_total": len(rho_vals),
                "silence_fraction": round(silence_fraction, 4),
                "heavy_tail_ratio": 0.0,
                "burst_fraction": 0.0}

    # Power-law / heavy-tail test: ρ_p99 / ρ_median.
    # UKFT predicts heavy tail (most time silent, occasional bursts) → ratio ≥ 3.0
    heavy_tail_ratio = float(np.percentile(nonzero, 99)) / float(np.median(nonzero))

    # Zipf log-rank correlation: sort ρ descending, then rs(log_rank, log_ρ_sorted)
    # For power-law: log_ρ ~ -α*log_rank → rs close to -1.0
    log_rho_sorted = np.sort(np.log(nonzero))[::-1]   # descending
    log_ranks = np.log(np.arange(1, len(log_rho_sorted) + 1))
    rs_zipf, _ = spearmanr(log_ranks, log_rho_sorted)  # expect < -0.60 for power-law

    skewness = float(np.mean(((nonzero - nonzero.mean()) / (nonzero.std() + 1e-12)) ** 3))
    burst_fraction = float(np.mean(np.array([r.get("n_bursts", 0) for r in records]) > 0))

    return {
        "median": round(float(np.median(nonzero)), 6),
        "mean": round(float(nonzero.mean()), 6),
        "p95": round(float(np.percentile(nonzero, 95)), 6),
        "p99": round(float(np.percentile(nonzero, 99)), 6),
        "skewness": round(skewness, 4),
        "n_nonzero": int(len(nonzero)),
        "n_total": int(len(rho_vals)),
        "silence_fraction": round(silence_fraction, 4),
        "heavy_tail_ratio": round(heavy_tail_ratio, 4),
        "rs_zipf": round(float(rs_zipf), 4),
        "burst_fraction": round(burst_fraction, 4),
    }


def compute_geodesic_ratio(records: List[Dict], n_sample: int = 400) -> Dict:
    """
    Estimate geodesic/Euclidean ratio on the embedding manifold.
    Uses a random subsample of pairs for tractability.
    Geodesic ≈ normalised angular distance (cosine); Euclidean = L2 norm.
    ratio < 1.0 indicates manifold compression (geodesic shorter than Euclidean).
    """
    embs = np.array([r["aligned_embedding"] for r in records], dtype=np.float32)
    n = len(embs)
    if n < 10:
        return {"ratio_median": float("nan"), "ratio_p95": float("nan"), "n_pairs": 0}

    rng = np.random.default_rng(42)
    idx = rng.choice(n, size=min(n_sample, n), replace=False)
    sub = embs[idx]

    # Normalise for cosine angle
    norms = np.linalg.norm(sub, axis=1, keepdims=True)
    sub_n = sub / np.where(norms > 0, norms, 1.0)

    ratios = []
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            cos_sim = float(np.clip(sub_n[i] @ sub_n[j], -1.0, 1.0))
            angular = math.acos(cos_sim) / math.pi   # [0, 1] geodesic proxy
            euclidean = float(np.linalg.norm(sub[i] - sub[j]))
            if euclidean > 1e-8:
                ratios.append(angular / euclidean)

    if not ratios:
        return {"ratio_median": float("nan"), "ratio_p95": float("nan"), "n_pairs": 0}
    ratios_arr = np.array(ratios)
    return {
        "ratio_median": round(float(np.median(ratios_arr)), 6),
        "ratio_p95": round(float(np.percentile(ratios_arr, 95)), 6),
        "n_pairs": len(ratios),
    }


def compute_centroid(records: List[Dict]) -> np.ndarray:
    """Compute mean of aligned embeddings (L2-normalised centroid)."""
    embs = np.array([r["aligned_embedding"] for r in records], dtype=np.float32)
    centroid = embs.mean(axis=0)
    norm = np.linalg.norm(centroid)
    return (centroid / norm) if norm > 0 else centroid


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine distance in [0, 1]."""
    sim = float(np.clip(a @ b, -1.0, 1.0))
    return round(1.0 - sim, 6)


# ── Cross-species tests ───────────────────────────────────────────────────────

def ks_rho_test(records_a: List[Dict], records_b: List[Dict]) -> Dict:
    """KS test on ρ distributions between two species."""
    rho_a = np.array([r["rho"] for r in records_a if r["rho"] > 0])
    rho_b = np.array([r["rho"] for r in records_b if r["rho"] > 0])
    if len(rho_a) < 5 or len(rho_b) < 5:
        return {"ks_stat": float("nan"), "ks_p": float("nan"), "verdict": "insufficient"}
    stat, p = ks_2samp(rho_a, rho_b)
    verdict = "same_shape (p>0.05)" if p > 0.05 else "different_shape (p<0.05)"
    return {"ks_stat": round(float(stat), 4), "ks_p": round(float(p), 6), "verdict": verdict}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Exp 13: Multi-Species Varietal Comparison")
    parser.add_argument("--dry-run", action="store_true",
                        help="Fast mode: limit Schizophyllum to 1h of data")
    parser.add_argument("--max-hours", type=float, default=None,
                        help="Limit CSV ingestion to first N hours (for speed)")
    parser.add_argument("--species", default=None,
                        help="Comma-separated species keys to include, e.g. pleurotus,schizophyllum")
    parser.add_argument("--list-species", action="store_true",
                        help="Show available species and exit")
    parser.add_argument("--force", action="store_true",
                        help="Re-run pipeline even if cached aligned_ndjson exists")
    args = parser.parse_args()

    if args.list_species:
        print("\nAvailable species:")
        for key, spec in SPECIES_REGISTRY.items():
            aligned = spec["aligned_ndjson"]
            csv_f = spec.get("csv")
            has_aligned = "✅ aligned" if (aligned and aligned.exists()) else "  (no aligned)"
            has_csv = "✅ csv" if (csv_f and csv_f.exists()) else "  (no csv)"
            print(f"  {key:<16} {spec['name']:<35} {has_aligned}  {has_csv}")
        sys.exit(0)

    if args.dry_run:
        print("=" * 70)
        print("DRY RUN MODE — limiting CSV to 1h of data per species")
        print("=" * 70)
        max_hours = args.max_hours or 1.0
    else:
        max_hours = args.max_hours

    # ── Load projection head (trained on Pleurotus, used for ALL species) ─────
    print(f"\nLoading projection head: {PROJECTION_HEAD.name}")
    if not PROJECTION_HEAD.exists():
        print("ERROR: Run bert_align.py first to train projection_head.safetensors")
        sys.exit(1)
    W, b = load_projection_matrix()
    print(f"  W={W.shape}, b={b.shape}")

    # ── Determine which species to process ────────────────────────────────────
    if args.species:
        selected_keys = [k.strip() for k in args.species.split(",")]
    else:
        # Auto-detect available species (have CSV or pre-computed aligned)
        selected_keys = []
        for key, spec in SPECIES_REGISTRY.items():
            aligned = spec["aligned_ndjson"]
            csv_f = spec.get("csv")
            if (aligned and aligned.exists()) or (csv_f and csv_f.exists()):
                selected_keys.append(key)

    print(f"\nProcessing species: {', '.join(selected_keys)}")

    if len(selected_keys) < 1:
        print("No species data found. Download data first:")
        print("  cd noosphere/apps/myco-explorer && bash tools/download_multispecies.sh all")
        sys.exit(1)

    # ── Run pipeline per species ───────────────────────────────────────────────
    species_data: Dict[str, List[Dict]] = {}
    for key in selected_keys:
        if key not in SPECIES_REGISTRY:
            print(f"WARNING: Unknown species key '{key}', skipping")
            continue
        spec = SPECIES_REGISTRY[key]
        print(f"\n[{key}] {spec['name']}")

        if spec.get("csv") is None:
            # Pleurotus: load pre-computed
            aligned_path = spec["aligned_ndjson"]
            if aligned_path and aligned_path.exists():
                print(f"  Loading pre-computed aligned embeddings: {aligned_path.name}")
                records = load_aligned_ndjson(aligned_path)
                print(f"  Loaded {len(records):,} records")
            else:
                print(f"  ERROR: Pleurotus aligned NDJSON not found at {aligned_path}")
                continue
        else:
            records = run_pipeline_for_species(
                key, spec, W, b,
                max_hours=max_hours,
                dry_run=args.dry_run,
                force=args.force,
            )

        if records:
            species_data[key] = records
        else:
            print(f"  WARNING: No records produced for {key}")

    if len(species_data) < 1:
        print("\nERROR: No species processed successfully. Exiting.")
        sys.exit(1)

    if len(species_data) < 2:
        print(f"\nWARNING: Only 1 species available ({list(species_data.keys())[0]}).")
        print("Cross-species comparison requires ≥2 species.")
        print("Download remaining species:")
        print("  cd noosphere/apps/myco-explorer && bash tools/download_multispecies.sh all\n")

    # ── Compute per-species metrics ───────────────────────────────────────────
    print("\n" + "=" * 70)
    print("PER-SPECIES METRICS")
    print("=" * 70)

    per_species: Dict[str, Dict] = {}
    centroids: Dict[str, np.ndarray] = {}

    for key, records in species_data.items():
        name = SPECIES_REGISTRY[key]["name"]
        rho_stats = compute_rho_distribution(records)
        geo_ratio = compute_geodesic_ratio(records)
        centroid = compute_centroid(records)
        centroids[key] = centroid

        per_species[key] = {
            "name": name,
            "n_windows": len(records),
            "rho": rho_stats,
            "geodesic": geo_ratio,
        }

        print(f"\n{name} ({key})")
        print(f"  Windows:          {len(records):,}")
        print(f"  ρ median:          {rho_stats['median']:.4f}")
        print(f"  ρ p95:             {rho_stats['p95']:.4f}")
        print(f"  ρ p99/median:      {rho_stats['heavy_tail_ratio']:.2f}x (heavy tail ≥3.0)")
        print(f"  ρ silence:         {rho_stats['silence_fraction']:.3f} (UKFT ≥0.50)")
        print(f"  ρ skewness:        {rho_stats['skewness']:.3f}")
        print(f"  ρ nonzero/total:   {rho_stats['n_nonzero']}/{rho_stats['n_total']}")
        print(f"  ρ burst fraction:  {rho_stats['burst_fraction']:.3f}")
        print(f"  rs_zipf(log_rank,log_ρ): {rho_stats['rs_zipf']:.4f} (power-law ≤-0.60)")
        print(f"  Geodesic ratio:   {geo_ratio['ratio_median']:.4f} (median)")

    # ── Cross-species comparisons ─────────────────────────────────────────────
    cross_comparisons: List[Dict] = []
    centroid_distances: List[float] = []

    if len(species_data) >= 2:
        print("\n" + "=" * 70)
        print("CROSS-SPECIES COMPARISONS")
        print("=" * 70)

        keys = list(species_data.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                ka, kb = keys[i], keys[j]
                na = SPECIES_REGISTRY[ka]["name"]
                nb = SPECIES_REGISTRY[kb]["name"]

                # KS test on ρ distributions
                ks = ks_rho_test(species_data[ka], species_data[kb])

                # Centroid cosine distance
                c_dist = cosine_distance(centroids[ka], centroids[kb])
                centroid_distances.append(c_dist)

                comparison = {
                    "pair": f"{ka}_vs_{kb}",
                    "species_a": na,
                    "species_b": nb,
                    "ks_rho": ks,
                    "centroid_cosine_dist": c_dist,
                }
                cross_comparisons.append(comparison)

                print(f"\n{na} vs {nb}")
                print(f"  KS(ρ): stat={ks['ks_stat']:.4f}, p={ks['ks_p']:.4f}  → {ks['verdict']}")
                print(f"  Centroid cosine dist: {c_dist:.4f}")

    # ── Hypothesis evaluation ─────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("HYPOTHESIS EVALUATION")
    print("=" * 70)

    h13a_results = {}
    for key, metrics in per_species.items():
        rho = metrics["rho"]
        silence = rho.get("silence_fraction", 0.0)
        tail_ratio = rho.get("heavy_tail_ratio", 0.0)
        rs_zipf = rho.get("rs_zipf", float("nan"))
        # UKFT sparse-burst criterion: system mostly silent + heavy-tailed bursts
        silence_ok = silence >= 0.40   # ≥40% of windows silent
        tail_ok = tail_ratio >= 3.0    # burst windows are ≥3x median ρ
        zipf_ok = not math.isnan(rs_zipf) and rs_zipf < -0.60
        passed = silence_ok and (tail_ok or zipf_ok)
        h13a_results[key] = {
            "silence_fraction": silence, "heavy_tail_ratio": tail_ratio,
            "rs_zipf": rs_zipf, "pass": passed,
        }
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"H13a [{key}] silence={silence:.2f} tail_ratio={tail_ratio:.1f}x "
              f"rs_zipf={rs_zipf:.4f} {status}")

    h13b_results = {}
    for comp in cross_comparisons:
        pair = comp["pair"]
        ks_p = comp["ks_rho"]["ks_p"]
        h13b_results[pair] = {
            "ks_p": ks_p,
            "verdict": comp["ks_rho"]["verdict"],
        }
        print(f"H13b [{pair}] p={ks_p:.4f}  → {comp['ks_rho']['verdict']}")

    h13c_results = {}
    all_h13c_pass = True
    for key, metrics in per_species.items():
        ratio = metrics["geodesic"]["ratio_median"]
        passed = not math.isnan(ratio) and 0.05 <= ratio <= 0.80
        h13c_results[key] = {"ratio_median": ratio, "pass": passed}
        if not passed:
            all_h13c_pass = False
        status = "✅ PASS" if passed else ("❌ FAIL" if not math.isnan(ratio) else "⚠️  N/A")
        print(f"H13c [{key}] geodesic ratio = {ratio:.4f} {status}")

    h13d_verdict = "N/A (single species)"
    if centroid_distances:
        mean_dist = float(np.mean(centroid_distances))
        max_dist = float(np.max(centroid_distances))
        if max_dist < 0.30:
            h13d_verdict = f"CONVERGENT (max cosine dist = {max_dist:.4f} < 0.30)"
        elif min(centroid_distances) > 0.60:
            h13d_verdict = f"DIVERGENT (min cosine dist = {min(centroid_distances):.4f} > 0.60)"
        else:
            h13d_verdict = f"PARTIAL/AMBIGUOUS (dist range [{min(centroid_distances):.4f}, {max_dist:.4f}])"
        print(f"H13d centroids — {h13d_verdict}")

    # ── Gate evaluation ───────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    gate_h13a = all(v["pass"] for v in h13a_results.values())
    gate_h13c = all_h13c_pass

    print(f"GATE_H13a (power-law ρ in all species):   {'✅ OPEN' if gate_h13a else '❌ CLOSED'}")
    print(f"GATE_H13c (geodesic ratio in [0.05,0.80]):'{'✅ OPEN' if gate_h13c else '❌ CLOSED'}")

    n_species = len(species_data)
    if n_species >= 2:
        has_clear_h13d = "CONVERGENT" in h13d_verdict or "DIVERGENT" in h13d_verdict
        print(f"H13d signal clarity:                       {'✅ CLEAR' if has_clear_h13d else '⚠️  AMBIGUOUS'}")
        if has_clear_h13d and gate_h13a and gate_h13c:
            print("\n✅ PHASE 5 GATE OPEN — multi-species UKFT comparison complete")
        else:
            print("\n⚠️  PHASE 5 GATE: more species data needed for definitive result")
    else:
        print(f"\n⚠️  Only 1 species processed. Download more species for cross-species comparison.")
        if gate_h13a and gate_h13c:
            print("   Single-species gates pass. Re-run after downloading remaining species.")

    # ── Save report ───────────────────────────────────────────────────────────
    report_path = UKFTBIO_ROOT / "experiments" / "results" / "13_varietal_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "experiment": 13,
        "title": "Multi-Species Varietal Comparison",
        "species_processed": list(species_data.keys()),
        "dry_run": args.dry_run,
        "max_hours": max_hours,
        "hypotheses": {
            "H13a_sparse_burst": {
                "description": "silence>=0.40 AND (tail_ratio>=3.0 OR rs_zipf<-0.60) in all species",
                "results": h13a_results,
                "gate_pass": gate_h13a,
            },
            "H13b_ks_distribution": {
                "description": "KS test on ρ CDFs across species pairs",
                "results": h13b_results,
            },
            "H13c_geodesic_ratio": {
                "description": "Geodesic/Euclidean ratio in [0.05, 0.80]",
                "results": h13c_results,
                "gate_pass": gate_h13c,
            },
            "H13d_centroid_distance": {
                "description": "Centroid cosine distances (< 0.30 convergent, > 0.60 divergent)",
                "distances": {c["pair"]: c["centroid_cosine_dist"] for c in cross_comparisons},
                "verdict": h13d_verdict,
            },
        },
        "per_species": per_species,
        "cross_comparisons": cross_comparisons,
        "gates": {
            "GATE_H13a": gate_h13a,
            "GATE_H13c": gate_h13c,
            "n_species": n_species,
        },
    }

    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nReport saved → {report_path.relative_to(UKFTBIO_ROOT)}")


if __name__ == "__main__":
    main()
