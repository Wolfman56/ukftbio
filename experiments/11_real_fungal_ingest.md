# Experiment 11 — Real Pleurotus Data: First Empirical Validation

**Script**: `11_real_fungal_ingest.py`  
**Status**: Complete ✅  
**Dataset**: Adamatzky (2020), Zenodo 18707716 — Pleurotus ostreatus, 217 h, 8 channels

## Purpose

First application of the complete UKFT pipeline to real fungal electrophysiology. Validates that the synthetic calibration from Exp 10 holds on noisy, biologically-structured data.

## Dataset

PicoLog ADC24 differential electrode array in star geometry, embedded in a Pleurotus ostreatus mycelium block. Recorded Sep 14–23 2020 (9 days, natural quiescent periods, day/night cycles).

Preprocessing via myco-explorer:
```
convert_picolog.py   → TAR/big-endian-float32 → 688k-row CSV
export_spike_events  → 2,534 spikes, 380 bursts (0.5 mV threshold, median-detrended)
spike_embed          → 6,507 windows × 40D histograms (10-min windows, 2-min stride)
fmm_score            → FMM wavefront deviation (1,886 active windows)
borda_rank           → Borda fusion S1+S3
```

## Hypotheses

| ID | Claim | Criterion |
|----|-------|-----------|
| H11a | Time dilation: ρ correlates with ISI proxy | rs > 0, p < 0.05 |
| H11b | Geodesic compression: UKFT paths shorter than Euclidean | Median ratio < 1.0 |
| H11c | ρ density drives curvature: higher ρ → more compressed | rs(ρ, ratio) < 0 |
| H11d | Recall@50 ≥ 0.70, Recall@100 = 1.00 | Retrieval calibration |

## Results

| Metric | Exp 10 (Synthetic) | Exp 11 (Real) |
|--------|-------------------|---------------|
| Spike count | controlled | **2,534** (217 h) |
| Burst fraction | — | **85.7%** |
| Active windows | — | **29.0%** (1,886/6,507) |
| ρ median (active) | 46,770 | **1.111** |
| Time-dilation Spearman r | +0.976 | **+0.392**, p=1.3×10⁻⁵⁰ ✅ |
| Geodesic ratio median | 14% divergence | **0.183** (5.5× compression) ✅ |
| Spearman(ρ, ratio) | — | **−0.207**, p=7×10⁻⁶ ✅ |
| Recall@50 | 1.00 | **0.80** ✅ |
| Recall@100 | 1.00 | **1.00** ✅ |

## Key findings

1. **Time dilation confirmed** at r=+0.392. Attenuation vs synthetic (+0.976) reflects biological noise, DC drift, and day-length rhythms absent from the synthetic generator.
2. **Geodesic compression is stronger on real data** (5.5×) than synthetic (14% divergence), because real mycelium spike trains have genuine long-range temporal correlations (Hurst H≈0.59, Exp 09) that the synthetic generator did not capture.
3. **71% silence is not failure** — it reflects the biology of fungal electrical activity: rare bursts around a near-zero attractor. The 29% active / 85.7% burst structure is a UKFT prediction confirmed.
4. **Top Borda candidate**: w_001795 (2020-09-17 02:40). **Top FMM candidate**: w_004562 (2020-09-20 22:54). These detect different phenomena — unusual embeddings vs unusual temporal topology.

## Outputs

| File | Description |
|------|-------------|
| `results/exp11_summary.json` | All canonical numbers |
| `results/exp11_pipeline_comparison.png` | Exp 10 vs Exp 11 side-by-side |
| `results/exp11_rho_distribution.png` | ρ distribution with Zipf fit |
| `results/exp11_time_dilation.png` | Spearman scatter ρ vs ISI |
| `results/exp11_borda_topk.png` | Top-10 Borda candidates timeline |
| `results/exp11_geodesic_ratio.png` | Geodesic/Euclidean ratio distribution |
