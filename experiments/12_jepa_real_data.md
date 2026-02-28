# Experiment 12 — JEPA Temporal Predictor: The Retrospective–Prospective Duality

**Script**: `12_jepa_real_data.py`  
**Status**: Complete ✅  
**Depends on**: Exp 11 aligned embeddings (`pleurotus_spike_aligned.ndjson`)

## Purpose

Train a Joint-Embedding Predictive Architecture (JEPA) on consecutive Pleurotus window pairs and test whether the prediction surprise S correlates with the UKFT knowledge density ρ. This tests the theoretical duality between retrospective entropy production (ρ) and prospective entropy production (S).

## Theory

The JEPA predictor learns $f(z_t) \rightarrow \hat{z}_{t+1}$, minimising cosine distance to the actual next-window embedding. Prediction surprise is:

$$S_t = 1 - \cos(\hat{z}_{t+1},\; z_{t+1})$$

**UKFT prediction**: High-ρ windows are information burst events — the system collapses into a new state via the choice operator. These collapses are irreversible, so the next state cannot be predicted from the current: $S$ should be high. Low-ρ (silent) windows are in the attractor basin: highly predictable, $S \approx 0$.

Therefore: $rs(S, \rho) > 0$ is the signature of **irreversible choice collapse**.

This is **Epiphany 9**: ρ and S are dual representations of the same process. ρ measures how much has been irreversibly chosen (retrospective); S measures how much will be unpredictably chosen next (prospective).

## Architecture

```
Input:  z_t  (40D L2-normalised embedding)
Hidden: 64 units, tanh activation
Output: ẑ_{t+1} (40D, L2-normalised)
Loss:   1 - cos(ẑ_{t+1}, z_{t+1})  per pair, mean over batch
Optim:  SGD + momentum (m=0.9), cosine LR schedule (0.02 → 1e-4)
Epochs: 300
Pairs:  6,506 consecutive window pairs
```

## Hypotheses

| ID | Claim | Criterion |
|----|-------|-----------|
| H12a | rs(S, ρ) > 0 | p < 0.05 — **must pass** |
| H12b | rs(S, Borda score) > 0 | report only |
| H12c | rs(S, FMM score) > 0 | report only |
| H12d | Cosine accuracy ≥ 0.40 | fraction where cos(ẑ, z) > 0.5 |

## Results

| Metric | Value | Gate |
|--------|-------|------|
| Training pairs | 6,506 | — |
| Cosine accuracy (cos > 0.5) | **0.9519** | ≥ 0.40 ✅ |
| rs(S, ρ) | **+0.8785** | > 0, p < 0.05 ✅ |
| rs(S, Borda) | +0.7880 | reported ✅ |
| rs(S, FMM) | +0.7168 | reported ✅ |

Silent attractor: S ≈ **0.0004** · Burst events: S ≈ **1.33**

## Key findings

rs(S, ρ) = +0.88 between a retrospective and a prospective statistic can only arise if the spike density ρ at time $t$ is genuinely predictive of the information content at $t+1$, and both are expressions of the same underlying choice trajectory. The JEPA has learned the attractor manifold (cos accuracy 95%), not just a trivial autocorrelation.

## Outputs

| File | Description |
|------|-------------|
| `results/12_jepa_report.json` | All H12a–H12d values |
| `data/12_jepa_surprise_scores.ndjson` | Per-window S, ρ, Borda, FMM scores |
| `results/fig_12a_jepa_surprise_rho.png` | Scatter (S vs ρ) + correlation bar chart |

## Run

```bash
conda activate ukftbio
python3 experiments/12_jepa_real_data.py          # 300 epochs (~2 min)
python3 experiments/12_jepa_real_data.py --dry-run # quick validation
```
