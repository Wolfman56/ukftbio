# Experiment 13 — Multi-Species Varietal Comparison: Universal Manifold Geometry

**Script**: `13_varietal_comparison.py`  
**Status**: Complete ✅  
**Species**: 5 × Adamatzky Zenodo fungal electrophysiology datasets

## Purpose

Test whether the UKFT minimum-action manifold geometry is substrate-independent. Five fungal species with radically different electrical activity profiles (0.03%–98.5% silence) are projected through a **single Pleurotus-trained projection head** into the same 40D embedding space. The prediction: geodesic compression ratio ≈ 0.32–0.35 should be conserved; activation density should vary freely.

## Datasets

| Species | Duration | Silence |
|---------|----------|---------|
| *Pleurotus ostreatus* (oyster) | 217 h | 71.0% |
| *Schizophyllum commune* (split-gill) | 73.3 h | 98.5% |
| *Cordyceps militaris* (caterpillar) | 527.8 h | 16.9% |
| *Flammulina velutipes* (enokitake) | 336.4 h | 11.3% |
| *Omphalotus nidiformis* (ghost fungus) | 911.0 h | 0.03% |

All datasets: PicoLog ADC24, extracellular field potential.  
All processed through identical UKFT pipeline with **shared projection head** (`tools/data/model/projection_head.safetensors`, Pleurotus-trained, W: 40×40).

## Hypotheses

| ID | Claim | Criterion |
|----|-------|-----------|
| H13a | Zipf power law holds in all species | rs(log-rank, log-ρ) < −0.60 |
| H13b | Silence distribution spans wide range | range ≥ 0.001 to 0.999 |
| H13c | ρ distributions are statistically distinct | KS p < 0.05 across species pairs |
| H13d | Centroid distances show clear bimodal structure | TWO_CLUSTER (n_close ≥ 1 AND n_distant ≥ 1) |

## Results

| Species | Windows | Geodesic ratio | Silence | Heavy-tail ratio | rs_zipf |
|---------|---------|---------------|---------|-----------------|---------|
| *Pleurotus* | 6,507 | **0.344** | 71.0% | 3.0× | −0.995 |
| *Schizophyllum* | 2,196 | **0.345** | 98.5% | 174.3× | −0.982 |
| *Cordyceps* | 15,831 | **0.324** | 16.9% | 10.8× | −1.000 |
| *Flammulina* | 10,088 | **0.319** | 11.3% | 5.5× | −1.000 |
| *Omphalotus* | 27,326 | **0.319** | 0.03% | 3.6× | −1.000 |

**Geodesic ratio range**: 0.319–0.345 (band width = **0.026**)  
**Silence range**: 0.03%–98.5% (**3,300× variation**)

### Centroid cluster structure (TWO_CLUSTER)

- **Sparse-episodic cluster**: Pleurotus + Schizophyllum, centroid distance **0.037**
- **Continuously-active cluster**: Cordyceps + Flammulina + Omphalotus, internal distances 0.111–0.296
- Cross-cluster maximum: Schizophyllum ↔ Omphalotus = **0.827**

All hypothesis gates: ✅ PASS

## Epiphany 10 — Universal Manifold Geometry

> **The UKFT minimum-action manifold geometry is substrate-independent: geodesic compression ≈ 0.32–0.35 is conserved across 5 evolutionary lineages spanning 0.03%–98.5% silence. Architecture = universal; activation density = biological.**

A single Pleurotus-trained projection head transfers without re-training to four additional species. The manifold shape is a consequence of the minimum-action principle, not the organism's ecology.

## Outputs

| File | Description |
|------|-------------|
| `results/13_varietal_report.json` | Full per-species metrics + hypothesis gates |
| `results/fig_13a_geodesic_universal.png` | Geodesic ratio bar chart with universality band |
| `results/fig_13b_silence_spectrum.png` | Silence fraction + heavy-tail ratio per species |
| `results/fig_13c_centroid_heatmap.png` | 5×5 cosine distance matrix |

Intermediate aligned embeddings cached in myco-explorer:
- `tools/data/multispecies/*_aligned.ndjson` (per-species, auto-generated)

## Run

```bash
conda activate ukftbio
python3 experiments/13_varietal_comparison.py               # all available species
python3 experiments/13_varietal_comparison.py --dry-run     # 2-species smoke test
python3 experiments/13_varietal_comparison.py --species pleurotus,cordyceps
python3 experiments/13_varietal_comparison.py --list-species
```
