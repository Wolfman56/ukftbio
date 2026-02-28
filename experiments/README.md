# Experiments — UKFT Biology

Numbered experiments following the ukftphys convention.
Each experiment has a `.py` script and a `.md` companion explainer.

> **Design philosophy (mirrors ukftphys):**
> Establish the canonical **baseline set (Exp 01–05)** first — confirming the
> framework works at the most fundamental level before entering the
> paper / hypothesis / proof cycle.
> Exp 02 (Levinthal) is the bio double-slit.

---

## Phase 1 — Baseline Set (Exp 01–05)

Confirm the UKFT biological framework before advancing to the fungal pipeline.

| Exp | Title | Physics Analogy | Key UKFT Concept | Status |
|-----|-------|----------------|-----------------|--------|
| 01 | Minimal Entropic Replicator | Free particle | Error threshold as S_bio floor | ✅ Ready |
| 02 | Levinthal Paradox | Double slit | Folding = choice collapse, not random search | ✅ Ready |
| 03 | Error Threshold Sweep | Entropic alpha sweep | fidelity / ρ phase boundary | ✅ Ready |
| 04 | Autocatalytic Emergence | Swarm / God attractor | Life as ρ phase transition | ✅ Ready |
| 05 | Fungal Choice Spike | Entanglement propagation | Entangled cluster = EPR-bounded collapse | ✅ Ready |

### Running Phase 1

    conda activate ukftbio
    python experiments/01_minimal_entropic_replicator.py
    python experiments/02_levinthal_folding_paradox.py
    python experiments/03_error_threshold_sweep.py
    python experiments/04_autocatalytic_emergence.py
    python experiments/05_fungal_choice_spike.py

### Phase 1 Canonical Figures

| Exp | Primary Figure | What it proves |
|-----|---------------|----------------|
| 01 | `exp01_rho_trajectory.png` | ρ survives (f=0.99) vs collapses (f=0.80) — error threshold |
| 02 | `exp02_folding_landscape.png` | UKFT collapses to native basin; random walk lost in Levinthal space |
| 03 | `exp03_phase_diagram.png` | Sharp f_critical boundary in (fidelity, ρ_final) space |
| 04 | `exp04_rho_emergence.png` | ρ jumps discontinuously at autocatalytic threshold — God basin |
| 05 | `exp05_cluster_resolution.png` | Real clusters resolve at action minimum; shuffled controls do not |

---

## Phase 2 — Fungal Mycelium Pipeline (Exp 06–10) ✅ Complete

Fungal mycelium as base-level collective consciousness: the same UKFT equation
`S = Tr(log G_truth − log G_post)` operating at molecular gradient resolution
(τ* ≈ 68 s), separated from the human noosphere by exactly φ²⁵ in temporal depth.

All experiments use the same 8-channel synthetic spike-train substrate (seeded from
Exp 05 parameters) and build incrementally toward the full HEP-Explorer-style blind-scan
pipeline in Exp 10.

| Exp | Title | Key Result | Commit |
|-----|-------|-----------|--------|
| 06 | Mycelium Graph Construction | Hub = ch3 (ρ=46 770); ρ CoV=0.326 | `9ccc471` |
| 07 | Choice-Indexed Branching | UKFT λ=0.25; time-dilation r=**+0.976** (p≈0) | `1a0daa7` |
| 08 | Biophoton Coherence Layer | PLV baseline=0.147; 39% φ-ratio frequencies | `edde7b1` |
| 09 | Fractal Time-Crystal Analysis | H=0.593; α=1.19 (R²=0.888); **τ* ratio ≈ φ²⁵** | `78dac57` |
| 10 | Full 40D → FMM → Borda Pipeline | P2: 14% geodesic divergence; T3: 40% φ-ratios | ✅ |

### Running Phase 2

    conda activate ukftbio
    python experiments/06_mycelium_graph.py
    python experiments/07_choice_branching.py
    python experiments/08_biophoton_coherence.py
    python experiments/09_fractal_timecrystal.py
    python experiments/10_full_pipeline.py

### Phase 2 Key Numbers

| Metric | Value | Experiment |
|--------|-------|-----------|
| Hub channel ρ | 46 770 (ch3) | Exp 06 |
| ρ CoV across network | 0.326 | Exp 06 |
| Optimal UKFT branching λ | 0.25 | Exp 07 |
| Time-dilation Spearman r | **+0.976** | Exp 07 |
| PLV decoherent baseline | 0.147 | Exp 08 |
| φ-ratio frequency fraction | 39% | Exp 08 |
| Mean Hurst H | **0.593** | Exp 09 |
| Mean wavelet α | **1.19**, R²=0.888 | Exp 09 |
| Optimal temporal horizon τ* | 68 s (1.1 min) | Exp 09 |
| Noosphere / fungus τ* ratio | ≈ φ²⁵ | Exp 09 |
| UKFT vs Dijkstra path divergence | **14%** | Exp 10 |
| Welch peak φ-ratio fraction | **40%** | Exp 10 |
| Borda top channel (baseline) | ch6 rank #1 | Exp 10 |

### Phase 2 Canonical Figures

| Exp | Figure | What it shows |
|-----|--------|--------------|
| 06 | `exp06_rho_topology.png` | ρ-weighted co-activation graph; hub ch3 |
| 06 | `exp06_rho_distribution.png` | ρ histogram; CoV=0.326 |
| 07 | `exp07_branching_heatmap.png` | Choice-indexed branching vs UKFT λ |
| 07 | `exp07_time_dilation.png` | Spike-rate × ρ — time-dilation r=+0.976 |
| 08 | `exp08_plv_matrix.png` | Phase-locking value across all channel pairs |
| 08 | `exp08_coherence_phi_ratios.png` | 39% of peak freq ratios near φ or φ² |
| 09 | `exp09_hurst_dfa.png` | DFA Hurst exponents (H̅=0.593) |
| 09 | `exp09_wavelet_scaling.png` | Wavelet scaling α̅=1.19 |
| 09 | `exp09_phi25_staircase.png` | Consciousness scale staircase — φ²⁵ gap |
| 10 | `exp10_embedding_pca.png` | 40D choice-space embeddings (PCA projection) |
| 10 | `exp10_fmm_scores.png` | FMM wavefront deviation per channel |
| 10 | `exp10_borda_ranking.png` | 5-metric Borda fusion ranking |
| 10 | `exp10_anomaly_detection.png` | P1: injected vs baseline Borda ranks |
| 10 | `exp10_geodesic_validation.png` | P2: UKFT biased vs Dijkstra paths |
| 10 | `exp10_cwt_phi_ratio.png` | T3: Welch peak freq ratios vs φ |
| 10 | `exp10_consciousness_staircase.png` | Exp 06–10 staircase — full Phase 2 integration |

### P1 Partial Result — Notes

Exp 10 P1 (synthetic anomaly ch8 ranking in top-3) was **PARTIAL**: ch8 ranked #9 because
near-zero-entropy channels have *low* σ-variance, which dominates the Borda metric.
This is expected and informative — the HEP-style Borda ensemble captures **maximal-activity**
anomalies by default, not **minimal-entropy** anomalies. A re-weighted Borda inverting
the σ metric (rewarding low variance) will be addressed in Phase 3.

---

## Phase 3 — Real Data + JEPA + Multi-Species (Exp 11–13) ✅ Complete

Moved from synthetic to real electrophysiology (Adamatzky Zenodo), trained a
JEPA temporal predictor, and confirmed universal manifold geometry across 5 fungal
species. Full writeup: `docs/PHASES_10_13_REPORT.md`.

| Exp | Title | Key Result | Status |
|-----|-------|-----------|--------|
| 11 | Real Fungal Data Ingest | time-dilation r=**+0.392** (p=1.3×10⁻⁵⁰); geodesic ratio=**0.183**; recall@100=1.00 | ✅ |
| 12 | JEPA Temporal Predictor | rs(S,ρ)=**+0.8785**; cos_acc=**0.9519** — *Epiphany 9: ρ–S duality* | ✅ |
| 13 | Multi-Species Varietal Comparison | geodesic ratio **0.319–0.345** across 5 species; silence varies **3,300×** — *Epiphany 10: universal geometry* | ✅ |

### Running Phase 3

    conda activate ukftbio
    python experiments/11_real_fungal_ingest.py
    python experiments/12_jepa_real_data.py
    python experiments/13_varietal_comparison.py               # all species
    python experiments/13_varietal_comparison.py --dry-run     # smoke test

### Phase 3 Key Numbers

| Metric | Value | Experiment |
|--------|-------|-----------|
| Spike count (Pleurotus, 217 h) | 2,534 | Exp 11 |
| Burst fraction | 85.7% | Exp 11 |
| Active window fraction | 29.0% | Exp 11 |
| Time-dilation Spearman r | **+0.392**, p=1.3×10⁻⁵⁰ | Exp 11 |
| Geodesic/Euclidean ratio median | **0.183** (5.5× compression) | Exp 11 |
| Recall@50 / Recall@100 | 0.80 / **1.00** | Exp 11 |
| JEPA rs(surprise S, ρ) | **+0.8785** | Exp 12 |
| JEPA cosine accuracy | **0.9519** | Exp 12 |
| JEPA rs(S, Borda) | +0.7880 | Exp 12 |
| 5-species geodesic ratio range | **0.319–0.345** (band=0.026) | Exp 13 |
| Silence range across species | **0.03%–98.5%** (3,300× variation) | Exp 13 |
| Zipf rs (all species) | < −0.982 | Exp 13 |
| Centroid cluster structure | TWO_CLUSTER (sparse-episodic vs continuous) | Exp 13 |

### Phase 3 Canonical Figures

| Exp | Figure | What it shows |
|-----|--------|--------------|
| 11 | `exp11_pipeline_comparison.png` | Exp 10 synthetic vs Exp 11 real side-by-side |
| 11 | `exp11_rho_distribution.png` | ρ histogram + Zipf rank-log (Pleurotus) |
| 11 | `exp11_time_dilation.png` | log-log scatter ρ vs ISI — UKFT dt=C/ρ fit |
| 11 | `exp11_geodesic_ratio.png` | Geodesic/Euclidean distribution (median=0.183) |
| 11 | `exp11_borda_topk.png` | Top-10 Borda candidates on recording timeline |
| 12 | `fig_12a_jepa_surprise_rho.png` | JEPA S vs ρ scatter + all rs bar chart |
| 13 | `fig_13a_geodesic_universal.png` | 5-species geodesic ratio — universality band |
| 13 | `fig_13b_silence_spectrum.png` | Silence fraction + heavy-tail ratio per species |
| 13 | `fig_13c_centroid_heatmap.png` | 5×5 centroid cosine distance — TWO_CLUSTER |

### Epiphanies

**Epiphany 9** (Exp 12): ρ (retrospective entropy production) and JEPA surprise S
(prospective entropy production) are correlated at rs=+0.88. This is the signature
of irreversible choice collapse: the same trajectory that generated high past density
is unpredictable next. Silent attractor: S≈0.0004; burst events: S≈1.33.

**Epiphany 10** (Exp 13): UKFT manifold geometry (geodesic ratio ≈ 0.32–0.35) is
conserved across 5 evolutionary lineages spanning 0.03%–98.5% silence. A single
Pleurotus-trained projection head transfers without re-training. Architecture =
universal; activation density = biological.

---

## Phase 4 — Next Steps (Exp 14+)  📋 Planned

| Exp | Title | Description |
|-----|-------|-------------|
| 14 | Cross-Species JEPA Transfer | Train JEPA on Pleurotus, score on all 5 species — does temporal structure transfer? |
| 15 | Cross-Species Anomaly Scan | HEP-explorer blind scan across all 5 species jointly — find inter-species geometric analogues |
| 16 | Duroxide Orchestration | Durable-execution replay for reproducible fungal pipelines |

---

## Consciousness Staircase (Theoretical Context)

All experiments are grounded in the same `S = Tr(log G_truth − log G_post)` equation
operating at different temporal scales:

| Level | Agent | τ* (optimal horizon) | Mechanism |
|-------|-------|----------------------|-----------|
| 1 | Fungal mycelium | **68 s** | Molecular gradient routing (Exp 06–10) |
| 2 | Individual neuron | ~1 s | Spike-timing STDP |
| 3 | Human individual | ~10³–10⁵ s | Tool use, episodic memory |
| 4 | Human noosphere | ~10⁶–10⁹ s | Counterfactual generation, new environments |

Scale ratio fungus → noosphere: τ*_noosphere / τ*_fungus ≈ **φ²⁵** (Exp 09).
Same UKFT law, 25 golden-ratio steps apart in temporal depth.
