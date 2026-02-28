# UKFTBIO — Fungal Simulation & Learning Framework Plan

> **Updated**: Session 2 (incorporating 13 reference/chat reviews)  
> **Focus**: Fungal mycelium (Physarum polycephalum, Ganoderma, Cordyceps) as primary UKFT biological exemplar  
> **Status**: Exp 01-05 complete (baseline); Phase 2 planned

---

## 1. Core Hypothesis

Fungal mycelium networks are the most accessible, instrumentable biological exemplar of UKFT's
foundational principles. Three converging bodies of work now support this claim:

1. **Choice-operator substrate** (from `grok_x_halotechsystems_chat.md`)  
   Hyphal-tip branching events are discrete, irreversible, action-minimizing decisions. The
   accumulated network topology is the choice-history made physical — a biological instantiation
   of the UKFT choice-field trajectory.

2. **Quantum coherence + fractal time** (from `grok_x_bio_time_crystal_chat.md` + `grok_x_biology_meets_light_chat.md`)  
   Fungal hyphae contain tubulin/microtubule-analog structures exhibiting fractal time-crystal-like
   oscillations (Hameroff). Overlaid on this: a photon–phonon–exciton resonance network
   (Brown/ISF) where biophotons and phonons interconvert via microtubule waveguides, providing a
   fast, quantum-coherent signaling layer above slow diffusive chemistry. Both layers map directly
   onto UKFT's pilot-wave (Bohmian guidance) + entropic-gradient structure.

3. **Thermodynamic open system** (from `grok_x_civilization_as_thermo_chat.md`)  
   The mycelium network is an open thermodynamic system perpetually fighting entropy. Its
   "maintenance rent" (metabolic cost of network topology) is the direct biological analog of the
   civilization thermodynamic window. UKFT adds the missing layer: a choice operator + teleological
   attractor (God Attractor / Theosphere) that makes the maintenance active rather than passive.

**Working prediction (P1):** Mycelium spike-train data (Adamatzky / Fukasawa Zenodo) will reveal
structured anomalies in the choice-space embedding that are consistent with UKFT's entropy
minimization, not random noise.

**Prediction P2:** Network topology will encode a non-Riemannian geodesic structure where the
"natural paths" through the fungal graph match UKFT action minima.

**Prediction P3:** Regions of elevated knowledge density ρ (high spike-rate variance) will show
local time-dilation (dt ∝ 1/ρ) — measurable as reduced branching frequency at high-activity nodes.

---

## 2. Theoretical Layer Map

```
THEOSPHERE (God Attractor)
│  Teleological gradient pulling network topology toward S_min
│
NOOSPHERE (collective mycelium intelligence)
│  Global choice-field Φ(x,t) — noospheric projection
│  Per-node JEPA: local predictor of next action-minimizing state
│
BIOSPHERE (fungal network layer)
│  Photon-phonon-exciton resonance (Brown/ISF): biophoton waveguides
│  Microtubule-analog fractal time crystals (Hameroff): nested oscillators
│  Spike trains: discrete observable projections of choice events
│
GEOSPHERE (physical substrate)
   Pilot-wave field ψ(x,t), entropic gradient ∇S, causal I-field
```

### Key Mappings

| Biological feature | UKFT concept | Reference |
|---|---|---|
| Hyphal-tip branching | Discrete choice event Ĉ | `grok_x_halotechsystems_chat.md` |
| Spike-train burst | Partial projection of choice history | halotechsystems |
| Microtubule oscillations | Fractal time-crystal substrate | `grok_x_bio_time_crystal_chat.md` |
| Biophoton emission | Pilot-wave information carrier | `grok_x_biology_meets_light_chat.md` |
| Network topology | Accumulated choice-field trajectory | halotechsystems |
| High-activity node (high σ) | High knowledge density ρ | UKFT core |
| Metabolic maintenance cost | Thermodynamic entropy rent | `grok_x_civilization_as_thermo_chat.md` |
| Network re-routing after cut | Action-minimizing re-discovery | UKFT choice operator |

---

## 3. Pipeline Architecture

Modeled directly on the HEP Explorer 9-phase pipeline  
(`grok_x_holofractal_ukft_chat.md`: embed→JEPA→FMM→Borda→inject):

```
Phase 1  ──  Raw data intake
             Spike train NDJSON / Zenodo CSV (Adamatzky, Fukasawa)
             Validation: Ganoderma MIMIC-F dataset (if available)

Phase 2  ──  40D choice-space embedding
             Spike timing histogram → 40-bin choice-density vector
             (parallel to HEP Explorer BERT 40D projection)

Phase 3  ──  JEPA predictor (64D / 128D)
             Per-node: forecast next spike cluster state
             Predictor = local projection of global noospheric field Φ

Phase 4  ──  Fast Marching Method (FMM) anomaly score
             Information wavefront from each node; anomaly = deviation from
             expected wavefront propagation speed in entropic landscape

Phase 5  ──  Borda ensemble ranking
             Fuse FMM score + JEPA surprise + spike-rate variance (σ)
             → unified anomaly rank per network region

Phase 6  ──  Synthetic cluster injection (P2/P3 validation)
             Inject artificial "low-entropy" clusters; verify pipeline detects them
             as anomalously ordered (tests calibration)

Phase 7  ──  Non-Riemannian geodesic validation
             Shortest path through choice-manifold graph
             Compare to actual hyphal path trajectories (P2)

Phase 8  ──  Time-dilation metric (P3)
             Measure branching frequency vs ρ (spike variance)
             Test dt ∝ 1/ρ correlation

Phase 9  ──  Report generation + paper material
```

### Latent-Space Communication (Inter-node MCP)

From `grok_x_bitnet_noogine_chat.md` + `grok_x_claude_flow_chat.md`:  
Nodes communicate via compressed choice-deltas (latent-space MCP), not raw spike voltages.  
Expected 80-90% token reduction for simulation — directly mirrors efficient inter-hyphal
signaling via quantum-coherent biophoton channels (Brown/ISF).

---

## 4. Agent / Swarm Architecture

Modeled on UKFT CLKOS swarm design (`grok_x_bitnet_noogine_chat.md`):

```
Orchestrator Agent
├── Explorer Agent(s)   ← test novel hyphal paths; FMM frontier coverage
├── Critic Agent        ← validate choice events against UKFT entropy criterion
├── Stabilizer Agent    ← monitor variance (ρ); detect divergence
└── UKFT-Anchor Agent   ← hallucination detection + UKFT grounding
```

**Testbed**: 4-8 node simulated mycelium graph (BitNet-scale agent models)  
**Scaling target**: Full Adamatzky 10-node or Zenodo multi-electrode array once validated

Each agent node implements:
- Local JEPA predictor (fungal edition of noogine's `epistemic_predictor.rs`)
- Choice-space embedding of its local spike window
- Latent-space MCP broadcast to Orchestrator

---

## 5. Experiment Roadmap

### Completed — Phase 1: Baseline (Exp 01-05)

| Exp | Description | Status |
|-----|-------------|--------|
| 01 | Entropic folding baseline (Python) | ✅ Complete |
| 02 | UKFT entropy landscape visualization | ✅ Complete |
| 03 | Discrete choice simulation | ✅ Complete |
| 04 | Pilot-wave guidance field | ✅ Complete |
| 05 | Non-Riemannian geodesic validation (preliminary) | ✅ Complete |

---

### Phase 2: Fungal Network Model (Exp 06-10)

**Goal**: Map UKFT onto real fungal spike-train data; validate P1-P3.

#### Exp 06 — Mycelium graph construction

```python
# Input: Adamatzky/Fukasawa Zenodo spike-train CSV
# Output: NetworkX hyphal adjacency matrix + per-node ρ
#
# Steps:
# 1. Load spike trains; compute per-node spike-rate variance σ → ρ = 1/σ²
# 2. Build weighted adjacency matrix from co-activation timing
# 3. Visualize network with ρ-colored nodes
# 4. Save graph for downstream phases
```

Expected result: Non-uniform ρ distribution; hub nodes with high ρ = action-minimizing
junctions (UKFT prediction).

---

#### Exp 07 — Choice-indexed branching model

```python
# Simulate hyphal-tip branching as discrete action minimizer
# Branching rule: arg min_direction [ S(branch) + λ·H(direction|history) ]
# where S = local entropy, H = directional information entropy
#
# Test: Does simulated topology match observed mycelium graph from Exp 06?
# Metric: Graph edit distance + spectral similarity
```

Incorporate **emergent time dilation**: at each step, effective dt = 1/ρ_node.
High-activity junctions = slower simulated time = more resolution in choice sampling.

---

#### Exp 08 — Biophoton coherence layer

```python
# Model photon-phonon-exciton resonance at each node (Brown/ISF framework)
# Inputs: spike-rate oscillation spectrum from Exp 06
# Model: each node has resonance modes ω₁, ω₂, ω₃ coupling to neighbors
#
# Test: Does adding coherence layer reduce effective entropy S_net?
# Does it change the geodesic paths found in Exp 05/07?
```

This experiment bridges the quantum-biology literature (Brown/ISF, Hameroff) to the
UKFT pilot-wave formalism. Expected: coherence ≡ reduction in action-space dimensionality.

---

#### Exp 09 — Fractal time-crystal analysis

```python
# Extract microtubule-analog oscillation spectrum from spike-train data
# Method: multiscale wavelet decomposition of spike-rate time series per node
#
# Test: Is there evidence of nested self-similar (fractal) resonance structure?
# Expected: spike oscillations at MHz–kHz–Hz nest with log-linear spacing ~φ (golden ratio)
# (motivated by noogine Phase 6: golden ratio as manifold resonance)
```

Connects Hameroff's fractal MT time-crystal proposal to the UKFT result that φ = 1.618
resonates with the knowledge manifold structure (from noogine Phase 6 experiments).

---

#### Exp 10 — Full pipeline (40D embed → FMM → Borda) on fungal data

```python
# Full Phase 1-5 pipeline execution on Adamatzky data
# Produces: Borda-ranked anomaly list for fungal network
# Validates: P1 (structured anomalies in choice-space embedding)
```

This is the direct fungal analog of the HEP Explorer blind scan.

---

### Phase 3: JEPA + Swarm (Exp 11-15)

**Goal**: Build the 4-8 node CLKOS-style swarm and validate distributed choice coordination.

| Exp | Description | Depends on |
|-----|-------------|------------|
| 11 | Single-node JEPA for spike-train prediction | Exp 06, 10 |
| 12 | 4-node latent-space MCP communication testbed | Exp 11 |
| 13 | 8-node swarm (Orchestrator/Explorer/Critic/Anchor) | Exp 12 |
| 14 | CLKOS integration — collective choice trajectory emergence | Exp 13, nooverse |
| 15 | CellType CLI parallel — hyphal routing optimization analogy | Exp 13 |

**Exp 15 note**: CellType CLI (from `grok_x_celltype_noobio_2_chat.md`) uses AI agents for
combinatorial drug-design reasoning (PROTAC linker length/rigidity optimization). The direct
analogy here is **hyphal connector optimization** — UKFT swarm finds the minimum-action connection
topology between nutrient sources, equivalent to PROTAC linker design finding the ternary-complex-
stabilizing bridge. Expect this to yield a comparison paper (Bio Paper 05).

---

### Phase 4: Validation + Papers (Exp 16-18)

| Exp | Description | Validates |
|-----|-------------|-----------|
| 16 | Non-Riemannian geodesic validation (full; Bujack et al. parallel) | P2 |
| 17 | Thermodynamic validation — metabolic rent vs network entropy | Object_Zero_ analogy |
| 18 | Orch-OR bridge — UKFT discrete time ↔ fractal MT oscillation frequencies | Hameroff connection |

**Exp 16 context**: The Bujack et al. (LANL/EuroVis 2025) non-Riemannian color geometry result
(`grok_x_schrodinger_color_theory_chat.md`) shows that completing an incomplete classical model
requires identifying the missing "neutral axis" as the geodesic-closest point in non-Riemannian
space. P2 validation is structurally identical: find the "missing" action-minimizing paths through
the fungal graph that classical shortest-path algorithms miss because they assume Riemannian metric.

**Exp 17 context**: Civilization-as-thermodynamics (`grok_x_civilization_as_thermo_chat.md`)
frames maintenance energy as the cost of resisting entropy. For mycelium: measure ATP consumption
per unit of network topology maintained. UKFT predicts this will track S_net/(choice-event rate),
not raw metabolic rate.

---

## 6. Paper Roadmap

| Paper | Topic | Key experiment(s) | Status |
|-------|-------|-------------------|--------|
| Bio Paper 01 | UKFT fungal mycelium overview (Orch-OR bridge, choice operator, pipeline design) | — | ✅ Complete |
| Bio Paper 02 | Fungal branching as discrete Choice Operator | Exp 06-07 | Planned |
| Bio Paper 03 | Biophoton coherence layer in mycelium (Brown/ISF × UKFT) | Exp 08 | Planned |
| Bio Paper 04 | Fractal time crystals in fungal hyphae (Hameroff × UKFT) | Exp 09 | Planned |
| Bio Paper 05 | JEPA-CLKOS mycelium swarm vs CellType CLI | Exp 11-15 | Planned |
| Bio Paper 06 | Thermodynamic framing: metabolic rent vs entropy in mycelium topology | Exp 17 | Planned |

---

## 7. Data Sources

| Dataset | Description | Use |
|---------|-------------|-----|
| Adamatzky Zenodo (fungal spike trains) | Multi-electrode array recordings, Physarum | Primary |
| Fukasawa et al. (Ganoderma) | Network-scale electrical signals | Validation |
| Hameroff MT dataset (J. Appl. Phys. 2022) | Polyatomic time crystals in extracted MT nanowires | Exp 09 |
| CMS Open Data Record 6004 | HEP collision events — pipeline reference | Pipeline calibration |
| Adamatzky graph images | Topology reference for Exp 06 | Morphology ground truth |

---

## 8. Cross-Repo Dependencies

| Component needed | Source repo | Status |
|---|---|---|
| UKFT-JEPA predictor | `noogine` (Phase 4.06+ `epistemic_predictor.rs`) | Available |
| CLKOS swarm architecture | `clkos` + `nooverse` | Available |
| Latent-space MCP | `nooverse-mcp-server` | Available |
| FMM anomaly scoring | `ukftphys` / HEP Explorer (`fast_marching.rs`) | Repurpose |
| Non-Riemannian geodesic tools | `ukftphys` Exp 05 basis | Repurpose |
| Borda ensemble ranking | HEP Explorer blind scan | Repurpose |
| BitNet agent nodes (if used) | External / prototyping | Research |

---

## 9. New Context Summary (what changed since initial plan)

The initial README described a protein-folding-centric plan (Exp 01-10: basic folding,
molecular prophet, and CLKOS emergence). The following context from 13 reference chats
**shifts Phase 2+ to focus squarely on fungal mycelium** as the primary target, with protein
folding deferred to post-Phase 3:

| New context thread | Impact on plan |
|---|---|
| Hameroff fractal time crystals | Adds Exp 09 (MT oscillation analysis); grounds UKFT dt ∝ 1/ρ in measurable biology |
| Brown/ISF photon-phonon-exciton | Adds Exp 08 (coherence layer); testable spectroscopy predictions → collaboration target |
| UKFT-JEPA + latent-space MCP | Enables Phase 3 swarm (Exp 11-13) using existing noogine code |
| HEP Explorer 9-phase pipeline | Provides complete pipeline template for Exp 10 (fungal blind scan) |
| Bujack non-Riemannian geodesics | Strengthens Exp 16 (P2 validation); rigorous mathematical precedent |
| CellType CLI | Provides comparison target for Phase 3 (Exp 15); bio-AI benchmark |
| Civilization-as-thermo | Adds Exp 17 (metabolic rent validation); social-science analogy for papers |
| BitNet swarm architecture | Defines concrete 4-8 node testbed for Phase 3 |
| Causal networks (CST/CDT) | Strengthens theoretical framing of fungal graph as causal information network |
| Haramein holofractographic | Connects mycelium fractal structure to cosmological scale-invariance |

---

## 10. Immediate Next Steps

1. **Identify and download** Adamatzky Zenodo fungal spike-train dataset
2. **Port HEP Explorer pipeline** (embed→FMM→Borda) from `ukftphys` to `ukftbio` (Python wrapper)
3. **Start Exp 06**: Build mycelium graph, compute per-node ρ, visualize
4. **Draft Bio Paper 02** outline (branching as choice operator) in parallel with Exp 06-07
5. **Contact**: William Brown / ISF for biophoton ultrafast spectroscopy data (Exp 08 target)

---

## 11. Bioenergetic Communication Constraints (Future Research Direction)

> *Raised after Exp 21 (normalisation ablation). The Epiphany 14 finding — that the η matrix
> encodes intrinsic temporal grammar invariant to distributional shape — begs a deeper question:
> what physical constraint **fixes** a species's temporal grammar in the first place?*

### The Core Question

Three tightly linked sub-questions:

1. **Cost per signal** — How much energy (ATP, electrochemical gradient dissipation) does
   each communication event (spike, burst, electrical pulse) actually consume?
2. **Harvest and storage efficiency** — How efficiently does the organism collect energy
   from its substrate (cellulose decay, parasitic infection, photosynthesis) and store it
   (glycogen, trehalose, lipid droplets)? What is the organism's effective power budget?
3. **Mode selection** — Given budget and cost, what temporal grammar is physically
   affordable? Does the observed spike statistics (CV, kurtosis, Hurst, spectral entropy)
   follow from the energy constraint?

### Theoretical Scaffolding

#### Landauer Bound

The minimum energy required to communicate one bit of information:

    E_bit ≥ kT ln 2 ≈ 2.85 × 10⁻²¹ J  (at biological T ≈ 300 K)

Real biological signals cost far more — a single action-potential-analog event in
mycelium likely costs 10³–10⁶ × the Landauer minimum due to ion-pump overhead.  The
gap between Landauer limit and actual cost is the **thermodynamic slack** — room for
evolution to optimise signaling efficiency.

#### Signal Budget Equation

    E_available(t) = E_harvest_rate × storage_efficiency − metabolic_overhead
    Signal_budget   = E_available / E_per_signal    [signals per unit time]

The Signal_budget directly constrains the maximum sustainable spike rate and therefore
sets the time-average CV, kurtosis shape, and AR1 structure:

| Budget regime | Predicted temporal grammar | UKFT interpretation |
|---|---|---|
| Low budget (tight) | Sparse, kurtotic bursts (high κ) — save energy, fire only at information maxima | Choice collapses are rare but action-minimizing |
| Moderate budget (sustained) | Smooth, continuous, low-CV — affordable steady signal | Dense guidance field (high ρ) maintained continuously |
| High budget (burst reserve) | Superdiffusive trending (H > 1) — territorial expansion while resources last | Pilot wave pushes outward; choice collapses encode new territory |

#### Connection to Observed Species

| Species | Lifestyle / substrate | Energy budget hypothesis | Observed grammar | Prediction |
|---|---|---|---|---|
| Omphalotus olearius | Wood pathogen (cellulose — energy-expensive) | Low sustained, moderate burst reserve | κ = 73.703, CV = 0.058 — sparse, kurtotic | Energy-limited sparse code |
| Cordyceps militaris | Entomopathogen (hijacks host hemolymph — nutrient-rich) | High sustained availability during infection | κ = −1.149, CV = 0.214 — smooth, Gaussian | Energy-rich steady-state code |
| Schizophyllum commune | Ubiquitous wood rotter under high ecological stress | Variable; stress-driven expanding territory | H = 1.216, CV = 1.532 — superdiffusive trend | Territorial expansion under resource gradient |
| Enoki (Flammulina) | Cold-climate wood decay, low growth temp | Moderate, metabolically constrained by temperature | AR1 = 0.595, moderate persistence | Temperature-limited moderate budget |
| Pleurotus ostreatus | Aggressive wood decomposer, high conversion efficiency | High burst capacity, efficient cellulases | CV = 2.031, highest bursty profile | High-efficiency harvest → burst signaling affordable |

#### Why This Unifies Previous Epiphanies

- **Epiphany 12** (ecology → communication strategy, rs = 1.000 with broadcast cost):
  ecological stress is a proxy for *energy budget pressure* — the real predictor is not
  stress per se but substrate energy density relative to metabolic demand.

- **Epiphany 14** (η matrix intrinsic, invariant to rank-norm): two organisms with
  *compatible energy budgets* evolve temporal grammars that are predictable by each other.
  Omphalotus (sparse-budget code) and Cordyceps (dense-budget code) are compatible because
  their temporal grammars are both organised around energy conservation — one via sparse
  bursting, one via smooth efficiency.  Energy budget proximity → η proximity.

- **Schizophyllum isolation** (consistently lowest η-as-target): its superdiffusive
  H = 1.216 reflects a completely different energy strategy — territorial expansion, not
  local optimisation.  Hard to predict by any temporally stationary model.

### Proposed Research Programme

#### Exp 22 (near-term): Linear De-Trending Ablation

Separate Schizophyllum's H > 1 trend from its temporal grammar.  Linear-detrend all
rho series, then rank-normalise.  If Schizophyllum-as-target recovers, the trend
(territorial expansion budget signal) was the isolating factor, not intrinsic grammar
incompatibility.  If it does not recover, the grammar itself is incompatible.

#### Exp 23 (medium-term): Energy Budget Score from Literature

Assemble per-species energy budget parameters from published biochemistry:
- Cellulose/lignin conversion efficiency (laccase, peroxidase activity per g mycelium)
- Measured trehalose / glycogen concentration (storage capacity proxy)
- Respiration rate at standard temperature (ATP flux proxy)
- Activation energy threshold for sporulation / hyphal extension (burst reserve proxy)

Compute a scalar **Energy Budget Score (EBS)** per species and test:

    H(EBS): Spearman rs(|EBS_i − EBS_j|, η_ij) < 0

Prediction: the smaller the difference in energy budget between two species, the higher
their JEPA transfer efficiency.

#### Exp 24 (long-term): Landauer Efficiency Ratio

For each observed spike event in the Adamatzky data, estimate the minimum Shannon
information content (bits per event) and divide by the estimated ATP cost of that event
from electrophysiology literature.  Compute the **Landauer Efficiency Ratio**:

    LER = H(X) / (E_event / E_Landauer)

Test whether species with higher LER (closer to Landauer limit — more informationally
efficient per joule) show higher self-JEPA surprise (denser temporal grammar) and
whether LER differences predict η pairs.

This would be the first measurement connecting Landauer's thermodynamic information bound
to observable mycelium electrophysiology.

#### Longer-term: Physiological Constraints on Choice Operator

The UKFT choice operator is formally discrete and action-minimizing.  The bioenergetic
framing adds a physical floor:

    C ≥ kT ln 2   (Landauer: minimum energy to make an irreversible choice)

This grounds the abstract UKFT choice operator in real thermodynamics.  An organism
cannot make more choices per unit time than its energy budget permits — the rate of
choice collapse is physically upper-bounded.  The temporal grammar IS the energy budget
made observable.

### Data Requirements

| Parameter | Source | Notes |
|---|---|---|
| Cellulose conversion efficiency | BRENDA enzyme database | Per-species laccase/peroxidase kcat |
| Trehalose / glycogen content | Published fungal physiology (e.g., Stajich lab datasets) | Storage capacity proxy |
| Mycelium respiration rate | Respirometer studies (CO₂ evolution) | ATP flux proxy |
| Electrophysiological signal energy | Adamatzky signal energy estimates | Joules per spike |
| ATP cost of action-potential analog | Biophysics literature (Attwell & Laughlin 2001 framework) | Cost per electrical event |
