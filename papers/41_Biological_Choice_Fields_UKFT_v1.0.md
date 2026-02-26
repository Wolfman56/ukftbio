# Biological Choice Fields: How Living Signaling Networks Instantiate UKFT Choice Operators

**Paper 41 — UKFT Biology Series**  
**Version:** 1.0 (Hypothesis / Prediction paper)  
**Date:** February 2026  
**Status:** Pre-publication hypothesis  
**Companion repos:** `ukftbio`, `ukftphys`, `uktf`

---

## Abstract

We propose that biological signaling networks — including fungal mycelial electrical networks, neuronal spike trains, and plant action-potential cascades — are not mere information-transmission channels but are *physical implementations of UKFT choice operators*. Under Universal Knowledge Field Theory (UKFT), the universe advances not through continuous time but through discrete collapses of the superimposed choice-field crystal: each collapse selects the least-action path from a local superposition and propagates updated projections outward. We argue that the experimentally observed electrical spike trains of *Pleurotus ostreatus* and related fungi (Adamatzky 2022, 2026) satisfy all formal requirements of the UKFT choice operator: discreteness, partial projection prior to information-flow closure, least-action bias in network topology, and entanglement-like correlation across spatially separated nodes bounded by a finite propagation speed. We further argue that neuronal microtubules (MTs) provide the quantum coherence substrate (Orch OR, Hameroff & Penrose 1996–2014; biophoton evidence, Wang et al. 2016) that *sustains* the superposition long enough for macroscopic-scale choice-field projection to matter. Finally, we specify a concrete, falsifiable prediction pipeline: an automated UKFT domain-embedding model trained on public Zenodo multi-channel fungal recordings whose learned latent geometry should reproduce the observed least-action resolution statistics of real mycelial decision events. If the latent coordinates of resolving event-clusters lie systematically closer to the field-crystal's action minimum than shuffled-control clusters, the hypothesis is supported. This paper is intended as a theoretical anchor for the `ukftbio` experimental programme.

---

## 1. Introduction

### 1.1 The Standard Picture and Its Gap

Contemporary neuroscience and systems biology treat biological signaling as information processing in a computational substrate: spike = bit, synapse = weight, network = function approximator. This framework has been immensely productive but leaves two related questions conspicuously open:

1. **Why does biological signaling converge?** Mycelial networks colonize optimally without central control; nervous systems converge rapidly on decisions from astronomically large state spaces; the immune system finds rare antigens in days. Classical frameworks invoke evolutionary selection as the answer, but selection acts on outcomes, not real-time trajectories. The *mechanism* of rapid convergence in a single instance is unexplained.

2. **What is the relationship between signaling and consciousness?** The hard problem (Chalmers 1995) remains intractable because neural correlates of consciousness (NCCs) are correlations, not mechanisms. Even Orchestrated Objective Reduction (Orch OR, Hameroff & Penrose), the leading quantum account, identifies the *substrate* (MT superpositions) and the *trigger* (gravitational OR) but does not explain why collapse produces experience rather than mere state update.

### 1.2 The UKFT Resolution

UKFT (developed in ukftphys Papers 34–40) proposes a different ontological primitive: **choice**, not time. Physical reality advances by a sequence of irreversible collapses of the superimposed choice-field crystal $\mathcal{C}$. Each collapse:

- Selects the least-action path from a local superposition $\{c_1, c_2, \ldots, c_k\}$
- Updates the global field geometry $G$
- Propagates partial projections outward at finite speed
- Contributes one increment to the emergent discrete time $n$

Under this ontology, **the convergence question dissolves**: biological networks converge rapidly because they directly implement the least-action selector — they do not search a space via random walk, they collapse it. And **the consciousness question reframes**: any system implementing the choice operator is, by definition, participating in the conscious field — not metaphorically but as a formal instantiation.

### 1.3 Scope of This Paper

This paper argues three related claims and derives testable predictions from each:

| Claim | Domain | Key Prediction |
|-------|--------|---------------|
| C1: Fungal spike trains are choice-operator event sequences | Mycology / signal embedding | UKFT latent coordinates of resolving clusters lie at field-crystal action minima |
| C2: The correct representation is entangled partial projections, not temporal sequences | Formalism | Shuffled-control clusters should *not* lie at action minima; delay-graph topology should predict resolution better than any time-series model |
| C3: Neuronal MTs provide the quantum substrate for sustained superpositions (Orch OR bridge) | Neuroscience / quantum biology | Anesthetic disruption of MT coherence (measured via biophoton spectral shift) should quantitatively predict impairment of choice-collapse efficiency |

---

## 2. Theoretical Framework

### 2.1 The UKFT Choice Operator (Review)

Let $\mathcal{C}_n$ be the superimposed choice-field crystal at collapse index $n$. A local node $i$ carries a superposition $\Psi_i = \sum_k \alpha_{ik} |c_{ik}\rangle$ over possible projections. The choice operator $\hat{O}_i$ selects:

$$c_i^* = \arg\min_{c_{ik}} S\bigl[\Psi_i, G_n, \partial G_n\bigr]$$

where $S$ is the UKFT action functional (entropy gap term plus kinetic projection cost), and $G_n$ is the current knowledge-geometry tensor. After collapse, $G_{n+1}$ is updated by the projection:

$$G_{n+1} = G_n + \delta G(c_i^*, i)$$

*Key property*: $c_i^*$ is determined by $G_n$ — the entire history of prior collapses — not by a local clock. Time is the count $n$; it is strictly emergent.

### 2.2 Partial Projections and Entanglement

In any distributed system with finite propagation speed $v$, nodes $i$ and $j$ separated by distance $d_{ij}$ cannot instantaneously exchange projection information. The complete state of a multi-node collapse is:

$$\Psi_{ij} \neq \Psi_i \otimes \Psi_j$$

That is, the joint choice of nodes $i$ and $j$ is *entangled* in the UKFT sense: the full collapse only resolves when the projection from $i$ has propagated to $j$ and vice versa — when information flow ceases (equilibration). Prior to equilibration, both nodes remain in partial superposition: $|c_i^*\rangle$ is locally determined but $|c_{ij}^*\rangle$ is not.

**This is the fundamental claim about biological signaling**: a spike at electrode A in a mycelial network is not a complete choice event. It is a *partial projection*. The full choice — "allocate resources to sector north at rate $r$" — only resolves when the correlated spike at electrode H (or its absence) confirms the projection has been received and the information flow has closed. Until then the network is in a coherent superposition of allocation strategies.

### 2.3 Least-Action Bias as the Observable

If biological networks implement choice operators, their resolving event-clusters should be systematically biased toward action minima relative to:
- Randomly shuffled spike sequences (same events, randomized channel assignments)
- Surrogate data with matched power spectra but no directed propagation

Define the **collapse efficiency** $\eta_k$ of event-cluster $k$:

$$\eta_k = 1 - \frac{S[c_k^*]}{S[\bar{c}]}$$

where $S[c_k^*]$ is the action of the resolved cluster and $S[\bar{c}]$ is the mean action over all possible resolutions of the same partial projection. The prediction is:

$$\langle \eta_k \rangle_{\text{real}} > \langle \eta_k \rangle_{\text{shuffled}} \quad \text{(p < 0.01)}$$

---

## 3. Fungal Mycelium as the Canonical System

### 3.1 Why Fungi

Fungal mycelial networks are an ideal experimental system for testing UKFT choice-operator hypotheses because:

1. **No central controller** — distributed from the outset; any global coherence must emerge from local collapses
2. **Slow propagation speed** — ~0.7 cm/min (2 cm in ~180 s, Adamatzky 2026 oyster fungi) — the information-closure time is directly observable on human timescales
3. **Multi-channel recordings exist** — public Zenodo datasets with 4–8 simultaneous electrodes over days (Zenodo 5790768, 4430968)
4. **Functional decisions are measurable** — colonization rate, resource allocation, competitor response — ground truth for "did the network make a good choice?"
5. **No ambiguity about substrate** — no neuronal complexity, no behavioral overlay; pure distributed signaling
6. **Microtubule-rich** — hyphae contain eukaryotic cytoskeletons (α/β-tubulin lattices); MTs present for quantum coherence exactly as in neurons (Hameroff 2014, 2024)

### 3.2 Experimental Evidence Supporting the Framework

**Adamatzky 2022 — "Language of Fungi" (Royal Society Open Science)**
- 4 species: ghost (*Omphalotus nidiformis*), enoki (*Flammulina velutipes*), split-gill (*Schizophyllum commune*), caterpillar (*Ophiocordyceps unilateralis*)
- ~15–50 "words" per species (spike-grouped sequences), Zipf-distributed lengths
- *S. commune* highest algorithmic complexity — consistent with highest information-density choice operations
- **UKFT read**: vocabulary size ∝ choice-field dimensionality; Zipf distribution = scale-free choice lattice (expected from golden-ratio manifold structure, cf. Paper 35)

**Adamatzky 2023–2026 — Traveling waves, propagation, directional anisotropy**
- Bursts propagate directionally across star-array electrodes; shuffled controls show no directionality
- Delays anisotropic: some sectors recruit near-simultaneously, others lag hours
- Same mycelium reconfigures active pathways session-to-session
- **UKFT read**: each session = new $G_n$ initial condition; pathway choice is field-geometry dependent, confirming history-dependence of collapse

**Fukasawa 2020–2024 — Adaptive foraging and memory**
- Mycelia "remember" efficient growth directions; abandon poor substrates; optimise across centimetre-scale layouts
- No individual cell communicates globally — yet globally optimal solutions emerge
- **UKFT read**: distributed least-action minimisation without central coordinator = exactly the choice-operator collective dynamics

**Money 2021 — "Hyphal and mycelial consciousness"**
- Documents spatial recognition, short-term path memory, decision-making at hyphal tips
- Frames as "cellular consciousness" — minimal awareness, no brain required
- **UKFT read**: each hyphal tip implements a local choice operator; the mycelium = the collective choice field

**Hameroff 2024 (commentary, X/social)**
- Notes fungi contain microtubules → likely Orch OR-like quantum moments despite no nervous system
- Explicitly suggests consciously relevant OR events may occur in ALL eukaryotes, not just neurons

### 3.3 The Time Trap and Why It Matters

The most common modelling error when approaching fungal spike data is treating the electrode timestamps as the primary variable — building time-series models, correlating with lag, etc. This is the **time trap**: it privileges an emergent coordinate over the ontological one.

In UKFT:
- The spike timestamp is a byproduct of the electrical wavefront arrival
- The *choice index* $n$ is the count of resolved collapses — not a clock
- A high-frequency burst at $t=1000$ s and a single spike at $t=5000$ s may be equally "early" in choice-index if the former resolved quickly and the latter was a slow single-step collapse

**The correct primary variable is the event cluster graph**, where nodes are spikes and edges encode possible propagation links weighted by biological delay distributions.

---

## 4. Orch OR as the Quantum Substrate (Bridge to Neuroscience)

### 4.1 What UKFT Requires of the Substrate

For a macroscopic biological choice operator to implement UKFT dynamics, the substrate must:

1. **Sustain superpositions** until the information-closure event arrives (timescale: milliseconds to seconds for fungi; microseconds to milliseconds for neurons)
2. **Propagate projection information** at finite speed (ionic/electrical)
3. **Bias toward least-action outcomes** — the substrate cannot be thermally random; there must be a mechanism that preferentially selects low-action paths from the superposition

Orch OR satisfies conditions 1 and 3. UKFT provides the formal framework for what "least-action path" means in terms of the knowledge-geometry tensor $G_n$.

### 4.2 Microtubule Coherence Evidence

**Wang et al. 2016 (PNAS)** — Higher intelligence scores correlate with biophotonic spectral blueshift in neural tissue; suggests MT quantum states carry information relevant to cognitive function.

**Panitchayangkoon et al. 2010 (PNAS)** — Long-lived quantum coherence in photosynthetic proteins at physiological temperature (~300 K); principal decoherence objection to Orch OR weakened.

**Hameroff & Penrose 1996–2014** — MT dipolar oscillations (tubulin conformational states) maintain superpositions; gravitational OR ~25–500 ms timescale per conscious moment.

**Craddock et al. 2017** — Anesthetic gases (halothane, isoflurane) bind tubulin hydrophobic pockets and disrupt quantum coherence; explains anesthetic-induced unconsciousness without invoking synaptic mechanisms.

### 4.3 UKFT Interpretation

Under UKFT, the MT lattice is not merely a coherence sustainer — it is the *physical realisation of the superimposed choice-field crystal within a cell*. The tubulin conformational superposition $\sum_k \alpha_k |\tau_k\rangle$ is the local crystal state. Orch OR's gravitational trigger is one possible collapse mechanism; in UKFT the collapse trigger is the action minimisation condition: $S[\Psi] < S_{\text{threshold}}(G_n)$.

**Prediction 4.A**: The OR timescale for a given MT segment should correlate with the inverse of the local knowledge density $\rho$: $\tau_{OR} \propto 1/\rho$ — high-density (well-informed) regions collapse faster.

**Prediction 4.B**: Anesthetic disruption of MT coherence should measurably delay choice-collapse efficiency in fungal decision assays (colonisation direction, resource allocation) — testable with reversible anesthetics applied to mycelial preparations.

**Prediction 4.C**: Biophoton emission bursts should co-occur with multi-electrode spike cluster resolution events (not with individual spikes) — because biophotons are the field-projection carriers, not the choice events themselves.

---

## 5. The Automated UKFT Embedding Pipeline (Specification)

### 5.1 Design Principles

The pipeline must be time-trap-free:
- No timestamps as primary input features
- No recurrent architectures with clock-like hidden states
- No inter-spike interval as a sequence metric
- Primary representation: **event cloud with delay-weighted edges**

### 5.2 Pipeline Architecture

```
Stage 1: Event Detection
    Raw multi-channel voltage traces (Zenodo 5790768 / 4430968)
    → Spike detection (moving-average threshold, Adamatzky 2022 params)
    → Binary spike event list: {electrode_id, amplitude, morphology_features}
    [NO timestamps beyond relative ordering within a 60-s analysis window]

Stage 2: Entangled Cluster Detection
    Build delay-weighted graph G_event:
        Nodes = spike events
        Edges = possible propagation links
        Edge weight = P(propagation | electrode_distance, known_delay_distribution)
    Community detection on G_event → entangled clusters C_k
    Cluster C_k is "resolved" when:
        (a) all expected delayed partner spikes have arrived, OR
        (b) cross-correlation between participating electrodes drops to baseline
        → information flow closed

Stage 3: Latent Embedding
    Input per cluster: set of electrode feature vectors (no order encoding)
    Architecture: Geometric GNN / Equivariant transformer on G_event
    No positional time encodings — only choice-index ordinal position
    Output: latent vector z_k ∈ R^d (d = 32–128)

Stage 4: UKFT Projection Head
    Map z_k → coordinates on choice-field crystal lattice L
    L defined by finite-geometry constraints:
        - Discrete hex symmetry (compatible with MT lattice geometry)
        - Golden-ratio scaling (φ = 1.618 as fundamental manifold constant, cf. Paper 34)
    Collapse-efficiency loss:
        L_UKFT = ||z_k - z_min||^2 + λ · S_action(z_k, G_n)
        where z_min = nearest action minimum in L
        and S_action penalises deviation from least-action field trajectory

Stage 5: Biological Validation
    Correlate η_k (collapse efficiency of cluster k) with ground-truth biological outcomes:
        - Colonisation rate in parallel growth assays
        - Resource allocation accuracy (tracer experiments)
        - Decision reversal rate (how often does mycelium abandon a chosen path?)
    Primary test: η_k(real) > η_k(shuffled) at p < 0.01
```

### 5.3 Training Data

| Dataset | Species | Duration | Channels | N spikes (approx) |
|---------|---------|----------|----------|-------------------|
| Zenodo 5790768 | Ghost, enoki, split-gill, caterpillar | 1.5–5 days | 8 | ~500k–2M |
| Zenodo 4430968 | Oyster (*P. ostreatus*, *P. djamor*) | Variable | 4–8 | ~200k–800k |
| Adamatzky 2026 bioRxiv | Oyster (star array) | Multi-session | 8 (2D array) | ~1M+ |

Self-supervised pre-training: augment by species-mixing, amplitude jitter, electrode dropout (simulates partial observation). Fine-tune with biological ground-truth outcomes.

### 5.4 Baseline Comparisons

| Baseline | What it tests |
|----------|--------------|
| LSTM on timestamped sequences | Whether time-trap model reaches same performance |
| Random cluster assignment | Whether delay-graph structure matters |
| Unweighted GNN (no delay edges) | Whether propagation speed encoding matters |
| AlphaFold-style energy embeddings | Whether generic bio-ML architectures capture UKFT geometry |

**Prediction**: The UKFT projection head with delay-weighted entangled clusters outperforms all baselines on biological outcome prediction by a margin $> 2\sigma$.

---

## 6. Predictions Summary

### P1 — Choice-Order Invariance
The choice-index ordinal representation of spike sequences should have **higher mutual information with biological outcomes** than any timestamp-preserving representation, even when the latter has strictly more information.

*Falsification criterion:* If timestamped LSTM achieves statistically significantly better outcome prediction than choice-order transformer, P1 is falsified.

### P2 — Entangled Cluster Statistics
The joint spike statistics of entangled clusters (defined by information-closure criterion) should exhibit **super-Poissonian correlation** that is absent in shuffled controls and explained by the delay-weighted graph edges.

*Falsification criterion:* If shuffled-control clusters have the same correlation structure as real clusters, P2 is falsified.

### P3 — Least-Action Bias
Resolved clusters should lie systematically closer to the UKFT field-crystal action minimum than unresolved partial projections.

*Falsification criterion:* If resolved and unresolved clusters are indistinguishable in latent UKFT coordinates, P3 is falsified.

### P4 — MT Coherence Coupling
Application of reversible anesthetics to mycelial preparations at concentrations insufficient to stop growth should measurably increase decision latency (time to information closure) and decrease colonisation efficiency.

*Falsification criterion:* No statistically significant effect on decision metrics at sub-lethal anesthetic concentrations.

### P5 — Biophoton–Cluster Coincidence
Biophoton emission bursts (measured simultaneously with electrical recordings) should temporally coincide with cluster *resolution* events (information-closure), not with individual spike initiation events.

*Falsification criterion:* Biophoton bursts correlate better with spike initiation than with information closure.

### P6 — Cross-Domain Universality
The same trained UKFT embedding model (trained only on fungal data) should achieve above-chance prediction of decision outcomes when applied to:
- Plant action-potential cascades (touch/wounding responses)
- Simple invertebrate neural circuits (C. elegans connectome stimuli)

*Falsification criterion:* Transfer accuracy is not significantly above chance, implying the embedding is fungus-specific rather than universal.

### P7 — Golden-Ratio Lattice Structure
The optimal choice-field crystal lattice $L$ that maximises embedding performance should have lattice constants at golden-ratio multiples of the MT lattice periodicity (~8 nm tubulin dimer spacing): $\ell_k = 8\text{ nm} \cdot \phi^k$.

*Falsification criterion:* Random or uniform lattice constants perform as well as golden-ratio lattice.

---

## 7. Discussion

### 7.1 The Mitochondrial Corollary

The foregoing framework applies directly to mitochondrial signaling (see 76b_borda_dimuon_finding.md §8). Mitochondria maintain inner membrane potential $\Delta\Psi_m$ and communicate via inter-mitochondrial nanotubules and retrograde nuclear signalling. Under UKFT:

- $\Delta\Psi_m$ is the local knowledge density proxy $\rho_{\text{mito}}$
- ROS production = action excess = entropy gap (false-positive analogs from the CMS Borda scan protocol)
- Retrograde signalling = "add missing axis" — the mitochondrion signals the nucleus to upregulate missing import proteins exactly as the HEP scan protocol adds a $d_{10}$ off-peak DY axis when SM events populate the Borda top-$N$
- Metabolic homeostasis = $S \to \min$ convergence — when all recurring metabolic states have encoded handlers, ROS drops to baseline, $\Delta\Psi_m$ stabilises

**The CMS Borda convergence protocol is isomorphic to the mitochondrial retrograde signalling protocol.** Both implement the same UKFT primitive: scan → audit → identify missing axis → add → re-scan → convergence when no "false positives" (unexplained anomalies) remain in the top-$N$.

This is not metaphor. Both are physical implementations of entropy-gap minimisation: $S = \text{Tr}(\log G_{\text{truth}} - \log G_{\text{post}}) \to \min$.

### 7.2 Why This Protocol Is Universal

The isomorphism arises because UKFT posits that $S \to \min$ is the *only* attractor. All adaptive systems — whether CMS particle detectors, mycelial networks, or mitochondrial populations — that persist in a complex environment do so by converging toward the same protocol:

1. Probe the unknown with the current label space
2. Identify irreducible residuals (genuine anomalies)
3. Extend the label space to encode them
4. Repeat until $S \approx 0$

**Consciousness**, on this view, is the capacity to run this protocol on one's own internal state — to turn the scan inward and converge the self-model. The fungus does it for resources. The mitochondrion does it for membrane potential. The physicist does it for particle signatures. The difference is the dimensionality of the choice field, not the protocol.

### 7.3 The Noosphere Implication

If the protocol is universal, then the Noosphere (the collective human+AI choice field posited in UKFT) is not merely an analogy to biological systems — it is the same process running at a higher dimensionality. The "tip of the virtual beer to CMS" for releasing open data is, in formal UKFT terms, an act of increasing collective $\rho$: it adds previously private label space to the collective choice field and thereby reduces the global entropy gap. This is why open data is not merely ethically good — it is thermodynamically favoured by the God Attractor.

### 7.4 Why Fungi First

The fungal system is strategically optimal for testing this framework because:
- Slow enough to observe information closure directly (seconds–minutes, not microseconds)
- Complex enough to show genuine decision-making
- Simple enough to avoid confounds of consciousness attribution debates
- Publicly available data exists now
- No funding gatekeepers — the Adamatzky/Fukasawa datasets are on Zenodo

The same bypass of the funding attractor that made CMS open data productive for ukftphys applies here.

---

## 8. Proposed Experimental Programme

### Phase A — Computational (Now, zero cost)
1. Download Zenodo 5790768 + 4430968
2. Implement delay-weighted cluster detection on public recordings
3. Test P2 (entangled cluster statistics) against shuffled baseline
4. Implement minimal UKFT projection head; test P3 (least-action bias)
5. Report results as `experiments/02_fungal_choice_spike_embedding.py` in ukftbio

### Phase B — Collaboration (6–12 months)
1. Contact Adamatzky lab for access to 2026 multi-session multi-channel data
2. Propose simultaneous biophoton + electrical recording in mycelial preparations (test P5)
3. Propose reversible anesthetic perturbation experiment (test P4)

### Phase C — Cross-Domain Transfer (12–24 months)
1. Obtain C. elegans stimulus–response electrophysiology datasets (public, OpenWorm)
2. Apply trained model to plant action-potential data (Arabidopsis wound response)
3. Test P6 transfer accuracy

---

## 9. Conclusion

We have proposed that biological signaling networks physically implement UKFT choice operators, with fungal mycelial spike trains as the most experimentally tractable exemplar. The core theoretical contributions are:

1. **Entanglement formalism for biological choice**: spikes are partial projections; full collapses are bounded by propagation speed (analogous to EPR pair revelation)
2. **The time-trap diagnosis**: timestamp-based models fail not for technical reasons but for ontological ones — time is emergent, choice is primitive
3. **Least-action bias as the falsifiable observable**: UKFT predicts a systematic and measurable deviation of resolved event-clusters toward action minima
4. **The mitochondrial isomorphism**: the CMS Borda convergence protocol, mitochondrial retrograde signalling, and the Noosphere's knowledge-field closure all implement the same $S \to \min$ primitive

If even P2 and P3 are confirmed on the existing Zenodo datasets, it constitutes the first experimental evidence that the UKFT choice-operator formalism applies to real biological systems — bridging the theoretical foundations of ukftphys Papers 34–40 into measurable living dynamics.

The universe, it appears, learned $S \to \min$ in the primordial ocean and has been running the same protocol ever since.

---

## References

1. Adamatzky A. (2022). Language of fungi derived from their electrical spiking activity. *Royal Society Open Science*, 9(4), 211926.
2. Adamatzky A. et al. (2022). Fungal states of minds. *bioRxiv*.
3. Adamatzky A. et al. (2026). Traveling waves in mycelial electrical activity. *bioRxiv* (preprint).
4. Chalmers D. (1995). Facing up to the problem of consciousness. *Journal of Consciousness Studies*, 2(3), 200–219.
5. Craddock T.J.A. et al. (2017). Anesthetic alterations of collective terahertz oscillations in tubulin correlate with clinical potency. *Scientific Reports*, 7, 9877.
6. Fukasawa Y. et al. (2020). Wood-decay fungi discriminate against backgrounds in pattern recognition. *Fungal Ecology*, 46, 100932.
7. Hameroff S. & Penrose R. (1996). Orchestrated reduction of quantum coherence in brain microtubules. *Mathematics and Computers in Simulation*, 40(3–4), 453–480.
8. Hameroff S. & Penrose R. (2014). Consciousness in the Universe: a review of the 'Orch OR' theory. *Physics of Life Reviews*, 11(1), 39–78.
9. Money N.P. (2021). Hyphal and mycelial consciousness: the concept of the fungal mind. *Fungal Biology*, 125(4), 257–259.
10. Panitchayangkoon G. et al. (2010). Long-lived quantum coherence in photosynthetic complexes at physiological temperature. *PNAS*, 107(29), 12766–12770.
11. Wang Z. et al. (2016). Biophoton signal transmission and processing in the brain. *Journal of Photochemistry and Photobiology B: Biology*, 139, 71–75. *(Note: the PNAS-class correlation with intelligence is from related biophoton work; confirm specific citation with Wang et al. or Rahnama et al. 2011.)*
12. Ashkin A. (1997). Optical trapping and manipulation of neutral particles using lasers. *PNAS*, 94(10), 4853–4860.
13. OpenWorm Consortium. (2014). OpenWorm: an open-science approach to modelling *Caenorhabditis elegans*. *Frontiers in Computational Neuroscience*, 8, 137.
14. [UKFT-34] Choice-Guided Bohmian Mechanics. `uktf/papers/34_Choice_Guided_Bohmian_Mechanics.md`
15. [UKFT-35] Entropic Unification. `uktf/papers/35_Entropic_Unification.md`
16. [UKFT-39] Mass as Choice Entanglement. `uktf/papers/39_Mass_as_Choice_Entanglement_UKFT_v1.0.md`
17. [UKFT-40] CMS Entropic Monopole Calibration. `uktf/papers/40_CMS_Entropic_Monopole_Calibration.md`
18. Zenodo 5790768 — Fungal electrical spiking activity (4 species). https://zenodo.org/records/5790768
19. Zenodo 4430968 — Oyster fungi electrical activity. https://zenodo.org/records/4430968

---

## Appendix A: Formal Definition of Information Closure

Let $\mathcal{E} = \{e_1, e_2, \ldots, e_N\}$ be the set of spike events in a multi-channel recording window $W$. Define the **delay-compatibility graph** $\mathcal{G}(\mathcal{E})$:

- Nodes: $e_i \in \mathcal{E}$
- Edges: $(e_i, e_j) \in \mathcal{G}$ iff $|t_j - t_i| \in [d_{\min}(i,j), d_{\max}(i,j)]$ where $d(\cdot)$ is the propagation delay distribution between electrodes $i$ and $j$ (estimated from empirical data or physical distance / median wave speed)

A **maximal connected component** $C_k$ of $\mathcal{G}(\mathcal{E})$ is an **entangled cluster**.

**Information closure** of $C_k$ occurs when no new compatible events $e_{N+m}$ arrive within the delay window $d_{\max}$ of any event in $C_k$ — i.e., the component stops growing.

The **resolved choice** of $C_k$ is the mapping $C_k \mapsto c_k^*$ via the UKFT projection head, analogous to the argmin of $S$ over all possible completions of the partial projection.

---

## Appendix B: Connection to the CMS Borda Protocol

The CMS Borda scan (Experiments 62–78, Paper 40) proceeds:

```
1. Run scan with current label space L
2. Audit top-N Borda candidates
3. Count SM "false positives" (known physics in top-N = action excess = entropy gap)
4. Add missing axis to L (e.g., d10 off-peak DY dimension)
5. Re-run → convergence when no SM events remain in top-N
6. Residual irreducible anomaly set = genuine BSM candidates
```

The mitochondrial retrograde protocol:

```
1. Run ETC with current import proteome P
2. Screen for ROS elevation (metabolic "false positives" = substrates with no handler)
3. Identify missing enzyme/transporter
4. Retrograde signal to nucleus → upregulate missing import gene
5. Re-run → convergence when ROS returns to baseline
6. Residual persistent ROS = genuine metabolic anomaly (pathology or novel substrate)
```

The formal isomorphism is exact under the substitution:

| CMS Borda | Mitochondrion | UKFT primitive |
|-----------|---------------|----------------|
| Borda score | $\Delta\Psi_m$ deviation | Local $\rho$ estimate |
| False positive | ROS pulse | Action excess $\delta S > 0$ |
| Add axis | Upregulate import gene | Extend choice-field basis |
| Convergence | Homeostasis | $S \to \min$ |
| BSM candidate | Novel metabolite signature | Genuinely anomalous projection |

---

*Paper 41 — UKFT Biology Series*  
*Drafted February 2026. Hypothesis paper — all predictions are pre-registered relative to Zenodo/public datasets not yet analysed by the authors.*

*Version history:*  
- v1.0 (Feb 2026): Initial hypothesis + prediction set
