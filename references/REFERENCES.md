# UKFTBIO — Primary Reference Library

> **Updated**: 2026-02-27  
> **Scope**: Papers and datasets directly used in myco-explorer pipeline  
> Classification labels: A = primary data source, B = theoretical foundation, C = methodology parallel, D = ancillary context

---

## 1. Fungal Electrical Activity (Primary Data Sources)

### [A1] Adamatzky 2026 — Directional Spiking, Oyster Mycelium
```
arXiv:2601.08099 [cs.ET]
"Directional Electrical Spiking, Bursting, and Information Propagation
 in Oyster Mycelium Recorded with a Star-Shaped Electrode Array"
Andrew Adamatzky
Submitted: 13 Jan 2026
DOI: 10.48550/arXiv.2601.08099
```
**Role**: PRIMARY dataset source. 8-channel star-shaped differential electrode array on
Pleurotus ostreatus; directional heterogeneity, clustered bursting, propagation delays
seconds→hours. Direct input for Exp 06-07 (graph construction, choice-indexed model).  
**UKFT mapping**: Each burst origin = choice event Ĉ; propagation delay = emergent dt field.

---

### [A2] Adamatzky 2021 — Language of Fungi
```
arXiv:2112.09907 [cs.ET]
"Language of fungi derived from electrical spiking activity"
Andrew Adamatzky
Submitted: 18 Dec 2021
Published: Royal Society Open Science, 2022
DOI: 10.48550/arXiv.2112.09907
```
**Role**: Establishes spike-train vocabulary framework. Omphalotus nidiformis, Flammulina
velutipes, Schizophyllum commune, Cordyceps militaris. Spikes → words → sentences. Lempel-Ziv
complexity hierarchy. S. commune generates most complex sentences.  
**UKFT mapping**: Word-length distribution ≈ Shannon entropy of choice-space; complexity
hierarchy → rho gradient map.  
**Data files**: Spike traces available from author supplementary / Zenodo (see A1 DOI).

---

### [A3] Adamatzky 2026 — Fungal Systems for Security & Resilience
```
arXiv:2602.10543 [cs.ET]
"Fungal systems for security and resilience"
Andrew Adamatzky
Submitted: 11 Feb 2026
DOI: 10.48550/arXiv.2602.10543
```
**Role**: Framing paper for mycelium as distributed sensing / anomaly-detection substrate.
Decentralised control, embodied memory, autonomous repair → infrastructure resilience.  
**UKFT mapping**: Anomaly-detection function = UKFT Borda anomaly scoring on spike-field;
self-repair = choice-space re-routing after perturbation (P2 test case).

---

### [A4] Adamatzky et al. 2021 — Fungal Electronics
```
arXiv:2111.11231 [cs.ET]
"Fungal electronics"
Andrew Adamatzky, Phil Ayres, Alexander E. Beasley, Alessandro Chiolerio,
Mohammad M. Dehshibi, Antoni Gandia, ..., Han A.B. Wösten
Submitted: 22 Nov 2021
DOI: 10.48550/arXiv.2111.11231
```
**Role**: Foundational survey of fungal electronic devices — impedance changes, potential
spikes, sensing capabilities. Establishes memristive and capacitative properties of mycelium.  
**UKFT mapping**: Memristive behaviour = pilot-wave memory; impedance = local rho.

---

### [A5] Adamatzky et al. 2021 — Logics in Fungal Mycelium Networks
```
arXiv:2112.07236 [cs.ET]
"Logics in fungal mycelium networks"
Andrew Adamatzky, Phil Ayres, Alexander E. Beasley, Nic Roberts,
Martin Tegelaar, Michail-Antisthenis Tsompanas, Han A.B. Wösten
Submitted: 14 Dec 2021
Journal: Logica Universalis (2022)
DOI: 10.48550/arXiv.2112.07236
```
**Role**: Demonstrates Boolean gate implementation in living mycelium via electrical signal
propagation. Provides formal basis for treating mycelium as a computing substrate.  
**UKFT mapping**: Gate logic = discrete choice outcomes; propagation patterns = I-field.

---

### [A6] Przyczyna, Szacilowski, Chiolerio, Adamatzky 2022 — Frequency Discrimination
```
arXiv:2210.01775 [cs.ET]
"Electrical frequency discrimination by fungi Pleurotus ostreatus"
Dawid Przyczyna, Konrad Szacilowski, Alessandro Chiolerio, Andrew Adamatzky
Submitted: 4 Oct 2022
DOI: 10.48550/arXiv.2210.01775
```
**Role**: Fungal mycelium networks discriminate input frequencies (fuzzy/threshold-based).
Key for Exp 08 (resonance layer) — fungal frequency response is the biological analogue of
the Brown/ISF phonon-photon interconversion modes.

---

### [A7] Mayne, Roberts, Phillips, Weerasekera, Adamatzky 2023 — Signal Propagation
```
arXiv:2304.10675 [cs.ET]
"Propagation of electrical signals by fungi"
Richard Mayne, Nic Roberts, Neil Phillips, Roshan Weerasekera, Andrew Adamatzky
Submitted: 20 Apr 2023
DOI: 10.48550/arXiv.2304.10675
```
**Role**: Frequency-modulated signal transmission through mycelium-bound composites. Establishes
propagation constants and attenuation — needed for FMM wavefront speed calibration (Exp 10).

---

### [A8] Gunji, Adamatzky, Mougkogiannis, Khrenikov 2025 — Quantum-like Coherence
```
arXiv:2509.01021 [cs.AI, nlin.AO]
"Quantum-like Coherence Derived from the Interaction between
 Chemical Reaction and Its Environment"
Yukio-Pegio Gunji, Andrew Adamatzky, Panagiotis Mougkogiannis, Andrei Khrenikov
Submitted: 31 Aug 2025
DOI: 10.48550/arXiv.2509.01021
```
**Role**: Defines "open computing" as contrast to classical closed AI. Closed vs open
computing in chemical reaction systems — quantum-like coherence emerges from environment
coupling, not isolation. Direct foundation for Exp 08 biophoton coherence layer.  
**UKFT mapping**: Open computing = UKFT choice operator acting on entropic environment.

---

## 2. Quantum Biology & Time Crystals (Theoretical Foundations)

### [B1] Tsompanas et al. 2022 — Polyatomic Time Crystals of Brain Neuron (J. Appl. Phys.)
```
Journal: Journal of Applied Physics 132, 194401 (2022)
"Polyatomic time crystals of the brain neuron cytoskeleton,
 and their relevance to the hard problem of consciousness"
M.A. Tsompanas, A. Adamatzky, J.A. Tuszynski, et al.
DOI: 10.1063/5.0102521
arXiv: NOT AVAILABLE — journal-only
```
**Role**: Demonstrates polyatomic time-crystal behaviour in extracted microtubule nanowires
via dielectric resonance. Nested MHz–GHz–THz oscillations (fractal structure). Basis for
Exp 09 (fractal time-crystal analysis in mycelium).  
**UKFT mapping**: Fractal resonance nesting ↔ φ=1.618 manifold structure (noogine Phase 6).
dt ∝ 1/ρ: time dilation at high-rho nodes parallels MT slow-time zones.

---

### [B2] Hameroff 2026 — Microtubules are Fractal Time Crystals (ResearchGate preprint)
```
ResearchGate preprint: "Microtubules are 'Fractal Time Crystals':
 Implications for Life and Consciousness"
Stuart Hameroff
Published: ~Feb 2026, ResearchGate
URL: https://www.researchgate.net/
arXiv: NOT AVAILABLE
```
**Role**: Argues MT time-crystal properties extend from nanowires to intact neurons.
Collective entanglement across neuronal networks → binding of "conscious now".
Source chat: grok_x_bio_time_crystal_chat.md  
**UKFT mapping**: Hameroff ladder (MT → neuron → brain → network) = UKFT geo/bio/noo/theo.

---

### [B3] Brown / ISF 2026 — Where Biology Meets Resonance (Non-arXiv)
```
Publication: Harmonic Science Perspectives, Vol 1, Issue 1 (2026)
"Where Biology Meets Resonance: Light, Vibration, and Living Order"
William Brown (International Space Federation)
URL: https://spacefed.com/biology/where-biology-meets-resonance-light-vibration-and-living-order
ResearchGate parallel: "The Quantum-Resonance Symphony of Life:
 The Triad of Photon-Phonon-Exciton Interconversion as an
 Engine of Cellular Organization" (~Jan 2026)
arXiv: NOT AVAILABLE — ISF/journal only
```
**Role**: Photon–phonon–exciton resonance network in cells. Microtubule waveguides for
energy routing. Phase-locking across cellular structures. Testable predictions: ultrafast
spectroscopy of coordinated biophoton timing. Basis for Exp 08.  
**UKFT mapping**: Resonance hierarchy = entropic landscape modulation; phase-locking ↔
pilot-wave coherence; biophoton emission = I-field carrier.

---

## 3. Quantum Gravity / Causal Networks (Theoretical Parallel)

### [C1] Loll 2019 — CDT Review
```
arXiv:1905.08669 [hep-th, gr-qc, hep-lat]
"Quantum Gravity from Causal Dynamical Triangulations: A Review"
R. Loll
Submitted: 21 May 2019
Published: Classical and Quantum Gravity (2020)
DOI: 10.1088/1361-6382/ab57c7
```
**Role**: Authoritative CDT review. Spacetime as sum over causal triangulations. Time and
geometry emerge from discrete causal structure. 64 pages, 28 figures.  
**UKFT mapping**: CDT causal ordering ↔ UKFT choice-indexed time. CDT "foliation" ↔
choice-ancestry tree. Key context for grok_x_causal_networks_chat.md.

---

### [C2] Loll, Ambjorn, Jurkiewicz 2006 — The Universe from Scratch
```
arXiv:hep-th/0509010
"The Universe from Scratch"
R. Loll, J. Ambjorn, J. Jurkiewicz
Submitted: 1 Sep 2005; v2 Oct 2006
Published: Contemporary Physics 47 (2006) 103-117
DOI: 10.1080/00107510600603344
```
**Role**: Accessible CDT overview. "Why spacetime has 4 dimensions" from first principles.
Good paper-citation anchor for UKFT vs CDT comparison discussion.

---

## 4. Non-Riemannian Geometry (Methodology Parallel)

### [C3] Bujack et al. 2025 — Non-Riemannian Color Geometry (CGF/EuroVis)
```
Journal: Computer Graphics Forum (EuroVis 2025)
"The Geometry of Color in the Light of a Non-Riemannian Space"
Roxana Bujack et al. (LANL)
Published: EuroVis 2025
arXiv: NOT AVAILABLE — CGF journal only
Note: Completing Schrödinger's 1920s color geometry; neutral axis =
      geodesic-closest point to black in perceptual metric space
```
**Role**: Structural parallel to UKFT P2 validation. The non-Riemannian completion of an
incomplete classical model by finding the "missing central structure" via geodesic distance —
exactly the operation for finding action-minimizing paths in the fungal graph that Euclidean
shortest-path misses. Methodology template for Exp 16.

---

## 5. Datasets (No arXiv — Direct Download)

### [D1] Adamatzky 2018 — Fungal Computer Supplementary Data (Zenodo)
```
Zenodo: https://zenodo.org/records/1451496
"Supplementary materials for: Adamatzky A. 2018 Towards fungal computer.
 Interface Focus 8: 20180029"
Andrew Adamatzky
DOI: 10.1098/rsfs.2018.0029
```
**Role**: Electrical trace data for multiple fungal species. Historical baseline for comparison
with newer 2026 star-array recordings (A1).

---

### [D2] Adamatzky 2026 — Oyster Spike-Train Data (to be confirmed)
```
Expected: Zenodo or author page (companion to arXiv:2601.08099)
Status: CHECK supplementary of 2601.08099 for data link
Note: 8-channel NDJSON or CSV spike-train files needed for Exp 06
```
**Action needed**: Email author or check ISTA/UWE repository for raw .csv spike traces.

---

## 6. AI / Agent Architecture (Ancillary Context)

### [D3] van Dijk et al. — CellType CLI (GitHub / preprint)
```
GitHub: https://github.com/celltypeinc/celltype-cli
Organization: CellType Inc. (YC-backed, Yale/van Dijk lab)
Context: grok_x_celltype_noobio_2_chat.md
arXiv: NOT AVAILABLE at time of writing
```
**Role**: "Claude Code for biology" — AI agent for drug discovery. PROTAC linker design
as benchmark task. Comparison target for Exp 15 (hyphal routing vs PROTAC linker optimization).

---

## 7. Summary Table

| ID | arXiv / DOI | Authors | Year | Classification |
|----|------------|---------|------|----------------|
| A1 | arXiv:2601.08099 | Adamatzky | 2026 | Primary data |
| A2 | arXiv:2112.09907 | Adamatzky | 2021 | Primary data |
| A3 | arXiv:2602.10543 | Adamatzky | 2026 | Primary data |
| A4 | arXiv:2111.11231 | Adamatzky et al. | 2021 | Primary data |
| A5 | arXiv:2112.07236 | Adamatzky et al. | 2021 | Primary data |
| A6 | arXiv:2210.01775 | Przyczyna et al. | 2022 | Primary data |
| A7 | arXiv:2304.10675 | Mayne et al. | 2023 | Primary data |
| A8 | arXiv:2509.01021 | Gunji et al. | 2025 | Primary data |
| B1 | DOI:10.1063/5.0102521 | Tsompanas et al. | 2022 | Theory (MT/time crystal) |
| B2 | ResearchGate preprint | Hameroff | 2026 | Theory (MT/time crystal) |
| B3 | ISF/spacefed.com | Brown | 2026 | Theory (biophoton coherence) |
| C1 | arXiv:1905.08669 | Loll | 2019 | Method parallel (CDT) |
| C2 | arXiv:hep-th/0509010 | Loll et al. | 2006 | Method parallel (CDT) |
| C3 | CGF/EuroVis 2025 | Bujack et al. | 2025 | Method parallel (non-Riemannian) |
| D1 | Zenodo:1451496 | Adamatzky | 2018 | Dataset |
| D2 | TBD | Adamatzky | 2026 | Dataset (needed) |
| D3 | GitHub/celltypeinc | van Dijk et al. | 2026 | Benchmark comparison |
