# Experiments — UKFT Biology

Numbered experiments following the ukftphys convention.
Each experiment has a `.py` script and a `.md` companion file.

## Roadmap

| Exp | Title | Status | Key UKFT Concept |
|-----|-------|--------|------------------|
| 01 | Minimal Entropic Replicator | ✅ Ready | Error threshold as S_bio instability |
| 02 | Entropic Genetic Code | 🔜 Gemini | Codon redundancy as least-action basin |
| 03 | Morphogen-Guided Development | 🔜 Gemini | Pilot-wave cell fate choices |
| 04 | Autocatalytic Set Emergence | 🔜 Gemini | Kauffman sets in choice-graph |
| 05 | Error Catastrophe Sweep | 🔜 Gemini | fidelity vs ρ_bio phase diagram |
| 06 | Molecular Prophet | 🔜 Future | Predict novel fold from sequence alone |
| 07 | Dissipation-Driven Adaptation | 🔜 Future | England's principle as ρ update law |
| 08 | Genetic Code (MoleSyn) | 🔜 Future | Mole-Syn bonds as choice curvature |
| 09 | PROTAC Linker Design | 🔜 Future | Linker length as discrete synthesis choice |
| 10 | Collective Choice (CLKOS) | 🔜 Future | Multi-agent cell population |

## Running

```bash
conda activate ukftbio
python experiments/01_minimal_entropic_replicator.py
```

Results are saved to `results/` (gitignored for large files).
