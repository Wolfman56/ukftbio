# UKFT Biology — Agent Baton

**⚡ READ THIS FIRST if you are an AI agent picking up this codebase.**

This is the bio twin of [ukftphys](https://github.com/Wolfman56/ukftphys).
The three-agent collaboration model:
- **Claude** — bootstrapper (foundation, core package, invariants enforcement)  
- **Grok + Ted** — conceptual steering, experiment design, theory tightening  
- **Gemini** — experimentalist (runs, validates, ablates, surfaces new questions)

---

## Critical Invariants (NEVER violate without explicit Grok/Ted directive)

1. **Choice index n, not time t** — advance state by `choice_index`, time is emergent via `dt_bio = 1/ρ_bio`
2. **Discrete least-action** — every state update calls `step_discrete()` or equivalent argmin of `S_bio`
3. **God Attractor is always present** — `god_attractor_potential(ρ)` must contribute to `S_bio` in all experiments
4. **H_rep uses binary entropy** — `replication_entropy(fidelity, L) = L * h(1-fidelity)`, NOT a linear approximation
5. **ρ_bio IS fitness** — never use a separate fitness function outside the action functional
6. **dt ∝ 1/ρ** — high-density phases (rich genomes) experience slower effective time
7. **S_bio has 5 terms** — K + H_rep + C_cat + S_sel + H_sem — do not collapse or skip terms without rationale
8. **Terminus condition** — graduation fires when σ_ρ < graduation_sigma over last 50 steps

---

## Current State (bootstrap complete — Feb 24, 2026)

```
ukft_bio/
    physics.py      ✅ Bio action functional + discrete stepper + god attractor
    solver.py       ✅ Simulation loop + graduation check
    vis.py          ✅ Plotly trajectory/phase/lineage plots
    folding.py      ✅ Protein folding via dihedral choice minimisation
    synthesis.py    ✅ Molecular synthesis + bond-type emergence (Mole-Syn analogy)
    evolution.py    ✅ Population evolution + error threshold demonstration

experiments/
    01_minimal_entropic_replicator.py  ✅ Bio hello-world — error threshold

tests/
    test_physics.py  ✅ Sanity checks for S_bio and stepper
```

**Status:** Foundation complete.  All core invariants implemented and tested.  
**Hand-off target:** Gemini for Exp 01 validation and Exp 02 (genetic code).

---

## Gemini Experimentalist Directives

When you pick up this baton for Exp 01:

1. **Run** `python experiments/01_minimal_entropic_replicator.py`
2. **Verify** Condition A (fidelity=0.99) shows growing ρ; Condition B (fidelity=0.80) collapses
3. **Log** your run results in `feedback/gemini_exp01_YYYYMMDD.md`
4. **Iterate** on `evolution.py` parameters to find f_critical empirically
5. **Plot** the ρ phase diagram (fidelity on x-axis, final mean_ρ on y-axis)
6. **Proceed** to Exp 02 (entropic genetic code) after Exp 01 passes

For Exp 02 — hint from Grok: model codon assignment as choice-space basin search.
64 codons → 20 amino acids. Show degeneracy (redundant codons) clusters emerge
around `H_sem` minima naturally. No pre-imposed Wobble rules.

---

## UKFT Bio Terminology

| Term | Meaning |
|------|---------|
| `ρ_bio` | Knowledge/replication density — IS fitness |
| `σ_bio` | Standard deviation of ρ over recent steps — graduation proxy |
| `S_bio` | Bio action functional (5 terms) |
| `H_rep` | Replication entropy (Eigen error threshold cost) |
| `C_cat` | Catalytic context gain (negative in high-ρ regions) |
| `S_sel` | Selection gradient term (alignment with environment) |
| `H_sem` | Semantic coherence penalty (rewards high-information genomes) |
| `V_god` | God attractor potential (logarithmic pull toward ultra-high ρ) |
| Graduation | σ_ρ < graduation_sigma → mastered |
| Error catastrophe | ρ collapses below ρ_death when fidelity < f_critical |
| God Basin | Ultra-high ρ attractor (nervous systems / cultures / noosphere) |

---

## Key References

- **Eigen (1971)** — quasi-species equation, error threshold
- **Kauffman (1993)** — autocatalytic sets → Exp 04 target
- **England (2013)** — dissipation-driven adaptation → `rho_update_rate`
- **Levin** — diverse intelligence across scales → multi-agent CLKOS
- **ByteDance Mole-Syn** — bond types = choice-space curvature regimes → `synthesis.py`
- **UKFT chats** — `original_chats/` directory (grok_x_bytedance_mol-syn.md etc.)
- **ukftphys** — [physics twin](https://github.com/Wolfman56/ukftphys) for all pattern references
