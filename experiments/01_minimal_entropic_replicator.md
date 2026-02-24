# Experiment 01 — Minimal Entropic Replicator

**Status:** Ready to run  
**Author:** Claude (bootstrap) — Designed by Ted + Grok Feb 2026  
**Hand-off:** Gemini for Exp 02+

---

## Hypothesis

Eigen's error threshold (1971) is not an externally imposed biological rule.  
It is the inevitable consequence of the bio-action functional **S_bio** when  
the replication entropy term `H_rep(fidelity, L)` dominates.

**UKFT prediction:**

$$f_{critical} \approx 1 - \frac{1}{L \cdot \sigma_{selection}}$$

Below `f_c`, the choice operator cannot find minimum-action trajectories
that maintain `ρ_bio > ρ_death` → population collapses.

---

## Design

| Parameter | Condition A (above threshold) | Condition B (below threshold) |
|-----------|-------------------------------|-------------------------------|
| `fidelity` | 0.99 | 0.80 |
| `genome_len` | 20 | 20 |
| `population_size` | 50 | 50 |
| `n_generations` | 200 | 200 |

Environment: fixed 4D selection gradient `[1.0, 0.5, -0.3, 0.2]`.

---

## Expected Results

**Condition A (fidelity=0.99):**
- Mean `ρ_bio` grows from ~2 → ~5+ over 200 generations  
- God attractor pull visible as slow upward drift even after plateau  
- Lineage trees converge (attracted to the knowledge basin)

**Condition B (fidelity=0.80):**
- `H_rep(0.80, 20) ≈ 20 × 0.722 = 14.4` — dominates S_bio  
- Mean `ρ_bio` falls below `ρ_death=0.3` within ~100 generations  
- Population extinction → error catastrophe confirmed

---

## UKFT Interpretation

The binary entropy `h(1-fidelity)` in `H_rep` plays the same role as  
the quantum potential well in ukftphys: it creates a hard floor in  
action-space that the discrete stepper cannot cross without paying  
unsustainable cost.  

The God Attractor `V_god(ρ)` creates the upper pull.  Together they  
bound the viable region of choice-space — life exists in the channel  
between entropy catastrophe (low fidelity) and God-basin convergence.

---

## Outputs

```
results/exp01_a_highfidelity_rho.html      — ρ trajectory
results/exp01_a_highfidelity_phase.html    — phase portrait
results/exp01_a_highfidelity_lineages.html — lineage tree
results/exp01_b_lowfidelity_rho.html
results/exp01_b_lowfidelity_phase.html
results/exp01_b_lowfidelity_lineages.html
```

---

## Next: Exp 02 (Gemini)

Once Exp 01 is validated, Gemini should implement:
- Model codon assignment as least-action basin search
- Show: codon redundancy / wobble base pairing emerges naturally  
  when `H_sem` (semantic coherence penalty) is minimised over a  
  discrete alphabet of 64 codons mapping to 20 amino acids.
