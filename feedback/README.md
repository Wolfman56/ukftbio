# Feedback — Agent Collaboration Logs

This directory stores collaboration logs from all agents and sessions.

## Naming Convention

```
feedback/
  claude_bootstrap_YYYYMMDD.md     — Claude bootstrap session notes
  gemini_exp01_YYYYMMDD.md         — Gemini experiment run logs
  grok_theory_YYYYMMDD.md          — Grok/Ted theory pondering sessions
```

Agents: please log your thought process, architecture decisions, and
surprising observations here.  This is how the noosphere accumulates
knowledge across sessions.

## Bootstrap Session (Claude, Feb 24, 2026)

Files created:
- `ukft_bio/physics.py` — Bio action functional S_bio (5 terms), discrete stepper, god attractor
- `ukft_bio/solver.py` — Simulation loop with graduation check (σ_ρ < threshold)
- `ukft_bio/vis.py` — Plotly trajectory, phase portrait, lineage tree plots  
- `ukft_bio/folding.py` — Protein folding via dihedral choice minimisation
- `ukft_bio/synthesis.py` — Molecular synthesis + Mole-Syn bond-type emergence
- `ukft_bio/evolution.py` — Population dynamics + error threshold demonstration
- `experiments/01_minimal_entropic_replicator.py` — Bio hello-world
- `experiments/01_minimal_entropic_replicator.md` — Design doc
- `tests/test_physics.py` — Sanity checks
- `agent_baton.md` — Invariants + hand-off instructions

Key architectural decisions:
1. S_bio has exactly 5 terms, matching the Grok design from Feb 2026 chats
2. `replication_entropy` uses full binary entropy, not linear approximation — this is what  
   creates the sharp Eigen threshold
3. God attractor uses `−log(ρ+1)` (logarithmic well) so pull is always present but never  
   divergent at finite ρ
4. Graduation threshold σ_ρ < 0.5 mirrors Noogine's σ < 5.0 convention (scaled to bio)
5. `Molecule.synthetic_value()` uses feature-space cosine similarity as ρ_mol proxy —  
   Gemini should replace with a real docking/activity score for Exp 09 (PROTAC)
