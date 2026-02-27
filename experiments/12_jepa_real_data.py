#!/usr/bin/env python3
"""
Experiment 12 — JEPA Temporal Surprise on Real Pleurotus Data
==============================================================
Joint-Embedding Predictive Architecture (JEPA) applied to myco-explorer
Phase-4 aligned embeddings from Pleurotus ostreatus Star_Lab03.

THEORY
------
A JEPA predictor f(z_t) → ẑ_{t+1} that operates in the geodesic-aligned
embedding space should encode the *expected* next state under minimum-action
dynamics.  Surprise S = 1 - cos(ẑ_{t+1}, z_{t+1}) is the *epistemic residual*:
the gap the predictor cannot close.

UKFT prediction:
  • High-ρ windows are INFORMATION BURST events: the system departs from its
    attractor, collapses into a new state (choice event), and produces entropy.
    These are prospectively SURPRISING — the next state cannot be predicted from
    the previous because the choice is irreversible.
  • Low-ρ windows (silent/baseline) form the geodesic attractor basin:
    highly predictable (the system returns to near-zero activity).
  Therefore: rs(S, ρ) > 0, p < 0.05   ← Gate H12a

  This is the UKFT duality:
    ρ measures information production retrospectively (spikes arrived)
    S measures information production prospectively (state was unpredictable)
    Their positive correlation is the signature of irreversible choice collapse.

Secondary analyses:
  • H12b: rs(S, borda_score) — borda scores high for windows with unusual
          multi-scale structure; expect positive (anomalies are surprising)
  • H12c: rs(S, fmm_score) > 0 (FMM scores high for frontier windows which are
          ambiguous/surprising)

ARCHITECTURE
------------
  JEPAPredictor: 40 → 64 (tanh) → 40 (L2-normalised)
  Loss: 1 - cosine(ẑ_{t+1}, z_{t+1})   (per-pair, mean over batch)
  Optimizer: SGD + momentum (m=0.9), cosine LR schedule
  Sequences: per-channel (8 channels), temporal order within channel
  k = 1 (predict next consecutive window in the same electrode pair)

GATE CRITERIA (Phase 4 ukftbio)
  H12a: rs(surprise, rho) > 0   AND p < 0.05  ← MUST PASS
        (information-burst events are irreversibly surprising — choice collapse)
  H12b: rs(surprise, borda)     result         ← REPORT
  H12c: rs(surprise, fmm)       result         ← REPORT
  H12d: test cosine accuracy (fraction where cosine(ẑ, z) > 0.5) ≥ 0.40
        (above chance, meaning the predictor learned the attractor manifold)

Usage:
    python3 experiments/12_jepa_real_data.py [--epochs N] [--dry-run]
"""

import argparse
import json
import math
import random
import sys
from pathlib import Path

import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).parent.parent
MYCO_ROOT = REPO_ROOT.parent / "noosphere" / "apps" / "myco-explorer"
DATA_DIR = MYCO_ROOT / "tools" / "data"
RESULTS_DIR = MYCO_ROOT / "results"

ALIGNED_NDJSON = DATA_DIR / "pleurotus_spike_aligned.ndjson"
BORDA_JSON = RESULTS_DIR / "borda_rankings.json"
FMM_JSON = RESULTS_DIR / "fmm_scores.json"

OUT_DIR = REPO_ROOT / "experiments" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_REPORT = OUT_DIR / "12_jepa_report.json"
OUT_SURPRISE = OUT_DIR / "12_jepa_surprise_scores.ndjson"

# ── Hyperparameters ──────────────────────────────────────────────────────────
HIDDEN_DIM = 64
K_AHEAD = 1           # predict k windows ahead within same channel
LR_INIT = 0.02
LR_MIN = 1e-4
MOMENTUM = 0.9
WEIGHT_DECAY = 1e-4
BATCH_SIZE = 64
VAL_FRAC = 0.20
SEED = 42


# ── Tiny pure-numpy JEPA predictor ───────────────────────────────────────────

class JEPAPredictor:
    """
    Single hidden-layer temporal predictor w/ L2-normalised output.
    Maps z_t (40D) → ẑ_{t+k} (40D, unit-norm).

    Parameters stored as numpy arrays for portability (no torch dependency).
    """

    def __init__(self, in_dim: int = 40, hidden_dim: int = HIDDEN_DIM):
        rng = np.random.default_rng(SEED)
        scale1 = math.sqrt(2.0 / in_dim)
        scale2 = math.sqrt(2.0 / hidden_dim)
        self.W1 = rng.normal(0, scale1, (hidden_dim, in_dim)).astype(np.float32)
        self.b1 = np.zeros(hidden_dim, dtype=np.float32)
        self.W2 = rng.normal(0, scale2, (in_dim, hidden_dim)).astype(np.float32)
        self.b2 = np.zeros(in_dim, dtype=np.float32)
        # Momentum buffers
        self.mW1 = np.zeros_like(self.W1)
        self.mb1 = np.zeros_like(self.b1)
        self.mW2 = np.zeros_like(self.W2)
        self.mb2 = np.zeros_like(self.b2)

    def forward(self, x: np.ndarray) -> tuple[np.ndarray, dict]:
        """x: (N, 40) → z_hat: (N, 40), cache for backward."""
        h = np.tanh(x @ self.W1.T + self.b1)        # (N, H)
        out = h @ self.W2.T + self.b2                # (N, 40)
        # L2-normalise
        norms = np.linalg.norm(out, axis=1, keepdims=True) + 1e-8
        z_hat = out / norms
        cache = {"x": x, "h": h, "out": out, "norms": norms}
        return z_hat, cache

    def loss_grad(self, x: np.ndarray, z_target: np.ndarray) -> tuple[float, dict]:
        """
        Cosine loss: L = mean(1 - dot(ẑ, z_target))
        Returns (loss, {dW1, db1, dW2, db2}).
        """
        z_hat, cache = self.forward(x)
        # Loss = 1 - mean cosine (targets already normalised)
        cos_vals = np.sum(z_hat * z_target, axis=1)   # (N,)
        loss = float(1.0 - np.mean(cos_vals))

        # Backward through normalisation
        N = x.shape[0]
        # dL/d(z_hat) = -z_target / N
        dz_hat = -z_target / N                         # (N, 40)
        # dL/d(out): chain through L2-norm
        h, out, norms = cache["h"], cache["out"], cache["norms"]
        # z_hat = out / norms  →  d(z_hat)/d(out) = (I - z_hat z_hat^T) / norms
        dot = np.sum(dz_hat * z_hat, axis=1, keepdims=True)
        dout = (dz_hat - dot * z_hat) / norms          # (N, 40)

        # dL/dW2, dL/db2
        dW2 = dout.T @ h                               # (40, H)
        db2 = dout.sum(axis=0)                         # (40,)
        dh = dout @ self.W2                            # (N, H)

        # tanh backward
        dh_pre = dh * (1.0 - h ** 2)                  # (N, H)
        dW1 = dh_pre.T @ x                             # (H, 40)
        db1 = dh_pre.sum(axis=0)                       # (H,)

        grads = {"dW1": dW1, "db1": db1, "dW2": dW2, "db2": db2}
        return loss, grads

    def step(self, grads: dict, lr: float):
        """SGD + momentum update with L2 weight decay."""
        for name, g in grads.items():
            param_name = name[1:]  # strip leading 'd'
            param = getattr(self, param_name)
            mom = getattr(self, "m" + param_name)
            # Add L2 decay to weights (not biases)
            if param_name.startswith("W"):
                g = g + WEIGHT_DECAY * param
            mom[:] = MOMENTUM * mom + g
            param -= lr * mom


# ── Data helpers ─────────────────────────────────────────────────────────────

def load_aligned(path: Path) -> list[dict]:
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def build_sequences(records: list[dict], k: int = K_AHEAD) -> dict[str, list[dict]]:
    """
    Group records by channel, sort by t_start_s.
    Returns {channel_id: [record, ...]} sorted temporally.
    """
    channels: dict[str, list] = {}
    for r in records:
        ch = str(r["channel"])
        channels.setdefault(ch, []).append(r)
    for ch in channels:
        channels[ch].sort(key=lambda r: r["t_start_s"])
    return channels


def make_pairs(sequences: dict, k: int = K_AHEAD) -> list[tuple[dict, dict]]:
    """
    Build (context_t, target_{t+k}) pairs across all channels.
    Each element is a (context_record, target_record) tuple.
    """
    pairs = []
    for ch, seq in sequences.items():
        for i in range(len(seq) - k):
            pairs.append((seq[i], seq[i + k]))
    return pairs


def cosine_lr(epoch: int, total: int) -> float:
    """Cosine annealing from LR_INIT to LR_MIN."""
    if total <= 1:
        return LR_INIT
    prog = epoch / (total - 1)
    return LR_MIN + 0.5 * (LR_INIT - LR_MIN) * (1.0 + math.cos(math.pi * prog))


def spearman_rs(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Spearman rank correlation + p-value (t-distribution approx)."""
    n = len(a)
    ra = np.argsort(np.argsort(a)).astype(np.float64)
    rb = np.argsort(np.argsort(b)).astype(np.float64)
    ra -= ra.mean(); rb -= rb.mean()
    denom = (np.sqrt(np.sum(ra**2)) * np.sqrt(np.sum(rb**2)))
    rs = float(np.sum(ra * rb) / denom) if denom > 0 else 0.0
    # t-test approximation
    if abs(rs) >= 1.0:
        p = 0.0
    else:
        t_stat = rs * math.sqrt((n - 2) / (1.0 - rs**2))
        # Two-tailed p using normal approx for large n
        from math import erfc, sqrt
        p = float(erfc(abs(t_stat) / sqrt(2.0)))
    return rs, p


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Exp 12: JEPA on real Pleurotus")
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--dry-run", action="store_true",
                        help="Run 5 epochs, skip save")
    args = parser.parse_args()

    epochs = 5 if args.dry_run else args.epochs
    print(f"{'[DRY-RUN] ' if args.dry_run else ''}Exp 12 JEPA — {epochs} epochs")
    print(f"  Aligned embeddings: {ALIGNED_NDJSON}")

    # ── Load data ─────────────────────────────────────────────────────────────
    if not ALIGNED_NDJSON.exists():
        sys.exit(f"ERROR: {ALIGNED_NDJSON} not found — run apply_projection.py first")

    print("Loading aligned embeddings...", end=" ", flush=True)
    records = load_aligned(ALIGNED_NDJSON)
    print(f"{len(records)} windows")

    # Build window→rho map
    wid_to_rho = {r["window_id"]: r["rho"] for r in records}
    wid_to_sigma = {r["window_id"]: r["sigma"] for r in records}

    # ── Build temporal sequences ──────────────────────────────────────────────
    sequences = build_sequences(records, K_AHEAD)
    print(f"  Channels: {sorted(sequences.keys())}")
    ch_lengths = {ch: len(seq) for ch, seq in sequences.items()}
    total_windows = sum(ch_lengths.values())
    print(f"  Windows per channel: {ch_lengths}")
    print(f"  Total: {total_windows} windows")

    pairs = make_pairs(sequences, K_AHEAD)
    print(f"  Training pairs (k={K_AHEAD}): {len(pairs)}")

    if len(pairs) < 20:
        sys.exit("ERROR: too few training pairs — check data")

    # Train/val split: last VAL_FRAC of each channel's sequence (temporal)
    train_pairs, val_pairs = [], []
    # Re-build splits per-channel to respect temporal ordering
    for ch, seq in sequences.items():
        n = len(seq)
        split = max(1, int(n * (1.0 - VAL_FRAC)))
        for i in range(n - K_AHEAD):
            if i < split - K_AHEAD:
                train_pairs.append((seq[i], seq[i + K_AHEAD]))
            else:
                val_pairs.append((seq[i], seq[i + K_AHEAD]))

    print(f"  Train pairs: {len(train_pairs)}  |  Val pairs: {len(val_pairs)}")

    def pairs_to_arrays(pairs):
        X = np.array([p[0]["aligned_embedding"] for p in pairs], dtype=np.float32)
        Y = np.array([p[1]["aligned_embedding"] for p in pairs], dtype=np.float32)
        return X, Y

    X_train, Y_train = pairs_to_arrays(train_pairs)
    X_val, Y_val = pairs_to_arrays(val_pairs)

    # ── Train JEPA predictor ──────────────────────────────────────────────────
    model = JEPAPredictor(in_dim=40, hidden_dim=HIDDEN_DIM)
    rng = random.Random(SEED)

    best_val_loss = float("inf")
    best_W1 = model.W1.copy()
    best_b1 = model.b1.copy()
    best_W2 = model.W2.copy()
    best_b2 = model.b2.copy()

    print(f"\nTraining JEPA (40→{HIDDEN_DIM}→40, cosine loss):")
    n_train = len(X_train)
    indices = list(range(n_train))

    for epoch in range(epochs):
        lr = cosine_lr(epoch, epochs)

        # Shuffle training pairs
        rng.shuffle(indices)

        train_losses = []
        for start in range(0, n_train, BATCH_SIZE):
            batch_idx = indices[start:start + BATCH_SIZE]
            x_b = X_train[batch_idx]
            y_b = Y_train[batch_idx]
            loss, grads = model.loss_grad(x_b, y_b)
            model.step(grads, lr)
            train_losses.append(loss)

        # Validation
        z_hat_val, _ = model.forward(X_val)
        val_cos = float(np.mean(np.sum(z_hat_val * Y_val, axis=1)))
        val_loss = 1.0 - val_cos

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_W1 = model.W1.copy()
            best_b1 = model.b1.copy()
            best_W2 = model.W2.copy()
            best_b2 = model.b2.copy()

        if epoch % max(1, epochs // 10) == 0 or epoch == epochs - 1:
            avg_train = float(np.mean(train_losses))
            print(f"  epoch {epoch:4d}/{epochs}  lr={lr:.5f}  "
                  f"train_loss={avg_train:.4f}  val_loss={val_loss:.4f}")

    # Restore best
    model.W1 = best_W1; model.b1 = best_b1
    model.W2 = best_W2; model.b2 = best_b2

    print(f"\nBest val loss: {best_val_loss:.4f}  (cosine accuracy threshold 0.5)")

    # ── Compute surprise on ALL pairs (train+val) ─────────────────────────────
    print("\nComputing per-window surprise scores...")
    all_pairs = train_pairs + val_pairs
    X_all, Y_all = pairs_to_arrays(all_pairs)

    # Batch inference
    z_hat_all, _ = model.forward(X_all)
    cos_all = np.sum(z_hat_all * Y_all, axis=1)           # (N,)
    surprise_all = 1.0 - cos_all                           # (N,)

    # Cosine accuracy (fraction where cos > 0.5)
    cos_acc = float(np.mean(cos_all > 0.5))
    print(f"  Cosine accuracy (cos>0.5): {cos_acc:.4f}  [gate ≥ 0.40]")

    # Map back to window IDs (use target window as the "surprised" window)
    surprise_records = []
    for pair, s_val in zip(all_pairs, surprise_all):
        ctx, tgt = pair
        surprise_records.append({
            "window_id": tgt["window_id"],
            "context_id": ctx["window_id"],
            "surprise": float(s_val),
            "rho": tgt["rho"],
            "sigma": tgt["sigma"],
            "n_spikes": tgt["n_spikes"],
            "t_start_s": tgt["t_start_s"],
            "channel": tgt["channel"],
        })

    # ── Load borda + fmm for cross-correlation ────────────────────────────────
    borda_data = borda_map = fmm_map = None
    if BORDA_JSON.exists():
        borda_list = json.loads(BORDA_JSON.read_text())
        borda_map = {r["window_id"]: r["borda_score"] for r in borda_list}
    if FMM_JSON.exists():
        fmm_raw = json.loads(FMM_JSON.read_text())
        if isinstance(fmm_raw, dict):
            fmm_map = {wid: v["fmm_score"] for wid, v in fmm_raw.items()}
        else:
            fmm_map = {r["window_id"]: r["fmm_score"] for r in fmm_raw}

    # ── Gate analysis ─────────────────────────────────────────────────────────
    print("\nGate analysis:")

    # H12a: surprise vs rho
    rho_vals = np.array([r["rho"] for r in surprise_records], dtype=np.float64)
    surp_vals = np.array([r["surprise"] for r in surprise_records], dtype=np.float64)

    rs_rho, p_rho = spearman_rs(surp_vals, rho_vals)
    h12a_pass = bool(rs_rho > 0 and p_rho < 0.05)
    print(f"  H12a  rs(surprise, ρ) = {rs_rho:+.4f}  p = {p_rho:.2e}  "
          f"→ {'✅ PASS' if h12a_pass else '❌ FAIL'}"
          f"  [expect > 0: burst=choice-collapse=surprise]")

    # H12b: surprise vs borda
    rs_borda = p_borda = None
    if borda_map:
        matched_b = [(r["surprise"], borda_map[r["window_id"]])
                     for r in surprise_records if r["window_id"] in borda_map]
        if matched_b:
            sb, bb = zip(*matched_b)
            rs_borda, p_borda = spearman_rs(np.array(sb), np.array(bb))
            print(f"  H12b  rs(surprise, borda) = {rs_borda:+.4f}  p = {p_borda:.2e}  "
                  f"[expect > 0: anomaly candidates are surprising]")

    # H12c: surprise vs fmm
    rs_fmm = p_fmm = None
    if fmm_map:
        matched_f = [(r["surprise"], fmm_map[r["window_id"]])
                     for r in surprise_records if r["window_id"] in fmm_map]
        if matched_f:
            sf, ff = zip(*matched_f)
            rs_fmm, p_fmm = spearman_rs(np.array(sf), np.array(ff))
            print(f"  H12c  rs(surprise, fmm)   = {rs_fmm:+.4f}  p = {p_fmm:.2e}  "
                  f"[expect > 0: frontier windows are surprising]")

    # H12d: cosine accuracy
    h12d_pass = bool(cos_acc >= 0.40)
    print(f"  H12d  cosine accuracy      = {cos_acc:.4f}  "
          f"→ {'✅ PASS' if h12d_pass else '❌ FAIL'}")

    # ── Descriptive statistics ────────────────────────────────────────────────
    print(f"\nSurprise distribution:")
    print(f"  mean = {float(np.mean(surp_vals)):.4f}")
    print(f"  median = {float(np.median(surp_vals)):.4f}")
    print(f"  std = {float(np.std(surp_vals)):.4f}")
    print(f"  p5 = {float(np.percentile(surp_vals, 5)):.4f}  "
          f"p95 = {float(np.percentile(surp_vals, 95)):.4f}")

    # Top-5 most surprising windows
    top5_idx = np.argsort(surp_vals)[-5:][::-1]
    print(f"\nTop-5 most surprising windows:")
    for i in top5_idx:
        r = surprise_records[i]
        print(f"  {r['window_id']}  surprise={r['surprise']:.4f}  "
              f"ρ={r['rho']:.3f}  ch={r['channel']}  "
              f"t={r['t_start_s']:.0f}s  n_spikes={r['n_spikes']}")

    # Top-5 most predictable (least surprising)
    bot5_idx = np.argsort(surp_vals)[:5]
    print(f"\nTop-5 most predictable windows (lowest surprise):")
    for i in bot5_idx:
        r = surprise_records[i]
        print(f"  {r['window_id']}  surprise={r['surprise']:.4f}  "
              f"ρ={r['rho']:.3f}  ch={r['channel']}  "
              f"t={r['t_start_s']:.0f}s  n_spikes={r['n_spikes']}")

    # ── GATE SUMMARY ──────────────────────────────────────────────────────────
    gate_open = h12a_pass and h12d_pass
    print(f"\n{'='*60}")
    print(f"GATE SUMMARY")
    print(f"  H12a  rs(surprise, ρ) > 0, p < 0.05  : {'✅ PASS' if h12a_pass else '❌ FAIL'}")
    print(f"  H12d  cosine accuracy ≥ 0.40          : {'✅ PASS' if h12d_pass else '❌ FAIL'}")
    print(f"\n{'✅ PHASE 4 GATE OPEN' if gate_open else '❌ GATE BLOCKED'}")
    print(f"{'='*60}")

    if gate_open:
        print("\nINSIGHT: UKFT duality between ρ and surprise confirmed.")
        print("  ρ = retrospective entropy production (spikes arrived, measured)")
        print("  S = prospective entropy production (state was unpredictable, JEPA)")
        print("  rs(S, ρ) > 0 is the UKFT signature of irreversible choice collapse:")
        print("  the same windows that ARE informative WERE also unpredictable.")
        print("  Silent attractor basin (ρ≈0): self-predicting, high cos-acc.")
        print("  Burst events (ρ>>0): irreducibly surprising = genuine choice events.")
        print("  This is the minimum-action principle operating in REVERSE:")
        print("  it tells you WHAT happened (ρ) but not WHAT WILL happen (S).")

    # ── Save ──────────────────────────────────────────────────────────────────
    if not args.dry_run:
        # Save per-window surprise scores as NDJSON
        with open(OUT_SURPRISE, "w") as f:
            for r in surprise_records:
                f.write(json.dumps(r) + "\n")
        print(f"\nSurprise scores: {OUT_SURPRISE}  ({len(surprise_records)} records)")

        # Save report
        report = {
            "experiment": "12_jepa_real_data",
            "description": "JEPA temporal predictor: surprise vs rho on aligned Pleurotus embeddings",
            "config": {
                "epochs": epochs,
                "k_ahead": K_AHEAD,
                "hidden_dim": HIDDEN_DIM,
                "lr_init": LR_INIT,
                "lr_min": LR_MIN,
                "momentum": MOMENTUM,
                "weight_decay": WEIGHT_DECAY,
                "batch_size": BATCH_SIZE,
                "val_frac": VAL_FRAC,
            },
            "data": {
                "total_windows": len(records),
                "paired_windows": len(surprise_records),
                "channels": list(sequences.keys()),
            },
            "results": {
                "best_val_loss": float(best_val_loss),
                "cosine_accuracy": float(cos_acc),
                "surprise_mean": float(np.mean(surp_vals)),
                "surprise_median": float(np.median(surp_vals)),
                "surprise_std": float(np.std(surp_vals)),
                "H12a_hypothesis": "rs(surprise, rho) > 0: burst=choice-collapse=surprise",
                "H12a_rs_surprise_rho": float(rs_rho),
                "H12a_p": float(p_rho),
                "H12a_pass": bool(h12a_pass),
                "H12b_rs_surprise_borda": float(rs_borda) if rs_borda is not None else None,
                "H12b_p": float(p_borda) if p_borda is not None else None,
                "H12c_rs_surprise_fmm": float(rs_fmm) if rs_fmm is not None else None,
                "H12c_p": float(p_fmm) if p_fmm is not None else None,
                "H12d_cos_acc": float(cos_acc),
                "H12d_pass": bool(h12d_pass),
                "gate_open": bool(gate_open),
            },
        }
        OUT_REPORT.write_text(json.dumps(report, indent=2))
        print(f"Report: {OUT_REPORT}")
    else:
        print("\n[DRY-RUN] Skipping save. Results:")
        print(f"  H12a rs={rs_rho:+.4f} p={p_rho:.2e} pass={h12a_pass}")
        print(f"  H12d cos_acc={cos_acc:.4f} pass={h12d_pass}")
        print(f"  gate_open={gate_open}")


if __name__ == "__main__":
    main()
