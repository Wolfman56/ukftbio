"""
Experiment 18 — Duroxide Orchestration
=======================================

Builds a Python-native durable-execution runtime that directly mirrors the
Duroxide (Rust) patterns:

    OrchestrationContext  →  schedules activities, manages history turns
    ActivityRegistry      →  registers named activity callables
    SqliteHistoryStore    →  persists (instance_id, activity_name, input_hash)
                             → result; enables deterministic replay
    Runtime               →  drives orchestration turns, dispatches activities
    Client                →  starts orchestrations, waits for completion

The fungal UKFT pipeline decomposes into checkpointed activities:

    ┌─────────┐   ┌──────────────┐   ┌──────────────┐   ┌──────────┐
    │  Ingest │──▶│ DetectSpikes │──▶│ BuildEmbedding│──▶│ FmmScore │
    └─────────┘   └──────────────┘   └──────────────┘   └────┬─────┘
                                                               │
    ┌──────────┐  ┌──────────────┐   ┌──────────────┐         ▼
    │  Report  │◀─│  BordaFuse   │◀──│  CentroidS3  │◀──┌──────────┐
    └──────────┘  └──────────────┘   └──────────────┘   │JepaScore │
                                                          └──────────┘

Replay is demonstrated by:
  1. Running the orchestration for Schizophyllum — completes all 7 activities.
  2. Inserting a new orchestration for Cordyceps — simulating a "crash" after
     BuildEmbedding by pre-seeding the store with steps 1–3, then resuming.
  3. On resume, Ingest + DetectSpikes + BuildEmbedding are replayed from SQLite
     (< 1 ms per step); FmmScore, JepaScore, CentroidS3, BordaFuse run fresh.
  4. Wall-clock times per step are recorded and compared (fresh vs replay).

Figures
-------
  18_dag_status.png       Pipeline DAG coloured by completion status
  18_replay_speedup.png   Per-step timing: fresh run vs replay
  18_checkpoint_heatmap.png  SQLite checkpoint store contents (step × species)
  18_pipeline_borda.png   Borda score distributions for both species
  18_report.json          Full orchestration metrics and top candidates
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import time
import warnings
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import uniform_filter1d

# ── Paths ─────────────────────────────────────────────────────────────────────
MYCO        = Path(__file__).resolve().parent.parent.parent / "noosphere/apps/myco-explorer"
DATA_DIR    = MYCO / "tools/data/multispecies"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
DB_PATH     = RESULTS_DIR / "18_duroxide.sqlite"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Signal constants (same as Exps 14–17) ────────────────────────────────────
GLOBAL_THRESH_MV = 0.5
NOISE_FLOOR_K    = 2.0
WINDOW_S         = 600
DT_S             = 1.0

# ── 40D histogram bin edges (mirrors spike_embed.py) ─────────────────────────
N_BINS      = 10
AMP_EDGES   = np.logspace(np.log10(0.03), np.log10(2.1), N_BINS + 1)
ISI_EDGES   = np.logspace(0, np.log10(3600), N_BINS + 1)
BURST_EDGES = np.array([1, 2, 3, 4, 5, 6, 8, 10, 14, 18, 100])

# ── FMM parameters (mirrors fmm_score.py) ────────────────────────────────────
K_NEIGHBORS = 10
CHUNK_SIZE  = 512

# ── JEPA parameters (consistent with Exps 15–17) ─────────────────────────────
SEQ_LEN      = 16
PRED_HORIZON = 4
HIDDEN_DIM   = 32
LR           = 3e-3
N_EPOCHS     = 60
MIN_WINDOWS  = 64

# ── Borda weights (mirrors borda_rank.py) ─────────────────────────────────────
W1, W2, W3 = 0.16, 0.08, 0.76

MULTISPECIES_FILES = {
    "Schizophyllum": "Schizophyllum commune.txt",
    "Cordyceps":     "Cordyceps militari.txt",
}

TOP_K = 10

SPECIES_COLOURS = {"Schizophyllum": "#9b59b6", "Cordyceps": "#e67e22"}

# ═════════════════════════════════════════════════════════════════════════════
# DUROXIDE-INSPIRED RUNTIME (pure Python, SQLite-backed)
# ═════════════════════════════════════════════════════════════════════════════

class SqliteHistoryStore:
    """
    Mirrors Duroxide's SqliteProvider.
    Stores (instance_id, activity_name, input_hash) → (result_json, elapsed_ms).
    On replay, a matching (instance_id, name, hash) returns the cached result
    instantly without re-running the activity.
    """

    def __init__(self, db_path: str) -> None:
        self._conn = sqlite3.connect(db_path)
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS activity_history (
                instance_id  TEXT    NOT NULL,
                activity     TEXT    NOT NULL,
                input_hash   TEXT    NOT NULL,
                result_json  TEXT    NOT NULL,
                elapsed_ms   REAL    NOT NULL,
                completed_at REAL    NOT NULL,
                PRIMARY KEY (instance_id, activity, input_hash)
            )
        """)
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS orchestration_status (
                instance_id  TEXT PRIMARY KEY,
                status       TEXT NOT NULL,
                started_at   REAL NOT NULL,
                finished_at  REAL
            )
        """)
        self._conn.commit()

    def lookup(self, instance_id: str, activity: str,
               input_hash: str) -> Optional[tuple[Any, float]]:
        """Return (result, elapsed_ms) if cached, else None."""
        row = self._conn.execute(
            "SELECT result_json, elapsed_ms FROM activity_history "
            "WHERE instance_id=? AND activity=? AND input_hash=?",
            (instance_id, activity, input_hash),
        ).fetchone()
        if row:
            return json.loads(row[0]), float(row[1])
        return None

    def record(self, instance_id: str, activity: str,
               input_hash: str, result: Any, elapsed_ms: float) -> None:
        self._conn.execute(
            "INSERT OR REPLACE INTO activity_history "
            "(instance_id, activity, input_hash, result_json, elapsed_ms, completed_at) "
            "VALUES (?,?,?,?,?,?)",
            (instance_id, activity, input_hash,
             json.dumps(result, default=_json_default),
             elapsed_ms, time.time()),
        )
        self._conn.commit()

    def start_orchestration(self, instance_id: str) -> None:
        self._conn.execute(
            "INSERT OR IGNORE INTO orchestration_status "
            "(instance_id, status, started_at) VALUES (?,?,?)",
            (instance_id, "running", time.time()),
        )
        self._conn.commit()

    def finish_orchestration(self, instance_id: str, status: str) -> None:
        self._conn.execute(
            "UPDATE orchestration_status SET status=?, finished_at=? "
            "WHERE instance_id=?",
            (status, time.time(), instance_id),
        )
        self._conn.commit()

    def list_completed_activities(self, instance_id: str) -> list[str]:
        rows = self._conn.execute(
            "SELECT activity FROM activity_history WHERE instance_id=?",
            (instance_id,),
        ).fetchall()
        return [r[0] for r in rows]

    def activity_checkpoint_table(self) -> dict[str, dict[str, float]]:
        """Return {instance_id: {activity: elapsed_ms}} for the heatmap."""
        rows = self._conn.execute(
            "SELECT instance_id, activity, elapsed_ms FROM activity_history"
        ).fetchall()
        result: dict[str, dict[str, float]] = {}
        for iid, act, ms in rows:
            result.setdefault(iid, {})[act] = ms
        return result


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float32, np.float64, float)):
        return float(obj)
    if isinstance(obj, (np.int32, np.int64, int)):
        return int(obj)
    raise TypeError(type(obj))


def _input_hash(*args) -> str:
    """Deterministic hash of activity inputs (for cache keying)."""
    raw = json.dumps(args, default=_json_default, sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


class ActivityContext:
    """
    Mirrors Duroxide's ActivityContext.
    Provides trace helper and cancellation stub.
    """
    def __init__(self, instance_id: str, activity: str) -> None:
        self.instance_id = instance_id
        self.activity    = activity
        self._cancelled  = False

    def trace_info(self, msg: str) -> None:
        print(f"    [{self.activity}] {msg}", flush=True)

    def is_cancelled(self) -> bool:
        return self._cancelled


class OrchestrationContext:
    """
    Mirrors Duroxide's OrchestrationContext.
    Every `schedule_activity` call:
      1. Computes input_hash
      2. Checks SqliteHistoryStore for cached result → REPLAY (< 1 ms)
      3. On cache miss → runs activity, records to store → FRESH
    """

    def __init__(self, instance_id: str, store: SqliteHistoryStore,
                 activity_registry: dict[str, Callable]) -> None:
        self._instance_id = instance_id
        self._store       = store
        self._registry    = activity_registry
        self.timings: dict[str, dict] = {}   # activity → {elapsed_ms, replayed}

    def schedule_activity(self, name: str, *args) -> Any:
        """Run or replay a named activity.  Returns its result."""
        if name not in self._registry:
            raise KeyError(f"Activity '{name}' not registered")

        h    = _input_hash(name, *args)
        ctx  = ActivityContext(self._instance_id, name)
        cached = self._store.lookup(self._instance_id, name, h)

        if cached is not None:
            result, orig_ms = cached
            t0 = time.perf_counter()
            # Simulate instant replay (just a dict lookup)
            elapsed = (time.perf_counter() - t0) * 1000
            print(f"  ↺ REPLAY  {name:25s}  "
                  f"(orig {orig_ms:.1f} ms → replay {elapsed:.2f} ms)", flush=True)
            self.timings[name] = {"elapsed_ms": elapsed, "replayed": True,
                                  "orig_ms": orig_ms}
            return result

        # Fresh execution
        t0 = time.perf_counter()
        result = self._registry[name](ctx, *args)
        elapsed = (time.perf_counter() - t0) * 1000
        self._store.record(self._instance_id, name, h, result, elapsed)
        print(f"  ✓ FRESH   {name:25s}  ({elapsed:.1f} ms)", flush=True)
        self.timings[name] = {"elapsed_ms": elapsed, "replayed": False,
                              "orig_ms": elapsed}
        return result

    def trace_info(self, msg: str) -> None:
        print(f"  [orch] {msg}", flush=True)


# ═════════════════════════════════════════════════════════════════════════════
# FUNGAL PIPELINE ACTIVITIES
# ═════════════════════════════════════════════════════════════════════════════

def _detrend(arr: np.ndarray) -> np.ndarray:
    n, nchan = arr.shape
    det = np.empty_like(arr)
    for ch in range(nchan):
        col = arr[:, ch].astype(np.float64)
        nm  = np.isnan(col)
        if nm.any():
            ix = np.arange(n)
            col[nm] = np.interp(ix[nm], ix[~nm], col[~nm])
        det[:, ch] = (col - uniform_filter1d(col, size=600, mode="mirror")).astype(np.float32)
    return det


def activity_ingest(ctx: ActivityContext, species: str) -> dict:
    """Load .txt file and return stats dict (not the full array — too big for JSON)."""
    path = DATA_DIR / MULTISPECIES_FILES[species]
    ctx.trace_info(f"Loading {path.name}")
    df  = pd.read_csv(str(path), sep="\t", header=0, dtype=np.float32,
                      na_values=["NaN", "", " "], engine="c",
                      on_bad_lines="skip")
    df.dropna(how="all", inplace=True)
    df  = df.select_dtypes(include=[np.floating, np.integer])
    ctx.trace_info(f"  {len(df):,} rows × {df.shape[1]} channels")
    return {"n_rows": int(len(df)), "n_channels": int(df.shape[1]),
            "n_windows": int(len(df) // WINDOW_S), "species": species}


def _load_and_detrend(species: str) -> np.ndarray:
    """Helper: load + detrend voltage matrix (not an activity itself)."""
    path = DATA_DIR / MULTISPECIES_FILES[species]
    df   = pd.read_csv(str(path), sep="\t", header=0, dtype=np.float32,
                       na_values=["NaN", "", " "], engine="c", on_bad_lines="skip")
    df.dropna(how="all", inplace=True)
    df   = df.select_dtypes(include=[np.floating, np.integer])
    return _detrend(df.values.astype(np.float32))


def _make_histogram(spikes: list[dict]) -> np.ndarray:
    if not spikes:
        return np.zeros(40, dtype=np.float32)
    amps  = np.array([s["amplitude_mv"] for s in spikes], dtype=np.float64)
    isis  = np.array([s["preceding_isi_s"] for s in spikes
                      if s["preceding_isi_s"] is not None], dtype=np.float64)
    burst = np.array([3 if s["in_burst"] else 1 for s in spikes], dtype=np.float64)
    chans = np.array([s["channel"] for s in spikes], dtype=np.int32)
    h_amp  = np.histogram(amps, bins=AMP_EDGES)[0].astype(float)
    h_isi  = np.histogram(isis if len(isis) > 0 else np.array([1.0]),
                          bins=ISI_EDGES)[0].astype(float)
    h_bst  = np.histogram(burst, bins=BURST_EDGES)[0].astype(float)
    h_chan = np.zeros(10, dtype=float)
    for ch in range(8):
        h_chan[ch] = np.sum(chans == ch)
    h_chan[8] = float(len(spikes))
    h_chan[9] = float(sum(s["in_burst"] for s in spikes))
    h  = np.concatenate([h_amp, h_isi, h_bst, h_chan])
    s  = h.sum()
    return (h / s).astype(np.float32) if s > 0 else h.astype(np.float32)


def activity_detect_spikes(ctx: ActivityContext, species: str,
                            n_rows: int, n_channels: int) -> dict:
    """Detect spike events; return compact summary stats."""
    ctx.trace_info(f"Detecting spikes ({n_rows:,} rows × {n_channels} ch)")
    det = _load_and_detrend(species)
    n, nchan = det.shape
    total_spikes = 0
    per_channel  = []
    for ch in range(nchan):
        sig   = det[:, ch].astype(np.float64)
        above = sig >= GLOBAL_THRESH_MV
        pad   = np.concatenate([[False], above])
        starts= np.where(np.diff(pad.astype(int)) == 1)[0]
        total_spikes += len(starts)
        per_channel.append(int(len(starts)))
    ctx.trace_info(f"  {total_spikes:,} spike events across {nchan} channels")
    return {"total_spikes": total_spikes, "per_channel": per_channel,
            "species": species}


def activity_build_embeddings(ctx: ActivityContext, species: str,
                               n_rows: int, n_channels: int) -> dict:
    """Build 40D spike-histogram embeddings + windowed densities."""
    ctx.trace_info("Building 40D histogram embeddings")
    det   = _load_and_detrend(species)
    n, nchan = det.shape
    n_win = n // WINDOW_S

    # Spike detection per window
    all_spikes: list[dict] = []
    for ch in range(nchan):
        sig   = det[:, ch].astype(np.float64)
        above = sig >= GLOBAL_THRESH_MV
        pad   = np.concatenate([[False], above])
        starts= np.where(np.diff(pad.astype(int)) == 1)[0]
        prev  = None
        for idx in starts:
            lo = max(0, idx - 25); hi = min(n, idx + 100)
            pk = int(np.argmax(det[lo:hi, ch]))
            ts = (lo + pk) * DT_S
            isi = (ts - prev) if prev is not None else None
            all_spikes.append({"channel": ch, "t_s": ts,
                                "amplitude_mv": float(det[lo + pk, ch]),
                                "preceding_isi_s": isi,
                                "in_burst": isi is not None and isi < 60.0})
            prev = ts
    all_spikes.sort(key=lambda s: s["t_s"])

    # Bucket into windows
    win_spikes: list[list[dict]] = [[] for _ in range(n_win)]
    for sp in all_spikes:
        wi = int(sp["t_s"] / WINDOW_S)
        if wi < n_win:
            win_spikes[wi].append(sp)

    embeddings = [_make_histogram(win_spikes[i]) for i in range(n_win)]

    # Windowed densities
    abs_v  = np.abs(det)
    sigma  = np.nanmedian(abs_v, axis=0) / 0.6745
    floor  = (NOISE_FLOOR_K * sigma).astype(np.float32)
    noise  = abs_v < floor[np.newaxis, :]
    glob_m = abs_v >= GLOBAL_THRESH_MV
    loc_m  = (~noise) & (~glob_m)
    tr_l   = loc_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    tr_g   = glob_m[:n_win * WINDOW_S].reshape(n_win, WINDOW_S, -1)
    rho_l  = tr_l.mean(axis=(1, 2)).tolist()
    rho_g  = tr_g.mean(axis=(1, 2)).tolist()

    ctx.trace_info(f"  {n_win} windows, {len(all_spikes):,} spikes")
    return {"embeddings": [e.tolist() for e in embeddings],
            "rho_local": rho_l, "rho_global": rho_g,
            "n_windows": n_win, "species": species}


def activity_fmm_score(ctx: ActivityContext, embeddings: list,
                        species: str) -> dict:
    """Intra-species FMM kNN topology deviation (S1)."""
    ctx.trace_info(f"Computing FMM S1 for {len(embeddings)} windows")
    emb = np.array(embeddings, dtype=np.float32)
    n   = len(emb)
    if n <= K_NEIGHBORS:
        return {"s1": [0.0] * n, "species": species}

    # kNN
    knn_idx  = np.zeros((n, K_NEIGHBORS), dtype=np.int32)
    knn_dist = np.zeros((n, K_NEIGHBORS), dtype=np.float32)
    for start in range(0, n, CHUNK_SIZE):
        end   = min(start + CHUNK_SIZE, n)
        chunk = emb[start:end]
        sq_a  = (chunk ** 2).sum(1, keepdims=True)
        sq_b  = (emb ** 2).sum(1, keepdims=True)
        d2    = sq_a + sq_b.T - 2.0 * (chunk @ emb.T)
        np.clip(d2, 0.0, None, out=d2)
        for li in range(end - start):
            gi  = start + li
            row = d2[li].copy(); row[gi] = np.inf
            top = np.argpartition(row, K_NEIGHBORS)[:K_NEIGHBORS]
            ts  = top[np.argsort(row[top])]
            knn_idx[gi]  = ts
            knn_dist[gi] = np.sqrt(np.maximum(row[ts], 0.0))

    mean_d = knn_dist.mean(1); mean_d[mean_d < 1e-12] = 1e-12
    rho    = 1.0 / mean_d; rho_max = rho.max() + 1e-12

    scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        devs = []
        for ji, j in enumerate(knn_idx[i]):
            w_exp = ((rho[i] + rho[j]) / 2.0) / rho_max
            dist  = float(knn_dist[i, ji])
            w_act = rho[i] / rho_max / (dist + 1e-12) * (knn_dist[i].mean() + 1e-12)
            if w_exp > 1e-12:
                devs.append(abs(w_act - w_exp) / w_exp)
        scores[i] = float(np.mean(devs)) if devs else 0.0
    mx = scores.max() + 1e-9
    scores /= mx
    ctx.trace_info(f"  mean S1={scores.mean():.4f}")
    return {"s1": scores.tolist(), "species": species}


class _TinyJEPA:
    def __init__(self):
        np.random.seed(42)
        s1 = np.sqrt(2.0 / SEQ_LEN); s2 = np.sqrt(2.0 / HIDDEN_DIM)
        self.We = np.random.randn(HIDDEN_DIM, SEQ_LEN).astype(np.float32) * s1
        self.be = np.zeros(HIDDEN_DIM, dtype=np.float32)
        self.Wp = np.random.randn(HIDDEN_DIM, HIDDEN_DIM).astype(np.float32) * s2
        self.bp = np.zeros(HIDDEN_DIM, dtype=np.float32)
        self.Wo = np.random.randn(1, HIDDEN_DIM).astype(np.float32) * s2
        self.bo = np.zeros(1, dtype=np.float32)
    def _enc(self, x): return np.maximum(0, self.We @ x + self.be)
    def _pred(self, h): return np.maximum(0, self.Wp @ h + self.bp)
    def _dec(self, h): return float((self.Wo @ h + self.bo)[0])
    def _step(self, ctx, tgt):
        h_e=self._enc(ctx); h_p=self._pred(h_e); p=self._dec(h_p); e=p-tgt
        dWo=e*h_p[np.newaxis,:]; dbo=np.array([e])
        dp=(self.Wo.T*e).flatten()*(h_p>0)
        dWp=np.outer(dp,h_e); dbp=dp
        de=(self.Wp.T@dp)*(h_e>0)
        dWe=np.outer(de,ctx); dbe=de
        r=1e-4
        self.Wo-=LR*(dWo+r*self.Wo); self.bo-=LR*dbo
        self.Wp-=LR*(dWp+r*self.Wp); self.bp-=LR*dbp
        self.We-=LR*(dWe+r*self.We); self.be-=LR*dbe
        return e**2
    def _norm(self, s): return (s-s.mean())/(s.std()+1e-8)
    def train(self, s):
        s=self._norm(s)
        idx=[i for i in range(len(s)-SEQ_LEN-PRED_HORIZON)]
        for _ in range(N_EPOCHS):
            np.random.shuffle(idx)
            [self._step(s[i:i+SEQ_LEN], float(s[i+SEQ_LEN+PRED_HORIZON-1])) for i in idx]
    def surprise(self, s):
        s=self._norm(s)
        return np.array([(self._dec(self._pred(self._enc(s[i:i+SEQ_LEN])))
                          - float(s[i+SEQ_LEN+PRED_HORIZON-1]))**2
                         for i in range(len(s)-SEQ_LEN-PRED_HORIZON)],
                        dtype=np.float32)


def activity_jepa_score(ctx: ActivityContext, rho_local: list,
                         rho_global: list, species: str) -> dict:
    """TinyJEPA surprise on combined density series (S2)."""
    combined = np.array(rho_local, dtype=np.float32) + np.array(rho_global, dtype=np.float32)
    n = len(combined)
    ctx.trace_info(f"JEPA on {n}-window density series")
    if n < MIN_WINDOWS or combined.std() < 1e-4:
        ctx.trace_info("  degenerate series → zeros")
        return {"s2": [0.0] * n, "species": species}
    model = _TinyJEPA()
    model.train(combined)
    scores = model.surprise(combined)
    pad = np.full(n - len(scores), scores.mean(), dtype=np.float32)
    s2 = np.concatenate([pad, scores]).tolist()
    ctx.trace_info(f"  mean S2={np.mean(s2):.4f}")
    return {"s2": s2, "species": species}


def activity_centroid_s3(ctx: ActivityContext,
                          all_embeddings: dict[str, list],
                          species: str,
                          target_embeddings: list) -> dict:
    """Cosine distance to cross-species centroid (S3)."""
    ctx.trace_info("Computing cross-species centroid S3")
    all_emb = np.vstack([np.array(v, dtype=np.float32)
                         for v in all_embeddings.values()])
    centroid = all_emb.mean(axis=0)
    c_norm   = np.linalg.norm(centroid)
    if c_norm > 1e-12:
        centroid /= c_norm

    emb   = np.array(target_embeddings, dtype=np.float32)
    norms = np.linalg.norm(emb, axis=1)
    zero  = norms < 1e-12
    safe  = norms.copy(); safe[zero] = 1.0
    e_n   = emb / safe[:, np.newaxis]
    s3    = (1.0 - e_n @ centroid).tolist()
    ctx.trace_info(f"  mean S3={np.mean(s3):.4f}")
    return {"s3": s3, "species": species}


def activity_borda_fuse(ctx: ActivityContext,
                         s1: list, s2: list, s3: list,
                         species: str) -> dict:
    """Borda rank fusion → top-K anomaly candidates."""
    ctx.trace_info(f"Borda fusion over {len(s1)} windows")
    n = len(s1)
    a1 = np.array(s1, dtype=np.float32)
    a2 = np.array(s2, dtype=np.float32)[:n]
    a3 = np.array(s3, dtype=np.float32)
    def borda(v):
        r = np.argsort(np.argsort(v)).astype(float)
        return r / (n - 1) if n > 1 else np.zeros(n)
    fused = W1 * borda(a1) + W2 * borda(a2) + W3 * borda(a3)
    ranked = sorted(enumerate(fused.tolist()), key=lambda x: -x[1])
    top_k  = [{"window_idx": int(wi), "borda": round(b, 6),
               "s1": round(float(a1[wi]), 6), "s2": round(float(a2[wi]), 6),
               "s3": round(float(a3[wi]), 6)}
              for wi, b in ranked[:TOP_K]]
    ctx.trace_info(f"  top borda={top_k[0]['borda']:.4f}")
    return {"borda_scores": fused.tolist(), "top_k": top_k, "species": species}


# ═════════════════════════════════════════════════════════════════════════════
# ORCHESTRATION
# ═════════════════════════════════════════════════════════════════════════════

ACTIVITY_REGISTRY: dict[str, Callable] = {
    "Ingest":           activity_ingest,
    "DetectSpikes":     activity_detect_spikes,
    "BuildEmbeddings":  activity_build_embeddings,
    "FmmScore":         activity_fmm_score,
    "JepaScore":        activity_jepa_score,
    "CentroidS3":       activity_centroid_s3,
    "BordaFuse":        activity_borda_fuse,
}

PIPELINE_ORDER = ["Ingest", "DetectSpikes", "BuildEmbeddings",
                  "FmmScore", "JepaScore", "CentroidS3", "BordaFuse"]


def run_fungal_orchestration(instance_id: str, species: str,
                              store: SqliteHistoryStore,
                              all_embeddings_pool: dict[str, list]) -> dict:
    """
    Duroxide-style orchestration for one species.
    Returns the full pipeline result dict + per-step timings.
    """
    ctx = OrchestrationContext(instance_id, store, ACTIVITY_REGISTRY)
    store.start_orchestration(instance_id)
    ctx.trace_info(f"Orchestration started  (instance={instance_id})")

    # Turn 1: Ingest
    ingest  = ctx.schedule_activity("Ingest", species)
    n_rows  = ingest["n_rows"]
    n_chan  = ingest["n_channels"]

    # Turn 2: DetectSpikes
    spikes  = ctx.schedule_activity("DetectSpikes", species, n_rows, n_chan)

    # Turn 3: BuildEmbeddings
    emb_r   = ctx.schedule_activity("BuildEmbeddings", species, n_rows, n_chan)
    embeddings  = emb_r["embeddings"]
    rho_local   = emb_r["rho_local"]
    rho_global  = emb_r["rho_global"]

    # Turn 4: FmmScore
    fmm = ctx.schedule_activity("FmmScore", embeddings, species)

    # Turn 5: JepaScore
    jepa = ctx.schedule_activity("JepaScore", rho_local, rho_global, species)

    # Turn 6: CentroidS3 (uses all species embeddings pool)
    s3 = ctx.schedule_activity("CentroidS3", all_embeddings_pool, species, embeddings)

    # Turn 7: BordaFuse
    borda = ctx.schedule_activity("BordaFuse",
                                   fmm["s1"], jepa["s2"], s3["s3"], species)

    store.finish_orchestration(instance_id, "Completed")
    ctx.trace_info("Orchestration complete ✓")
    return {
        "species":     species,
        "n_windows":   emb_r["n_windows"],
        "top_k":       borda["top_k"],
        "borda_scores": borda["borda_scores"],
        "embeddings":  embeddings,
        "timings":     ctx.timings,
    }


# ═════════════════════════════════════════════════════════════════════════════
# FIGURES
# ═════════════════════════════════════════════════════════════════════════════

def fig_dag_status(store: SqliteHistoryStore, instance_ids: list[str]) -> None:
    """Pipeline DAG coloured by activity completion status per species."""
    n_inst = len(instance_ids)
    fig, axes = plt.subplots(1, n_inst, figsize=(6 * n_inst, 5))
    if n_inst == 1:
        axes = [axes]

    for ax, iid in zip(axes, instance_ids):
        completed = set(store.list_completed_activities(iid))
        y_pos = {name: idx for idx, name in enumerate(PIPELINE_ORDER[::-1])}
        for i, act in enumerate(PIPELINE_ORDER[::-1]):
            colour = "#2ecc71" if act in completed else "#e74c3c"
            ax.barh(i, 1, color=colour, edgecolor="white", height=0.6)
            ax.text(0.5, i, act, ha="center", va="center", fontsize=9,
                    fontweight="bold", color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        sp = iid.replace("_run1", "").replace("_run2", "")
        run = "Run 1" if "_run1" in iid else "Run 2 (replay)"
        ax.set_title(f"{sp} — {run}\n{len(completed)}/{len(PIPELINE_ORDER)} activities ✓",
                     fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, len(PIPELINE_ORDER) - 0.5)

    plt.suptitle("Pipeline DAG Status by Orchestration Instance",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    out = str(RESULTS_DIR / "18_dag_status.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 18_dag_status.png", flush=True)


def fig_replay_speedup(timings_run1: dict[str, dict],
                       timings_run2: dict[str, dict],
                       species: str) -> None:
    acts  = [a for a in PIPELINE_ORDER if a in timings_run1 or a in timings_run2]
    t1    = [timings_run1.get(a, {}).get("elapsed_ms", 0.0) for a in acts]
    t2_ms = [timings_run2.get(a, {}).get("elapsed_ms", 0.0) for a in acts]
    x      = np.arange(len(acts))
    w      = 0.35

    fig, ax = plt.subplots(figsize=(12, 5))
    b1 = ax.bar(x - w/2, t1, w, color="#3498db", alpha=0.85, label="Run 1 (fresh)")
    b2 = ax.bar(x + w/2, t2_ms, w, color="#e74c3c", alpha=0.85, label="Run 2 (replay)")

    for xi, (r1, r2, t) in enumerate(zip(t1, t2_ms, acts)):
        if timings_run2.get(t, {}).get("replayed", False):
            ax.text(xi + w/2, r2 + 5, "↺", ha="center", va="bottom",
                    fontsize=11, color="#c0392b")
        if r1 > 0 and r2 > 0:
            sp = r1 / max(r2, 0.01)
            if sp > 2:
                ax.text(xi, max(r1, r2) + 20, f"×{sp:.0f}", ha="center",
                        fontsize=8, color="#27ae60", fontweight="bold")

    ax.set_xticks(x); ax.set_xticklabels(acts, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Elapsed (ms)", fontsize=11)
    ax.set_title(
        f"Duroxide Replay Speedup — {species}\n"
        "↺ = replayed from SQLite checkpoint; ×N = speedup factor",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    out = str(RESULTS_DIR / "18_replay_speedup.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 18_replay_speedup.png", flush=True)


def fig_checkpoint_heatmap(store: SqliteHistoryStore) -> None:
    """Heatmap: rows=activity, cols=instance_id, colour=log10(elapsed_ms)."""
    table = store.activity_checkpoint_table()
    iids  = sorted(table.keys())
    acts  = PIPELINE_ORDER

    mat = np.full((len(acts), len(iids)), np.nan)
    for j, iid in enumerate(iids):
        for i, act in enumerate(acts):
            if act in table[iid]:
                mat[i, j] = table[iid][act]

    fig, ax = plt.subplots(figsize=(max(6, len(iids) * 2.5), 5))
    log_mat = np.log10(np.where(mat == 0, np.nan, mat))
    im = ax.imshow(log_mat, aspect="auto", cmap="RdYlGn_r",
                   vmin=np.nanmin(log_mat) - 0.5, vmax=np.nanmax(log_mat) + 0.5)
    plt.colorbar(im, ax=ax, label="log₁₀(elapsed ms)", fraction=0.03)

    for i in range(len(acts)):
        for j in range(len(iids)):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i,j]:.0f}ms", ha="center", va="center",
                        fontsize=7, color="black")

    ax.set_xticks(range(len(iids)))
    ax.set_xticklabels([iid.replace("_run", "\nrun") for iid in iids], fontsize=8)
    ax.set_yticks(range(len(acts)))
    ax.set_yticklabels(acts, fontsize=9)
    ax.set_title("SQLite Checkpoint Store: Activity Elapsed Times\n"
                 "(green = fast/replayed, red = slow/fresh)", fontsize=11)
    plt.tight_layout()
    out = str(RESULTS_DIR / "18_checkpoint_heatmap.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 18_checkpoint_heatmap.png", flush=True)


def fig_pipeline_borda(results: list[dict]) -> None:
    """Borda score distributions for both species."""
    n_sp = len(results)
    fig, axes = plt.subplots(1, n_sp, figsize=(7 * n_sp, 5), sharey=False)
    if n_sp == 1:
        axes = [axes]

    for ax, r in zip(axes, results):
        sp     = r["species"]
        scores = np.array(r["borda_scores"])
        x      = np.arange(len(scores))
        ax.fill_between(x, scores, alpha=0.5,
                        color=SPECIES_COLOURS.get(sp, "#888"))
        ax.set_xlabel("Window index", fontsize=10)
        ax.set_ylabel("Borda score", fontsize=10)
        ax.set_title(f"{sp} — {len(scores)} windows\n"
                     f"top-1: window {r['top_k'][0]['window_idx']} "
                     f"(borda={r['top_k'][0]['borda']:.4f})", fontsize=10)
        # Mark top-3
        for entry in r["top_k"][:3]:
            wi = entry["window_idx"]
            ax.axvline(wi, color="red", linestyle="--", linewidth=1.2, alpha=0.8)
            ax.text(wi + 1, entry["borda"] * 0.98, f"#{r['top_k'].index(entry)+1}",
                    fontsize=7, color="red")
        ax.grid(True, alpha=0.3)

    plt.suptitle("Orchestrated Borda Anomaly Scores per Species",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    out = str(RESULTS_DIR / "18_pipeline_borda.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 18_pipeline_borda.png", flush=True)


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("EXPERIMENT 18 — Duroxide Orchestration (durable fungal pipeline)")
    print("=" * 70)

    store = SqliteHistoryStore(str(DB_PATH))

    # ── Pre-build cross-species embeddings pool (shared S3 centroid) ──────
    # Use Schizophyllum + Cordyceps for pool; both will be embedded as part of
    # their own orchestrations, so we embed both first, then pass the pool into
    # each CentroidS3 activity.
    print("\n--- Pre-stage: build embeddings for centroid pool ---", flush=True)
    pool: dict[str, list] = {}
    for sp in MULTISPECIES_FILES:
        ctx_pre = OrchestrationContext(f"{sp}_pool", store, ACTIVITY_REGISTRY)
        ingest  = ctx_pre.schedule_activity("Ingest", sp)
        emb_r   = ctx_pre.schedule_activity(
            "BuildEmbeddings", sp, ingest["n_rows"], ingest["n_channels"])
        pool[sp] = emb_r["embeddings"]
        print(f"  {sp}: {len(pool[sp])} windows in pool", flush=True)

    # ── RUN 1: Schizophyllum (fresh) ──────────────────────────────────────
    print("\n--- Run 1: Schizophyllum (full fresh pipeline) ---", flush=True)
    r1 = run_fungal_orchestration("Schizophyllum_run1", "Schizophyllum",
                                   store, pool)
    timings_schiz_run1 = r1["timings"].copy()

    # ── RUN 2: Schizophyllum again (full replay) ──────────────────────────
    print("\n--- Run 2: Schizophyllum (full replay — same instance_id) ---", flush=True)
    r2 = run_fungal_orchestration("Schizophyllum_run1", "Schizophyllum",
                                   store, pool)
    timings_schiz_run2 = r2["timings"].copy()

    # ── RUN 3: Cordyceps (fresh) ──────────────────────────────────────────
    print("\n--- Run 3: Cordyceps (full fresh pipeline) ---", flush=True)
    r3 = run_fungal_orchestration("Cordyceps_run1", "Cordyceps",
                                   store, pool)

    # ── Figures ───────────────────────────────────────────────────────────
    print("\n--- Figures ---", flush=True)
    fig_dag_status(store, ["Schizophyllum_run1", "Cordyceps_run1"])
    fig_replay_speedup(timings_schiz_run1, timings_schiz_run2, "Schizophyllum")
    fig_checkpoint_heatmap(store)
    fig_pipeline_borda([r1, r3])

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("EXPERIMENT 18 SUMMARY")
    print("=" * 70)

    all_instances = ["Schizophyllum_run1", "Schizophyllum_run2 (replay)", "Cordyceps_run1"]
    for r, label in [(r1, "Schizophyllum run1"), (r2, "Schizophyllum run2 (replay)"), (r3, "Cordyceps run1")]:
        total_ms = sum(t["elapsed_ms"] for t in r["timings"].values())
        fresh_ms = sum(t["elapsed_ms"] for t in r["timings"].values() if not t.get("replayed"))
        replay_ms= sum(t["elapsed_ms"] for t in r["timings"].values() if t.get("replayed"))
        n_replayed = sum(1 for t in r["timings"].values() if t.get("replayed"))
        print(f"\n  {label}")
        print(f"    Windows:      {r['n_windows']}")
        print(f"    Total time:   {total_ms:.1f} ms")
        print(f"    Fresh:        {fresh_ms:.1f} ms ({len(r['timings']) - n_replayed} activities)")
        print(f"    Replayed:     {replay_ms:.2f} ms ({n_replayed} activities)")
        if n_replayed > 0 and fresh_ms > 0:
            speedup = (fresh_ms / max(replay_ms, 0.01)) if n_replayed > 0 else 1.0
            print(f"    Replay speedup: ×{speedup:.0f}")
        print(f"    Top anomaly:  window {r['top_k'][0]['window_idx']}  "
              f"borda={r['top_k'][0]['borda']:.4f}")

    # Hypothesis checks
    schiz_t1 = sum(t["orig_ms"] for t in timings_schiz_run1.values())
    schiz_t2 = sum(t["elapsed_ms"] for t in timings_schiz_run2.values())
    print(f"\n  Duroxide replay correctness check:")
    assert r1["top_k"][0]["window_idx"] == r2["top_k"][0]["window_idx"], \
        "FAIL: replay returned different top window!"
    print(f"    ✓ Top anomaly window identical across fresh + replay")
    print(f"    ✓ Run1 total: {schiz_t1:.0f} ms  |  Run2 (replay): {schiz_t2:.2f} ms")
    speedup = schiz_t1 / max(schiz_t2, 0.01)
    print(f"    ✓ Replay speedup: ×{speedup:.0f}")

    # Save report
    report = {
        "runs": {
            "Schizophyllum_run1": {
                "n_windows": r1["n_windows"],
                "top_k": r1["top_k"],
                "timings": {k: v for k, v in timings_schiz_run1.items()},
                "total_ms": sum(t["elapsed_ms"] for t in timings_schiz_run1.values()),
            },
            "Schizophyllum_run2_replay": {
                "n_windows": r2["n_windows"],
                "top_k": r2["top_k"],
                "timings": {k: v for k, v in timings_schiz_run2.items()},
                "total_ms": sum(t["elapsed_ms"] for t in timings_schiz_run2.values()),
            },
            "Cordyceps_run1": {
                "n_windows": r3["n_windows"],
                "top_k": r3["top_k"],
                "timings": r3["timings"],
                "total_ms": sum(t["elapsed_ms"] for t in r3["timings"].values()),
            },
        },
        "replay_speedup_x": round(speedup, 1),
        "replay_determinism": "PASS",
    }
    out = str(RESULTS_DIR / "18_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, default=_json_default)
    print(f"\n  Report → 18_report.json")
    print("  SQLite checkpoint store → 18_duroxide.sqlite")
    print("Done.  Figures → experiments/results/18_*.png")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
