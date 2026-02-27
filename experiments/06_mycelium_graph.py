"""
Experiment 06 -- Mycelium Graph Construction
Build weighted hyphal network from 8-channel spike-train data.
Compute per-node knowledge density ρ and co-activation adjacency matrix.

Run:
    python experiments/06_mycelium_graph.py

Data:
    Real:       noosphere/apps/myco-explorer/tools/data/pleurotus_spike_events.ndjson
                (run export_spike_events.py first, requires Adamatzky 2026 supplementary)
    Synthetic:  auto-generated if real data absent (8-channel star array, ~Adamatzky 2026)

Outputs (results/):
    exp06_network_rho.png         -- graph coloured by knowledge density ρ
    exp06_coactivation_heatmap.png -- channel co-activation matrix
    exp06_rho_distribution.png    -- per-channel ρ histogram
    exp06_graph.json              -- serialised NetworkX graph (for Exp 07+)
    exp06_node_stats.csv          -- per-channel stats table
"""

import sys
import json
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from scipy import stats

RESULTS = Path(__file__).resolve().parent.parent / "results"
RESULTS.mkdir(exist_ok=True)

# Real data path (noosphere myco-explorer pipeline output)
NOOSPHERE = Path(__file__).resolve().parent.parent.parent / "noosphere"
REAL_DATA = NOOSPHERE / "apps" / "myco-explorer" / "tools" / "data" / "pleurotus_spike_events.ndjson"

# ── Geometry ──────────────────────────────────────────────────────────────────
# 8-channel star-shaped differential electrode array (Adamatzky 2026 arXiv:2601.08099)
N_CHANNELS = 8
CHANNEL_ANGLES_DEG = np.arange(0, 360, 45)   # 0,45,90,...,315
CHANNEL_XY = {
    ch: (np.cos(np.radians(a)), np.sin(np.radians(a)))
    for ch, a in enumerate(CHANNEL_ANGLES_DEG)
}

# ── Analysis parameters ───────────────────────────────────────────────────────
BIN_S             = 60.0   # 1-minute spike-rate bins for ρ computation
COACT_WINDOW_S    = 30.0   # co-activation window (same burst → adjacent)
MIN_COACT_WEIGHT  = 0.01   # minimum edge weight to include in graph
np.random.seed(42)


# ── Synthetic data generator (Adamatzky 2026 profile) ────────────────────────

def generate_synthetic_spike_trains(duration_s=7200.0):
    """
    8-channel synthetic spike train mimicking Adamatzky 2026 star-array.

    Biological parameters calibrated to arXiv:2601.08099:
    - Spike duration: 2–21 hours (observed range); here we use a 2-hour window
      to get sufficient statistics. Scale as needed.
    - Spike amplitude: 0.3–2.1 mV log-normal
    - Background ISI: log-normal centred ~280 s (Adamatzky 2021 arXiv:2112.09907)
    - Burst structure: 3-8 spikes, propagation delay ~180 s / 2 cm = 90 s/cm

    Returns list of dicts matching export_spike_events.py NDJSON format.
    """
    events = []

    # Channel baseline rates (spikes/hour) — non-uniform to produce interesting ρ
    base_rates = {
        0: 4.5,   # hub-like (high ρ)
        1: 2.1,
        2: 5.8,   # hub
        3: 1.3,
        4: 4.0,
        5: 1.8,
        6: 6.2,   # hub
        7: 2.7,
    }

    # Burst events: correlated across channel pairs
    burst_schedule = []
    t = 120.0
    while t < duration_s - 300:
        root_ch = np.random.choice(N_CHANNELS)
        burst_size = np.random.randint(2, 6)
        # Propagation: fires nearby channels with 60-180s delay
        partners = []
        for _ in range(burst_size - 1):
            other = (root_ch + np.random.randint(1, 4)) % N_CHANNELS
            delay = np.random.exponential(90.0)
            partners.append((other, delay))
        burst_schedule.append((t, root_ch, partners))
        t += np.random.exponential(600.0)

    event_id = 0
    for t_root, root_ch, partners in burst_schedule:
        # Root spike
        amp = np.random.lognormal(np.log(0.8), 0.4)
        amp = float(np.clip(amp, 0.03, 2.1))
        dur = float(np.random.exponential(60.0) + 5.0)
        events.append({
            "event_id": f"syn_ch{root_ch}_{event_id:05d}",
            "channel": root_ch,
            "t_start_s": t_root,
            "t_peak_s": t_root + dur * 0.3,
            "t_end_s": t_root + dur,
            "amplitude_mv": amp,
            "duration_s": dur,
            "in_burst": True,
            "burst_id": f"b{event_id:04d}",
            "preceding_isi_s": None,
        })
        # Partner spikes
        for p_ch, delay in partners:
            p_amp = amp * np.random.uniform(0.6, 1.0)
            p_dur = float(np.random.exponential(40.0) + 5.0)
            events.append({
                "event_id": f"syn_ch{p_ch}_{event_id:05d}p",
                "channel": p_ch,
                "t_start_s": t_root + delay,
                "t_peak_s": t_root + delay + p_dur * 0.3,
                "t_end_s": t_root + delay + p_dur,
                "amplitude_mv": float(p_amp),
                "duration_s": p_dur,
                "in_burst": True,
                "burst_id": f"b{event_id:04d}",
                "preceding_isi_s": delay,
            })
        event_id += 1

    # Background isolated spikes
    for ch in range(N_CHANNELS):
        rate = base_rates[ch] / 3600.0   # per second
        t = np.random.exponential(1.0 / rate)
        while t < duration_s:
            amp = float(np.clip(np.random.lognormal(np.log(0.3), 0.5), 0.03, 1.0))
            dur = float(np.random.exponential(20.0) + 2.0)
            events.append({
                "event_id": f"syn_ch{ch}_{event_id:05d}bg",
                "channel": ch,
                "t_start_s": t,
                "t_peak_s": t + dur * 0.3,
                "t_end_s": t + dur,
                "amplitude_mv": amp,
                "duration_s": dur,
                "in_burst": False,
                "burst_id": None,
                "preceding_isi_s": None,
            })
            event_id += 1
            t += np.random.exponential(1.0 / rate)

    events.sort(key=lambda e: e["t_start_s"])
    return events


# ── Core analysis ─────────────────────────────────────────────────────────────

def compute_spike_rate_series(events, ch, duration_s):
    """Bin spike counts into BIN_S windows → spike-rate time series per channel."""
    n_bins = max(1, int(duration_s / BIN_S))
    counts = np.zeros(n_bins)
    for e in events:
        if e["channel"] == ch:
            idx = min(int(e["t_start_s"] / BIN_S), n_bins - 1)
            counts[idx] += 1
    return counts / BIN_S   # spikes/second


def compute_rho(rate_series):
    """
    Knowledge density ρ = 1 / σ²  (σ = spike-rate variance).
    High ρ → low variance → predictable, high-information node.
    UKFT: high-ρ nodes are action-minimizing junctions (Exp 06 UKFT prediction).
    """
    sigma2 = float(np.var(rate_series))
    return 1.0 / sigma2 if sigma2 > 1e-9 else 1e6


def build_coactivation_matrix(events, duration_s):
    """
    Build N_CHANNELS × N_CHANNELS weighted co-activation matrix.
    Entry [i,j] = fraction of spikes on channel i that have a co-occurring
    spike on channel j within COACT_WINDOW_S.
    Normalised row-wise → [0, 1] edge weights.
    """
    mat = np.zeros((N_CHANNELS, N_CHANNELS))
    # Index events by channel
    ch_events = {ch: [] for ch in range(N_CHANNELS)}
    for e in events:
        ch_events[e["channel"]].append(e["t_peak_s"])

    for i in range(N_CHANNELS):
        if not ch_events[i]:
            continue
        for j in range(N_CHANNELS):
            if i == j or not ch_events[j]:
                continue
            t_i = np.array(ch_events[i])
            t_j = np.array(ch_events[j])
            n_coact = 0
            for t in t_i:
                diffs = np.abs(t_j - t)
                if np.any(diffs <= COACT_WINDOW_S):
                    n_coact += 1
            mat[i, j] = n_coact / len(t_i)

    return mat


def build_graph(coact_mat, rho_dict):
    """Build NetworkX directed weighted graph from co-activation matrix."""
    G = nx.DiGraph()
    for ch in range(N_CHANNELS):
        x, y = CHANNEL_XY[ch]
        G.add_node(ch,
                   rho=rho_dict[ch],
                   pos=(x, y),
                   angle_deg=int(CHANNEL_ANGLES_DEG[ch]))
    for i in range(N_CHANNELS):
        for j in range(N_CHANNELS):
            w = coact_mat[i, j]
            if w >= MIN_COACT_WEIGHT:
                G.add_edge(i, j, weight=float(w))
    return G


# ── Plotting ──────────────────────────────────────────────────────────────────

def fig1_network_rho(G, node_stats):
    """Figure 1: Network coloured by ρ with edge width ∝ co-activation weight."""
    fig, ax = plt.subplots(figsize=(8, 8))

    rhos = np.array([G.nodes[n]["rho"] for n in G.nodes()])
    rho_norm = (rhos - rhos.min()) / (rhos.max() - rhos.min() + 1e-9)
    cmap = plt.cm.plasma
    node_colors = [cmap(r) for r in rho_norm]

    pos = {n: G.nodes[n]["pos"] for n in G.nodes()}
    edges = list(G.edges(data=True))
    edge_weights = [d["weight"] for _, _, d in edges]
    max_w = max(edge_weights) if edge_weights else 1.0
    edge_widths = [3.0 * w / max_w for w in edge_weights]
    edge_alphas = [0.3 + 0.7 * w / max_w for w in edge_weights]

    # Draw edges
    for (u, v, d), lw, alpha in zip(edges, edge_widths, edge_alphas):
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                               width=lw, alpha=alpha, ax=ax,
                               arrows=True, arrowsize=15,
                               edge_color="#aaaaaa", connectionstyle="arc3,rad=0.1")

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=800, ax=ax, alpha=0.92)
    labels = {n: f"ch{n}\nρ={G.nodes[n]['rho']:.1f}" for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, ax=ax)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=mcolors.Normalize(rhos.min(), rhos.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label("Knowledge density ρ = 1/σ²", fontsize=10)

    ax.set_title("Exp 06 — Pleurotus ostreatus 8-channel star array\n"
                 "Node colour = UKFT knowledge density ρ  |  Edge = co-activation strength",
                 fontsize=11)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.axis("off")

    out = RESULTS / "exp06_network_rho.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")
    return out


def fig2_coactivation_heatmap(coact_mat):
    """Figure 2: Co-activation matrix heatmap."""
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(coact_mat, cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=coact_mat.max())
    fig.colorbar(im, ax=ax, label="Co-activation fraction")
    ax.set_xticks(range(N_CHANNELS))
    ax.set_yticks(range(N_CHANNELS))
    ax.set_xticklabels([f"ch{i}" for i in range(N_CHANNELS)])
    ax.set_yticklabels([f"ch{i}" for i in range(N_CHANNELS)])
    ax.set_xlabel("Target channel")
    ax.set_ylabel("Source channel")
    ax.set_title("Exp 06 — Inter-channel co-activation matrix\n"
                 f"(window ±{COACT_WINDOW_S:.0f} s)", fontsize=11)
    # Annotate with values
    for i in range(N_CHANNELS):
        for j in range(N_CHANNELS):
            ax.text(j, i, f"{coact_mat[i,j]:.2f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if coact_mat[i, j] > coact_mat.max() * 0.6 else "black")
    out = RESULTS / "exp06_coactivation_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")
    return out


def fig3_rho_distribution(node_stats):
    """Figure 3: Per-channel ρ bar chart with UKFT hub threshold."""
    channels = [s["channel"] for s in node_stats]
    rhos = np.array([s["rho"] for s in node_stats])
    spike_counts = [s["n_spikes"] for s in node_stats]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # Panel 1: ρ per channel
    rho_median = float(np.median(rhos))
    colors = ["#e63946" if r > rho_median * 1.5 else "#457b9d" for r in rhos]
    bars = ax1.bar(channels, rhos, color=colors, edgecolor="white")
    ax1.axhline(rho_median * 1.5, ls="--", color="#e63946", lw=1.5,
                label=f"Hub threshold = 1.5×median = {rho_median*1.5:.1f}")
    ax1.set_ylabel("ρ = 1/σ²  (knowledge density)", fontsize=10)
    ax1.set_title("Exp 06 — Per-channel knowledge density ρ", fontsize=11)
    ax1.legend(fontsize=9)

    # Annotate hub channels
    for bar, r, ch in zip(bars, rhos, channels):
        if r > rho_median * 1.5:
            ax1.text(bar.get_x() + bar.get_width() / 2, r + rhos.max() * 0.02,
                     "HUB", ha="center", va="bottom", fontsize=8,
                     color="#e63946", fontweight="bold")

    # Panel 2: spike count per channel
    ax2.bar(channels, spike_counts, color="#a8dadc", edgecolor="white")
    ax2.set_xlabel("Channel")
    ax2.set_ylabel("Spike count")
    ax2.set_title("Spike count per channel", fontsize=10)
    ax2.set_xticks(channels)
    ax2.set_xticklabels([f"ch{c}" for c in channels])

    fig.tight_layout()
    out = RESULTS / "exp06_rho_distribution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")
    return out


# ── Serialise graph ───────────────────────────────────────────────────────────

def save_graph_json(G, out_path):
    """Save graph as JSON for Exp 07+ consumption."""
    data = {
        "nodes": [{
            "id": n,
            "rho": G.nodes[n]["rho"],
            "pos": list(G.nodes[n]["pos"]),
            "angle_deg": G.nodes[n]["angle_deg"],
        } for n in G.nodes()],
        "edges": [{
            "source": u,
            "target": v,
            "weight": d["weight"],
        } for u, v, d in G.edges(data=True)],
        "meta": {
            "n_channels": N_CHANNELS,
            "coact_window_s": COACT_WINDOW_S,
            "bin_s": BIN_S,
            "experiment": "06_mycelium_graph",
            "ukft_prediction": "Hub nodes (high rho) = action-minimizing junctions (P1)",
        }
    }
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Saved: {Path(out_path).name}")


def save_node_stats_csv(node_stats, out_path):
    """Save per-channel statistics table."""
    fields = ["channel", "angle_deg", "n_spikes", "n_bursts", "spike_rate_mean",
              "spike_rate_std", "rho", "sigma", "is_hub"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(node_stats)
    print(f"  Saved: {Path(out_path).name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Exp 06 — Mycelium Graph Construction")
    print("=" * 60)

    # ── Load data ──
    if REAL_DATA.exists():
        print(f"\nLoading REAL spike data: {REAL_DATA}")
        events = []
        with open(REAL_DATA) as f:
            for line in f:
                events.append(json.loads(line.strip()))
        data_source = "Adamatzky 2026 (arXiv:2601.08099) — real"
    else:
        print(f"\nReal data not found: {REAL_DATA}")
        print("Generating synthetic 8-channel star-array data")
        print("(calibrated to Adamatzky 2026 spike parameters)")
        events = generate_synthetic_spike_trains(duration_s=7200.0)
        data_source = "Synthetic (Adamatzky 2026 profile, 7200s)"

    t_max = max(e["t_end_s"] for e in events)
    print(f"\nData source  : {data_source}")
    print(f"Total spikes : {len(events):,}")
    print(f"Duration     : {t_max:.0f} s ({t_max/3600:.2f} h)")

    # ── Per-channel stats & ρ ──
    print("\nComputing per-channel knowledge density ρ...")
    node_stats = []
    rho_dict = {}
    rho_median_ref = None

    for ch in range(N_CHANNELS):
        ch_events = [e for e in events if e["channel"] == ch]
        rate_series = compute_spike_rate_series(events, ch, t_max)
        rho = compute_rho(rate_series)
        sigma = float(np.std(rate_series))
        rho_dict[ch] = rho
        node_stats.append({
            "channel": ch,
            "angle_deg": int(CHANNEL_ANGLES_DEG[ch]),
            "n_spikes": len(ch_events),
            "n_bursts": len(set(e["burst_id"] for e in ch_events if e.get("burst_id"))),
            "spike_rate_mean": round(float(np.mean(rate_series)), 6),
            "spike_rate_std": round(sigma, 6),
            "rho": round(rho, 4),
            "sigma": round(sigma, 6),
            "is_hub": False,  # filled below
        })
        print(f"  ch{ch} ({int(CHANNEL_ANGLES_DEG[ch]):3d}°): "
              f"{len(ch_events):4d} spikes  ρ={rho:.2f}  σ={sigma:.4f}")

    # Mark hub channels (ρ > 1.5 × median)
    rhos = np.array([s["rho"] for s in node_stats])
    rho_median = float(np.median(rhos))
    hub_threshold = rho_median * 1.5
    hub_channels = []
    for s in node_stats:
        s["is_hub"] = s["rho"] > hub_threshold
        if s["is_hub"]:
            hub_channels.append(s["channel"])

    print(f"\nρ median: {rho_median:.2f}  |  Hub threshold: {hub_threshold:.2f}")
    print(f"Hub channels (ρ > {hub_threshold:.2f}): {hub_channels}")

    # UKFT interpretation
    print("\n[UKFT] Hub channels = high knowledge density = action-minimizing junctions")
    print("[UKFT] Prediction P1: structured burst propagation routes THROUGH hubs")

    # ── Co-activation matrix ──
    print("\nBuilding co-activation adjacency matrix...")
    coact_mat = build_coactivation_matrix(events, t_max)
    print(f"  Mean co-activation: {coact_mat[coact_mat > 0].mean():.3f}")
    print(f"  Non-zero edges: {(coact_mat >= MIN_COACT_WEIGHT).sum()}")

    # ── NetworkX graph ──
    G = build_graph(coact_mat, rho_dict)
    print(f"\nGraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Graph metrics
    if G.number_of_edges() > 0:
        in_deg = dict(G.in_degree(weight="weight"))
        out_deg = dict(G.out_degree(weight="weight"))
        top_hub = max(in_deg, key=in_deg.get)
        print(f"  Highest in-degree hub: ch{top_hub} "
              f"(in={in_deg[top_hub]:.3f}, ρ={rho_dict[top_hub]:.2f})")

    # ── Figures ──
    print("\nGenerating figures...")
    fig1_network_rho(G, node_stats)
    fig2_coactivation_heatmap(coact_mat)
    fig3_rho_distribution(node_stats)

    # ── Serialise ──
    print("\nSerialising outputs...")
    save_graph_json(G, RESULTS / "exp06_graph.json")
    save_node_stats_csv(node_stats, RESULTS / "exp06_node_stats.csv")

    # ── Summary ──
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Hub channels      : {hub_channels}")
    print(f"ρ range           : {rhos.min():.2f} – {rhos.max():.2f}")
    print(f"ρ coefficient-of-variation: {rhos.std()/rhos.mean():.3f}")
    print(f"Graph edges       : {G.number_of_edges()}")
    print(f"Mean edge weight  : {np.mean([d['weight'] for _,_,d in G.edges(data=True)]):.3f}"
          f" (co-activation fraction)")
    print()

    # Hypothesis test: do hub channels have higher in-degree? (Spearman)
    in_degrees = np.array([G.in_degree(ch, weight="weight") for ch in range(N_CHANNELS)])
    rho_arr = np.array([rho_dict[ch] for ch in range(N_CHANNELS)])
    if np.std(in_degrees) > 0:
        r, p = stats.spearmanr(rho_arr, in_degrees)
        print(f"Spearman ρ vs in-degree: r={r:.3f}, p={p:.4f}")
        if p < 0.05:
            print("[UKFT] P1 SUPPORTED: High-ρ channels ARE the preferential signal hubs")
        elif p < 0.10:
            print("[UKFT] P1 marginal (p<0.10): weak hub structure — more data needed")
        else:
            print("[UKFT] P1 not supported on this dataset — revisit with real Zenodo data")

    print("\nExp 06 complete.")
    print("Next: Exp 07 (choice-indexed branching model)")


if __name__ == "__main__":
    main()
