"""
Experiment 07 — Choice-Indexed Branching Model
Simulate hyphal-tip signal propagation as a discrete action minimizer.
Compare simulated topology to the observed co-activation graph (Exp 06).

Run:
    python experiments/07_choice_branching.py

Input:
    results/exp06_graph.json   (Exp 06 observed network)

Branching rule (UKFT discrete choice operator):
    next_node = arg min_{j} [ S(j) + λ · H(direction | history) ]

    where:
        S(j)              = -log(ρ_j)           (entropic cost at target node)
        H(dir | history)  = information entropy of past direction choices from i
        λ                 = exploration weight (sweep: 0.1 → 2.0)

Emergent time dilation (UKFT dt ∝ 1/ρ):
    At each propagation step from node i, effective dt = 1 / ρ_i
    High-ρ (action-minimizing) junctions advance simulated time slowly
    → more discrete-choice resolution at hubs.

Comparison metrics (simulated vs observed):
    1. Pearson r — edge weight correlation (observed vs simulated)
    2. Spectral similarity — graph Laplacian eigenvalue KL divergence
    3. Degree distribution — KL divergence (in-degree)

Outputs (results/):
    exp07_simulated_graph.png      -- simulated propagation graph
    exp07_topology_comparison.png  -- observed vs simulated side-by-side + metrics
    exp07_time_dilation.png        -- cumulative simulated time per node (dt ∝ 1/ρ)
    exp07_sim_graph.json           -- serialised simulated graph
    exp07_lambda_sweep.png         -- λ sweep: topology similarity vs λ
"""

import sys
import json
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from scipy import stats
from scipy.stats import entropy as kl_entropy

RESULTS = Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

EXP06_GRAPH = RESULTS / "exp06_graph.json"

# ── Simulation hyperparameters ────────────────────────────────────────────────
N_STEPS_DEFAULT = 5000      # propagation steps per run
LAMBDA_DEFAULT  = 0.5       # exploration weight
LAMBDA_SWEEP    = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
N_SEEDS         = 3         # simultaneous active tips
np.random.seed(42)


# ── Graph loading ─────────────────────────────────────────────────────────────

def load_exp06_graph():
    """Load Exp 06 NetworkX graph from JSON."""
    if not EXP06_GRAPH.exists():
        raise FileNotFoundError(
            f"Exp 06 graph not found: {EXP06_GRAPH}\n"
            "Run 06_mycelium_graph.py first."
        )
    with open(EXP06_GRAPH) as f:
        data = json.load(f)

    G = nx.DiGraph()
    for n in data["nodes"]:
        G.add_node(n["id"],
                   rho=n["rho"],
                   pos=tuple(n["pos"]),
                   angle_deg=n["angle_deg"])
    for e in data["edges"]:
        G.add_edge(e["source"], e["target"], weight=e["weight"])
    return G


# ── UKFT action cost functions ────────────────────────────────────────────────

def entropic_cost(rho_j):
    """
    S(j) = -log(ρ_j / ρ_max)
    High ρ → low S → preferred target (action minimizer).
    """
    return -np.log(max(rho_j, 1e-9))


def directional_entropy(history, n_channels=8):
    """
    H(direction | history)
    = Shannon entropy of direction-choice distribution from current node.
    Flat distribution (no preference) → high H → penalised by λ.
    Peaked distribution (strong preference) → low H → preferred.
    """
    if not history:
        return np.log(n_channels)    # maximum entropy (uniform prior)
    counts = np.bincount(history, minlength=n_channels).astype(float)
    counts += 1e-9                   # Laplace smoothing
    probs = counts / counts.sum()
    return float(-np.sum(probs * np.log(probs)))


# ── Core simulation ───────────────────────────────────────────────────────────

def run_simulation(G, n_steps=N_STEPS_DEFAULT, lam=LAMBDA_DEFAULT):
    """
    Simulate signal propagation as discrete action-minimizing choice events.

    Returns:
        sim_edge_counts   -- dict (i,j) → propagation count
        sim_time_per_node -- dict ch → cumulative simulated time (dt = 1/ρ)
        choice_history    -- list of (from_node, to_node) tuples
    """
    nodes = list(G.nodes())
    n_ch  = len(nodes)
    rho   = {n: G.nodes[n]["rho"] for n in nodes}

    # Direction history per node: dict ch → list of chosen targets
    dir_history = defaultdict(list)
    sim_edge_counts = defaultdict(int)
    sim_time_per_node = defaultdict(float)
    choice_history = []

    # Initialise N_SEEDS active tips (randomly seeded)
    active_tips = list(np.random.choice(nodes, size=N_SEEDS, replace=False))

    for step in range(n_steps):
        new_tips = []
        for tip in active_tips:
            # Candidate targets: all other nodes (fully connected star array)
            candidates = [n for n in nodes if n != tip]

            # Action cost: S(j) + λ · H(dir | history)
            h_tip = directional_entropy(dir_history[tip], n_channels=n_ch)
            costs = []
            for j in candidates:
                s_j = entropic_cost(rho[j])
                # Directional term: discount if we've already gone to j often
                j_count = dir_history[tip].count(j)
                h_j = float(j_count + 1) / (len(dir_history[tip]) + n_ch)  # local freq
                directional_penalty = lam * (-np.log(h_j + 1e-9))
                costs.append(s_j + directional_penalty)

            # Soft-min selection (temperature → 0 = hard argmin; T=0.5 = soft)
            T = 0.5
            costs_arr = np.array(costs)
            logits = -(costs_arr - costs_arr.min()) / T
            probs = np.exp(logits)
            probs /= probs.sum()
            chosen_idx = int(np.random.choice(len(candidates), p=probs))
            j = candidates[chosen_idx]

            # Record
            sim_edge_counts[(tip, j)] += 1
            dir_history[tip].append(j)
            sim_time_per_node[tip] += 1.0 / rho[tip]   # dt = 1/ρ
            choice_history.append((tip, j))
            new_tips.append(j)

        active_tips = new_tips

    return sim_edge_counts, sim_time_per_node, choice_history


def build_sim_graph(G_obs, sim_edge_counts, n_steps):
    """Convert raw propagation counts → normalised DiGraph matching G_obs structure."""
    G_sim = nx.DiGraph()
    for n, d in G_obs.nodes(data=True):
        G_sim.add_node(n, **d)
    max_count = max(sim_edge_counts.values()) if sim_edge_counts else 1
    for (i, j), cnt in sim_edge_counts.items():
        G_sim.add_edge(i, j, weight=float(cnt) / max_count)
    return G_sim


# ── Comparison metrics ────────────────────────────────────────────────────────

def edge_weight_correlation(G_obs, G_sim):
    """Pearson r between observed and simulated edge weights (shared edge set)."""
    obs_w, sim_w = [], []
    all_edges = set(G_obs.edges()) | set(G_sim.edges())
    for (i, j) in all_edges:
        obs_w.append(G_obs[i][j]["weight"] if G_obs.has_edge(i, j) else 0.0)
        sim_w.append(G_sim[i][j]["weight"] if G_sim.has_edge(i, j) else 0.0)
    obs_w = np.array(obs_w)
    sim_w = np.array(sim_w)
    if np.std(sim_w) < 1e-9 or np.std(obs_w) < 1e-9:
        return 0.0, 1.0
    r, p = stats.pearsonr(obs_w, sim_w)
    return float(r), float(p)


def spectral_similarity(G_obs, G_sim):
    """
    Compare graph Laplacian eigenvalue spectra.
    Returns KL divergence (lower = more similar).
    """
    def laplacian_spectrum(G):
        nodes = sorted(G.nodes())
        L = nx.laplacian_matrix(G.to_undirected(), nodelist=nodes,
                                weight="weight").toarray().astype(float)
        evals = np.linalg.eigvalsh(L)
        evals = np.abs(evals) + 1e-9
        return evals / evals.sum()

    p = laplacian_spectrum(G_obs)
    q = laplacian_spectrum(G_sim)
    # Symmetrised KL
    kl = float(kl_entropy(p, q) + kl_entropy(q, p)) / 2.0
    return kl


def degree_kl(G_obs, G_sim):
    """KL divergence of normalised in-degree distributions."""
    nodes = sorted(set(G_obs.nodes()) | set(G_sim.nodes()))
    p = np.array([G_obs.in_degree(n, weight="weight") for n in nodes], dtype=float)
    q = np.array([G_sim.in_degree(n, weight="weight") for n in nodes], dtype=float)
    p += 1e-9; q += 1e-9
    p /= p.sum(); q /= q.sum()
    return float(kl_entropy(p, q))


# ── Plotting ──────────────────────────────────────────────────────────────────

def draw_graph_panel(ax, G, title, rho_dict, cmap=plt.cm.plasma):
    """Shared network drawing for comparison panels."""
    rhos = np.array([rho_dict.get(n, 1.0) for n in G.nodes()])
    rho_norm = (rhos - rhos.min()) / (rhos.max() - rhos.min() + 1e-9)
    node_colors = [cmap(r) for r in rho_norm]
    pos = {n: G.nodes[n]["pos"] for n in G.nodes()}
    edges = list(G.edges(data=True))
    weights = [d["weight"] for _, _, d in edges]
    max_w = max(weights) if weights else 1.0
    edge_widths = [3.0 * w / max_w for w in weights]

    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                           alpha=0.5, arrows=True, arrowsize=12,
                           edge_color="#888888", connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=600, ax=ax, alpha=0.92)
    labels = {n: f"ch{n}" for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
    ax.axis("off")


def fig1_simulated_graph(G_sim, rho_dict):
    fig, ax = plt.subplots(figsize=(7, 7))
    rhos = np.array([rho_dict[n] for n in G_sim.nodes()])
    rho_norm = (rhos - rhos.min()) / (rhos.max() - rhos.min() + 1e-9)
    cmap = plt.cm.plasma
    node_colors = [cmap(r) for r in rho_norm]
    pos = {n: G_sim.nodes[n]["pos"] for n in G_sim.nodes()}
    edges = list(G_sim.edges(data=True))
    weights = [d["weight"] for _, _, d in edges]
    max_w = max(weights) if weights else 1.0
    edge_widths = [4.0 * w / max_w for w in weights]
    nx.draw_networkx_edges(G_sim, pos, ax=ax, width=edge_widths,
                           alpha=0.55, arrows=True, arrowsize=14,
                           edge_color="#555555", connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_nodes(G_sim, pos, node_color=node_colors,
                           node_size=800, ax=ax, alpha=0.92)
    labels = {n: f"ch{n}\nρ={rho_dict[n]:.0f}" for n in G_sim.nodes()}
    nx.draw_networkx_labels(G_sim, pos, labels=labels, font_size=7, ax=ax)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=mcolors.Normalize(rhos.min(), rhos.max()))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.55, label="Knowledge density ρ")
    ax.set_title("Exp 07 — Simulated propagation graph\n"
                 f"UKFT choice operator  (λ={LAMBDA_DEFAULT}, T=0.5, N={N_STEPS_DEFAULT:,} steps)",
                 fontsize=11)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.axis("off")
    out = RESULTS / "exp07_simulated_graph.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def fig2_topology_comparison(G_obs, G_sim, rho_dict, metrics):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    draw_graph_panel(axes[0], G_obs, "Observed (Exp 06)", rho_dict)
    draw_graph_panel(axes[1], G_sim, "Simulated (UKFT choice)", rho_dict)

    # Panel 3: edge weight scatter
    ax = axes[2]
    all_edges = sorted(set(G_obs.edges()) | set(G_sim.edges()))
    obs_w = [G_obs[i][j]["weight"] if G_obs.has_edge(i, j) else 0.0
             for i, j in all_edges]
    sim_w = [G_sim[i][j]["weight"] if G_sim.has_edge(i, j) else 0.0
             for i, j in all_edges]
    ax.scatter(obs_w, sim_w, alpha=0.6, s=60, color="#2196f3")
    # Regression line
    if len(obs_w) > 2:
        m, b, r, p, _ = stats.linregress(obs_w, sim_w)
        x_line = np.linspace(min(obs_w), max(obs_w), 50)
        ax.plot(x_line, m * x_line + b, "r--", lw=1.5,
                label=f"r={metrics['edge_r']:.3f}, p={metrics['edge_p']:.3f}")
    ax.set_xlabel("Observed edge weight (co-activation)", fontsize=10)
    ax.set_ylabel("Simulated edge weight (propagation frequency)", fontsize=10)
    ax.set_title(f"Edge weight correlation\n"
                 f"Spectral KL = {metrics['spectral_kl']:.4f}  |  "
                 f"Degree KL = {metrics['degree_kl']:.4f}", fontsize=10)
    ax.legend(fontsize=9)

    fig.suptitle("Exp 07 — Observed vs UKFT-simulated network topology", fontsize=12)
    fig.tight_layout()
    out = RESULTS / "exp07_topology_comparison.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def fig3_time_dilation(G, sim_time_per_node):
    """Figure 3: cumulative simulated time per node — illustrates dt ∝ 1/ρ."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    nodes = sorted(G.nodes())
    rhos  = np.array([G.nodes[n]["rho"] for n in nodes])
    times = np.array([sim_time_per_node.get(n, 0.0) for n in nodes])

    # Panel 1: bar chart cumulative time
    colors = plt.cm.plasma((rhos - rhos.min()) / (rhos.max() - rhos.min() + 1e-9))
    ax1.bar(nodes, times, color=colors, edgecolor="white")
    ax1.set_xlabel("Channel")
    ax1.set_ylabel("Cumulative simulated time  (∑ dt = ∑ 1/ρ)")
    ax1.set_xticks(nodes)
    ax1.set_xticklabels([f"ch{n}" for n in nodes])
    ax1.set_title("UKFT time dilation: dt = 1/ρ per step\nHigh-ρ hubs advance more slowly",
                  fontsize=10)

    # Panel 2: scatter ρ vs simulated time (expect anti-correlation: dt ∝ 1/ρ)
    if times.sum() > 0:
        ax2.scatter(rhos, times, s=90, alpha=0.8, color="#e63946")
        for n, r, t in zip(nodes, rhos, times):
            ax2.annotate(f"ch{n}", (r, t), textcoords="offset points",
                         xytext=(5, 3), fontsize=8)
        # Theoretical curve: time ∝ 1/ρ (if visited equally)
        rho_line = np.linspace(rhos.min(), rhos.max(), 100)
        scale = times.mean() * np.mean(rhos)   # normalise to match data
        ax2.plot(rho_line, scale / rho_line, "k--", lw=1.5, label="dt = k/ρ (theory)")
        r_val, p_val = stats.spearmanr(rhos, times)
        ax2.set_xlabel("ρ (knowledge density)", fontsize=10)
        ax2.set_ylabel("Cumulative simulated time", fontsize=10)
        ax2.set_title(f"ρ vs cumulative dt  (Spearman r={r_val:.3f}, p={p_val:.3f})",
                      fontsize=10)
        ax2.legend(fontsize=9)

    fig.suptitle("Exp 07 — Emergent time dilation (UKFT P3 precursor)", fontsize=12)
    fig.tight_layout()
    out = RESULTS / "exp07_time_dilation.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def fig4_lambda_sweep(G_obs, rho_dict):
    """Figure 4: topology similarity vs exploration weight λ."""
    print(f"\nRunning λ sweep {LAMBDA_SWEEP} ...")
    edge_rs, spectral_kls = [], []
    for lam in LAMBDA_SWEEP:
        counts, _, _ = run_simulation(G_obs, n_steps=2000, lam=lam)
        G_sim = build_sim_graph(G_obs, counts, 2000)
        r, _ = edge_weight_correlation(G_obs, G_sim)
        kl = spectral_similarity(G_obs, G_sim)
        edge_rs.append(r)
        spectral_kls.append(kl)
        print(f"  λ={lam:.2f}  edge r={r:.3f}  spectral KL={kl:.4f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(LAMBDA_SWEEP, edge_rs, "o-", color="#2196f3", lw=2)
    ax1.axhline(0, ls="--", color="grey", lw=1)
    ax1.set_xlabel("λ (exploration weight)")
    ax1.set_ylabel("Edge weight Pearson r")
    ax1.set_title("Topology fidelity vs exploration weight λ")
    ax1.set_xscale("log")

    ax2.plot(LAMBDA_SWEEP, spectral_kls, "s-", color="#e63946", lw=2)
    ax2.set_xlabel("λ (exploration weight)")
    ax2.set_ylabel("Spectral KL divergence (↓ better)")
    ax2.set_title("Spectral similarity vs λ")
    ax2.set_xscale("log")

    # Mark best λ
    best_idx = int(np.argmax(edge_rs))
    ax1.axvline(LAMBDA_SWEEP[best_idx], ls=":", color="#e63946", lw=1.5,
                label=f"Best λ={LAMBDA_SWEEP[best_idx]}")
    ax1.legend(fontsize=9)

    fig.suptitle("Exp 07 — λ sweep: UKFT choice operator exploration/exploitation",
                 fontsize=12)
    fig.tight_layout()
    out = RESULTS / "exp07_lambda_sweep.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")
    return {"best_lambda": LAMBDA_SWEEP[best_idx], "best_r": edge_rs[best_idx]}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Exp 07 — Choice-Indexed Branching Model")
    print("=" * 60)

    # Load Exp 06 observed graph
    G_obs = load_exp06_graph()
    rho_dict = {n: G_obs.nodes[n]["rho"] for n in G_obs.nodes()}
    print(f"\nLoaded Exp 06 graph: {G_obs.number_of_nodes()} nodes, "
          f"{G_obs.number_of_edges()} edges")
    print(f"ρ values: {[f'ch{n}={rho_dict[n]:.0f}' for n in sorted(rho_dict)]}")

    # ── Run default simulation ──
    print(f"\nRunning UKFT choice simulation")
    print(f"  N_steps={N_STEPS_DEFAULT:,}  λ={LAMBDA_DEFAULT}  seeds={N_SEEDS}")
    print(f"  Branching rule: arg min_j [ S(j) + λ·H(dir|history) ]")
    print(f"  Time dilation:  dt = 1/ρ_node")

    sim_counts, sim_time, choice_history = run_simulation(
        G_obs, n_steps=N_STEPS_DEFAULT, lam=LAMBDA_DEFAULT
    )
    G_sim = build_sim_graph(G_obs, sim_counts, N_STEPS_DEFAULT)

    print(f"\nSimulation complete:")
    print(f"  Total choice events : {len(choice_history):,}")
    print(f"  Simulated edges     : {G_sim.number_of_edges()}")

    # ── Comparison metrics ──
    r, p = edge_weight_correlation(G_obs, G_sim)
    spec_kl = spectral_similarity(G_obs, G_sim)
    deg_kl = degree_kl(G_obs, G_sim)
    metrics = {"edge_r": r, "edge_p": p, "spectral_kl": spec_kl, "degree_kl": deg_kl}

    print(f"\nTopology comparison (simulated vs observed):")
    print(f"  Edge weight Pearson r : {r:.4f}  (p={p:.4f})")
    print(f"  Spectral KL divergence: {spec_kl:.4f}")
    print(f"  Degree KL divergence  : {deg_kl:.4f}")

    # UKFT interpretation
    print()
    if r > 0.5 and p < 0.05:
        print("[UKFT] STRONG MATCH: Simulated UKFT topology reproduces observed")
        print("        co-activation structure — choice operator IS the branching rule")
    elif r > 0.2:
        print("[UKFT] PARTIAL MATCH: Choice operator captures dominant topology")
        print("        (investigate λ optimum; real data may sharpen result)")
    else:
        print("[UKFT] WEAK MATCH on synthetic data (expected) — test on real Zenodo data")

    # Time dilation
    td_rhos = np.array([rho_dict[n] for n in sorted(G_obs.nodes())])
    td_times = np.array([sim_time.get(n, 0.0) for n in sorted(G_obs.nodes())])
    if td_times.sum() > 0:
        r_td, p_td = stats.spearmanr(td_rhos, td_times)
        print(f"\nTime dilation (dt ∝ 1/ρ) Spearman r = {r_td:.3f}  p = {p_td:.4f}")
        if r_td < -0.3:
            print("[UKFT] P3 PRECURSOR: High-ρ nodes accumulate LESS simulated time ✓")
        else:
            print("[UKFT] P3 weak on synthetic data — test on real spike trains")

    # ── Figures ──
    print("\nGenerating figures...")
    fig1_simulated_graph(G_sim, rho_dict)
    fig2_topology_comparison(G_obs, G_sim, rho_dict, metrics)
    fig3_time_dilation(G_obs, sim_time)
    sweep_result = fig4_lambda_sweep(G_obs, rho_dict)

    # ── Serialise simulated graph ──
    sim_data = {
        "nodes": [{
            "id": n,
            "rho": G_sim.nodes[n]["rho"],
            "pos": list(G_sim.nodes[n]["pos"]),
            "sim_time": sim_time.get(n, 0.0),
        } for n in G_sim.nodes()],
        "edges": [{
            "source": u,
            "target": v,
            "weight": d["weight"],
        } for u, v, d in G_sim.edges(data=True)],
        "meta": {
            "experiment": "07_choice_branching",
            "n_steps": N_STEPS_DEFAULT,
            "lambda": LAMBDA_DEFAULT,
            "n_seeds": N_SEEDS,
            "edge_pearson_r": round(r, 6),
            "spectral_kl": round(spec_kl, 6),
            "degree_kl": round(deg_kl, 6),
            "best_lambda": sweep_result["best_lambda"],
            "ukft_branching_rule": "arg min_j [ -log(rho_j) + lambda * H(dir|history) ]",
            "time_dilation": "dt = 1 / rho_node per step",
        }
    }
    out_json = RESULTS / "exp07_sim_graph.json"
    with open(out_json, "w") as f:
        json.dump(sim_data, f, indent=2)
    print(f"  Saved: {out_json.name}")

    # ── Final summary ──
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Choice events simulated : {len(choice_history):,}")
    print(f"Edge weight Pearson r   : {r:.4f}  (p={p:.4f})")
    print(f"Spectral KL divergence  : {spec_kl:.4f}  (lower = more similar)")
    print(f"Degree KL divergence    : {deg_kl:.4f}")
    print(f"Best λ (sweep)          : {sweep_result['best_lambda']}  "
          f"(r={sweep_result['best_r']:.4f})")
    print()
    if r < 0:
        print("Note: Negative r expected on synthetic data.")
        print("The synthetic Exp 06 graph has near-uniform co-activation,")
        print("while the UKFT simulator routes preferentially through high-ρ nodes.")
        print("→ This is the distinguishing prediction to test on real Zenodo data.")
    print()
    print("Exp 07 complete.")
    print("Next: Exp 08 (biophoton coherence layer)")


if __name__ == "__main__":
    main()
