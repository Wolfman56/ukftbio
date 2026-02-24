"""
ukft_bio.vis — Plotly Visualizations for Bio Trajectories
===========================================================

Bio equivalent of ukftphys/ukft_sim/vis.py.
All plots use Plotly for interactive HTML outputs (matches ukftphys convention).
"""

from __future__ import annotations

from typing import List, Optional
import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("WARNING: plotly not installed. Visualizations unavailable.")


def plot_rho_trajectory(
    trajectory,
    title: str = "ρ_bio Trajectory",
    output_html: Optional[str] = None,
) -> None:
    """Plot ρ_bio over discrete choice steps n."""
    if not PLOTLY_AVAILABLE:
        return
    ns = [r.n for r in trajectory]
    rhos = [r.rho for r in trajectory]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ns, y=rhos, mode="lines", name="ρ_bio",
                             line=dict(color="#00cc99", width=2)))
    fig.update_layout(
        title=title,
        xaxis_title="Choice step n",
        yaxis_title="Knowledge density ρ_bio",
        template="plotly_dark",
    )
    if output_html:
        fig.write_html(output_html)
        print(f"  Saved: {output_html}")
    else:
        fig.show()


def plot_phase_portrait(
    rho_vals: List[float],
    fidelity_vals: List[float],
    title: str = "Phase Portrait: ρ_bio vs Fidelity",
    output_html: Optional[str] = None,
) -> None:
    """Phase portrait — replication density vs fidelity (shows Eigen threshold)."""
    if not PLOTLY_AVAILABLE:
        return
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=fidelity_vals, y=rho_vals,
        mode="lines+markers",
        marker=dict(size=4, color=np.arange(len(rho_vals)),
                    colorscale="Viridis", showscale=True,
                    colorbar=dict(title="Step n")),
        line=dict(width=1, color="rgba(255,255,255,0.2)"),
        name="trajectory",
    ))
    fig.update_layout(
        title=title,
        xaxis_title="Replication fidelity",
        yaxis_title="ρ_bio",
        template="plotly_dark",
    )
    if output_html:
        fig.write_html(output_html)
        print(f"  Saved: {output_html}")
    else:
        fig.show()


def plot_lineage_tree(
    lineages: List[List[float]],
    title: str = "Lineage Trajectories in ρ-Space",
    output_html: Optional[str] = None,
) -> None:
    """
    Each lineage is a list of ρ values over n steps.
    Shows how different replicator lineages diverge/converge in choice-space.
    """
    if not PLOTLY_AVAILABLE:
        return
    fig = go.Figure()
    colors = ["#00cc99", "#ff6b6b", "#4ecdc4", "#ffe66d", "#a8e063",
              "#ff9f43", "#c44dff", "#2dd4bf"]
    for i, lineage in enumerate(lineages):
        fig.add_trace(go.Scatter(
            y=lineage,
            mode="lines",
            name=f"Lineage {i + 1}",
            line=dict(color=colors[i % len(colors)], width=1.5),
            opacity=0.8,
        ))
    fig.update_layout(
        title=title,
        xaxis_title="Choice step n",
        yaxis_title="ρ_bio",
        template="plotly_dark",
    )
    if output_html:
        fig.write_html(output_html)
        print(f"  Saved: {output_html}")
    else:
        fig.show()
