import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Dict

def violin_helper(
    df: pd.DataFrame,
    x: str,
    y: str,
    group: str,
    color_map: dict,
    show_points: bool = False,
    jitter_strength: float = 0.01,
):
    """
    Create one violin per x value, colored by a group, without side-by-side splitting.
    """
    fig = go.Figure()
    rng = np.random.default_rng(42)

    for xi in df[x].unique():
        sub_df = df[df[x] == xi]
        if sub_df.empty:
            continue
        group_val = sub_df[group].iloc[0]
        color = color_map.get(group_val, "gray")

        fig.add_trace(go.Violin(
            y=sub_df[y],
            x=[xi] * len(sub_df),
            line_color=color,
            name=str(xi),
            showlegend=False
        ))

        # Optional scatter (jittered) points
        if show_points:
            jitter_x = rng.normal(loc=0, scale=jitter_strength, size=len(sub_df))
            x_jittered = [str(xi)] * len(sub_df)  # Keep as str for Plotly alignment
            fig.add_trace(go.Scatter(
                x=[str(xi)] * len(sub_df),
                y=sub_df[y],
                mode="markers",
                marker=dict(color=color, size=6, opacity=0.7),
                showlegend=False,
                hoverinfo="x+y"
            ))

    fig.update_layout(
        xaxis_title=x,
        yaxis_title=y
    )

    return fig

def rename_legend(fig: go.Figure, new_names: Dict[str, str]):
    fig.for_each_trace(lambda t: t.update(name = new_names[t.name],
                                        legendgroup = new_names[t.name],
                                        hovertemplate = t.hovertemplate.replace(t.name, new_names[t.name])
                                        )
                    )
    return fig