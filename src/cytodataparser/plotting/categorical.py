import plotly.express as px
import plotly.graph_objects as go
import polars as pl
import pandas as pd
import numpy as np
from typing import Optional, Dict, Union, Callable, Any, List, Tuple
import warnings

from .themes import get_color_map
from .helpers import violin_helper, rename_legend
from cytodataparser import CytoGateParser
from cytodataparser.utils import helpers
from cytodataparser.analysis import run_ttest, run_anova


def pval_to_star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'


def extract_pairwise_pvals(anova_output: dict) -> Dict[Tuple[str, str], float]:
    """
    Extracts pairwise p-values from the posthoc section of a run_anova() result.
    Output: {("group1", "group2"): p_value}
    """
    pvals = {}

    if "posthoc" in anova_output and "comparisons" in anova_output["posthoc"]:
        for comparison in anova_output["posthoc"]["comparisons"]:
            g1 = comparison["group1"]
            g2 = comparison["group2"]
            p = comparison["p_value"]
            pvals[(g1, g2)] = p
    else:
        print("Warning: No posthoc comparisons found in ANOVA output.")

    return pvals

def extract_ttest_pval(ttest_result: dict) -> Dict[Tuple[str, str], float]:
    """
    Extracts a single pairwise p-value from run_ttest() result.
    Returns: {("group1", "group2"): p_value} or empty dict if invalid
    """
    try:
        g1, g2 = map(str, ttest_result["groups"])
        p = ttest_result["p_value"]
        if p is not None and not np.isnan(p):
            return {(g1, g2): p}
        else:
            print("Warning: t-test returned NaN p-value.")
            return {}
    except Exception as e:
        print(f"Error extracting t-test p-value: {e}")
        return {}

def add_significance(
    fig: go.Figure,
    group_labels: List[str],
    p_values: Dict[Tuple[str, str], float],
    y_max: float,
    y_pad: float = 0.05,
    star_func=None
):
    """
    Add significance lines/stars to a grouped bar or categorical figure.

    Parameters:
    - fig: Plotly Figure (assumed to be bar, violin, or box)
    - group_labels: Ordered labels on x-axis (e.g., ["Control", "Treatment A", "Treatment B"])
    - p_values: Dict of (group1, group2): p-value
    - y_max: Current max y value in data (used as base for placing annotations)
    - y_pad: Padding above max y value per annotation layer
    - star_func: Optional function to convert p-value to annotation (e.g., "***")
    """

    if star_func is None:
        star_func = pval_to_star

    # Build mapping from group label -> x positions
    x_lookup = {label: [] for label in group_labels}

    for trace in fig.data:
        if hasattr(trace, "x") and hasattr(trace, "y") and trace.x is not None:
            for x in trace.x:
                x_str = str(x)
                if x_str in x_lookup:
                    x_lookup[x_str].append(x)
                elif isinstance(x, (int, float)) and str(int(x)) in x_lookup:
                    x_lookup[str(int(x))].append(x)

    # Get representative x positions (e.g., bar center per category)
    x_coord = {}
    for label, positions in x_lookup.items():
        if positions:
            x_coord[label] = np.mean(positions)
        else:
            x_coord[label] = group_labels.index(label)  # fallback to category index

    used_levels = 0
    for (g1, g2), p in p_values.items():
        if g1 not in x_coord or g2 not in x_coord:
            continue

        x1, x2 = x_coord[g1], x_coord[g2]
        if x1 > x2:
            x1, x2 = x2, x1

        level = used_levels
        y = y_max + y_pad * (level + 1)
        y_star = y + y_pad / 2

        # Horizontal line
        fig.add_shape(
            type="line", x0=x1, x1=x2, y0=y, y1=y,
            line=dict(color="black", width=1), xref="x", yref="y"
        )

        # Vertical ticks
        for x in [x1, x2]:
            fig.add_shape(
                type="line", x0=x, x1=x, y0=y, y1=y - y_pad * 0.2,
                line=dict(color="black", width=1), xref="x", yref="y"
            )

        # Annotation
        fig.add_annotation(
            x=(x1 + x2) / 2, y=y_star,
            text=star_func(p), showarrow=False,
            yanchor="bottom", font=dict(size=12),
            xref="x", yref="y"
        )

        used_levels += 1


def categorical_plot(
    cgp: CytoGateParser,
    node: Union[List[str], str],
    x: str,
    y: str = "pct_parent",
    sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]] = None,
    legend_names: Optional[Dict[str, str]] = None,
    color: Optional[str] = None,
    facet_row: Optional[str] = None,
    facet_col: Optional[str] = None,
    label_sig: bool = True,
    plot_type: str = "box",
    show_points: bool = False,
    jitter_strength: float = 0.01,  # currently unused
    color_map: Optional[dict] = None,
    agg_func: str = "mean",
    error: str = "sem",
    verbose: bool = True,
    **kwargs
):
    """
    Plot a specific measure from a specific GateNode across samples in a CytoGateParser.
    Smart defaults ensure clean and interpretable plots for bench scientists.
    """

    if isinstance(node, str):
        node = [node]

    # Get and flatten data
    samples = cgp.get_nodes(terms=node, sample_criteria=sample_criteria)
    if not samples:
        raise ValueError(f"No samples found for node: {node}")
    message = ""
    nodes_found = ""
    first_node_name = ""
    seen_node_name_sets = set()

    if color != "node" and facet_row != "node" and facet_col != "node":
        for sample in samples:
            node_list = sample["nodes"]
            if len(node_list) > 1:
                if not message:
                    message = f"Multiple nodes found when using terms: {node}. Make your search terms more specific.\n"

            # Extract names and use them to track uniqueness
            node_names = tuple(n.name for n in node_list)
            if node_names and node_names not in seen_node_name_sets:
                if not first_node_name:
                    first_node_name = node_names[0]
                    message += f"Using first node: {first_node_name}\n"
                seen_node_name_sets.add(node_names)
                nodes_found += f"{[n.name for n in node_list]}\n"

    if message and verbose:
        message += f"Nodes found: {nodes_found.strip()}\n"
        warnings.warn(message, category=RuntimeWarning)

    df = helpers.flatten_samples(samples)
    if message:
        df = df.filter(pl.col("node") == first_node_name)
    pdf = df.to_pandas()

    if color not in pdf.columns:
        color = None
        warnings.warn(message="Color column not found in DataFrame", category=RuntimeWarning)
    # Clean categorical ordering
    pdf[x] = pd.Categorical(pdf[x], categories=sorted(pdf[x].unique()), ordered=True)
    if color:
        pdf[color] = pd.Categorical(pdf[color], categories=sorted(pdf[color].unique()), ordered=True)
        if color_map is None:
            color_map = get_color_map(color, pdf[color].unique().tolist())
            
    plot_x = x
    if "category_orders" not in kwargs and color:
        kwargs["category_orders"] = {
            plot_x: pdf[plot_x].cat.categories.tolist(),
            color: pdf[color].cat.categories.tolist()
        }

    # Main plot
    if plot_type == "box":
        fig = px.box(
            pdf, x=plot_x, y=y, color=color,
            facet_row=facet_row, facet_col=facet_col,
            points="all" if show_points else "outliers",
            color_discrete_map=color_map,
            **kwargs
        )

    elif plot_type == "violin":
        fig = px.violin(
            pdf, x=plot_x, y=y, color=color,
            facet_row=facet_row, facet_col=facet_col,
            points="all" if show_points else "outliers",
            color_discrete_map=color_map,
            **kwargs
        )

    # TODO: implement grouped_bar better
    elif plot_type == "grouped_bar" or plot_type == "bar":
        group_cols = [x]
        if color:
            group_cols.append(color)

        grouped = pdf.groupby(group_cols, observed=True)[y].agg([agg_func]).reset_index()

        if error in {"sem", "std"}:
            err_func = lambda s: s.std() / np.sqrt(len(s)) if error == "sem" else s.std()
            grouped["error"] = pdf.groupby(group_cols, observed=True)[y].agg(err_func).reset_index(drop=True)
        metadata_cols = [col for col in pdf.columns if col not in group_cols + [y]]
        meta_df = pdf[group_cols + metadata_cols].drop_duplicates(subset=group_cols)

        # Merge back the metadata
        final_df = pd.merge(grouped, meta_df, on=group_cols, how="left")

        fig = px.bar(
            final_df,
            x=x,
            y=agg_func,
            color=color,
            error_y="error" if "error" in grouped else None,
            barmode="group" if plot_type == "grouped_bar" else "relative",
            facet_row=facet_row,
            facet_col=facet_col,
            color_discrete_map=color_map,
            **kwargs
        )
        y_data = final_df[agg_func]
        y_err = final_df["error"] if "error" in final_df else 0
        y_max = (y_data + y_err).max()
        fig.update_yaxes(range=[0, y_max * 1.25])
    else:
        raise ValueError(f"Unsupported plot_type: {plot_type}")

    # Axis labels
    fig.update_layout(
        xaxis_title=x.replace("_", " ").title(),
        yaxis_title=y.replace("_", " ").title(),
        margin=dict(t=40, b=40, l=40, r=40),
    )

    # Auto significance annotation
    if label_sig and isinstance(fig, go.Figure):
        group_labels = pdf[plot_x].cat.categories.tolist()
        y_max = pdf[y].max()

        try:
            p_map = {}

            if pdf[x].nunique() == 2:
                if verbose:
                    print("Running t-test")
                t_result = run_ttest(cgp, node, groupby=x, metric=y)
                p_map = extract_ttest_pval(t_result)

            elif pdf[x].nunique() > 2:
                if verbose:
                    print("Running ANOVA")
                a_result = run_anova(cgp, node=node, groupby=x, sample_criteria=sample_criteria, metric=y)
                p_map = extract_pairwise_pvals(a_result)

            if not p_map or all(p >= 0.05 for p in p_map.values()):
                fig.add_annotation(
                    text="<b>No significant differences</b>",
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.98,
                    y=0.98,
                    xanchor="right",
                    yanchor="top",
                    font=dict(size=13, color="dimgray"),
                    align="right",
                    bgcolor="rgba(255,255,255,0.85)",
                    bordercolor="lightgray",
                    borderwidth=1,
                    borderpad=6,
                    opacity=0.95
                )

            else:
                add_significance(fig, group_labels, p_map, y_max)

        except Exception as e:
            print(f"Warning: significance annotation failed — {e}")

    n_groups = pdf[plot_x].nunique()

    fig.update_layout(
        width=max(400, min(200 + 150 * n_groups, 1000)),
        margin=dict(t=100, b=40, l=40, r=40)
    )

    # Uncommenting causes boxes and violins to stack instead of group
    # fig.update_traces(width=0.5)        
    
    if legend_names is not None:
        fig = rename_legend(fig, legend_names)
    return fig
