from typing import List, Dict, Union, Callable, Optional, Tuple, Any
from cytodataparser import CytoGateParser
import polars as pl
import numpy as np

#TODO: Allow for multiple groupby columns
def describe_metric(
    cgp: CytoGateParser,
    node: List[str],
    sample_criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]] = None,
    metric: str = "pct_parent",
    groupby: Optional[str] = None,
    quantiles: Tuple[float, float] = (0.25, 0.75),
    return_dataframe: bool = False,
    return_nodes: bool = False
) -> Union[dict, pl.DataFrame]:
    """
    Compute descriptive statistics for a given node metric across samples.
    Optionally groups results by metadata field.
    """
    matches = cgp.get_nodes(node, sample_criteria)

    records = []
    matched_samples = []

    for match in matches:
        meta = match["metadata"]
        for n in match["nodes"]:
            value = n.measures.get(metric)
            if value is not None:

                record = {"value": value, "node_path": n.name}
                record.update(meta)
                records.append(record)
                if return_nodes:
                    matched_samples.append({
                        "metadata": meta,
                        "nodes": match["nodes"]
                    })

    if not records:
        raise ValueError("No valid values found for the specified node(s) and metric.")

    df = pl.DataFrame(records)

    q1, q3 = quantiles
    stats = {
        "gates": sorted(df["node_path"].unique().to_list()),
        "metric": metric,
        "n_total": df.shape[0]
    }

    if groupby and groupby in df.columns:
        grouped = df.group_by(groupby)
        result = []
        group_stats = {}
        for group_name, subdf in grouped:
            vals = subdf["value"].to_numpy()
            group_result = {
                "mean": float(np.mean(vals)),
                "std": float(np.std(vals, ddof=1)),
                "median": float(np.median(vals)),
                f"q{int(100*q1)}": float(np.quantile(vals, q1)),
                f"q{int(100*q3)}": float(np.quantile(vals, q3))
            }
            result.append({
                groupby: group_name,
                "n": len(vals),
                **group_result
            })
            group_stats[group_name] = group_result

        if return_dataframe:
            return pl.DataFrame(result)
        else:
            stats.update({
                "groupby": groupby,
                "n_per_group": {row[groupby]: row["n"] for row in result},
                "group_stats": group_stats
            })
            if return_nodes:
                stats["samples"] = matched_samples
            return stats

    else:
        vals = df["value"].to_numpy()
        stats.update({
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals, ddof=1)),
            "median": float(np.median(vals)),
            f"q{int(100*q1)}": float(np.quantile(vals, q1)),
            f"q{int(100*q3)}": float(np.quantile(vals, q3))
        })
        if return_dataframe:
            return pl.DataFrame([stats])
        else:
            if return_nodes:
                stats["samples"] = matched_samples
            return stats
