from typing import Union, Tuple, List, Dict, Any, Callable, Optional
from cytodataparser import CytoGateParser
import cytodataparser.utils.helpers as helpers
from .helpers import run_posthoc
from scipy.stats import ttest_ind, ttest_rel, f_oneway, kruskal, shapiro, levene, chi2_contingency, pearsonr, spearmanr
import warnings
from collections import defaultdict
import numpy as np

#TODO: Include alpha in function definition and add significance level to results
# TODO: Test passing node as str, List[str] and List[List[str]]
def run_ttest(
    cgp: CytoGateParser,
    node: Union[str, List[str], List[List[str]]],
    sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]] = None,
    metric: str = "pct_parent",
    groupby: Optional[str] = None,
    alpha: float = 0.05,
    flavor: str = "auto",  # Options: "auto", "student", "welch", "paired"
    return_values: bool = False,
    return_nodes: bool = False
) -> dict:
    """
    Run a t-test between groups (unpaired) or between nodes (paired) on a given metric.

    Parameters:
        cgp: CytoGateParser object
        node: A single node (str) for group-wise comparison, or a list of two nodes for within-sample paired test
        sample_criteria: Filter criteria for sample metadata
        metric: Metric to compare (e.g., "Count", "proliferation index")
        groupby: Metadata field to use for grouping (required for unpaired tests)
        flavor: "auto", "student", "welch", or "paired"
        return_values: If True, include raw group or paired values in output

    Returns:
        dict containing test results and metadata
    """
    paired = isinstance(node, list) and len(node) == 2 and all(isinstance(n, list) for n in node)
    results = {
        "test": "t_test",
        "metric": metric,
        "paired": paired,
    }

    if paired:
        results["test"] = "paired_t_test"
        # Paired t-test: compare two nodes within each sample
        node_a = node[0]
        node_b = node[1]
        if isinstance(node_a, str):
            node_a = [node_a]
        if isinstance(node_b, str):
            node_b = [node_b]
        matches_a = cgp.get_nodes(node_a, sample_criteria)
        matches_b = cgp.get_nodes(node_b, sample_criteria)

        indexed_a = {helpers.stringify_metadata(m["metadata"]): m for m in matches_a}
        indexed_b = {helpers.stringify_metadata(m["metadata"]): m for m in matches_b}

        common_keys = set(indexed_a) & set(indexed_b)
        if len(common_keys) < 2:
            raise ValueError("Not enough overlapping samples with both nodes for paired test.")

        values_a, values_b = [], []
        for key in sorted(common_keys):  # sort for deterministic order
            node_a = indexed_a[key]["nodes"][0]
            node_b = indexed_b[key]["nodes"][0]
            val_a = node_a.measures.get(metric)
            val_b = node_b.measures.get(metric)
            if val_a is not None and val_b is not None:
                values_a.append(val_a)
                values_b.append(val_b)

        if len(values_a) < 2:
            raise ValueError("Not enough paired samples to run paired t-test.")

        stat, pval = ttest_rel(values_a, values_b)

        if return_nodes:
            results.update({
                "nodes": [node_a, node_b],
            })
        results.update({
            "n_samples": len(values_a),
            "statistic": float(stat),
            "p_value": float(pval)
        })

        # TODO: Test to make sure node_a[0] really returns the name - used to be node_a.name
        if return_values:
            results["values"] = {
                node_a[0]: values_a,
                node_b[0]: values_b
            }

    else:
        # Unpaired test: compare metric for same node across groups
        if groupby is None:
            raise ValueError("Unpaired t-tests require a `groupby` field for metadata grouping.")
        if isinstance(node, str):
            node = [node]

        matches = cgp.get_nodes(terms=node, sample_criteria=sample_criteria)
        group_vals: Dict[Any, List[float]] = {}

        for match in matches:
            group = match["metadata"].get(groupby)
            for n in match["nodes"]:
                val = n.measures.get(metric)
                if val is not None:
                    group_vals.setdefault(group, []).append(val)

        if len(group_vals) != 2:
            raise ValueError(f"T-test requires exactly two groups; found: {list(group_vals.keys())}")

        (group1, vals1), (group2, vals2) = group_vals.items()

        if flavor == "student":
            results["test"] = "student_t_test"
            stat, pval = ttest_ind(vals1, vals2, equal_var=True)
        elif flavor in ["auto", "welch"]:
            results["test"] = "welch_t_test"
            stat, pval = ttest_ind(vals1, vals2, equal_var=False)
        else:
            raise ValueError(f"Unsupported flavor '{flavor}' for unpaired test.")
        
        if return_nodes:
            results.update({
                "nodes": node,
            })
        results.update({
            "groupby": groupby,
            "groups": [group1, group2],
            "n_per_group": {group1: len(vals1), group2: len(vals2)},
            "statistic": float(stat), # type: ignore
            "p_value": float(pval), # type: ignore
            "group_means": {
                group1: float(np.mean(vals1)),
                group2: float(np.mean(vals2))
            },
            "group_stds": {
                group1: float(np.std(vals1, ddof=1)),
                group2: float(np.std(vals2, ddof=1))
            }
        })

        if return_values:
            results["values"] = {
                group1: vals1,
                group2: vals2
            }
        
    if results["p_value"] < alpha:
        results["significant"] = True
    return results

def run_anova(
    cgp: CytoGateParser,
    node: List[str],
    groupby: Optional[str],
    sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]] = None,
    metric: str = "pct_parent",
    flavor: str = "auto",  # "auto", "anova", "kruskal"
    posthoc: str = "auto",  # "auto", "tukey", "games-howell", "dunn"
    return_values: bool = False,
    verbose: bool = True
) -> dict:
    """
    Run a one-way ANOVA or Kruskal-Wallis test across multiple groups on a node metric.
    Includes assumption checks for normality and equal variances if flavor="auto".
    Optionally performs post-hoc tests if specified.

    Parameters:
        cgp: CytoGateParser object
        node: Node path (list of terms) to analyze
        sample_criteria: Optional metadata filtering
        metric: Metric to extract
        groupby: Metadata field to group by
        flavor: "auto", "anova", or "kruskal"
        posthoc: "auto", "tukey", "games-howell", or "dunn"
        return_values: If True, include raw values in output

    Returns:
        dict with test results
    """
    if groupby is None:
        raise ValueError("ANOVA requires a `groupby` metadata field.")

    matches = cgp.get_nodes(node, sample_criteria)
    group_vals: Dict[str, List[float]] = {}

    for match in matches:
        group = match["metadata"].get(groupby)
        for n in match["nodes"]:
            val = n.measures.get(metric)
            if val is not None:
                group_vals.setdefault(group, []).append(val)

    if len(group_vals) < 3:
        raise ValueError(f"ANOVA requires ≥3 groups; found: {list(group_vals.keys())}")

    groups, values = zip(*group_vals.items())

    # Check assumptions if flavor is "auto"
    test_type = flavor
    assumptions = {"normality": True, "equal_variance": True}
    if flavor == "auto":
        # Normality: require each group to pass Shapiro
        normality_pvals = [shapiro(vals)[1] for vals in values if len(vals) >= 3]
        if any(p < 0.05 for p in normality_pvals):
            assumptions["normality"] = False

        # Equal variances: Levene’s test
        if len(values) >= 2:
            _, p_levene = levene(*values)
            if p_levene < 0.05:
                assumptions["equal_variance"] = False

        if not all(assumptions.values()):
            if not assumptions["normality"]:
                if verbose:
                    warnings.warn(
                        "ANOVA assumption of normality not met. "
                        "Falling back to Kruskal–Wallis test."
                    )
            elif not assumptions["equal_variance"]:
                if verbose:
                    warnings.warn(
                        "ANOVA assumption of equal variances not met. "
                        "Falling back to Kruskal–Wallis test."
                    )
            else:
                if verbose:
                    warnings.warn(
                        "ANOVA assumptions not met (normality and equal variance). "
                        "Falling back to Kruskal–Wallis test."
                    )
            test_type = "kruskal"
        else:
            test_type = "anova"

    # Perform the selected test
    if test_type == "anova":
        stat, pval = f_oneway(*values)
    elif test_type == "kruskal":
        stat, pval = kruskal(*values)
    else:
        raise ValueError(f"Unsupported flavor: '{flavor}'")

    result = {
        "test": test_type,
        "gate": sorted({n.name for match in matches for n in match["nodes"]}),
        "metric": metric,
        "groupby": groupby,
        "groups": list(groups),
        "n_per_group": {g: len(v) for g, v in group_vals.items()},
        "statistic": float(stat),
        "p_value": float(pval),
        "group_means": {g: float(np.mean(v)) for g, v in group_vals.items()},
        "group_stds": {g: float(np.std(v, ddof=1)) for g, v in group_vals.items()},
    }

    if flavor == "auto":
        result["assumption_check"] = assumptions

    if return_values:
        result["values"] = group_vals

    # Run post-hoc analysis if requested
    if posthoc:
        if pval > 0.05:
            if verbose:
                warnings.warn(
                    f"ANOVA p-value was not statistically significant (p = {pval:.3f}). "
                    "Post-hoc results may not be meaningful."
                )
        result["posthoc"] = run_posthoc(group_vals, test_type, posthoc, assumptions)

    return result

def run_chi2_test(
    cgp: CytoGateParser,
    row_field: str,
    col_field: str,
    sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]] = None,
    correction: bool = True
) -> dict:
    """
    Perform chi-squared test of independence between two metadata fields.

    Parameters:
        cgp: CytoGateParser
        row_field: Metadata field to use for rows
        col_field: Metadata field to use for columns
        sample_criteria: Optional filtering of samples by metadata
        correction: Whether to apply Yates' correction (default: True)

    Returns:
        Dict summarizing chi-squared test result and contingency table.
    """
    # Filter samples
    samples = [
        s["metadata"] for s in cgp.samples
        if sample_criteria is None or all(
            f(v) if callable(f := sample_criteria[k]) else (v in f if isinstance(f, range) else v == f)
            for k, v in s["metadata"].items() if k in sample_criteria
        )
    ]

    if not samples:
        raise ValueError("No samples matched the sample_criteria.")

    # Build contingency counts
    contingency = defaultdict(lambda: defaultdict(int))
    row_categories = set()
    col_categories = set()

    for meta in samples:
        row_val = meta.get(row_field)
        col_val = meta.get(col_field)
        if row_val is not None and col_val is not None:
            contingency[row_val][col_val] += 1
            row_categories.add(row_val)
            col_categories.add(col_val)

    row_categories = sorted(row_categories)
    col_categories = sorted(col_categories)

    # Convert to observed matrix
    observed = np.array([
        [contingency[r][c] for c in col_categories]
        for r in row_categories
    ])

    if observed.shape[0] < 2 or observed.shape[1] < 2:
        raise ValueError("Chi-squared test requires at least a 2x2 table.")

    # Perform the test
    chi2, p, dof, expected = chi2_contingency(observed, correction=correction)

    # Evaluate expected counts
    low_expected = expected < 5 # type: ignore
    num_low = np.sum(low_expected)
    total = expected.size # type: ignore
    percent_low = num_low / total

    if num_low > 0:
        warnings.warn(
            f"Chi-squared test assumption warning: "
            f"{num_low} of {total} cells ({percent_low:.1%}) have expected counts < 5. "
            "Interpret the p-value with caution."
        )

    return {
        "test": "chi2",
        "row_field": row_field,
        "col_field": col_field,
        "row_labels": row_categories,
        "col_labels": col_categories,
        "chi2": float(chi2), # type: ignore
        "p_value": float(p), # type: ignore
        "degrees_of_freedom": int(dof), # type: ignore
        "observed": observed.tolist(),
        "expected": expected.tolist() # type: ignore
    }

def run_correlation(
    cgp: CytoGateParser,
    node_a: List[str],
    node_b: Optional[List[str]] = None,
    metric: Union[str, List[str]] = "Count",
    sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]] = None,
    method: str = "pearson",
    return_values: bool = False
) -> dict:
    """
    Compute correlation between two node metrics across samples.

    Supports two modes:
    1. Compare one metric across two nodes (provide node_a, node_b, and metric:str)
    2. Compare two metrics within one node (provide node_a and metric:Tuple[str, str])

    Parameters:
        cgp: CytoGateParser
        node_a: Path to first node (or the shared node if comparing two metrics)
        node_b: Optional path to second node
        metric: Metric name (str) or a tuple of two metric names
        sample_criteria: Optional filter for samples
        method: "pearson" or "spearman"
        return_values: If True, include raw values in output

    Returns:
        dict with correlation result and raw values
    """

    if method not in {"pearson", "spearman"}:
        raise ValueError("Method must be 'pearson' or 'spearman'.")

    is_metric_pair = isinstance(metric, List) and len(metric) == 2
    is_node_pair = node_b is not None and isinstance(metric, str)

    label_a = ""
    label_b = ""

    if not (is_metric_pair or is_node_pair):
        raise ValueError(
            "You must provide either two nodes and one metric, or one node and two metrics."
        )

    x, y = [], []

    if is_metric_pair:
        metric_a, metric_b = metric[0], metric[1]
        matches = cgp.get_nodes(node_a, sample_criteria)

        for match in matches:
            node = match["nodes"][0]
            val1 = node.measures.get(metric_a)
            val2 = node.measures.get(metric_b)
            if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                x.append(val1)
                y.append(val2)

        label_a = f"{'/'.join(node_a)}:{metric_a}"
        label_b = f"{'/'.join(node_a)}:{metric_b}"

    elif is_node_pair:
        matches_a = cgp.get_nodes(node_a, sample_criteria)
        matches_b = []
        if node_b is not None:
            matches_b = cgp.get_nodes(node_b, sample_criteria)

        index_a = {helpers.stringify_metadata(m["metadata"]): m for m in matches_a}
        index_b = {helpers.stringify_metadata(m["metadata"]): m for m in matches_b}
        shared_keys = sorted(set(index_a) & set(index_b))

        for key in shared_keys:
            node1 = index_a[key]["nodes"][0]
            node2 = index_b[key]["nodes"][0]
            val1 = node1.measures.get(metric)
            val2 = node2.measures.get(metric)
            if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                x.append(val1)
                y.append(val2)

        label_a = "/".join(node_a)
        if node_b is not None:
            label_b = "/".join(node_b)

    if len(x) < 3:
        raise ValueError("Too few valid paired samples to compute correlation.")

    corr_func = pearsonr if method == "pearson" else spearmanr
    corr, pval = corr_func(x, y)

    results = {
        "test": "correlation",
        "method": method,
        "gates": [label_a, label_b],
        "metrics": metric if isinstance(metric, str) else f"{metric[0]} vs {metric[1]}",
        "n": len(x),
        "correlation": float(corr), # type: ignore
        "p_value": float(pval), # type: ignore
        
    }
    if return_values:
        results["values"] = {
            "x": x,
            "y": y
        }
    return results