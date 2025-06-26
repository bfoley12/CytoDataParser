from typing import Dict, List
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import shapiro, levene, f_oneway, kruskal
import numpy as np
import pingouin as pg
import scikit_posthocs as sp

def run_posthoc(group_vals: Dict[str, List[float]], test_type: str, posthoc: str, assumptions: Dict[str, bool]) -> dict:
    import pandas as pd

    # Flatten into long-format DataFrame
    data = [(val, group) for group, vals in group_vals.items() for val in vals]
    df = pd.DataFrame(data, columns=["value", "group"])

    # Select test automatically
    if posthoc == "auto":
        if test_type == "kruskal":
            posthoc = "dunn"
        elif not assumptions["equal_variance"]:
            posthoc = "games-howell"
        else:
            posthoc = "tukey"

    comparisons = []

    if posthoc == "tukey":
        tukey = pairwise_tukeyhsd(df["value"], df["group"])
        for i in range(len(tukey.meandiffs)):
            comparisons.append({
                "group1": tukey.groupsunique[tukey._multicomp.pairindices[i][0]],
                "group2": tukey.groupsunique[tukey._multicomp.pairindices[i][1]],
                "p_value": float(tukey.pvalues[i]),
                "mean_diff": float(tukey.meandiffs[i])
            })

    elif posthoc == "games-howell":
        if pg is None:
            raise ImportError("'pingouin' must be installed to run Games-Howell posthoc test.")
        gh = pg.pairwise_gameshowell(dv="value", between="group", data=df)
        for _, row in gh.iterrows():
            comparisons.append({
                "group1": row["A"],
                "group2": row["B"],
                "p_value": float(row["pval"]),
                "mean_diff": float(row["mean(A)-mean(B)"])
            })

    elif posthoc == "dunn":
        if sp is None:
            raise ImportError("'scikit-posthocs' must be installed to run Dunn's test.")
        dunn = sp.posthoc_dunn(df, val_col="value", group_col="group", p_adjust="bonferroni")
        for i in range(len(dunn.columns)):
            for j in range(i+1, len(dunn.columns)):
                group1 = dunn.columns[i]
                group2 = dunn.columns[j]
                comparisons.append({
                    "group1": group1,
                    "group2": group2,
                    "p_value": float(dunn.iloc[i, j]),
                    "mean_diff": float(np.mean(group_vals[group1]) - np.mean(group_vals[group2]))
                })

    else:
        raise ValueError(f"Unsupported posthoc test: {posthoc}")

    return {
        "test": posthoc,
        "comparisons": comparisons
    }