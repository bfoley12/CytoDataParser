import polars as pl
from collections import defaultdict
from typing import List
from cytodataparser.structures import GateNode

def flatten_samples(samples: List[GateNode]) -> pl.DataFrame:
    """
    Flattens a list of samples into a Polars DataFrame.

    Each sample must be a dict with keys:
    - "metadata": a dict of sample-level metadata
    - "nodes": a list of GateNode objects with .name and .measures

    Returns:
        A polars DataFrame where each row represents a GateNode
        annotated with the sample's metadata.
    """
    rows = []

    for sample in samples:
        metadata = sample["metadata"]
        nodes = sample["nodes"]
        for node in nodes:
            row = {
                **node.flatten(),
                **metadata
            }
            rows.append(row)

    return pl.DataFrame(rows)

def stringify_metadata(meta: dict) -> str:
    """Stable way to create a unique key for metadata dict."""
    return "|".join(f"{k}={meta[k]}" for k in sorted(meta))