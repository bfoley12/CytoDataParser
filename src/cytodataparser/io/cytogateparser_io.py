from cytodataparser.core import CytoGateParser
import polars as pl
from typing import Optional
from pathlib import Path
import json
from datetime import date, datetime
from cytodataparser.structures.tree import GateTree


def load_from_xlsx(file_path: str, sheet_name: Optional[str] = None) -> CytoGateParser:
    """
    Load a CytoGateParser from an Excel file.

    Parameters:
        file_path (str): Path to the Excel file.
        sheet_name (str): Name of the sheet containing gate data.

    Returns:
        CytoGateParser: An instance of CytoGateParser loaded with data from the specified sheet.
    """
    df = pl.read_excel(file_path, sheet_name=sheet_name)
    return CytoGateParser(df)

def load_from_csv(file_path: str) -> CytoGateParser:
    """
    Load a CytoGateParser from a CSV file.

    Parameters:
        file_path (str): Path to the CSV file.

    Returns:
        CytoGateParser: An instance of CytoGateParser loaded with data from the CSV file.
    """
    df = pl.read_csv(file_path)
    return CytoGateParser(df)

def load_from_json(path: str | Path) -> CytoGateParser:
    with open(path, "r") as f:
        data = json.load(f)

    restored_samples = []
    for sample in data["samples"]:
        restored_samples.append({
            "metadata": sample["metadata"],
            "tree": GateTree.from_dict(sample["tree"])
        })

    return CytoGateParser.from_samples(
        samples=restored_samples,
        metadata_cols=data.get("metadata_cols")
    )

def _sanitize(obj):
    """Recursively convert unsupported types to JSON-safe formats."""
    if isinstance(obj, (date, datetime)):
        return obj.isoformat()
    elif isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    else:
        return obj  # leave as is (str, int, float, etc.)

def save_to_json(cgp: CytoGateParser, path: str | Path):
    sanitized_samples = []

    for sample in cgp.samples:
        tree = sample["tree"]
        if hasattr(tree, "to_dict"):
            tree = tree.to_dict()  # Convert GateTree -> dict
        sanitized_samples.append({
            "metadata": _sanitize(sample["metadata"]),
            "tree": _sanitize(tree)
        })

    data = {
        "samples": sanitized_samples,
        "metadata_cols": cgp.metadata_cols,
        "version": 1
    }

    with open(path, "w") as f:
        json.dump(data, f, indent=2)