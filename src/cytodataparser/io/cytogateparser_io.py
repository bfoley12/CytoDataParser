import polars as pl
from typing import Optional, List, Union
from pathlib import Path
import json
from datetime import date, datetime
from cytodataparser.structures import GateTree, Sample

# TODO: Implement json loading
def load_file(path: str | Path, sheet_name: Optional[Union[str, None]]=None) -> List[Sample]:
    """
    Load a CytoGateParser from any of: xlsx, xls, csv, or json

    Parameters:
        file_path (str): Path to the Excel file.
        sheet_name (str): Name of the sheet containing gate data.

    Returns:
        List[Sample]: A list of samples to be loaded into a CytoGateParser.
    """

    df = pl.DataFrame()
    if isinstance(path, str):
        path = Path(path)
    ending = path.suffix
    if ending == ".xlsx" or ending == ".xls":
        df = pl.read_excel(path, sheet_name=sheet_name)
    elif ending == ".csv":
        df = pl.read_csv(path)
    elif ending == ".json":
        return load_json(path)
    else:
        raise ValueError("Unexpected filetype encountered. Please use one of xlsx, xls, csv, or json")

    samples = samples_from_polars(df)
    return samples

def samples_from_polars(data: pl.DataFrame) -> List[Sample]:
    """
    Load a CytoGateParser from an Excel file.

    Parameters:
        file_path (str): Path to the Excel file.
        sheet_name (str): Name of the sheet containing gate data.

    Returns:
        CytoGateParser: An instance of CytoGateParser loaded with data from the specified sheet.
    """

    metadata = [col for col in data.columns if '|' not in col]
    metadata = data[metadata]
    samples_prep = [
        {
            "metadata": {k: v[0] if isinstance(v, list) and len(v) == 1 else v
                            for k, v in metadata[row_idx].to_dict(as_series=False).items()},
            "tree": GateTree(row)
        }
        for row_idx, row in enumerate(data.iter_rows(named=True))
    ]
    samples = []
    for sample in samples_prep:
        samples.append(Sample(sample["metadata"], sample["tree"]))

    return samples

# TODO: Actually implement
def json_to_polars(path: str | Path) -> pl.DataFrame:
    # TODO: raise error if json is malformed
    return pl.DataFrame()

def load_json(path: str | Path) -> List[Sample]:
    with open(path, "r") as f:
        data = json.load(f)

    restored_samples = []
    for sample in data["samples"]:
        restored_samples.append(
            Sample(
                metadata=sample["metadata"],
                tree=GateTree.from_dict(sample["tree"])
            )
        )

    return restored_samples

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

# TODO: Move to core or helpers
def save_to_json(cgp, path: str | Path):
    sanitized_samples = []

    for sample in cgp.samples:
        tree = sample.tree
        if hasattr(tree, "to_dict"):
            tree = tree.to_dict()  # Convert GateTree -> dict
        sanitized_samples.append({
            "metadata": _sanitize(sample.metadata),
            "tree": _sanitize(tree)
        })

    data = {
        "samples": sanitized_samples,
    }

    with open(path, "w") as f:
        json.dump(data, f, indent=2)