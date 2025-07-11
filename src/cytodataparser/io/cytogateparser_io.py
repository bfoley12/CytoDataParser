import polars as pl
from typing import Optional, List, Union
from pathlib import Path
import json
from datetime import date, datetime
from cytodataparser.structures import GateTree, Sample
import flowkit as fk
import sys
import warnings

# TODO: for load_wsp, path should point directly to the .wsp - need to make this more clear
def load_file(path: Union[str, Path], 
              sheet_name: Optional[Union[str, None]]=None, 
              fcs_files: Union[str, List[str], Path, List[Path]]="./Unmixed",
              verbose: Optional[bool]=True) -> List[Sample]:
    """
    Load a CytoGateParser from any of: xlsx, xls, csv, or json

    Parameters:
        path (str): Path to the Excel file.
        sheet_name (str): Name of the sheet containing gate data.
        fcs_files (Union[str, List[str], Path, List[Path]]): If str or path, this must be the top level folder containing all fcs files. load_wsp will search all subfolders (excluding "Reference Group") for fcs files
            If List[str] or List[Path], depending on if the entries end in .fcs will be treated as directories to search or as files to load
    Returns:
        List[Sample]: A list of samples to be loaded into a CytoGateParser.
    """

    df = pl.DataFrame()
    if isinstance(path, str):
        path = Path(path)

    ending = path.suffix

    # TODO: make load_excel and load_csv methods to unify interface
    if ending == ".xlsx" or ending == ".xls":
        df = pl.read_excel(path, sheet_name=sheet_name)
    elif ending == ".csv":
        df = pl.read_csv(path)
    elif ending == ".json":
        return load_json(path)
    elif ending == ".wsp":
        return load_wsp(path, fcs_files, verbose=verbose)
    else:
        raise ValueError("Unexpected filetype encountered. Please use one of xlsx, xls, csv, or json")

    samples = samples_from_polars(df)
    return samples

def load_wsp(
    path: Union[str, Path],
    fcs_files: Union[str, List[str], Path, List[Path]] = "./Unmixed",
    verbose: Optional[bool]=True
) -> List[Sample]:
    """
    Load Samples from a FlowJo Workspace file and associated FCS files.

    Parameters:
        path (Union[str, Path]): Path to the FlowJo .wsp file.
        fcs_files (Union[str, List[str], Path, List[Path]]): If str or Path, it is a top-level folder
            where all FCS files are stored. If a list, each entry is either a directory (searched recursively,
            skipping 'Reference Group') or a specific .fcs file.

    Returns:
        List[dict]: A list of dicts with 'metadata' and 'tree' (gated counts) for each sample.
    """
    if not verbose:
        warnings.filterwarnings("ignore", category=UserWarning)

    wsp_path = Path(path)

    if isinstance(fcs_files, (str, Path)):
        fcs_files = [fcs_files] # type: ignore

    # Resolve all paths relative to the user's script
    resolved_fcs_paths = [resolve_relative_to_entrypoint(p) for p in fcs_files] # type: ignore

    # Find FCS files
    fcs_paths = []
    for item in resolved_fcs_paths:
        item_path = Path(item)
        if item_path.is_dir():
            fcs_paths.extend([
                str(fcs) for fcs in item_path.rglob("*.fcs")
                if "Reference Group" not in fcs.parts
            ])
        elif item_path.is_file() and item_path.suffix.lower() == ".fcs":
            fcs_paths.append(str(item_path))

    # Map FCS base name to full path (case-insensitive match)
    fcs_map = {fcs.lower(): fcs for fcs in fcs_paths}

    output = []

    for sample_path in fcs_map.values():
        # TODO: Implement disk checking
        #fcs_name = Path(sample.fcs_file).name.lower()
        #if fcs_name not in fcs_map:
        #    print(f"FCS file {fcs_name} not found on disk. Skipping.")
        #    continue
        print(sample_path)
        sample_wsp = fk.Workspace(wsp_path, sample_path)
        sample = sample_wsp.get_samples()
        if sample == []:
            print(f"Unable to load fcs file: {sample_path}")
            continue
        sample = sample[0]
        
        print(sample)

        # Extract metadata
        metadata = sample.get_metadata()
        # Extract gating tree (counts per gate)
        sample_wsp.analyze_samples()
        result = sample_wsp.get_analysis_report()

        output.append(*samples_from_wsp(sample_wsp))

    return output

def samples_from_wsp(wsp: fk.Workspace) -> List[Sample]:
    # Prep table to make GateTree
    metadata = pl.DataFrame([wsp.get_keywords(sample.id) for sample in wsp.get_samples()])
    wsp.analyze_samples()
    summary = wsp.get_analysis_report()

    # Format gate_path to work with GateNode
    summary['full_path'] = summary.apply(
        lambda row: "/".join(row['gate_path'] + (row['gate_name'],)),
        axis=1
    )

    # Pivot and format to work with GateNode/GateTree
    summary = summary.pivot(index = "sample", columns="full_path", values=["count", "absolute_percent", "relative_percent"])

    summary.columns = [f"{col[1]} | {col[0]}" for col in summary.columns]

    # Optional: reset index if you want 'sample' as a column
    summary = summary.reset_index()
    #summary['sample'] = summary['sample'].str.replace(r'\.fcs$', '', case=False, regex=True)
    summary = pl.DataFrame(summary)

    # Add on metadata
    summary = summary.join(metadata, left_on="sample", right_on = "$FIL")

    return samples_from_polars(summary)


def samples_from_polars(data: pl.DataFrame) -> List[Sample]:
    """
    Format Samples from a polars dataframe

    Parameters:
        data (pl.DataFrame): a polars dataframe, with data columns denoted by *gate_path* | *metric* and metadata columns lacking "|"

    Returns:
        List[Sample]: A list of samples to be loaded into a CytoGateParser
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

def load_json(path: Union[str, Path]) -> List[Sample]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    try:
        with open(path, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Malformed JSON file: {e}") from e

    if not isinstance(data, dict) or "samples" not in data:
        raise ValueError("JSON must contain a top-level 'samples' key.")

    if not isinstance(data["samples"], list):
        raise ValueError("'samples' must be a list.")

    restored_samples = []
    for i, sample in enumerate(data["samples"]):
        if not isinstance(sample, dict):
            raise ValueError(f"Sample at index {i} is not a dictionary.")

        if "metadata" not in sample or "tree" not in sample:
            raise ValueError(f"Sample at index {i} must contain 'metadata' and 'tree' keys.")

        restored_samples.append(
            Sample(
                metadata=_sanitize(sample["metadata"]),
                tree=_sanitize(GateTree.from_dict(sample["tree"]))
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

def resolve_relative_to_entrypoint(path: Union[str, Path]) -> Path:
    """Resolve a path relative to the user's script entrypoint. In interactive mode, use cwd."""
    path = Path(path)
    if path.is_absolute():
        return path

    try:
        main_module = sys.modules["__main__"]
        entry_dir = Path(getattr(main_module, "__file__", Path.cwd())).resolve()
        return (entry_dir / path).resolve()
    except Exception:
        # Absolute fallback: current working directory
        return path.resolve()