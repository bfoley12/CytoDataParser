from __future__ import annotations
from typing import Optional, List, Dict, Any, Callable, Union
import polars as pl
from .structures import GateTree, GateNode, Sample
from .utils.predicates import parse_string_condition, from_range
from datetime import date, datetime

#TODO: Allow for loading with a .wsp file instead of polars dataframe
class CytoGateParser:
    """
    Main class for managing cytometry gating data.
    Parses a Polars DataFrame into metadata and gating trees, one per sample.
    """

    def __init__(self, samples: List[Sample], original_df: pl.DataFrame = None, metadata_cols=None):
        """
        Initialize the CytoGateParser.

        Parameters:
            df (pl.DataFrame): Full dataset containing both metadata and gating results.
            metadata_cols (Optional[List[str]]): List of columns to treat as metadata.
                                             If None, inferred by excluding columns with '|' in name.
        """
        self.original_df = original_df
        self.metadata_cols = metadata_cols if metadata_cols else self._gather_metadata(samples)
        self.samples: List[Sample] = samples

    # TODO: Implement loading from xlsx, csv, and polars dataframe
    @classmethod
    def from_xlsx(cls, path: str) -> CytoGateParser:
        """
        Construct a CytoGateParser from an xlsx file
        """
        data = pl.read_excel(path)
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
        return cls(samples, data, metadata.columns)
    @classmethod
    def from_samples(cls, samples: List[Dict[str, Any]], original_df: Optional[pl.DataFrame] = None) -> CytoGateParser:
        """
        Construct a CytoGateParser from a list of samples.
        """

        return cls(samples, original_df)

    @staticmethod
    def _flatten_tree_to_row(tree: Any) -> dict:
        """
        Flatten a GateTree into a flat row of {path | metric: value}.
        """
        result = {}

        def recurse(node: Any, path: List[str]):
            full_path = "/".join(path + [node.name]) if path else node.name
            for measure_name, val in node.measures.items():
                result[f"{full_path} | {measure_name}"] = val
            for child in node.children:
                recurse(child, path + [node.name])

        recurse(tree.root, [])
        return result
    
    def _gather_metadata(self, samples: List[Sample]):
        all_metadata_cols = set()
        for sample in samples:
            all_metadata_cols.update(sample.metadata.keys())

        return list(all_metadata_cols)

    # Out of date: Used when loading from polars dataframe. Maybe keep?
    def _infer_metadata_cols(self, df: pl.DataFrame) -> List[str]:
        """Infer metadata columns as those that do not contain '|' in their names."""
        return [col for col in df.columns if '|' not in col]

    def _build_samples(self) -> List[Dict[str, Any]]:
        """
        Build a list of samples, each represented as a dictionary with metadata and a GateTree.

        Returns:
            List[Dict[str, Any]]: Each sample contains {'metadata': ..., 'tree': ...}
        """
        return [
            {
                "metadata": {k: v[0] if isinstance(v, list) and len(v) == 1 else v
                              for k, v in self.metadata[row_idx].to_dict(as_series=False).items()},
                "tree": GateTree(row)
            }
            for row_idx, row in enumerate(self.data.iter_rows(named=True))
        ]

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, index: int) -> Dict[str, Any]:
        return self.samples[index]

    def get_tree(self, index: int) -> GateTree:
        return self.samples[index].tree

    def get_metadata(self, index: int) -> Dict[str, Any]:
        return self.samples[index].metadata

    #TODO: Allow for exact values like {"Strain": "B6"}, currently it has to be e.g., {"Strain": "== B6"}
    def find_samples(self, criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]) -> List[int]:
        """
        Find sample indices matching specified metadata criteria.

        Parameters:
            criteria (Dict[str, Union[Any, str, range, Callable[[Any], bool]]]):
                Metadata fields and either:
                    - exact values (e.g., {"Strain": "B6"})
                    - string expressions (e.g., {"Age": "> 10"})
                    - range (e.g., {"Age": range(10, 20)})
                    - callable predicates (e.g., {"Age": lambda x: x > 10})

        Returns:
            List[int]: Indices of samples where metadata matches all criteria.
        """
        matched_samples = []
        for sample in self.samples:
            metadata = sample.metadata
            match = True
            for key, condition in criteria.items():
                value = metadata.get(key)

                if isinstance(condition, str):
                    try:
                        condition = parse_string_condition(condition)
                    except Exception:
                        match = False
                        break
                elif isinstance(condition, range):
                    condition = from_range(condition)

                if callable(condition):
                    if not condition(value):
                        match = False
                        break
                elif value != condition:
                    match = False
                    break

            if match:
                matched_samples.append(sample)
        return matched_samples

    #TODO: If sample_criteria isn't specified, get all samples (currently returns [])
    def get_nodes(self, terms: List[str], sample_criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]= None, exclude_children: bool=True, sample_idx: int=None) -> List[GateNode]:
        """
        Find all nodes across all samples that match the given terms.

        Parameters:
            terms (List[str]): Each string may be a full path or a partial sub-path component.
            sample_idx (int, optional): Index of the sample to search in. If None, search all samples.
                Defaults to None.

        Returns:
            List[GateNode]: All matching nodes across desired trees.
        """
        matched = []
        samples = self.samples
        if sample_idx:
            return self.samples[sample_idx].get_nodes(terms, exclude_children=exclude_children)
        if sample_criteria:
            samples = self.find_samples(sample_criteria)
            for sample in samples:
                matched.append({
                        "metadata": sample.metadata,
                        "nodes": sample.get_nodes(terms, exclude_children=exclude_children)
                    })
            return matched
        for sample in samples:
            matched.append({
                        "metadata": sample.metadata,
                        "nodes": sample.get_nodes(terms, exclude_children=exclude_children)
                    })
        return matched

    #TODO: Make more accurate to tree structure instead of the original DataFrame
    def to_polars(self) -> pl.DataFrame:
        """
        Convert the original DataFrame to a Polars DataFrame.

        Returns:
            pl.DataFrame: The original DataFrame.
        """
        return self.original_df

    def __repr__(self):
        return f"CytoGateParser(num_samples={len(self)})"
    
    def _find_samples_index(self, criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]) -> List[int]:
        """
        Find sample indices matching specified metadata criteria.

        Parameters:
            criteria (Dict[str, Union[Any, str, range, Callable[[Any], bool]]]):
                Metadata fields and either:
                    - exact values (e.g., {"Strain": "B6"})
                    - string expressions (e.g., {"Age": "> 10"})
                    - range (e.g., {"Age": range(10, 20)})
                    - callable predicates (e.g., {"Age": lambda x: x > 10})

        Returns:
            List[int]: Indices of samples where metadata matches all criteria.
        """
        matched_indices = []
        for sample in self.samples:
            metadata = sample["metadata"]
            match = True
            for key, condition in criteria.items():
                value = metadata.get(key)

                if isinstance(condition, str):
                    try:
                        condition = parse_string_condition(condition)
                    except Exception:
                        match = False
                        break
                elif isinstance(condition, range):
                    condition = from_range(condition)

                if callable(condition):
                    if not condition(value):
                        match = False
                        break
                elif value != condition:
                    match = False
                    break

            if match:
                matched_indices.append(idx)
        return matched_indices
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, CytoGateParser):
            return NotImplemented
        return (
            self._as_serializable() == other._as_serializable()
        )

    def _as_serializable(self) -> dict:
        """Return a fully comparable representation of the object."""
        def sanitize(obj):
            if isinstance(obj, (datetime, date)):
                return obj.isoformat()
            elif isinstance(obj, dict):
                return {k: sanitize(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [sanitize(v) for v in obj]
            elif hasattr(obj, "to_dict"):
                return sanitize(obj.to_dict())
            else:
                return obj

        return {
            "samples": [
                {
                    "metadata": sanitize(sample["metadata"]),
                    "tree": sanitize(sample["tree"])
                } for sample in self.samples
            ],
            "metadata_cols": self.metadata_cols
        }