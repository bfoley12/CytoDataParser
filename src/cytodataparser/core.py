from __future__ import annotations
from typing import Optional, List, Dict, Any, Callable, Union
import polars as pl
from .structures import GateTree, GateNode
from .utils.predicates import parse_string_condition, from_range
from datetime import date, datetime

#TODO: Allow for loading with a .wsp file instead of polars dataframe
class CytoGateParser:
    """
    Main class for managing cytometry gating data.
    Parses a Polars DataFrame into metadata and gating trees, one per sample.
    """

    def __init__(self, df: pl.DataFrame, metadata_cols: Optional[List[str]] = None):
        """
        Initialize the CytoGateParser.

        Parameters:
            df (pl.DataFrame): Full dataset containing both metadata and gating results.
            metadata_cols (Optional[List[str]]): List of columns to treat as metadata.
                                             If None, inferred by excluding columns with '|' in name.
        """
        self.original_df = df
        self.metadata_cols = metadata_cols or self._infer_metadata_cols(df)
        self.metadata = df.select(self.metadata_cols)
        self.data = df.drop(self.metadata_cols)
        self.samples: List[Dict[str, Any]] = self._build_samples()

    @classmethod
    def from_samples(cls, samples: List[Dict[str, Any]], metadata_cols: Optional[List[str]] = None) -> "CytoGateParser":
        """
        Construct a CytoGateParser from a serialized list of samples.
        """
        import polars as pl

        metadata_rows = []
        data_rows = []
        for sample in samples:
            metadata_rows.append(sample["metadata"])
            flat_data = cls._flatten_tree_to_row(sample["tree"])
            
            data_rows.append(flat_data)

        df = pl.DataFrame([dict(**m, **d) for m, d in zip(metadata_rows, data_rows)])
        return cls(df, metadata_cols=metadata_cols)

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
        return self.samples[index]["tree"]

    def get_metadata(self, index: int) -> Dict[str, Any]:
        return self.samples[index]["metadata"]

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
        matched_indices = []
        for idx, sample in enumerate(self.samples):
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
        return [self[i] for i in matched_indices]

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
        if sample_idx:
            return self.get_tree(sample_idx).get_nodes(terms, exclude_children=exclude_children)
        if sample_criteria:
            sample_indices = self._find_samples_index(sample_criteria)
            for idx in sample_indices:
                matched.append({
                        "metadata": self.samples[idx]["metadata"],
                        "nodes": self.samples[idx]["tree"].get_nodes(terms, exclude_children=exclude_children)
                    })
            return matched
        for sample in self.samples:
            matched.append({
                        "metadata": sample["metadata"],
                        "nodes": sample["tree"].get_nodes(terms, exclude_children=exclude_children)
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
        for idx, sample in enumerate(self.samples):
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