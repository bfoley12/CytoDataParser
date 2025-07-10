from __future__ import annotations
from typing import Optional, List, Dict, Any, Callable, Union, Iterator
import polars as pl
from .structures import GateTree, Sample
from .utils.predicates import parse_string_condition, from_range
from datetime import date, datetime
from cytodataparser.structures import NodeResult
from cytodataparser.io import load_file
from pathlib import Path
#import flowkit as fk

class CytoGateParser:
    """
    Main class for managing cytometry gating data.
    Parses a Polars DataFrame into metadata and gating trees, one per sample.
    """

    def __init__(self, samples: List[Sample]):
        """
        Initialize the CytoGateParser.

        Parameters:
            df (pl.DataFrame): Full dataset containing both metadata and gating results.
            metadata_cols (Optional[List[str]]): List of columns to treat as metadata.
                                             If None, inferred by excluding columns with '|' in name.
        """
        self.samples: List[Sample] = samples

    # TODO: Implement loading from wsp
    '''
    @classmethod
    def from_wsp(cls, path: str) -> CytoGateParser:
        
        Construct a CytoGateParser from a FlowJo WSP file

        Parameters:
            path (str): The path (relative or absolute) to the wsp file (including filename)

        Returns:
            CytoGateParser: A CytoGateParser Object
        
        ws = fk.parse_wsp(path)
    '''

    @classmethod
    def from_file(cls, path: str, sheet_name: Optional[Union[str, None]]=None, fcs_files: Union[str, List[str], Path, List[Path]]="./Unmixed") -> CytoGateParser:
        """
        Construct a CytoGateParser from a file

        Parameters:
            path (str): the path to the file containing the sample information (one of xlsx, xls, csv, or json)
            sheetname (str, optional): sheet_name, if loading from specific xlsx sheet
        """
        return cls(load_file(path, sheet_name, fcs_files))
    
    @classmethod
    def from_samples(cls, samples: List[Sample]) -> CytoGateParser:
        """
        Construct a CytoGateParser from a list of samples.
        """

        return cls(samples)

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
    
    def add_samples(self, path: str, sheet_name: Optional[str] = None) -> CytoGateParser:
        """
        Adds samples to the current CytoGateParser instance.

        Parameters:
            path (str): the path to the file containing the sample information (one of xlsx, xls, csv, or json)
            sheetname (str, optional): sheet_name, if loading from specific xlsx sheet
        """
        self.samples = self.samples + load_file(path, sheet_name)

        return self
    
    def _gather_metadata(self, samples: List[Sample]):
        all_metadata_cols = set()
        for sample in samples:
            all_metadata_cols.update(sample.metadata.keys())

        return list(all_metadata_cols)

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, index: int) -> Sample:
        return self.samples[index]

    def get_tree(self, index: int) -> GateTree:
        return self.samples[index].tree

    def get_metadata(self, index: int) -> Dict[str, Any]:
        return self.samples[index].metadata
    
    #TODO: Allow for exact values like {"Strain": "B6"}, currently it has to be e.g., {"Strain": "== B6"}
    def filter(self, criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]) -> CytoGateParser:
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
            CytoGateParser: A new CytoGateParser with filters applied.
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
        return CytoGateParser.from_samples(matched_samples)

    #TODO: Allow for exact values like {"Strain": "B6"}, currently it has to be e.g., {"Strain": "== B6"}
    # DEPRECATED
    def find_samples(self, criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]) -> List[Sample]:
        """
        Deprecated: Use filter instead
        Find sample indices matching specified metadata criteria.

        Parameters:
            criteria (Dict[str, Union[Any, str, range, Callable[[Any], bool]]]):
                Metadata fields and either:
                    - exact values (e.g., {"Strain": "B6"})
                    - string expressions (e.g., {"Age": "> 10"})
                    - range (e.g., {"Age": range(10, 20)})
                    - callable predicates (e.g., {"Age": lambda x: x > 10})

        Returns:
            List[Sample]: List of samples meeting criteria.
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
    def get_nodes(self, 
                  terms: Optional[Union[List[List[str]], List[str], None]]=None, 
                  sample_criteria: Optional[Dict[str, Union[Any, str, range, Callable[[Any], bool]]]]= None, 
                  exclude_children: bool=True, 
                  sample_idx: Optional[int]=None
                  ) -> List[NodeResult]:
        """
        Find all nodes across all samples that match the given terms.

        Parameters:
            terms (Union[List[List[str]], List[str]], optional): Each string may be a full path or a partial sub-path component. If None, get all nodes for requested samples
            sample_criteria (Dict[str, Union[Any, str, range, Callable[[Any], bool]]], optional): Metadata describing samples to query across.
                Defaults to None (query all samples).
            sample_index (int, optional): Allows specifying index to query sample directly

        Returns:
            List[NodeResult]: All matching nodes across desired samples, with metadata.
                Always of form {"metadata": Dict[str, Any], "nodes": List[GateNode]}
        """
        matched = []
        samples = self.samples
        if sample_idx:
            return [NodeResult(
                metadata=self.samples[sample_idx].metadata,
                nodes=self.samples[sample_idx].get_nodes(terms, exclude_children=exclude_children)
            )]
        if sample_criteria:
            samples = self.find_samples(sample_criteria)
            for sample in samples:
                matched.append(NodeResult(
                        metadata = sample.metadata,
                        nodes = sample.get_nodes(terms, exclude_children=exclude_children)
                ))
            return matched
        for sample in samples:
            matched.append(NodeResult(
                        metadata = sample.metadata,
                        nodes = sample.get_nodes(terms, exclude_children=exclude_children)
            ))
        return matched

    #TODO: Implement
    '''
    def to_polars(self) -> pl.DataFrame:
        """
        Convert the original DataFrame to a Polars DataFrame.

        Returns:
            pl.DataFrame: The original DataFrame.
        """
        return self.original_df
    '''

    def __repr__(self):
        return f"CytoGateParser(num_samples={len(self)})"
    
    def _find_samples_index(self, criteria: Dict[str, Union[Any, str, range, Callable[[Any], bool]]]) -> List[int]:
        """
        Deprecated in favor of find_samples()

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
        for sample, idx in zip(self.samples, range(len(self.samples))):
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
            ]
        }