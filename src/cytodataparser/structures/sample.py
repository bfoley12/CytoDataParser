from __future__ import annotations
from typing import Optional, Dict, List, Any, Union
from .tree import GateTree
from .node import GateNode
from pathlib import Path

class Sample:
    def __init__(self, metadata: Dict[str, Any], tree: GateTree, fcs_file: Optional[Path]=None):
        self.metadata = metadata
        self.tree = tree
        self.fcs_file = fcs_file

    def get_nodes(self, terms, exclude_children) -> List[GateNode]:
        return self.tree.get_nodes(terms=terms, exclude_children=exclude_children)
    
    def __getitem__(self, key: Union[List[str], str]) -> Any:
        if key in self.metadata:
            return self.metadata[key]
        if isinstance(key, str):
            key = [key]
        result = self.tree.get_nodes(key, exclude_children=False)
        if result:
            return result
        raise KeyError(f"{key} not found in metadata or tree.")