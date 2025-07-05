from __future__ import annotations
from typing import Optional, Dict, List, Any, Union
from .tree import GateTree
from .node import GateNode

# TODO: Define __getitem__
class Sample:
    def __init__(self, metadata: Dict[str, Any], tree: GateTree):
        self.metadata = metadata
        self.tree = tree

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