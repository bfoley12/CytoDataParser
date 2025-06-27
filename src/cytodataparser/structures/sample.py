from __future__ import annotations
from typing import Optional, Dict, List, Any
from .tree import GateTree

class Sample:
    def __init__(self, metadata: Dict[str, Any], tree: GateTree) -> Sample:
        self.metadata = metadata
        self.tree = tree

    def get_nodes(self, terms, exclude_children):
        return self.tree.get_nodes(terms=terms, exclude_children=exclude_children)