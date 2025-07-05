from dataclasses import dataclass
from collections.abc import Mapping
from typing import Dict, Any, List, Iterator
from cytodataparser.structures.node import GateNode


@dataclass
class NodeResult(Mapping):
    metadata: Dict[str, Any]
    nodes: List[GateNode]

    def __getitem__(self, key: str):
        return {"metadata": self.metadata, "nodes": self.nodes}[key]

    def __iter__(self) -> Iterator[str]:
        return iter(["metadata", "nodes"])

    def __len__(self) -> int:
        return 2