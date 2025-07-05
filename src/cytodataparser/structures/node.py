from __future__ import annotations
from typing import Optional, Dict, List
import polars as pl

# TODO: create __getitem__()
class GateNode:
    """
    Represents a single gate in the cytometry gating hierarchy.
    Stores information about its parent, children, and measured attributes.
    """

    def __init__(
        self,
        name: str = "",
        parent: Optional[GateNode] = None,
        measures: Optional[Dict[str, float]] = {},
    ):
        """
        Initialize a GateNode.

        Parameters:
            name (str): Full path name of the gate (e.g., 'Cells/Singlets/CD3+').
            parent (Optional[GateNode]): Reference to the parent GateNode, or None if root.
        """
        self.name: str = name
        self.parent: Optional[GateNode] = parent
        self.children: List[GateNode] = []
        
        # If measures isn't provided, initialize as an empty dictionary
        self.measures: Dict[str, float] = dict(measures) if measures else {}  # e.g., {"Count": 1000, "MFI": 432.1}

    def add_child(self, child: GateNode):
        """Add a child node to this gate."""
        self.children.append(child)

    def is_root(self) -> bool:
        """Return True if the node has no parent."""
        return self.parent is None

    def is_leaf(self) -> bool:
        """Return True if the node has no children."""
        return len(self.children) == 0

    def degree(self) -> int:
        """Return the number of children of this node."""
        return len(self.children)

    def depth(self) -> int:
        """Return the depth of this node (distance from root)."""
        depth = 0
        node = self
        while node.parent:
            depth += 1
            node = node.parent
        return depth

    def height(self) -> int:
        """Return the height of this node (longest path to a leaf)."""
        if self.is_leaf():
            return 0
        return 1 + max(child.height() for child in self.children)

    def num_descendants(self) -> int:
        """Return the number of descendant nodes under this node."""
        return sum(child.num_descendants() + 1 for child in self.children)

    def is_balanced(self) -> bool:
        """Return True if the heights of all children differ by at most 1."""
        if self.is_leaf():
            return True
        heights = [child.height() for child in self.children]
        return max(heights) - min(heights) <= 1

    def summary(self) -> Dict[str, float | int | bool]:
        """Return a dictionary summarizing key structural metrics of the node."""
        return {
            "depth": self.depth(),
            "height": self.height(),
            "degree": self.degree(),
            "num_descendants": self.num_descendants(),
            "is_leaf": self.is_leaf(),
            "is_balanced": self.is_balanced(),
        }
    
    def flatten(self):
        return {"node": self.name, **self.measures}
    
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "measures": self.measures,
            "children": [child.to_dict() for child in self.children]
        }
    
    def update_pct(self):
        if self.parent and self.parent.measures["Count"] > 0:
            self.measures["pct_parent"] = self.measures["Count"] / self.parent.measures["Count"] if self.parent.measures["Count"] > 0 else 0


    @staticmethod
    def from_dict(data: dict) -> "GateNode":
        node = GateNode(name=data["name"], measures=data["measures"])
        node.children = [GateNode.from_dict(c) for c in data["children"]]
        return node
    
    def __str__(self, level=0):
        prefix = "    " * level + ("└── " if level > 0 else "")
        s = prefix + self.name + "\n"
        for child in self.children:
            s += child.__str__(level + 1)
        return s

    def __repr__(self):
        return self.__str__()

    # TODO: Finish working on this
    def __getitem__(self, key: str) -> float | GateNode | List[GateNode] | int | bool:
        if key in self.measures:
            return self.measures[key]
        
        summary = self.summary()
        if key in summary:
            return summary[key]
        
        for child in self.children:
            if child.name.split("/")[-1] == key:
                return child
        
        raise KeyError(f"'{key}' not found in measures, summary, or children of node '{self.name}'")
        
