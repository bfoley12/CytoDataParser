from __future__ import annotations
from typing import Optional, Dict, List, Any
import polars as pl
from .node import GateNode

# TODO: Make sure from_dict and constructor still work
class GateTree:
    """
    Represents a tree of gates for a single sample.
    Automatically builds its structure from a single row of cytometry data.
    """

    # TODO: Consider changing row to be a dict
    def __init__(self, row: Dict[str, Any]):
        """
        Initialize the GateTree using a row of cytometry data.

        Parameters:
            row (pl.Series): A single sample's data with columns in the format 'path/to/gate | Measurement'.
        """
        self.nodes: Dict[str, GateNode] = {}
        self.root: Optional[GateNode] = None  # Defer actual root assignment
        self._build_tree_from_row(row)
        self._build_tree_from_row(row)
        self._insert_ungated_nodes()
        self._compute_percentages()

    def _build_tree_from_row(self, row: Dict[str, Any]):
        """Build the tree structure from a single sample row."""
        for col, value in row.items():
            if "|" not in col:
                continue  # Skip metadata or malformed columns
            path_part, measure_name = [x.strip() for x in col.split("|", 1)]
            node = self._ensure_node_exists(path_part)
            node.measures[measure_name] = value
            #if "Count" in node.measures.keys():
            #    node.update_pct()

    def _ensure_node_exists(self, path: str) -> GateNode:
        parts = path.split("/")
        for i in range(1, len(parts) + 1):
            sub_path = "/".join(parts[:i])
            if sub_path not in self.nodes:
                parent_path = "/".join(parts[:i - 1]) if i > 1 else ""
                parent_node = self.nodes.get(parent_path)
                node = GateNode(name=sub_path, parent=parent_node)
                self.nodes[sub_path] = node

                if parent_node:
                    parent_node.add_child(node)
                else:
                    if self.root is not None:
                        raise ValueError(f"Multiple root nodes found. Current root: {self.root.name}, new root: {sub_path}")
                    self.root = node
        return self.nodes[path]


    def _compute_percentages(self):
        """
        Compute and store percent of parent and root for all nodes with a 'Count' measurement.
        """
        for node in self.nodes.values():
            count = node.measures.get("Count")
            if count is None:
                continue

            if node.parent and "Count" in node.parent.measures:
                parent_count = node.parent.measures["Count"]
                if parent_count > 0:
                    node.measures["pct_parent"] = count / parent_count * 100

            if self.root and "Count" in self.root.measures:
                root_count = self.root.measures["Count"]
                if root_count > 0:
                    node.measures["pct_root"] = count / root_count * 100

    def _insert_ungated_nodes(self):
        """
        For each node, if the sum of its children counts is less than its own count,
        create a new 'ungated' node to account for the difference.
        """
        for node in list(self.nodes.values()):
            if "Count" not in node.measures or not node.children:
                continue

            node_count = node.measures["Count"]
            total_child_count = sum(
                child.measures.get("Count", 0) for child in node.children
            )

            if node_count > total_child_count:
                ungated_count = node_count - total_child_count
                ungated_path = f"{node.name}/Ungated"
                if ungated_path not in self.nodes:
                    ungated_node = GateNode(name=ungated_path, parent=node)
                    ungated_node.measures["Count"] = ungated_count
                    node.add_child(ungated_node)
                    self.nodes[ungated_path] = ungated_node

    def get_nodes(self, terms: List[str], exclude_children: bool = True) -> List[GateNode]:
        """
        Retrieve nodes by either full path or path component match.
        If a single term is given, perform full path or partial match.
        If multiple terms are given, return nodes that contain all parts in order.
        If exclude_children is True, only return the deepest matched nodes whose name ends in the last term.

        Parameters:
            terms (List[str]): One or more full paths or sub-path segments.
            exclude_children (bool): Whether to exclude children of matched nodes (default True).

        Returns:
            List[GateNode]: Matching nodes.
        """
        if isinstance(terms, str):
            terms = [terms.strip()]
        matches = set()
        if len(terms) == 1 and terms[0] in self.nodes:
            matches.add(self.nodes[terms[0]])
        else:
            for path, node in self.nodes.items():
                segments = path.split("/")
                match_idx = 0
                for seg in segments:
                    if match_idx < len(terms) and seg == terms[match_idx]:
                        match_idx += 1
                    if match_idx == len(terms):
                        if exclude_children:
                            if segments[-1] == terms[-1]:
                                matches.add(node)
                        else:
                            matches.add(node)
                        break
        return list(matches)

    def max_depth(self) -> int:
        """Return the maximum depth in the tree."""
        return max((node.depth() for node in self.nodes.values()), default=0)

    def num_leaves(self) -> int:
        """Return the number of leaf (terminal) nodes in the tree."""
        return sum(1 for node in self.nodes.values() if node.is_leaf())

    def avg_branching_factor(self) -> float:
        """Return the average number of children across all non-leaf nodes."""
        non_leaf_nodes = [node for node in self.nodes.values() if not node.is_leaf()]
        if not non_leaf_nodes:
            return 0.0
        return sum(node.degree() for node in non_leaf_nodes) / len(non_leaf_nodes)

    def terminal_gate_density(self) -> float:
        """
        Return the ratio of leaf nodes to total nodes.
        Indicates how shallow or deep the gating tree is overall.
        """
        if not self.nodes:
            return 0.0
        return self.num_leaves() / len(self.nodes)

    def summary(self) -> Dict[str, float | int]:
        """
        Return a summary of overall tree structure metrics.

        Returns:
            Dict[str, float | int]: Tree-level metrics.
        """
        return {
            "num_nodes": len(self.nodes),
            "num_leaves": self.num_leaves(),
            "max_depth": self.max_depth(),
            "avg_branching_factor": self.avg_branching_factor(),
            "terminal_gate_density": self.terminal_gate_density(),
        }
    
    def to_dict(self) -> dict:
        if self.root:
            return self.root.to_dict()
        return {"Error": "root not set"}
    
    @staticmethod
    def from_dict(data: dict) -> "GateTree":
        tree = GateTree.__new__(GateTree)  # bypass __init__
        tree.root = GateNode.from_dict(data)
        tree.nodes = GateTree._make_nodes_dfs(tree.root)
        return tree
    
    @staticmethod
    def _make_nodes_dfs(root: GateNode) -> dict:
        nodes = {root.name: root}
        for child in root.children:
            nodes.update(GateTree._make_nodes_dfs(child))
        return nodes


    def __str__(self):
        return str(self.root)

    def __repr__(self):
        return self.__str__()
