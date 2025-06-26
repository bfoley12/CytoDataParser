import polars as pl
from cytodataparser.structures import GateTree

def test_tree_construction():
    row = pl.DataFrame({
        "Cells | Count": [1000, 10],
        "Cells/CD4 | Count": [600, 1],
        "Cells/CD4/Living | Count": [600, 1],
        "Cells/CD8 | Count": [400, 20],
    }).row(0, named=True)
    tree = GateTree(row)
    assert tree.root.name == "Cells"
    assert "Cells/CD4" in tree.nodes
    assert tree.get_node("Cells/CD8").measures["Count"] == 400
    assert tree.get_node("Cells/CD8").measures["pct_parent"] == 400 / 1000 * 100
    assert tree.max_depth() == 2
    assert tree.get_node("Cells/CD4/Ungated").measures["Count"] == 100
    