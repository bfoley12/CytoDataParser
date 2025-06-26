from cytodataparser.structures import GateNode

def test_node_structure():
    root = GateNode("Root")
    child = GateNode("Root/Child", parent=root)
    root.add_child(child)

    assert root.is_root()
    assert not root.is_leaf()
    assert child.is_leaf()
    assert child.depth() == 1
    assert root.height() == 1
    assert child.parent == root
    assert root.children == [child]