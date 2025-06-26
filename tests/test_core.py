import polars as pl
from cytodataparser import CytoGateParser

def test_parser_metadata_split_and_find():
    df = pl.DataFrame({
        "MouseID": ["A", "B"],
        "Strain": ["B6", "F1"],
        "Cells | Count": [1000, 900],
        "Cells/CD4 | Count": [500, 300],
    })
    parser = CytoGateParser(df)
    assert len(parser) == 2
    assert parser.get_metadata(0)["Strain"] == "B6"
    match = parser.find_samples({"Strain": "== F1"})
    assert match == [1]
