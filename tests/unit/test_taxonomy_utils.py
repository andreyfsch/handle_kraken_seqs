from bigtree import Node
from dataset import taxonomy_utils


def test_extract_max_subseqs_set_simple():
    root = Node("root", rank="superkingdom")
    leaf = Node("leaf", rank="sequence", seq="ACGTACGT")
    root.children = [leaf]
    result = taxonomy_utils.extract_max_subseqs_set(root, 4)
    assert result == {"ACGT", "CGTA", "GTAC", "TACG"}


def test_extract_max_subseqs_set_short_sequence():
    root = Node("root", rank="superkingdom")
    leaf = Node("short", rank="sequence", seq="ACG")
    root.children = [leaf]
    result = taxonomy_utils.extract_max_subseqs_set(root, 4)
    assert result == set()


def test_generate_seqs_by_taxon_tree(tmp_path, monkeypatch):
    # Patch paths to mock test data
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    tree = taxonomy_utils.generate_seqs_by_taxon_tree()
    # Should find at least one sequence node
    leaf_nodes = [n for n in tree.descendants if getattr(
        n, "rank", "") == "sequence"]
    assert len(leaf_nodes) > 0


def test_generate_seqs_by_taxon_tree_multisegment(monkeypatch):
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    tree = taxonomy_utils.generate_seqs_by_taxon_tree()
    segment_names = [
        "NC_024503.1", "NC_024504.1", "NC_024505.1", "NC_024506.1",
        "NC_024507.1", "NC_024508.1", "NC_024509.1", "NC_024510.1",
        "NC_024499.1", "NC_024500.1"
    ]
    found = set()
    for seg in segment_names:
        if any(
            seg in getattr(node, "node_name", "")
            for node in tree.descendants
        ):
            found.add(seg)
    assert len(found) == len(
        segment_names), f"Missing segments: {set(segment_names) - found}"


def test_get_subseqs_from_final_node_handles_empty():
    leaf = Node("empty", rank="sequence", seq="")
    result = taxonomy_utils.get_subseqs_from_final_node(leaf, n=3)
    assert isinstance(result, list)
    assert len(result) == 0
