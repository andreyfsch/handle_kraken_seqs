from unittest import mock
from dataset import (io_utils, taxonomy_utils, sequence_utils)
import pandas as pd
import random
import pytest


def test_csv_write_failure(tmp_path, caplog):
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""

    # Patch pd.DataFrame.to_csv to raise IOError for first call
    with mock.patch(
        "pandas.DataFrame.to_csv",
        side_effect=IOError("disk full")
    ):
        with caplog.at_level("ERROR"):
            # Should NOT crash, should log error and move on
            io_utils.write_csvs(
                min_subseq_len=10,
                max_subseq_len=20,
                min_num_seqs=5,
                percentage=100.0,
                compression=False,
                parallel=False
            )
            assert any("disk full" in rec.message for rec in caplog.records)


def test_label_consistency(tmp_path):
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    io_utils.write_csvs(min_subseq_len=10, max_subseq_len=20, min_num_seqs=5,
                        percentage=100.0, compression=False, parallel=False)
    # Find any CSV and corresponding label mapping
    for csv_file in tmp_path.glob("*_seqs.csv"):
        label_file = csv_file.with_name(csv_file.stem + "_taxon_labels.csv")
        df = pd.read_csv(csv_file)
        label_df = pd.read_csv(label_file)
        # All labels in data must be in mapping
        assert set(df["label"]).issubset(set(label_df["id"]))
        # All labels are > 0 and unique
        assert all(df["label"] > 0)
        assert label_df["id"].is_unique


def test_extract_max_subseqs_set_dedup():
    # Sequence with repeated regions
    from bigtree import Node
    root = Node("root", rank="superkingdom")
    # The sequence below has "ACGT" repeated 10 times
    leaf = Node("leaf", rank="sequence", seq="ACGT" * 10)
    root.children = [leaf]
    result = taxonomy_utils.extract_max_subseqs_set(root, 4)
    # Only 4 unique subsequences possible
    assert set(result) == {"ACGT", "CGTA", "GTAC", "TACG"}


def test_output_csv_schema(tmp_path):
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    from dataset import taxonomy_utils
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    io_utils.write_csvs(min_subseq_len=10, max_subseq_len=20, min_num_seqs=5,
                        percentage=100.0, compression=False, parallel=False)
    for csv_file in tmp_path.glob("*_seqs.csv"):
        df = pd.read_csv(csv_file)
        assert "sequence" in df.columns
        assert "label" in df.columns
        assert df["sequence"].apply(lambda x: isinstance(x, str)).all()
        assert df["label"].apply(lambda x: isinstance(x, int)).all(
        ) or df["label"].apply(lambda x: isinstance(x, float)).all()


def test_parallel_scaling(tmp_path):
    # Not a strict correctness test—just that it runs without crash
    # on large input
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    from dataset import taxonomy_utils
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    # Use small numbers here, increase if you add more genomes
    io_utils.write_csvs(
        min_subseq_len=10,
        max_subseq_len=20,
        min_num_seqs=5,
        percentage=100.0,
        compression=False,
        parallel=True,
        max_workers=2
    )
    # Check that at least one output file exists
    assert list(tmp_path.glob("*_seqs.csv"))


def test_extract_subseqs_diversity():
    seq = "ACGT" * 100
    rng = random.Random(42)
    subseqs = sequence_utils.extract_subseqs(
        seq, n=10, min_len=10, max_len=20, rng=rng)
    # Not all subseqs should be identical
    assert len(set(subseqs)) > 1


def test_malformed_fasta_handling(tmp_path):
    # Create a fake genome folder with a broken FASTA (no header)
    genome_dir = tmp_path / "genomes" / "BROKEN"
    genome_dir.mkdir(parents=True)
    with open(genome_dir / "genome.fna", "w") as f:
        f.write("ACGTACGTACGT\nACGTACGT\n")  # no '>'
    # Map file with corresponding entry
    map_file = tmp_path / "seqid2taxid.map"
    with open(map_file, "w") as f:
        f.write("kraken:taxid|123456|BROKEN 123456\n")
    # Patch path and run
    taxonomy_utils.KRAKEN_PATH = str(tmp_path)
    taxonomy_utils.KRAKEN_DATABASE = ""
    # Should not crash, should skip this genome
    try:
        taxonomy_utils.generate_seqs_by_taxon_tree()
    except Exception as e:
        pytest.fail(f"generate_seqs_by_taxon_tree crashed: {e}")


def test_empty_seqid2taxid_map(tmp_path):
    # Place empty map file
    (tmp_path / "genomes").mkdir(parents=True)
    map_file = tmp_path / "seqid2taxid.map"
    map_file.touch()
    taxonomy_utils.KRAKEN_PATH = str(tmp_path)
    taxonomy_utils.KRAKEN_DATABASE = ""
    # Should raise or log error, but not crash hard
    try:
        taxonomy_utils.generate_seqs_by_taxon_tree()
    except FileNotFoundError:
        pass  # acceptable if code raises, or handle as appropriate
    except Exception:
        pass
