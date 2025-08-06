import pandas as pd
import os
import shutil
import pytest

ROOT = "/home/andrey/generate_kraken_dataset/handle_kraken_seqs/"


def get_small_tsv(original_path, tmp_path, n=5):
    """
    Copy the first n rows from original_path (a Kraken2 TSV) to tmp_path
    as a new file.
    Returns the new file's path.
    """
    df = pd.read_csv(original_path, sep="\t", nrows=n)
    small_path = os.path.join(tmp_path, "small_library_report.tsv")
    df.to_csv(small_path, sep="\t", index=False)
    return small_path


@pytest.fixture(autouse=True)
def clean_tmp_path(tmp_path):
    # Clean up tmp_path before each test to avoid stale files
    yield
    shutil.rmtree(tmp_path, ignore_errors=True)
