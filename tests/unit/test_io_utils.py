from dataset import io_utils
import pandas as pd


def test_write_csvs_creates_output(tmp_path):
    # Patch output path to temp directory
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    # Patch taxonomy_utils KRAKEN_PATH
    from dataset import taxonomy_utils
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""

    # Run with test params (small for speed)
    io_utils.write_csvs(
        min_subseq_len=10,
        max_subseq_len=20,
        min_num_seqs=5,
        percentage=100.0,
        compression=False,
        parallel=False
    )
    # Check output files exist and are not empty
    found = False
    for fname in tmp_path.iterdir():
        if fname.suffix == ".csv":
            found = True
            df = pd.read_csv(fname)
            assert not df.empty
    assert found


def test_write_csvs_label_mapping(tmp_path):
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    from dataset import taxonomy_utils
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    io_utils.write_csvs(min_subseq_len=10, max_subseq_len=20, min_num_seqs=5,
                        percentage=100.0, compression=False, parallel=False)
    mapping_files = [f for f in tmp_path.iterdir() if "taxon_labels" in f.name]
    assert mapping_files, "No taxon_labels CSV found"
    for f in mapping_files:
        df = pd.read_csv(f)
        assert not df.empty


def test_write_csvs_parallel_vs_serial(tmp_path):
    io_utils.CSV_OUTPUT_PATH = str(tmp_path)
    from dataset import taxonomy_utils
    taxonomy_utils.KRAKEN_PATH = "tests/data"
    taxonomy_utils.KRAKEN_DATABASE = ""
    io_utils.write_csvs(min_subseq_len=10, max_subseq_len=20, min_num_seqs=5,
                        percentage=100.0, compression=False, parallel=True)
    files_parallel = {f.name for f in tmp_path.iterdir()}
    # Clear tmp_path
    for f in tmp_path.iterdir():
        f.unlink()
    io_utils.write_csvs(min_subseq_len=10, max_subseq_len=20, min_num_seqs=5,
                        percentage=100.0, compression=False, parallel=False)
    files_serial = {f.name for f in tmp_path.iterdir()}
    assert files_parallel == files_serial
