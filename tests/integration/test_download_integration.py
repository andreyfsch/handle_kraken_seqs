import pytest
import os
import pandas as pd
from download import download_seqs, webscrap
from tests.integration.conftest import get_small_tsv, ROOT


@pytest.mark.slow
def test_seqid2taxid_map_present_and_valid(tmp_path):
    """
    Integration test: Actually downloads and extracts the viral Kraken2 DB
    tarball, then checks that seqid2taxid.map exists within the extracted
    structure.
    """
    df = webscrap.fetch_kraken2_collections()
    viral = df[df["Collection"].str.lower().str.contains("viral")]
    assert not viral.empty, "No viral row found in kraken2 collections table"

    tarball_url = viral.iloc[0]["tarball_url"]
    assert tarball_url and tarball_url.startswith(
        "http"), f"No valid viral tarball_url found: {tarball_url}"

    kraken_db_dir = tmp_path / "kraken_viral"
    ok = download_seqs.download_and_extract_kraken_db_tarball(
        tarball_url, str(kraken_db_dir)
    )
    assert ok, f"Failed to fetch/extract viral Kraken2 tarball: {tarball_url}"

    # Find seqid2taxid.map anywhere in the output
    found = list(kraken_db_dir.rglob("seqid2taxid.map"))
    assert found, "No seqid2taxid.map found in extracted kraken DB"
    print(f"seqid2taxid.map found at: {found[0]}")


@pytest.mark.slow
def test_downloaded_fastas_are_valid(tmp_path):
    tsv_path = os.path.join(ROOT, "tests", "data", "library_report.tsv")
    out_dir = tmp_path / "fasta"
    # <--- FIX: Actually call the downloader!
    count = download_seqs.download_fasta_batch(
        str(tsv_path), str(out_dir), num_threads=2
    )
    assert count > 0, "No files were downloaded!"
    df = pd.read_csv(tsv_path, sep='\t')
    url_count = len(df["URL"])
    if not out_dir.exists():
        pytest.fail(f"Download directory {out_dir} was not created!")
    accession_dirs = [
        d for d in out_dir.iterdir()
        if (
            d.is_dir()
            and (
                download_seqs.extract_accession_from_seqname(d.name)
                is not None
            )
        )
    ]
    assert len(accession_dirs) == url_count
    for d in accession_dirs:
        fna = d / "genome.fna"
        assert fna.exists()
        assert fna.stat().st_size > 0
        with open(fna) as f:
            line = f.readline()
            assert line.startswith(">")


@pytest.mark.slow
def test_download_fasta_batch_real_data(tmp_path):
    original_tsv = os.path.join(
        ROOT, "tests", "data", "library_report.tsv")
    small_tsv = get_small_tsv(original_tsv, tmp_path, n=50)
    out_dir = tmp_path / "fasta"
    count = download_seqs.download_fasta_batch(
        str(small_tsv), str(out_dir), num_threads=2
    )
    assert count > 0, "No files were downloaded!"
    df = pd.read_csv(small_tsv, sep='\t')
    url_count = len(df["URL"])
    if not out_dir.exists():
        pytest.fail(f"Download directory {out_dir} was not created!")
    accession_dirs = [
        d for d in out_dir.iterdir()
        if d.is_dir() and download_seqs.extract_accession_from_seqname(d.name)
    ]
    assert len(accession_dirs) == url_count
    for d in accession_dirs:
        fna = d / "genome.fna"
        assert fna.exists()
        assert fna.stat().st_size > 0
