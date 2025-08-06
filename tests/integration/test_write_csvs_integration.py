import pytest
import os
from dataset import io_utils
from download import download_seqs, webscrap
import shutil
from pathlib import Path
from tests.integration.conftest import get_small_tsv, ROOT


def restrict_seqid2taxid_map(kraken_db_dir, genomes_dir):
    """
    Keep only the entries in seqid2taxid.map for which a genome exists in
    genomes_dir.
    """
    seqid2taxid = Path(kraken_db_dir) / "seqid2taxid.map"
    if not seqid2taxid.exists():
        raise RuntimeError(f"seqid2taxid.map not found at {seqid2taxid}")

    # List genome accessions actually downloaded
    available = {d.name for d in Path(genomes_dir).iterdir() if d.is_dir()}
    # Read and filter map
    rows = []
    with open(seqid2taxid) as f:
        for line in f:
            acc = line.split("\t")[0].strip()
            # Extract the accession after the last "|"
            if "|" in acc:
                acc_name = acc.split("|")[-1]
            else:
                acc_name = acc
            if acc_name in available:
                rows.append(line)
    # Overwrite with filtered map
    with open(seqid2taxid, "w") as f:
        f.writelines(rows)


def ensure_viral_kraken_db(tmp_path):
    # Fetch the viral tarball URL from webscrap.py
    df = webscrap.fetch_kraken2_collections()
    viral_row = df[df['Collection'].str.lower().str.contains("viral")].iloc[0]
    tarball_url = viral_row['tarball_url']
    # Download tarball to a file, then extract
    kraken_db_dir = tmp_path / "kraken_viral"
    kraken_db_dir.mkdir(parents=True, exist_ok=True)
    tarball_file = kraken_db_dir / "viral.tar.gz"
    # Download the tarball
    import requests
    with requests.get(tarball_url, stream=True) as r:
        r.raise_for_status()
        with open(tarball_file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    # Extract tarball
    import tarfile
    with tarfile.open(tarball_file, "r:gz") as tar:
        tar.extractall(path=kraken_db_dir)
    assert (
        kraken_db_dir / "seqid2taxid.map"
    ).exists(), "seqid2taxid.map missing"
    return kraken_db_dir


@pytest.mark.slow
def test_write_csvs_runs_with_50_real_downloaded(tmp_path, monkeypatch):
    """
    Integration test: use top 50 rows of real TSV, download actual FASTAs,
    filter seqid2taxid.map, and verify write_csvs works and produces output.
    """
    original_tsv = os.path.join(ROOT, "tests", "data", "library_report.tsv")
    small_tsv = get_small_tsv(original_tsv, tmp_path, n=50)
    out_dir = tmp_path / "fasta"
    count = download_seqs.download_fasta_batch(
        str(small_tsv), str(out_dir), num_threads=2)
    assert count > 0, "No files were downloaded!"

    kraken_db_dir = ensure_viral_kraken_db(tmp_path)

    genomes_dir = out_dir
    restrict_seqid2taxid_map(kraken_db_dir, genomes_dir)

    target_genomes_dir = kraken_db_dir / "genomes"
    target_genomes_dir.mkdir(parents=True, exist_ok=True)
    for d in genomes_dir.iterdir():
        if d.is_dir():
            shutil.move(str(d), str(target_genomes_dir / d.name))

    downloaded_accessions = set(d.name for d in (
        kraken_db_dir / "genomes").iterdir() if d.is_dir())
    map_path = kraken_db_dir / "seqid2taxid.map"
    with open(map_path) as f:
        mapped_accessions = set(
            line.split('\t')[0]
            for line in f
            if line.strip()
        )
    print("Downloaded accessions:", downloaded_accessions)
    print("Mapped accessions:", mapped_accessions)
    print("Missing in map:", downloaded_accessions - mapped_accessions)

    monkeypatch.setattr("dataset.taxonomy_utils.KRAKEN_PATH", str(tmp_path))
    monkeypatch.setattr(
        "dataset.taxonomy_utils.KRAKEN_DATABASE", "kraken_viral")
    monkeypatch.setattr("dataset.io_utils.CSV_OUTPUT_PATH",
                        str(tmp_path / "csvs"))

    io_utils.write_csvs(
        min_subseq_len=4,
        max_subseq_len=8,
        min_num_seqs=1,
        percentage=100.0,
        compression=False,
        parallel=False,
        max_workers=2,
        rng=None,
    )
    # Check output CSVs exist
    csvs = list((tmp_path / "csvs").glob("*.csv"))
    assert csvs, "No CSV output produced by write_csvs"
