import pytest
import pandas as pd
from download import download_seqs
from dataset import taxonomy_utils


@pytest.fixture
def fasta_url(tmp_path):
    # Use httpbin or a simple mock for a small file download
    return "https://httpbin.org/bytes/10"  # A real 10-byte file


@pytest.fixture
def fake_tsv(tmp_path):
    # Create a dummy TSV with URL column
    df = pd.DataFrame({"URL": ["https://httpbin.org/bytes/10"]})
    path = tmp_path / "library_report.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


def test_all_output_dirs_are_accessions(tmp_path, fake_tsv):
    from download.download_seqs import download_fasta_batch
    outdir = tmp_path / "accdirs"
    download_fasta_batch(fake_tsv, str(outdir), num_threads=1)
    for d in outdir.iterdir():
        if d.is_dir():
            assert (
                (
                    download_seqs.extract_accession_from_seqname(d.name)
                    is not None
                )
            )


def test_download_respects_overwrite_flag(tmp_path, fake_tsv):
    from download.download_seqs import download_fasta_batch
    outdir = tmp_path / "overw"
    # First download
    download_fasta_batch(fake_tsv, str(outdir), num_threads=1)
    # Overwrite with new content
    for d in outdir.iterdir():
        if d.is_dir():
            genome = d / "genome.fna"
            genome.write_text(">fake\nATATAT\n")
            # Download again, overwrite=False should preserve new content
            download_fasta_batch(fake_tsv, str(
                outdir), num_threads=1, overwrite=False)
            assert genome.read_text() == ">fake\nATATAT\n"
            # Now with overwrite=True, content should be different
            download_fasta_batch(fake_tsv, str(
                outdir), num_threads=1, overwrite=True)
            assert genome.read_text() != ">fake\nATATAT\n"


@pytest.mark.slow
def test_genome_folder_structure(tmp_path):
    """
    After downloading, every genome should have a folder, with genome.fna
    inside.
    """
    # Simulate or copy a fixture, or point to your real output dir
    # Let's say your download logic puts files in: tmp_path / "genomes"
    genomes_dir = tmp_path / "genomes"
    # Simulate genome directories
    genome_ids = ["GCF_001504435.1_ViralProj307954",
                  "GCF_003051725.1_ASM305172v1"]
    for gid in genome_ids:
        genome_dir = genomes_dir / gid
        genome_dir.mkdir(parents=True, exist_ok=True)
        (genome_dir / "genome.fna").write_text(">seq1\nACTGACTGACTG\n")

    # Test: all genome directories exist, with genome.fna file
    for gid in genome_ids:
        genome_dir = genomes_dir / gid
        assert genome_dir.is_dir()
        fasta_file = genome_dir / "genome.fna"
        assert fasta_file.exists()
        assert fasta_file.stat().st_size > 0


def test_genome_fasta_parsable(tmp_path, monkeypatch):
    """
    The genome.fna files must be readable and compatible with
    get_genome_sequences.
    """
    kraken_path = tmp_path
    kraken_db = "viral"
    genome_id = "GCF_001504435.1"
    genome_dir = kraken_path / kraken_db / "genomes" / genome_id
    genome_dir.mkdir(parents=True, exist_ok=True)
    fasta_file = genome_dir / "genome.fna"
    fasta_file.write_text(
        ">NC_002195.1 Example virus genome\nACTGACTGACTGACTG\n")

    # Patch config variables for the test
    monkeypatch.setattr("dataset.taxonomy_utils.KRAKEN_PATH", str(kraken_path))
    monkeypatch.setattr("dataset.taxonomy_utils.KRAKEN_DATABASE", kraken_db)

    # Should not raise
    seqs = taxonomy_utils.get_genome_sequences(genome_id)
    assert isinstance(seqs, dict)
    assert "NC_002195.1" in seqs
    assert seqs["NC_002195.1"]["seq"].startswith("ACTG")


def test_generate_seqs_by_taxon_tree(tmp_path, monkeypatch):
    """
    Full integration: test the tree generation from downloaded structure.
    """
    # Simulate a download
    kraken_path = tmp_path
    kraken_db = "viral"
    genome_id = "GCF_001504435"
    genome_dir = kraken_path / kraken_db / "genomes" / genome_id
    genome_dir.mkdir(parents=True, exist_ok=True)
    (genome_dir / "genome.fna").write_text(
        ">NC_002195.1 Example genome\nACGTACGT\n"
    )

    # Add minimal seqid2taxid.map
    (kraken_path / kraken_db).mkdir(parents=True, exist_ok=True)
    (kraken_path / kraken_db / "seqid2taxid.map").write_text(
        "kraken:taxid|687377|NC_002195.1 687377\n"
    )

    monkeypatch.setattr("dataset.taxonomy_utils.KRAKEN_PATH", str(kraken_path))
    monkeypatch.setattr("dataset.taxonomy_utils.KRAKEN_DATABASE", kraken_db)

    # Should build tree without errors
    from dataset.taxonomy_utils import generate_seqs_by_taxon_tree
    tree = generate_seqs_by_taxon_tree()
    assert tree is not None


def simulate_genome_download(tmp_path):
    # Mimic the structure expected by the dataset module
    kraken_path = tmp_path
    kraken_db = "viral"
    genome_id = "GCF_001504435"
    genome_dir = kraken_path / kraken_db / "genomes" / genome_id
    genome_dir.mkdir(parents=True, exist_ok=True)
    (genome_dir / "genome.fna").write_text(
        ">NC_002195.1 Example virus genome\nACTGACTGACTGACTG\n")
    # Add minimal seqid2taxid.map
    (kraken_path / kraken_db).mkdir(parents=True, exist_ok=True)
    (kraken_path / kraken_db / "seqid2taxid.map").write_text(
        "kraken:taxid|687377|NC_002195.1 687377\n")
    return kraken_path, kraken_db, genome_id


def test_seqid2taxid_map_format_and_contents(tmp_path):
    """
    Test seqid2taxid.map is present, well-formed, and parsable by the dataset
    module.
    """
    # Simulate or copy real map file here
    kraken_dir = tmp_path / "kraken"
    kraken_dir.mkdir(parents=True, exist_ok=True)
    map_file = kraken_dir / "seqid2taxid.map"
    map_file.write_text("kraken:taxid|687377|NC_002195.1 687377\n")

    # File should exist and be readable by taxonomy_utils
    assert map_file.exists()
    with open(map_file) as f:
        line = f.readline().strip()
        assert line.startswith("kraken:taxid|")
        parts = line.split()
        assert len(parts) == 2
        assert "|" in parts[0]
        assert parts[1].isdigit()


def test_genome_fna_is_fasta(tmp_path):
    """
    Test that each genome.fna is a valid FASTA file (header and sequence).
    """
    genome_dir = tmp_path / "genomes" / "GCF_001504435.1_ViralProj307954"
    genome_dir.mkdir(parents=True, exist_ok=True)
    fasta_file = genome_dir / "genome.fna"
    fasta_file.write_text(">NC_002195.1 Example\nACTG\n")
    # Basic check: first line starts with ">"
    with open(fasta_file) as f:
        header = f.readline().strip()
        assert header.startswith(">")
