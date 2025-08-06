import os
import logging
import requests
import gzip
import shutil
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)


def extract_accession_from_seqname(seqname: str) -> str:
    # Return the first token (before any whitespace)
    return seqname.split(' ')[0]


def extract_accession_from_filename(filename: str) -> str:
    """
    (For test compatibility) Extract accession from typical file names,
    fallback to using seqname logic.
    """
    # Try to extract before first _ if pattern like ACC_XXX,
    # else before first dot, else whole base
    if not filename or not isinstance(filename, str):
        return None
    # E.g., NC_123456.1_SomeLabel_genomic.fna.gz -> NC_123456.1
    base = filename.split("_")[0]
    if '.' in base:
        return base
    return filename.split('.')[0]


def download_and_prepare_fasta(url, accession, output_root, overwrite=False):
    """
    Download gzipped FASTA from url, extract to genome.fna under
    output_root/accession/.
    Returns True if successful.
    """
    target_dir = os.path.join(output_root, accession)
    os.makedirs(target_dir, exist_ok=True)
    fasta_path = os.path.join(target_dir, "genome.fna")
    tmp_gz = os.path.join(target_dir, "tmp_download.fna.gz")

    if os.path.exists(fasta_path) and not overwrite:
        logger.info(f"File {fasta_path} exists, skipping.")
        return True

    try:
        # Download
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(tmp_gz, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        # Decompress
        with gzip.open(tmp_gz, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(tmp_gz)
        return True
    except Exception as e:
        logger.error(f"Failed: {url} ({accession}): {e}")
        try:
            if os.path.exists(tmp_gz):
                os.remove(tmp_gz)
        except Exception:
            pass
        return False


def download_fasta_batch(
    library_report_tsv,
    output_dir,
    num_threads=None,
    overwrite=False
):
    """
    Download all FASTA files listed in a Kraken2 library_report.tsv,
    with a tqdm progress bar over all downloads.

    For each row, use 'Sequence Name' column to extract accession.

    Returns
    -------
    num_downloaded : int
        The number of files successfully downloaded and extracted.
    """
    logger.info(f"Reading library report: {library_report_tsv}")
    df = pd.read_csv(library_report_tsv, sep='\t')
    if 'URL' not in df.columns:
        raise KeyError(
            (
                f"'URL' column not found in "
                f"{library_report_tsv}. Columns are: {df.columns.tolist()}"
            )
        )
    if 'Sequence Name' not in df.columns:
        logger.warning(
            "'Sequence Name' column not found in TSV. Will try to extract "
            "accessions from file names. This is deprecated and should only "
            "be used for old-format TSVs."
        )

    tasks = []
    for i, row in tqdm(
        df.iterrows(),
        total=len(df),
        desc="All FASTA downloads"
    ):
        url = row["URL"]
        if 'Sequence Name' in df.columns:
            accession = extract_accession_from_seqname(row["Sequence Name"])
        else:
            accession = extract_accession_from_filename(os.path.basename(url))
        if not accession or not isinstance(accession, str):
            logger.warning(f"Skipping row {i}: no accession in Sequence Name.")
            continue
        tasks.append((url, accession))

    # ---- PARALLEL download here ----
    successes = 0
    num_threads = num_threads or os.cpu_count() or 4
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(download_and_prepare_fasta, url,
                            accession, output_dir, overwrite)
            for url, accession in tasks
        ]
        for f in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="All FASTA downloads (parallel)"
        ):
            if f.result():
                successes += 1
    return successes


def download_and_extract_kraken_db_tarball(tarball_url_or_path, output_dir):
    """
    Download (if url) or extract (if local file) a Kraken2 tarball and
    extract it.
    Returns True if seqid2taxid.map is found after extraction.
    """
    import tarfile

    os.makedirs(output_dir, exist_ok=True)
    # If it's already a file, just extract
    if os.path.exists(tarball_url_or_path):
        tarball_path = tarball_url_or_path
    else:
        # Download logic as before
        tarball_path = os.path.join(output_dir, "kraken_db.tar.gz")
        with requests.get(tarball_url_or_path, stream=True, timeout=180) as r:
            r.raise_for_status()
            with open(tarball_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
    try:
        with tarfile.open(tarball_path, "r:gz") as tar:
            tar.extractall(path=output_dir)
        # Look for seqid2taxid.map anywhere under output_dir
        for root, dirs, files in os.walk(output_dir):
            if "seqid2taxid.map" in files:
                return True
        logger.error("seqid2taxid.map not found in extracted tarball!")
        return False
    except Exception as e:
        logger.error(f"Failed to extract tarball: {tarball_url_or_path}: {e}")
        return False
