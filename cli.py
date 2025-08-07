# cli.py

"""
Kraken2 Viral Genomes Dataset Utility CLI
=========================================

Command-line interface for managing Kraken2 viral genomes datasets.
Allows users to list available Kraken2/Bracken collections, download
associated databases and genome FASTA files, and generate CSV datasets
suitable for machine learning workflows.

This CLI is intended to streamline the process of retrieving, verifying,
and preparing Kraken2 viral genome data for bioinformatics and ML tasks.

Main Features
-------------
- List all available Kraken2/Bracken collections from the AWS reference
  site.
- Download all files (library report, tarball, and checksum) associated
  with a selected collection.
- Download all genome FASTA files listed in a collection's
  `library_report.tsv`.
- Generate labeled CSV datasets from downloaded genomes with customizable
  parameters.

Usage Examples
--------------
List available collections:
    $ python cli.py list-collections

Download all files for the "Viral" collection:
    $ python cli.py download-db --collection Viral --outdir ./kraken_downloads

Download all genomes from a library report:
    $ python cli.py download-genomes \
        --library-report ./kraken_downloads/library_report.tsv \
        --output-dir ./kraken_downloads/genomes

Generate CSVs from downloaded genomes:
    $ python cli.py generate-csv \
        --output-dir ./csvs \
        --min-subseq-len 100 \
        --max-subseq-len 400

Functions
---------
download_file(url, dest, show_progress=True)
    Download a file from a URL to a local destination with optional
    progress bar.

verify_md5(file_path, md5_file)
    Verify the MD5 checksum of a file against a .md5 file.

list_collections(args)
    Print a table of all available Kraken2/Bracken collections.

download_db(args)
    Download the files (library report, tarball, MD5) for a chosen
    collection.

download_genomes(args)
    Download all genome FASTA files listed in a library_report.tsv file.

generate_csv(args)
    Generate labeled CSV datasets from downloaded genomes.

get_parser()
    Construct and return an argument parser for the CLI.

run_cli(args=None)
    Run the CLI with the provided arguments (or from sys.argv).

Notes
-----
- Requires an internet connection for downloading files and listing
  collections.
- CSV generation requires the genome files to be organized in the
  expected structure and a properly set up Kraken2 database.

See Also
--------
download.webscrap.fetch_kraken2_collections : For fetching collection info
download.download_seqs.download_fasta_batch : For downloading FASTA files
dataset.io_utils.write_csvs : For generating CSV datasets

"""

import argparse
import logging
import os
import sys
import hashlib
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse

import requests
from tqdm import tqdm

from download.webscrap import fetch_kraken2_collections
from download.download_seqs import download_fasta_batch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s:%(name)s: %(message)s"
)
logger = logging.getLogger("kraken_cli")


def verify_all_md5(md5_file, root_dir):
    """
    Verifies all files in an md5 file exist under root_dir and match their md5.
    """
    failed = []
    with open(md5_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            expected = parts[0]
            fname = parts[-1].lstrip('*')
            # Locate file relative to root_dir
            fpath = os.path.join(root_dir, fname)
            if not os.path.exists(fpath):
                print(f"Missing file for md5 check: {fpath}")
                failed.append((fname, 'missing'))
                continue
            # Compute md5
            hasher = hashlib.md5()
            with open(fpath, "rb") as file:
                for chunk in iter(lambda: file.read(8192), b""):
                    hasher.update(chunk)
            found = hasher.hexdigest()
            if found != expected:
                print(
                    (
                        f"MD5 mismatch for {fname}:\n"
                        f"  Expected: {expected}\n"
                        f"  Found:    {found}"
                    )
                )
                failed.append((fname, 'mismatch'))
            else:
                print(f"MD5 ok for {fname}")
    return failed


def download_file(url: str, dest: str, show_progress: bool = True) -> str:
    logger.info(f"Starting download: {url}")
    try:
        r = requests.get(url, stream=True, timeout=60)
        r.raise_for_status()
        total = int(r.headers.get('content-length', 0))
        pbar = tqdm(total=total, unit='B', unit_scale=True,
                    desc=os.path.basename(dest)) if show_progress else None
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    if pbar:
                        pbar.update(len(chunk))
        if pbar:
            pbar.close()
        logger.info(f"Download complete: {dest}")
        return dest
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        raise


def get_md5_for_file(md5_file, target_filename):
    """
    Scan an md5 file and extract the hash for a specific filename.
    """
    with open(md5_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts or len(parts) < 2:
                continue
            # Filename is last part
            filename = parts[-1]
            # Remove possible leading '*' (for binary mode md5sum)
            filename = filename.lstrip('*')
            if filename == target_filename:
                return parts[0]
    raise ValueError(f"MD5 for {target_filename} not found in {md5_file}")


def verify_md5(target_file, md5_file):
    """Checks the MD5 hash of target_file against the md5_file content."""
    import hashlib
    # Read the hash value from the md5 file
    hash_str = get_md5_for_file(md5_file, target_file.split('/')[-1])

    # Compute hash of downloaded file
    with open(target_file, "rb") as f:
        file_hash = hashlib.md5()
        for chunk in iter(lambda: f.read(8192), b""):
            file_hash.update(chunk)
    md5sum = file_hash.hexdigest()
    # The bug: sometimes md5_file has a different file name!
    if md5sum != hash_str:
        print("MD5 MISMATCH")
        print(f"Computed: {md5sum} (for {target_file})")
        return False
    return True


def list_collections(args):
    df = fetch_kraken2_collections()
    print(df[["Collection", "Contains", "Date", "Archive size (GB)",
          "Index size (GB)"]].to_string(index=False))


def download_db(args):
    df = fetch_kraken2_collections()
    collection = args.collection
    row = df[df['Collection'].str.lower() == collection.lower()]
    if row.empty:
        logger.error(
            (
                f"Collection '{collection}' not found. "
                "Use list-collections to see available options."
            )
        )
        sys.exit(1)
    row = row.iloc[0]
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    to_download = []
    if row['library_report_url']:
        to_download.append(
            (
                row['library_report_url'],
                os.path.join(outdir, "library_report.tsv")
            )
        )
    if row['tarball_url']:
        tarball_name = os.path.basename(urlparse(row['tarball_url']).path)
        to_download.append(
            (row['tarball_url'], os.path.join(outdir, tarball_name)))
    if row['md5_url']:
        md5_name = os.path.basename(urlparse(row['md5_url']).path)
        to_download.append((row['md5_url'], os.path.join(outdir, md5_name)))

    if not to_download:
        logger.error(
            f"No downloadable files found for collection {collection}.")
        sys.exit(1)

    logger.info(f"Downloading files for collection '{collection}' to {outdir}")
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        list(tqdm(pool.map(lambda args: download_file(*args),
             to_download), total=len(to_download), desc="Files"))

    if row['md5_url'] and row['tarball_url'] and not args.skip_md5:
        tarball_path = os.path.join(outdir, os.path.basename(
            urlparse(row['tarball_url']).path))
        md5_path = os.path.join(outdir, os.path.basename(
            urlparse(row['md5_url']).path))

        # --- Extract tarball using your existing function ---
        from download.download_seqs import (
            download_and_extract_kraken_db_tarball
        )

        logger.info(f"Extracting tarball {tarball_path} ...")
        # We just want to extract a *local* tarball, so let's call the function
        if not os.path.exists(tarball_path):
            logger.error(f"Tarball file not found at {tarball_path}!")
            sys.exit(2)

        # Call your function but patch:
        # it expects a URL but we give it a local file
        ok = download_and_extract_kraken_db_tarball(tarball_path, outdir)
        if not ok:
            logger.error(
                "Failed to extract Kraken tarball for MD5 verification!")
            sys.exit(2)

        # --- Verify ALL files listed in the md5 file ---
        logger.info(
            (
                f"Verifying MD5 for all files listed in {md5_path} "
                f"under {outdir} ..."
            )
        )
        failures = verify_all_md5(md5_path, outdir)
        if failures:
            logger.error(f"MD5 failed for files: {failures}")
            sys.exit(2)
        else:
            logger.info("All MD5 checks passed.")


def download_genomes(args):
    threads = args.threads or os.cpu_count()
    logger.info(
        (
            f"Downloading genomes listed in {args.library_report} "
            f"to {args.output_dir} using {threads} threads."
        )
    )
    download_fasta_batch(
        args.library_report,
        args.output_dir,
        num_threads=threads
    )


def generate_csv(args):
    from dataset import io_utils
    import random

    rng = None
    if args.seed is not None:
        rng = random.Random(args.seed)

    logger.info("Generating CSVs with parameters:")
    for key, val in vars(args).items():
        logger.info(f"  {key}: {val}")

    io_utils.write_csvs(
        min_subseq_len=args.min_subseq_len,
        max_subseq_len=args.max_subseq_len,
        rng=rng,
        compression=args.compression,
        min_num_seqs=args.min_num_seqs,
        percentage=args.percentage,
        parallel=args.parallel,
        max_workers=args.max_workers
    )
    logger.info("CSV generation complete.")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Kraken2 Viral Genomes Dataset CLI")
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Command to run")

    parser_list = subparsers.add_parser(
        "list-collections", help="List available Kraken2/Bracken collections")
    parser_list.set_defaults(func=list_collections)

    parser_dl = subparsers.add_parser(
        "download-db",
        help="Download the files for a specific Kraken2 collection"
    )
    parser_dl.add_argument("--collection", required=True,
                           help="Collection name (see list-collections)")
    parser_dl.add_argument(
        "--outdir",
        default="./kraken_downloads",
        help="Directory to save downloads"
    )
    parser_dl.add_argument("--skip-md5", action="store_true",
                           help="Skip MD5 check after download")
    parser_dl.set_defaults(func=download_db)

    parser_genomes = subparsers.add_parser(
        "download-genomes",
        help="Download all FASTA sequences in library_report.tsv"
    )
    parser_genomes.add_argument(
        "--library-report", required=True, help="Path to library_report.tsv")
    parser_genomes.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for genome .fna files"
    )
    parser_genomes.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "Number of parallel download threads "
            "(default: all CPUs)"
        )
    )
    parser_genomes.set_defaults(func=download_genomes)

    parser_csv = subparsers.add_parser(
        "generate-csv",
        help="Generate CSV datasets from the downloaded genomes"
    )
    parser_csv.add_argument("--min-subseq-len", type=int, default=100,
                            help="Minimum length of subsequences to extract")
    parser_csv.add_argument("--max-subseq-len", type=int, default=512,
                            help="Maximum length of subsequences to extract")
    parser_csv.add_argument("--seed", type=int, default=None,
                            help="Random seed for reproducibility")
    parser_csv.add_argument("--compression", action="store_true",
                            default=False, help="Compress CSV output files")
    parser_csv.add_argument(
        "--min-num-seqs",
        type=int,
        default=5000,
        help="Minimum number of unique subsequences per node"
    )
    parser_csv.add_argument(
        "--percentage",
        type=float,
        default=98.0,
        help=(
            "Percent of nodes with at least given sequence count"
        )
    )
    parser_csv.add_argument(
        "--parallel",
        action="store_true",
        default=True,
        help="Extract/write per node in parallel"
    )
    parser_csv.add_argument("--max-workers", type=int,
                            default=10, help="Max parallel worker processes")
    parser_csv.set_defaults(func=generate_csv)

    return parser


def run_cli(args=None):
    parser = get_parser()
    parse_args = parser.parse_args(args)
    try:
        parse_args.func(parse_args)
    except Exception as e:
        logger.error(f"CLI encountered an error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    run_cli()
