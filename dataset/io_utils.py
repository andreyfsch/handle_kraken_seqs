"""
io_utils.py

Summary
-------
Tools for dataset output and CSV writing for downstream ML training and
evaluation.
Handles taxonomic-level dataset construction, deduplication, and robust
CSV export.

Functions
    Main pipeline to generate CSVs of sequences/labels for each taxonomic
    level.
write_csvs
    Main pipeline to generate CSVs of sequences/labels for each taxonomic
    level.

Notes
-----
Designed for reproducibility, robustness, and safe concurrent data writing.
"""
import pandas as pd
import numpy as np
import os
import random
import logging
import signal
import concurrent.futures
import multiprocessing as mp
import glob
from tqdm.auto import tqdm
from bigtree import find_attrs
from dataset.taxonomy_utils import (
    generate_seqs_by_taxon_tree, extract_max_subseqs_set,
    get_subseqs_from_final_node
)
from config import CSV_OUTPUT_PATH

DEFAULT_NODE_TIMEOUT = int(os.environ.get("CSV_NODE_TIMEOUT_SEC", "120"))


def _init_worker():
    # keep BLAS single-threaded per process to avoid oversubscription
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    try:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    except Exception:
        pass


class _Timeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _Timeout()


# Set up module logger
logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format=(
            "%(asctime)s %(levelname)s:%(name)s: %(message)s"
        )
    )


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = []
        try:
            for _ in range(size):
                chunk.append(next(it))
        except StopIteration:
            if chunk:
                yield chunk
            break
        yield chunk


def _extract_and_write_node_guarded(job, timeout_s=DEFAULT_NODE_TIMEOUT):
    old = signal.signal(signal.SIGALRM, _alarm_handler)
    # Use setitimer for sub-second timers and to avoid drift
    signal.setitimer(signal.ITIMER_REAL, timeout_s)
    try:
        return _extract_and_write_node(job)
    except _Timeout:
        node = job[0]
        idx = job[4]
        logger.error(
            f"Timeout for node '{getattr(node, 'node_name', '?')}' "
            f"({idx}) after {timeout_s}s; skipping."
        )
        return (node.node_name, idx, None, False, f"timeout {timeout_s}s")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old)


def process_chunk(jobs_chunk):
    return [_extract_and_write_node_guarded(job) for job in jobs_chunk]


def _extract_and_write_node(args):
    (
        node,
        num_seqs_extraction,
        min_subseq_len,
        max_subseq_len,
        idx,
        output_path,
        rng,
        compression
    ) = args
    """
    Helper for parallel node processing: extract subseqs, write DataFrame,
    and return label info.
    """
    try:
        virus_dict = {"sequence": [], "label": []}
        subseqs = get_subseqs_from_final_node(
            node,
            num_seqs_extraction,
            min_subseq_len,
            max_subseq_len,
            rng,
            parallel=False
        )
        for subseq in subseqs:
            virus_dict["sequence"].append(subseq)
            virus_dict["label"].append(idx)
        part_path = f"{output_path}.part{idx:05d}"
        df = pd.DataFrame.from_dict(virus_dict)
        if compression:
            df.to_csv(part_path, index=False, compression='gzip')
        else:
            df.to_csv(part_path, index=False)
        return (node.node_name, idx, part_path, True, None)
    except Exception as e:
        logger.error(f"Failed to process node '{node.node_name}' ({idx}): {e}")
        return (node.node_name, idx, None, False, str(e))


def write_csvs(
        min_subseq_len: int = 100,
        max_subseq_len: int = 512,
        rng: random.Random | None = None,
        compression: bool = True,
        min_num_seqs: int = 5000,
        percentage: float = 98.0,
        parallel: bool = True,
        max_workers: int | None = None,
) -> None:
    """
    Generate and write CSV files containing DNA subsequences and taxon labels
    for ML workflows, with robust logging, error handling, and optional
    parallel execution.

    For each taxonomic level (species, genus, etc.), the function determines
    the number of unique possible subsequences available for each taxon node.
    A selection threshold and per-node sequence extraction count are set as
    follows:

      - If all nodes have at least `min_num_seqs` possible subsequences, the
        smallest maximum among nodes is chosen as the number of sequences to
        extract per node, and all nodes are kept.
      - If some nodes have fewer than `min_num_seqs`, the number of sequences
        extracted per node is set as the largest value for which at least
        `percentage`% of nodes have at least that many sequences. Nodes below
        this threshold are skipped.

    Each qualifying node has the determined number of subsequences extracted
    (with labeling), and written as a temporary CSV. All part files are then
    concatenated into a single output CSV per taxonomic level. Optionally,
    extraction and CSV writing per node can be parallelized across CPUs.

    Parameters
    ----------
    min_subseq_len : int, optional
        Minimum length of subsequences to extract (default: 100).
    max_subseq_len : int, optional
        Maximum length of subsequences to extract (default: 512).
    rng : random.Random, optional
        Random number generator instance for reproducibility. If None,
        uses global randomness.
    compression : bool, optional
    min_num_seqs : int, optional
        Minimum number of unique subsequences required per node
        (default: 5000).
    percentage : float, optional
        Percentage of nodes that must have at least the extracted sequence
        count (default: 98.0). Nodes below this threshold are excluded from
        the output.
    parallel : bool, optional
        If True, extract and write each node in parallel using all available
        CPUs (default: True). If False, run in serial for debugging or
        low-memory.
    max_workers : int or None, optional
        The maximum number of parallel worker processes to use
        (default: all CPUs).

    Returns
    -------
    None
        Writes output files to disk.

    Output Files
    ------------
    {CSV_OUTPUT_PATH}/{taxonomic_level}_seqs.csv[.gz] :
        Contains columns: 'sequence' (the DNA subsequence) and 'label'
        (numeric taxon label). May be compressed if `compression=True`.
    - Extraction and CSV writing for each node is parallelized if
      `parallel=True`.
    - Each node’s output is written to a temporary part file, then all parts
      are merged.
    - Any node failing extraction or writing is logged and skipped; the
      pipeline continues.
    - Only nodes meeting the threshold are included for each level; rare
      nodes may be dropped.
    - Label indices are 1-based and consistent within each level.
    - Compression uses `.csv.gz` extension when enabled.
    - Taxon-label mapping is always written to a separate file and not
      merged with sequence data.
    - For strict reproducibility, pass a seeded `random.Random` to `rng`.
    - Only nodes meeting the threshold are included for each level; rare nodes
        may be dropped.
    - Label indices are 1-based and consistent within each level.
    - Compression uses `.csv.gz` extension when enabled.
    - Taxon-label mapping is always written to a separate file and not merged
        with sequence data.
    - For strict reproducibility, pass a seeded `random.Random` to `rng`.

    Examples
    --------
    >>> from io_utils import write_csvs
    >>> import random
    >>> write_csvs(
    ...     min_subseq_len=150,
    ...     max_subseq_len=300,
    ...     min_num_seqs=5000,
    ...     percentage=98.0,
    ...     rng=random.Random(42),
    ...     compression=True,
    ...     parallel=True,
    ...     max_workers=8
    ... )
    >>> # Output: species_seqs.csv.gz, genus_seqs.csv.gz, ...,
    ... # plus species_seqs_taxon_labels.csv, etc.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s:%(name)s: %(message)s"
    )

    os.makedirs(CSV_OUTPUT_PATH, exist_ok=True)
    logger.info("Starting tree construction.")
    tree = generate_seqs_by_taxon_tree()
    logger.info("Tree built. Starting dataset extraction.")
    taxonomic_levels = ["species", "genus", "family",
                        "order", "class", "phylum", "kingdom"]

    for taxonomic_level in taxonomic_levels:
        logger.info(f"Processing taxonomic level: {taxonomic_level}")
        taxon_names2ids = {taxonomic_level: [], "id": []}
        all_taxonomic_nodes = find_attrs(tree, "rank", taxonomic_level)
        all_max_len_subseqs = []
        node_name_to_max = {}

        for node in all_taxonomic_nodes:
            max_subseqs = extract_max_subseqs_set(node, min_subseq_len)
            n_seqs = len(max_subseqs)
            all_max_len_subseqs.append(n_seqs)
            node_name_to_max[node.path_name + '/' + node.node_name] = n_seqs

        if not all_max_len_subseqs:
            logger.warning(
                f"No nodes found at level {taxonomic_level}; skipping.")
            continue

        min_max = min(all_max_len_subseqs)

        # Logic: If every node has at least min_num_seqs, use min_max
        if min_max >= min_num_seqs:
            num_seqs_extraction = min_max
            selected_nodes = [
                node for node in all_taxonomic_nodes
                if node_name_to_max[
                    node.path_name + '/' + node.node_name
                ] >= num_seqs_extraction
            ]
        else:
            cutoff = int(np.percentile(all_max_len_subseqs, 100 - percentage))
            num_seqs_extraction = cutoff
            selected_nodes = [
                node for node in all_taxonomic_nodes
                if (
                    node_name_to_max[node.path_name + '/' + node.node_name]
                    >= num_seqs_extraction
                )
            ]
        if compression:
            output_path = f"{CSV_OUTPUT_PATH}/{taxonomic_level}_seqs.csv.gz"
        else:
            output_path = f"{CSV_OUTPUT_PATH}/{taxonomic_level}_seqs.csv"

        # Parallel extraction/writing per node
        jobs = []
        for idx, node in enumerate(selected_nodes, start=1):
            taxon_names2ids[taxonomic_level].append(node.node_name)
            taxon_names2ids["id"].append(idx)
            jobs.append((
                node, num_seqs_extraction, min_subseq_len, max_subseq_len, idx,
                output_path, rng, compression))
        results = []
        if parallel and jobs:
            if max_workers is None:
                max_workers = min(32, os.cpu_count() or 1)

            CHUNK_SIZE = int(os.environ.get("CSV_CHUNK_SIZE", "50"))
            job_chunks = list(chunked_iterable(jobs, CHUNK_SIZE))

            # process chunks in waves so workers get respawned periodically
            # ~10k nodes per wave with CHUNK_SIZE=50
            WAVE_SIZE = int(os.environ.get("CSV_WAVE_CHUNKS", "200"))

            remaining = list(range(len(job_chunks)))
            pbar = tqdm(
                total=len(jobs),
                desc="Extracting & writing",
                unit="node"
            )

            while remaining:
                wave_idx = remaining[:WAVE_SIZE]
                remaining = remaining[WAVE_SIZE:]

                ctx = mp.get_context("spawn")
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=max_workers,
                    mp_context=ctx,
                    initializer=_init_worker,
                ) as executor:
                    future_map = {executor.submit(
                        process_chunk, job_chunks[i]): i for i in wave_idx}
                    completed_this_wave = []

                    try:
                        for fut in concurrent.futures.as_completed(future_map):
                            i = future_map[fut]
                            chunk_results = fut.result()
                            results.extend(chunk_results)
                            pbar.update(len(chunk_results))
                            completed_this_wave.append(i)
                    except Exception as e:
                        # A worker died (BrokenProcessPool or similar).
                        # Drop back the unfinished chunk indices into
                        # 'remaining' to retry next wave.
                        failed = set(wave_idx) - set(completed_this_wave)
                        remaining = list(failed) + remaining
                        # Optional: log the reason; continue to next wave
                        logger.error(
                            (
                                f"Wave aborted due to worker failure: {e}. "
                                f"Retrying {len(failed)} chunks in next wave."
                            )
                        )

            pbar.close()
        else:
            for job in tqdm(jobs, desc="Extracting & writing", unit="node"):
                results.append(_extract_and_write_node(job))

        logger.info(
            f"Extracted {num_seqs_extraction} sequences from "
            f"{len(jobs)} nodes at level {taxonomic_level} "
            f"with {max_workers} workers."
        )

        # Concatenate all part files into final output
        part_files = sorted(glob.glob(f"{output_path}.part*"))
        frames = []
        for part_file in part_files:
            if os.path.exists(part_file):
                try:
                    frames.append(
                        pd.read_csv(
                            part_file,
                            compression='gzip' if compression else None
                        )
                    )
                except Exception as e:
                    logger.error(f"Could not read {part_file}: {e}")
        if frames:
            final_df = pd.concat(frames, ignore_index=True)
            final_df.to_csv(output_path, index=False,
                            compression='gzip' if compression else None)
            logger.info(
                f"Wrote final CSV: {output_path} ({len(final_df)} rows)")
            # Clean up temp part files
            for part_file in part_files:
                try:
                    os.remove(part_file)
                except Exception as e:
                    logger.warning(
                        f"Could not delete temp file {part_file}: {e}")
        else:
            logger.warning(
                f"No valid data frames to write for level {taxonomic_level}.")
        mapping_path = output_path.replace(
            '.csv.gz', '_taxon_labels.csv'
        ).replace(
            '.csv', '_taxon_labels.csv'
        )
        taxon_df = pd.DataFrame.from_dict(taxon_names2ids)
        mapping_path = (
            output_path.replace('.csv.gz', '_taxon_labels.csv')
            .replace('.csv', '_taxon_labels.csv')
        )
        try:
            taxon_df.to_csv(mapping_path, index=False)
        except Exception as e:
            logger.error(
                f"Failed to write taxon label mapping to {mapping_path}: {e}"
            )
        logger.info(f"Wrote label mapping CSV: {mapping_path}")

    logger.info("All taxonomic levels processed.")
