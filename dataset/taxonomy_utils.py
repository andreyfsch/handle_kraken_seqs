"""
taxonomy_utils.py

Summary
-------
This module provides core utilities for extracting and organizing viral
reference genome sequences into a hierarchical taxonomic tree using
Kraken2 and NCBI taxonomy.

Extended Summary
----------------
Key functionalities include:
    - Mapping RefSeq genome identifiers to NCBI taxonomy IDs using
      Kraken2 outputs.
    - Loading genome sequences from FASTA files organized by Kraken2.
    - Building an in-memory taxonomic tree (with the `bigtree` library)
      where each node represents a taxonomic rank, and leaf nodes contain
      actual genome sequences.
    - Traversing and extracting data from the taxonomic tree, including:
        - Extracting representative subsequences from final nodes for
          downstream tasks.
        - Generating sets of all unique subsequences of a given window size
          for deduplication or sequence selection.
    - Robust error handling and logging throughout the genome loading and
      taxonomy resolution process.

The resulting data structures enable downstream dataset creation,
sequence sampling, and taxonomic label assignment for machine learning
workflows (e.g., foundation model fine-tuning).

Dependencies
------------
- Kraken2 database with standard folder structure and map files.
- `taxoniq` for NCBI taxonomy queries.
- `bigtree` for tree representation.
- Standard Python libraries.

Examples
--------
Typical usage:

    >>> from taxonomy_utils import generate_seqs_by_taxon_tree
    >>> tree = generate_seqs_by_taxon_tree()
    >>> # tree is a bigtree.Node. Use traversal or extraction utilities
    >>> # as needed.

Notes
-----
Designed for robust, large-scale, reproducible dataset generation in
computational genomics.

"""
# Standard library imports
import os
import random
import logging
import concurrent.futures
from typing import Dict, Set

# Third-party imports
import tqdm
from bigtree import Node, add_path_to_tree, find_attrs
from progress.bar import Bar

# Local application imports
from dataset import sequence_utils
from config import KRAKEN_PATH, KRAKEN_DATABASE


def get_tax_ids() -> Dict[str, int]:
    """
    Load mapping from RefSeq IDs to NCBI Taxonomy IDs using the
    Kraken2 map file.

    Returns
    -------
    dict of str to int
        Dictionary where keys are RefSeq sequence names (e.g., "NC_002195.1")
        and values are corresponding NCBI Taxonomy IDs.

    Raises
    ------
    FileNotFoundError
        If the Kraken2 seqid2taxid.map file does not exist.

    Notes
    -----
    This function expects the Kraken2 map file to exist at:
    '{KRAKEN_PATH}/{KRAKEN_DATABASE}/seqid2taxid.map'

    Example
    -------
    >>> tax_ids = get_tax_ids()
    >>> tax_ids["NC_002195.1"]
    687377
    """
    tax_ids = {}
    with open(f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/seqid2taxid.map") as f:
        for line in f.readlines():
            # Line example:
            # kraken:taxid|687377|NC_002195.1 687377
            ref_seq = line.rstrip().split("|")[-1].split()[0]
            tax_id = int(line.rstrip().split("|")[-1].split()[-1])
            tax_ids[ref_seq] = tax_id

    return tax_ids


def get_num_refseq_files() -> int:
    """
    Return the number of RefSeq files (genome folders) in the Kraken2
    database structure.

    Returns
    -------
    int
        Number of reference sequence genome directories found under
        '{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes'.

    Notes
    -----
    This function does not check the contents of each folder;
    only their existence.
    """
    logger = logging.getLogger(__name__)
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    if not os.path.isdir(gen_dir):
        logger.error(f"Genome directory {gen_dir} does not exist.")
        raise FileNotFoundError(f"Genome directory {gen_dir} does not exist.")
    return len(
        [
            name
            for name in os.listdir(gen_dir)
            if os.path.isdir(os.path.join(gen_dir, name))
        ]
    )


def parse_refseq_folder(args):
    """
    Helper for multiprocessing: parses genome and taxonomy for one RefSeq.

    Parameters
    ----------
    args : tuple
        (ref_seq, tax_ids, KRAKEN_PATH, KRAKEN_DATABASE)
    Returns
    -------
    dict
        {
            "ref_seq": ...,
            "tax_id": ...,
            "sequences": ...,
            "ranked_taxons": ...,
            "error": None or str
        }
    """
    ref_seq, tax_ids, KRAKEN_PATH, KRAKEN_DATABASE = args
    try:
        import taxoniq  # Local import for worker
        tax_id = tax_ids[ref_seq]
        # Inline genome parsing as before
        sequences = {}
        fasta_path = (
            f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/"
            f"{ref_seq}/genome.fna"
        )
        with open(fasta_path) as f:
            lines = f.readlines()
            subseq = None
            subseq_ref = ""
            subseq_name = ""
            for i in range(len(lines)):
                if lines[i].startswith(">"):
                    if subseq:
                        sequences[subseq_ref] = {
                            "seq_name": subseq_name, "seq": subseq}
                    line_split = lines[i].split()
                    subseq_ref = line_split[0][1:]
                    subseq_name = " ".join(line_split[1:])
                    subseq = ""
                elif i == len(lines) - 1:
                    subseq += lines[i].rstrip()
                    sequences[subseq_ref] = {
                        "seq_name": subseq_name, "seq": subseq}
                else:
                    subseq += lines[i].rstrip()
        t = taxoniq.Taxon(tax_id)
        ranked_taxons = t.ranked_lineage
        return {
            "ref_seq": ref_seq,
            "tax_id": tax_id,
            "sequences": sequences,
            "ranked_taxons": ranked_taxons,
            "error": None
        }
    except Exception as e:
        return {"ref_seq": ref_seq, "error": str(e)}


def generate_seqs_by_taxon_tree() -> Node:
    """
    Build a hierarchical taxonomic tree populated with sequences from
    Kraken2 RefSeq genome files.

    This function:
        - Scans the Kraken2 genomes directory (as specified in config) for
          all reference sequence folders.
        - For each reference sequence:
            - Loads all FASTA sequences for that genome.
            - Resolves its full NCBI taxonomy lineage using the
              `taxoniq` library.
            - Builds a hierarchical tree (using bigtree's `Node`)
              representing the taxonomy, with each taxonomic rank as a node
              and sequences as leaf nodes.
            - Attaches each sequence under the appropriate taxon in the tree,
              branching down to the reference sequence and sequence leaves.
        - Avoids duplicate nodes and sequence entries via internal tracking.
        - Reports progress using a terminal progress bar.

    Returns
    -------
    Node
        The root node of the constructed taxonomic tree.
        Each node includes attributes for taxonomic rank, scientific name,
        and (for sequence nodes) sequence data and description.

    Raises
    ------
    FileNotFoundError
        If the Kraken2 reference files or required taxonomy maps are missing.
    KeyError
        If a reference sequence cannot be mapped to a taxonomy ID.
    Exception
        For other errors in taxonomy lookup or file reading.

    See Also
    --------
    get_tax_ids : Loads RefSeq-to-taxonomy mapping.
    get_genome_sequences : Loads genome FASTA files into dicts.
    bigtree.Node : Tree structure used for the hierarchy.

    Notes
    -----
    - Taxonomy is resolved using the `taxoniq` package.
    - The function depends on the configuration of `KRAKEN_PATH` and
      `KRAKEN_DATABASE`.
    - Sequences are attached to leaves with node attributes: 'seq_name',
      'seq'.
    - Uses a progress bar to report status during tree construction.
    - Suitable for viral RefSeq collections, but can be adapted for other
      taxonomic databases.
    - Designed for robust, large-scale, reproducible dataset generation in
      computational genomics.
    - Duplicates and missing data are logged and skipped, so the function is
      robust to partial failures in the underlying data structure.

    Examples
    --------
    >>> from taxonomy_utils import generate_seqs_by_taxon_tree
    >>> tree = generate_seqs_by_taxon_tree()
    >>> # tree is a bigtree.Node object. Use bigtree utilities to traverse
    >>> # or extract data.
    """
    logger = logging.getLogger(__name__)
    taxon_tree = Node("root")
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    try:
        tax_ids = get_tax_ids()
    except FileNotFoundError as e:
        logger.error(
            f"Kraken2 seqid2taxid.map file not found at expected location: {e}"
        )
        raise
    except Exception as e:
        logger.error(f"Unexpected error loading tax IDs: {e}")
        raise

    if not os.path.isdir(gen_dir):
        logger.error(f"Genome directory {gen_dir} does not exist.")
        raise FileNotFoundError(f"Genome directory {gen_dir} does not exist.")

    refseq_dirs = [
        name for name in os.listdir(gen_dir)
        if os.path.isdir(os.path.join(gen_dir, name))
    ]

    num_refseq_files = len(refseq_dirs)
    visited_nodes = set()
    registered_subseq_refs = set()

    # Multiprocessing: Use all available CPUs
    max_workers = os.cpu_count() or 1

    results = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers
    ) as executor:
        # Prepare job args (pass all dependencies explicitly to avoid
        # fork/import bugs)
        jobs = [
            (ref_seq, tax_ids, KRAKEN_PATH, KRAKEN_DATABASE)
            for ref_seq in refseq_dirs
        ]
        # Show progress bar for parallel work
        for result in tqdm.tqdm(
            executor.map(parse_refseq_folder, jobs),
            total=len(jobs),
            desc="Loading genomes+taxonomy in parallel"
        ):
            results.append(result)

    bar = Bar("Building taxon tree", max=num_refseq_files)
    for item in results:
        if item.get("error"):
            logger.warning(f"Error for {item['ref_seq']}: {item['error']}")
            bar.next()
            continue
        ref_seq = item["ref_seq"]
        ranked_taxons = item["ranked_taxons"]
        sequences = item["sequences"]
        path_parent_rank = ""
        for idx, ranked_taxon in enumerate(reversed(ranked_taxons)):
            slash = "" if path_parent_rank == "" else "/"
            path_parent_rank += slash + ranked_taxon.scientific_name
            try:
                if path_parent_rank in visited_nodes:
                    already_added = True
                    for subseq_ref in sequences.keys():
                        if subseq_ref not in registered_subseq_refs:
                            already_added = False

                    if idx == len(ranked_taxons) - 1 and not already_added:
                        add_path_to_tree(
                            taxon_tree,
                            f"{path_parent_rank}/{ref_seq}",
                            node_attrs={"rank": "ref_seq"},
                        )
                        for subseq_ref, seq_dict in sequences.items():
                            add_path_to_tree(
                                taxon_tree,
                                (
                                    f"{path_parent_rank}/{ref_seq}/"
                                    f"{subseq_ref}"
                                ),
                                node_attrs={
                                    "rank": "sequence",
                                    "seq_name": seq_dict["seq_name"],
                                    "seq": seq_dict["seq"]
                                },
                            )
                            registered_subseq_refs.add(subseq_ref)
                    else:
                        continue
                else:
                    if ranked_taxon.rank.name == "superkingdom":
                        visited_nodes |= {path_parent_rank}
                        taxon_tree = Node.from_dict({
                            "name": ranked_taxon.scientific_name,
                            "rank": ranked_taxon.rank.name
                        })
                    else:
                        visited_nodes |= {path_parent_rank}
                        add_path_to_tree(
                            taxon_tree,
                            path_parent_rank,
                            node_attrs={"rank": ranked_taxon.rank.name},
                        )
                        already_added = True
                        for subseq_ref in sequences.keys():
                            if subseq_ref not in registered_subseq_refs:
                                already_added = False
                        if idx == len(ranked_taxons) - 1 and not already_added:
                            add_path_to_tree(
                                taxon_tree,
                                f"{path_parent_rank}/{ref_seq}",
                                node_attrs={
                                    "rank": "ref_seq",
                                },
                            )
                            for subseq_ref, seq_dict in sequences.items():
                                add_path_to_tree(
                                    taxon_tree,
                                    (
                                        f"{path_parent_rank}/{ref_seq}/"
                                        f"{subseq_ref}"
                                    ),
                                    node_attrs={
                                        "rank": "sequence",
                                        "seq": seq_dict["seq"],
                                        "seq_name": seq_dict["seq_name"]
                                    },
                                )
                                registered_subseq_refs.add(subseq_ref)
            except Exception as e:
                logger.error(
                    f"Error updating taxon tree at {path_parent_rank}: {e}")
                continue
        bar.next()
    bar.finish()
    return taxon_tree


def get_genome_sequences(ref_seq: str) -> Dict[str, Dict[str, str]]:
    """
    Load all sequences for a specific RefSeq identifier from its genome FASTA
    file.

    Parameters
    ----------
    ref_seq : str
        Reference sequence identifier (e.g., "NC_019947.1").

    Returns
    -------
    dict of str to dict
        Dictionary mapping subsequence IDs to dictionaries with keys:
        - 'seq_name': str, the full description from the FASTA header.
        - 'seq': str, the DNA sequence.

    Raises
    ------
    FileNotFoundError
        If the genome FASTA file does not exist.

    Notes
    -----
    This function expects FASTA files named 'genome.fna' within each
    genome directory:
    '{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/{ref_seq}/genome.fna'

    Example
    -------
    >>> genome = get_genome_sequences("NC_019947.1")
    >>> list(genome.keys())
    ['NC_019947.1']
    >>> genome['NC_019947.1']['seq_name']
    'Tomato yellow mottle virus segment DNA-B, complete sequence'
    """
    sequences = {}
    # Start of fasta sequence example:
    # >NC_019947.1 Tomato yellow mottle virus segment DNA-B, complete sequence
    with open(
        f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/genomes/{ ref_seq }/genome.fna"
    ) as f:
        lines = f.readlines()
        subseq = None
        subseq_ref = ""
        subseq_name = ""
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                if subseq:
                    sequences[subseq_ref] = {"seq_name": subseq_name,
                                             "seq": subseq}
                line_split = lines[i].split()
                subseq_ref = line_split[0][1:]
                subseq_name = " ".join(line_split[1:])
                subseq = ""
            elif i == len(lines) - 1:
                subseq += lines[i].rstrip()
                sequences[subseq_ref] = {"seq_name": subseq_name,
                                         "seq": subseq}
            else:
                subseq += lines[i].rstrip()

    return sequences


def get_subseqs_from_final_node(
        final_node: Node, n: int, min_len: int = 100,
        max_len: int = 512, rng: random.Random | None = None,
        parallel: bool = True, max_workers: int | None = None) -> list[str]:
    """
    Extract representative subsequences from a taxon node, in parallel if
    requested.

    For each sequence leaf in the final node, determines the proportional
    number of subsequences to sample based on sequence length, and samples
    subsequences (using a random window) for each.

    When `parallel=True`, the extraction for each sequence is performed in
    parallel across available CPU cores, using
    `concurrent.futures.ProcessPoolExecutor`.

    Parameters
    ----------
    final_node : Node
        The final taxonomic node from which to extract sequences.
    n : int
        Number of subsequences to extract (total across all sequences).
    min_len : int, optional
        Minimum subsequence length (default is 100).
    max_len : int, optional
        Maximum subsequence length (default is 512).
    rng : random.Random, optional
        Random number generator instance for reproducibility. If None,
        uses global randomness (note: for parallel mode, each worker
        receives its own default RNG unless you provide a seed per
        sequence).
    parallel : bool, optional
        If True, extract subsequences from each sequence in parallel using
        all available CPU cores.
        If False, extraction is performed serially (default is True).
    max_workers : int, optional
        The maximum number of parallel workers to use (default: all CPUs).

    Returns
    -------
    list of str
        List of extracted DNA subsequences from all sequences under
        `final_node`.

    Notes
    -----
    - The total number of subsequences is always exactly `n`, distributed
      across all sequence leaves proportionally to their length. Counts are
      rounded and compensated as needed.
    - When `parallel=True`, each sequence is processed independently in a
      separate process.
    - If the number of sequence leaves is small, parallelization may
      provide little speedup.
    - Reproducibility: For strict reproducibility in parallel mode, you
      must manage random seeds yourself per sequence. Otherwise, worker
      RNGs are independent.
    - Suitable for large taxonomic nodes with many constituent sequences.

    Examples
    --------
    Extract 10 subsequences from all sequences under `node`, in parallel:

    >>> subseqs = get_subseqs_from_final_node(
    ...     node, 10, min_len=150, max_len=400, parallel=True)
    >>> print(subseqs[0])
    'AGGCTT...'

    Or extract serially (for debugging):

    >>> subseqs = get_subseqs_from_final_node(
    ...     node, 10, min_len=150, max_len=400, parallel=False)
    >>> print(subseqs[-1])
    'TTGGCA...'
    """
    sequence_nodes = find_attrs(final_node, "rank", "sequence")

    sequences = {}
    n_seqs = {}
    total_len = 0
    for sequence_node in sequence_nodes:
        ref_seq = sequence_node.node_name
        sequence = sequence_node.get_attr("seq")
        sequences[ref_seq] = sequence
        total_len += len(sequence)
    if total_len == 0:
        return []
    for ref_seq, sequence in sequences.items():
        fraction = len(sequence) / total_len
        n_seqs[ref_seq] = int(round(fraction * n, 0))
    # compensate rounding errors
    if sum(n_seqs.values()) != n:
        while sum(n_seqs.values()) != n:
            if sum(n_seqs.values()) < n:
                n_seqs[min(n_seqs, key=n_seqs.get)] += 1
            elif sum(n_seqs.values()) > n:
                n_seqs[max(n_seqs, key=n_seqs.get)] -= 1
    for ref_seq, sequence in sequences.items():
        if n_seqs[ref_seq] == 0:
            if len(sequence) >= min_len:
                n_seqs[ref_seq] = 1
                n_seqs[max(n_seqs, key=n_seqs.get)] -= 1
            else:
                continue

    # Assign a unique deterministic seed for each sequence/job
    if rng is not None and hasattr(rng, 'randrange'):
        base_seed = rng.randrange(1 << 30)
    else:
        base_seed = random.randrange(1 << 30)

    # Prepare arguments for each parallel call
    jobs = [
        # or pass in rng/seed as desired
        (sequence, n_seqs[ref_seq], min_len, max_len, base_seed + i)
        for i, (ref_seq, sequence) in enumerate(sequences.items())
        if n_seqs[ref_seq] > 0
    ]

    subseqs = []

    if parallel:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            results = list(executor.map(extract_subseqs_worker, jobs))
        for res in results:
            subseqs.extend(res)
    else:
        for args in jobs:
            subseqs.extend(sequence_utils.extract_subseqs(*args))

    return subseqs


def extract_subseqs_worker(args: tuple) -> list[str]:
    """
    Parallel worker function for extracting random DNA subsequences
    from a sequence.

    This function is designed to be used as a worker target for
    `concurrent.futures.ProcessPoolExecutor` or similar parallelization
    frameworks. It unpacks arguments and calls `extract_subseqs`, returning
    a list of subsequences extracted from the input DNA sequence.

    Parameters
    ----------
    args : tuple
        A tuple containing:
            - sequence (str): The full DNA sequence to sample from.
            - n (int): Number of subsequences to extract.
            - min_len (int): Minimum length for each subsequence.
            - rng (random.Random or None): Optional random number generator
              instance for reproducibility. If None, uses the global RNG.
              for reproducibility. If None, uses the global RNG.

    Returns
    -------
    list of str
        List of extracted DNA subsequences.

    Notes
    - This function is required for compatibility with Python's
      multiprocessing, which cannot pickle lambda functions or
      locally-defined functions.
      which cannot pickle lambda functions or locally-defined functions.
    - Each worker receives its own arguments tuple; random seed reproducibility
      is only guaranteed if a seeded RNG or seed is passed per job.
    - Intended for use with `executor.map(extract_subseqs_worker, jobs)` where
      `jobs` is a list of argument tuples as described above.
    - Handles exceptions by propagating them up to the parent process;
      any errors will terminate the corresponding parallel job.

    Examples
    --------
    >>> args = ("ACGTACGTACGT", 5, 3, 6, None)
    >>> result = extract_subseqs_worker(args)
    >>> print(result)
    ['CGTACG', 'TACGTA', 'ACGTAC', 'GTACGT', 'TACGTA']
    """
    return sequence_utils.extract_subseqs(*args)


def extract_max_subseqs_set(
    final_taxon_node: Node, window_size: int
) -> Set[str]:
    """
    Generate a set of all unique possible subsequences of a given window size
    from all sequences under the specified taxonomic node.

    Parameters
    ----------
    final_taxon_node : Node
        Taxonomic node containing sequence leaves.
    window_size : int
        Length of the sliding window for subsequence extraction.

    Returns
    -------
    set of str
        Set of unique subsequences of the specified length.

    Notes
    -----
    This is used to deduplicate and analyze sequence content for a taxonomic
    group.

    Example
    -------
    >>> unique_seqs = extract_max_subseqs_set(node, 100)
    >>> len(unique_seqs)
    12345
    """
    max_subseqs = set()
    seqs = find_attrs(final_taxon_node, "rank", "sequence")

    for sequence_node in seqs:
        sequence = sequence_node.get_attr("seq")
        for i in range(0, len(sequence) - window_size + 1):
            max_subseqs.add(sequence[i:i+window_size])

    return max_subseqs
