import math
import os
import random
from typing import Dict, List, Set
import numpy as np
import pandas as pd

# from datetime import datetime
# import pandas as pd
# from zoneinfo import ZoneInfo
import taxoniq
from bigtree import Node, add_path_to_tree, find_attrs
from progress.bar import Bar

KRAKEN_PATH = "/home/andrey/generate_kraken_dataset/kraken2"
KRAKEN_DATABASE = "viral"


def extract_subseqs(seq, n, min_len, max_len):
    subseqs = []
    # extract non overlapping subsequences of max_len randomly
    if len(seq) >= 2 * n * max_len:
        blacklist = []
        while len(subseqs) < n:
            idx = random.randrange(0, len(seq) - max_len + 1)
            try_again = False
            for pos in blacklist:
                if idx >= pos and idx <= pos + max_len:
                    try_again = True
                    break
            if try_again:
                continue
            else:
                blacklist.append(idx)
                subseqs.append(seq[idx: idx + max_len])
    # extract non adjacent overlapping subsequences of max_len (window between)
    elif len(seq) < 2 * n * max_len and len(seq) >= n * max_len:
        rest = (len(seq) / max_len) - n
        rest_bases = int(math.floor(len(seq) / rest))
        window_bases = int((rest_bases / n) - 1)
        left_start = 0
        right_start = len(seq) - max_len
        if n % 2 == 0:
            operations = int(n / 2)
        else:
            operations = int((n - 1) / 2)
        for i in range(operations):
            subseqs.append(seq[left_start: left_start + max_len])
            left_start += max_len + window_bases
            subseqs.append(seq[right_start: right_start + max_len])
            right_start -= max_len + window_bases
        if n % 2 != 0:
            mid_seq = int(math.floor(len(seq) / 2))
            mid_max_len = int(math.floor(max_len / 2))
            subseqs.append(seq[mid_seq - mid_max_len: mid_seq + mid_max_len])
    # extract overlapping subsequences of length between min_len and max_len
    else:
        iter = 0
        iter_limit = n * 2
        while len(subseqs) < n:
            iter += 1
            if iter >= iter_limit:
                max_len_subseq = max(subseqs, key=len)
                max_len_subseqs = (i for i, x in enumerate(subseqs)
                                   if len(x) == max_len_subseq)
                subseqs.pop(max_len_subseqs.pop())
            subseq_len = random.randrange(min_len, min(max_len, len(seq)))
            idx = random.randrange(0, len(seq) - subseq_len + 1)
            subseq = seq[idx: (idx + subseq_len)]
            complement = get_complement(subseq)
            if subseq not in subseqs:
                subseqs.append(subseq)
            if complement not in subseqs:
                subseqs.append(complement)
        # compensate double append when n is odd
        if len(subseqs) > n:
            subseqs.pop()
    return subseqs


def get_complement(seq: str) -> str:
    """Return the complement sequence of a given seq."""
    # IUPAC notation https://doi.org/10.1021/bi00822a023
    complementary_bases = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "Y": "R",
        "R": "Y",
        "W": "W",
        "S": "S",
        "K": "M",
        "M": "K",
        "D": "H",
        "H": "D",
        "V": "B",
        "B": "V",
        "N": "N",
    }
    return "".join([complementary_bases[i] for i in seq])[::-1]


def get_tax_ids() -> Dict:
    """Return a dict with RefSeq as key and NCBI Taxonomy ID as value obtained
    by accessing the map file provided by Kraken2."""
    tax_ids = {}
    with open(f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/seqid2taxid.map") as f:
        for line in f.readlines():
            # kraken:taxid|687377|NC_002195.1 687377
            ref_seq = line.rstrip().split("|")[-1].split()[0]
            tax_id = int(line.rstrip().split("|")[-1].split()[-1])
            tax_ids[ref_seq] = tax_id
    f.close()

    return tax_ids


def get_num_refseq_files() -> int:
    """Return the number of RefSeq files contained in
    the kraken2 file structure."""
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    return len(
        [
            name
            for name in os.listdir(gen_dir)
            if os.path.isdir(os.path.join(gen_dir, name))
        ]
    )


def generate_seqs_by_taxon_tree() -> Node:
    """Return hierarchical taxon tree according to Taxoniq NCBI taxonomy
    containing RefSeq kraken2 files' sequences, branching them on
    SubRefSeqs leaves."""
    taxon_tree = Node("root")
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    tax_ids = get_tax_ids()

    file_tree = os.walk(gen_dir)
    file_tree.__next__()

    num_refseq_files = get_num_refseq_files()
    visited_nodes = set()
    registered_subseq_refs = set()
    bar = Bar("Loading taxon tree", max=num_refseq_files)

    for subdir, _, _ in file_tree:
        ref_seq = subdir.split("/")[-1]
        tax_id = tax_ids[ref_seq]

        sequences = get_genome_sequences(ref_seq)
        t = taxoniq.Taxon(tax_id)
        ranked_taxons = t.ranked_lineage
        path_parent_rank = ""

        for idx, ranked_taxon in enumerate(reversed(ranked_taxons)):
            slash = "" if path_parent_rank == "" else "/"
            path_parent_rank += slash + ranked_taxon.scientific_name

            if path_parent_rank in visited_nodes:
                already_added = True
                for subseq_ref in sequences.keys():
                    if subseq_ref not in registered_subseq_refs:
                        already_added = False

                if idx == len(ranked_taxons) - 1 and not already_added:
                    add_path_to_tree(
                        taxon_tree,
                        f"{ path_parent_rank }/{ ref_seq }",
                        node_attrs={"rank": "ref_seq"},
                    )
                    for subseq_ref, seq_dict in sequences.items():
                        add_path_to_tree(
                            taxon_tree,
                            f"{ path_parent_rank }/{ ref_seq }/{ subseq_ref }",
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
                            f"{ path_parent_rank }/{ ref_seq }",
                            node_attrs={
                                "rank": "ref_seq",
                            },
                        )
                        for subseq_ref, seq_dict in sequences.items():
                            add_path_to_tree(
                                taxon_tree,
                                f"{ path_parent_rank }/{ ref_seq }/" +
                                f"{ subseq_ref }",
                                node_attrs={
                                    "rank": "sequence",
                                    "seq": seq_dict["seq"],
                                    "seq_name": seq_dict["seq_name"]
                                },
                            )
                            registered_subseq_refs.add(subseq_ref)
        bar.next()
    bar.finish()
    return taxon_tree


def get_genome_sequences(ref_seq: str) -> Dict:
    """Return dict of subseq_ref sequences retrieved from
    kraken2 genome fasta file identified by ref_seq."""
    sequences = {}
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
    f.close()

    return sequences


def get_subseqs_from_final_node(
        final_node: Node, n: int, min_len: int = 100,
        max_len: int = 512) -> List:
    """Return list of n subsequences extracted from final_node ranging from
    min_len to max_len length."""
    sequence_nodes = find_attrs(final_node, "rank", "sequence")

    sequences = {}
    subseqs = []
    n_seqs = {}
    total_len = 0
    for sequence_node in sequence_nodes:
        ref_seq = sequence_node.node_name
        sequence = sequence_node.get_attr("seq")
        sequences[ref_seq] = sequence
        total_len += len(sequence)
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
    for ref_seq, sequence in sequences.items():
        subseqs.extend(extract_subseqs(
            sequence, n_seqs[ref_seq], min_len, max_len
        ))

    return subseqs


def extract_max_subseqs_set(final_taxon_node: Node, window_size: int) -> Set:
    """Return a set with all unique possible sequences of a given window size
    length of a given final taxon node."""
    max_subseqs = set()
    seqs = find_attrs(final_taxon_node, "rank", "sequence")

    for sequence_node in seqs:
        sequence = sequence_node.get_attr("seq")
        for i in range(0, len(sequence) - window_size):
            max_subseqs.add(sequence[i:i+window_size-1])

    return max_subseqs


def write_csvs(min_subseq_len=100, max_subseq_len=512):
    taxon_names2ids = {}
    virus_dict = {"sequence": [], "label": []}

    tree = generate_seqs_by_taxon_tree()

    taxonomic_levels = []
    generic_taxon = taxoniq.Taxon(100)
    for generic_taxon_level in generic_taxon.ranked_lineage:
        taxonomic_levels.append(generic_taxon_level.rank.name)

    for taxonomic_level in taxonomic_levels:
        taxon_names2ids = {taxonomic_level: [], "id": []}
        all_taxonomic_nodes = find_attrs(tree, "rank", taxonomic_level)
        all_max_len_subseqs = []
        target_node_paths_names = []

        for node in all_taxonomic_nodes:
            max_subseqs = extract_max_subseqs_set(node, min_subseq_len)
            all_max_len_subseqs.append(len(max_subseqs))

        mean = np.mean(all_max_len_subseqs)
        std = np.std(all_max_len_subseqs)
        num_seqs_extraction = int(math.ceil((mean + std) * 0.01))

        for node in all_taxonomic_nodes:
            max_subseqs = extract_max_subseqs_set(node, min_subseq_len)
            if len(max_subseqs) >= num_seqs_extraction:
                target_node_paths_names.append(
                    node.path_name+'/'+node.node_name)

        output_path = f"{ taxonomic_level }_seqs.csv"
        taxon_names2ids_path = f"{ taxonomic_level }_seqs.csv"
        idx = 1
        for node in all_taxonomic_nodes:
            if node.path_name+'/'+node.node_name in target_node_paths_names:
                taxon_names2ids[taxonomic_level].append(node.node_name)
                taxon_names2ids["id"].append(idx)

                virus_dict = {"sequence": [], "label": []}
                subseqs = get_subseqs_from_final_node(
                    node, num_seqs_extraction)
                for subseq in subseqs:
                    print(subseq)
                    virus_dict["sequence"].append(subseq)
                    virus_dict["label"].append(idx)
                idx += 1
                virus_df = pd.DataFrame.from_dict(virus_dict)
                virus_df.to_csv(
                    output_path, mode='a',
                    header=not os.path.exists(output_path))
        taxon_df = pd.DataFrame.from_dict(taxon_names2ids)
        taxon_df.to_csv(taxon_names2ids_path, index=False)


write_csvs()
