import math
import os
import random
import subprocess
import sys
import multiprocessing
from progress.bar import Bar
from typing import Dict, List, Set, Tuple
# from datetime import datetime
# import pandas as pd
# from zoneinfo import ZoneInfo
import taxoniq
from bigtree import Node, add_path_to_tree, \
    find_attrs, clone_tree, \
    copy_and_replace_nodes_from_tree_to_tree


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
        while len(subseqs) < n:
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


def get_tax_id_taxoniq(ranked_taxon):
    str_taxon = str(ranked_taxon)
    str_tax_id = str_taxon[str_taxon.find("(") + 1: str_taxon.find(")")]

    return int(str_tax_id)


def get_complement(seq):
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


def get_tax_id(ref_seq):
    with open(f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/seqid2taxid.map") as f:
        for line in f.readlines():
            # kraken:taxid|687377|NC_002195.1 687377
            if line.rstrip().split("|")[-1].split()[0] == ref_seq:
                tax_id = int(line.split("|")[-1].split()[-1])
                break
    f.close()

    return tax_id


def get_tax_ids():
    tax_ids = {}
    with open(f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/seqid2taxid.map") as f:
        for line in f.readlines():
            ref_seq = line.rstrip().split("|")[-1].split()[0]
            tax_id = int(line.rstrip().split("|")[-1].split()[-1])
            tax_ids[ref_seq] = tax_id
    f.close()

    return tax_ids


def get_num_refseq_files():
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
    containing ref_seq kraken2 files' sequences, branching them on
    subref_seqs leaves."""
    taxon_tree = Node("root")
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    tax_ids = get_tax_ids()

    file_tree = os.walk(gen_dir)
    file_tree.__next__()

    num_refseq_files = get_num_refseq_files()
    visited_nodes = set()
    registered_subseq_refs = set()
    bar = Bar("Loading taxon tree", max=num_refseq_files)
    # print(f"Loading taxon tree", flush=True)
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
            # print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
            # print(f"Loading { path_parent_rank }", flush=True)
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
                    for subseq_ref, seq_list in sequences.items():
                        add_path_to_tree(
                            taxon_tree,
                            f"{ path_parent_rank }/{ ref_seq }/{ subseq_ref }",
                            node_attrs={
                                "rank": "sequence",
                                "seq_name": seq_list[0],
                                "seq": seq_list[1]
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
                        for subseq_ref, seq_list in sequences.items():
                            add_path_to_tree(
                                taxon_tree,
                                f"{ path_parent_rank }/{ ref_seq }/" +
                                f"{ subseq_ref }",
                                node_attrs={
                                    "rank": "sequence",
                                    "seq": seq_list[1],
                                    "seq_name": seq_list[0]
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
        f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/{ref_seq}/genome.fna"
    ) as f:
        lines = f.readlines()
        subseq = None
        subseq_ref = ""
        subseq_name = ""
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                if subseq:
                    sequences[subseq_ref] = [subseq_name, subseq]
                line_split = lines[i].split()
                subseq_ref = line_split[0][1:]
                subseq_name = " ".join(line_split[1:])
                subseq = ""
            elif i == len(lines) - 1:
                subseq += lines[i].rstrip()
                sequences[subseq_ref] = [subseq_name, subseq]
            else:
                subseq += lines[i].rstrip()
    f.close()

    return sequences


def get_subseqs(rank, seqs, seq_len_spec={}, n=100, min_len=100, max_len=512):
    len_seqs = len(seqs)
    print(
        f"Generating { n } subsequences of length between \
        { min_len }bp and { max_len }bp from { len_seqs } { rank }..."
    )
    subseqs = {}
    num_load = 0
    num_too_small = 0
    for tax_id, part_seqs in seqs.items():
        total_len = len("".join([seq for seq in part_seqs.values()]))
        if total_len > max_len:
            sys.stdout.write(
                f"\033[KRead { num_load } of { len_seqs } species\r")
            seq_len_spec[tax_id] = total_len
            num_load += 1
            subseqs[tax_id] = {}
            n_seqs = {}
            # populate number of sequences extracted from each subsequence
            for part_seq_ref, seq in part_seqs.items():
                seq_fraction = len(seq) / total_len
                n_seqs[part_seq_ref] = int(round(seq_fraction * n, 0))
            # compensate rounding errors
            if sum(n_seqs.values()) != n:
                while sum(n_seqs.values()) != n:
                    if sum(n_seqs.values()) < n:
                        n_seqs[min(n_seqs, key=n_seqs.get)] += 1
                    elif sum(n_seqs.values()) > 100:
                        n_seqs[max(n_seqs, key=n_seqs.get)] -= 1
            for part_seq_ref, seq in part_seqs.items():
                # include smallest subsequences if size >= min_len
                if n_seqs[part_seq_ref] == 0:
                    if len(seq) >= min_len:
                        n_seqs[part_seq_ref] = 1
                        n_seqs[max(n_seqs, key=n_seqs.get)] -= 1
                    else:
                        # discard subsequences that are too small
                        continue
                part_subseqs = extract_subseqs(
                    seq, n_seqs[part_seq_ref], min_len, max_len
                )
                subseqs[tax_id][part_seq_ref] = part_subseqs
        else:
            num_too_small += 1
    sys.stdout.write("\n")
    print("Done!")
    print(
        f"{ num_too_small } { rank } had sequences smaller \
        than { max_len }bp and were left out of pool"
    )
    return subseqs, num_too_small


# def write_csv(n_seqs=0, n_subseqs=100,
#   min_subseq_len=100, max_subseq_len=512):
#     virus_dict = {"sequence": [], "label": []}

#     seq_len_spec = {}
#     viral_seqs = get_sequences(n_seqs)
#     viral_subseqs, num_too_small = get_subseqs(viral_seqs, seq_len_spec)

#     if n_seqs != 0:
#         while num_too_small > 0:
#             viral_seqs_rest, num_too_small = get_sequences(num_too_small)
#             for tax_id, seqs in viral_seqs_rest.items():
#                 viral_subseqs[tax_id] = seqs

#     taxids2ids = {}
#     ids2seqlen = {}
#     idx = 1
#     print("Preparing data to write csv...")
#     for tax_id, subseqs_dict in viral_subseqs.items():
#         sys.stdout.write(
#             f"\033[KPrepared { idx } of { len(viral_subseqs) } classes\r")
#         taxids2ids[tax_id] = idx
#         ids2seqlen[idx] = seq_len_spec[tax_id]
#         for part_ref_seq, seqs in subseqs_dict.items():
#             for seq in seqs:
#                 virus_dict["sequence"].append(seq)
#                 virus_dict["label"].append(idx)
#         idx += 1
#     sys.stdout.write("\n")
#     print(f"Writing { idx - 1 } classes to csv...")
#     virus_df = pd.DataFrame.from_dict(virus_dict)
#     virus_df.to_csv(
#         f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/virus_seqs.csv", index=False)
#     print("Done!")


def get_final_taxon_nodes_multiple_refseq(tree: Node) -> List:
    """Return list of final taxon nodes from a taxon tree containing
    more than one ref_seq."""
    all_final_taxon_nodes = []
    all_ref_seq_nodes = find_attrs(tree, "rank", "ref_seq")
    for ref_seq_node in all_ref_seq_nodes:
        all_final_taxon_nodes.append(ref_seq_node.parent)
    all_final_taxon_nodes = set(all_final_taxon_nodes)
    final_taxon_nodes_mult_refseq = []
    for leaf in all_final_taxon_nodes:
        if len(find_attrs(leaf, "rank", "ref_seq")) > 1:
            final_taxon_nodes_mult_refseq.append(leaf)
    return final_taxon_nodes_mult_refseq


def final_taxon_node_to_dict_subseqs(final_taxon_node: Node) -> Dict:
    """Return the dict: {
    ref_seq:{subseq_ref_seq:
        {seq_name: subseq_name, seq: subseq}, ...}, ...}
    of a given final taxon node.
    """
    sequences = {}
    for ref_seq_node in find_attrs(final_taxon_node, "rank", "ref_seq"):
        ref_seq = ref_seq_node.node_name
        subseq_dict = {}
        for subseq_node in find_attrs(ref_seq_node, "rank", "sequence"):
            subseq_ref_seq = subseq_node.node_name
            subseq_dict[subseq_ref_seq] = {
                "seq_name": subseq_node.get_attr("seq_name"),
                "seq": subseq_node.get_attr("seq")}
        sequences[ref_seq] = subseq_dict
    return sequences


def get_dict_fasta_subseqs(final_taxon_node: Node) -> Dict:
    """Return a dict grouping multiple ref_seq sequences in a final taxon node
    by subseq name.
    """
    subseqs_dict = final_taxon_node_to_dict_subseqs(final_taxon_node)
    fasta_subseqs = {}
    for ref_seq, subseqs in subseqs_dict.items():
        for subseq_val in subseqs.values():
            if subseq_val["seq_name"] not in fasta_subseqs.keys():
                fasta_subseqs[subseq_val["seq_name"]] = []
            fasta_subseqs[subseq_val["seq_name"]].append(
                {ref_seq: subseq_val["seq"]})
    return fasta_subseqs


def handle_merge_mult_refseqs(final_taxon_node: Node) -> Dict:
    """Return a dict merging groups of multiple ref_seq sequences
    in a final taxon node that have only one sequence and contain
    keywords in their name.
    """
    final_taxon_node_name = final_taxon_node.node_name
    subseqs_dict = final_taxon_node_to_dict_subseqs(final_taxon_node)
    fasta_subseqs = get_dict_fasta_subseqs(final_taxon_node)

    merge_subseqs = False
    keywords_merge = ["strain", "phage",
                      "isolate", "clone", "complete", "partial"]
    compare_subseq_ref = None
    for subseq_name, subseqs in fasta_subseqs.items():
        for subseq in subseqs:
            ref_seq = next(iter(subseq))
            if merge_subseqs:
                if compare_subseq_ref == ref_seq:
                    merge_subseqs = False
            if (len(subseqs) == 1
                    and any(keyword in subseq_name
                            for keyword in keywords_merge)):
                compare_subseq_ref = ref_seq
                merge_subseqs = True

    if merge_subseqs:
        new_fasta_subseqs = {}
        new_fasta_subseqs[final_taxon_node_name] = []
        for ref_seq, subseqs in subseqs_dict.items():
            for subseq_dict in subseqs.values():
                new_fasta_subseqs[final_taxon_node_name].append(
                    {ref_seq: subseq_dict["seq"]})
        fasta_subseqs = new_fasta_subseqs
    return fasta_subseqs


def write_final_taxon_node_mult_refseqs_fastas(final_taxon_node: Node) -> None:
    """Write fasta files for multiple ref_seq final taxon node."""
    final_taxon_node_name = final_taxon_node.node_name
    if not os.path.exists(f"fasta/{ final_taxon_node_name }"):
        os.makedirs(f"fasta/{ final_taxon_node_name }")

    fasta_subseqs = handle_merge_mult_refseqs(final_taxon_node)

    name_too_long = 0
    for subseq_name, subseqs in fasta_subseqs.items():
        new_subseq_name = subseq_name.replace('/', '-')
        if len(new_subseq_name) > 200:
            name_too_long += 1
            with open(
                f"fasta/{ final_taxon_node_name }/name_too_long"
                f"{ name_too_long }",
                "w"
            ) as f:
                f.write(new_subseq_name)
            f.close()
            with open(
                f"fasta/{ final_taxon_node_name }/name_too_long"
                f"{ name_too_long }.fasta",
                "w"
            ) as f:
                for subseq in subseqs:
                    for ref_seq, seq in subseq.items():
                        f.write(f">{ ref_seq }\n{ seq }\n")
            f.close()
        else:
            with open(
                f"fasta/{ final_taxon_node_name }/{ new_subseq_name }.fasta",
                "w"
            ) as f:
                for subseq in subseqs:
                    for ref_seq, seq in subseq.items():
                        f.write(f">{ ref_seq }\n{ seq }\n")
            f.close()


def get_refseqs_multiple_differentiative(tree: Node) -> Dict:
    """Return dict grouping sequences leaves that belong to the same ref_seq
    in a taxon tree but are actually variants of the same taxon by keyword."""
    differentiative_words = ['DNA', 'RNA', 'VP', 'NSP', 'VD', 'PB', '(PB',
                             'B_', 'II_', '(M.CviS', 'contig', 'Contig',
                             'Circle', 'dsRNA', 'Maloyas', 'Diaz', 'Limbo',
                             'Poko', 'Forest', 'Juan']
    keywords_merge = ["strain", "phage", "isolate", "clone", "variant"]
    multiple_differentiative_leaves = {}
    for leaf in find_attrs(tree, "rank", "ref_seq"):
        if len(find_attrs(leaf, "rank", "sequence")) > 1:
            seq_names = {}
            for sequence_node in find_attrs(leaf, "rank", "sequence"):
                seq_names[sequence_node.node_name] = sequence_node.get_attr("seq_name")
            checked_refs = []
            for seq_ref_a, seq_name_a in seq_names.items():
                for seq_ref_b, seq_name_b in seq_names.items():
                    if {seq_ref_a, seq_ref_b} in checked_refs:
                        continue
                    checked_refs.append({seq_ref_a, seq_ref_b})
                    diff = set(seq_name_a.split()) ^ set(seq_name_b.split())
                    if len(diff) > 2:
                        continue
                    for difference in diff:
                        differentiative = False
                        for word in differentiative_words:
                            if difference.startswith(word):
                                differentiative = True
                        if difference in seq_name_b.split():
                            prev_diff = seq_name_b.split(
                            )[seq_name_b.split().index(difference) - 1]
                            if prev_diff not in seq_name_a:
                                continue
                        elif difference in seq_name_a.split():
                            prev_diff = seq_name_a.split(
                            )[seq_name_a.split().index(difference) - 1]
                            if prev_diff not in seq_name_b:
                                continue
                        if prev_diff in keywords_merge and not differentiative:
                            if leaf.node_name not in \
                                    multiple_differentiative_leaves.keys():
                                multiple_differentiative_leaves[
                                    leaf.node_name] = {}

                            for sequence_node in find_attrs(leaf,
                                                            "rank", "sequence"):
                                if sequence_node.node_name in [seq_ref_a,
                                                               seq_ref_b]:
                                    add_set_to_dict(
                                        multiple_differentiative_leaves[
                                            leaf.node_name],
                                        prev_diff, sequence_node)
    return multiple_differentiative_leaves


def add_set_to_dict(dictionary, key, value):
    if key not in dictionary.keys():
        dictionary[key] = {value}
    else:
        dictionary[key].add(value)


def write_fasta_final_taxon_nodes_multiple_differentiative(tree: Node):
    """Write fasta files of ref_seq nodes in a taxon tree containing variants
    of the same taxon."""
    multiple_differentiative_leaves = get_refseqs_multiple_differentiative(
        tree)
    for ref_seq, diffs in multiple_differentiative_leaves.items():
        for prev_diff, seqs in diffs.items():
            write_final_taxon_node_multiple_differentiative(
                ref_seq, prev_diff, seqs)


def align_ref_seq_nodes_multiple_differentiative(
        tree: Node, num_cores: int = None) -> None:
    """Align fasta files of ref_seq nodes in a taxon tree containing variants
    of the same taxon using num_cores processor cpu cores. All available
    cpu cores are selected by default."""
    multiple_differentiative_leaves = get_refseqs_multiple_differentiative(
        tree)
    for ref_seq, diffs in multiple_differentiative_leaves.items():
        for prev_diff in diffs.keys():
            if num_cores:
                align_sequences_multiple_differentiative(ref_seq, prev_diff,
                                                         num_cores)
            else:
                align_sequences_multiple_differentiative(ref_seq, prev_diff)


def write_and_align_multiple_differentiative(
        tree: Node, num_cores: int = None) -> None:
    """Write fasta files of ref_seq nodes in a taxon tree containing variants
    of the same taxon and and align them using num_cores cpu processor cores.
    All available cpu cores are selected by default."""
    write_fasta_final_taxon_nodes_multiple_differentiative(tree)
    align_ref_seq_nodes_multiple_differentiative(tree, num_cores)


def write_final_taxon_node_multiple_differentiative(
        ref_seq: str, seqs: List) -> None:
    """Write fasta file containing seqs of ref_seq node in a taxon tree"""
    fasta_dir = f"fasta/{ ref_seq }"
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)
    with open(f"{ fasta_dir }/{ ref_seq }.fasta", "w") as f:
        for seq in seqs:
            f.write(f">{ seq.node_name }\n{ seq.get_attr('seq') }\n")
    f.close()


def align_sequences_mult_refseqs(
        final_taxon_node: Node,
        num_cores: int = multiprocessing.cpu_count()) -> None:
    """Align fasta files of final_taxon_node in a taxon tree containing
    more than one ref_seq using num_cores processor cpu cores. All available
    cpu cores are selected by default."""
    fasta_dir = f"fasta/{ final_taxon_node.node_name }"
    aln_dir = f"aln/{ final_taxon_node.node_name }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    for fasta_file in os.listdir(fasta_dir):
        aln_file = f"{ aln_dir }/{ fasta_file.replace('.fasta', '.aln') }"
        subprocess.run([
            "clustalo", "-i", f"{ fasta_dir }/{ fasta_file }",
            "-o", aln_file,
            "--threads", f"{ num_cores }",
            "--log", "clustalo.log",
            "--verbose",
            "--force"
        ])


def align_sequences_multiple_differentiative(
        ref_seq: str,
        prev_diff: str,
        num_cores: int = multiprocessing.cpu_count()
        ) -> None:
    """Align fasta file of ref_seq node in a taxon tree variants
    of the same taxon using num_cores processor cpu cores. All available
    cpu cores are selected by default."""
    fasta_dir = f"fasta/{ ref_seq }/{ prev_diff}"
    aln_dir = f"aln/{ ref_seq }/{ prev_diff }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    aln_file = f"{ aln_dir }/{ ref_seq }/{ prev_diff }.aln"
    subprocess.run([
        "clustalo", "-i", f"{ fasta_dir }/{ ref_seq }/{ prev_diff }.fasta",
        "-o", aln_file,
        "--threads", f"{ num_cores }",
        "--log", "clustalo.log",
        "--verbose",
        "--force"
    ])


def align_all_final_final_taxon_nodes_mult_refseqs(
        tree: Node, num_cores: int = None) -> None:
    """Write fasta files of final_taxon_nodes in a taxon tree containing
    more than one ref_seq and align them using num_cores processor cpu cores.
    All available cpu cores are selected by default."""
    final_taxon_nodes_mult_refseq = get_final_taxon_nodes_multiple_refseq(tree)
    bar = Bar("Align multiple RefSeq final taxon nodes",
              max=len(final_taxon_nodes_mult_refseq))
    for final_taxon_node in final_taxon_nodes_mult_refseq:
        write_final_taxon_node_mult_refseqs_fastas(final_taxon_node)
        if num_cores:
            align_sequences_mult_refseqs(final_taxon_node, num_cores)
        else:
            align_sequences_mult_refseqs(final_taxon_node)
        bar.next()
    bar.finish()


def find_mutation_positions(
        reference_sequence: str, compare_sequence: str) -> Dict:
    """Return dict containing mutation position as key
    and mutation as value."""
    mutations = {}
    for i in range(len(reference_sequence)):
        if reference_sequence[i] != compare_sequence[i]:
            mutations[i] = compare_sequence[i]

    point_indel_mutations = {}
    indel_start = None
    for pos, mut in mutations.items():
        if (pos - 1 not in mutations.keys() and
                pos + 1 not in mutations.keys()):
            point_indel_mutations[pos] = mut
        else:
            if pos - 1 not in mutations.keys():
                indel_start = pos
            if indel_start:
                point_indel_mutations[indel_start] += mut
            if pos + 1 not in mutations.keys():
                indel_start = None

    return point_indel_mutations


def extract_mutations_from_aln(folder: str, filename: str) -> Tuple[str, Dict]:
    """Return reference sequence and a mutation dict containing ref_seq as key
    and mutations dict as value."""
    with open(f"aln/{ folder }/{ filename }.aln") as f:
        lines = f.readlines()
        subseqs = {}
        subseq = ""
        subseq_ref = ""
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                if subseq:
                    subseqs[subseq_ref] = subseq
                subseq_ref = lines[i].split()[0][1:]
                subseq = ""
            elif i == len(lines) - 1:
                subseq += lines[i].rstrip()
                subseqs[subseq_ref] = subseq
            else:
                subseq += lines[i].rstrip()
    f.close()

    mutations = {}
    min_dash = 999999999
    for ref_seq, seq in subseqs.items():
        if seq.count("-") < min_dash:
            min_dash = seq.count("-")
            reference_refsq = ref_seq
            reference_sequence = subseqs[reference_refsq]
    for ref_seq, compare_sequence in subseqs.items():
        if ref_seq == reference_refsq:
            continue
        if compare_sequence != reference_sequence:
            seq_mutations = find_mutation_positions(
                reference_sequence, compare_sequence)
            mutations[ref_seq] = seq_mutations
    return reference_sequence, mutations


def add_mutations_multiple_refseq_to_tree(tree: Node) -> Node:
    new_tree = clone_tree(tree, Node)
    """Return taxon tree substituting """
    for final_taxon_node in get_final_taxon_nodes_multiple_refseq(tree):
        add_path_to_tree(
            new_tree,
            f"{ final_taxon_node.path_name }/ref_seq",
            node_attrs={
                "rank": "ref_seq",
            },
        )
        for aln_file in [
            os.path.splitext(f)[0]
            for f in os.listdir(f"aln/{ final_taxon_node.node_name }")
            if os.path.isfile(f)
        ]:
            reference_sequence, mutations = extract_mutations_from_aln(
                final_taxon_node.node_name, aln_file)
            copy_and_replace_nodes_from_tree_to_tree(
                from_tree=tree,
                to_tree=new_tree,
                from_paths=[final_taxon_node.path_name],
                to_paths=[final_taxon_node.path_name],
                delete_children=True,
                with_full_path=True,
            )
            if not mutations:
                add_path_to_tree(
                    new_tree,
                    f"{ final_taxon_node.path_name }/ref_seq/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_sequence
                    },
                )
            else:
                add_path_to_tree(
                    new_tree,
                    f"{ final_taxon_node.path_name }/ref_seq/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_sequence,
                        "mutations": mutations
                    },
                )
    return new_tree


def add_mutations_multiple_differentiative_to_tree(tree: Node) -> Node:
    new_tree = clone_tree(tree, Node)
    all_ref_seq_nodes = find_attrs(tree, "rank", "ref_seq")
    for ref_seq in get_refseqs_multiple_differentiative(tree).keys():
        for ref_seq_node in all_ref_seq_nodes:
            if ref_seq_node.node_name == ref_seq:
                target_node = ref_seq_node
                parent = ref_seq_node.parent
        add_path_to_tree(
            new_tree,
            f"{ parent.path_name }/{ ref_seq }",
            node_attrs={
                "rank": "ref_seq",
            },
        )
        for aln_file in [
            os.path.splitext(f)[0]
            for f in os.listdir(f"aln/{ ref_seq }")
            if os.path.isfile(f)
        ]:
            reference_sequence, mutations = extract_mutations_from_aln(
                ref_seq, aln_file)
            copy_and_replace_nodes_from_tree_to_tree(
                from_tree=tree,
                to_tree=new_tree,
                from_paths=[target_node.path_name],
                to_paths=[target_node.path_name],
                delete_children=True,
                with_full_path=True,
            )
            if not mutations:
                add_path_to_tree(
                    new_tree,
                    f"{ target_node.path_name }/ref_seq/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_sequence
                    },
                )
            else:
                add_path_to_tree(
                    new_tree,
                    f"{ target_node.path_name }/ref_seq/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_sequence,
                        "mutations": mutations
                    },
                )
    return new_tree


def extract_max_subseqs(final_taxon_node: Node, window_size: int) -> Set:
    max_subseqs = set()
    seqs = find_attrs(final_taxon_node, "rank", "sequence")

    for sequence_node in seqs:
        sequence = sequence_node.get_attr("seq")
        mutations = sequence_node.get_attr("mutations")
        if mutations:
            add_seqs = get_mutation_subseqs(sequence, window_size, mutations)
            for seq in add_seqs:
                max_subseqs.add(seq)
        else:
            for i in range(0, len(sequence) - window_size):
                max_subseqs.add(sequence[i:i+window_size])

    return max_subseqs


def get_mutation_subseqs(sequence, window_size, mutations):
    add_seqs = []
    for i in range(0, len(sequence) - window_size - 1):
        add_subseq = sequence[i:i+window_size-1]
        window_end = i + window_size - 1
        if "-" in add_subseq:
            add_subseq = get_seq_next_to_window_subseq_gaps(
                sequence, add_subseq, window_size, window_end)
        if add_subseq:
            add_seqs.append(add_subseq)
        new_seqs = get_sequences_mutations(
            sequence, i, window_size, mutations)
        for seq in new_seqs:
            add_seqs.append(seq)
    return add_seqs


def get_sequences_mutations(sequence, window_start, window_size, mutations):
    new_seqs = []
    add_subseq = sequence[window_start:window_start+window_size-1]
    mutations_in_range = get_mutations_in_range(
        window_start, window_size, mutations)
    if num_gaps_empty(mutations_in_range):
        for seq in get_subseq_mutations(add_subseq,
                                        mutations_in_range,
                                        window_size):
            new_seqs.append(seq)
    else:
        for seq in get_mutation_seqs_gaps(sequence, window_start,
                                          window_size,
                                          mutations):
            new_seqs.append(seq)
    return new_seqs


def num_gaps_empty(mutations_in_range):
    num_gaps = {}
    for ref_seq, muts in mutations_in_range.items():
        num_gaps[ref_seq] = 0
        for mut in muts:
            mut_start = next(iter(mut))
            if "-" in mut[mut_start]:
                num_gaps[ref_seq] += mut[mut_start].count("-")
    for gaps in num_gaps.values():
        if gaps > 0:
            return False
    return True


def get_mutation_seqs_gaps(
        sequence, window_start, window_size, mutations):
    mutations_in_range = get_mutations_in_range(
        window_start, window_size, mutations)
    new_seqs = []
    subseq = sequence[window_start:window_start+window_size-1]
    for ref_seq, muts in mutations_in_range.items():
        for mut_start, mut in muts.items():
            new_subseq = subseq
            new_subseq[mut_start:mut_start+len(mut)-1] = mut
            new_subseq = new_subseq.replace("-", "")
        if not new_subseq:
            continue
        elif len(new_subseq) < window_size:
            window_end = window_start + window_size - 1
            seq_next_gap = get_seq_next_to_window_mut_gaps(
                sequence, new_subseq, window_size,
                mutations, window_end, ref_seq)
            if seq_next_gap:
                new_seqs.append(seq_next_gap)
        else:
            new_seqs.append(new_subseq)
    return new_seqs


def get_seq_next_to_window_mut_gaps(
        sequence, subseq, window_size, mutations, window_end, ref_seq):
    diff = window_size - len(subseq)
    new_window_end = window_end + diff - 1
    if new_window_end > len(sequence) - 1:
        return None
    seq_diff = sequence[window_end:window_end+diff-1]
    new_subseq = subseq + seq_diff
    for mut_start, mut in mutations[ref_seq].items():
        new_subseq = subseq
        new_subseq[mut_start:mut_start+len(mut)-1] = mut
        new_subseq = new_subseq.replace("-", "")
    if len(new_subseq) < window_size:
        return get_seq_next_to_window_mut_gaps(
            sequence, new_subseq, window_size, mutations,
            new_window_end, ref_seq)
    else:
        return new_subseq


def get_seq_next_to_window_subseq_gaps(
        sequence, subseq, window_size, window_end):
    new_subseq = subseq.replace("-", "")
    diff = window_size - len(new_subseq)
    new_window_end = window_end + diff - 1
    if new_window_end > len(sequence) - 1:
        return None
    seq_diff = sequence[window_end:window_end+diff-1]
    new_subseq = new_subseq + seq_diff
    new_subseq = subseq.replace("-", "")
    if len(new_subseq) < window_size:
        return get_seq_next_to_window_subseq_gaps(
            sequence, new_subseq, window_size, new_window_end)
    else:
        return new_subseq


def get_subseq_mutations(subseq, mutations_in_range, window_size):
    new_seqs = []
    for muts in mutations_in_range.values():
        new_seq = subseq
        for mut in muts:
            mut_start = next(iter(mut))
            mut_len = len(mut[mut_start])
            mut_end = mut_start + mut_len - 1
            if mut_end > window_size - 1:
                mut_end = window_size - 1
            new_seq[mut_start:mut_end] = mut[mut_start]
        new_seqs.append(new_seq)


def get_mutations_in_range(window_start, window_size, mutations):
    mutations_in_range = {}
    for ref_seq, muts in mutations.items():
        for mut_start, mut in muts.items():
            if mut_start in range(window_start, window_start+window_size) \
                    or mut_start + len(mut) - 1 \
                    in range(window_start, window_start+window_size):
                if ref_seq in mutations_in_range.keys():
                    mutations_in_range[ref_seq] = \
                        {mut_start-window_start: mut}
                else:
                    mutations_in_range[ref_seq][
                        mut_start - window_start] = mut
    return mutations_in_range


tree = generate_seqs_by_taxon_tree()
