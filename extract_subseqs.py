import math
import re
import os
import random
import subprocess
import sys
import pathlib
from progress.bar import Bar
from Bio import Align
from datetime import datetime
try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo
# import pandas as pd
import taxoniq
from bigtree import Node, add_path_to_tree, find_name, tree_to_newick, \
    newick_to_tree, find_attrs, find_attr, clone_tree, \
    copy_and_replace_nodes_from_tree_to_tree

from itertools import groupby


def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


KRAKEN_PATH = "/home/aschoier/workspace/kraken2"
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


def get_tax_id(seq_ref):
    with open(f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/seqid2taxid.map") as f:
        for line in f.readlines():
            # kraken:taxid|687377|NC_002195.1 687377
            if line.rstrip().split("|")[-1].split()[0] == seq_ref:
                tax_id = int(line.split("|")[-1].split()[-1])
                break
    f.close()

    return tax_id


def get_tax_ids():
    tax_ids = {}
    with open(f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/seqid2taxid.map") as f:
        for line in f.readlines():
            seq_ref = line.rstrip().split("|")[-1].split()[0]
            tax_id = int(line.rstrip().split("|")[-1].split()[-1])
            tax_ids[seq_ref] = tax_id
    f.close()

    return tax_ids


def get_num_genomes():
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    return len(
        [
            name
            for name in os.listdir(gen_dir)
            if os.path.isdir(os.path.join(gen_dir, name))
        ]
    )


def generate_seqs_by_taxon_tree():
    taxon_tree = Node("root")
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    tax_ids = get_tax_ids()

    file_tree = os.walk(gen_dir)
    file_tree.__next__()

    num_genomes = get_num_genomes()
    visited_nodes = set()
    registered_subseq_refs = set()
    bar = Bar("Loading taxon tree", max=num_genomes)
    #print(f"Loading taxon tree", flush=True)
    for subdir, _, _ in file_tree:
        seq_ref = subdir.split("/")[-1]
        tax_id = tax_ids[seq_ref]

        sequences = get_genome_sequences(seq_ref)
        t = taxoniq.Taxon(tax_id)
        ranked_taxons = t.ranked_lineage
        path_parent_rank = ""
        for idx, ranked_taxon in enumerate(reversed(ranked_taxons)):
            slash = "" if path_parent_rank == "" else "/"
            path_parent_rank += slash + ranked_taxon.scientific_name
            #print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
            #print(f"Loading { path_parent_rank }", flush=True)
            if path_parent_rank in visited_nodes:
                already_added = True
                for subseq_ref in sequences.keys():
                    if subseq_ref not in registered_subseq_refs:
                        already_added = False
                if idx == len(ranked_taxons) - 1 and not already_added:
                    add_path_to_tree(
                        taxon_tree,
                        f"{ path_parent_rank }/{ seq_ref }",
                        node_attrs={"rank": "genome"},
                    )
                    for subseq_ref, seq_list in sequences.items():
                        add_path_to_tree(
                            taxon_tree,
                            f"{ path_parent_rank }/{ seq_ref }/{ subseq_ref }",
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
                            f"{ path_parent_rank }/{ seq_ref }",
                            node_attrs={
                                "rank": "genome",
                            },
                        )
                        for subseq_ref, seq_list in sequences.items():
                            add_path_to_tree(
                                taxon_tree,
                                f"{ path_parent_rank }/{ seq_ref }/" +
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


def get_genome_size(sequences):
    total_len = 0
    for seq in sequences.values():
        total_len += len(seq)
    return total_len


def get_genome_sequences(seq_ref):
    sequences = {}
    with open(
        f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/{seq_ref}/genome.fna"
    ) as f:
        lines = f.readlines()
        subseq = None
        subseq_ref = ""
        subseq_name = ""
        for i in range(len(lines)):
            # >NC_019947.1 Tomato yellow mottle virus segment DNA-B, complete sequence
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


def get_sequences(n=0):
    gen_dir = f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes"
    num_genomes = len(
        [
            name
            for name in os.listdir(gen_dir)
            if os.path.isdir(os.path.join(gen_dir, name))
        ]
    )
    taxons = ["kingdom", "phylum", "order",
              "class", "family", "genus", "species"]
    sequences = {}
    # sample of n genomes
    if n != 0:
        genome_names = []
        sample = random.sample(range(0, num_genomes), n)
        first = True
        i = 0
        for subdir, _, _ in os.walk(gen_dir):
            if first:
                first = False
                continue
            if i in sample:
                genome_names.append(subdir)
            i += 1
    if n == 0:
        total_genomes = num_genomes
        print(f"Loading SeqRef sequences from { num_genomes } files...")
    else:
        total_genomes = n
        print(f"Loading SeqRef sequences from { n } of {num_genomes} files...")
    first = True
    num_load = 0

    tax_ids = get_tax_ids()

    for subdir, _, _ in os.walk(gen_dir):
        if first:  # skip /genomes
            first = False
            continue

        # sys.stdout.write(f'\033[KLoaded { num_load } of { total_genomes } SeqRef sequences\r')
        num_load += 1
        seq_ref = subdir.split("/")[-1]
        tax_id = tax_ids[seq_ref]

        t = taxoniq.Taxon(tax_id)
        ranked_taxons = t.ranked_lineage
        for ranked_taxon in reversed(ranked_taxons):
            if ranked_taxon.rank.name in taxons:
                print(ranked_taxon.rank, ranked_taxon.tax_id)
                ranked_tax_id = get_tax_id_taxoniq(ranked_taxon)
        os._exit(os.EX_OK)

        if n != 0 and subdir not in genome_names:
            continue
        # walk over fasta
        with open(subdir + "/genome.fna") as f:
            lines = f.readlines()
            subseq = None
            for i in range(len(lines)):
                if lines[i].startswith(">"):
                    # >NC_019947.1 Tomato yellow mottle virus segment DNA-B, complete sequence
                    if subseq:
                        for rank in sequences.keys():
                            if seq_ref not in rank_tax_ids[rank].keys():
                                continue
                            tax_id = rank_tax_ids[rank][seq_ref]
                            sequences[rank][tax_id][seq_ref][subseq_ref] = subseq
                    subseq_ref = lines[i].split()[0][1:]
                    subseq = ""
                elif i == len(lines) - 1:
                    subseq += lines[i].rstrip()
                    for rank in sequences.keys():
                        if seq_ref not in rank_tax_ids[rank].keys():
                            continue
                        tax_id = rank_tax_ids[rank][seq_ref]
                        sequences[rank][tax_id][seq_ref][subseq_ref] = subseq
                else:
                    subseq += lines[i].rstrip()
        f.close()
        if n != 0 and len(sequences) == n:
            break
    sys.stdout.write("\n")
    print("Done!")
    return sequences


def get_subseqs(rank, seqs, seq_len_spec={}, n=100, min_len=100, max_len=512):
    len_seqs = len(seqs)
    print(
        f"Generating { n } subsequences of length between { min_len }bp and { max_len }bp from { len_seqs } { rank }..."
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
        f"{ num_too_small } { rank } had sequences smaller than { max_len }bp and were left out of pool"
    )
    return subseqs, num_too_small


def write_csv(n_seqs=0, n_subseqs=100, min_subseq_len=100, max_subseq_len=512):
    virus_dict = {"sequence": [], "label": []}

    seq_len_spec = {}
    viral_seqs = get_sequences(n_seqs)
    viral_subseqs, num_too_small = get_subseqs(viral_seqs, seq_len_spec)

    if n_seqs != 0:
        while num_too_small > 0:
            viral_seqs_rest, num_too_small = get_sequences(num_too_small)
            for tax_id, seqs in viral_seqs_rest.items():
                viral_subseqs[tax_id] = seqs

    tax_ids = get_tax_ids()
    taxids2ids = {}
    ids2seqlen = {}
    idx = 1
    print("Preparing data to write csv...")
    for tax_id, subseqs_dict in viral_subseqs.items():
        sys.stdout.write(
            f"\033[KPrepared { idx } of { len(viral_subseqs) } classes\r")
        taxids2ids[tax_id] = idx
        ids2seqlen[idx] = seq_len_spec[tax_id]
        for part_ref_seq, seqs in subseqs_dict.items():
            for seq in seqs:
                virus_dict["sequence"].append(seq)
                virus_dict["label"].append(idx)
        idx += 1
    sys.stdout.write("\n")
    print(f"Writing { idx - 1 } classes to csv...")
    write_map_file(taxids2ids, "taxid2specid")
    write_map_file(ids2seqlen, "specid2seqlen")
    virus_df = pd.DataFrame.from_dict(virus_dict)
    virus_df.to_csv(
        f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/virus_seqs.csv", index=False)
    print("Done!")


def get_leaves_multiple_refseq(tree):
    all_leaves = []
    all_genomes = find_attrs(tree, "rank", "genome")
    for genome in all_genomes:
        all_leaves.append(genome.parent)
    all_leaves = set(all_leaves)
    leaves_mult_refseq = []
    for leaf in all_leaves:
        if len(find_attrs(leaf, "rank", "genome")) > 1:
            print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
            print(f"{ leaf.node_name } has more than one RefSeq", flush=True)
            leaves_mult_refseq.append(leaf)
    return leaves_mult_refseq


def get_sequences_leaf_mult_refseq(leaf):
    sequences = {}
    for genome in find_attrs(leaf, "rank", "genome"):
        sequence = ""
        for seq in find_attrs(genome, "rank", "sequence"):
            sequence += seq.get_attr("seq")
        sequences[genome.node_name] = sequence
    return sequences


def get_subsequences_leaf_mult_refseq(leaf):
    sequences = {}
    for refseq in find_attrs(leaf, "rank", "genome"):
        sequence = {}
        for subseq in find_attrs(refseq, "rank", "sequence"):
            sequence[subseq.node_name] = [
                subseq.get_attr("seq_name"), subseq.get_attr("seq")]
        sequences[refseq.node_name] = sequence
    return sequences


def write_leaf_mult_refseqs_fasta(leaf):
    if os.path.isfile(f"fasta/{ leaf.node_name }.fasta"):
        return
    sequences = get_sequences_leaf_mult_refseq(leaf)
    with open(f"fasta/{ leaf.node_name }.fasta", "w") as f:
        for ref_seq, sequence in sequences.items():
            f.write(f">{ ref_seq }\n{ sequence }\n")
    f.close()


def write_leaf_mult_refseqs_fastas(leaf):
    if not os.path.exists(f"fasta/{leaf.node_name}"):
        os.makedirs(f"fasta/{ leaf.node_name }")

    file_too_long = 0
    sequences = get_subsequences_leaf_mult_refseq(leaf)
    fasta_subseqs = {}
    for refseq, subseqs in sequences.items():
        for subseq_ref, subseq_list in subseqs.items():
            if subseq_list[0] not in fasta_subseqs.keys():
                fasta_subseqs[subseq_list[0]] = []
            fasta_subseqs[subseq_list[0]].append({refseq: subseq_list[1]})

    merge_subseqs = False
    keywords_merge = ["strain", "phage", "isolate", "clone", "complete", "partial"]
    compare_subseqref = None
    for subseq_name, subseqs in fasta_subseqs.items():
        if merge_subseqs:
            if compare_subseqref == list(subseqs[0].keys())[0]:
                merge_subseqs = False
        if (len(subseqs) == 1
                and any(ext.lower() in subseq_name for ext in keywords_merge)):
            compare_subseqref = list(subseqs[0].keys())[0]
            merge_subseqs = True

    if merge_subseqs:
        new_fasta_subseqs = {}
        new_fasta_subseqs[leaf.node_name] = []
        for refseq, subseqs in sequences.items():
            for subseq_ref, subseq_list in subseqs.items():
                new_fasta_subseqs[leaf.node_name].append(
                    {refseq: subseq_list[1]})
        fasta_subseqs = new_fasta_subseqs

    for subseq_name, subseqs in fasta_subseqs.items():
        new_subseq_name = subseq_name.replace('/', '-')
        print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
        print(f"writing fasta for { leaf.node_name }/{ new_subseq_name }", flush=True)
        if len(new_subseq_name) > 200:
            file_too_long += 1
            with open(f"fasta/{ leaf.node_name }/name_too_long{ file_too_long }", "w") as f:
                f.write(new_subseq_name)
            f.close()
            with open(f"fasta/{ leaf.node_name }/name_too_long{ file_too_long }.fasta", "w") as f:
                for subseq in subseqs:
                    for refseq, seq in subseq.items():
                        f.write(f">{ refseq }\n{ seq }\n")
            f.close()
        else:
            with open(f"fasta/{ leaf.node_name }/{ new_subseq_name }.fasta", "w") as f:
                for subseq in subseqs:
                    for refseq, seq in subseq.items():
                        f.write(f">{ refseq }\n{ seq }\n")
            f.close()


def get_leaves_multiple_differentiative(tree):
    differentiative_words = ['DNA', 'RNA', 'VP', 'NSP', 'VD', 'PB', '(PB', 'B_', 'II_', '(M.CviS',
        'contig', 'Contig', 'Circle', 'dsRNA', 'Maloyas', 'Diaz', 'Limbo', 'Poko', 'Forest', 'Juan']
    keywords_merge = ["strain", "phage", "isolate", "clone", "variant"]
    multiple_differentiative_leaves = []
    for leaf in find_attrs(tree, "rank", "genome"):
        if len(find_attrs(leaf, "rank", "sequence")) > 1:
            seq_names = {}
            for seq in find_attrs(leaf, "rank", "sequence"):
                seq_names[seq.node_name] = seq.get_attr("seq_name")
            checked_refs = set()
            diff_words = {}
            for seq_ref_a, seq_name_a in seq_names.items():
                for seq_ref_b, seq_name_b in seq_names.items():
                    if seq_ref_a in checked_refs and seq_ref_b in checked_refs:
                        continue
                    checked_refs.add(seq_ref_b)
                    diff = set(seq_name_a.split()) ^ set(seq_name_b.split())
                    for difference in diff:
                        differentiative = False
                        for word in differentiative_words:
                            if difference.startswith(word):
                                differentiative = True
                        if difference in seq_name_b.split():
                            prev_diff = seq_name_b.split()[seq_name_b.split().index(difference) - 1]
                        elif difference in seq_name_a.split():
                            prev_diff = seq_name_a.split()[seq_name_a.split().index(difference) - 1]
                        if prev_diff in keywords_merge and not differentiative:
                            if prev_diff in diff_words.keys():
                                diff_words[prev_diff].append(difference)
                            else:
                                diff_words[prev_diff] = [difference]
            if len(diff_words) == 1 and len(next(iter(diff_words.values()))) == len(seq_names):
                multiple_differentiative_leaves.append(leaf)
    return multiple_differentiative_leaves


def write_leaves_multiple_differentiative(tree):
    multiple_differentiative_leaves = get_leaves_multiple_differentiative(tree)
    for leaf in multiple_differentiative_leaves:
        write_leaf_multiple_differentiative(leaf)


def align_leaves_multiple_differentiative(tree):
    multiple_differentiative_leaves = get_leaves_multiple_differentiative(tree)
    for leaf in multiple_differentiative_leaves:
        align_sequences_multiple_differentiative(leaf)


def write_leaf_multiple_differentiative(leaf):
    fasta_dir = f"fasta/{ leaf.node_name }"
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)
    seqs = {}
    for seq in find_attrs(leaf, "rank", "sequence"):
        seqs[seq.node_name] = seq.get_attr("seq")
        with open(f"{ fasta_dir }/{ leaf.node_name }.fasta", "w") as f:
            for refseq, seq in seqs.items():
                f.write(f">{ refseq }\n{ seq }\n")
        f.close()


def align_sequences(leaf):
    fasta_file = f"fasta/{ leaf.node_name }.fasta"
    aln_file = f"aln/{ leaf.node_name }.aln"
    subprocess.run([
        "clustalo", "-i", fasta_file,
        "-o", aln_file,
        "--threads", "50",
        "--log", "clustalo.log",
        "--verbose",
        "--force"
    ])


def align_sequences_mult_refseqs(leaf):
    fasta_dir = f"fasta/{ leaf.node_name }"
    aln_dir = f"aln/{ leaf.node_name }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    for fasta_file in os.listdir(fasta_dir):
        aln_file = f"{ aln_dir }/{ fasta_file.replace('.fasta', '.aln') }"
        print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
        print(f"aligning { aln_dir }/{ fasta_file.replace('.fasta', '.aln') }", flush=True)
        subprocess.run([
            "clustalo", "-i", f"{ fasta_dir }/{ fasta_file }",
            "-o", aln_file,
            "--threads", "50",
            "--log", "clustalo.log",
            "--verbose",
            "--force"
        ])


def align_sequences_multiple_differentiative(leaf):
    fasta_dir = f"fasta/{ leaf.node_name }"
    aln_dir = f"aln/{ leaf.node_name }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    aln_file = f"{ aln_dir }/{ leaf.node_name }.aln"
    subprocess.run([
        "clustalo", "-i", f"{ fasta_dir }/{ leaf.node_name }.fasta",
        "-o", aln_file,
        "--threads", "50",
        "--log", "clustalo.log",
        "--verbose",
        "--force"
    ])


def align_all_leaves_mult_refseq(tree):
    leaves_mult_refseq = get_leaves_multiple_refseq(tree)
    bar = Bar("Align multiple RefSeq species", max=len(leaves_mult_refseq))
    for leaf in leaves_mult_refseq:
        write_leaf_mult_refseqs_fasta(leaf)
        align_sequences(leaf)
        bar.next()
    bar.finish()


def align_all_leaves_mult_refseqs(tree):
    leaves_mult_refseq = get_leaves_multiple_refseq(tree)
    bar = Bar("Align multiple RefSeq species", max=len(leaves_mult_refseq))
    for leaf in leaves_mult_refseq:
        write_leaf_mult_refseqs_fastas(leaf)
        align_sequences_mult_refseqs(leaf)
        bar.next()
    bar.finish()


def find_mutation_positions(reference_seq, seq):
    mutations = {}
    for i in range(len(reference_seq)):
        if reference_seq[i] != seq[i]:
            mutations[i] = seq[i]
    return mutations


def extract_mutations_from_aln(folder, filename):
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
    reference_refsq = list(subseqs.keys())[0]
    reference_seq = subseqs[reference_refsq]
    for seq_ref, seq in subseqs.items():
        if seq_ref == reference_refsq:
            continue
        if seq != reference_seq:
            seq_mutations = find_mutation_positions(reference_seq, seq)
            mutations[seq_ref] = seq_mutations
    return reference_seq, mutations


def add_mutations_to_tree(tree):
    new_tree = clone_tree(tree, Node)
    for leaf in get_leaves_multiple_refseq(tree):
        reference_seq, mutations = extract_mutations_from_aln(leaf.node_name)
        copy_and_replace_nodes_from_tree_to_tree(
            from_tree=tree,
            to_tree=new_tree,
            from_paths=[leaf.path_name],
            to_paths=[leaf.path_name],
            delete_children=True,
            with_full_path=True,
        )
        if not mutations:
            add_path_to_tree(
                new_tree,
                f"{ leaf.path_name }/genome",
                node_attrs={
                    "rank": "sequence",
                    "seq": reference_seq,
                },
            )
        else:
            add_path_to_tree(
                new_tree,
                f"{ leaf.path_name }/genome",
                node_attrs={
                    "rank": "sequence",
                    "seq": reference_seq,
                    "mutations": mutations
                },
            )
    return new_tree


def replace_subseqs(tree):
    new_tree = clone_tree(tree, Node)
    all_leaves = []
    all_genomes = find_attrs(tree, "rank", "genome")
    for genome in all_genomes:
        all_leaves.append(genome.parent)
    all_leaves = set(all_leaves)
    for leaf in all_leaves:
        sequences = get_sequences_leaf_mult_refseq(leaf)
        copy_and_replace_nodes_from_tree_to_tree(
            from_tree=tree,
            to_tree=new_tree,
            from_paths=[leaf.path_name],
            to_paths=[leaf.path_name],
            delete_children=True,
            with_full_path=True,
        )
        for seq_ref, seq in sequences.items():
            add_path_to_tree(
                new_tree,
                f"{ leaf.path_name }/genome",
                node_attrs={
                    "rank": "sequence",
                    "seq": seq,
                },
            )

    return new_tree


def get_multiple_subseqs_in_multiple_refseq(tree):
    leaves_mult_refseq = get_leaves_multiple_refseq(tree)
    multiple_subseqs = []
    for leaf in leaves_mult_refseq:
        genomes = find_attrs(leaf, "rank", "genome")
        for genome in genomes:
            if len(find_attrs(genome, "rank", "sequence")) > 1:
                multiple_subseqs.append(leaf)
                break
    return multiple_subseqs


tree = generate_seqs_by_taxon_tree()
align_leaves_multiple_differentiative(tree)
