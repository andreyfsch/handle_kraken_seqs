import math
import os
import random
import subprocess
import sys
import pathlib
from progress.bar import Bar

import pandas as pd
import taxoniq
from bigtree import Node, add_path_to_tree, find_name, tree_to_newick, \
    newick_to_tree, find_attrs

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


def get_taxon_rank(taxon_id):
    t = taxoniq.Taxon(taxon_id)
    return t.rank.name


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


def get_species_id(seq_ref):
    tax_id = get_tax_id(seq_ref)
    with open(f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/taxid2specid.map") as f:
        for line in f.readlines():
            if line.rstrip().split()[0] == tax_id:
                species_id = int(line.split()[-1])
                break
    f.close()

    return species_id


def write_map_file(mapdict, filename):
    with open(
        f"{ KRAKEN_PATH }/{ KRAKEN_DATABASE }/{ filename }.map", "w"
    ) as f:
        for key, val in mapdict.items():
            f.write(str(key) + " " + str(val) + "\n")
    f.close()


def create_taxon_structure(sequences, ranked_taxon_ids):
    pass


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

    bar = Bar("Loading taxon tree", max=num_genomes)
    for subdir, _, _ in file_tree:
        seq_ref = subdir.split("/")[-1]
        tax_id = tax_ids[seq_ref]

        visited_nodes = set()
        sequences = get_genome_sequences(seq_ref)
        t = taxoniq.Taxon(tax_id)
        ranked_taxons = t.ranked_lineage
        path_parent_rank = ""
        for idx, ranked_taxon in enumerate(reversed(ranked_taxons)):
            slash = "" if path_parent_rank == "" else "/"
            path_parent_rank += slash + ranked_taxon.scientific_name
            if path_parent_rank in visited_nodes:
                if idx == len(ranked_taxons) - 1:
                    add_path_to_tree(
                        taxon_tree,
                        f"{ path_parent_rank }/{ seq_ref }",
                        node_attrs={"rank": "genome"},
                    )
                    for subseq_ref, seq in sequences.items():
                        add_path_to_tree(
                            taxon_tree,
                            f"{ path_parent_rank }/{ seq_ref }/{ subseq_ref }",
                            node_attrs={"rank": "sequence", "seq": seq},
                        )
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
                    if idx == len(ranked_taxons) - 1:
                        add_path_to_tree(
                            taxon_tree,
                            f"{ path_parent_rank }/{ seq_ref }",
                            node_attrs={"rank": "genome"},
                        )
                        for subseq_ref, seq in sequences.items():
                            add_path_to_tree(
                                taxon_tree,
                                f"{ path_parent_rank }/{ seq_ref }/" +
                                f"{ subseq_ref }",
                                node_attrs={"rank": "sequence", "seq": seq},
                            )
        bar.next()
    bar.finish()
    return taxon_tree


def get_genome_sequences(seq_ref):
    sequences = {}
    with open(
        f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/{seq_ref}/genome.fna"
    ) as f:
        lines = f.readlines()
        subseq = None
        subseq_ref = ""
        for i in range(len(lines)):
            # >NC_019947.1 Tomato yellow mottle virus segment DNA-B, complete sequence
            if lines[i].startswith(">"):
                if subseq:
                    sequences[subseq_ref] = subseq
                subseq_ref = lines[i].split()[0][1:]
                subseq = ""
            elif i == len(lines) - 1:
                subseq += lines[i].rstrip()
                sequences[subseq_ref] = subseq
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


def seqs_by_taxon_tree_to_newick():
    tree = generate_seqs_by_taxon_tree()
    with open("seqs_by_taxon_tree.newick", "w") as f:
        f.write(tree_to_newick(tree, attr_list=["rank", "seq"]))
    f.close()


def newick_seqs_by_taxon_to_tree():
    with open("seqs_by_taxon_tree.newick", "r") as f:
        newick = f.read()
    f.close()
    return newick_to_tree(newick)


def get_leaves_multiple_refseq(tree):
    # tree = newick_seqs_by_taxon_to_tree()
    all_leaves = []
    all_genomes = find_attrs(tree, "rank", "genome")
    for genome in all_genomes:
        all_leaves.append(genome.parent)
    all_leaves = set(all_leaves)
    leaves_mult_refseq = []
    for leaf in all_leaves:
        if len(find_attrs(leaf, "rank", "genome")) > 1:
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


def write_leaf_mult_refseqs_fasta(leaf):
    if os.path.isfile(f"fasta/{ leaf.node_name }.fasta"):
        return
    sequences = get_sequences_leaf_mult_refseq(leaf)
    with open(f"fasta/{ leaf.node_name }.fasta", "w") as f:
        for ref_seq, sequence in sequences.items():
            f.write(f">{ ref_seq }\n{ sequence }\n")
    f.close()


def align_sequences(leaf):
    if os.path.isfile(f"aln/{ leaf.node_name }.aln"):
        return
    subprocess.run([
        "clustalo", "-i", f"fasta/{ leaf.node_name }.fasta",
        "-o", f"aln/{ leaf.node_name }.aln",
        "--threads", "40",
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
        # remove_fasta_files(leaf)
    bar.finish()


def remove_fasta_files(leaf):
    pathlib.Path.unlink(f"{ leaf.node_name }.fasta")
    # pathlib.Path.unlink(f"{ species.node_name }.aln")


def extract_subseqs_from_aln(filename):
    with open(f"{ filename }.aln") as f:
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
    has_gap = False
    for subseq_ref, seq in subseqs.items():
        if "-" in seq:
            has_gap = True
            break
    if not has_gap:
        return {subseq_ref: subseqs[subseq_ref]}
    else:
        subseq_gaps = {}
        for subseq_ref, seq in subseqs.items():
            for i in range(0, len(seq), 512):
                if "-" not in seq[i: i + 512]:
                    consensus = True
                    for subseq_ref_search, seq_search in subseqs.items():
                        if subseq_ref_search == subseq_ref:
                            continue
                        if (seq[i: i + 512] != seq_search[i: i + 512]):
                            consensus = False
                            break
                    if consensus:
                        subseq_gaps[f"consensus.{i}"] = seq[i: i + 512]
                    elif seq[i: i + 512] not in subseq_gaps.values():
                        subseq_gaps[f"{subseq_ref}.{i}.no-consensus"] = seq[i: i + 512]
                else:
                    for subseq_ref_search, seq_search in subseqs.items():
                        if subseq_ref_search == subseq_ref:
                            continue
                        if (
                                "-" not in seq_search[i: i + 512] and
                                seq_search[i: i + 512]
                                not in subseq_gaps.values()
                        ):
                            subseq_gaps[
                                f"{subseq_ref}.{i}.no-gap"
                            ] = seq_search[i: i + 512]
                            break
                        if seq_search[i: i + 512] not in subseq_gaps.values():
                            subseq_gaps[f"{subseq_ref}.{i}.no-gap_not-found"] = seq_search[i: i +
                                                                          512].replace("-", "")
        return subseq_gaps


tree = generate_seqs_by_taxon_tree()
align_all_leaves_mult_refseq(tree)

# def get_consensus_seq_species_mult_refseq(species):
#   sequences = get_sequences_species_mult_refseq(species)
#   TACCACAGGTTACGCTGAGTTATTTT
#   TACCACAGGTTACGCTGAGTTATTTT


# viral_seqs = get_sequences()

# seq_ref_lens = {}
# for taxid, seq_refs in viral_seqs["phylum"].items():
#     num_seq_refs = len(seq_refs)
#     total_len = 0
#     for subseq_refs in seq_refs.values():
#         for subseq in subseq_refs.values():
#             total_len += len(subseq)
#     print(taxid, num_seq_refs, total_len)


# for rank, taxids in viral_seqs.items():
# print(rank)
# print(len(taxids))
# biggest = 0
# biggest_taxid = None
# smallest = 99999999999
# smallest_taxid = None
# for taxid, seq_refs in taxids.items():
# total_len = sum(len(val) for val in [subseq for subseq in seq_refs.values()])
# if total_len > biggest:
# biggest = total_len
# biggest_taxid = taxid
# if total_len < smallest:
# smallest = total_len
# smallest_taxid = taxid
# print(biggest_taxid, biggest)
# print(smallest_taxid, smallest)

# write_csv()
