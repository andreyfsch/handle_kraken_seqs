import math
import os
import random
import subprocess
import sys
from progress.bar import Bar
from datetime import datetime
# import pandas as pd
from zoneinfo import ZoneInfo
import taxoniq
from bigtree import Node, add_path_to_tree, \
    find_attrs, clone_tree, \
    copy_and_replace_nodes_from_tree_to_tree
import collections


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


def get_num_refseq_files():
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

    num_refseq_files = get_num_refseq_files()
    visited_nodes = set()
    registered_subseq_refs = set()
    bar = Bar("Loading taxon tree", max=num_refseq_files)
    # print(f"Loading taxon tree", flush=True)
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
    # >NC_019947.1 Tomato yellow mottle virus segment DNA-B, complete sequence
    with open(
        f"{KRAKEN_PATH}/{KRAKEN_DATABASE}/genomes/{seq_ref}/genome.fna"
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


def get_subsequences_leaf_mult_refseq(leaf):
    sequences = {}
    for refseq in find_attrs(leaf, "rank", "genome"):
        sequence = {}
        for subseq in find_attrs(refseq, "rank", "sequence"):
            sequence[subseq.node_name] = [
                subseq.get_attr("seq_name"), subseq.get_attr("seq")]
        sequences[refseq.node_name] = sequence
    return sequences


def write_leaf_mult_refseqs_fastas(leaf):
    leaf_name = leaf.node_name
    if not os.path.exists(f"fasta/{ leaf_name }"):
        os.makedirs(f"fasta/{ leaf_name }")

    file_too_long = 0
    sequences = get_subsequences_leaf_mult_refseq(leaf)
    fasta_subseqs = {}
    for refseq, subseqs in sequences.items():
        for subseq_ref, subseq_list in subseqs.items():
            if subseq_list[0] not in fasta_subseqs.keys():
                fasta_subseqs[subseq_list[0]] = []
            fasta_subseqs[subseq_list[0]].append({refseq: subseq_list[1]})

    merge_subseqs = False
    keywords_merge = ["strain", "phage",
                      "isolate", "clone", "complete", "partial"]
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
        new_fasta_subseqs[leaf_name] = []
        for refseq, subseqs in sequences.items():
            for subseq_ref, subseq_list in subseqs.items():
                new_fasta_subseqs[leaf_name].append(
                    {refseq: subseq_list[1]})
        fasta_subseqs = new_fasta_subseqs

    for subseq_name, subseqs in fasta_subseqs.items():
        new_subseq_name = subseq_name.replace('/', '-')
        print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
        print(
            f"writing fasta for { leaf_name }/{ new_subseq_name }",
            flush=True
        )
        if len(new_subseq_name) > 200:
            file_too_long += 1
            with open(
                f"fasta/{ leaf_name }/name_too_long{ file_too_long }",
                "w"
            ) as f:
                f.write(new_subseq_name)
            f.close()
            with open(
                f"fasta/{ leaf_name }/name_too_long{ file_too_long }.fasta",
                "w"
            ) as f:
                for subseq in subseqs:
                    for refseq, seq in subseq.items():
                        f.write(f">{ refseq }\n{ seq }\n")
            f.close()
        else:
            with open(
                f"fasta/{ leaf_name }/{ new_subseq_name }.fasta",
                "w"
            ) as f:
                for subseq in subseqs:
                    for refseq, seq in subseq.items():
                        f.write(f">{ refseq }\n{ seq }\n")
            f.close()


def get_leaves_multiple_differentiative(tree):
    differentiative_words = ['DNA', 'RNA', 'VP', 'NSP', 'VD', 'PB', '(PB',
                             'B_', 'II_', '(M.CviS', 'contig', 'Contig',
                             'Circle', 'dsRNA', 'Maloyas', 'Diaz', 'Limbo',
                             'Poko', 'Forest', 'Juan']
    keywords_merge = ["strain", "phage", "isolate", "clone", "variant"]
    multiple_differentiative_leaves = {}
    for leaf in find_attrs(tree, "rank", "genome"):
        if len(find_attrs(leaf, "rank", "sequence")) > 1:
            seq_names = {}
            for seq in find_attrs(leaf, "rank", "sequence"):
                seq_names[seq.node_name] = seq.get_attr("seq_name")
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
                            for seq in find_attrs(leaf, "rank", "sequence"):
                                if seq.node_name == seq_ref_a:
                                    add_set_to_dict(
                                        multiple_differentiative_leaves,
                                        leaf.node_name, seq)
                                elif seq.node_name == seq_ref_b:
                                    add_set_to_dict(
                                        multiple_differentiative_leaves,
                                        leaf.node_name, seq)
    return multiple_differentiative_leaves


def add_set_to_dict(dictionary, key, value):
    if key not in dictionary.keys():
        dictionary[key] = {value}
    else:
        dictionary[key].add(value)


def write_leaves_multiple_differentiative(tree):
    multiple_differentiative_leaves = get_leaves_multiple_differentiative(tree)
    for ref_seq, seqs in multiple_differentiative_leaves.items():
        write_leaf_multiple_differentiative(ref_seq, seqs)


def align_leaves_multiple_differentiative(tree):
    multiple_differentiative_leaves = get_leaves_multiple_differentiative(tree)
    for ref_seq in multiple_differentiative_leaves.keys():
        align_sequences_multiple_differentiative(ref_seq)


def write_and_align_multiple_differentiative(tree):
    write_leaves_multiple_differentiative(tree)
    align_leaves_multiple_differentiative(tree)


def write_leaf_multiple_differentiative(ref_seq, seqs):
    fasta_dir = f"fasta/{ ref_seq }"
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)
    print(f"writing fasta for { fasta_dir }/{ ref_seq }")
    with open(f"{ fasta_dir }/{ ref_seq }.fasta", "w") as f:
        for seq in seqs:
            f.write(f">{ seq.node_name }\n{ seq.get_attr('seq') }\n")
    f.close()


def align_sequences_mult_refseqs(leaf):
    fasta_dir = f"fasta/{ leaf.node_name }"
    aln_dir = f"aln/{ leaf.node_name }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    for fasta_file in os.listdir(fasta_dir):
        aln_file = f"{ aln_dir }/{ fasta_file.replace('.fasta', '.aln') }"
        print(datetime.now(ZoneInfo("America/Sao_Paulo")), flush=True)
        print(
            f"aligning { aln_dir }/{ fasta_file.replace('.fasta', '.aln') }",
            flush=True
        )
        subprocess.run([
            "clustalo", "-i", f"{ fasta_dir }/{ fasta_file }",
            "-o", aln_file,
            "--threads", "50",
            "--log", "clustalo.log",
            "--verbose",
            "--force"
        ])


def align_sequences_multiple_differentiative(ref_seq):
    fasta_dir = f"fasta/{ ref_seq }"
    aln_dir = f"aln/{ ref_seq }"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    aln_file = f"{ aln_dir }/{ ref_seq }.aln"
    print(f"aligning { aln_dir }/{ ref_seq }")
    subprocess.run([
        "clustalo", "-i", f"{ fasta_dir }/{ ref_seq }.fasta",
        "-o", aln_file,
        "--threads", "50",
        "--log", "clustalo.log",
        "--verbose",
        "--force"
    ])


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
    min_dash = 999999999
    for seq_ref, seq in subseqs.items():
        if seq.count("-") < min_dash:
            min_dash = seq.count("-")
            reference_refsq = seq_ref
            reference_seq = subseqs[reference_refsq]
    for seq_ref, seq in subseqs.items():
        if seq_ref == reference_refsq:
            continue
        if seq != reference_seq:
            seq_mutations = find_mutation_positions(reference_seq, seq)
            mutations[seq_ref] = seq_mutations
    return reference_seq, mutations


def add_mutations_multiple_refseq_to_tree(tree):
    new_tree = clone_tree(tree, Node)
    for leaf in get_leaves_multiple_refseq(tree):
        add_path_to_tree(
            new_tree,
            f"{ leaf.path_name }/genome",
            node_attrs={
                "rank": "genome",
            },
        )
        for aln_file in [
            os.path.splitext(f)[0]
            for f in os.listdir(f"aln/{ leaf.node_name }")
            if os.path.isfile(f)
        ]:
            reference_seq, mutations = extract_mutations_from_aln(
                leaf.node_name, aln_file)
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
                    f"{ leaf.path_name }/genome/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_seq
                    },
                )
            else:
                add_path_to_tree(
                    new_tree,
                    f"{ leaf.path_name }/genome/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_seq,
                        "mutations": mutations
                    },
                )
    return new_tree


def add_mutations_multiple_differentiative_to_tree(tree):
    new_tree = clone_tree(tree, Node)
    genomes = find_attrs(tree, "rank", "genome")
    for ref_seq in get_leaves_multiple_differentiative(tree).keys():
        for gen in genomes:
            if gen.node_name == ref_seq:
                leaf = gen
                parent = gen.parent
        add_path_to_tree(
            new_tree,
            f"{ parent.path_name }/{ ref_seq }",
            node_attrs={
                "rank": "genome",
            },
        )
        for aln_file in [
            os.path.splitext(f)[0]
            for f in os.listdir(f"aln/{ ref_seq }")
            if os.path.isfile(f)
        ]:
            reference_seq, mutations = extract_mutations_from_aln(
                ref_seq, aln_file)
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
                    f"{ leaf.path_name }/genome/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_seq
                    },
                )
            else:
                add_path_to_tree(
                    new_tree,
                    f"{ leaf.path_name }/genome/{ aln_file }",
                    node_attrs={
                        "rank": "sequence",
                        "seq": reference_seq,
                        "mutations": mutations
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


def get_maximum_extract_subseqs(leaf, min_subseq_len):
    sequences = find_attrs(leaf, "rank", "sequence")
    max_subseqs = 0
    for seq in sequences:
        seq_len = len(seq.get_attr("seq"))
        if seq_len >= min_subseq_len:
            max_subseqs += seq_len - min_subseq_len
        if seq.get_attr("mutations") is None:
            seq.set_attrs({"max_subseqs": max_subseqs})
        else:
            point_mutations = {}
            indels = {}
            seq_ref_mutations = seq.get_attr("mutations")
            for seq_ref, mutations in seq_ref_mutations.items():
                seq_ref_point_mutations = {}
                seq_ref_indels = {}
                for pos, mut in mutations.items():
                    if (not pos - 1 in mutations.keys() and
                            not pos + 1 in mutations.keys()):
                        seq_ref_point_mutations[pos] = mut
                    else:
                        seq_ref_indels[pos] = mut
                point_mutations[seq_ref] = seq_ref_point_mutations
                indels[seq_ref] = seq_ref_indels
            max_subseqs += extract_max_subseqs_mutations(
                point_mutations, indels, min_subseq_len)


def get_mutations_inside_window(seq_ref_mutations, window_size):
    window_groups = {}
    for seq_ref, mutations in seq_ref_mutations.items():
        checked_pos = []
        window_groups[seq_ref] = {}
        for pos_a, mut_a in mutations.items():
            for pos_b, mut_b in mutations.items():
                if pos_a == pos_b:
                    continue
                if {pos_a, pos_b} in checked_pos:
                    continue
                else:
                    checked_pos.append({pos_a, pos_b})
                if abs(pos_a - pos_b) <= window_size:
                    if pos_a not in window_groups[seq_ref].keys():
                        window_groups[seq_ref][pos_a] = {pos_a: mut_a,
                                                         pos_b: mut_b}
                    else:
                        window_groups[seq_ref][pos_a][pos_b] = mut_b
    return window_groups


def get_window_group_score(window_group, window_size):
    all_gap = True
    ordered_mutation_keys = sorted(window_group.keys())
    for mutation in window_group.values():
        if mutation != '-':
            all_gap = False
            break
    if all_gap:
        gap_groups = {}
        for index_pos in range(len(ordered_mutation_keys)):
            if index_pos <= len(ordered_mutation_keys) - 1:
                if ordered_mutation_keys[index_pos + 1] - \
                   ordered_mutation_keys[index_pos] == 1:
                    if ordered_mutation_keys[index_pos] \
                            not in gap_groups.keys():
                        gap_groups[ordered_mutation_keys[index_pos]] = [
                            ordered_mutation_keys[index_pos]]
                    gap_groups[ordered_mutation_keys[index_pos]].append(
                        ordered_mutation_keys[index_pos + 1])
                else:
                    gap_groups[ordered_mutation_keys[index_pos]] = [
                        ordered_mutation_keys[index_pos]]
        gaps = {}
        num_gaps = 0
        for gap_start, gap_list in gap_groups.items():
            gap_end = max(gap_list)
            gaps[gap_start] = gap_end
            num_gaps += 1
        discount = 0
        prev_end = None
        for gap_start, gap_end in gaps.items():
            if prev_end:
                discount += gap_start - prev_end
            prev_end = gap_end

        return ((window_size - 1) * num_gaps) - discount

    else:
        biggest_index = -1
        mutation = '-'
        while (mutation == '-' and
                biggest_index > -len(window_group)):
            biggest_mutation_key = ordered_mutation_keys[biggest_index]
            mutation = window_group[biggest_mutation_key]
            biggest_index -= 1
        smallest_index = 0
        mutation = '-'
        while (mutation == '-' and
                smallest_index < len(window_group)):
            smallest_mutation_key = ordered_mutation_keys[smallest_index]
            mutation = window_group[smallest_mutation_key]
            smallest_index += 1
        smallest_mutation_key = ordered_mutation_keys[0]

        biggest_mutation_distance = biggest_mutation_key - \
            smallest_mutation_key
    

def get_point_mutation_score(pos, point_mutations, indels, window_size):
    max_subseqs = 0
    checked_seq_refs = []
    checked_point_mutations = {}
    for seq_ref_a, seq_ref_point_mutations_a \
            in point_mutations.items():
        for seq_ref_b, seq_ref_point_mutations_b \
                in point_mutations.items():
            if seq_ref_a == seq_ref_b:
                continue
            if {seq_ref_a, seq_ref_b} in checked_seq_refs:
                continue
            else:
                checked_seq_refs.append({seq_ref_a, seq_ref_b})
                if seq_ref_a not in checked_point_mutations.keys():
                    checked_point_mutations[seq_ref_a] = []
                if seq_ref_b not in checked_point_mutations.keys():
                    checked_point_mutations[seq_ref_b] = []
            comon_positions = set(
                seq_ref_point_mutations_a.keys()
            ) & set(seq_ref_point_mutations_b.keys())
            for pos in comon_positions:
                if (seq_ref_point_mutations_a[pos] == '-' and
                        seq_ref_point_mutations_b[pos] == '-'):
                    continue
                if (pos in checked_point_mutations[seq_ref_a] and
                        pos in checked_point_mutations[seq_ref_b]):
                    continue
                elif (pos not in checked_point_mutations[seq_ref_a] and
                        pos not in checked_point_mutations[seq_ref_b]):
                    if (seq_ref_point_mutations_a[pos] !=
                            seq_ref_point_mutations_b[pos]):
                        if (seq_ref_point_mutations_a[pos] == '-' or
                                seq_ref_point_mutations_b[pos] == '-'):
                            max_subseqs += 1
                        else:
                            max_subseqs += 2
                    else:
                        max_subseqs += 1
                    checked_point_mutations[seq_ref_a].append(pos)
                    checked_point_mutations[seq_ref_b].append(pos)
                else:
                    if seq_ref_point_mutations_a[pos] != \
                            seq_ref_point_mutations_b[pos]:
                        if (pos not in
                                checked_point_mutations[seq_ref_a] and
                                seq_ref_point_mutations_a[pos] != '-'):
                            max_subseqs += 1
                        if (pos not in
                                checked_point_mutations[seq_ref_b] and
                                seq_ref_point_mutations_b[pos] != '-'):
                            max_subseqs += 1
                    if pos not in \
                            checked_point_mutations[seq_ref_a]:
                        checked_point_mutations[seq_ref_a] \
                            .append(pos)
                    if pos not in \
                            checked_point_mutations[seq_ref_b]:
                        checked_point_mutations[seq_ref_b] \
                            .append(pos)
            for pos, mut in seq_ref_point_mutations_a.items():
                if (pos not in comon_positions and
                        pos not in checked_point_mutations[seq_ref_a]):
                    if mut != '-':
                        max_subseqs += 1
                    checked_point_mutations[seq_ref_a].append(pos)
            for pos, mut in seq_ref_point_mutations_b.items():
                if (pos not in comon_positions and
                        pos not in checked_point_mutations[seq_ref_b]):
                    if mut != '-':
                        max_subseqs += 1
                    checked_point_mutations[seq_ref_b].append(pos)

def extract_max_subseqs_mutations(
        point_mutations, indels, window_size):
    max_subseqs = 0
    checked_indels = set()
    checked_seq_refs = []
    for seq_ref_a, seq_ref_indels_a in indels.items():
        for seq_ref_b, seq_ref_indels_b in indels.items():
            if seq_ref_a == seq_ref_b:
                continue
            if {seq_ref_a, seq_ref_b} in checked_seq_refs:
                continue
            else:
                checked_seq_refs.append({seq_ref_a, seq_ref_b})
            ordered_a = collections.OrderedDict(
                sorted(seq_ref_indels_a.items()))
            indels_a = []
            ordered_b = collections.OrderedDict(
                sorted(seq_ref_indels_b.items()))
            indels_b = []
            prev_pos_a = None
            current_indel = {}
            for pos_a, mut_a in ordered_a.items():
                if pos_a - prev_pos_a == 1:
                    current_indel[pos_a] = mut_a
                else:
                    indels_a.append(current_indel)
                    current_indel = {pos_a: mut_a}
                prev_pos_a = pos_a
            for pos_b, mut_b in ordered_b.items():
                if pos_b - prev_pos_a == 1:
                    current_indel[pos_b] = mut_b
                else:
                    indels_b.append(current_indel)
                    current_indel = {pos_b: mut_b}
                prev_pos_a = pos_b
            checked_indel_groups = []
            for indel_a in indels_a:
                for indel_b in indels_b:
                    if {indel_a, indel_b} in checked_indel_groups:
                        continue
                    else:
                        checked_seq_refs.append({indel_a, indel_b})
                    if indel_a == indel_b:
                        if indel_a not in checked_indels:
                            checked_indels.add(indel_a)
                            max_subseqs += len(indel_a)
                    else:
                        for pos_a, mut_a in indel_a.items():
                            if pos_a in indel_b.keys():
                                if mut_a != indel_b[pos_a]:
                                    max_subseqs += 1
    
    

    return max_subseqs
