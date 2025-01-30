import taxoniq
import matplotlib.pyplot as plt

def get_seqref_tax_ids():
    tax_ids = {}
    with open('/home/andrey/kraken2/viral/seqid2taxid.map') as f:
        for line in f.readlines():
            seq_ref = line.rstrip().split('|')[-1].split()[0]
            tax_id = int(line.rstrip().split('|')[-1].split()[-1])
            if tax_id in tax_ids.keys():
                tax_ids[tax_id].append(seq_ref)
            else:
                tax_ids[tax_id] = [seq_ref]
    f.close()

    return tax_ids

def get_spec_tax_ids():
    tax_ids = {}
    with open('/home/andrey/kraken2/viral/taxid2specid.map') as f:
        for line in f.readlines():
            tax_id = line.rstrip().split()[0]
            idx = line.rstrip().split()[-1]
            tax_ids[int(tax_id)] = int(idx)
    f.close()

    return tax_ids

def get_spec_seq_lens():
    seq_lens = {}
    with open('/home/andrey/kraken2/viral/specid2seqlen.map') as f:
        for line in f.readlines():
            idx = line.rstrip().split()[0]
            seq_len = line.rstrip().split()[-1]
            seq_lens[int(idx)] = int(seq_len)
    f.close()

    return seq_lens

def get_tax_id_taxoniq(taxon):
    taxon_str = str(taxon)
    str_tax_id = taxon_str[taxon_str.find("(")+1:taxon_str.find(")")]
    tax_id = int(str_tax_id)

    return tax_id

def populate_dict(target_dict, key, val):
    if key in target_dict.keys():
        target_dict[key].append(val)
    else:
        target_dict[key] = [val]

spec_tax_ids = get_seqref_tax_ids()
spec_seq_lens = get_spec_seq_lens()

genus_tax_id2species_id = {}
family_tax_id2species_id = {}
order_tax_id2species_id = {}
class_tax_id2species_id = {}
phylum_tax_id2species_id = {}
kingdom_tax_id2species_id = {}

taxon_dicts = {'genus': genus_tax_id2species_id, 'family': family_tax_id2species_id, 'order': order_tax_id2species_id,
               'class': class_tax_id2species_id, 'phylum': phylum_tax_id2species_id, 'kingdom': kingdom_tax_id2species_id}

tax_ranks_mapfiles = {'genus': 'taxid2genid', 'family': 'taxid2famid', 'order': 'taxid2ordid', 
                  'class': 'taxid2clasid', 'phylum': 'taxid2phylid', 'kingdom': 'taxid2kingid'}

id_ranks_mapfiles = {'genus': 'seqref2genid', 'family': 'seqref2famid', 'order': 'seqref2ordid',
                     'class': 'seqref2clasid', 'phylum': 'seqref2phylid', 'kingdom': 'seqref2kingid'}

taxids2rankids = {}

for tax_id, seqrefs in spec_tax_ids.items():
    t = taxoniq.Taxon(tax_id)
    ranked_taxons = t.ranked_lineage
    for ranked_taxon in ranked_taxons:
        ranked_tax_id = get_tax_id_taxoniq(ranked_taxon)
        if ranked_taxon.rank.name in taxon_dicts.keys():
            for seqref in seqrefs:
                populate_dict(taxon_dicts[ranked_taxon.rank.name], ranked_tax_id, seqref)
            taxids2rankids[ranked_taxon.rank.name] = {}

for rank, ranktaxid2seqrefs in taxon_dicts.items():
    with open('/home/andrey/kraken2/viral/'+tax_ranks_mapfiles[rank]+'.map', 'w') as f:
        idx = 1
        for ranked_tax_id in ranktaxid2seqrefs.keys():
            f.write(str(ranked_tax_id)+' '+str(idx)+'\n')
            idx += 1
            taxids2rankids[rank][ranked_tax_id] = idx
    f.close()

    with open('/home/andrey/kraken2/viral/'+id_ranks_mapfiles[rank]+'.map', 'w') as f:
        for ranked_tax_id, seqrefs in ranktaxid2seqrefs.items():
            for seqref in seqrefs:
                f.write(str(seqref)+' '+str(ranked_tax_id)+'\n')
    f.close()

