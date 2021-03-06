import sys
import os
from collections import defaultdict, namedtuple

ncbi_folder = '/home/hawkjo/dbs/NCBI/DATA/'


def build_geneid_given_refseq_rna_acc_dict():
    gene2refseq_file = os.path.join(ncbi_folder, 'gene2refseq')
    geneid_given_refseq = {}
    with open(gene2refseq_file) as f:
        next(f)
        for line in f:
            var = line.strip().split('\t')
            geneid_given_refseq[var[3]] = var[1]
            geneid_given_refseq[var[3].split('.')[0]] = var[1]
    return geneid_given_refseq


def build_refseq_rna_acc_given_geneid_dict():
    gene2refseq_file = os.path.join(ncbi_folder, 'gene2refseq')
    refseq_given_geneid = {}
    with open(gene2refseq_file) as f:
        next(f)
        for line in f:
            var = line.strip().split('\t')
            refseq_given_geneid[var[1]] = [var[3], var[3].split('.')[0]]
    return refseq_given_geneid


def build_GOid_given_gene_dict():
    gene2go_file = os.path.join(ncbi_folder, 'gene2go')
    GOid_given_gene = defaultdict(list)
    with open(gene2go_file) as f:
        next(f)
        for line in f:
            var = line.strip().split('\t')
            GOid_given_gene[var[1]].append(var[2])
    return GOid_given_gene


def build_go_entries_given_gene_id_dict():
    gene2go_entry = namedtuple('gene2go_entry',
                               'tax_id GeneID GO_ID Evidence Qualifier GO_term PubMed Category')
    gene2go_file = os.path.join(ncbi_folder, 'gene2go')
    go_entries_given_gene_id = defaultdict(list)
    with open(gene2go_file) as f:
        next(f)
        for line in f:
            # Make and clean entry
            # Qualifier is an optional argument. Should be blank rather than dash.
            # PubMed will be used as reference argument. If not given, replace with all zeros.
            # Category will be used as aspect, which must be 'F', 'P', or 'C'.
            var = line.strip().split('\t')
            if var[4] == '-':
                var[4] = ''           # Qualifier
            if var[6] == '-':
                var[6] = '0000000'    # PubMed
            if var[7][0] in ['F', 'P', 'C']:
                var[7] = var[7][0]
            else:
                sys.exit('Bad Category:\n' + line)
            entry = gene2go_entry(*var)
            go_entries_given_gene_id[var[1]].append(entry)
    return go_entries_given_gene_id


def build_symbol_given_gene_id_dict():
    gene_info_file = os.path.join(ncbi_folder, 'gene_info')
    symbol_given_gene_id = {}
    with open(gene_info_file) as f:
        next(f)
        for line in f:
            var = line.strip().split('\t')
            symbol_given_gene_id[var[1]] = var[2]
    return symbol_given_gene_id


def output_GO_GAF2_file_from_ncbi_geneids(lst_of_gene_ids, filename):
    print 'Reading GO entries...'
    go_entries_given_gene_id = build_go_entries_given_gene_id_dict()
    print 'Reading NCBI symbols...'
    symbol_given_gene_id = build_symbol_given_gene_id_dict()

    print 'Writing output...'
    missing_symbols = 0
    missing_go_entries = 0
    with open(filename, 'w') as out:
        out.write('!gaf-version: 2.0\n')
        for gene_id in lst_of_gene_ids:
            if gene_id not in symbol_given_gene_id:
                missing_symbols += 1
                continue
            elif gene_id not in go_entries_given_gene_id:
                missing_go_entries += 1
                continue

            gene_symbol = symbol_given_gene_id[gene_id]
            for go_entry in go_entries_given_gene_id[gene_id]:
                out.write('\t'.join(['NCBI_Gene', gene_id, gene_symbol, go_entry.Qualifier,
                                     go_entry.GO_ID, 'PMID:'+go_entry.PubMed, go_entry.Evidence,
                                     '', '', go_entry.Category, '', '', 'gene_product',
                                     'taxon:' + go_entry.tax_id, '20140212', 'NCBI_Gene', '', ''])
                          + '\n')
    print
    print 'Summary:'
    print 'Input gene IDs:', len(lst_of_gene_ids)
    print 'Missing symbols:', missing_symbols
    print 'IDs with zero GO entries:', missing_go_entries


def build_gene_id_given_protein_gi_dict():
    gene_id_given_protein_gi = {}
    filename = os.path.join(ncbi_folder, 'gene2accession')
    with open(filename) as f:
        next(f)
        for line in f:
            var = line.strip().split()
            gene_id_given_protein_gi[var[6]] = var[1]
    return gene_id_given_protein_gi


def build_protein_gi_given_gene_id_dict():
    protein_gi_given_gene_id = {}
    filename = os.path.join(ncbi_folder, 'gene2accession')
    with open(filename) as f:
        next(f)
        for line in f:
            var = line.strip().split()
            protein_gi_given_gene_id[var[1]] = var[6]
    return protein_gi_given_gene_id


def build_ncbi_symbol_given_id_dict():
    symbol_given_id = {}
    fname = os.path.join(ncbi_folder, 'gene_info')
    with open(fname) as f:
        next(f)
        for line in f:
            var = line.strip().split()
            symbol_given_id[var[1]] = var[2]
    return symbol_given_id
