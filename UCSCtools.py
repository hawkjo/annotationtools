import sys


def reduce_to_longest_isoforms(lst_of_gene_records):
    """reduce_to_longest_isoforms(lst_of_gene_records):

    Given a list of gene records with a UCSC name as gene.name, returns a list of the longest
    isoforms.
    """
    knownIsoformfile = '/home/hawkjo/genomes/hg19/annotation/knownIsoforms.txt'
    ucscname2isoform_id = {}
    for line in open(knownIsoformfile):
        var = line.strip().split()
        ucscname2isoform_id[var[1]] = int(var[0])

    longest_isoforms_by_id = {}
    for gene in lst_of_gene_records:
        iso_id = ucscname2isoform_id[gene.name]
        if iso_id not in longest_isoforms_by_id or len(gene) > longest_isoforms_by_id[iso_id]:
            longest_isoforms_by_id[iso_id] = gene

    longest_isoform_records = longest_isoforms_by_id.values()
    longest_isoform_records.sort(key=lambda gene: lst_of_gene_records.index(gene))
    return longest_isoform_records


def build_refseq_given_ucscname_dict():
    knownToRefSeq_file = '/home/hawkjo/genomes/hg19/annotation/knownToRefSeq.txt'
    refseq_given_ucscname = {}
    for line in open(knownToRefSeq_file):
        var = line.strip().split()
        refseq_given_ucscname[var[0]] = var[1]
    return refseq_given_ucscname


def build_ucscname_given_refseq_dict():
    knownToRefSeq_file = '/home/hawkjo/genomes/hg19/annotation/knownToRefSeq.txt'
    ucscname_given_refseq = {}
    for line in open(knownToRefSeq_file):
        var = line.strip().split()
        ucscname_given_refseq[var[1]] = var[0]
    return ucscname_given_refseq


def build_ucscname_given_ncbi_gene_id_dict():
    from NCBItools import build_refseq_rna_acc_given_geneid_dict
    refseq_list_given_ncbi_id = build_refseq_rna_acc_given_geneid_dict()
    ucscname_given_refseq = build_ucscname_given_refseq_dict()
    ucscname_given_ncbi_id = {}
    for ncbi_id, refseq_list in refseq_list_given_ncbi_id.items():
        for refseq in refseq_list:
            try:
                # If either ncbi_id not assigned yet or no ucscname, not a problem.
                if ucscname_given_ncbi_id[ncbi_id] != ucscname_given_refseq[refseq]:
                    sys.exit('Multiple ucscnames for single ncbi id')
            except KeyError:
                pass
            try:
                # If no ucscname, not a problem.
                ucscname_given_ncbi_id[ncbi_id] = ucscname_given_refseq[refseq]
            except KeyError:
                pass
    return ucscname_given_ncbi_id


def build_ncbi_gene_id_given_ucscname_dict():
    ncbi_gene_id_given_ucscname = {}
    for gene_id, ucscname in build_ucscname_given_ncbi_gene_id_dict().items():
        assert ucscname not in ncbi_gene_id_given_ucscname
        ncbi_gene_id_given_ucscname[ucscname] = gene_id
    return ncbi_gene_id_given_ucscname


def build_ncbi_symbol_given_ucscname_dict():
    from NCBItools import build_ncbi_symbol_given_id_dict
    ncbi_symbol_given_id = build_ncbi_symbol_given_id_dict()
    ncbi_symbol_given_ucscname = {}
    for ucscname, gene_id in build_ncbi_gene_id_given_ucscname_dict().items():
        try:
            ncbi_symbol_given_ucscname[ucscname] = ncbi_symbol_given_id[gene_id]
        except KeyError:
            pass
    return ncbi_symbol_given_ucscname
