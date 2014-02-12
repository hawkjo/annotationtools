import sys

def reduce_to_longest_isoforms( lst_of_gene_records ):
    """Given a list of gene records with a UCSC name as gene.name, returns a list of the longest isoforms."""
    knownIsoformfile = '/home/hawkjo/scratch/genomes/hg19/annotation/knownIsoforms.txt'
    ucscname2isoform_id = {}
    for line in open(knownIsoformfile):
        var = line.strip().split()
        ucscname2isoform_id[var[1]] = int(var[0])

    longest_isoforms_by_id = {}
    for gene in lst_of_gene_records:
        iso_id = ucscname2isoform_id[ gene.name ]
        if iso_id not in longest_isoforms_by_id or len(gene) > longest_isoforms_by_id[iso_id]:
            longest_isoforms_by_id[iso_id] = gene
    
    longest_isoform_records = longest_isoforms.values()
    longest_isoform_records.sort( key=lambda gene: lst_of_gene_records.index( gene ) )
    return longest_isoform_records

def ucscname2refseq( lst_of_ucsc_names ):
    knownToRefSeq_file = '/home/hawkjo/scratch/genomes/hg19/annotation/knownToRefSeq.txt'
    ucscname2refseq = {}
    for line in open(knownToRefSeq_file):
        var = line.strip().split()
        ucscname2refseq[var[0]] = var[1]
    return [ ucscname2refseq[ ucscname ] for ucscname in lst_of_ucsc_names ]
