import sys
import os
from glob import glob
import itertools
from collections import defaultdict, Counter
from Bio import SeqIO


def get_qualifier_force_single(feature, qual_name, allow_zero=False):
    try:
        qual = feature.qualifiers[qual_name]
    except KeyError:
        if allow_zero:
            return ''
        else:
            raise
    assert len(qual) == 1, (feature, qual_name, qual)
    return qual[0]


class GenBankCDS(object):
    def __init__(self, rec, cds_feature):
        self.gene_id = get_qualifier_force_single(cds_feature, 'gene')
        self.chrm = rec.id
        self.strand = cds_feature.location.strand
        assert self.strand in [-1, 1]
        self.get_parts_and_boundaries(cds_feature)
        self.length = sum(abs(part[0] - part[1]) + 1 for part in self.parts)
        self.codon_start = get_qualifier_force_single(cds_feature, 'codon_start')
        try:
            self.protein_id = get_qualifier_force_single(cds_feature, 'protein_id')
            self.translation = get_qualifier_force_single(cds_feature, 'translation')
        except KeyError:
            assert 'exception' in cds_feature.qualifiers, cds
            self.protein_id = None
            self.translation = None
        
    def get_parts_and_boundaries(self, cds_feature):
        self.parts = []
        self.boundaries = set()
        for loc in cds_feature.location.parts:
            self.parts.append((loc.nofuzzy_start, loc.nofuzzy_end, loc.strand))
            self.boundaries.update((loc.nofuzzy_start, loc.nofuzzy_end))

    def __str__(self):
        return '%s: length %d, %s' % (self.gene_id, self.length, self.parts)

class GenBankGene(object):
    def __init__(self, rec, gene_feature):
        self.gene_id = get_qualifier_force_single(gene_feature, 'gene')
        syn_str = get_qualifier_force_single(gene_feature, 'gene_synonym', allow_zero=True)
        syn_str.replace(' ', '')
        self.gene_synonyms = set(syn_str.split(';'))
        note = get_qualifier_force_single(gene_feature, 'note', allow_zero=True)
        self.is_readthrough = bool('readthrough' in note.lower())

        self.chrm = rec.id
        self.gene_start = gene_feature.location.nofuzzy_start
        self.gene_end = gene_feature.location.nofuzzy_end

        self.cdss = []
        self.cds_parts = set()
        self.cds_boundaries = set()
        self.longest_cds = None
        self.longest_cds_length = None

    def contains_cds(self, cds):
        for part in cds.parts:
            if not (self.gene_start <= part[0] <= self.gene_end 
                    and self.gene_start <= part[1] <= self.gene_end):
                return False
        return True

    def add_cds(self, cds):
        # First assert that gene id is the same and all parts are contained in the gene.
        assert cds.gene_id == self.gene_id or cds.gene_id in self.gene_synonyms, '%s\n%s' % (self, cds)
        assert self.contains_cds(cds), '%s\n%s' % (self, cds)

        self.cdss.append(cds)
        self.cds_parts.update(cds.parts)
        self.cds_boundaries.update(cds.boundaries)
        if cds.length > self.longest_cds_length:
            self.longest_cds = cds
            self.longest_cds_length = cds.length

    def overlaps_region(self, chrm, start, end):
        if self.chrm == chrm \
                and (self.gene_start <= start <= self.gene_end
                     or self.gene_start <= end <= self.gene_end
                     or start <= self.gene_start <= end
                     or start <= self.gene_end <= end):
            return True
        else:
            return False

    def overlaps_gene(self, other_gene):
        return self.overlaps_region(other_gene.chrm, other_gene.gene_start, other_gene.gene_end)

    def __str__(self):
        if self.longest_cds:
            longest_prot_id = self.longest_cds.protein_id
            all_prot_ids = ','.join([cds.protein_id for cds in self.cdss])
        else:
            longest_prot_id = None
            all_prot_ids = None
        return '%s\t%s:%d-%d\tnum_CDSs=%s\tlongest_CDS_protein=%s\tlongest_CDS_len=%s\tall_CDS_proteins=%s' \
                % (self.gene_id,
                   self.chrm,
                   self.gene_start,
                   self.gene_end,
                   len(self.cdss),
                   longest_prot_id,
                   self.longest_cds_length,
                   all_prot_ids)


class GenBankGenomeGenes(object):
    def __init__(self, fpath):
        assert fpath.endswith('genomic.gbff'), 'Incorrect file type.'
        self.load_genes(fpath)

    def get_genes_overlapping_region(self, chrm, start, end):
        overlapping_genes = []
        for gene in self.gene_list_given_chrm[chrm]:
            if gene.overlaps_region(chrm, start, end):
                overlapping_genes.append(gene)
        return overlapping_genes

    def load_genes(self, fpath):
        """
        This method reads through refseq *_genomic.gbff files and extracts CDS isoform information.
        """
        # Read in all the genes
        print 'Reading genes from genomic file'
        genes_given_id = defaultdict(list)
        readthrough_genes = set()
        for rec in SeqIO.parse(open(fpath), 'gb'):
            if ('FIX_PATCH' in rec.annotations['keywords'] 
                    or 'NOVEL_PATCH' in rec.annotations['keywords'] 
                    or 'ALTERNATE_LOCUS' in rec.annotations['keywords']):
                continue
            sys.stdout.write('.')
            sys.stdout.flush()
            for feature in rec.features:
                if feature.type == 'gene':
                    gene = GenBankGene(rec, feature)
                    if gene.is_readthrough:
                        # We reject any readthrough genes.
                        readthrough_genes.add(gene.gene_id)
                    else:
                        genes_given_id[gene.gene_id].append(gene)
                elif feature.type == 'CDS':
                    cds = GenBankCDS(rec, feature)
                    if cds.protein_id is None:
                        # Some CDSs have exceptions indicating no protein is directly coded. Skip those
                        continue
                    for gene in genes_given_id[cds.gene_id]:
                        if gene.contains_cds(cds):
                            gene.add_cds(cds)
                            break
                    else:
                        if cds.gene_id not in readthrough_genes:
                            sys.exit('Error: Did not find gene for cds:\n\t%s' % cds)
        print
    
        # In my observations, all large gene groups which use the same exons should be called isoforms
        # of the same gene.  They appear to simply be the invention of over-ambitious scientists,
        # potentially for more genes by their name or something like that.
        print 'Combining isoform "gene families"'
        # Find groups which share cds 'exons'
        genes_given_part = defaultdict(list)
        for gene in itertools.chain(*genes_given_id.values()):
            for part in gene.cds_parts:
                genes_given_part[part].append(gene)
        overlapping_proteins = set()
        for genes in genes_given_part.values():
            if len(genes) > 1:
                overlapping_proteins.add(', '.join(sorted([gene.gene_id for gene in genes])))
        # Keep only the longest member of the family.
        for family_str in overlapping_proteins:
            family = family_str.split(', ')
            if len(family) <= 2:
                continue
            family_gene_iter = itertools.chain(*[genes_given_id[gene_id] for gene_id in family])
            max_len_gene = max(family_gene_iter, key=lambda gene: gene.longest_cds_length)
            for gene_id in family:
                if gene_id != max_len_gene.gene_id:
                    del genes_given_id[gene_id]
    
        # Remove genes with no cdss
        print 'Removing genes with no CDSs'
        for gene_id in genes_given_id.keys():
            genes = genes_given_id[gene_id]  # Handle with care due to deletions
            genes_to_keep = []
            for i in xrange(len(genes)):
                if genes[i].longest_cds is not None:
                    genes_to_keep.append(i)
            if not genes_to_keep:
                del genes_given_id[gene_id]
            else:
                genes_given_id[gene_id] = [genes[i] for i in genes_to_keep]
    
        # Verify good results
        print 'Checking results'
        print 'Checking gene names and overlaps'
        dist_num_genes_for_gene_id = Counter()
        genes_with_dups = set()
        for gene_id, genes in genes_given_id.items():
            dist_num_genes_for_gene_id[len(genes)] += 1
            if len(genes) > 1:
                genes_with_dups.add(gene_id)
            for i, g1 in enumerate(genes):
                for g2 in genes[i+1:]:
                    assert not g1.overlaps_gene(g2), '%s\n%s\n%s' % (gene_id, g1, g2)
        print 'Dist of number of genes per gene_id:'
        for k, v in dist_num_genes_for_gene_id.items():
            print '\t%s: %s' % (k, v)
        print 'Genes with duplicates:', ', '.join(sorted(genes_with_dups))
    
        # Organize loaded genes by chrm
        gene_list_given_chrm = defaultdict(list)
        for gene in itertools.chain(*genes_given_id.values()):
            gene_list_given_chrm[gene.chrm].append(gene)
        for chrm in gene_list_given_chrm.keys():
            gene_list_given_chrm[chrm].sort(key=lambda gene: gene.gene_start)

        self.genes_given_id = dict(genes_given_id)
        self.gene_list_given_chrm = dict(gene_list_given_chrm)
