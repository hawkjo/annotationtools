import sys
import time
from collections import defaultdict


class GenomicRegion(object):
    """GenomicRegion

    A contiguous region of the genome, along with capability to store the names of any genes it may
    overlap
    """
    def __init__(self, gene_names, chrom, start, end):
        if isinstance(gene_names, str):
            self.geneNames = set([gene_names])
        elif isinstance(gene_names, list):
            self.geneNames = set(gene_names)
        elif isinstance(gene_names, set):
            self.geneNames = gene_names
        else:
            sys.exit("Invalid input type for geneNames")
        self.chrom = chrom
        self.start = start
        self.end = end

    def merge(self, region):
        if self.chrom != region.chrom: return
        self.geneNames.union(region.geneNames)
        self.start = min(self.start, region.start)
        self.end = max(self.end, region.end)

    def can_merge(self, region):
        if self.chrom != region.chrom:
            return False
        elif self.start <= region.start <= self.end or \
                self.start <= region.end <= self.end or \
                region.start <= self.start <= region.end:
            # Definition of can_merge here includes the possibility of two regions
            # lying immediately next to each other but not having any bases in common.
            # The math above is due to the half-open convention of start and end.
            return True
        else:
            return False


class Gene(object):
    def __init__(self, line=None, filetype='knownGene'):
        assert filetype in ['knownGene', 'ensGene'], \
            'Unknown filetype: %s' % filetype
        self.filetype = filetype
        if line:
            self._parse(line)
        else:
            sys.exit("No line specified for gene.")

    def _parse(self, line):
        var = line.strip().split('\t')

        if self.filetype == 'knownGene':
            self.name, self.chrom, self.strand = var[0:3]  # strand is '+' or '-'
            self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCount = map(int, var[3:8])
            self.exonStarts = [int(s) for s in var[8].split(',')[:self.exonCount]]
            self.exonEnds = [int(s) for s in var[9].split(',')[:self.exonCount]]
        elif self.filetype == 'ensGene':
            self.bin = int(var[0])
            self.name, self.chrom, self.strand = var[1:4]  # strand is '+' or '-'
            self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCount = map(int, var[4:9])
            self.exonStarts = [int(s) for s in var[9].split(',')[:self.exonCount]]
            self.exonEnds = [int(s) for s in var[10].split(',')[:self.exonCount]]
            self.score = int(var[11])
            self.name2, self.cdsStartStat, self.cdsEndStat = var[12:15]
            self.exonFrams = [int(s) for s in var[15].split(',')[:self.exonCount]]

        self.introns = []
        for intron_start, intron_end in zip(self.exonEnds, self.exonStarts[1:]):
            self.introns.append(GenomicRegion(self.name, self.chrom, intron_start, intron_end))
        self.exons = []
        self.CDSs = []
        self.UTRs = []
        for start, end in zip(self.exonStarts, self.exonEnds):
            self.exons.append(GenomicRegion(self.name, self.chrom, start, end))

            # Now to break the exons up into coding and non-coding regions
            if self.cdsStart == self.cdsEnd:
                # Non-coding gene. All exons are UTRs
                self.UTRs.append(GenomicRegion(self.name, self.chrom, start, end))
            elif start >= self.cdsStart and end <= self.cdsEnd:
                # Exon entirely within coding region
                self.CDSs.append(GenomicRegion(self.name, self.chrom, start, end))
            elif start < self.cdsStart and self.cdsStart < end <= self.cdsEnd:
                # Exon has coding region and 5' UTR
                self.UTRs.append(GenomicRegion(self.name, self.chrom, start, self.cdsStart))
                self.CDSs.append(GenomicRegion(self.name, self.chrom, self.cdsStart, end))
            elif self.cdsStart <= start < self.cdsEnd and end > self.cdsEnd:
                # Exon has coding reigon and 3' UTR
                self.CDSs.append(GenomicRegion(self.name, self.chrom, start, self.cdsEnd))
                self.UTRs.append(GenomicRegion(self.name, self.chrom, self.cdsEnd, end))
            elif start < self.cdsStart < self.cdsEnd < end:
                # Exon has coding region and both 5' and 3' UTRs
                self.UTRs.append(GenomicRegion(self.name, self.chrom, start, self.cdsStart))
                self.CDSs.append(GenomicRegion(self.name, self.chrom, self.cdsStart, self.cdsEnd))
                self.UTRs.append(GenomicRegion(self.name, self.chrom, self.cdsEnd, end))
            else:
                # Coding region exists but does not overlap exon
                self.UTRs.append(GenomicRegion(self.name, self.chrom, start, end))

    def __str__(self):
        return "%s %s %s %d %d %d %d %d %s %s" % \
            (self.name, self.chrom, self.strand, self.txStart, self.txEnd,
             self.cdsStart, self.cdsEnd, self.exonCount,
             ','.join(str(start) for start in self.exonStarts),
             ','.join(str(end) for end in self.exonEnds))


def knownGene_as_dict(knownGene_file, filetype='knownGene'):
    genes = {}
    for line in open(knownGene_file):
        gene = Gene(line, filetype)
        genes[gene.name] = gene
    return genes


def knownGene_as_list(knownGene_file, filetype='knownGene'):
    return [Gene(line, filetype) for line in open(knownGene_file)]


def get_all_sorted_regions(knownGene_file, region_type):
    regions = defaultdict(list)
    for gene in knownGene_as_list(knownGene_file):
        regions[gene.chrom].extend(getattr(gene, region_type))
    for chrom in regions.keys():
        regions[chrom].sort(key=lambda region: region.start)
    return regions


def get_sorted_disjoint_regions(knownGene_file, region_type):
    sorted_disjoint_regions = defaultdict(list)
    for chrom, regions in get_all_sorted_regions(knownGene_file, region_type).items():
        for region in regions:
            if (sorted_disjoint_regions[chrom]
                    and sorted_disjoint_regions[chrom][-1].can_merge(region)):
                sorted_disjoint_regions[chrom][-1].merge(region)
            else:
                sorted_disjoint_regions[chrom].append(region)
    return sorted_disjoint_regions


class RegionNode(object):
    def __init__(self, sorted_regions_list):
        mid = int(len(sorted_regions_list)/2)
        exon_region = sorted_regions_list[mid]
        self.geneNames = exon_region.geneNames
        self.chrom = exon_region.chrom
        self.start = exon_region.start
        self.end = exon_region.end
        if len(sorted_regions_list[:mid]) > 0:
            self.left = RegionNode(sorted_regions_list[:mid])
        else:
            self.left = None
        if len(sorted_regions_list[mid+1:]) > 0:
            self.right = RegionNode(sorted_regions_list[mid+1:])
        else:
            self.right = None

    def local_overlap(self, start, end):
        if self.start <= start < self.end or \
                self.start < end <= self.end or \
                start <= self.start < end:
            # Convention here says that start and end are half-open as well
            return True
        else:
            return False

    def overlaps(self, start, end):
        if self.local_overlap(start, end):
            return True
        elif self.start >= end:
            if self.left:
                return self.left.overlaps(start, end)
            else:
                return False
        elif self.end <= start:
            if self.right:
                return self.right.overlaps(start, end)
            else:
                return False
        else:
            sys.exit("overlaps error. Impossible overlap comparison.")

    def overlapping_regions(self, start, end):
        if self.local_overlap(start, end):
            overlapping_regions = []
            if self.left:
                overlapping_regions.extend(self.left.overlapping_regions(start, end))
            overlapping_regions.append((self.start, self.end))
            if self.right:
                overlapping_regions.extend(self.right.overlapping_regions(start, end))
            return overlapping_regions
        elif self.start >= end:
            if self.left:
                return self.left.overlapping_regions(start, end)
            else:
                return []  # Equivalent to False, but compatible with list extension
        elif self.end <= start:
            if self.right:
                return self.right.overlapping_regions(start, end)
            else:
                return []  # Equivalent to False, but compatible with list extension
        else:
            sys.exit("overlapping regions error. Impossible overlap comparison.")

    def overlapping_genes(self, start, end):
        if self.local_overlap(start, end):
            overlapping_genes = []
            if self.left:
                overlapping_genes.extend(self.left.overlapping_genes(start, end))
            overlapping_genes.extend(self.geneNames)
            if self.right:
                overlapping_genes.extend(self.right.overlapping_genes(start, end))
            return overlapping_genes
        elif self.start >= end:
            if self.left:
                return self.left.overlapping_genes(start, end)
            else:
                return []  # Equivalent to False, but compatible with list extension
        elif self.end <= start:
            if self.right:
                return self.right.overlapping_genes(start, end)
            else:
                return []  # Equivalent to False, but compatible with list extension
        else:
            sys.exit("overlapping genes error. Impossible overlap comparison.")

    def dist_to_closest(self, start, end):
        if self.local_overlap(start, end):
            return 0
        elif self.start >= end:
            dist = self.start - end + 1
            if self.left:
                return min(dist, self.left.dist_to_closest(start, end))
            else:
                return dist
        elif self.end <= start:
            dist = start - self.end + 1
            if self.right:
                return min(dist, self.right.overlaps(start, end))
            else:
                return dist
        else:
            sys.exit("dist_to_closest error. Impossible overlap comparison.")


class RegionBST(object):
    def __init__(self, knownGene_file, region_type):
        self.region_bst = {}
        for chrom, sorted_region_list \
                in get_sorted_disjoint_regions(knownGene_file, region_type).items():
            self.region_bst[chrom] = RegionNode(sorted_region_list)

    def overlaps(self, chrom, start, end):
        if chrom not in self.region_bst:
            return False
        else:
            return self.region_bst[chrom].overlaps(start, end)

    def overlapping_regions(self, chrom, start, end):
        if chrom not in self.region_bst:
            return []
        else:
            return self.region_bst[chrom].overlapping_regions(start, end)

    def overlapping_genes(self, chrom, start, end):
        if chrom not in self.region_bst:
            return []
        else:
            return self.region_bst[chrom].overlapping_genes(start, end)

    def dist_to_closest(self, chrom, start, end):
        if chrom not in self.region_bst:
            return -1
        else:
            dist = self.region_bst[chrom].overlaps(start, end)
            assert dist >= 0, "Negative distance on same chromosome reported"
            return dist

if __name__ == "__main__":
    print "Starting..."
    knownGene_file = '/home/hawkjo/genomes/hg19/knownGene.txt'
    start_time = time.time()
    cds_bst = RegionBST(knownGene_file, 'CDSs')
    print "Build time: %f sec" % (time.time() - start_time)

    if False:
        import pickle
        with open('cds_bst.pkl', 'w') as out:
            pickle.dump(cds_bst, out, pickle.HIGHEST_PROTOCOL)

    rand_cds = []
    for line in open('test_knownGene/rand_cds.txt'):
        chrom, start, end = line.strip().split()
        rand_cds.append((chrom, int(start), int(end)))
    rand_intron = []
    for line in open('test_knownGene/rand_introns.txt'):
        chrom, start, end = line.strip().split()
        rand_intron.append((chrom, int(start), int(end)))

    cds_found = 0
    start_time = time.time()
    for cds in rand_cds:
        if cds_bst.overlaps(*cds):
            cds_found += 1

    # 100 CDSs are tested, but 24 of them are for RNA genes, which should be ignored.
    print "%d of 76 CDSs found" % (cds_found)
    print "CDS testing time: %f sec" % (time.time() - start_time)

    start_time = time.time()
    knownGene_dict = knownGene_as_dict(knownGene_file)
    print "knownGene_as_dict time: %f sec" % (time.time() - start_time)

    start_time = time.time()
    intron_found = 0
    for intron in rand_intron:
        if cds_bst.overlaps(*intron):
            print
            print intron
            for gene in cds_bst.overlapping_genes(*intron):
                print knownGene_dict[gene]
            intron_found += 1
    print "%d of %d introns found" % (intron_found, len(rand_intron))
    print "Intron testing time: %f sec" % (time.time() - start_time)
