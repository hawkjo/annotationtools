import sys,time
from collections import defaultdict

knownGene_file = '/home/hawkjo/genomes/hg19/knownGene.txt'

class Gene(object):
    def __init__(self,line=None):
        if line:
            self._parse(line)
        else:
            sys.exit("No line specified for gene.")

    def _parse(self,line):
        var = line.strip().split('\t')
        self.name = var[0]
        self.chrom = var[1]
        self.strand = var[2]
        self.txStart = int(var[3])
        self.txEnd = int(var[4])
        self.cdsStart = int(var[5])
        self.cdsEnd = int(var[6])
        self.exonCount = int(var[7])
        self.exonStarts = [int(s) for s in var[8].split(',')[:self.exonCount]]
        self.exonEnds = [int(s) for s in var[9].split(',')[:self.exonCount]]
        self.exons = []
        self.introns = []
        self.CDSs = []
        self.UTRs = []
        for start, end in zip(self.exonStarts, self.exonEnds):
            self.exons.append( Exon(self.name, self.chrom, start, end ) )

            # Now to break the exons up into coding and non-coding regions
            if self.cdsStart == self.cdsEnd:
                # Non-coding gene. All exons are UTRs
                self.UTRs.append( Exon(self.name, self.chrom, start, end ) )
            elif start >= self.cdsStart and end <= self.cdsEnd:
                # Exon entirely within coding region
                self.CDSs.append( Exon(self.name, self.chrom, start, end ) )
            elif start < self.cdsStart and self.cdsStart < end <= self.cdsEnd:
                # Exon has coding region and 5' UTR
                self.UTRs.append( Exon(self.name, self.chrom, start, self.cdsStart ) )
                self.CDSs.append( Exon(self.name, self.chrom, self.cdsStart, end ) )
            elif self.cdsStart <= start < self.cdsEnd and end > self.cdsEnd:
                # Exon has coding reigon and 3' UTR
                self.CDSs.append( Exon(self.name, self.chrom, start, self.cdsEnd ) )
                self.UTRs.append( Exon(self.name, self.chrom, self.cdsEnd, end ) )
            elif start < self.cdsStart < self.cdsEnd < end:
                # Exon has coding region and both 5' and 3' UTRs
                self.UTRs.append( Exon(self.name, self.chrom, start, self.cdsStart ) )
                self.CDSs.append( Exon(self.name, self.chrom, self.cdsStart, self.cdsEnd ) )
                self.UTRs.append( Exon(self.name, self.chrom, self.cdsEnd, end ) )
            else:
                # Coding region exists but does not overlap exon
                self.UTRs.append( Exon(self.name, self.chrom, start, end ) )

    def __str__(self):
        return "%s %s %s %d %d %d %d %d %s %s" % \
                (   self.name, self.chrom, self.strand, self.txStart, self.txEnd,
                    self.cdsStart, self.cdsEnd, self.exonCount,
                    ','.join(str(start) for start in self.exonStarts),
                    ','.join(str(end) for end in self.exonEnds) )

def knownGene_as_dict():
    genes = {}
    for line in open(knownGene_file):
        gene = Gene(line)
        genes[gene.name] = gene
    return genes

def knownGene_as_list():
    return [Gene(line) for line in open(knownGene_file)]

class Exon(object):
    """Simple exon object with geneName, chrom, start, and end."""
    def __init__(self, gene_name, chrom, start, end):
        self.geneName = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end

def get_all_sorted_exons():
    exons = defaultdict(list)
    for gene in knownGene_as_list():
        exons[gene.chrom].extend(gene.exons)
    for chrom in exons.keys():
        exons[chrom].sort(key = lambda exon: exon.start)
    return exons

class ExonRegion(object):
    """A contiguous region where every base is covered by an exon"""
    def __init__(self, exon):
        self.geneNames = set(exon.geneName)
        self.chrom = exon.chrom
        self.start = exon.start
        self.end = exon.end

    def add(self, exon):
        if self.chrom != exon.chrom: return 
        self.geneNames.add( exon.geneName )
        self.start = min( self.start, exon.start )
        self.end = max( self.end, exon.end )

    def should_add(self, exon):
        if self.chrom != exon.chrom:
            return False
        elif self.start <= exon.start <= self.end or \
                self.start <= exon.end <= self.end or \
                exon.start <= self.start <= exon.end:
            # Definition of should_add here includes the possibility of two exons
            # lying immediately next to each other but not having any bases in common.
            # The math above is due to the half-open convention of start and end.
            return True
        else:
            return False

def get_sorted_disjoint_exon_regions():
    exon_regions = defaultdict(list)
    for chrom, exons in get_all_sorted_exons().items():
        for exon in exons:
            if exon_regions[chrom] and exon_regions[chrom][-1].should_add( exon ):
                exon_regions[chrom][-1].add( exon )
            else:
                exon_regions[chrom].append( ExonRegion(exon) )
    return exon_regions

class ExonRegionNode(object):
    def __init__(self, sorted_exon_regions_list):
        mid = int(len(sorted_exon_regions_list)/2)
        exon_region = sorted_exon_regions_list[mid]
        self.geneNames = exon_region.geneNames
        self.chrom = exon_region.chrom
        self.start = exon_region.start
        self.end = exon_region.end
        if len( sorted_exon_regions_list[:mid] ) > 0:
            self.left = ExonRegionNode( sorted_exon_regions_list[:mid] )
        else: self.left = None
        if len( sorted_exon_regions_list[mid+1:] ) > 0:
            self.right = ExonRegionNode( sorted_exon_regions_list[mid+1:] )
        else: self.right = None

    def local_overlap(self, start, end):
        if self.start <= start < self.end or \
                self.start < end <= self.end or \
                start <= self.start < end:
            # Convention here says that start and end are half-open as well
            return True
        else:
            return False

    def overlaps(self, start, end):
        if self.local_overlap( start, end ):
            return True
        elif self.start > end:
            if self.left: return self.left.overlaps( start, end )
            else: return False
        elif self.end < start:
            if self.right: return self.right.overlaps( start, end )
            else: return False
        else:
            sys.exit("overlaps error. Impossible overlap comparison.")

    def overlapping_regions(self, start, end):
        if self.local_overlap( start, end ):
            overlapping_regions = []
            if self.left:
                overlapping_regions.extend( self.left.overlaps( start, end ) )
            overlapping_regions.append( (self.start, self.end) )
            if self.right:
                overlapping_regions.extend( self.right.overlaps( start, end ) )
            return overlapping_regions
        elif self.start > end:
            if self.left: return self.left.overlaps( start, end )
            else: return [] # Equivalent to False, but compatible with list extension
        elif self.end < start:
            if self.right: return self.right.overlaps( start, end )
            else: return [] # Equivalent to False, but compatible with list extension
        else:
            sys.exit("overlaps error. Impossible overlap comparison.")

    def inorderwalk(self):
        if self.left:
            for x in self.left.inorderwalk(): yield x
        yield self
        if self.right:
            for x in self.right.inorderwalk(): yield x

class ExonRegionBST(object):
    def __init__(self):
        self.exon_region_bst = {}
        for chrom, sorted_exon_region_list in get_sorted_disjoint_exon_regions().items():
            self.exon_region_bst[chrom] = ExonRegionNode( sorted_exon_region_list )

    def overlaps(self, chrom, start, end):
        if chrom not in self.exon_region_bst: return False
        else: return self.exon_region_bst[chrom].overlaps( start, end )

if __name__ == "__main__":
    print "Starting..."
    start_time = time.time()
    exon_bst = ExonRegionBST()
    print "Build time: %f sec" % (time.time() - start_time)
    print exon_bst.exon_region_bst.keys()
    print exon_bst.overlaps( 'chr1', 2000000,20000000)
    print exon_bst.overlaps( 'chr1', 2000,2100)
