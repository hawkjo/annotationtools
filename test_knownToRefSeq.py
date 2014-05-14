import sys, os
from itertools import izip
from collections import defaultdict

NCBI_DATA_dir  = '/home/hawkjo/scratch/dbs/NCBI/DATA/'
gene2refseq_file = os.path.join(NCBI_DATA_dir, 'gene2refseq') 

print 'Reading field names...'
with open(gene2refseq_file) as f:
    line = f.readline().strip()[9:334]
    gene2refseq_fields = line.split()

print 'Reading gene2refseq...'
field_given_value = {}
for line in open(gene2refseq_file):
    values = line.strip().split('\t')
    for field_name, value in izip( gene2refseq_fields, values ):
        field_given_value[value.split('.')[0]] = field_name

print 'Reading kg2refseq...'
kg2rs_file = '/home/hawkjo/scratch/genomes/hg19/annotation/knownToRefSeq.txt'
refseqs_in_kg = [line.strip().split()[1] for line in open(kg2rs_file)]

print 'Counting...'
counts = defaultdict(int)
failures = 0
print 'Failures:'
for rs in refseqs_in_kg:
    try:
        counts[ field_given_value[rs] ] += 1
    except KeyError:
        print rs
        failures += 1

print
print '%7d  KeyErrors' % failures
for k,v in sorted(counts.items(), key=lambda tup:-tup[1]):
    print '%7d  %s' % (v, k)
