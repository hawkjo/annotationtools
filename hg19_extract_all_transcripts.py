import sys, knownGene
from Bio import SeqIO

print 'Reading in hg19...'
hg19_file = '/home/hawkjo/genomes/hg19/hg19.fa'
hg19 = {}
for record in SeqIO.parse(open(hg19_file), 'fasta'):
    hg19[record.name] = record

print 'Reading knownGene...'
knownGene_file = '/home/hawkjo/genomes/hg19/annotation/knownGene.txt'
knownGene_list = knownGene.knownGene_as_list(knownGene_file)

print 'Extracting transcripts...'
transcript_list = []
stop_codons = ['TAA','TAG','TGA']
genes_with_start_and_stop_codons = 0
genes_with_start_codons = 0
genes_with_stop_codons = 0
bad_genes = 0
for gene in knownGene_list:
    transcript = hg19['chr1'][0:0] # Get sequence record with empty sequence
    for exon in gene.exons:
        transcript += hg19[exon.chrom][exon.start:exon.end]
    if gene.strand == '-':
        transcript = transcript.reverse_complement()
    elif gene.strand != '+': sys.exit('Unexpected strand')
    transcript = transcript.upper()
    transcript.id = gene.name
    transcript.name = gene.name
    transcript.description = ''
    transcript_list.append(transcript)

    # Testing
    if gene.CDSs:
        start = 0
        for exon in gene.exons:
            if exon.start > gene.cdsStart:
                break
            elif exon.end < gene.cdsStart:
                start += exon.end - exon.start
            else:
                start += gene.cdsStart - exon.start

        end = len(transcript)
        for exon in reversed(gene.exons):
            if exon.end < gene.cdsEnd:
                break
            elif exon.start > gene.cdsEnd:
                end -= exon.end - exon.start
            else:
                end -= exon.end - gene.cdsEnd

        if gene.strand == '-':
            tmp = start
            start = len(transcript) - end
            end = len(transcript) - tmp

        if str(transcript.seq)[start:start+3] == 'ATG':
            has_start_codon = True
        else:
            has_start_codon = False
        if str(transcript.seq)[end-3:end] in stop_codons or \
                str(transcript.seq)[end:end+3] in stop_codons:
            has_stop_codon = True
        else:
            has_stop_codon = False

        if has_start_codon and has_stop_codon:
            genes_with_start_and_stop_codons += 1
        elif has_start_codon:
            genes_with_start_codons += 1
        elif has_stop_codon:
            genes_with_stop_codons += 1
        else:
            #print gene
            #print start, end, len(transcript)
            #SeqIO.write([transcript], sys.stdout, 'fasta')
            #raw_input('Press enter...')
            bad_genes += 1

print 'Writing output...'
outfile = 'hg19_transcripts.fa'
SeqIO.write(transcript_list, open(outfile,'w'), 'fasta')

print
print 'Summary:'
print 'genes_with_start_and_stop_codons:', genes_with_start_and_stop_codons 
print 'genes_with_start_codons:', genes_with_start_codons 
print 'genes_with_stop_codons:', genes_with_stop_codons 
print 'bad_genes:', bad_genes 
