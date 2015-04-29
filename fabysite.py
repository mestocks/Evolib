import sys

from evolib.SequenceFormats import FastaFormat

fileName = sys.argv[1]
fileObject = open(fileName, 'r')

F = FastaFormat(fileObject)

chrom = "Contig" + fileName[2:6]

dna = set(['A', 'T', 'G', 'C'])
bp = 1
for site in F.iter_sites():
    if len(list(set(site.upper()) - dna)) > 0:
        print chrom, bp, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
    else:
        ref = site[0].upper()
        print chrom, bp, ref,
        for i in site[1:]:
            if i.upper() == ref:
                print '0',
            else:
                print '1',
        print ""
    bp += 1
