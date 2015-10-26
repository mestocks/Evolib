#! /usr/bin/python2

"""
...| vcfINFO <col1>,<col2>
"""

import sys

args = sys.argv[1:]

from evolib.iterators import vcf_iter

colstr = args[0]
cols = colstr.split(",")

print "CHROM", "POS", ' '.join(cols)

for row in vcf_iter(sys.stdin):
    chrom, pos = row["CHROM"], row["POS"]

    print chrom, pos,

    for col in cols:
        colobj = row['INFO'][col]
        if colobj is None:
            print "NA",
        else:
            print colobj,
    print "\n",

    

    
