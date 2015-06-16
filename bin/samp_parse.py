"""
cat <aln.vcf> | python samp_parse.py [ <INFO>,<INFO> ] <FORMAT>
"""
import sys

args = sys.argv[1:]

if len(args) == 1:
    info_arg = None
    format_arg = args[0]
elif len(args) > 1:
    info_arg = args[0].split(",")
    format_arg = args[1]
else:
    raise ValueError, "No arguments"
    

from evolib.iterators import vcf_iter

for row in vcf_iter(sys.stdin):
    chrom, pos = row['CHROM'], row['POS']
    print chrom, pos,
    
    if info_arg is not None:
        for i in info_arg:
            info = row['INFO'][i]
            print info,

    frmt = (d[format_arg] for d in row.iter_samples())
    
    for f in frmt:
        print f,
    print ""
