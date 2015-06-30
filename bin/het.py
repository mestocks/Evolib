import sys

from evolib.NGSFormats import VariantCallFormat

vcf = VariantCallFormat(sys.stdin)

args = sys.argv[1:]
names = args[0].split(",")

for row in vcf:
    chrom, pos = row['CHROM'], row['POS']
    try:
        GTs = list((smp['GT'] for smp in row.iter_samples() if smp.name in names and smp['GT'] != "./." and smp['DP'] != "." and smp['DP'] is not None and int(smp['DP']) > 7))
    except KeyError:
        GTs = []

    ngts = len(GTs)
    nhet = GTs.count("0/1")

    if ngts == 0:
        print chrom, pos, "NA"
    else:
        print chrom, pos, nhet / float(ngts)
