import sys

from evolib.NGSFormats import VariantCallFormat

vcf = VariantCallFormat(sys.stdin)

args = sys.argv[1:]
names = args[0].split(",")

def gtiter(row, names):
    for samp in row.iter_samples():
        if samp.name in names:
            dp = samp['DP']
            if dp != "." and dp is not None and int(dp) > 7:
                gt = samp['GT']
                if gt is None or gt == "./.":
                    yield '0/0'
                else:
                    yield gt

for row in vcf:
    chrom, pos = row['CHROM'], row['POS']
    GTs = list(gtiter(row, names))
    #GTs = list((smp['GT'] for smp in row.iter_samples() if smp.name in names and smp['GT'] is not None and smp['GT'] != "./." and smp['DP'] != "." and smp['DP'] is not None and int(smp['DP']) > 7))
    
    ngts = len(GTs)
    nhet = GTs.count("0/1")

    if ngts == 0:
        print chrom, pos, "NA"
    else:
        print chrom, pos, nhet / float(ngts)
