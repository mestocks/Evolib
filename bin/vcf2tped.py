import sys

from evolib.iterators import vcf_iter
from evolib.tools.GeneralMethods import create_ids

n = 1
for row in vcf_iter(sys.stdin):
    
    chrom, pos = row["CHROM"].replace("Contig",""), row["POS"]
    
    gpos, snpid = "0", "rs" + str(n)
    
    gts = (smp['GT'] for smp in row.iter_samples())
    alleles = []

    sep = ""
    gtstr = ""
    for gt in gts:
        if gt == "./.":
            tmpstr = sep + "0 0"
        else:
            gtsplit = [a + 1 for a in map(int, gt.split("/"))]
            alleles += gtsplit
            tmpstr = sep + ' '.join(map(str, gtsplit))

        gtstr += tmpstr
        sep = " "

    ualleles = set(alleles)

    if len(ualleles) > 2:
        pass
    elif alleles == []:
        print chrom, snpid, gpos, pos, gtstr
    else:
        if len(ualleles - set([1, 2])) == 0:
            print chrom, snpid, gpos, pos, gtstr
        else:
            allist = map(str, list(ualleles))
            gtstr = gtstr.replace(allist[0], '1')
            gtstr = gtstr.replace(allist[1], '2')
            print chrom, snpid, gpos, pos, gtstr

    n += 1
