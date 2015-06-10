import sys

from evolib.iterators import vcf_iter

for row in vcf_iter(sys.stdin):

    chrom, pos = row["CHROM"].replace("Contig",""), row["POS"]
    gpos, snpid = "0", "rs" + str(pos)
    
    gts = (smp['GT'] for smp in row.iter_samples())

    sep = ""
    gtstr = ""
    for gt in gts:
        if gt == "./.":
            tmpstr = sep + "0 0"
        else:
            gtsplit = (a + 1 for a in map(int, gt.split("/")))
            tmpstr = sep + ' '.join(map(str, gtsplit))

        gtstr += tmpstr
        sep = " "
        
    print chrom, snpid, gpos, pos, gtstr
