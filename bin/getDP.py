import sys

from evolib.iterators import vcf_iter

for row in vcf_iter(sys.stdin):
            
    dp = map(str, (smp['DP'] for smp in row.iter_samples()))

    chrom, pos = row['CHROM'], row['POS']

    print ' '.join([chrom, str(pos)]), ' '.join(dp) 
        
        #for i in pop1:
        #    printstr += i,
        #print hom1, het1, alt1,
        #print "|",
        #for i in pop2:
        #    print i,
        #print hom2, het2, alt2,
        
        #print ""
    #alleles = row.alleles()
    #print chrom, pos, alleles
    #pop1 = ''.join([alleles[2 * f] + alleles[(2 * f) + 1] for f in Findices if int(row[f + 9]['DP']) > 7]).replace('N','')
    #pop2 = ''.join([alleles[2 * i] + alleles[(2 * i) + 1] for i in Iindices if int(row[i + 9]['DP']) > 7]).replace('N','')
    
    #uall = list(set(pop1 + pop2))

    #if len(pop1) > 0 and len(pop2) > 0:
    #    if len(uall) == 2:
    #        a = uall[0]
    #        b = uall[1]
    #        xa = pop1.count(a)
    #        xb = pop1.count(b)
    #        ya = pop2.count(a)
    #        yb = pop2.count(b)
    #        dxy = Dxy(xa, xb, ya, yb)
    #        print chrom, pos, dxy
        

    
