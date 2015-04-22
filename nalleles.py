import sys

#names = ['2014_SOL_MB8687', '2014_SOL_MB8688']
names = ['2014_SOL_MB8685', '2014_SOL_MB8686']

from evolib.iterators import vcf_iter

for row in vcf_iter(sys.stdin):
    
    gts = list((smp['GT'] for smp in row.iter_samples() if smp.name in names and smp['DP'] > 7 and smp['GT'] != './.' and smp['GT'] is not None))

    nsam = len(gts)
    ref = gts.count('0/0')
    het = gts.count('0/1')
    alt = gts.count('1/1')

    # ind conditions
    if nsam == 2 and nsam == sum([ref, het, alt]):
        
        chrom, pos = row['CHROM'], row['POS']
        
        printlst = [chrom,
                    str(pos),
                    str(2 * ref + het),
                    str(2 * alt + het)]
        
        print ' '.join(printlst)
        
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
        

    
