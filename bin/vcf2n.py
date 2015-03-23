import sys
import time

sys.path.append('/home/bo1mesx/lib/Evolib')

from evolib.iterators import vcf_iter

######

#args = sys.argv[1:]

#subsample_file = open(args[0], 'r')
#subsamples = [s.rstrip() for s in subsample_file]

#print subsamples
######

#chr_num = [1270, 2063, 3047, 3208, 3357, 3750, 3913, 4111, 4944, 9013]
#chr_num = [0]
#inc_chrom = ['Contig' + str(i) for i in chr_num]

######
# list of morph indices
#faeder = []
#satellite = []
#indepedent = []

#class TestA():
#    def __init__(self, init_string):
#        self.istring = init_string
        
#class TestB():
#    pass

#startA = time.time()

#for i in xrange(10000):
#    clA = TestA(str(i))

#print time.time() - startA
    
#startB = time.time()

#clB = TestB()
#for j in xrange(10000):
#    clB.init_string = str(j)
    
#print time.time() - startB
######
# iter population sample

# need direct access to name of each sample from 

######

vcfstream = vcf_iter(sys.stdin)
#vcfstream = vcf_iter4(sys.stdin)

#for line in sys.stdin:
#    #if line.startswith('#') is False:
#    if line[0] != '#':
#        values = line.rstrip().split('\t')
#        chrom = values[0]
#    elif line[:2] == '##':
#        pass
#    else:
#        pass

for row in vcfstream:
    #pass
    chrom, pos = row['CHROM'], row['POS']
    
    gts = (smp['GT'].split('/') for smp in row.iter_samples() if smp['GT'] != './.')
    #print ' '.join(map(str, gts))
    alls = [int(gt[0]) + int(gt[1]) for gt in gts]
    
    nsam = 2 * sum(1 for g in alls)
    nder = sum(alls)
    print ' '.join(map(str, alls)), nsam, nder
    if nsam == 0:
        print chrom, pos, 'NA', 'NA'
    else:
        print chrom, pos, nsam - nder, nder
    #print row
    #smp = row.smap(str, 'DP')
    #print chrom, pos, ' '.join(smp)
    #chrom, pos = row['CHROM'], row['POS']
    #print chrom, pos
#    if faeder == []:
        # add sample names from header to faeder, satellite and independent
    
#    CHROM = row['CHROM']
#    if CHROM in inc_chrom:
        
        
