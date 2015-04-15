import sys

args = sys.argv[1:]

from evolib.iterators import vcf_iter
from evolib.stats.PopGenStats import Dxy

if args[0] == "F_SI":
    str1 = ['f']
    str2 = ['s', 'i']


dataset = args[1]

#2014_SOL_TB8858__S05318_1  2014_SOL_TB8860__S05328_1  2014_SOL_TB8862__1509615_S  2014_SOL_TB8864__1560702_S  2014_SOL_TB8866__1458775_F  2014_SOL_TB8868__1415786_F

#2014_SOL_TB8859__S05325_1  2014_SOL_TB8861__S05331_1  2014_SOL_TB8863__1518438_S  2014_SOL_TB8865__1507057_S  2014_SOL_TB8867__1463052_F  2014_SOL_TB8869__1460666_F



if dataset == "wgs11":
    morph = {'2014_SOL_TB8858': 'i', '2014_SOL_TB8859': 'i',
             '2014_SOL_TB8860': 'i', '2014_SOL_TB8861': 'i',
             '2014_SOL_TB8862': 's', '2014_SOL_TB8863': 's',
             '2014_SOL_TB8864': 's', '2014_SOL_TB8865': 's',
             '2014_SOL_TB8866': 'f', '2014_SOL_TB8867': 'f',
             '2014_SOL_TB8868': 'f'}
elif dataset == "wgs":
    morph = {'2014_SOL_MB8684': 'i',
             '2014_SOL_MB8685': 's',
             '2014_SOL_MB8686': 's',
             '2014_SOL_MB8687': 'f',
             '2014_SOL_MB8688': 'f'}


#morph = dict([(m.split("\t")[0], m.split("\t")[4]) for m in open('/home/mist/current/ruffProject/pugnax_gwas/gwas/phe_r1bm_RUFFv0.2.txt', 'r')])



names1 = None
names2 = None

for row in vcf_iter(sys.stdin):
    
    if names1 is None:
        names1 = [value for value in row.header.names[9:] if morph[value] in str1]
        names2 = [value for value in row.header.names[9:] if morph[value] in str2]
        
    pop1 = list((row[ind1]['GT'] for ind1 in names1 if row[ind1]['GT'] != './.' and row[ind1]['DP'] > 7))
    pop2 = list((row[ind2]['GT'] for ind2 in names2 if row[ind2]['GT'] != './.' and row[ind2]['DP'] > 7))

    n1 = len(pop1)
    n2 = len(pop2)

    hom1 = pop1.count('0/0')
    het1 = pop1.count('0/1')
    alt1 = pop1.count('1/1')

    hom2 = pop2.count('0/0')
    het2 = pop2.count('0/1')
    alt2 = pop2.count('1/1')

    if n1 == len(names1) and n2 == len(names2) and (n1 + n2) == sum([hom1, het1, alt1, hom2, het2, alt2]):
        # fae/non-fae conditions
        # ((0/0 + 0/1 > 0) and (1/1 + 0/1 == 0)) or ((1/1 + 0/1 > 0) and (0/0 + 0/1 == 0))
        
        if (alt1 == 0 and (hom1 + het1 > 0) and (hom2 + het2 == 0)) or (hom1 == 0 and (alt1 + het1 > 0) and (alt2 + het2 == 0)):
            
            chrom, pos = row['CHROM'], row['POS']

            printlst = [chrom,
                        str(pos),
                        str(pos + 1),
                        str(2 * hom1 + het1),
                        str(2 * hom2 + het2),
                        str(2 * alt1 + het1),
                        str(2 * alt2 + het2)]

            #print hom1, het1, alt1
            #print hom2, het2, alt2
            #print row
            print '\t'.join(printlst)
            
    
