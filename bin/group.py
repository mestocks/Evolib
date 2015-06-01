"""
... | python casecontrol.py <key> <file>

where <file> takes the form:

<key>:case:<id1>,<id2>,...
<key>:control:<id3>,<id4>,...

"""

import sys

args = sys.argv[1:]

from evolib.iterators import vcf_iter

key = args[0]
filename = args[1]

names1 = []
names2 = []

for k in open(filename, 'r'):
    ksplit = k.rstrip().split(":")
    if ksplit[0] == key:
        if ksplit[1] == "case":
            cases += ksplit[2].split(",")
        elif ksplit[1] == "control":
            controls += ksplit[2].split(",")

for row in vcf_iter(sys.stdin):
            
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
        
        chrom, pos = row['CHROM'], row['POS']
        # case1, case2, control1, control2
        printlst = [chrom,
                    str(pos),
                    str(2 * hom1 + het1),
                    str(2 * alt1 + het1),
                    str(2 * hom2 + het2),
                    str(2 * alt2 + het2)]

        print '\t'.join(printlst)
            
    
