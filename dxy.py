import sys

args = sys.argv[1:]

ind5 = ['2014_SOL_MB8684']
sat5 = ['2014_SOL_MB8685', '2014_SOL_MB8686']
fae5 = ['2014_SOL_MB8687', '2014_SOL_MB8688']

ind12 = ['2014_SOL_TB8858', '2014_SOL_TB8859',
         '2014_SOL_TB8860', '2014_SOL_TB8861']

sat12 = ['2014_SOL_TB8862', '2014_SOL_TB8863',
         '2014_SOL_TB8864', '2014_SOL_TB8865']

fae12 =  ['2014_SOL_TB8866', '2014_SOL_TB8867',
          '2014_SOL_TB8868', '2014_SOL_TB8869']

ind17 = ind5 + ind12
sat17 = sat5 + sat12
fae17 = fae5 + fae12

if args[1] == "wgs5":
    if args[0] == "F_I":
        ref = ['2014_SOL_MB8684']
        alt = ['2014_SOL_MB8687', '2014_SOL_MB8688']
    elif args[0] == "F_S":
        ref = ['2014_SOL_MB8685', '2014_SOL_MB8686']
        alt = ['2014_SOL_MB8687', '2014_SOL_MB8688']
    elif args[0] == "S_I":
        ref = ['2014_SOL_MB8684']
        alt = ['2014_SOL_MB8685', '2014_SOL_MB8686']
    elif args[0] == "F_SI":
        ref = ['2014_SOL_MB8684', '2014_SOL_MB8685', '2014_SOL_MB8686']
        alt = ['2014_SOL_MB8687', '2014_SOL_MB8688']
elif args[1] == "wgs17":
    if args[0] == "F_I":
        ref = ind17
        alt = fae17
    elif args[0] == "F_S":
        ref = sat17
        alt = fae17
    elif args[0] == "S_I":
        ref = ind17
        alt = sat17
    elif args[0] == "F_SI":
        ref = sat17 + ind17
        alt = fae17
    
    
#print ref, alt

#ref = ['2014_SOL_MB8684']
#ref = ['2014_SOL_MB8685', '2014_SOL_MB8686']
#alt = ['2014_SOL_MB8685', '2014_SOL_MB8686']
#alt = ['2014_SOL_MB8687', '2014_SOL_MB8688']
#names = ['2014_SOL_MB8685', '2014_SOL_MB8686']

nref = len(ref)
nalt = len(alt)

from evolib.iterators import vcf_iter
from evolib.stats.PopGenStats import Dxy

for row in vcf_iter(sys.stdin):

    refgts = list((r['GT'] for r in row.iter_samples() if r.name in ref and r['DP'] > 7 and r['GT'] != './.' and r['GT'] is not None))
    
    altgts = list((a['GT'] for a in row.iter_samples() if a.name in alt and a['DP'] > 7 and a['GT'] != './.' and a['GT'] is not None))

    nsamref = len(refgts)
    refref = refgts.count('0/0')
    refhet = refgts.count('0/1')
    refalt = refgts.count('1/1')
    
    nsamalt = len(altgts)
    altref = altgts.count('0/0')
    althet = altgts.count('0/1')
    altalt = altgts.count('1/1')

    chrom, pos = row['CHROM'], row['POS']

    # ind conditions
    if nsamref == nref and nsamalt == nalt and nsamref == sum([refref, refhet, refalt]) and nsamalt == sum([altref, althet, altalt]):
        
        printlst = [chrom,
                    str(pos),
                    str(2 * refref + refhet),
                    str(2 * refalt + refhet),
                    str(2 * altref + althet),
                    str(2 * altalt + althet)]
        #print row
        dxy = Dxy(2 * refref + refhet,
                  2 * refalt + refhet,
                  2 * altref + althet,
                  2 * altalt + althet)
                  
        print ' '.join(printlst + [str(dxy)])

    else:

        print chrom, pos, "NA", "NA", "NA", "NA", "NA"
