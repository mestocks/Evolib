###### ######

from evolib.SequenceFormats import FastaFormat
from evolib.NGSFormats import VariantCallFormat

from evolib.lib.DNAmethods import synNonsynProbs

from evolib.lib.FileIterators import gff_iter

###### ######

openRef = open("/home/mist/Documents/Projects/spruce/papo/manuscripts/20140107papo/data/nxtgen/blast_Pa.txt", 'r')
refdb = {}
Pa2Pg = {}
for line in openRef:
    info = line.rstrip().split()
    PgName = info[0]
    PaName = info[1]
    regions = map(int, info[2:])
    refdb.update({PgName: regions})
    Pa2Pg.update({PaName: PgName})

openGFF = open("/home/mist/Documents/Projects/spruce/papo/manuscripts/20140107papo/data/reference/High-Confidence.gff3", 'r')
myGFF = gff_iter(openGFF)

gffdb = {}
for i in myGFF:
    if i.name in refdb.keys():
        gffdb.update({i.name: i.exons})

annFile = open("/home/mist/Documents/Projects/spruce/papo/manuscripts/20140107papo/data/nxtgen/abies_ann.txt", 'r')

for line in annFile:
    line_split = line.rstrip().split()
    if len(line_split) == 1:
        locus, exons = line_split[0], []
    else:
        locus, exons = line_split[0], map(int, line_split[1:])
    openFile = open("/home/mist/Documents/Projects/spruce/papo/manuscripts/20140107papo/data/nxtgen/Pa_" + locus + ".fsa", 'r')
    myData = FastaFormat(openFile)
    myData.annotate(exons)
    myData.justSynonymous()
    print myData
    #print locus, myData.validSites, myData.seg_sites()




#print 'Syn', myData.validSites, myData.seg_sites(), myData.thetaW() / myData.validSites, myData.thetaPi() / myData.validSites, myData.tajD()

#coding_iter = myData.codingBySite()

#while True:
#    try:
#        Bp = coding_iter.next()
#        print Bp.alleles(), Bp.siteClassNarrow, Bp.siteClassBroad
#    except StopIteration:
#        break

