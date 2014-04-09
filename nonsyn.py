###### ######

from evolib.SequenceFormats import FastaFormat
from evolib.NGSFormats import VariantCallFormat

from evolib.lib.DNAmethods import synNonsynProbs

###### ######

#myData = VariantCallFormat("test/pavy2012_A.vcf")
openFile = open("test/Pa_can008.fsa", 'r')
myData = FastaFormat(openFile)

myData.annotate([1, 1106])
myData.justSynonymous()

print 'Syn', myData.validSites, myData.seg_sites(), myData.thetaW() / myData.validSites, myData.thetaPi() / myData.validSites, myData.tajD()

#coding_iter = myData.codingBySite()

#while True:
#    try:
#        Bp = coding_iter.next()
#        print Bp.alleles(), Bp.siteClassNarrow, Bp.siteClassBroad
#    except StopIteration:
#        break

