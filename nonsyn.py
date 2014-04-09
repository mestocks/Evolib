###### ######

from evolib.SequenceFormats import FastaFormat
from evolib.NGSFormats import VariantCallFormat

from evolib.lib.DNAmethods import synNonsynProbs

###### ######

#myData = VariantCallFormat("test/pavy2012_A.vcf")
openFile = open("test/Pa_can008.fsa", 'r')
myData = FastaFormat(openFile)

print 'All', myData.validSites, myData.seg_sites(), myData.thetaW() / myData.validSites, myData.thetaPi() / myData.validSites, myData.tajD()

myData.justSynonymous()

print 'Syn', myData.valid_syn, myData.seg_sites(), myData.thetaW() / myData.valid_syn, myData.thetaPi() / myData.valid_syn, myData.tajD()

#coding_iter = myData.codingBySite()

#while True:
#    try:
#        Bp = coding_iter.next()
#        print Bp.alleles(), Bp.siteClassNarrow, Bp.siteClassBroad
#    except StopIteration:
#        break

