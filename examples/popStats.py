import sys

from evolib.SequenceFormats import FastaFormat

files = sys.argv[1:]

for f in files:
    file = open(f, 'r')
    F = FastaFormat(file)
    print f, F.nsamples(), F.length(), F.seg_sites(), F.thetaW() / F.length(), F.thetaPi() / F.length(), F.tajD()
