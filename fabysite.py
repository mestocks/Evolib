import sys

from evolib.SequenceFormats import FastaFormat

fileObject = open(sys.argv[1], 'r')

F = FastaFormat(fileObject)

for site in F.iter_sites():
    print site
