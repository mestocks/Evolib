import sys

args = sys.argv[1:]

from evolib.iterators import vcf_iter
from evolib.formats.config_parse import idsets

key = args[0]
filename = args[1]

rdict = idsets(filename)
inds = rdict[key]

for row in vcf_iter:
    
