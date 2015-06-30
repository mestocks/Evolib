import sys

from evolib.NGSFormats import VariantCallFormat

vcf = VariantCallFormat(sys.stdin)

for row in vcf:
    print row
