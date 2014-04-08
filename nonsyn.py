from evolib.lib.DNAmethods import sites2codons
from evolib.lib.GeneralMethods import block_iter

from evolib.NGSFormats import VariantCallFormat

###### ######

myVCF = VariantCallFormat("test/pavy2012_A.vcf")

for bp in myVCF.codingBySite():
    print bp.siteClassNarrow, bp.siteClassBroad

"""
locus = 


"""
