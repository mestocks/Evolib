###### ######

from evolib.NGSFormats import VariantCallFormat

###### ######

myVCF = VariantCallFormat("test/pavy2012_A.vcf")

for bp in myVCF.codingBySite():
    print bp.siteClassNarrow, bp.siteClassBroad

