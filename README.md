Evolib
======

**Evolib** is a python library for analysing DNA sequence data. It consists of a number of objects that abstract common sequence formats so that cleaner, simpler code can written.

Quick guide
-------------

VCF files can be efficiently iterated:
```python
from evolib.NGSFormats import VariantCallFormat

openFile = open("myData.vcf", 'r')
myVCF = VariantCallFormat(openFile)

for bp in myVCF:
  print bp['CHROM'], bp['POS'], bp.genotypes(), bp.heterozygosity()
```

Similarily, for data stored in fasta format:
```python
from evolib.SequenceFormats import FastaFormat

f = open("myData.fsa", 'r')
F = FastaFormat(f)

print F.seg_sites(), F.tajD()
```

Populations can be easily defined and statistics calculated:
```python
F = FastaFormat(f)
F.define_pops([12, 12])

for pop in F:
  print pop.seg_sites(), pop.tajD()
```
