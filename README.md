Evolib
======

**Evolib** is a python library for analysing DNA sequence data. It consists of a number of objects that abstract common sequence formats so that cleaner, simpler code can written.

Quick guide
-------------

VCF files can be efficiently iterated:
```python
import sys

from evolib.NGSFormats import VariantCallFormat

myVCF = VariantCallFormat(sys.stdin)

for row in myVCF:
    print row['CHROM'], row['POS']
```

Similarily, for data stored in fasta format:
```python
from evolib.SequenceFormats import FastaFormat

f = open("myData.fsa", 'r')
F = FastaFormat(f)

print F.seg_sites(), F.tajD()
```
