Evolib
======

**Evolib** is a python library for analysing DNA sequence data. It consists of a number of objects that abstract common sequence formats so that cleaner, simpler code can written.

Note that the library should be considered an early alpha release. 

Quick guide
-------------

Iterating VCF files:
```python
import sys

from evolib.NGSFormats import VariantCallFormat

myVCF = VariantCallFormat(sys.stdin)

for row in myVCF:
    chrom = row['CHROM']
    pos = int(row['POS'])
    dps = map(int, (smp['DP'] for smp in row.iter_samples()))
    print chrom, pos, len((dp for dp in dps if dp > 8))
    
```

Similarily, for alignments stored in fasta format:
```python
from evolib.SequenceFormats import FastaFormat

f = open("myData.fsa", 'r')
F = FastaFormat(f)

print F.seg_sites(), F.tajD()
```
