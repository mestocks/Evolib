## Object structure and inheritance


### Sequence data

It is often of interest to compare summary statistics calculated from observed data to the same statistics calculated from simulated datasets. The class ```SequenceData``` was therefore designed to be inherited by classes representing both real and simulated sequence data.
```
        SequenceData
        /          \
   FastaFormat    msFormat
```
This means that any statistical methods defined in ```SequenceData``` are available to both ```FastaFormat``` and ```msFormat```, reducing code redundancy and ensuring that any calculations are indentical between observed and simulated datasets. For example, in the following code the same ```.tajD()``` method is inherited and called by both the ```FastaFormat``` and ```msFormat``` classes:
```python
import sys
from evolib.SequenceFormats import FastaFormat
from evolib.SimulationMethods import ms_iter

# Observed data: print Tajima's D
fsaFile = open("myData.fsa", 'r')
Fsa = FastaFormat(fsaFile)
print Fsa.tajD()

# Simulated data: print Tajima's D
msData = sys.stdin
for ms in ms_iter(msData):
  print ms.tajD()
```
