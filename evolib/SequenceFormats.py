#import random

# Classes
#from lib.DNAobjects import FastaSequence
#from lib.DataObjects import SequenceData
#from lib.DataObjects import BinaryTable, SeqTable

###### ######

from evolib.format.IteratorObjects import FastaAlignment

class FastaFormat(FastaAlignment):
    """
    Example usage:
    
       >>> from evolib.SequenceFormats import FastaFormat
       
       Example 1:
       >>> fileObject = open('example.fsa', 'r')
       >>> F = FastaFormat(fileObject)
       
       Example 2:
       >>> import sys
       >>> stdinObject = sys.stdin
       >>> F = FastaFormat(stdinObject)
       
    """
    pass
