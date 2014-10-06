import random

# Classes
from lib.DNAobjects import FastaSequence
from lib.DataObjects import SequenceData
from lib.DataObjects import BinaryTable, SeqTable

###### ######

class FastaFormat(SequenceData):
    """
    *FastaFormat* - class representation of DNA sequence data in fasta format.
        evolib.SequenceFormats::FastaFormat
    
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

    def __init__(self, fileObject):
        self._fromFile(fileObject)
        
            
    def __getitem__(self, item):
        return FastaSequence(self.sequences[item], seqID = self.ids[item])

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        nseq = len(self.sequences)
        
        for i in range(nseq):
            sequence = FastaSequence(self.sequences[i], seqID = self.ids[i])
            
            yield sequence
            
            
    def __str__(self):
        
        stringList = []
        for s in range(len(self.sequences)):
            seq = FastaSequence(self.sequences[s], seqID = self.ids[s])
            stringList.append(seq)
        theString = '\n'.join(map(str, stringList))
            
        return theString
    

    def _fromFile(self, fileObject):
        """
        Returns a matrix where M[i][j] refers to the jth site of the ith individual.
        """        
        step1 = ''.join(map(lambda line: line, fileObject))
        step2 = step1.split('\n>')
        seq_table = [part.partition('\n')[2].replace('\n','') for part in step2]
        seq_names = [part.partition('\n')[0].replace('>', '') for part in step2]
        
        self.sequences = seq_table
        self.ids = seq_names
        
        self.Seqs = SeqTable(seq_table)
        self.IO = self._getBinaryTable(self.Seqs)
        
        
    def length(self):
        return self.validSites
    
    
    def nsamples(self):
        return len(self.sequences)

    
