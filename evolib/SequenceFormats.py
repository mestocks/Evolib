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
    
    On the creation of a FastaFormat instance, two attributes are created from 
    the fasta file and assigned to the object:
    
        FastaFormat.Sequences - is an object of type `lib.DNAmethods::SeqTable' 
            and 
        
        FastaFormat.IOtable - is an object of type `lib.DNAmethods::IOPolyTable' 
            and 
    
    
    
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

###### ######

class msFormat(SequenceData):
    
    
    def __init__(self, text):
        self.text = text
        lines = [i for i in text.split('\n')[2:] if i != '']
        self.IO = self._getBinaryTable(lines)
        
            
    def __str__(self):
        return self.text
    
    
    def _getBinaryTable(self, seqs):
        
        IO = BinaryTable()
        for line in seqs:
            IO.add_sample(line)
            
        return IO
            
    
    def nsamples(self):
        
        try:
            n = len(self.IO[0])
        except IndexError:
            n = None
            
        return n
    
    def sample_sites(self, p):
        
        newIO = BinaryTable()
        for i in self.IO:
            pick = random.random()
            if pick < p:
                newIO.append(i)
                
        self.IO = newIO
