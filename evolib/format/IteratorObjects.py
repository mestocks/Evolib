class Site():
    
    def __init__(self, alleles):
        self._alleles = alleles
        
    
    def __str__(self):
        return self.alleles()

    def alleles(self):
        return self._alleles
    
    
    def hasMissingData(self, dna = ['A', 'T', 'C', 'G']):
        
        if set(self.alleles()) <= set(dna):
            answer = False
        else:
            answer = True
            
        return answer
    
    def numberOfAlleles(self):
        alleles = set(self.alleles())
        return len(alleles)

###### ######

class FastqRead():
    def __init__(self, seq, seqid, qual):
        self.seq = seq
        self.seqid = seqid
        self.qual = qual
        
    def __str__(self):
        return '\n'.join([self.seqid, self.seq, '+', self.qual])
    
    def quality(self):
        qualities = map(ord, self.qual)
        return qualities

    def __len__(self):
        return len(self.seq)

###### ######

import numpy
from scipy import stats

#from evolib.lib.DataObjects import Site
from evolib.lib.DNAobjects import Genotypes
from evolib.lib.StatMethods import chisquared

class VCFrow(list, Site):
    
    def __init__(self, values, classes, header):
         
        self.NA = '.'
        self.classes = classes
        self.header = header
        self.values = values
        self.FormatClass = None
        self.lookupindex = dict([(key, index) for index, key in enumerate(self.header)])
        
        return list.__init__(self, values)
    
    def __getitem__(self, index):
        
        if isinstance(index, str):
            index = self.lookupindex[index]
            
        if isinstance(index, slice):
            value = []
            for i in range(len(self.header)).__getitem__(index):
                v = self._get_colvalue(i)
                value.append(v)
        else:
            value = self._get_colvalue(index)
        
        return value
    
    def _get_colvalue(self, index):
        
        if index > 8:
            if self.FormatClass is None:
                self.FormatClass = self.classes[8](self.values[8])
            value = self.classes[index](self.values[index], self.FormatClass.value)
        else:
            value = self.classes[index](self.values[index])
            
        return value
        
        
    def iter_samples(self):
        for Sample in self.__getitem__(slice(9, None)):
            yield Sample


###### ######

class GFFRecord():
    def __init__(self):
        self.exons = []
        self.name = None
        
    def add(self, items):
        if items[2] == "exon":
            self._exon(items)
        
        if self.name is None:
            self.name = items[0]
            
        try:
            self.records.append(items)
        except AttributeError:
            self.records = [items]
            
    def _exon(self, record):
        self.exons.extend([int(record[3]), int(record[4])])

###### ######

from evolib.lib.DataObjects import SequenceData

class msFormat(SequenceData):
    
    
    def __init__(self, text):
        self.text = text
        lines = [i for i in text.split('\n')[2:] if i != '']
        self.Seqs = lines
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



###### ######

from evolib.lib.DNAobjects import FastaSequence
from evolib.data.AlignmentObjects import DnaPopulationData
from evolib.lib.DataObjects import BinaryTable, SeqTable

class FastaAlignment(DnaPopulationData):
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
        
        self.DNAdata = SeqTable(seq_table)
        self.IOdata = self._get_IOdata(self.DNAdata)
        
        
    def length(self):
        return self.validSites
    
    
    def nsamples(self):
        return len(self.sequences)

