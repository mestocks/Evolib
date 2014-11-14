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

from evolib.generic.AlignmentSite import VCFSite

class VCFrow(VCFSite):
    
    def __init__(self, values, classes, header, nsamples):
        
        self.NA = '.'
        self.classes = classes
        self.header = header
        self.values = values
        self.nsamples = nsamples
        
        self.FormatClass = None
        self.lookupindex = None
        
    
    def __getitem__(self, index):

        if self.lookupindex is None:
            self.lookupindex = dict([(key, index) for index, key in enumerate(self.header)])
        
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

    def __str__(self):
        return '\t'.join(self.values)
    
    def _get_colvalue(self, index):
        
        if index > 8:
            if self.FormatClass is None:
                self.FormatClass = self.classes[8](self.values[8])
            value = self.classes[index](self.values[index], self.FormatClass.value)
        else:
            value = self.classes[index](self.values[index])
            
        return value
        
        
    def iter_samples(self):
        
        self.FormatClass = self.classes[8](self.values[8])
        maxRange = self.nsamples + 9
        
        FormatValue = self.FormatClass.value
        SampleClass = self.classes[9]
        
        return (SampleClass(self.values[i], FormatValue) for i in xrange(9, maxRange))


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

from evolib.data.AlignmentObjects import IOPopulationData

class msFormat(IOPopulationData):
    
    
    def __init__(self, text):
        self.text = text
        lines = [i for i in text.split('\n')[2:] if i != '']
        self.Seqs = lines
        #self.IOdata = self._alt_get_IOdata(lines)
        self.IOdata = self._get_IOdata(lines)
        
            
    def __str__(self):
        return self.text
    
    #def nsamples(self):
    #    
    #    try:
    #        n = len(self.IOdata[0])
    #    except IndexError:
    #        n = None
    #        
    #    return n
    
    def sample_sites(self, p):
        
        newIO = BinaryTable()
        for i in self.IOdata:
            pick = random.random()
            if pick < p:
                newIO.append(i)
                
        self.IOdata = newIO



###### ######

from evolib.generic.AlignmentSite import FastaSite
from evolib.generic.GeneticSequence import FastaSequence
from evolib.data.AlignmentObjects import DnaPopulationData
from evolib.data.DataObjects import BinaryTable, SeqTable

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
        
        for i in self.iter_seqs():
            yield i

    def __str__(self):

        stringlist = (str(fseq) for fseq in self.iter_seqs())

        return '\n'.join(stringlist)


    def iter_seqs(self):
        
        nseq = len(self.sequences)
        
        for i in range(nseq):
            sequence = FastaSequence(self.sequences[i], seqID = self.ids[i])
            
            yield sequence 
    

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
        self.IOdata = self._alt_get_IOdata(self.DNAdata)
        
        
    def length(self):
        return self.validSites
    
    
    #def nsamples(self):
    #    return len(self.sequences)

