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
