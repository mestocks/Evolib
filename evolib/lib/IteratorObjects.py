class FastqRead():
    def __init__(self, seq, seqid, quality):
        self.seq = seq
        self.seqid = seqid
        self.quality = quality
        
    def __str__(self):
        return '\n'.join([self.seqid, self.seq, '+', self.quality])

###### ######

class Site():
    
    def __init__(self, alleles):
        self._alleles = alleles
        
    
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
