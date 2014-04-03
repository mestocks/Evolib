import random

###### ######

class DNAsequence():
    
    def __init__(self, seq, seqID = None):
        
        if seqID is None:
            seqID = 'r' + str(random.randint(1, 1000000))
            
        self.name = seqID
        self.sequence = seq
        
    def __str__(self):
        return self.name, self.sequence
    
    def complement(self):
        
        compDict = {'A': 'T', 'a': 't',
                    'T': 'A', 't': 'a',
                    'G': 'C', 'g': 'c',
                    'C': 'G', 'c': 'g',
                    'N': 'N', 'n': 'n',
                    '-': '-'}
        
        self.sequence = ''.join([compDict[bp] for bp in self.sequence])
        
    def reverse(self):
        self.sequence = self.sequence[::-1]
        

###### ######

class FastaSequence(DNAsequence):
    
    def __str__(self):
        
        n = 70
        IDstring = '>' + self.name + '\n'
        sequence = '\n'.join([self.sequence[i: i + n] for i in range(0, len(self.sequence), n)])
        
        return IDstring + sequence
