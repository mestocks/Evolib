import random

###### ######

class DNAsequence():
    
    def __init__(self, seq, seqID = None):
        
        if seqID is None:
            seqID = 'r' + str(random.randint(1, 1000000))
            
        self.name = seqID
        self.sequence = seq
        
    def __getitem__(self, index):
        return self.sequence[index]
        
    def __len__(self):
        return len(self.sequence)
        
    def __str__(self):
        return self.sequence
    
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

###### ######

class Genotypes():
    
    
    def __init__(self, benotypes, ref, alt):
        self.benotypes = benotypes
        self.ref, self.alt = ref, alt
    
        
    def allele_numbers(self, ualleles):
        
        alleles = ''.join([b[0] + b[1] for b in self.benotypes])
        nalleles = [alleles.count(a) for a in ualleles]
        
        return nalleles
    
    
    def number_of_alleles(self):
        return len(set(''.join([b[0] + b[1] for b in self.benotypes])))
    
    
    def number_of_genotypes(self):
        return len(set([b[0] + b[1] for b in self.benotypes]))
    
    
    def subset(self, indices):
        self.benotypes = [self.benotypes[b] for b in indices]
        
        
    def unique_alleles(self):
        
        ualleles =  list(set(''.join([b[0] + b[1] for b in self.benotypes])))
        ualleles.sort()
        
        return ualleles
