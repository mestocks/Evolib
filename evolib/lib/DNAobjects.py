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
    """
    Class for dealing with the genotype values given in each row of 
    a VCF file. evolib.lib.VCFrow::ROW_BASECLASS provides a wrapper 
    for this.
    """
    def __init__(self, raw_genotypes, possible_alleles):
        self.raw_genotypes = raw_genotypes
        self.possible_alleles = possible_alleles
        
    def iter_benotypes(self):
        
        for Sample in raw_genotypes:
            b = Sample.binary_call()
            
            yield b[0] + b[1]

class Genotypes_old():
    
    
    def __init__(self, benotypes, possible_alleles):
        self.benotypes = benotypes
        self.possible_alleles = possible_alleles
    
    
    def __iter__(self):
        for g in self.iter_genotypes():
            yield g
            
    def __str__(self):
        return ' '.join([g[0] + "," + g[1] for g in self.iter_genotypes()])
        
    def allele_numbers(self, ualleles):
        
        alleles = ''.join([b[0] + b[1] for b in self.benotypes])
        nalleles = [alleles.count(a) for a in ualleles]
        
        return nalleles
    
    def iter_benotypes(self):
        
        for b in self.benotypes:
            yield (b[0], b[1])
    
    def iter_genotypes(self):
        
        for b in self.iter_benotypes():
            
            if b == ('.', '.'):
                yield ('N', 'N')
                
            else:
                one = self.possible_alleles[int(b[0])]
                two = self.possible_alleles[int(b[1])]
                
                yield (one, two)
    
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


