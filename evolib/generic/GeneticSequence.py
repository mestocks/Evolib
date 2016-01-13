import random

from evolib.tools.GeneralMethods import block_iter, non_overlapping_iter
from evolib.tools.DNAmethods import complement, reverse, dna_to_amino

###### ######

class GeneticSequence(object):

    def __init__(self, seq, seqID = None):
        
        if seqID is None:
            seqID = 'r' + str(random.randint(1, 1000000))
            
        self.name = seqID
        self.sequence = seq
        
    def __getitem__(self, index):
        
        if isinstance(index, int):
            item = self.sequence[index]
        elif isinstance(index, slice):
            seq = self.sequence[index]
            item = type(self)(seq, seqID = self.name)

        return item
        
    def __len__(self):
        return len(self.sequence)
        
    def __str__(self):
        return self.sequence

######

class AminoAcidSequence(GeneticSequence):
    pass

######

class FastaAminoAcidSequence(AminoAcidSequence):
    
    def __str__(self):
        
        n = 60
        IDstring = '>' + self.name + '\n'
        sequence = '\n'.join([self.sequence[i: i + n] for i in range(0, len(self.sequence), n)])
        
        return IDstring + sequence

######

class DNAsequence(GeneticSequence):
    
    def complement(self):
        cseq = complement(self.sequence)

        return type(self)(cseq, self.name)
        
    def reverse(self):
        rsequence = reverse(self.sequence)

        return type(self)(rsequence, self.name)

    def reverse_complement(self):
        cseq = complement(self.sequence)
        rcseq = reverse(cseq)

        return type(self)(rcseq, self.name)

    def _translate(self, frame = 0):

        codons = (codon for codon in non_overlapping_iter(self.sequence, size = 3, start = frame))
        aas = ''
        for codon in codons:
            if '-' in codon or 'N' in codon:
                aas += 'X'
            elif len(codon) != 3:
                aas += 'X'
            else:
                aas += dna_to_amino[codon]

        return aas

    def translate(self, frame = 0):
        
        aas = self._translate(frame)

        return AminoAcidSequence(aas, self.name)
        

###### ######

class FastaSequence(DNAsequence):
    
    def __str__(self):
        
        n = 60
        IDstring = '>' + self.name + '\n'
        sequence = '\n'.join([self.sequence[i: i + n] for i in range(0, len(self.sequence), n)])
        
        return IDstring + sequence

    def translate(self, frame = 0):
        
        aas = self._translate(frame)

        return FastaAminoAcidSequence(aas, self.name)



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
    
    
    def __init__(self, iter_benotypes, possible_alleles):
        self.benotypes = [b for b in iter_benotypes]
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
        
    def is_heterozygote(self):
        ishet = []
        for b in self.iter_benotypes():
            if b != ('.', '.'):
                if b[0] != b[1]:
                    ishet.append(True)
                else:
                    ishet.append(False)
        
        return ishet
        
        
    def genotype_numbers(self, binary_genotypes, ualleles):
        """
        Assumes mono and bi-allelic sites only.
        
        binary_genotypes = ['00', '01', ...'10']
        ualleles = ['0', '1']
        
        returns [num_AA, num_Aa, num_aa]
        """
        derived = ualleles[1]
        ngens = [0 for i in range(3)]
        for g in binary_genotypes:
            c = g.count(derived)
            ngens[c] += 1
            
        return ngens

    
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


