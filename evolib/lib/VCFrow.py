import numpy
from scipy import stats

from DataObjects import Site
from DNAobjects import Genotypes
from StatMethods import chisquared

class ROW_BASECLASS(list, Site):
    
    def __init__(self, values, classes, header):
         
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
            
"""
    def benotypes(self):
        
        btypes = []
        for sample in self.__getitem__(slice(9, None)):
            btypes.append(sample.binary_call())
        
        return btypes
    
    def get_genotypes(self):
        possible_alleles = self.possible_alleles()
        benotypes = self.benotypes()
        
        return Genotypes(benotypes, possible_alleles)

    def genotypes(self):
        
        REF = self['REF']
        ALT = self['ALT']
        
        ref = REF.value
        alt = ALT.value

        #if alt == '.':
        #    genotypes = [(ref, ref) for i in self.__getitem__(slice(9, None))]
        #else:
        genotypes = []
        for sample in self.__getitem__(slice(9, None)):
            genotype = sample.genotype_call(ref, alt)
            genotypes.append(genotype)
        
        return genotypes
    
    def iter_samples(self):
        for Sample in self.__getitem__(slice(9, None)):
            yield Sample
            
    
    def number_of_alleles(self, include = None):
        
        btypes = self.benotypes()
        alls = []
        index = 0
        for gen in btypes:
            if include is None or index in include:
                if gen[0] != None and gen[1] != None:
                    bps = gen[0] + gen[1]
                    if bps != 'NN':
                        alls.append(bps)
            index += 1
        allstr = ''.join(alls)
        
        return len(set(allstr))
    
    def number_of_genotypes(self, include = None):
        
        btypes = self.benotypes()
        gens = []
        index = 0
        for b in btypes:
            if include is None or index in include:
                if b[0] != None and b[1] != None:
                    gens.append(b[0] + b[1])
            index += 1
        
        return len(set(gens))
    
    def possible_alleles(self):
        
        ref, alt = self['REF'].value, self['ALT'].value
        print ref, alt
        alleles = [ref] + alt.split(',')
        
        return alleles
    
    def variationIsSegregating(self):
        
        benotypes = []
        for sample in self.__getitem__(slice(9, None)):
            benotype = sample['GT']
            if benotype != [None, None]:
                benotypes.append(benotype)
        benotypes = [b[0] + b[1] for b in benotypes]
            
        if len(set(benotypes)) > 1:
            isSegregating = True
        else:
            isSegregating = False
        
        return isSegregating
    
    def alleles(self):
        genotypes = self.genotypes()
        gametes = ''
        for gamete in genotypes:
            if gamete == (None, None):
                gametes += 'NN'
            else:
                gametes += gamete[0] + gamete[1]
        #gametes = ''.join([gamete[0] + gamete[1] for gamete in genotypes])
        #site = Site(gametes)
        
        return gametes

    
    def heterozygosity(self, minGQ = 0, minDP = 0):
        """
"""
        Calculate the hetorozygosity for this site. Counts the number of individuals that are
        heterozygotes and then divides this by the number of individuals. 
"""
        """
        REF = self['REF']
        ALT = self['ALT']
        
        ref = REF.value
        alt = ALT.value
        
        if alt == '.':
            htzgsty = 0.0
        else:
            ishet = [sample.is_heterozygote(ref, alt) for sample in self.__getitem__(slice(9, None)) if sample['GQ'] >= minGQ and sample['DP'] >= minDP]
            nsamples = len(ishet)
            nhet = sum(ishet)
            
            if nhet == 0:
                htzgsty = 0.0
            else:
                htzgsty = nhet / float(nsamples)
            
        return htzgsty
    
    def chi_squared(self, case_indices, control_indices):
        
        btypes = self.benotypes()
        
        case_genotypes = [btypes[i][0] + btypes[i][1] for i in case_indices]
        control_genotypes = [btypes[i][0] + btypes[i][1] for i in control_indices]
        
        ualleles = list(set(''.join(case_genotypes + control_genotypes)))
        ualleles.sort()
        
        gens = [ualleles[0] + ualleles[0], ualleles[0] + ualleles[1], ualleles[1] + ualleles[1]]
        
        case_count = [case_genotypes.count(gen) for gen in gens]
        control_count = [control_genotypes.count(gen) for gen in gens]
        total_count = [case_count[0] + control_count[0], 
                       case_count[1] + control_count[1],
                       case_count[2] + control_count[2]]
        
        
        R, S = sum(case_count), sum(control_count)
        N = sum(total_count)
    
        obsCase_A = (2 * case_count[0]) + case_count[1]
        obsCase_a = case_count[1] + (2 * case_count[2])
        obsControl_A = (2 * control_count[0]) + control_count[1]
        obsControl_a = control_count[1] + (2 * control_count[2])
        obsR, obsS = 2 * R, 2 * S
        obsN_A = obsCase_A + obsControl_A
        obsN_a = obsCase_a + obsControl_a
        obsN = obsN_A + obsN_a
    
        expCase_A = ( (2 * R) * ( (2 * total_count[0]) + total_count[1] ) ) / (2.0 * N)
        expCase_a = ((2 * R) * (total_count[1] + (2 * total_count[2]))) / (2.0 * N)
        expControl_A = ((2 * S) * ((2 * total_count[0]) + total_count[1])) / (2.0 * N)
        expControl_a = ((2 * S) * (total_count[1] + (2 * total_count[2]))) / (2.0 * N)
        
        obs = numpy.array([obsCase_A, obsCase_a, obsControl_A, obsControl_a])
        exp = numpy.array([expCase_A, expCase_a, expControl_A, expControl_a])
        chi = stats.chisquare(obs, exp, 1)
        p_value = chi[1]
        """
"""
        chis = []
        for i in range(len(obs)):
            x2 = chisquared(obs[i], exp[i])
            if x2 is not None:
                chis.append(x2)
        
        chi = sum(chis)
        
        p_value = None
        
        return chi
"""
        """
        return p_value
"""
####################

class ROW_BASECLASS_OLD(list):
    
    def __str__(self):
        return '\t'.join(map(str, self))
    
    def __getitem__(self, index):
        
        if isinstance(index, str):
            newitem = [i for i in list.__getitem__(self, slice(None, None)) if i.col_name == index][0]
        else:
            newitem = list.__getitem__(self, index)
        
        return newitem
    
    def genotypes(self):
        
        if self['ALT'].value == '.':
            genotypes = [(self['REF'].value, self['REF'].value) for i in self.__getitem__(slice(9, None))]
        else:
            genotypes = []
            for sample in self.__getitem__(slice(9, None)):
                genotype = sample.genotype_call(self['REF'].value, self['ALT'].value)
                genotypes.append(genotype)
        
        return genotypes
    
    def variationIsSegregating(self):
        genotypes = self.genotypes()
        gametes = ''.join([gamete[0] + gamete[1] for gamete in genotypes])
        if len(set(gametes)) > 1:
            isSegregating = True
        else:
            isSegregating = False
        
        return isSegregating
    
    def alleles(self):
        genotypes = self.genotypes()
        gametes = ''.join([gamete[0] + gamete[1] for gamete in genotypes])
        
        return gametes

    
    def heterozygosity(self, minGQ = 0, minDP = 0):
        """
        Calculate the hetorozygosity for this site. Counts the number of individuals that are
        heterozygotes and then divides this by the number of individuals. 
        """
        if self['ALT'].value == '.':
            htzgsty = 0.0
        else:
            ishet = [sample.is_heterozygote(self['REF'].value, self['ALT'].value) for sample in self.__getitem__(slice(9, None)) if sample['GQ'] >= minGQ and sample['DP'] >= minDP]
            nsamples = len(ishet)
            nhet = sum(ishet)
            
            if nhet == 0:
                htzgsty = 0.0
            else:
                htzgsty = nhet / float(nsamples)
            
        return htzgsty
        
