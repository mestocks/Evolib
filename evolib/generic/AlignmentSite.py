class Site(object):
    
    def __init__(self, alleles):
        self._alleles = alleles
        
    
    def __str__(self):
        return self.alleles()

    def alleles(self):
        return self._alleles
    
    
    def has_missing_data(self, dna = ['A', 'T', 'C', 'G']):
        
        if set(self.alleles()) <= set(dna):
            answer = False
        else:
            answer = True
            
        return answer
    
    def number_of_alleles(self):
        alleles = set(self.alleles())
        return len(alleles)


###### ######


class FastaSite(Site):
    pass


###### ######


class VCFSite(Site):

    def alleles(self):
        if str(self['ALT']) == ".":
            alleles = ''.join([str(self['REF']) + str(self['REF']) for i in self.iter_samples()])
        else:
            alleles = ''.join([s.genotype_str(str(self['REF']), str(self['ALT'])) for s in self.iter_samples()])

        return alleles
