class Site(object):
    
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


class FastaSite(Site):
    pass


###### ######


class VCFSite(Site):

    def alleles(self):
        raise IndexError, "self.alleles() not implemented yet"
