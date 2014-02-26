class ROW_BASECLASS(list):
    
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

    def genotypes(self):
        
        REF = self['REF']
        ALT = self['ALT']
        
        ref = REF.value
        alt = ALT.value

        if alt == '.':
            genotypes = [(ref, ref) for i in self.__getitem__(slice(9, None))]
        else:
            genotypes = []
            for sample in self.__getitem__(slice(9, None)):
                genotype = sample.genotype_call(ref, alt)
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
        
