class COL_BASE_STR(str):
    
    def __init__(self, value):
        self.chr_value = value
        
        return str.__init__(self)


class COL_BASE_INT(int):
    
    def __init__(self, value):
        self.chr_value = value
        
        return int.__init__(self)

class COL_BASECLASS(object):
    def __init__(self, value):
        """
        The original string version (as it appeared in the .vcf file) 
        of 'value' is preserved, with type 'string', as 'self.chr_value'. 
        The "working" version of 'value' is stored as 'self.value'. This allows 
        the type to change depending on the data type whilst maintaining the 
        original formatting.
        """
        self.chr_value, self.value = value, value
    
    def __str__(self):
        return self.chr_value

#######################

class CHROM(COL_BASE_STR):
    #__slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'CHROM'

class POS(COL_BASE_INT):
    __slots__ = ['col_name', 'chr_value']
    col_name = 'POS'
    
    def __init__(self, value):
        self.chr_value, self.value = value, int(value)

class ID(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'ID'

class REF(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'REF'

class ALT(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'ALT'

class QUAL(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'QUAL'
    
    def __init__(self, value):
        self.chr_value, self.value = value, float(value)

class FILTER(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'FILTER'

class INFO(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'INFO'
    
    def __init__(self, value):
        self.chr_value = value
        self.value = self._parse(value)
    
    def __getitem__(self, index):
        
        return self.value[index]
    
    def _parse(self, chr_value):
        """ 
        Each INFO entry is delimited by ';'. Each entry consist of the descriptor and 
        value separated by '='. E.g:
        INFO
        DP=1;AF=0;AC1=0;DP4=1,0,0,0;MQ=25;FQ=-24.3
        """
        value = dict([tuple(i.split('=')) for i in chr_value.split(';') if '=' in i])
        value = self._DP(value)
        value = self._EFF(value)
        value = self._MQ(value)
        
        return value
    
    def _DP(self, value):
        if 'DP' in value.keys():
            try:
                value['DP'] = int(value['DP'])
            except TypeError:
                pass
        
        return value
    
    def _EFF(self, value):
        
        if 'EFF' in value.keys():
            value['EFF'] = EFF(value['EFF'])
        
        return value
    
    def _MQ(self, value):
        
        if 'MQ' in value.keys():
            value['MQ'] = int(value['MQ'])
        
        return value

class FORMAT(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    col_name = 'FORMAT'
    
    def __init__(self, value):
        self.chr_value = value
        self.value = self._parse(value)
    
    def _parse(self, chr_value):
        
        value = chr_value.split(':')
        
        return value

class SAMPLE(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    
    def __init__(self, value, format):
        self.chr_value = value
        self.value = self._parse(value, format)
    
    def __getitem__(self, key):
        return self.value[key]
        
    def _parse(self, chr_value, format):
        
        list_value = chr_value.split(':')
        value = dict(zip(format, list_value))
        value = self._GQ(value)
        value = self._GT(value)
        value = self._PL(value)
        value = self._DP(value)
        
        return value
    
    def _DP(self, value):
        
        if 'DP' in value.keys():
            value['DP'] = int(value['DP'])
        
        return value
    
    def _GQ(self, value):
        
        if 'GQ' in value.keys():
            value['GQ'] = int(value['GQ'])
        
        return value
    
    def _GT(self, value):
        
        sep = '/'
        
        if 'GT' in value.keys():
            if '|' in value['GT']:
                sep = '|'
            value['GT'] = map(int, value['GT'].split(sep))
        
        return value
    
    def _PL(self, value):
        
        if 'PL' in value.keys():
            try:
                value['PL'] = int(value['PL'])
            except ValueError:
                value['PL'] = map(int, value['PL'].split(','))
    
        return value
    
    
    
    def genotype_call(self, ref, alt):
        
        possibleAlleles = [ref] + alt.split(',')
        #alleleCombinations = list(itertools.combinations_with_replacement(possibleAlleles, 2))
        #binaryCombinations = itertools.combinations_with_replacement(range(len(possibleAlleles)), 2)
        #indexOrder = [((k * (k + 1)) / 2) + j for (j, k) in list(binaryCombinations)]
        #genotypeCombinations = [alleleCombinations[i] for i in indexOrder]
        #print genotypeCombinations
        #PL = self.value['PL']
        #findZero = PL.index(0)
        #gCall = genotypeCombinations[findZero]
        #print self.value
        gCall = (possibleAlleles[self.value['GT'][0]], possibleAlleles[self.value['GT'][1]])
        
        return gCall

    def is_heterozygote(self, ref, alt):
        
        genotype = self.genotype_call(ref, alt)
        
        if genotype[0] == genotype[1]:
            ishet = False
        else:
            ishet = True
            
        return ishet

