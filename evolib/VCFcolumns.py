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
    """
    
    """
    
    def __init__(self, value):
        self.chr_value = value
        self.value = value

    def __getitem__(self, index):
        """
        The actual parsing of the INFO line only occurs when 
        __getitem__ is called. This is to save on unnecesary 
        processing of the VCF file that occurs when not all 
        columns are used. The same also applies to each item 
        within INFO.
        """
        if isinstance(self.value, str):
            self.value = self._parse(self.chr_value)

        INFO_parse = {'DP': self._DP, 
                      'MQ': self._MQ}

        if index in INFO_parse.keys():
            self.value = INFO_parse[index](self.value)

        value = self.value[index]

        return value
    
    
    def _parse(self, chr_value):
        """ 
        Each INFO entry is delimited by ';'. Each entry consist of the descriptor and 
        value separated by '='. E.g:
        INFO
        DP=1;AF=0;AC1=0;DP4=1,0,0,0;MQ=25;FQ=-24.3
        """
        value = dict([tuple(i.split('=')) for i in chr_value.split(';') if '=' in i])
        
        return value
    
    def _DP(self, value):
        if 'DP' in value.keys():
            try:
                value['DP'] = int(value['DP'])
            except TypeError:
                pass
        
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
    
    def __init__(self, value, format_value):
        self.chr_value = value
        self.value = value
        self.format_value = format_value
    
    def __getitem__(self, key):

        if isinstance(self.value, str):
            self.value = self._parse(self.chr_value, self.format_value)
            
        SAMPLE_parse = {'DP': self._DP, 
                        'GT': self._GT, 
                        'GQ': self._GQ, 
                        'PL': self._PL}
        
        if key in SAMPLE_parse.keys():
            item = self.value[key]
            if isinstance(item, str):
                self.value[key] = SAMPLE_parse[key](item)

        return self.value[key]
        
    def _parse(self, chr_value, format):
        
        list_value = chr_value.split(':')
        value = dict(zip(format, list_value))
        
        return value
    
    def _DP(self, item):
        return int(item)
    
    def _GT(self, item):
        
        sep = '/'
        if '|' in item:
            sep = '|'
        item = map(int, item.split(sep))
        
        return item
    
    def _GQ(self, item):        
        return int(item)
    

    def _PL(self, item):

        try:
            item = int(item)
        except ValueError:
            item = map(int, item.split(','))
    
        return item
    
    
    
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
        GT = self['GT']
        gCall = (possibleAlleles[GT[0]], possibleAlleles[GT[1]])
        
        return gCall

    def is_heterozygote(self, ref, alt):
        
        genotype = self.genotype_call(ref, alt)
        
        if genotype[0] == genotype[1]:
            ishet = False
        else:
            ishet = True
            
        return ishet

