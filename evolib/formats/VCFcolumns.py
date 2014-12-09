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

        INFO_parse = {'DP': self._DP}#, 
                     # 'MQ': self._MQ}

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

class FORMAT2(COL_BASECLASS):

    def __init__(self):
        self.format_dict = {}

    def __getitem__(self, value):
        
        if value not in self.format_dict:
            new_value = dict([(j, i) for (i, j) in enumerate(value.split(":"))])
            self.format_dict.update({value: new_value})

        return self.format_dict[value]



###### ######

class SAMPLE2(COL_BASECLASS):

    __slots__ = ['chr_value', 'value', 'split_value', 'Format']#, 'SAMPLE_parse']
        
    def __init__(self, value, Format):
        
        self.chr_value = value
        self.value = value
        self.split_value = None
        self.Format = Format
        #self.SAMPLE_parse = {'DP': self._DP,
        #                     'GT': self._GT, 
        #                     'GQ': self._GQ, 
        #                     'PL': self._PL}
    

    def __getitem__(self, key):
        
        #self.SAMPLE_parse = {'DP': self._DP,
        #                     'GT': self._GT, 
        #                     'GQ': self._GQ, 
        #                     'PL': self._PL}
    
        if self.chr_value == "./.":
            item = None
        else:
            self.split_value = self.chr_value.split(":")
            item = self.split_value[self.Format[self.Format.value][key]]

        return item

    def getitem__(self, key):

        if self.split_value is None:
            self.split_value = self.chr_value.split(":")
            
        if self.chr_value == "./.":
            item = None
        else:
            #split_value = self.chr_value.split(":")
            item = self.split_value[self.Format[self.Format.value][key]]
            
        if key in self.SAMPLE_parse:
            item = self.SAMPLE_parse[key](item)
            
        return item
    

    def __str__(self):
        return self.chr_value

    def str_fetch(self, key):        

        if self.chr_value == "./.":
            item = '0'
        else:
            #if self.split_value is None:
            self.split_value = self.chr_value.split(":")
                
            item = self.split_value[self.Format[self.Format.value][key]]

        return item
    
    def _DP(self, item):

        if item is None:
            item = 0

        return int(item)
        
        #new_item = 0
        #if item is not None:
        #    new_item = int(item)
            
        #return new_item
    
    
    def _GT(self, item):

        new_item = './.'
        if item is not None:
            sep = '/'
            if '|' in item:
                sep = '|'
                
            new_item = item.split(sep)
        
        return new_item
    
    def _GQ(self, item):

        new_item = './.'
        if item is not None:
            new_item = int(item)
            
        return new_item
    

    def _PL(self, item):
        
        if item is not None:
            try:
                item = int(item)
            except ValueError:
                item = map(int, item.split(','))
    
        return item
        
    def is_het(self):
        
        ishet = None

        if self.split_value is None:
            self.split_value = self.chr_value.split(":")
            
        if self.chr_value == "./.":
            item = None
        else:
            item = self.split_value[self.Format[self.Format.value]['GT']]

        if item == '0/0' or item == '1/1':
            ishet = False
        elif item == '0/1' or item == '1/0':
            ishet = True

        return ishet


###### ######

class SAMPLE(COL_BASECLASS):
    __slots__ = ['col_name', 'chr_value', 'value']
    
    """
    SAMPLE["GQ"]
    :: self._get_formatted_item() [1]
    :: :: parse [1]
    :: :: add key to dictionary (few if statements)
    
    value['GT']
    {'DL:DP': {'DL': 0, 'DP': 1}}
    
    
    """
    
    def __init__(self, value, format_value):
        self.chr_value = value
        self.value = value
        self.format_value = format_value
        
    
    def __getitem__(self, key):
        
        if isinstance(key, str) is False:
            raise TypeError, 'Only str indexing is supported.'

        #self.value = self._parse(self.chr_value, self.format_value)
        #return self._GT(self.value[key])
        
        return self._get_formatted_item(key)

    def __str__(self):
        return self.chr_value
    
    
    def _get_formatted_item(self, key):
        """
        Formats certain items that need to be converted. 
        For example, 'DP' is converted to an integer. 
        However, it only does this if it has not been done
        before. This minimises repetition and results in 
        considerable increases in speed as only those 
        items that are explicitly called are parsed.
        """
        # Check if the column has been parsed.
        if isinstance(self.value, str):
            self.value = self._parse(self.chr_value, self.format_value)
            
        # Some values need to be converted. If it's 
        # not been done already then re-format by 
        # calling the relevant function.
        SAMPLE_parse = {'DP': self._DP, 
                        'GT': self._GT, 
                        'GQ': self._GQ, 
                        'PL': self._PL}
        
        if key in SAMPLE_parse.keys():
            item = self.value[key]
            if item is not None:
                if isinstance(item, str):
                    self.value[key] = SAMPLE_parse[key](item)
            elif item is None:
                # Missing DP values should be 0 rather than None.
                if key == 'DP':
                    self.value[key] = 0
                
        return self.value[key]
    
        
    def _parse(self, chr_value, format_col):
        """
        Splits and parses the string by ':' to give a 
        dictionary. For example, '0/1:23:0,45,0', returns 
        {'GT': '0/1', 'DP': '23', 'PL': '0,45,0'}. Some vcf files do not 
        contain all format fields (when there are no 
        reads mapping to the position for a certain 
        individual). In such scenarios, None is entered 
        into the relevant position in the list. For example, 
        './.' would parse to {'GT': './.', 'DP': None, 'PL': None}. 
        """
        value = []
        list_value = chr_value.split(":")
        
        for i in range(len(format_col)):
            f = format_col[i]
            
            try:
                v = list_value[i]
            except IndexError:
                v = None
            
            value.append((f, v))
        
        return dict(value)
    
    
    def _DP(self, item):
        return int(item)
    
    
    def _GT(self, item):
        """
        
        """
        sep = '/'
        if '|' in item:
            sep = '|'
            
        #new_item = []
        #item_split = item.split(sep)
        #for i in item_split:
        #    new_item.append(i)
        new_item = item.split(sep)
        
        return new_item
    
    def _GQ(self, item):        
        return int(item)
    

    def _PL(self, item):

        try:
            item = int(item)
        except ValueError:
            item = map(int, item.split(','))
    
        return item
    
    def binary_call(self):
        # This will throw an error if self['GT'] 
        # is None (i.e. non-iterable). Do I need to 
        # take this into account?
        return tuple(self['GT'])
    
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
        gCall = (possibleAlleles[int(GT[0])], possibleAlleles[int(GT[1])])
        
        return gCall

    def is_heterozygote(self, ref, alt):
        
        genotype = self.genotype_call(ref, alt)
        
        if genotype[0] == genotype[1]:
            ishet = False
        else:
            ishet = True
            
        return ishet

