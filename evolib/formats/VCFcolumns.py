
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
    
class INFO(COL_BASECLASS):
    __slots__ = ['chr_value', 'value']
    
    def __getitem__(self, index):
        """
        The actual parsing of the INFO line only occurs when 
        __getitem__ is called. This is to save on unnecesary 
        processing of the VCF file that occurs when not all 
        columns are used. The same also applies to each item 
        within INFO.
        
        Each INFO entry is delimited by ';'. Each entry consist of the descriptor and 
        value separated by '='. E.g:
        INFO
        DP=1;AF=0;AC1=0;DP4=1,0,0,0;MQ=25;FQ=-24.3
        """
        if isinstance(self.value, str):
            self.value = dict([tuple(i.split('=')) for i in chr_value.split(';') if '=' in i])
        
        INFO_parse = {'DP': int, 'MQ': int}
        
        value = self.value[index]
        
        if index in INFO_parse:
            value = INFO_parse[index](value)
        
        return value
        
        
###### ######
    
class FORMAT(COL_BASECLASS):

    def __init__(self):
        self.format_dict = {}

    def __getitem__(self, value):
        
        if value not in self.format_dict:
            new_value = dict([(j, i) for (i, j) in enumerate(value.split(":"))])
            self.format_dict.update({value: new_value})

        return self.format_dict[value]

###### ######

class SAMPLE(COL_BASECLASS):
    __slots__ = ['chr_value', 'value', 'Format']
    
    def __init__(self, value, Format):
        self.chr_value = value
        self.value = value
        self.Format = Format
        
    def __getitem__(self, key):
        
        if self.chr_value == './.':
            item = None
        else:
            if isinstance(self.value, str):
                self.value = self.chr_value.split(':')
                # problems may occur if the number of items in the
                # FORMAT column != number items in the sample
                self.SAMPLE_parse = {'DP': int}
                
            item = self.value[self.Format[self.Format.value][key]]
            
            if key in self.SAMPLE_parse:
                print key, self.value
                item = self.SAMPLE_parse[key](item)
            
        return item

    def __str__(self):
        return self.chr_value
    
    def is_het(self):
        
        gt = self['GT']
        
        if gt is None:
            het = None
        else:
            if gt == '0/1' or gt == '1/0':
                het = True
            else: 
                het = False
            
        return het
    
    def str_fetch(self, key):        

        if self.chr_value == "./.":
            item = 'NA'
        else:
            #if self.split_value is None:
            self.split_value = self.chr_value.split(":")
                
            item = self.split_value[self.Format[self.Format.value][key]]

        return item

    def genotype_str(self, ref, alt):
        possibleAlleles = [ref] + alt.split(',')
        GT = self['GT'].split('/')
        
        one, two = possibleAlleles[int(GT[0])], possibleAlleles[int(GT[1])]
        
        return one + two
