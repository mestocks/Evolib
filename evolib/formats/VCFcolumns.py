
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

        * There are many vcf files that may have an INFO field 
        that is present on some rows but not others. There are 
        two ways to deal with this:
          1) allow __getitem__ to throw an IndexError when a field
             is absent. The problem with this is that there is not 
             a satisfying, efficient way to either check for this 
             absence, or to safely Except the IndexError. Things get
             messy.
          2) make __getitem__ return None when a field is absent.
        Probably best to go with option 2.
        """
        if isinstance(self.value, str):
            self.value = dict((tuple(i.split('=')) for i in self.chr_value.split(';') if '=' in i))

        # Return None if field is not present in INFO
        if index not in self.value:
            item = None
        else:
            item = self.value[index]
        
        return item
        
        
###### ######
    
class FORMAT(COL_BASECLASS):

    def __init__(self):
        self.format_dict = {}

    def __getitem__(self, value):
        
        if value not in self.format_dict:
            new_value = dict(((j, i) for (i, j) in enumerate(value.split(":"))))
            self.format_dict.update({value: new_value})

        return self.format_dict[value]

###### ######

class SAMPLE(COL_BASECLASS):
    __slots__ = ['chr_value', 'value', 'Format', 'name', 'nvalues']
    
    def __init__(self, value, Format, name = None):
        self.chr_value = value
        self.value = value
        self.Format = Format
        self.name = name
        
    def __getitem__(self, key):
        
        if self.chr_value == './.':
            item = None
        else:
            if isinstance(self.value, str):
                self.value = self.chr_value.split(':')
                self.nvalues = len(self.value)

            fdict = self.Format[self.Format.value]
            if key in fdict:
                findex = fdict[key]
                if findex >= self.nvalues:
                    item = None
                else:
                    item = self.value[findex]
            else:
                item = None
            
        return item

    def __str__(self):
        return self.chr_value

