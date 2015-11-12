import re
import sys

#from libc.string cimport const_char
#typedef const char specialChar

from cython cimport array

def split_iter2(rowstring):#, ncols, delim = '\t'):
    #str crowstring = rowstring
    for i in rowstring:
        print i

def split_iter(string, delim = '\t'):
    for col in re.finditer(r'[^' + delim + ']+', string.rstrip()):
        yield col.group(0)
        
def type_iter(string, rowtypes):
    i = 0
    for col in split_iter(string):
        yield rowtypes[i](col)
        i+=1

class INFO(object):
    def __init__(self, string):
        self.items = ((i.split('=')[0], i.split('=')[1]) for i in split_iter(string, ';'))
        self.itemdict = None
    def __getitem__(self, value):
        if self.itemdict is None:
            self.itemdict = dict(self.items)
        return self.itemdict[value]

class FORMAT(object):
    def __init__(self, string):
        self.items = split_iter(string, ':')
        self.itemlist = None
    def __getitem__(self, value):
        if self.itemlist is None:
            self.itemlist = list(self.items)
        return self.itemlist[value]
    
class SAMPLE(object):
    def __init__(self, string):
        self.items = split_iter(string, ':')
        self.itemlist = None
    def __getitem__(self, value):
        if self.itemlist is None:
            self.itemlist = list(self.items)
        return self.itemlist[value]
    
class VCFrow(object):
    def __init__(self, string, rowtypes, headerDict):
        self.row = type_iter(string, rowtypes)
        self.headerDict = headerDict
        self.rowlist = None

    def __getitem__(self, value):
        if self.rowlist is None:
                self.rowlist = list(self.row)
        if isinstance(value, str):
            value = self.headerDict[value]
        return self.rowlist[value]
            
class VariantCallFormat(object):
    def __init__(self, FileObject):
        self.FileObject = FileObject

    def __iter__(self):
        for row in self.FileObject:
            if row.startswith("##"):
                pass
            elif row.startswith("#CHROM"):
                headers = row[1:].rstrip().split('\t')
                headerDict = dict(((j, i) for i, j in enumerate(headers)))
                rowtypes = [str, int, str,
                            str, str, float,
                            str, str, str] + [str] * (len(headers) - 9)
                            #str, INFO, FORMAT] + [SAMPLE] * (len(headers) - 9)
            else:
                yield VCFrow(row, rowtypes, headerDict)
    
if __name__ == "__main__":

    for line in sys.stdin:
        split_iter2(line)
    """
    myVCF = VariantCallFormat(sys.stdin)

    for row in myVCF:
        chrom = row['CHROM']
        pos = row['POS']
        ID = row['ID']
        ref = row['REF']
        alt = row['ALT']
        info = row['INFO']
        frm = row['FORMAT']
        s1 = row[12]
        #info1 = info['AC']
        #frm1 = frm[0]
        
        #print chrom, pos, ID, ref, s1[0]#, frm1, info1
    """
