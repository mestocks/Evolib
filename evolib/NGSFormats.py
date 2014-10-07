from evolib.lib.DataObjects import SequenceData
from format.Iterators import fastq_iter, vcf_iter

###### ######

class VariantCallFormat(SequenceData):
    """
    
    Usage:
    
    VariantCallFormat(filename)
    
    Returns a VCF class. Each row is represented by a <ROWCLASS>, 
    with each column represented by a <COLCLASS>.
    
    Example:
    
        >>> from SequenceFormats import VariantCallFormat
        >>> openFile = open("myData.vcf", 'r')
        >>> myVCF = VariantCallFormat(openFile)
        >>> for bp in myVCF:
        >>>     print bp['CHROM'], bp['POS'], bp['REF'], bp.genotypes(), bp.heterozygosity()
    """
    
    def __init__(self, FileObject):
        
        self.FileObject = FileObject
        
        
    def __iter__(self):
        for row in self.bySite():
            yield row
            
    
    def bySite(self):
        
        for row in vcf_iter(self.FileObject):
            yield row

###### ######

class FastqFormat():
    
    def __init__(self, FileObject):
        self.FileObject = FileObject
        
    def __iter__(self):
        for row in fastq_iter(self.FileObject):
            yield row
            
    def __len__(self):
        return sum(1 for row in fastq_iter(self.FileObject))
