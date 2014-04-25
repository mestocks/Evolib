from lib.DataObjects import SequenceData
from lib.FileIterators import fastq_iter, vcf_iter

###### ######

class VariantCallFormat(SequenceData):
    """
    
    Usage:
    
    VariantCallFormat(filename)
    
    Returns a VCF class. Each row is represented by a <ROWCLASS>, 
    with each column represented by a <COLCLASS>.
    
    Example:
    
        >>> from SequenceFormats import VariantCallFormat
        >>> myVCF = VariantCallFormat("data.vcf")
        >>> for bp in myVCF:
        >>>     print bp['CHROM'], bp['POS'], bp['REF'], bp.genotypes(), bp.heterozygosity()
    """
    
    def __init__(self, FileObject):
        
        self.FileObject = FileObject
        
        #self.chromosome_info = dict(self.chr_info())
        
        
    def __iter__(self):
        for row in self.bySite():
            yield row
            
    
    def bySite(self):
        
        for row in vcf_iter(self.FileObject):
            yield row
                
    def chr_info(self):
        ### This no longer works (as it is dependent on 
        ### a file name rather than a file object. This 
        ### has downstream consequences for interating 
        ### through multiple vcf files.
        chromosomes = []
        past_chrom = None
        reverse = None
        
        file = open(self.file_name, 'r')
        
        for row in file:
            
            if row.startswith('#') is False:
                
                value_list = row.rstrip().split('\t')
                chrom = value_list[0]
                
                if chrom != past_chrom:
                    past_chrom = chrom
                    past_pos = int(value_list[1])
                    reverse = None
                    add_info = (chrom, )
                    
                elif reverse is None:
                    
                    if past_pos < int(value_list[1]):
                        reverse = False
                        add_info += (reverse, )
                        chromosomes.append(add_info)
                        
                    elif past_pos > int(value_list[1]):
                        reverse = True
                        add_info += (reverse, )
                        chromosomes.append(add_info)
                    
        return chromosomes

###### ######

class FastqFormat():
    
    def __init__(self, FileObject):
        self.FileObject = FileObject
        
    def __iter__(self):
        for row in self.bySite():
            yield row
            
    def bySite(self):
        for row in fastq_iter(self.FileObject):
            yield row
