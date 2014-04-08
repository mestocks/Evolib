from lib.DataObjects import SequenceData

from lib.VCFrow import ROW_BASECLASS, ROW_BASECLASS_OLD
from lib.VCFcolumns import CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE


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
    
    def __init__(self, file_name):
        
        self.file_name = file_name
        self.header = self.get_header(file_name)
        
        nsamples = len(self.header) - 9
        self.col_classes = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT] + [SAMPLE] * nsamples
        
        self.chromosome_info = dict(self.chr_info())
        
        
    def __iter__(self):
        for row in self.bySite():
            yield row
            
    
    def bySite(self):
        
        file = open(self.file_name, 'r')
        
        for line in file:
            
            if line.startswith('#') is False:
                
                value_list = line.rstrip().split('\t')
                v_col_row = ROW_BASECLASS(value_list, self.col_classes, self.header)
                
                yield v_col_row
                
                
    def chr_info(self):
        
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
    
    
    def get_header(self, file_name):
        
        file = open(file_name, 'r')
        
        for line in file:
            if line.startswith('#') is False:
                break
            elif line.startswith('##') is True:
                pass
            else:
                header = line[1:].rstrip().split('\t')
                
        file.close()
        
        return header
