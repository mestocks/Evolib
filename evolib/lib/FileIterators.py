from DataObjects import FastqRead
from VCFrow import ROW_BASECLASS
from VCFcolumns import CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE

def fastq_iter(FileObject):
    
    linenum = 0
    for line in FileObject:
        
        cline = line.rstrip()
        remainder = linenum % 4
        
        if remainder == 0:
            assert cline.startswith('@')
            seqid = cline
        elif remainder == 1:
            seq = cline
        elif remainder == 2:
            assert cline.startswith('+')
        elif remainder == 3:
            quality = cline
            FsqR = FastqRead(seq, seqid, quality)
            yield FsqR
        else:
            raise IndexError, "Problem with counter. remainder should not be > 3."
        
        linenum += 1

def vcf_iter(FileObject):
    
    for line in FileObject:
        
        if line.startswith('##'):
            pass
        
        elif line.startswith('#'):
            header = line[1:].rstrip().split('\t')
            
        else:
            value_list = line.rstrip().split('\t')
            nsamples = len(header) - 9
            col_classes = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT] + [SAMPLE] * nsamples
            Row = ROW_BASECLASS(value_list, col_classes, header)
            
            yield Row
