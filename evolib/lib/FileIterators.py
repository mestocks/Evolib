from DataObjects import FastqRead

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
