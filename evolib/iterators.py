from evolib.formats.IteratorObjects import FastqRead, GFFRecord, msFormat, VCFrow, FastaAlignment

from evolib.formats.VCFcolumns import INFO, FORMAT, SAMPLE

import os
import collections


###### ######


def fasta_iter(file_paths):

    if isinstance(file_paths, str):
        file_paths = [file_paths]

    assert isinstance(file_paths, collections.Iterable), "file_paths is not iterable: %r" % file_paths
    
    for fpath in file_paths:

        assert os.path.isfile(fpath), "fpath does not exist: %r" % fpath
        
        FileObject = open(fpath, 'r')
        FsaFormat = FastaAlignment(FileObject)
        
        yield FsaFormat
        

###### ######


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


###### ######


def gff_iter(FileObject):
    
    name = 'start'
    for line in FileObject:
        if line.startswith('#') is False:
            items = line.rstrip().split('\t')
            if name == 'start':
                name = items[0]
                Gff = GFFRecord()
                Gff.add(items)
                
            elif items[0] != name:
                yield Gff
                name = items[0]
                Gff = GFFRecord()
                Gff.add(items)
                
            else:
                Gff.add(items)
            
    yield Gff



###### ######

class Header(object):
    
    def __init__(self, header):
        self.names = header
        self.nsamples = len(header) - 9
        self.classes = [str, int, str, str, str, float, str, INFO, None] + [SAMPLE] * self.nsamples
        self.str_item = dict([(key, (self.classes[index], index)) for index, key in enumerate(header)])
        self.int_item = dict([(index, (self.classes[index], index)) for index, key in enumerate(header)])
        self.preamble = ''

    def __str__(self):
        return self.preamble

###### ######

def vcf_iter(FileObject):
    
    preamble = ''
    for line in FileObject:
        
        if line[0] != '#':
            # rstrip costs 0.5s per 1 million rows
            values = line.rstrip().split('\t')
            Row.values = values
            Row.Format.value = values[8]
            
            yield Row
            
        elif line[:2] == '##':
            preamble += line
        
        else:
            preamble += line.rstrip()
            header = line[1:].rstrip().split('\t')
            headerClass = Header(header)
            Format = FORMAT()
            headerClass.preamble = preamble
            Row = VCFrow([None] * len(header), headerClass, Format)
            
###### ######


def ms_iter(fileObject):
    
    header = True
    iteration = ''
    
    for line in fileObject:
        
        if header is False:
            if line.startswith('//'):
                
                msClass = msFormat(iteration.rstrip())
                yield msClass
                iteration = ''
                
            else:
                iteration += line
                
        elif header is True:
            if line.startswith('//'):
                header = False
                
    msClass = msFormat(iteration.rstrip())
    yield msClass


###### ######


def fwdpp_iter(fileObject):

    header = True
    iteration = ''
    ns = 0
    
    for line in fileObject:
        
        if header is False:
            if line.startswith('//'):
                if ns % 2 == 0:
                    selClass = msFormat(iteration.rstrip())
                    yield neuClass, selClass
                else:
                    neuClass = msFormat(iteration.rstrip())
                iteration = ''
                ns += 1
            else:
                iteration += line
                
        elif header is True:
            if line.startswith('//'):
                ns += 1
                header = False

    if ns % 2 == 0:
        selClass = msFormat(iteration.rstrip())
        yield neuClass, selClass
