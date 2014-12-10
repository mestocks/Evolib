from evolib.formats.IteratorObjects import FastqRead, GFFRecord, msFormat, VCFrow, FastaAlignment, VCFrow3

from evolib.formats.VCFcolumns import CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, FORMAT2, SAMPLE, SAMPLE2

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


def vcf_iter(FileObject):
    
    for line in FileObject:
        
        if line.startswith('##'):
            pass
        
        elif line.startswith('#'):
            header = line[1:].rstrip().split('\t')
            nsamples = len(header) - 9
            col_classes = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT] + [SAMPLE] * nsamples
            
        else:
            value_list = line.rstrip().split('\t')
            Row = VCFrow(value_list, col_classes, header, nsamples)
            
            yield Row

###### ######

class Header(object):
    
    def __init__(self, header):
        self.names = header
        self.nsamples = len(header) - 9
        self.classes = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, None] + [SAMPLE2] * self.nsamples
        self.str_item = dict([(key, (self.classes[index], index)) for index, key in enumerate(header)])
        self.int_item = dict([(index, (self.classes[index], index)) for index, key in enumerate(header)])

def vcf_iter3(FileObject):
    
    for line in FileObject:
        
        if line.startswith('##'):
            pass
        
        elif line.startswith('#'):
            header = line[1:].rstrip().split('\t')
            headerClass = Header(header)
            Format = FORMAT2()
            
        else:
            values = line.rstrip().split('\t')
            Row = VCFrow3(values, headerClass, Format)
            
            yield Row
        

def vcf_iter2(FileObject):
    
    for line in FileObject:
        
        if line.startswith('##'):
            pass
        
        elif line.startswith('#'):
            header = line[1:].rstrip().split('\t')
            nsamples = len(header) - 9
            Format = FORMAT2()
            col_classes = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, Format] + [SAMPLE2] * nsamples
            
        else:
            value_list = line.rstrip().split('\t')
            Row = VCFrow(value_list, col_classes, header, nsamples)
            
            yield Row


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
