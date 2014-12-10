###### ######

class FastqRead():
    def __init__(self, seq, seqid, qual):
        self.seq = seq
        self.seqid = seqid
        self.qual = qual
        
    def __str__(self):
        return '\n'.join([self.seqid, self.seq, '+', self.qual])
    
    def quality(self):
        qualities = map(ord, self.qual)
        return qualities

    def __len__(self):
        return len(self.seq)


###### ######

from evolib.generic.AlignmentSite import VCFSite

# 


class VCFrow3(object):
    # row["CHROM"]
    # row[:2] -> refers only to samples

    # header.names = ['CHROM']
    # header.str_item = {'CHROM': (CHROM, 0), ...}
    # header.int_item = {0: (CHROM, 0), ...}

    def __init__(self, values, header, Format):
        self.values = values
        self.header = header
        self.Format = Format
        self.Format.value = self.values[8]
        self.nsamples = header.nsamples

    def __str__(self):
        return '\t'.join(self.values)

    def __getitem__(self, item):

        #if isinstance(item, slice):
        #    value = self.values[slice]
        #else:
        if isinstance(item, str):
            itemClass, index = self.header.str_item[item]
        elif isinstance(item, int):
            itemClass, index = self.header.int_item[item]
            
        if index > 8:
            value = itemClass(self.values[index], self.Format)
        else:
            value = itemClass(self.values[index])

        return value

    def smap(self, mode, key):
        # currently ~ 8,500 rows/second (inc. CHROM and POS)
        return map(mode, (col.str_fetch(key) for col in self.iter_samples()))

    def iter_samples(self):
        
        maxRange = self.nsamples + 9
        SampleClass, index = self.header.int_item[9]

        for i in xrange(9, maxRange):
            yield SampleClass(self.values[i], self.Format)


class VCFrow(VCFSite):

    #__slots__ = ['NA', 'classes', 'header', 'values', 'nsamples', 'FormatClass', 'lookupindex']
    
    def __init__(self, values, classes, header, nsamples):
        
        self.NA = '.'
        self.classes = classes
        self.header = header
        self.values = values
        self.nsamples = nsamples
        
        self.FormatClass = None
        self.lookupindex = None
        
    
    def __getitem__(self, item):

        if self.lookupindex is None:
            self.lookupindex = dict([(key, index) for index, key in enumerate(self.header)])
            
        if isinstance(item, str):
            item = self.lookupindex[item]

        if isinstance(item, slice):
            value = [v for v in self.classes[item]]

        if isinstance(item, slice):
            value = []
            for i in range(len(self.header)).__getitem__(item):
                v = self._get_colvalue(i)
                value.append(v)
        else:
            value = self._get_colvalue(item)

        
        return value

    def __str__(self):
        return '\t'.join(self.values)
    
    def _old_get_colvalue(self, index):
        
        if index > 8:
            if self.FormatClass is None:
                self.FormatClass = self.classes[8](self.values[8])
            value = self.classes[index](self.values[index], self.FormatClass.value)
        else:
            value = self.classes[index](self.values[index])
            
        return value

    def _get_colvalue(self, index):

        Format = self.classes[8]
        Format.value = self.values[8]

        if index > 8:
            value = self.classes[index](self.values[index], Format)
        else:
            value = self.classes[index](self.values[index])
            
        return value

    def iter_samples(self):

        Format = self.classes[8]
        Format.value = self.values[8]
        
        maxRange = self.nsamples + 9
        
        #FormatValue = self.FormatClass.value
        SampleClass = self.classes[9]

        for i in xrange(9, maxRange):
            yield SampleClass(self.values[i], Format)
        
        #return (SampleClass(self.values[i], Format) for i in xrange(9, maxRange))      
        
    def old_iter_samples(self):
        
        self.FormatClass = self.classes[8](self.values[8])
        maxRange = self.nsamples + 9
        
        FormatValue = self.FormatClass.value
        SampleClass = self.classes[9]
        
        return (SampleClass(self.values[i], FormatValue) for i in xrange(9, maxRange))


###### ######

class GFFRecord():
    def __init__(self):
        self.exons = []
        self.name = None
        
    def add(self, items):
        if items[2] == "exon":
            self._exon(items)
        
        if self.name is None:
            self.name = items[0]
            
        try:
            self.records.append(items)
        except AttributeError:
            self.records = [items]
            
    def _exon(self, record):
        self.exons.extend([int(record[3]), int(record[4])])

###### ######

from evolib.data.AlignmentObjects import IOPopulationData

class msFormat(IOPopulationData):
    
    
    def __init__(self, text):
        self.text = text
        lines = [i for i in text.split('\n')[2:] if i != '']
        self.Seqs = lines
        #self.IOdata = self._alt_get_IOdata(lines)
        self.IOdata = self._get_IOdata(lines)
        
            
    def __str__(self):
        return self.text
    
    #def nsamples(self):
    #    
    #    try:
    #        n = len(self.IOdata[0])
    #    except IndexError:
    #        n = None
    #        
    #    return n
    
    def sample_sites(self, p):
        
        newIO = BinaryTable()
        for i in self.IOdata:
            pick = random.random()
            if pick < p:
                newIO.append(i)
                
        self.IOdata = newIO



###### ######

from evolib.formats.ParseMethods import parse_fasta_alignment
from evolib.generic.AlignmentSite import FastaSite
from evolib.generic.GeneticSequence import FastaSequence
from evolib.data.AlignmentObjects import DnaPopulationData
from evolib.data.DataObjects import BinaryTable, SeqTable, SeqTable2

class FastaAlignment(DnaPopulationData):
    """
    *FastaFormat* - class representation of DNA sequence data in fasta format.
        evolib.SequenceFormats::FastaFormat
    
    Example usage:
    
       >>> from evolib.SequenceFormats import FastaFormat
       
       Example 1:
       >>> fileObject = open('example.fsa', 'r')
       >>> F = FastaFormat(fileObject)
       
       Example 2:
       >>> import sys
       >>> stdinObject = sys.stdin
       >>> F = FastaFormat(stdinObject)
       
    """
    
    def __init__(self, fileObject):
        
        self._from_file(fileObject)


    def _from_file(self, fileObject):
        
        seqs, ids = parse_fasta_alignment(fileObject)
        self._attach_data(seqs, ids)

    ######
        
    def __getitem__(self, item):

        if isinstance(item, str):
            Fseq = FastaSequence(self.DNAdata2[item], item)
        elif isinstance(item, int):
            Fseq = FastaSequence(self.DNAdata2[item], self.DNAdata2.ids[item])
        else:
            raise TypeError, "String or integer required"
        
        return Fseq
    

    def __iter__(self):
        
        for i in self.iter_seqs():
            yield i
                

    def __str__(self):

        stringlist = (str(fseq) for fseq in self.iter_seqs())

        return '\n'.join(stringlist)

    ######

    def iter_seqs(self):
        
        nseq = self.nsamples()
        
        for i in xrange(nseq):
            sequence = FastaSequence(self.DNAdata2[i], seqID = self.DNAdata2.ids[i])
            
            yield sequence

    ######
            
    def ids(self):
        return self.DNAdata2.ids


    def nsamples(self):
        return self.__len__()


    def sequences(self):
        return self.DNAdata2.sequences

    
    #def length(self):
    #    return self.validSites
    
