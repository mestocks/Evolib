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

# Site object inheritance not compatible
#from evolib.generic.AlignmentSite import VCFSite

class VCFrow(object):
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
        
        # 1 item costs 0.7s per 1 million rows (str lookup)
        # if statement along costs 0.5s per 1 million row (str lookup)
        if isinstance(item, str):
            itemClass, index = self.header.str_item[item]
        elif isinstance(item, int):
            itemClass, index = self.header.int_item[item]
        
        if index <= 8:
            # initialising an empty str class costs ~0.6s per 1 million rows (CHROM)
            value = itemClass(self.values[index])
            #value = self.values[index]
        else:
            value = itemClass(self.values[index], self.Format)

        return value

    def smap(self, mode, key):
        # currently ~ 8,500 rows/second (inc. CHROM and POS)
        return map(mode, (col.str_fetch(key) for col in self.iter_samples()))

    def iter_samples(self):
        
        maxRange = self.nsamples + 9
        SampleClass, index = self.header.int_item[9]

        for i in xrange(9, maxRange):
            yield SampleClass(self.values[i], self.Format, self.header.names[i])
        


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
from evolib.data.DataObjects import SeqTable

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
    
    def __init__(self, *args):
        
        nargs = len(args)
        if nargs == 1:
            if isinstance(args[0], list):
                seqs, ids = self._from_sequence(args[0])
            else:
                seqs, ids = parse_fasta_alignment(args[0])
        elif nargs == 2:
                seqs, ids = self._from_sequence(args[0], args[1])
        else:
            raise TypeError, "Wrong number of arguments"
        
        self._attach_data(seqs, ids)


    ######
        
    def __getitem__(self, item):

        if isinstance(item, str):
            Fseq = FastaSequence(self.DNAdata[item], item)
        elif isinstance(item, int):
            Fseq = FastaSequence(self.DNAdata[item], self.DNAdata.ids[item])
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

    def pop(self, index = None):
        
        seq, seqid = self.DNAdata.pop(index)
        self.IOdata = self._get_IOdata(self.DNAdata)

        return FastaSequence(seq, seqid)

    ######

    def iter_seqs(self):
        
        nseq = self.nsamples()
        
        for i in xrange(nseq):
            sequence = FastaSequence(self.DNAdata[i], seqID = self.DNAdata.ids[i])
            
            yield sequence


    
