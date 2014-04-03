# Classes
from lib.Stats import PopStats
from lib.DNAmethods import IOPolyTable, SeqTable
from lib.DNAobjects import FastaSequence

# Methods
from lib.DNAmethods import minorMajorAllele, binarizeDNA

###### ######

class FastaFormat(PopStats):
    """
    Example usage:
    
       >>> from evolib.SequenceFormats import FastaFormat
       
       Example 1:
       >>> fileObject = open('example.fsa', 'r')
       >>> F = FastaFormat(fileObject)
       
       Example 2:
       >>> import sys
       >>> stdinObject = sys.stdin
       >>> F = FastaFormat(stdinObject)
       
       Example 3:
       >>> seq = 'ATCTGATGCTGAC'
       >>> F = FastaFormat(seq)
       
       Example 4:
       >>> id = 'seq1'
       >>> seq = 'ATCGATGTCGTGAC'
       >>> F = FastaFormat(seq, id)
       
       Example 5:
       >>> ids = ['seq1', 'seq2', 'seq3']
       >>> seqs = ['ATCGATGTCGTGAC', 'ATCGATGTCGTGAC', 'ATCGATGTCGTGAC']
       >>> F = FastaFormat(seqs, ids)
    """

    def __init__(self, *args):
        
        if len(args) == 1:
            arg1 = args[0]
            if isinstance(arg1, file):
                self._fromFile(arg1)
            elif isinstance(arg1, str):
                self._fromSequence(None, [arg1])
            elif isinstance(arg1, list):
                self._fromSequence(None, arg1)
            else:
                raise TypeError, 'Wrong arg type. File object or list of sequences.'
        elif len(args) == 2:
            arg1 = args[0]
            arg2 = args[1]
            if isinstance(arg1, str):
                assert isinstance(arg2, str)
                self._fromSequence([arg2], [arg1])
            elif isinstance(arg1, list):
                assert isinstance(arg2, list)
                assert len(arg1) == len(arg2)
                self._fromSequence(arg2, arg1)
            else:
                raise TypeError, 'Wrong arg type. File object or list of sequences.'
            
    def __iter__(self):
        
        nseq = len(self.sequences)
        for i in range(nseq):
            sequence = FastaSequence(self.sequences[i], seqID = self.ids[i])
            yield sequence
            
    def __str__(self):
        
        stringList = []
        for s in range(len(self.sequences)):
            seq = FastaSequence(self.sequences[s], seqID = self.ids[s])
            stringList.append(seq)
        theString = '\n'.join(map(str, stringList))
            
        return theString

    def _fromFile(self, fileObject):
        """
        Returns a matrix where M[i][j] refers to the jth site of the ith individual.
        """        
        step1 = ''.join(map(lambda line: line, fileObject))
        step2 = step1.split('\n>')
        seq_table = [part.partition('\n')[2].replace('\n','') for part in step2]
        seq_names = [part.partition('\n')[0].replace('>', '') for part in step2]
        self.ids = seq_names
        
        self.sequences = SeqTable(seq_table)
        self.IOlen = 0
        self.IOtable = IOPolyTable()
        for site in self.sequences.seqsBySite():
            if site.hasMissingData():
                pass
            elif site.numberOfAlleles() > 2:
                pass
            elif site.numberOfAlleles() == 1:
                self.IOlen += 1
            else:
                self.IOlen += 1
                siteIO = binarizeDNA(site)
                self.IOtable.append(siteIO)

    def _fromSequence(self, ids, seqs):
        
        if ids is None:
            nseqs = len(seqs)
            ids = self._create_ids(nseqs)
        
        self.ids = ids
        self.sequences = SeqTable(seqs)

    def _create_ids(self, nseqs):
        """
        Creates a list of sequence ids equal to the number 
        of sequences as shown in the examples below:
           ['seq1', 'seq2', 'seq3']
           ['seq01', 'seq02', ...'seq24']
        """
        ids = []
        for i in range(nseqs):
            seqnum = i + 1
            num_zeros = len(str(nseqs)) - len(str(seqnum))
            zeros = '0' * num_zeros
            id1 = 'seq' + zeros + str(seqnum)
            ids.append(id1)
            
        return ids
    
    def length(self):
        return self.IOlen
    
    def nsamples(self):
        return len(self.sequences)

### ###

class msFormat(PopStats):
    
    def __init__(self, text):
        self.text = text
        lines = [i for i in text.split('\n')[2:] if i != '']
        
        self.IOtable = IOPolyTable()
        for line in lines:
            self.IOtable.add_sample(line)
            
    def __str__(self):
        return self.text
    
    def nsamples(self):
        
        try:
            nsamples = len(self.IOtable[0])
        except IndexError:
            nsamples = None
            
        return nsamples
