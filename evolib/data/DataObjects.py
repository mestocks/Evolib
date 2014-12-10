import numpy

class IOtable(numpy.ndarray):
    """
    Boolean representation of segregating sites in a population.
    IOtable[i][j] represents the ith segregating sites in the 
    jth individual.
    """
    def __new__(cls, input_array):
        return numpy.asarray(input_array).view(cls)

    def __array_finalize__(self, obj):
        if obj is None: return

    def minor_allele_frequencies(self):
        return numpy.sum(self, axis = 1)

    def site_frequency_spectrum(self):

        maf = self.minor_allele_frequencies()
        n = len(self[0])
        
        if n % 2 == 0:
            minlen = (n / 2) + 1
        else:
            minlen = (n + 1) / 2
        
        return numpy.bincount(maf, minlength = minlen)

###### ######

class BinaryTable(list):
    
    def add_sample(self, item):
        s = len(item)
        
        #if self.nsamples() is not None:
        #    assert s == self.seg_sites(), "Cannot add a different number of segregating sites."
        
        for i in range(s):
            try:
                self[i] += item[i]
                
            except IndexError:
                self.append(item[i])


###### ######

from evolib.generic.AlignmentSite import Site
from evolib.tools.GeneralMethods import loopByColumn

class SeqTable(object):
    
    def __init__(self, sequences, ids):
        self.sequences, self.ids = sequences, ids
        self.seqlookup = dict([(ids[i], sequences[i]) for i in xrange(len(ids))])

    def __getitem__(self, index):

        if isinstance(index, str):
            item = self.seqlookup[index]
        elif isinstance(index, int):
            item = self.seqlookup[self.ids[index]]
        else:
            raise TypeError, "String or integer required"

        return item

    def __len__(self):
        return len(self.sequences)

    def seqsBySite(self):
        
        for site in loopByColumn(self.sequences):
            SiteClass = Site(site)
            yield SiteClass

            
###### ######
#, block_iter, member_iter
#from evolib.tools.DNAmethods import minorMajorAllele, binarizeDNA, sites2codons, synNonsynProbs

"""
class SequenceData():

    Class representation of DNA sequences from multiple 
    samples. 

    def __init__(self, seqs, ids = None):
        
        if isinstance(seqs, list) is False:
            raise TypeError, 'List expected.'
        
        n = len(seqs)
        
        if isinstance(ids, list) is False:
            if ids is None:
                ids = self._create_ids(n)
            else:
                raise TypeError, 'List expected.'
            
        self._from_sequence(self, seqs, ids)
        self.pop_nsam = None


    def _create_ids(self, nseqs):

        Creates a list of sequence ids equal to the number 
        of sequences as shown in the examples below:
           ['seq1', 'seq2', 'seq3']
           ['seq01', 'seq02', ...'seq24']

        ids = []
        for i in range(nseqs):
            seqnum = i + 1
            num_zeros = len(str(nseqs)) - len(str(seqnum))
            zeros = '0' * num_zeros
            id1 = 'seq' + zeros + str(seqnum)
            ids.append(id1)
            
        return ids
    
    
    def _from_sequence(self, seqs, ids):
        
        self.seqs = seqs
        self.ids = ids
        
        self.Seqs = SeqTable(seqs)
        self.IO = self._getBinaryTable(self.Seqs)
    
    
    def _getBinaryTable(self, seqs):
        
        self.validSites = 0
        IO = BinaryTable()
        self.sites = []
        for site in seqs.seqsBySite():
            self.sites.append(site)
            if site.hasMissingData():
                pass
            elif site.numberOfAlleles() > 2:
                pass
            elif site.numberOfAlleles() == 1:
                self.validSites += 1
            else:
                self.validSites += 1
                siteIO = binarizeDNA(site.alleles())
                IO.append(siteIO)
                
        return IO
    
    
    def annotate(self, exons):
        self.exons = exons
    
    
    def codingBySite(self, start = 1):
        
        amino_acids = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 
                       'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 
                       'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 
                       'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
                       'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 
                       'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 
                       'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 
                       'TAG':'*', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
                       'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 
                       'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C', 
                       'TGA':'*', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 
                       'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 
                       'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
        
        for i in block_iter(self.bySite(), start = start):
            SiteObject = i[0]
            SiteObject.siteClassBroad = 'Coding'
            
            if len(i[1]) == 3:
                Pos1 = i[1][0]
                Pos2 = i[1][1]
                Pos3 = i[1][2]
                codon_sites = [Pos1, Pos2, Pos3]
                codons = sites2codons(Pos1.alleles(), Pos2.alleles(), Pos3.alleles())
                ucodons = set(codons)
                SiteObject.ucodons = ucodons
                SiteObject.uaminos = None
                
                if sum([b.hasMissingData() for b in codon_sites]) == 0:
                    aminos = [amino_acids[j] for j in ucodons]
                    uaminos = set(aminos)
                    SiteObject.uaminos = uaminos
                    if len(ucodons) == 1:
                        
                        SiteObject.siteClassNarrow = 'MONO'
                        yield SiteObject
                        
                    elif len(ucodons) > 2:
                        
                        SiteObject.siteClassNarrow = 'MULTI'
                        yield SiteObject
                        
                    else:
                        
                        if SiteObject.numberOfAlleles() > 1:
                            
                            aminos = [amino_acids[j] for j in ucodons]
                            uaminos = set(aminos)
                            
                            if '*' in uaminos:
                                
                                SiteObject.siteClassNarrow = 'STOP'
                                yield SiteObject
                                
                            elif len(uaminos) == 1:
                                
                                SiteObject.siteClassNarrow = 'SYN'
                                yield SiteObject
                                
                            elif len(uaminos) == 2:
                                
                                SiteObject.siteClassNarrow = 'NONSYN'
                                yield SiteObject
                                
                            else:
                                raise IndexError, 'This should not happen (1).'
                            
                        else:
                            SiteObject.siteClassNarrow = 'MONO'
                            yield SiteObject
                            
                else:
                    SiteObject.siteClassNarrow = None
                    yield SiteObject
            else:
                SiteObject.siteClassNarrow = None
                #yield SiteObject

    
    def bySite(self):
        
        for site in self.Seqs.seqsBySite():
            yield site
    
    
    def define_pops(self, pop_nsam):
        self.pop_nsam = pop_nsam
        
    
    def getSite(self, index):
        return self.sites[index]
        
        
    def justSynonymous(self):
        
        member = member_iter(self.exons, end = 10000)
        coding = self.codingBySite()
        
        IO = BinaryTable()
        self.validSites = 0
        
        while True:
            try:
                isCoding = member.next()
                if isCoding is True:
                    
                    try:
                        Bp = coding.next()
                        if Bp.siteClassNarrow == 'SYN':
                            siteIO = binarizeDNA(Bp.alleles())
                            IO.append(siteIO)
                            
                        if Bp.siteClassNarrow is not None:
                            probs = [synNonsynProbs(codon)[0] for codon in Bp.ucodons]
                            prob = sum(probs) / len(probs)
                            self.validSites += prob
                                
                    except StopIteration:
                        break
                    
            except StopIteration:
                break
                
        self.IO = IO
        
    
    def populations(self):
        pops = self.pop_nsam
        
        start = 0
        for i in pops:
            finish = start + i
            IO = BinaryTable()
            for io in self.IO:
                IO.append(io[start: finish])
            start = finish
            
            yield IO
            
        
    def seg_sites(self):
        return self.IO.seg_sites()
    
    def thetaW(self):
        return self.IO.thetaW()
    
    def thetaPi(self):
        return self.IO.thetaPi()
    
    def tajD(self):
        return self.IO.tajD()
        
    def wh97(self):
        assert self.pop_nsam is not None
        return self.IO.wh97(self.pop_nsam)
        
"""
