import PopGenStats

# Methods
from GeneralMethods import loopByColumn
from DNAmethods import minorMajorAllele, binarizeDNA

###### ######

class Site(str):
    
    def hasMissingData(self, dna = ['A', 'T', 'C', 'G']):
        
        if set(self) <= set(dna):
            answer = False
        else:
            answer = True
            
        return answer
    
    def numberOfAlleles(self):
        alleles = set(self)
        return len(alleles)

###### ######

class BinaryTable(list):
    
    def add_sample(self, item):
        s = len(item)
        
        if self.nsamples() is not None:
            assert s == self.segregating_sites(), "Cannot add a different number of segregating sites."
        
        for i in range(s):
            try:
                self[i] += item[i]
                
            except IndexError:
                self.append(item[i])
    
    def nsamples(self):
        
        if self == []:
            n = None
        else:
            n = len(self[0])
            
        return n
                
    def segregating_sites(self):
        
        if self == []:
            s = 0
        else:
            s = len(self)
            
        return s
    
    def thetaw(self):
        n = self.nsamples()
        s = self.segregating_sites()
        
        return PopGenStats.WattersonsTheta(n, s)
        
    def thetapi(self):
        pi = 0.0
        
        if self.nsamples() is not None:
            n = self.nsamples()
            pi = PopGenStats.TajimasTheta(n, self)
        
        return pi
    
    def tajimasd(self):
        n = self.nsamples()
        s = self.segregating_sites()
        tw = self.thetaw()
        pi = self.thetapi()
        
        return PopGenStats.TajimasD(n, s, tw, pi)

###### ######

class SeqTable(list):
    

    def seqsBySite(self):
        
        for site in loopByColumn(self):
            SiteClass = Site(site)
            yield SiteClass

###### ######

class SequenceData():
    
    
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
    
    
    def _from_sequence(self, seqs, ids):
        
        self.seqs = seqs
        self.ids = ids
        
        self._getSeqTable(seqs)
        self._getBinaryTable(self.Seqs)
    
    
    def _getBinaryTable(self, seqs):
        
        self.validSites = 0
        self.IO = BinaryTable()
        for site in seqs.seqsBySite():
            
            if site.hasMissingData():
                pass
            elif site.numberOfAlleles() > 2:
                pass
            elif site.numberOfAlleles() == 1:
                self.validSites += 1
            else:
                self.validSites += 1
                siteIO = binarizeDNA(site)
                self.IO.append(siteIO)
                
    
    def _getSeqTable(self, seqs):
        
        self.Seqs = SeqTable(seqs)
        

    def define_pops(self, npops):
        
