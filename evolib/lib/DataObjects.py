import PopGenStats

# Methods
from GeneralMethods import loopByColumn, block_iter
from DNAmethods import minorMajorAllele, binarizeDNA, sites2codons

###### ######

class Site():
    
    def __init__(self, alleles):
        self._alleles = alleles
        
    
    def alleles(self):
        return self._alleles
    
    
    def hasMissingData(self, dna = ['A', 'T', 'C', 'G']):
        
        if set(self.alleles()) <= set(dna):
            answer = False
        else:
            answer = True
            
        return answer
    
    
    def numberOfAlleles(self):
        alleles = set(self.alleles())
        return len(alleles)
    
    
    def nonsyn_site(self):
        return
    
    
    def site_class(self):
        return 'Exon'


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
                
    def seg_sites(self):
        
        if self == []:
            s = 0
        else:
            s = len(self)
            
        return s
    
    def thetaW(self):
        n = self.nsamples()
        s = self.seg_sites()
        
        return PopGenStats.WattersonsTheta(n, s)
        
    def thetaPi(self):
        pi = 0.0
        
        if self.nsamples() is not None:
            n = self.nsamples()
            pi = PopGenStats.TajimasTheta(n, self)
        
        return pi
    
    def tajD(self):
        n = self.nsamples()
        s = self.seg_sites()
        tw = self.thetaW()
        pi = self.thetaPi()
        
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
        
        self.Seqs = SeqTable(seqs)
        self.IO = self._getBinaryTable(self.Seqs)
    
    
    def _getBinaryTable(self, seqs):
        
        self.validSites = 0
        IO = BinaryTable()
        for site in seqs.seqsBySite():
            
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
    
    
    def codingBySite(self):
        
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
        
        for i in block_iter(self.bySite()):
            SiteObject = i[0]
            SiteObject.siteClassBroad = 'Coding'
            Pos1 = i[1][0]
            Pos2 = i[1][1]
            Pos3 = i[1][2]
            codon_sites = [Pos1, Pos2, Pos3]
            codons = sites2codons(Pos1.alleles(), Pos2.alleles(), Pos3.alleles())
            ucodons = set(codons)
            
            if len(ucodons) == 1:
                
                SiteObject.siteClassNarrow = None
                yield SiteObject
                
            elif len(ucodons) > 2:
                
                SiteObject.siteClassNarrow = None
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
                    SiteObject.siteClassNarrow = None
                    yield SiteObject

    
    def bySite(self):
        
        for site in self.Seqs.seqsBySite():
            yield site
    
    
    def define_pops(self, pop_nsam):
        self.pop_nsam = pop_nsam
        
    
    def populations(self):
        pops = self.pop_nsam
        
        start = 0
        for i in pops:
            finish = start + i
            newTable = SeqTable(self.Seqs[start: finish])
            newBinary = self._getBinaryTable(newTable)
            
            start = finish
            
            yield newBinary
            
        
    def seg_sites(self):
        return self.IO.seg_sites()
    
    def thetaW(self):
        return self.IO.thetaW()
    
    def thetaPi(self):
        return self.IO.thetaPi()
    
    def tajD(self):
        return self.IO.tajD()
