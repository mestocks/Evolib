from evolib.tools.GeneralMethods import create_ids

from evolib.stats.StatObjects import IOstats, altIOstats

from evolib.data.DataObjects import BinaryTable, SeqTable, IOtable
from evolib.tools.DNAmethods import binarizeDNA, booleanDNA, booleanIO

class DnaPopulationData(altIOstats):
    """
    Class representation of DNA sequences from multiple 
    samples. 
    """
    def __init__(self, seqs, ids = None):
        
        if isinstance(seqs, list) is False:
            raise TypeError, 'List expected.'
        
        n = len(seqs)
        
        if isinstance(ids, list) is False:
            if ids is None:
                ids = create_ids(n, "seq")
            else:
                raise TypeError, 'List expected.'

        self._attach_data(seqs, ids)
        self.pop_nsam = None


    def _attach_data(self, sequences, ids):
        
        self.DNAdata = SeqTable(sequences, ids)
        self.IOdata = self._alt_get_IOdata(self.DNAdata)

    ######
        
    def __len__(self):
        return len(self.DNAdata)
        
    ######
    
    def _get_IOdata(self, seqs):
        
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

    def _alt_get_IOdata(self, seqs):

        self.validSites = 0
        io = []
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
                siteIO = booleanDNA(site.alleles())
                io.append(siteIO)

        IO = IOtable(io)
                
        return IO


    

###### ######

class IOPopulationData(IOstats):

    def __init__(self, seqs):

        self.IOdata = self._get_IOdata(seqs)
        #self.IOdata = self._alt_get_IOdata(seqs)
        
    def _get_IOdata(self, seqs):
        
        IO = BinaryTable()
        for line in seqs:
            IO.add_sample(line)
            
        return IO

    def _alt_get_IOdata(self, seqs):
        
        io = []
        if seqs != []:
            for s in range(len(seqs[0])):
                site = ''.join([f[s] for f in seqs[:]])
                bio = booleanIO(site)
                io.append(bio)
                
        IO = IOtable(io)
        
        return IO
