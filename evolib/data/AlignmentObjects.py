from evolib.tools.GeneralMethods import create_ids

from evolib.stats.StatObjects import IOstats
from evolib.generic.AlignmentSite import Site

from evolib.data.DataObjects import SeqTable, IOtable
from evolib.tools.DNAmethods import booleanDNA, booleanIO
from evolib.generic.GeneticSequence import DNAsequence

from evolib.tools.GeneralMethods import loopByColumn

from evolib.tools.DNAmethods import dna_to_amino, synNonsynProbs

class DnaPopulationData(IOstats):
    """
    Class representation of DNA sequences from multiple 
    samples. 
    """
    def __init__(self, *args):

        if len(args) == 1:
            seqs, ids = self._from_sequence(args[0])
        elif len(args) == 2:
            seqs, ids = self._from_sequence(args[0], args[1])
        else:
            raise TypeError, "Wrong number of arguments"
        
        self._attach_data(seqs, ids)


    def _from_sequence(self, seqs, ids = None):
    
        if isinstance(seqs, list) is False:
            raise TypeError, 'List expected.'
        
        n = len(seqs)
        
        if isinstance(ids, list) is False:
            if ids is None:
                ids = create_ids(n, "seq")
            else:
                raise TypeError, 'List expected.'

        return seqs, ids
        

    def _attach_data(self, sequences, ids):
        
        self.DNAdata = SeqTable(sequences, ids)
        self.IOdata = self._get_IOdata(self.DNAdata)

    ######
        
    def __len__(self):
        return len(self.DNAdata)
        
    ######
    
    def _get_IOdata(self, seqs):

        self.validSites = 0
        io = []
        
        for site in seqs.iter_sites():

            SiteClass = Site(site)
            
            if SiteClass.has_missing_data():
                pass
            elif SiteClass.number_of_alleles() > 2:
                pass
            elif SiteClass.number_of_alleles() == 1:
                self.validSites += 1
            else:
                self.validSites += 1
                siteIO = booleanDNA(SiteClass.alleles())
                io.append(siteIO)

        IO = IOtable(io)
                
        return IO

    ######
    
    def iter_sites(self):
        for site in self.DNAdata.iter_sites():
            yield site
    
    ######
            
    def ids(self):
        return self.DNAdata.ids


    def nsamples(self):
        return self.__len__()


    def sequences(self):
        return self.DNAdata.sequences

    
    def length(self):
        return self.validSites

    ######

    def index(self, key):
        return self.DNAdata.index(key)

    def pop(self, index = None):
        seq, seqid = self.DNAdata.pop(index)
        self.IOdata = self._get_IOdata(self.DNAdata)
        
        return DNAsequence(seq, seqid)

    ######

    def coding(self, refseq):

        dna = ['A', 'T', 'G', 'C', 'a', 't', 'g', 'c']
        nsam = self.nsamples()
        
        inc = (i for i in xrange(len(refseq)) if refseq[i] in dna)
        cds_seqs = ['' for j in xrange(nsam)]

        for site in inc:
            for ind in xrange(nsam):
                cds_seqs[ind] += self.DNAdata.sequences[ind][site]
            
        return type(self)(cds_seqs, self.DNAdata.ids)

    ######

    def nonsyn(self, frame):

        nsyn, nnon = 0, 0
        nsam = len(self.DNAdata.sequences)
        syn_seqs, nonsyn_seqs = ['' for n in range(nsam)], ['' for n in range(nsam)]
        for codons in loopByColumn(self.DNAdata.sequences, start = frame, size = 3):

            nmiss = sum([len(set(i) - set('ATGCatgc')) for i in codons])
            if nmiss > 0:
                # contains non-ATGC data
                pass
            
            else:
                
                ucodons = list(set(codons))
                nucodons = len(ucodons)
                if nucodons > 2:
                    # > 1 segregating site in codon
                    pass
                
                elif nucodons == 1:
                    # monomorphic site
                    nsp, nnp = synNonsynProbs(ucodons[0])
                    nsyn += 3 * nsp
                    nnon += 3 * nnp
                    
                else:
                    codon1, codon2 = ucodons[0], ucodons[1]
                    codon_count = [(codons.count(codon1), codon1), (codons.count(codon2), codon2)]
                    codon_count.sort(reverse = True)
                    major = codon_count[0][1]
                    nsp, nnp = synNonsynProbs(major)
                    nsyn += 3 * nsp
                    nnon += 3 * nnp
                    
                    sindex = [i for i in range(len(codon1)) if codon1[i] != codon2[i]][0]
                    for s in range(len(codons)):
                        aa1, aa2 = dna_to_amino[ucodons[0]], dna_to_amino[ucodons[1]]
                        if aa1 == aa2:
                            syn_seqs[s] += codons[s][sindex]
                        else:
                            nonsyn_seqs[s] += codons[s][sindex]

        SynClass = type(self)(syn_seqs, self.DNAdata.ids)
        NonSynClass = type(self)(nonsyn_seqs, self.DNAdata.ids)

        SynClass.validSites = nsyn
        NonSynClass.validSites = nnon

        return SynClass, NonSynClass
            

###### ######

import random

class IOPopulationData(IOstats):

    def __init__(self, seqs):

        if isinstance(seqs, IOtable):
            self.IOdata = seqs
        else:
            self.IOdata = self._get_IOdata(seqs)


    def _get_IOdata(self, seqs):
        
        io = []
        if seqs != []:
            for s in range(len(seqs[0])):
                site = ''.join([f[s] for f in seqs[:]])
                bio = booleanIO(site)
                io.append(bio)
                
        IO = IOtable(io)
        
        return IO

    def nonsyn_sample_sites(self, p):

        n = len(self.IOdata)
        one, two = [], []
        
        for i in xrange(n):
            rint = random.random()

            if rint < p:
                one.append(list(self.IOdata[i]))
            else:
                two.append(list(self.IOdata[i]))

        table1 = IOtable(one)
        table2 = IOtable(two)

        OneClass, TwoClass = IOPopulationData(table1), IOPopulationData(table2)

        return OneClass, TwoClass

        
