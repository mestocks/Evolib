from GeneralMethods import loopByColumn

###### ######

def minorMajorAllele(seq):
        
        useq = list(set(seq))
        counts = [(seq.count(i), i) for i in useq]
        counts.sort()
        minor = counts[0][1]
        major = counts[1][1]
        
        return minor, major

def binarizeDNA(DNA):
    
    bases = set(['A', 'T', 'C', 'G'])
    uDNA = set(DNA)
    
    assert uDNA <= bases
    
    if len(uDNA) == 1:
        robotDNA = DNA.replace(DNA[0], '0')
    else:
        minor, major = minorMajorAllele(DNA)
        robotDNA = DNA.replace(minor, '1')
        robotDNA = robotDNA.replace(major, '0')
    
    return robotDNA

###### ######

class IOPolyTable(list):
    
    def add_sample(self, item):
        nsamples = len(item)
        for i in range(nsamples):
            try:
                self[i] += item[i]
            except IndexError:
                self.append(item[i])

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

class SeqTable(list):
    
    def seqsBySite(self):
        
        for site in loopByColumn(self):
            SiteClass = Site(site)
            yield SiteClass
