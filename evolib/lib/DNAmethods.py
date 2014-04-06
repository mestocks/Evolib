

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
