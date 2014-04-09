
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


def sites2codons(site1, site2, site3):
	
	nsamples = len(site1)
	codons = []
	for i in range(nsamples):
		codon = site1[i] + site2[i] + site3[i]
		codons.append(codon)
		
	return codons


def synNonsynProbs(codon):
	
	dna = set(['A', 'T', 'G', 'C'])
	
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
	
	s, n = 0, 0
	AA = amino_acids[codon]
	
	for i in range(3):
		base_to_replace = codon[i]
		diff = dna - set(base_to_replace)
		
		for j in diff:
			new_key = [k for k in codon]
			new_key[i] = j
			aa = amino_acids[''.join(new_key)]
			
			if aa != AA:
				n += 1
			else:
				s += 1
	
	return s / 9.0, n / 9.0
	

def codonSynNonsyn():
	dna = set(['A', 'T', 'G', 'C'])
	
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
	total = 0
	S, N = 0, 0
	for key in amino_acids.keys():
		AA = amino_acids[key]
		print key, AA,
		s, n = 0, 0
		for i in range(3):
			base_to_replace = key[i]
			diff = dna - set(base_to_replace)
			for j in diff:
				new_key = [k for k in key]
				new_key[i] = j
				aa = amino_acids[''.join(new_key)]
				if aa != AA:
					n += 1
					N += 1
				else:
					s += 1
					S += 1
				total += 1
		print s, n, s / 9.0, n / 9.0
	print '--------------'
	print 'Total', '-', S, N, S / float(total), N / float(total)
