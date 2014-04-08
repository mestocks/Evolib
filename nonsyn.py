from evolib.lib.DNAmethods import sites2codons
from evolib.lib.GeneralMethods import block_iter

from evolib.NGSFormats import VariantCallFormat

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

myVCF = VariantCallFormat("test/pavy2012_A.vcf")

rowIter = myVCF.bySite()

for i in block_iter(rowIter):
    site = i[0]
    pos1 = i[1][0]
    pos2 = i[1][1]
    pos3 = i[1][2]
    codon_sites = [pos1, pos2, pos3]
    codons = sites2codons(pos1.alleles(), pos2.alleles(), pos3.alleles())
    ucodons = set(codons)
    
    if len(ucodons) == 1:
        pass
    elif len(ucodons) > 2:
        pass
    else:
        if site.numberOfAlleles() > 1:
            aminos = [amino_acids[j] for j in ucodons]
            if len(set(aminos)) == 1:
                print site['POS'], ucodons, aminos, 'SYN'
            elif len(set(aminos)) == 2:
                print site['POS'], ucodons, aminos, 'NONSYN'
    # assess whether there are two many mutations in the codon
    # 
    


"""

x.


"""
