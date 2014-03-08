# 

BEGIN {
    nloci = 0
    lastLocus = "first"   
}
{
    nsamples = NF - 9
    CHROM = $1
    REF = $4
    ALT = $5
    split($9,FORMAT,":")

    if (CHROM != lastLocus) {
	nloci++
	loci[nloci] = CHROM
	lastLocus = CHROM    
    }
    # Loop through the columns/fields corresponding to each 
    # of the samples. 
    for (s = 1; s <= nsamples; s++) {

	i = s + 9
	split($i,ind,":")

	if (s > 9) {
	    sample = prefix s
	} else {
	    sample = prefix "0" s
	}

	for (j in ind) {

	    if (FORMAT[j] == "GT") split(ind[j],GT,"/")
	    if (FORMAT[j] == "DP") DP = ind[j]
	}
	if (DP < minDP) {
	    A = "N"
	    B = "N"
	    seqArray[CHROM,sample"a"] = seqArray[CHROM,sample"a"] A
	    seqArray[CHROM,sample"b"] = seqArray[CHROM,sample"b"] B
	} else if (ALT == ".") {
	    A = REF
	    B = REF
	    seqArray[CHROM,sample"a"] = seqArray[CHROM,sample"a"] A
	    seqArray[CHROM,sample"b"] = seqArray[CHROM,sample"b"] B
	} else {
	    genos = REF","ALT
	    split(genos,genoArray,",")
	    
	    A = genoArray[GT[1] + 1]
	    B = genoArray[GT[2] + 1]

	    seqArray[CHROM,sample"a"] = seqArray[CHROM,sample"a"] A
	    seqArray[CHROM,sample"b"] = seqArray[CHROM,sample"b"] B
	}
    }
}
END {
    for (comb in seqArray) {

	split(comb,sep,SUBSEP)
	print ">Pa_"sep[1]"_"sep[2],seqArray[sep[1],sep[2]]
    }
}
