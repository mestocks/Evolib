{
    N[$e]++
    S[$e] += $d
}
END {
    for (factor in N) {
	if (f == "count") {
	    print factor,N[factor]
	}
	else if (f == "sum") {
	    print factor,S[factor]
	}
	else if (f == "mean") {
	    print factor,S[factor]/N[factor]
	}
    }
}
