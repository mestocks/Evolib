{
    A[$e] += $d
}
END {
    for (factor in A) {
	print factor,A[factor]
    }
}
