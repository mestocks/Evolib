BEGIN {

    # A[i * win] = i
    # B[i] = i * win

    #min = 0
    #max = 1
    range = max - min
    win = range / nbreaks
    A[min] = win
    
    for (i = 1; i <= nbreaks; i++) {
	A[i * win] = i
	B[i] = i * win
    }
  }
{

    d = $c / win
    split(d,e,".")
    C[e[1]+1]++
}
END {

    for (j = 1; j <= nbreaks; j++) {
	if (j in C) {
	    print j,B[j],C[j]
	} else {
	    print j,B[j],"0"
	}
    }

}
