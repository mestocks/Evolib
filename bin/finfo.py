
import sys

from evolib.iterators import vcf_iter

cols = ["FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum"]
cutoff = [60, 40, -12.5, 2, -8]
direction = ["lt", "gt", "gt", "gt", "gt"]

n = 0
for row in vcf_iter(sys.stdin):

    if n == 0:
        print row.header.preamble
        n = 1

    allow = True

    for i in xrange(len(cols)):
        c, co, d = cols[i], cutoff[i], direction[i]
        col = row['INFO'][c]

        if col is None:
            allow = False
            break
        else:
            if d == "gt":
                if float(col) >= co:
                    pass
                else:
                    allow = False
                    break
            elif d == "lt":
                if float(col) <= co:
                    pass
                else:
                    allow = False
                    break
    
    if allow is True:
        print row
