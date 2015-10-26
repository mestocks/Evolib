import sys
from evolib.iterators import vcf_iter

config_file = open(sys.argv[1], 'r')
cols, direction, cutoff = [], [], []
for line in config_file:
    cfg_list = line.rstrip().split(",")
    cols.append(cfg_list[0])
    direction.append(cfg_list[1])
    cutoff.append(float(cfg_list[2]))

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
