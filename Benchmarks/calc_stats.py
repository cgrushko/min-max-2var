#!/usr/bin/python
import numpy as np
import sys
if (len(sys.argv) < 2):
    print "Usage: calc_stats.py FILE_WITH_STATS"
    sys.exit()
filename = sys.argv[1]
print "count\tmean\tmedian\tstd"
with open(filename) as f:
    odd = False
    for line in f:
        odd = not odd
        if odd:
            a=np.array(line.split(),dtype=int)
            if len(a) > 0:             
                print '{0}\t{2:.3f}\t{3}\t{1:.3f}'.format(len(a), a.std(), a.mean(), int(np.median(a)))
