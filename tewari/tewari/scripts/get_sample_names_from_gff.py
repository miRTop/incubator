import sys
import os

fns = list(sys.argv[1:])


for fn in fns:
    with open(fn) as inh:
        for line in inh:
            if line.find("COLDATA") > 0:
                sample = line.strip().split()[-1]
                print "%s,%s" % (os.path.basename(fn), sample)
