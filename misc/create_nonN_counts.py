from __future__ import print_function
import sys
import collections

ref = sys.argv[1]

counts = collections.OrderedDict()

with open(ref, "r") as f:
    for line in f:
        l = line.rstrip()
        if l[0]==">":
            chrom = l.split()[0][1:]
            counts[chrom] = collections.defaultdict(int)
        else:
            for letter in l:
                counts[chrom][letter] += 1

for chrom in counts:
    print("\t".join([chrom, str(sum(counts[chrom].values())), str(sum([counts[chrom][k] for k in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']]))]))
