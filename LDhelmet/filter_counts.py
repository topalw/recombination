#!/usr/bin/env python3

import sys

counts = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for line in counts:
    if line[0:3] == 'CHR':
        continue
    else:
        tmp = line.split('\t')
        total = int(tmp[3])
        countREF = int(tmp[4].split(':')[1])
        countALT = int(tmp[5].split(':')[1])
        if countREF > 1 or countREF < int(total-1):
            outfile.write('\t'.join([tmp[0],tmp[1],'\n']))
        else: 
            continue
