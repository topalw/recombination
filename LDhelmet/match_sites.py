#!/usr/bin/env python3

import sys

print('Usage is match_sites.py outgroup1 outgroup2 output_outgroup1 output_outgroup2 samples')

outgroup1 = open(sys.argv[1], 'r')
outgroup2 = open(sys.argv[2], 'r')
out_outg1 = open(sys.argv[3], 'w')
out_outg2 = open(sys.argv[4], 'w')
samples = open(sys.argv[5], 'r')
common = open('common.match', 'w')

outgroup1_ls = []
outgroup2_ls = []
samples_ls = []

# get chr:pos from outgroups

for line in outgroup1:
    tmp = line.rstrip()
    tmp = tmp.split('\t')
    outgroup1_ls.append(':'.join(tmp))
print('Read outgroup 1')

for line in outgroup2:
    tmp = line.rstrip()
    tmp = tmp.split('\t')
    outgroup2_ls.append(':'.join(tmp))
print('Read outgroup 2')
# now go through the sample file and match positions is too slow
# write file in list and check lists might be faster

for line in samples:
    tmp = line.split('\t')
    tmp2 = ':'.join([tmp[0],tmp[1]])
    samples_ls.append(tmp2)
# now each list has elements in form chr:position
# find matching elements

# we need to split and write them in tab delimited
# this matching step is faster than matching in each line of samples file

### STILL TAKING TOO MUCH TIME ###
### SORT LISTS AND SEARCH BISECT OR USE SETS ###
set1 = set(outgroup1_ls)
set2 = set(outgroup2_ls)
samples_set = set(samples_ls)

print('Writing file 1')

for element in set1:
    if element in samples_set:
        out_outg1.write(''.join(['\t'.join(element.split(':')), '\n']))
    else:
        continue

print('Writing file 2')

for element in set2:
    if element in samples_set:
        out_outg2.write(''.join(['\t'.join(element.split(':')), '\n']))
    else:
        continue

print('DONE')
