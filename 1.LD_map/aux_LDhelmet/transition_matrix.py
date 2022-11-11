import sys 
import numpy as np


print('Usage is: transition_matrix.py outgroup')
print('Script also assumes that the outgroup file has only overlapping variants between the outgroup and the samples')

outgroupfile = open(sys.argv[1], 'r')

transmat = np.zeros( (4,4) ) # transition matrix 

order = ['A', 'C', 'G', 'T'] # will help fill in matrix without many if statements!

for line in outgroupfile:
    tmp = line.rstrip().split('\t') # first splitter
    allele1 = tmp[4].split(':')[0] # nucleotide 1
    counts1 = int(tmp[4].split(':')[1]) # counts 
    allele2 = tmp[5].split(':')[0] # nucleotide 2 
    counts2 = int(tmp[5].split(':')[1]) # counts 
    if counts1 > counts2: # define ancestral - derived based on frequency 
        ancestral = allele1
        derived = allele2
    else:
        ancestral = allele2
        derived = allele1
    transmat[order.index(ancestral)][order.index(derived)] += 1  # add a +1  at the correct position

# after all is done move to frequencies 


print(transmat)

