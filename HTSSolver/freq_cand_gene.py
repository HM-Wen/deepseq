#!/usr/bin/env python 
"""
Find the frequency of the candidate genes in different replicates 
TODO fix the input plate format 
"""
import sys, re 

def read_file(fname):
    """
    read a single column file
    """
    gname = dict()
    fh = open(fname, 'rU')
    for line in fh:
        line = line.strip('\n\r')
        if re.match(r'#', line):
            continue     
        gname[line]=1 
    fh.close()
    return gname

try:
    rep_file_1 = sys.argv[1]
    rep_file_2 = sys.argv[2]
    rep_file_3 = sys.argv[3]
except:
    print __doc__
    sys.exit(-1)

genes_1 = read_file(rep_file_1)
genes_2 = read_file(rep_file_2)
genes_3 = read_file(rep_file_3)

#print len(genes_1), len(genes_2), len(genes_3)
master_dict = dict()
for ele in genes_1:
    master_dict[ele]=1 
for ele in genes_2:
    master_dict[ele]=1 
for ele in genes_3:
    master_dict[ele]=1 
#print len(master_dict)

for master in master_dict:
    sys.stdout.write(master+'\t')

    if master in genes_1:
        sys.stdout.write('1\t')
    else:
        sys.stdout.write('0\t')

    if master in genes_2:
        sys.stdout.write('1\t')
    else:
        sys.stdout.write('0\t')

    if master in genes_3:
        sys.stdout.write('1\n')
    else:
        sys.stdout.write('0\n')

