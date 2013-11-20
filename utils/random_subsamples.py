#!/usr/bin/env python 
"""
Take a sub sample of reads from FASTQ file. 

Usage: python random_subsamples.py in.fastq(.bz2|.gz) out.fastq.bz2 number-of_sub_sample_entries

Requirement:
    Biopython:- http://biopython.org
    https://github.com/vipints/genomeutils/blob/master/gfftools/helper.py
"""
import os 
import re
import sys
import bz2 
import random
from Bio import SeqIO
from gfftools import helper 

def __main__():

    try:
        fastqFile = sys.argv[1]
        outFile = sys.argv[2]
        subCnt = int(sys.argv[3]) # sub sample count
    except:
        print __doc__
        sys.exit(-1)
    
    ## open the file
    ffh = helper._open_file(fastqFile)

    ## counting number of reads in fastq file 
    read_cnt = 0
    for rec in SeqIO.parse(ffh, 'fastq'):
        read_cnt +=1 
    ffh.close()
    print 'Number of reads in FASTQ: ', read_cnt

    assert subCnt <= read_cnt, str(subCnt) + ' (sub-sample count) should be less than total read count ' + str(read_cnt)
    try:
        accept_prob = (1.0*subCnt)/read_cnt
    except:
        accept_prob = 1

    ## outfile directory check for creating the new file 
    try:
        subFile = bz2.BZ2File(outFile, 'wb')
    except Exception as error:
        sys.exit(error)

    cnt, sub_cnt = 0, 0
    print 'Writing compressed file...'

    ffh = helper._open_file(fastqFile)
    for rec in SeqIO.parse(ffh, 'fastq'):
        rnb = random.random()
        cnt += 1
        if rnb <= accept_prob:
            sub_cnt += 1 
            subFile.write(rec.format("fastq"))
        if subCnt == sub_cnt:
            print '...done'
            break
    ffh.close()
    subFile.close()

    print 'Number of reads scanned: ', cnt
    print 'Number of reads in: ', sub_cnt

if __name__ == "__main__":
    __main__()
