#!/usr/bin/env python 
"""
Take a sub sample of reads from paired FASTQ file. 

Usage: python sub_sample_reads.py base_path_with_paired_fastq_files

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
        fastq_path = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    fastq_1_file = None 
    fastq_2_file = None

    # TODO expecting the left and right reads in the base folder with .fastq ending. Make it is common general form

    ## get the files from base path 
    for filename in os.listdir(fastq_path):
        if re.search(r'_1.fastq', filename):
            fastq_1_file = filename 
        if re.search(r'_2.fastq', filename):
            fastq_2_file = filename 
        
    print fastq_1_file, fastq_2_file
    print 
    
    ## count the number of reads and calculate the sub sample 
    fqh = helper.open_file('%s/%s' % (fastq_path, fastq_1_file)) 
    read_cnt = 0
    for rec in SeqIO.parse(fqh, 'fastq'):
        read_cnt +=1 
    fqh.close() 

    print '%d Number of reads in FASTQ' % read_cnt
    print 

    ## what percentage sub-sample 
    percentage = 1
    sub_count = int(round((percentage*read_cnt)/100.0))
    assert sub_count <= read_cnt, ' %d (sub-sample count) should be less than total read count %d'  % (sub_count, read_cnt)
    
    print "%d Sub sample count" % sub_count  
    print 

    try:
        accept_prob = (1.0*sub_count)/read_cnt
    except:
        accept_prob = 1

    print accept_prob

    sub_fastq_1_file = "%d_percentage_%s.bz2" % (percentage, fastq_1_file)
    sub_fastq_2_file = "%d_percentage_%s.bz2" % (percentage, fastq_2_file)

    ## writing out sub sample files
    try:
        subFile_1 = bz2.BZ2File("%s/%s" % (fastq_path, sub_fastq_1_file), 'wb')
        subFile_2 = bz2.BZ2File("%s/%s" % (fastq_path, sub_fastq_2_file), 'wb')
    except Exception as error:
        sys.exit(error)

    total_cnt = 0 
    sample_cnt = 0 
    left_reads = dict() 

    fqh = helper.open_file('%s/%s' % (fastq_path, fastq_1_file)) 
    for rec in SeqIO.parse(fqh, 'fastq'):
        rnb = random.random()
        total_cnt += 1

        if rnb <= accept_prob:
            ## @UNC15-SN850_63:4:1101:1103:2151/1 @UNC15-SN850_63:4:1101:1103:2151/2
            read_id = rec.id.split('/')
            if len(read_id) > 1:
                left_reads[read_id[0]] = 0
            else:    
                ## @UNC11-SN627:294:C236MACXX:5:1101:1430:2218 1:N:0:GGNTAC  @UNC11-SN627:294:C236MACXX:5:1101:1430:2218 2:N:0:GGNTAC 
                left_reads[rec.id] = 0

            sample_cnt += 1 
            subFile_1.write(rec.format("fastq"))

        if sub_count == sample_cnt:
            break 

    fqh.close() 
    subFile_1.close() 

    fqh = helper.open_file('%s/%s' % (fastq_path, fastq_2_file)) 
    for rec in SeqIO.parse(fqh, 'fastq'):
        read_id = rec.id.split('/')

        if len(read_id) > 1:
            ## @UNC15-SN850_63:4:1101:1103:2151/1
            if read_id[0] in left_reads:
                subFile_2.write(rec.format("fastq"))
        else:
            ## @UNC11-SN627:294:C236MACXX:5:1101:1430:2218 1:N:0:GGNTAC
            if rec.id in left_reads:
                subFile_2.write(rec.format("fastq"))
    fqh.close() 
    subFile_2.close() 

    print "%s/%s" % (fastq_path, sub_fastq_1_file)
    print "%s/%s" % (fastq_path, sub_fastq_2_file)
    print 

    print '%d Number of reads scanned' % total_cnt
    print '%d Number of reads in' % sample_cnt
    print 


if __name__ == "__main__":
    __main__()
