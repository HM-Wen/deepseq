#!/usr/bin/env python 
"""
fetch the unique reads from aligment file
"""

import os 
import sys 
import pysam 

def index_bam(bam_file):
    """
    """
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)


def __main__():
    try:
        in_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1) 
    
    index_bam(in_file) 
    sam_reader = pysam.Samfile(in_file, "rb")
    sam_reader.close()


if __name__=="__main__":
    __main__()
