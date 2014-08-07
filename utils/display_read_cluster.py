#!/usr/bin/env python 
"""
Finding and displaying short reads clustered in the genome.

Usage: display_read_cluster.py align.bam I:1000-5000 10 
    alignment file in BAM format 
    genomic location where the cluster of reads 
    distance in basepairs for two reads to be in the same cluster

Reuirements: 
"""

import os 
import re 
import sys
import pylab 
import collections
from bx.intervals.cluster import ClusterTree

def alignment_parse(alignment_file, cluster_name):
    
    infh = os.popen('/fml/ag-raetsch/share/software/samtools/samtools view ' + alignment_file + ' ' + cluster_name)
    for line in infh:
        line = line.strip().split('\t')
        yield line[0], line[2], int(line[3]), int(line[3])+41
    infh.close()

def plot_cluster(cluster_name, read_ids, location_db, freq_db):

    read_info = []
    for read_id in read_ids:
        start, end = location_db[read_id]
        freq = freq_db[read_id]
        read_info.append((start, end, freq))
    read_info.sort(reverse=True)
    min_freq = min(l[2] for l in read_info)
    max_freq = max(l[2] for l in read_info)
    print 'Minimum frequency: ' + str(min_freq)
    print 'Maximum frequency: ' + str(max_freq)
    freq_cmap = pylab.get_cmap('Blues')
    pylab.clf()
    labels_done = []
    for rindex, (start, end, freq) in enumerate(read_info):
        if min_freq == max_freq:
            color = 'gray'
        else:
            freq_percent = float(freq - min_freq) / float(max_freq - min_freq)
            color = freq_cmap(freq_percent)
        if freq in [min_freq, max_freq]:
            label = "Frequency: %s" % (freq)
        else:
            label = None
        if label and label not in labels_done:
            labels_done.append(label)
        else:
            label = None
        pylab.barh(rindex, end - start, left=start, height=0.8, color=color, label=label)
    pylab.title(cluster_name)
    pylab.xlabel("Coordinates")
    pylab.ylabel("Reads")
    pylab.legend(loc='upper right')
    pylab.yticks(())
    out_file = "%s.png" % (cluster_name)
    pylab.savefig(out_file)

if __name__ == '__main__':
    
    try:
        alignment_file = sys.argv[1]
        cluster_name = sys.argv[2]
        cluster_distance = int(sys.argv[3])
    except:
        print __doc__
        sys.exit(-1)

    sys.exit(-1)
    cluster_distance = 10  # change according to the location you wish to look
    cluster_trees = collections.defaultdict(lambda:ClusterTree(cluster_distance, 2))
    # Parse alignment file 
    read_id_map, cnt, location_db, freq_db = dict(), 0, dict(), dict()
    align_generator = alignment_parse(alignment_file, cluster_name)
    for read_id, match_id, start, end in align_generator:
        if not read_id in read_id_map: # get rid of long read id with small numbers 
            cnt +=1 
            read_id_map[read_id] = cnt
        cluster_trees[match_id].insert(start, end, read_id_map[read_id])
        location_db[read_id_map[read_id]] = start, end  # location of reads in genome 
        if read_id_map[read_id] in freq_db: # occurence of each read 
            freq_db[read_id_map[read_id]] += 1
        else:
            freq_db[read_id_map[read_id]] = 1
    print 
    print 'Number of clusters: ' + str(len(cluster_trees))

    for chrom, cluster_tree in cluster_trees.items():
        for start, end, read_ids in cluster_tree.getregions():
            if len(read_ids) < 100:continue # just remove the cluster having less number of reads.
            print 
            print 'number of reads ' + str(len(read_ids))
            plot_cluster(cluster_name, read_ids, location_db, freq_db)
            break # only one at this time 
    print 
