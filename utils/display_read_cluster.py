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
import pysam 
import collections
from bx.intervals.cluster import ClusterTree


def alignment_parse(bam_file, region_sltd):
    """
    Extract contents from BAM file
    """
    
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")

    chr_name, start_s, stop_s = re.search(r'^(.+):(\d+)-(\d+)$', region_sltd).group(1,2,3)

    for rec in sam_reader.fetch(chr_name, int(start_s), int(stop_s)):
        yield rec.qname, sam_reader.getrname(rec.rname), rec.pos, rec.pos+rec.qlen

    sam_reader.close()

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

    # - Distance in basepairs for two reads to be in the same cluster;
    #   for instance 20 would group all reads with 20bp of each other
    # - Number of reads necessary for a group to be considered a cluster;
    #   2 returns all groups with 2 or more overlapping reads
    cluster_trees = collections.defaultdict(lambda:ClusterTree(cluster_distance, 2))

    read_id_map, cnt, location_db, freq_db = dict(), 0, dict(), dict()
    align_generator = alignment_parse(alignment_file, cluster_name)

    for read_id, match_id, start, end in align_generator:
        #print read_id, match_id, start, end 
        if not read_id in read_id_map: #make read id compact 
            cnt +=1 
            read_id_map[read_id] = cnt

        cluster_trees[match_id].insert(start, end, read_id_map[read_id])
        break 

    sys.exit(-1)

    # Parse alignment file 
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
