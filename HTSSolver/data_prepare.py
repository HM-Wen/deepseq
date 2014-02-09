#!/usr/bin/env python
"""
Program to prepare the dataset for RNAi screens 

Requirements:
    Numpy: 
    xlutils: http://www.python-excel.org/ deals xls files 
    openpyxl: http://packages.python.org/openpyxl/ deals xlsx files 
"""
import os 
import re 
import sys
import csv 
try:
    import xlutils
    import openpyxl
except:
    pass 
import collections 
import numpy as NP
import scipy.io as SIO 

def _open_compx_file(fnam):
    """
    open a file and returns the file handler
    currently supporting csv, tsv
    """
    ## TODO give support for xls, xlsx
    try:
        if os.path.splitext(fnam)[1] in [".csv", ".tsv"]:
            FH = csv.reader(open(fnam, 'rb'))
    except Exception as error:
        sys.exit(error)

    return FH 

def _open_txt_file(sfile):
    """
    open a file and returns the generator 
    supports the txt format 
    """

    try:
        sfh = open(sfile, 'rU')
    except Exception as error:
        sys.exit(error)

    for line in sfh:
        line = line.strip('\n\r')
        if not line:
            continue
        yield line 

    sfh.close()

def get_file_contents(htsfile, plate_sz):
    """
    Parse the assay experiment generated file
    """
    htsfh = _open_compx_file(htsfile)
    
    ## depends on the well type this should change 
    hts_score = NP.zeros(plate_sz)
    filtered_val = 0 

    for rec in htsfh:
        if not rec:
            continue 
        #if re.search(r'^Results', rec[0]):
        if re.search(r'^Nexp', rec[0]):
            filtered_val = 1
            continue 
        if filtered_val:
            if not rec[0]:
                continue
            ## check for the well rows  
            if rec[0] in ['%c' % x for x in range(65, 65+plate_sz[0])]: 
                if rec[-1] == '':
                    rec = rec[:-1]
                ## find the index for array 
                ## TODO fix the file parsing 
                ## 13613 1213 1213 #file-1
                ## 67563 39438 59835
                ## A 13652 434 378 #file-2
                ## B 8754 245 674

                hts_score[ord(rec[0])-65] = [float(score) for score in rec[1:]]

                ## TODO check for the header line to remove may be a additional variable from use will help 
                if rec[0]=='P':
                    filtered_val = 0

    return hts_score

def get_plate_reference(plate_file):
    """
    Parse the plate ID informations
    example: cat plate_id.txt 
    Plate_# Assay_Plate_ID siRNA_Source_PlateID
    1 S119_HA_01 S1-HM-01   
    2 S119_HA_02 S1-HM-02   
    3 S119_HA_03 S1-HM-03 
    """

    plate_id = [] 
    fn_prefix = dict()

    pfh = _open_txt_file(plate_file)
    for pln in pfh:
        ## TODO regexp just match to a non-digit character 
        if re.search(r'^Plate', pln, re.IGNORECASE):
            continue 
        ## TODO split based on tab also, check for the required fields after splitting 
        pln = pln.split(' ')
        plate_id.append((pln[1], pln[2]))
        fn_prefix[pln[1]]=1 
    ## convert to np object TODO how to create a numpy array from generator  
    barcode_raw=NP.zeros((len(plate_id),2), dtype=NP.object)
    for i in range(len(plate_id)):
        barcode_raw[i,0]=plate_id[i][0]
        barcode_raw[i,1]=plate_id[i][1]
    return barcode_raw, fn_prefix 

def get_genes_info(smfname):
    """
    Parse siRNA annotation from the provided file
    """
    ## TODO expecting xlsx xls csv formatted gene reference list

    sfh = _open_txt_file(smfname)

    ## single line file count 
    num_lines = sum(1 for line in open(smfname))

    ## TODO predefined # of colums from annotation file 
    anno_desc = NP.zeros((num_lines, 9), dtype=NP.object)
    for idx, sln in enumerate(sfh):
        parts = sln.split(' ')
        
        ## don't need informations from the blank, negative and positive control wells   
        if len(parts)>9:
            anno_desc[idx,0] = parts[0]
            anno_desc[idx,1] = parts[1]
            anno_desc[idx,2] = parts[9]
            anno_desc[idx,3] = parts[8]
            anno_desc[idx,4] = parts[7]
            anno_desc[idx,5] = parts[9]
            anno_desc[idx,6] = parts[4]
            anno_desc[idx,7] = parts[6]
            anno_desc[idx,8] = parts[5]

    ## process the elements  
    ## predefined # of columns from reference annotation
    header = NP.zeros((1,9), dtype=NP.object)
    for ele in range(len(anno_desc[0])):
        header[0,ele] = anno_desc[0][ele]

    ## delete the header entry from the main array 
    anno_desc = NP.delete(anno_desc, 0, 0) 

    ## delete empty element row
    remove_idx = [] 
    for edx in range(len(anno_desc)):
        if not anno_desc[edx,:].any():
            remove_idx.append(edx) 
    anno_desc = NP.delete(anno_desc, remove_idx, 0)

    return header, anno_desc

def logical_plate_dist(plate_ctrl_ref, well_checks,plate_sz):
    """
    Make the mark on plates based on the relevance of the experiment control setup
    plate_ctrl_ref: a logical numpy array with false depends on the plate size  
    well_checks: selected wells value should change to True 
                 A the complete row
                 A1 first well (1,1)
                 1 first complete column 
                 A2,C2 row A and C of column 2
    """
    if well_checks == 'None':
        return plate_ctrl_ref
        
    well_checks = well_checks.split(',')

    alpha = ['%c' % x for x in range(65, 65+plate_sz[0])]
    for rc in well_checks:
        if len(rc) == 1 and rc.upper() in alpha:
            idx = ord(rc.upper())-65 # 0-based index  
            plate_ctrl_ref[idx] = [1 for j in range(plate_sz[1])] # true for entire row 
        elif len(rc) == 1 and int(rc) in [k+1 for k in range(0,plate_sz[1])]: 
            plate_ctrl_ref[:, int(rc)-1]=[1 for j in range(plate_sz[0])] # true for entire colomun 
        elif len(rc)>1:
            i,j = rc[0], rc[1:]
            if i.upper() in alpha and int(j) in [k+1 for k in range(0,plate_sz[1])]:
                idx = ord(i.upper())-65
                plate_ctrl_ref[idx][int(j)-1]=1 # true for specific well a12,c12 ...

    return plate_ctrl_ref

def __main__():

    ## replicates, plates and control and test well(s) input data handling in general way TODO 
    ## Usage python data_prepare.py plate_ID.txt CSV_file_path siRNA_reference plate_size  
    ##                              pos_cl_1 pos_cl_2 pos_cl_3 
    ##                              minus_1 minus_2 blank_1 
    try:
        plate_ref = sys.argv[1]
        assayfile_path = sys.argv[2]
        smRNA_fname = sys.argv[3]
        plate_size_nr = int(sys.argv[4]) # 0 for 8*12 1 for 16*24
        plus_ctrl_1 = sys.argv[5] # 3 C A2,C2 
        plus_ctrl_2 = sys.argv[6] # 3 C A2,C2 
        plus_ctrl_3 = sys.argv[7] # 3 C A2,C2 
        minus_ctrl_1 = sys.argv[8] # negative controls 
        minus_ctrl_2 = sys.argv[9] # negative reference controls
        blank_ctrl = sys.argv[10] # blank wells 

        #plate_stats_cl = int(sys.argv[11]) # Plate Quality Metrics calculations, accept values 0,1
        #log_option_cl = int(sys.argv[12]) # 1 for log2, 2 for log10, 3 for raw data in the QC and scoring calculations 
        #kn_cl = int(sys.argv[13]) ## 1 Automatically knockout outlier control points. accepts 0,1
        #cutof_outlier = int(sys.argv[14]) ## outlier threshold 1 for 0.5*stdev, 2 for 1*stdev, 3 for 2*stdev, 4 for 3*stdev, 5 for 4*stdev 
        #max_knout = int(sys.argv[15]) ## max % of points to be knocked from each control set in each plate if they are outside the selected outlier threshold 
        #delete_lines_nr = sys.argv[] 
    except:
        print __doc__
        sys.exit(-1)

    ## control plate selection validation 
    cross_checks = [] 
    for xk in [plus_ctrl_1, plus_ctrl_2, plus_ctrl_3, minus_ctrl_1, minus_ctrl_2, blank_ctrl]:
        for wl in xk.split(','):
            cross_checks.append(wl.upper())
    assert len(cross_checks) == len(set(cross_checks)), 'please select different control/blank reference sets, found common controls, cannot continue.' 

    ## 96=8*12 or 384=16*24 well 
    plate_size = [(8, 12), (16,24)]
    
    ## get plate records
    barcode, prefix = get_plate_reference(plate_ref)

    ## fetch the data files from data storage path 
    plate_config = [] 
    for files in os.listdir(os.path.abspath(assayfile_path)):
        ## check for validating true files from the repository 
        fl_check = False
        for prefix_term in prefix.keys():
            if files.startswith(prefix_term):
                fl_check = True 
                break
        if not fl_check:
            continue 
        ## plate type 
        exp_score = get_file_contents(os.path.abspath(assayfile_path)+'/'+files, plate_size[plate_size_nr])
        ## TODO plate normalization FL and control 

        plate_config.append((files, exp_score))

    plate_config.sort()
    ## check file processing order for double checking
    col_data = [xd[1] for xd in plate_config]
    ## 3d array 16X24X#plates
    plate_data = NP.dstack(col_data)
    
    ## get reference annotation of small RNA from the file
    ref_header, ref_source = get_genes_info(smRNA_fname)
    #print len(ref_header[0])

    ## plate control references : define the reference plate check according to the selection 
    pos_ctrl_1 = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    pos_ctrl_1 = logical_plate_dist(pos_ctrl_1, plus_ctrl_1, plate_size[plate_size_nr])

    pos_ctrl_2 = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    pos_ctrl_2 = logical_plate_dist(pos_ctrl_2, plus_ctrl_2, plate_size[plate_size_nr])

    pos_ctrl_3 = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    pos_ctrl_3 = logical_plate_dist(pos_ctrl_3, plus_ctrl_3, plate_size[plate_size_nr])

    neg_ctrl_1 = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    neg_ctrl_1 = logical_plate_dist(neg_ctrl_1, minus_ctrl_1, plate_size[plate_size_nr])

    neg_ctrl_2 = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    neg_ctrl_2 = logical_plate_dist(neg_ctrl_2, minus_ctrl_2, plate_size[plate_size_nr])

    blank_wl = NP.zeros(plate_size[plate_size_nr], dtype=NP.bool)
    blank_wl = logical_plate_dist(blank_wl, blank_ctrl, plate_size[plate_size_nr])

    ## TODO out path need to set 
    #output_file_name_cl = "customHTS_out.txt" 
    #output_destination_path_input_cl = "/Users/vipin/Desktop/" 
    ## get more input variables which includes the statistic method 
    mat_chart = dict(barcode_raw = barcode,
                    data = plate_data,
                    vendor_header = ref_header,
                    vendor_source = ref_source,
                    source_col_number_cl = len(ref_header[0]),
                    pos_controls_cl = pos_ctrl_1,
                    posNumber2_controls_cl = pos_ctrl_2,
                    posNumber3_controls_cl = pos_ctrl_3,
                    neg_controls_cl = neg_ctrl_1,
                    negref_controls_cl = neg_ctrl_2,
                    blanks_cl = blank_wl
                    #plate_stats_cl = plate_stats_cl,
                    #log_option_selection_cl = log_option_cl,
                    #knockout_points_cl = kn_cl, 
                    #outlier_threshold_cl = cutof_outlier,
                    #max_knockout_cl = max_knout,
                    #plate_size_cl = plate_size_nr+1;
                    #output_file_name_cl = output_file_name_cl,
                    #output_destination_path_input_cl = output_destination_path_input_cl
                    )

    SIO.savemat('hts_data.mat', 
                mat_chart,
                oned_as='row', 
                format='5')
    
    ## TODO:
    ## RL POL-2 promoter house keeping 
    ## WNT signal effects the FL expression, signal ... 
    ## how the test control works in both cells, This has to validate the presebnc
    ## 1 -ve and 2 positive controls which one make up and other make down 
    ## 23 well, 24 well -black list unwanted well for analysis also you will find those candidates from the inside
    ## normalize the controls to the individual wells. 
    ## Normalize the plate intensity between different plates  
    ## Individual well over -ve control, visualize the result based on up/down positive along with the control
    ## Normalize between plates 
    ## black listing the candidates
    ## then distribution between replicates, and look for hte candidates, which are interested to look for the further analysis
    ## the ... plate -> replicate -> 1,2,3

if __name__ == "__main__":
    __main__()
