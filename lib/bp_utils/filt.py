#! /usr/bin/env python

import sys, gzip, math, numpy
import pysam
from scipy import stats

def filter_by_control(tumor_bp_file, output_file, matched_control_bp_file, merged_control_file,
                      min_support_num, min_median_mapq, min_max_clip_size, max_control_num_thres,
                      ignore_seq_consist, permissible_range):

    use_matched_control = True if matched_control_bp_file != "" else False
    if use_matched_control: matched_control_db = pysam.TabixFile(matched_control_bp_file)

    use_merged_control = True if merged_control_file != "" else False
    if use_merged_control: merged_control_db = pysam.TabixFile(merged_control_file)

    hout = open(output_file, 'w')
    with gzip.open(tumor_bp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mapqs = [int(x) for x in F[6].split(';')]
            clip_sizes = [int(x) for x in F[7].split(';')]

            if len(mapqs) < min_support_num: continue
            if numpy.median(mapqs) < min_median_mapq: continue
            if max(clip_sizes) < min_max_clip_size: continue

            # filtering using merged control file
            merged_control_filt_flag = False 
            if use_merged_control:
                tabixErrorFlag = 0
                try:
                    # records = merged_control_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                    start_pos = int(F[1]) - 1 - permissible_range
                    if start_pos < 0: start_pos = 0
                    records = merged_control_db.fetch(F[0], start_pos , int(F[1]) + 1 + permissible_range)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        # if record[0] == F[0] and record[1] == F[1] and record[2] == F[2]:
                        if record[0] == F[0] and (int(F[1]) - permissible_range) <= int(record[1]) and int(record[1]) <= (int(F[1]) + permissible_range) and record[3] == F[3]:
                            if ignore_seq_consist or record[4] == F[4]:
                                merged_control_filt_flag = True

            if merged_control_filt_flag: continue

            # get readnum from matched control file
            if use_matched_control:
                num_matched_control = 0
                tabixErrorFlag = 0
                try:
                    # records = matched_control_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                    start_pos = int(F[1]) - 1 - permissible_range
                    if start_pos < 0: start_pos = 0
                    records = matched_control_db.fetch(F[0], start_pos, int(F[1]) + 1 + permissible_range)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1
            
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        # if record[0] == F[0] and record[1] == F[1] and record[2] == F[2]:
                        if record[0] == F[0] and (int(F[1]) - permissible_range) <= int(record[1]) and int(record[1]) <= (int(F[1]) + permissible_range) and record[3] == F[3]:
                            if ignore_seq_consist or record[4] == F[4]:
                                num_matched_control = len(record[6].split(';'))
            
            else:
                num_matched_control = "---"

            if use_matched_control and num_matched_control > max_control_num_thres: continue

            print >> hout, F[0] + '\t' + str(int(F[1])+1) + '\t' + F[3] + '\t' + F[4] + '\t' + str(len(mapqs)) + '\t' + str(num_matched_control)


    hout.close()


def filter_by_allele_freq(input_file, output_file, tumor_bam, matched_control_bam, tumor_AF_thres, control_AF_thres, max_fisher_pvalue):

    hout = open(output_file, 'w')

    print >> hout, '\t'.join(["Chr", "Pos", "Dir", "Junc_Seq", 
                              "Num_Tumor_Total_Read", "Num_Tumor_Var_Read", "Num_Control_Total_Read", "Num_Control_Var_Read",
                              "Minus_Log_Fisher_P_value"])

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tumor_num = int(F[4])
            control_num = int(F[5])
            region = F[0] + ':' + F[1] + '-' + F[1]

            depth_tumor_info = pysam.depth(tumor_bam, "-r", region)
            #depth_tumor = int(depth_tumor_info[0].rstrip('\n').split('\t')[2])
            depth_tumor = int(depth_tumor_info.rstrip('\n').split('\t')[2])
            # depth_tumor = int(depth_tumor_info.split('\t')[2])
            # depth_tumor = int(depth_tumor_info[0].split('\t')[2])
            AF_tumor = float(tumor_num) / depth_tumor
            if AF_tumor < tumor_AF_thres: continue

            # print '\t'.join(F)
            if matched_control_bam != "":
                depth_control_info = pysam.depth(matched_control_bam, "-r", region)
                #depth_control = int(depth_control_info[0].rstrip('\n').split('\t')[2]) if len(depth_control_info) != 0 else 0
                depth_control = int(depth_control_info.rstrip('\n').split('\t')[2]) if len(depth_control_info) != 0 else 0
                # depth_control = int(depth_control_info.split('\t')[2]) if len(depth_control_info) != 0 else 0
                # depth_control = int(depth_control_info[0].split('\t')[2]) if len(depth_control_info) != 0 else 0
                control_AF = float(control_num) / depth_control if depth_control > 0 else 1.0

            else:
                depth_control = "---"
                control_AF = "---"

            if control_AF != "---" and control_AF > control_AF_thres: continue 

            lpvalue = "---"
            if control_AF != "":
                oddsratio, pvalue = stats.fisher_exact([[depth_tumor - tumor_num, tumor_num], [depth_control - control_num, control_num]], 'less')
                if pvalue < 1e-100: pvalue = 1e-100
                lpvalue = (- math.log(pvalue, 10) if pvalue < 1 else 0)
                lpvalue = round(lpvalue, 4) 

                if 10**(-lpvalue) > float(max_fisher_pvalue): continue

            print >> hout, '\t'.join(F[:4]) + '\t' + str(depth_tumor) + '\t' + str(tumor_num) + '\t' + \
                           str(depth_control) + '\t' + str(control_num) + '\t' + str(lpvalue)

    hout.close()



