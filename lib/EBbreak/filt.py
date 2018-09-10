#! /usr/bin/env python

import sys, gzip, math, numpy
import pysam
from scipy import stats
import my_seq

def filter_by_merged_control(tumor_bp_file, output_file, merged_control_file,
                             min_median_mapq, min_max_clip_size, permissible_range):

    """
    filtering with merged control
    """

    use_merged_control = True if merged_control_file != "" else False
    if use_merged_control: merged_control_db = pysam.TabixFile(merged_control_file)

    hout = open(output_file, 'w')
    with gzip.open(tumor_bp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mapqs = [int(x) for x in F[6].split(';')]
            clip_sizes = [int(x) for x in F[7].split(';')]

            # remove breakpoint if supporting read does not meet the criteria below
            if numpy.median(mapqs) < min_median_mapq: continue
            if max(clip_sizes) < min_max_clip_size: continue
            #if numpy.median(base_qualities) < 20: continue
            #if numpy.sort(base_qualities)[-2] < 30: continue

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
                            #if ignore_seq_consist or record[4] == F[4]:
                            merged_control_filt_flag = True

            if merged_control_filt_flag: continue

            print >> hout, F[0] + '\t' + str(int(F[1])+1) + '\t' + F[3] + '\t' + F[4]

    hout.close()



def filter_by_allele_freq(input_file, output_file, tumor_bam, matched_control_bam, tumor_AF_thres, control_AF_thres, max_fisher_pvalue):

    """
    filtering by allele frequency
    """

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
            depth_tumor = int(depth_tumor_info.rstrip('\n').split('\t')[2])
            AF_tumor = float(tumor_num) / depth_tumor
            if AF_tumor < tumor_AF_thres: continue

            if matched_control_bam != "":
                depth_control_info = pysam.depth(matched_control_bam, "-r", region)
                depth_control = int(depth_control_info.rstrip('\n').split('\t')[2]) if len(depth_control_info) != 0 else 0
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



def filter_by_base_quality(input_file, output_file, tumor_bam, min_support_num, permissible_range):

    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tumor_bamfile = pysam.Samfile(tumor_bam, "rb")
            key_baseq = []
            tumor_var_read = 0
            tabixErrorFlag = 0
            try:
                records = tumor_bamfile.fetch(F[0], max(int(F[1]) - 1, 0) , int(F[1]) + 1)
            
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1

            if tabixErrorFlag == 0:
                for read in records:
                    flags = format(int(read.flag), "#014b")[:1:-1]

                    # skip if not aligned
                    if flags[2] == "1": continue

                    # skip supplementary alignment
                    if flags[8] == "1" or flags[11] == "1": continue

                    # skip duplicated reads
                    if flags[10] == "1": continue

                    # no clipping
                    if len(read.cigar) == 1: continue

                    left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
                    right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)
                    #print right_clipping

                    if F[2] == "+":
                        if right_clipping < 2: continue

                        juncPos_current = int(read.pos + 1) + read.alen - 1

                        if int(juncPos_current) != int(F[1]): continue

                        juncseq_start = read.rlen - right_clipping
                        juncseq_end = min(juncseq_start + 8, read.rlen)
                        juncseq = read.seq[juncseq_start:juncseq_end]

                        #if numpy.mean(read.query_qualities[juncseq_start:juncseq_end]) < 10:
                        #    continue

                        if juncseq in F[3]: 
                            tumor_var_read += 1
                            key_baseq.append(str(numpy.mean(read.query_qualities[juncseq_start:juncseq_end])))

                    if F[2] == "-":
                        if left_clipping < 2: continue
                        
                        juncPos_current = int(read.pos + 1)

                        if int(juncPos_current) != int(F[1]): continue

                        juncseq_start = left_clipping
                        juncseq_end = max(left_clipping - 8, 0)
                        juncseq = my_seq.reverse_complement(read.seq[juncseq_start:juncseq_end])
                        
                        #if numpy.mean(read.query_qualities[juncseq_start:juncseq_end]) < 10:
                        #    continue
                        
                        if juncseq in F[3]: 
                            tumor_var_read += 1
                            key_baseq.append(str(numpy.mean(read.query_qualities[juncseq_end:juncseq_start])))


            if tumor_var_read < min_support_num: continue
            if numpy.median(list(map(float, key_baseq))) < 20: continue
            if numpy.sort(list(map(float, key_baseq)))[-2] < 30: continue

            print >> hout, F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]  + '\t' + str(tumor_var_read)

    hin.close()
    hout.close()



def filter_by_matched_control(input_file, output_file, matched_control_bp_file, normal_bam, max_control_num_thres, permissible_range):

    """
    filtering with matched control
    """

    use_matched_control = True if matched_control_bp_file != "" else False
    if use_matched_control: normal_bamfile = pysam.Samfile(normal_bam, "rb")


    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        for line in hin:
            
            F = line.rstrip('\n').split('\t')
            normal_var_read = 0

            if use_matched_control:

                tabixErrorFlag = 0
                try:
                    records = normal_bamfile.fetch(F[0], max(int(F[1]) - 1, 0) , int(F[1]) + 1)
                
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for read in records:
                            flags = format(int(read.flag), "#014b")[:1:-1]

                            # skip if not aligned
                            if flags[2] == "1": continue

                            # skip supplementary alignment
                            if flags[8] == "1" or flags[11] == "1": continue

                            # skip duplicated reads
                            if flags[10] == "1": continue

                            # no clipping
                            if len(read.cigar) == 1: continue

                            if F[2] == "+":
                                right_clipping2 = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)
                                
                                if right_clipping2 < 2: continue

                                juncPos_current2 = int(read.pos + 1) + read.alen - 1
                                if max(int(F[1]) - permissible_range, 0) <= int(juncPos_current2) and int(juncPos_current2) <= (int(F[1]) + permissible_range):

                                    #if juncseq2 in F[3]: 
                                    normal_var_read += 1

                            if F[2] == "-":
                                left_clipping2 = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)

                                if left_clipping2 < 2: continue
                                juncPos_current2 = int(read.pos + 1)
                                if max(int(F[1]) - permissible_range, 0) <= int(juncPos_current2) and int(juncPos_current2) <= (int(F[1]) + permissible_range):

                                   #if juncseq2 in F[3]: 
                                    normal_var_read += 1

            else:
                normal_var_read == "---"

            if use_matched_control and normal_var_read > max_control_num_thres: continue

            print >> hout,  F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]  + '\t' + F[4] + '\t' + str(normal_var_read)



