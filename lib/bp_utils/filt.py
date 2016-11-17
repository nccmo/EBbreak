#! /usr/bin/env python

import sys, gzip, math, re, subprocess, numpy
import pysam
from scipy import stats
import my_seq
import swalign

def filter_by_control(tumor_bp_file, output_file, matched_control_bp_file, merged_control_file,
                      min_support_num, min_median_mapq, min_max_clip_size, max_control_num_thres):

    use_matched_control = True if matched_control_bp_file != "" else False
    if use_matched_control: matched_control_db = pysam.TabixFile(matched_control_bp_file)

    use_merged_control = True if merged_control_file != "" else False
    if use_merged_control: merged_control_db = pysam.TabixFile(merged_control_file)

    hout = open(output_file, 'w')
    with gzip.open(tumor_bp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mapqs = [int(x) for x in F[5].split(';')]
            clip_sizes = [int(x) for x in F[6].split(';')]

            if len(mapqs) < min_support_num: continue
            if numpy.median(mapqs) < min_median_mapq: continue
            if max(clip_sizes) < min_max_clip_size: continue

            # filtering using merged control file
            merged_control_filt_flag = False 
            if use_merged_control:
                tabixErrorFlag = 0
                try:
                    records = merged_control_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if record[0] == F[0] and record[1] == F[1] and record[2] == F[2] and record[3] == F[3]:
                            merged_control_filt_flag = True

            if merged_control_filt_flag: continue

            # get readnum from matched control file
            if use_matched_control:
                num_matched_control = 0
                tabixErrorFlag = 0
                try:
                    records = matched_control_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1
            
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if record[0] == F[0] and record[1] == F[1] and record[2] == F[2] and record[3] == F[3]:
                            num_matched_control = len(record[5].split(';'))
            
            else:
                num_matched_control = "---"

            if use_matched_control and num_matched_control > max_control_num_thres: continue

            print >> hout, '\t'.join(F[:4]) + '\t' + str(len(mapqs)) + '\t' + str(num_matched_control)


    hout.close()


def filter_by_allele_freq(input_file, output_file, tumor_bam, matched_control_bam, tumor_AF_thres, control_AF_thres, max_fisher_pvalue):

    hout = open(output_file, 'w')
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

            # print '\t'.join(F)
            if matched_control_bam != "":
                depth_control_info = pysam.depth(matched_control_bam, "-r", region)
                depth_control = int(depth_control_info.rstrip('\n').split('\t')[2]) if depth_control_info != "" else 0
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

            print >> hout, '\t'.join(F) + '\t' + str(depth_tumor) + '\t' + str(depth_control) + '\t' + str(lpvalue)

    hout.close()



def assemble_seq(readid2seq, junc_seq, tmp_file_path):

    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print >> hout, '>' + tid
        print >> hout, readid2seq[tid]
    hout.close()
    
    hout = open(tmp_file_path + ".tmp3.assemble_output.fq", 'w')
    subprocess.call(["fml-asm", tmp_file_path + ".tmp3.assemble_input.fa"], stdout = hout) 
    hout.close()

    line_num = 0
    temp_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_output.fq", 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')

                aln_1 = sw.align(tseq, junc_seq)
                if aln_1.score >= 35:
                    ttcontig = tseq[aln_1.r_end:]
                    if len(ttcontig) > len(temp_contig): temp_contig = ttcontig
                
                aln_2 = sw.align(tseq, my_seq.reverse_complement(junc_seq))
                if aln_2.score >= 35:
                    ttcontig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                    if len(ttcontig) > len(temp_contig): temp_contig = ttcontig

    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.fa"])
    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_output.fq"])
    return temp_contig

 
def generate_contig(input_file, output_file, tumor_bp_file, tumor_bam, reference_genome):

    """
    tumor_bp_db = pysam.TabixFile(tumor_bp_file)

    readid2key = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')               
            tabixErrorFlag = 0
            try:
                records = tumor_bp_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1

            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    if record[0] == F[0] and record[1] == F[1] and record[2] == F[2] and record[3] == F[3]:
                        for readid in record[4].split(';'):
                            readid2key[re.sub(r'/\d$', '', readid)] = ','.join(F[:4])

 
    bamfile = pysam.Samfile(tumor_bam, "rb")

    hout = open(output_file + ".tmp2.contig.unsorted", 'w')
    for read in bamfile.fetch():
       
        if read.qname in readid2key:
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            print >> hout, readid2key[read.qname] + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence

    hout.close()

    hout = open(output_file + ".tmp2.contig.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig.unsorted"], stdout = hout)
    hout.close()

    """

    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    key2contig = {}
    with open(output_file + ".tmp2.contig.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if temp_key != F[0]:
                if len(temp_id2seq) > 0:
                    key2contig[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, output_file)

                temp_key = F[0]
                temp_id2seq = {}
                FF = temp_key.split(',')
                if FF[2] == "+":
                    temp_junc_seq = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - 20, int(FF[1]))
                else:
                    temp_junc_seq = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + 20))

            temp_id2seq[F[1]] = F[2]

        if len(temp_id2seq) > 0: 
            key2contig[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, output_file)


    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            key = ','.join(F[:4])

            if key not in key2contig: continue
            contig = key2contig[key]
            if len(contig) < 50: continue
            if contig[:8] != F[3][:8]: continue

            
            print >> hout, '\t'.join(F) + '\t' + contig

    hout.close()
