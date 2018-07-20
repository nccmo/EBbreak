#!/usr/bin/env python

import sys, pysam, re, subprocess
import my_seq
# import utils

SAre = re.compile('([^ \t\n\r\f\v,]+),(\d+),([\-\+]),(\w+),(\d+),(\d+);')
cigarMDRe = re.compile('(\d+)([MD])')
cigarHIMSRe = re.compile('(\d+)([HIMS])')
cigarHSRe_right = re.compile('(\d+)([HS])$')
cigarHSRe_left = re.compile('^(\d+)([HS])')

def parse_bp_from_bam(input_bam, output_file, key_seq_size, min_major_clip_size, max_minor_clip_size):

    """
    function for getting breakpoints from BAM file
    input:bam file
    output:output_file + ".bp.tmp.txt"
    breakpoint info (juncChr, juncPos-1, juncPos, Dir (right clipping:+, left clipping:-), Juncseq, ID + ("/1" or "/2"), MAPQ, Clipping size, AlignmentSize)
    """

    bamfile = pysam.Samfile(input_bam, "rb")
    hout = open(output_file, "w")
 
    # maybe add the regional extraction of bam files
    for read in bamfile.fetch():

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if not aligned
        if flags[2] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # no clipping
        if len(read.cigar) == 1: continue

        # get the clipping size in the both side
        left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
        right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

        if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size: continue

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")

        # when the right side is clipped...
        if right_clipping >= min_major_clip_size:

            clipLen_current = right_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen

            juncChr_current = chr_current
            juncPos_current = pos_current + alignmentSize_current - 1
            juncDir_current = "+"
    
            juncseq_start = readLength_current - clipLen_current
            juncseq_end = readLength_current - clipLen_current + key_seq_size 
            juncseq = read.seq[juncseq_start:juncseq_end]

            print >> hout, '\t'.join([juncChr_current, str(juncPos_current-1), str(juncPos_current), juncDir_current, juncseq, 
                                      read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(right_clipping), str(alignmentSize_current)])

        if left_clipping >= min_major_clip_size:

            clipLen_current = left_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
     
            juncChr_current = chr_current
            juncPos_current = pos_current
            juncDir_current = "-"

            juncseq_end = clipLen_current
            juncseq_start = clipLen_current - key_seq_size 
            juncseq = my_seq.reverse_complement(read.seq[juncseq_start:juncseq_end])

            print >> hout, '\t'.join([juncChr_current, str(juncPos_current-1), str(juncPos_current), juncDir_current, juncseq, 
                                      read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(left_clipping), str(alignmentSize_current)])


    bamfile.close()
    hout.close()



def cluster_breakpoint(input_file, output_file, check_interval):

    """
    function for clustering breakpoints
    input:output_file + ".bp.tmp.txt"               result of "parse_bp_from_bam" function 
    output:output_file + ".bp.clustered.tmp.txt"    breakpoint info (juncChr, juncPos-1, juncPos, Dir (right clipping:+, left clipping:-), Juncseq, IDs + ("/1" or "/2"), MAPQs, Clipping sizes, AlignmentSizes)
    """

    hout = open(output_file, 'w')

    tmp_chr = ""
    tmp_pos = 0

    key2read = {}
    key2mapq = {}
    key2clipsize = {}
    key2alnsize = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if tmp_chr != F[0] or int(F[1]) - tmp_pos > check_interval:
                del_list = []
                for key in key2read:
                    print >> hout, key + '\t' + ';'.join(key2read[key]) + '\t' + ';'.join(key2mapq[key]) + '\t' + \
                                   ';'.join(key2clipsize[key]) + '\t' + ';'.join(key2alnsize[key])
                    del_list.append(key)

                for key in del_list:
                    del key2read[key]
                    del key2mapq[key]
                    del key2clipsize[key]
                    del key2alnsize[key]                

                tmp_chr = F[0]  
                tmp_pos = int(F[1])

             
            key = F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3] + '\t' + F[4]
            if key not in key2read: key2read[key] = []
            if key not in key2mapq: key2mapq[key] = []
            if key not in key2clipsize: key2clipsize[key] = []
            if key not in key2alnsize: key2alnsize[key] = []

            key2read[key].append(F[5])
            key2mapq[key].append(F[6])
            key2clipsize[key].append(F[7])
            key2alnsize[key].append(F[8])

    # final flush
    for key in key2read:
        print >> hout, key + '\t' + ';'.join(key2read[key]) + '\t' + ';'.join(key2mapq[key]) + '\t' + \
                       ';'.join(key2clipsize[key]) + '\t' + ';'.join(key2alnsize[key])


    hout.close()


