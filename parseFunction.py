#!/usr/bin/env python

"""
    functions for parsing breakpoint containing read pairs and improperly aligned read pairs
"""

import sys, pysam, re, subprocess, collections
# import utils

def parseJunctionFromBam(inputBAM, outputFilePath, min_mapping_qual, min_major_clip_size, max_minor_clip_size):

    """
    This function utilizes the SA tags (SA:Z:rname, pos, strand, CIGAR, mapQ, number of mismatch).
    The strand of the supplementary alignment information in the SA tag is determined by the orignal sequence (before taking complement).
    Therefore, please not that when the primary alignment is in the reverse direction, the sequence shown in the bam file does not match
    to the SA tags..
    """

    bamfile = pysam.Samfile(inputBAM, "rb")
    hOUT = open(outputFilePath, "w")
 
    SAre = re.compile('([^ \t\n\r\f\v,]+),(\d+),([\-\+]),(\w+),(\d+),(\d+);')
    cigarMDRe = re.compile('(\d+)([MD])')
    cigarHIMSRe = re.compile('(\d+)([HIMS])')
    cigarHSRe_right = re.compile('(\d+)([HS])$')
    cigarHSRe_left = re.compile('^(\d+)([HS])')


    # maybe add the regional extraction of bam files
    for read in bamfile.fetch():

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # no clipping
        if len(read.cigar) == 1: continue

        # skip if below the minimum mapping quality
        if (read.mapq < min_mapping_qual): continue

        # skip if the read aligned to hs37d5"
        # (in the future, this step will be replaced to some more sophisticated way;
        # (e.g., the user can input the target chromosomes and ignore if the read is aligned to non-target chromosomes, and so on..
        if bamfile.getrname(read.tid) == "hs37d5" or bamfile.getrname(read.rnext) == "hs37d5": continue


        # get the clipping size in the both side
        left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
        right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

        if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size: continue

        """
        # skip if there is no SA tags
        SA_str = None 
        for item in read.tags:
            if item[0] == "SA":
                SA_str = SAre.match(item[1])
                
        if SA_str is None: continue
        """

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")

        # get the soft clipping information on the supplementary alignment
        # right_clipping_SA = 0
        # tempMatch = cigarHSRe_right.search(cigar_SA)
        # if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        # left_clipping_SA = 0
        # tempMatch = cigarHSRe_left.search(cigar_SA)
        # if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))

        # skip if the both sides of the supplementary alignemt is clipped than the specified threshould
        # if left_clipping_SA > max_minor_clip_size and right_clipping_SA > max_minor_clip_size: continue


        # when the right side is clipped...
        if right_clipping >= min_major_clip_size:

            clipLen_current = right_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen

            juncChr_current = chr_current
            juncPos_current = pos_current + alignmentSize_current - 1
            juncDir_current = "+"
    
            juncseq_start = readLength_current - clipLen_current
            juncseq_end = readLength_current - clipLen_current + 8
            juncseq = read.seq[juncseq_start:juncseq_end]

            print >> hOUT, '\t'.join([juncChr_current, str(juncPos_current), juncDir_current, juncseq, read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq)])

        if left_clipping >= min_major_clip_size:

            clipLen_current = left_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
     
            juncChr_current = chr_current
            juncPos_current = pos_current
            juncDir_current = "-"

            juncseq_end = clipLen_current
            juncseq_start = clipLen_current - 8
            juncseq = read.seq[juncseq_start:juncseq_end]

            print >> hOUT, '\t'.join([juncChr_current, str(juncPos_current), juncDir_current, juncseq, read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq)])


    bamfile.close()
    hOUT.close()



def getPairStartPos(inputFilePath, outputFilePath):

    """
        script for obtaining the position information about the pair read from the junction file 

    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath + ".tmp", 'w')
    num = 1

    reChrPos = re.compile('^([^ \t\n\r\f\v,]+):(\d+)')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        ID = F[6]
        chr = ""
        pos = ""

        # obtain the information about the start site of the pair read 
        chrpos = reChrPos.search(F[12])
        if chrpos is not None:
            chr = chrpos.group(1)
            pos = chrpos.group(2)
        else:
            print '\t'.join(F)
            print "the 13th column did not match to (chr):(start) pattern"
            sys.exit()

        # change the pair read num
        if ID[-1:] == "1":
            ID = ID[:-1] + "2"
        else:
            ID = ID[:-1] + "1"

        print >> hOUT, chr + '\t' + str(int(pos) - 1) + '\t' + pos + '\t' + ID + '\t' + str(num)

        num = num + 1

    hIN.close()
    hOUT.close()


    hOUT = open(outputFilePath, 'w')
    subprocess.call(["sort", "-k1,1", "-k3,3n", outputFilePath + ".tmp"], stdout = hOUT)
    hOUT.close()


    ####################
    # delete intermediate file
    subprocess.call(["rm", outputFilePath + '.tmp'])



def cluster_breakpoint(input_file, output_file, check_interval):

    """
        script for merging and summarizing junction read pairs
    """
 
    hout = open(output_file, 'w')

    tmp_chr = ""
    tmp_pos = 0

    key2read = {}
    key2mapq = {}

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if tmp_chr != F[0] or int(F[1]) - tmp_pos > check_interval:
                del_list = []
                for key in key2read:
                    print >> hout, key + '\t' + ';'.join(key2read[key]) + '\t' + ';'.join(key2mapq[key])
                    del_list.append(key)

                for key in del_list:
                    del key2read[key]
                    del key2mapq[key]

                tmp_chr = F[0]  
                tmp_pos = int(F[1])

             
            key = F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]
            if key not in key2read: key2read[key] = []
            if key not in key2mapq: key2mapq[key] = []

            key2read[key].append(F[4])
            key2mapq[key].append(F[5])


    # final flush
    for key in key2read:
        print >> hout, key + '\t' + ';'.join(key2read[key]) + '\t' + ';'.join(key2mapq[key])


    hout.close()


