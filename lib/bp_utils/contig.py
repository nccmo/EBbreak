#! /usr/bin/env python

import sys, os, math, re, subprocess
import pysam
import my_seq
import swalign


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
    sret = subprocess.call(["fml-asm", tmp_file_path + ".tmp3.assemble_input.fa"], stdout = hout) 
    hout.close()

    if sret != 0:
        print >> sys.stderr, "fml-asm error, error code: " + str(sret)
        sys.exit()
 
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

 
def generate_contig(input_file, output_file, tumor_bp_file, tumor_bam, reference_genome, min_contig_length):

    tumor_bp_db = pysam.TabixFile(tumor_bp_file)

    readid2key = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')               
            if F[0] == "Chr": continue

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
            if len(contig) < min_contig_length: continue
            if contig[:8] != F[3][:8]: continue

            
            print >> hout, '\t'.join(F) + '\t' + contig

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])



def psl_check(psl_file, key2seq, align_margin): 

    tempID = ""
    temp_align2score = {}
    key2info = {}
    with open(psl_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue

            if tempID != F[9]:
                if tempID != "":
                    for k, v in sorted(temp_align2score.items(), key=lambda x: x[1], reverse = False):
                        key2info[tempID].append(k)
                        if len(key2info[tempID]) >= 10: break

                tempID = F[9]
                temp_align2score = {}
                key2info[tempID] = []

            inseq = key2seq[tempID][0:int(F[11])]
            talign = ','.join([F[13], F[15], F[16], F[8], inseq, str(int(F[10]) - int(F[0]))])
            if int(F[10]) - int(F[0]) < align_margin:
                temp_align2score[talign] = int(F[10]) - int(F[0])

        for k, v in sorted(temp_align2score.items(), key=lambda x: x[1], reverse = False):
            key2info[tempID].append(k)
            if len(key2info[tempID]) >= 10: break

    return key2info


def alignment_contig(input_file, contig_file, output_file, reference_genome, blat_option, virus_db, repeat_db):

    
    blat_cmds = ("blat " + blat_option).split(' ')

    key2seq = {}
    hout = open(output_file + ".tmp4.contig.alignment_check.fa", 'w')
    with open(contig_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            key2seq[key] = F[9]
            print >> hout, '>' + key
            print >> hout, F[9]

    hout.close()

    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check.fa",
                           output_file + ".tmp4.contig.alignment_check.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2info_human = psl_check(output_file + ".tmp4.contig.alignment_check.psl", key2seq, 10)

    if virus_db != "":
    
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [virus_db, output_file + ".tmp4.contig.alignment_check.fa",
                               output_file + ".tmp4.contig.alignment_check_virus.psl"],
                               stdout = FNULL, stderr = subprocess.STDOUT)
    
        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()
   
        key2info_virus = psl_check(output_file + ".tmp4.contig.alignment_check_virus.psl", key2seq, 20)
    else:
        key2info_virus = {}

    if repeat_db != "":
    
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [repeat_db, output_file + ".tmp4.contig.alignment_check.fa",
                        output_file + ".tmp4.contig.alignment_check_repeat.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2info_repeat = psl_check(output_file + ".tmp4.contig.alignment_check_repeat.psl", key2seq, 20)
    else:
        key2info_repeat = {}


    hout = open(output_file, 'w')    
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n')
        print >> hout, header + '\t' + '\t'.join(["Contig", "Human_Alignment", "Virus_Alignment", "RepBase_Alignment"]) 
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])

            seq = key2seq[key] if key in key2seq else "---"
            human = ';'.join(key2info_human[key]) if key in key2info_human and len(key2info_human[key]) > 0 else "---"
            virus = ';'.join(key2info_virus[key]) if key in key2info_virus and len(key2info_virus[key]) > 0 else "---"
            repeat = ';'.join(key2info_repeat[key]) if key in key2info_repeat and len(key2info_repeat[key]) > 0 else "---"
            print >> hout, '\t'.join(F) + '\t' + seq + '\t' + human + '\t' + virus + '\t' + repeat

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat.psl"])


     
