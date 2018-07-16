#! /usr/bin/env python

import sys, os, math, re, subprocess
import pysam
import my_seq
import swalign
import annot_utils.gene, annot_utils.exon


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

    cap3output = tmp_file_path + ".tmp3.assemble_input.fa.cap3.screen"
    cap3input = tmp_file_path + ".tmp3.assemble_input.fa"
    sret = subprocess.Popen("cap3 %s -i 25 -j 31 -o 16 -s 251 -z 1 -c 10 > %s" % (cap3input,cap3output), shell=True).wait()

    if sret != 0:
        print >> sys.stderr, "cap3 error, error code: " + str(sret)
        #sys.exit()
 
    line_num = 0
    temp_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_input.fa.cap.contigs", 'r') as hin:
        Allhin = hin.read()
        hin2 = Allhin.replace('\n', '').split(">")
        tseq = []
        for i in hin2:
            if i != "":
                seq =  re.sub(r'\d+','',i[6:])
                tseq.append(seq)

        if len(tseq) == 1:
            ttseq = tseq[0]
            aln_1 = sw.align(ttseq, junc_seq)
            if aln_1.score >= 35:
                ttcontig = ttseq[aln_1.r_end:]
                if len(ttcontig) > len(temp_contig): temp_contig = ttcontig

            aln_2 = sw.align(ttseq, my_seq.reverse_complement(junc_seq))
            if aln_2.score >= 35:
                ttcontig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                if len(ttcontig) > len(temp_contig): temp_contig = ttcontig

            return temp_contig

        else:
            return ""


def assemble_seq2(readid2seq, tmp_file_path):

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print >> hout, '>' + tid
        print >> hout, readid2seq[tid]
    hout.close()


    cap3output = tmp_file_path + ".tmp3.assemble_input.fa.cap3.screen"
    cap3input = tmp_file_path + ".tmp3.assemble_input.fa"
    sret = subprocess.Popen("cap3 %s -i 25 -j 31 -o 16 -s 251 -z 1 -c 10 > %s" % (cap3input,cap3output), shell=True).wait()


    if sret != 0:
        print >> sys.stderr, "cap3 error, error code: " + str(sret)
        #sys.exit()
 
    line_num = 0
    temp_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_input.fa.cap.contigs", 'r') as hin:
        Allhin = hin.read()
        hin2 = Allhin.replace('\n', '').split(">")
        tseq = []
        for i in hin2:
            if i != "":
                seq =  re.sub(r'\d+','',i[6:])
                tseq.append(seq)

        if len(tseq) == 1:
            ttseq = tseq[0]
            return ttseq

        else:
            return ""

 
def generate_contig(input_file, output_file, tumor_bp_file, tumor_bam, reference_genome, min_contig_length):

    tumor_bp_db = pysam.TabixFile(tumor_bp_file)

    #readid2key gets paired-reads if either of the pair contains breakpoint
    readid2key = {}
    #readidkey2 gets reads if the read contains breakpoint (pair read is not included)
    readid2key2 = {}
    readid2key3 = {}
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
                    if record[0] == F[0] and (int(record[1])+1) == int(F[1]) and record[3] == F[2] and record[4] == F[3]:
                        for readid in record[5].split(';'):
                            readid2key[re.sub(r'/\d$', '', readid)] = ','.join(F[:4])
                            readid2key2[readid] =  ','.join(F[:4])


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

    #get reads if the read contains breakpoint (pair read is not included)
    hout = open(output_file + ".tmp2.contig2.unsorted", 'w')
    for read in bamfile.fetch():

        flags = format(int(read.flag), "#014b")[:1:-1]

        ID = read.qname + ("/1" if flags[6] == "1" else "/2")

        if ID in readid2key2:
            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            print >> hout, readid2key2[ID] + '\t' + ID + '\t' + read.query_sequence

    hout.close()

    #get reads if partner read contains breakpoint
    hout = open(output_file + ".tmp2.contig3.unsorted", 'w')
    for read in bamfile.fetch():

        flags = format(int(read.flag), "#014b")[:1:-1]

        pair_ID = read.qname + ("/2" if flags[6] == "1" else "/1")

        if pair_ID in readid2key2:
            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            print >> hout, readid2key2[pair_ID] + '\t' + ("/1" if flags[6] == "1" else "/2") + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence


    hout.close()


    hout = open(output_file + ".tmp2.contig.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig.unsorted"], stdout = hout)
    hout.close()

    hout = open(output_file + ".tmp2.contig2.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig2.unsorted"], stdout = hout)
    hout.close()

    hout = open(output_file + ".tmp2.contig3.sorted", 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", output_file + ".tmp2.contig3.unsorted"], stdout = hout)
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

    temp_key2 = ""
    temp_id2seq2 = {}
    temp_junc_seq2 = ""
    key2contig2 = {}
    with open(output_file + ".tmp2.contig2.sorted") as hin2:
        for line in hin2:
            F = line.rstrip('\n').split('\t')
            if temp_key2 != F[0]:
                if len(temp_id2seq2) > 0:
                    key2contig2[temp_key2] = assemble_seq(temp_id2seq2, temp_junc_seq2, output_file)

                temp_key2 = F[0]
                temp_id2seq2 = {}
                FF = temp_key2.split(',')
                if FF[2] == "+":
                    temp_junc_seq2 = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - 20, int(FF[1]))
                else:
                    temp_junc_seq2 = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + 20))

            temp_id2seq2[F[1]] = F[2]

        if len(temp_id2seq2) > 0: 
            key2contig2[temp_key2] = assemble_seq(temp_id2seq2, temp_junc_seq2, output_file)

    #readid2key2[pair_ID] + '\t' + ("/1" if flags[6] == "1" else "/2") + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence
    temp_key3 = ""
    temp_id2seq3 = {}
    key2contig3 = {}

    temp_key4 = ""
    temp_id2seq4 = {}
    key2contig4 = {}

    with open(output_file + ".tmp2.contig3.sorted") as hin3:
        for line in hin3:
            F = line.rstrip('\n').split('\t')
            if F[1] == "/1":
                if temp_key3 != F[0]:
                    if len(temp_id2seq3) > 0:
                        key2contig3[temp_key3] = assemble_seq2(temp_id2seq3, output_file)
                    temp_key3 = F[0]
                    temp_id2seq3 = {}
                temp_id2seq3[F[2]] = F[3]
            if F[1] == "/2":
                if temp_key4 != F[0]:
                    if len(temp_id2seq4) > 0:
                        key2contig4[temp_key4] = assemble_seq2(temp_id2seq4, output_file)
                    temp_key4 = F[0]
                    temp_id2seq4 = {}
                temp_id2seq4[F[2]] = F[3]

        if len(temp_id2seq3) > 0: 
            key2contig3[temp_key3] = assemble_seq2(temp_id2seq3, output_file)

        if len(temp_id2seq4) > 0: 
            key2contig4[temp_key4] = assemble_seq2(temp_id2seq4, output_file)


    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            key = ','.join(F[:4])

            if key not in key2contig: continue
            long_contig = key2contig[key]
            if key not in key2contig2: continue
            short_contig = key2contig2[key]
            if key not in key2contig3: continue
            pair_contig1 = key2contig3[key]
            if key not in key2contig4: continue
            pair_contig2 = key2contig4[key]

            if len(short_contig) > len(long_contig):
                contig = short_contig
            else:
                contig = long_contig
            if len(contig) < min_contig_length: continue
            # if contig[:8] != F[3][:8]: continue

            
            print >> hout, '\t'.join(F) + '\t' + contig + '\t' + long_contig + '\t' + short_contig + '\t' + pair_contig1 + '\t' + pair_contig2

    hout.close()


    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig2.sorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig2.unsorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig3.sorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig3.unsorted"])

    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.fa*"])




def psl_check(psl_file, key2seq, align_margin = 10000): 

    tempID = ""
    temp_align2score = {}
    key2align = {}
    key2best_score = {}
    key2margin = {}
    with open(psl_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue

            if tempID != F[9]:
                if tempID != "":
                    for k, v in sorted(temp_align2score.items(), key=lambda x: x[1], reverse = False):
                        key2align[tempID].append(k)
                        if len(key2align[tempID]) >= 10: break
                        if key2best_score[tempID] != float("inf") and key2margin[tempID] == float("inf"): # second key
                            key2margin[tempID] = temp_align2score[k] - key2best_score[tempID] 
                        if key2best_score[tempID] == float("inf"): # first key
                            key2best_score[tempID] = temp_align2score[k]

                tempID = F[9]
                temp_align2score = {}
                key2align[tempID] = []
                key2best_score[tempID] = float("inf") 
                key2margin[tempID] = float("inf")

            inseq = key2seq[tempID][0:int(F[11])]
            talign = ','.join([F[13], F[15], F[16], F[8], inseq, str(int(F[10]) - int(F[0]))])
            if int(F[10]) - int(F[0]) < align_margin:
                temp_align2score[talign] = int(F[10]) - int(F[0])

        for k, v in sorted(temp_align2score.items(), key=lambda x: x[1], reverse = False):
            key2align[tempID].append(k)
            if len(key2align[tempID]) >= 10: break
            if key2best_score[tempID] != float("inf") and key2margin[tempID] == float("inf"): # second key
                key2margin[tempID] = temp_align2score[k] - key2best_score[tempID]
            if key2best_score[tempID] == float("inf"): # first key
                key2best_score[tempID] = temp_align2score[k]

    return [key2align, key2best_score, key2margin]


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

    #partner1
    key2seq2 = {}
    hout = open(output_file + ".tmp4.contig.alignment_check2.fa", 'w')
    with open(contig_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            key2seq2[key] = F[12]
            print >> hout, '>' + key
            print >> hout, F[12]

    hout.close()

    #partner2
    key2seq3 = {}
    hout = open(output_file + ".tmp4.contig.alignment_check3.fa", 'w')
    with open(contig_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            key2seq3[key] = F[13]
            print >> hout, '>' + key
            print >> hout, F[13]

    hout.close()

    # reference genome ##############################
    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check.fa",
                           output_file + ".tmp4.contig.alignment_check.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human, key2bscore_human, key2margin_human = psl_check(output_file + ".tmp4.contig.alignment_check.psl", key2seq)


    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check2.fa",
                           output_file + ".tmp4.contig.alignment_check2.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human2, key2bscore_human2, key2margin_human2 = psl_check(output_file + ".tmp4.contig.alignment_check2.psl", key2seq)


    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check3.fa",
                           output_file + ".tmp4.contig.alignment_check3.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human3, key2bscore_human3, key2margin_human3 = psl_check(output_file + ".tmp4.contig.alignment_check3.psl", key2seq)
    ################################################


    # virus genome #################################
    if virus_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [virus_db, output_file + ".tmp4.contig.alignment_check.fa",
                               output_file + ".tmp4.contig.alignment_check_virus.psl"],
                               stdout = FNULL, stderr = subprocess.STDOUT)
    
        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()
   
        key2align_virus, key2bscore_virus, key2margin_virus = psl_check(output_file + ".tmp4.contig.alignment_check_virus.psl", key2seq)
    else:
        key2info_virus, key2bscore_virus, key2margin_virus = {}, {}, {}

    if virus_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [virus_db, output_file + ".tmp4.contig.alignment_check2.fa",
                               output_file + ".tmp4.contig.alignment_check_virus2.psl"],
                               stdout = FNULL, stderr = subprocess.STDOUT)
    
        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()
   
        key2align_virus2, key2bscore_virus2, key2margin_virus2 = psl_check(output_file + ".tmp4.contig.alignment_check_virus2.psl", key2seq)
    else:
        key2info_virus2, key2bscore_virus2, key2margin_virus2 = {}, {}, {}

    if virus_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [virus_db, output_file + ".tmp4.contig.alignment_check3.fa",
                               output_file + ".tmp4.contig.alignment_check_virus3.psl"],
                               stdout = FNULL, stderr = subprocess.STDOUT)
    
        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()
   
        key2align_virus3, key2bscore_virus3, key2margin_virus3 = psl_check(output_file + ".tmp4.contig.alignment_check_virus2.psl", key2seq)
    else:
        key2info_virus3, key2bscore_virus3, key2margin_virus3 = {}, {}, {}
    ################################################


    # repeat genome #################################
    if repeat_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [repeat_db, output_file + ".tmp4.contig.alignment_check.fa",
                        output_file + ".tmp4.contig.alignment_check_repeat.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_repeat, key2bscore_repeat, key2margin_repeat = psl_check(output_file + ".tmp4.contig.alignment_check_repeat.psl", key2seq)
    else:
        key2info_repeat, key2bscore_repeat, key2margin_repeat = {}, {}, {} 

    if repeat_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [repeat_db, output_file + ".tmp4.contig.alignment_check2.fa",
                        output_file + ".tmp4.contig.alignment_check_repeat2.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_repeat2, key2bscore_repeat2, key2margin_repeat2 = psl_check(output_file + ".tmp4.contig.alignment_check_repeat2.psl", key2seq)
    else:
        key2info_repeat2, key2bscore_repeat2, key2margin_repeat2 = {}, {}, {} 

    if repeat_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [repeat_db, output_file + ".tmp4.contig.alignment_check3.fa",
                        output_file + ".tmp4.contig.alignment_check_repeat3.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_repeat3, key2bscore_repeat3, key2margin_repeat3 = psl_check(output_file + ".tmp4.contig.alignment_check_repeat3.psl", key2seq)
    else:
        key2info_repeat3, key2bscore_repeat3, key2margin_repeat3 = {}, {}, {} 
    ################################################

    # mitochondria_genome #################################
    if mitochondria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [mitochondria_db, output_file + ".tmp4.contig.alignment_check.fa",
                        output_file + ".tmp4.contig.alignment_check_mitochondria.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_mitochondria, key2bscore_mitochondria, key2margin_mitochondria = psl_check(output_file + ".tmp4.contig.alignment_check_mitochondria.psl", key2seq)
    else:
        key2info_mitochondria, key2bscore_mitochondria, key2margin_mitochondria = {}, {}, {} 

    if mitochondria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [mitochondria_db, output_file + ".tmp4.contig.alignment_check2.fa",
                        output_file + ".tmp4.contig.alignment_check_mitochondria2.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_mitochondria2, key2bscore_mitochondria2, key2margin_mitochondria2 = psl_check(output_file + ".tmp4.contig.alignment_check_mitochondria2.psl", key2seq)
    else:
        key2info_mitochondria2, key2bscore_mitochondria2, key2margin_mitochondria2 = {}, {}, {} 

    if mitochondria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [mitochondria_db, output_file + ".tmp4.contig.alignment_check3.fa",
                        output_file + ".tmp4.contig.alignment_check_mitochondria3.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_mitochondria3, key2bscore_mitochondria3, key2margin_mitochondria3 = psl_check(output_file + ".tmp4.contig.alignment_check_mitochondria3.psl", key2seq)
    else:
        key2info_mitochondria3, key2bscore_mitochondria3, key2margin_mitochondria3 = {}, {}, {} 
    ################################################



    # bacteria_genome #################################
    if bacteria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [bacteria_db, output_file + ".tmp4.contig.alignment_check.fa",
                        output_file + ".tmp4.contig.alignment_check_bacteria.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_bacteria, key2bscore_bacteria, key2margin_bacteria = psl_check(output_file + ".tmp4.contig.alignment_check_bacteria.psl", key2seq)
    else:
        key2info_bacteria, key2bscore_bacteria, key2margin_bacteria = {}, {}, {} 

    if bacteria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [bacteria_db, output_file + ".tmp4.contig.alignment_check2.fa",
                        output_file + ".tmp4.contig.alignment_check_bacteria2.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_bacteria2, key2bscore_bacteria2, key2margin_bacteria2 = psl_check(output_file + ".tmp4.contig.alignment_check_bacteria2.psl", key2seq)
    else:
        key2info_bacteria2, key2bscore_bacteria2, key2margin_bacteria2 = {}, {}, {} 

    if bacteria_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [bacteria_db, output_file + ".tmp4.contig.alignment_check3.fa",
                        output_file + ".tmp4.contig.alignment_check_bacteria3.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_bacteria3, key2bscore_bacteria3, key2margin_bacteria3 = psl_check(output_file + ".tmp4.contig.alignment_check_bacteria3.psl", key2seq)
    else:
        key2info_bacteria3, key2bscore_bacteria3, key2margin_bacteria3 = {}, {}, {} 
    ################################################


    hout = open(output_file, 'w')    
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n')
        print >> hout, header + '\t' + '\t'.join(["Contig", "Junc_Seq_Consistency", "Human_Alignment", "Human_Mismatch", "Human_Margin",
                                                  "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                                  "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                                  "Bacteria_Alignment", "Bacteria_Mismatch", "Bacteria_Margin",
                                                  "Contig_pair1", "Human_Alignment_pair1", "Human_Mismatch_pair1", "Human_Margin_pair1",
                                                  "Virus_Alignment_pair1", "Virus_Mismatch_pair1", "Virus_Margin_pair1",
                                                  "Repeat_Alignment_pair1", "Repeat_Mismatch_pair1", "Repeat_Margin_pair1",
                                                  "Mitochondria_Alignment_pair1", "Mitochondria_Mismatch_pair1", "Mitochondria_Margin_pair1",
                                                  "Bacteria_Alignment_pair1", "Bacteria_Mismatch_pair1", "Bacteria_Margin_pair1",
                                                  "Contig_pair2", "Human_Alignment_pair2", "Human_Mismatch_pair2", "Human_Margin_pair2",
                                                  "Virus_Alignment_pair2", "Virus_Mismatch_pair2", "Virus_Margin_pair2",
                                                  "Repeat_Alignment_pair2", "Repeat_Mismatch_pair2", "Repeat_Margin_pair2",
                                                  "Mitochondria_Alignment_pair2", "Mitochondria_Mismatch_pair2", "Mitochondria_Margin_pair2",
                                                  "Bacteria_Alignment_pair2", "Bacteria_Mismatch_pair2", "Bacteria_Margin_pair2",
                                                  ])

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])

            seq = key2seq[key] if key in key2seq and len(key2seq[key]) > 0 else "---"
            pair_seq1 = key2seq2[key] if key in key2seq2 and len(key2seq2[key]) > 0 else "---"
            pair_seq2 = key2seq3[key] if key in key2seq3 and len(key2seq3[key]) > 0 else "---"

            junc_seq_consistency = "TRUE" if seq[:8] == F[3][:8] else "FALSE"

            align_human = ';'.join(key2align_human[key]) if key in key2align_human and len(key2align_human[key]) > 0 else "---"
            bscore_human = str(key2bscore_human[key]) if key in key2bscore_human and key2bscore_human[key] != float("inf") else "---"
            margin_human = str(key2margin_human[key]) if key in key2margin_human and key2margin_human[key] != float("inf") else "---"
 
            align_virus = ';'.join(key2align_virus[key]) if key in key2align_virus and len(key2align_virus[key]) > 0 else "---"
            bscore_virus = str(key2bscore_virus[key]) if key in key2bscore_virus and key2bscore_virus[key] != float("inf") else "---"
            margin_virus = str(key2margin_virus[key]) if key in key2margin_virus and key2margin_virus[key] != float("inf") else "---"

            align_repeat = ';'.join(key2align_repeat[key]) if key in key2align_repeat and len(key2align_repeat[key]) > 0 else "---"
            bscore_repeat = str(key2bscore_repeat[key]) if key in key2bscore_repeat and key2bscore_repeat[key] != float("inf") else "---"
            margin_repeat = str(key2margin_repeat[key]) if key in key2margin_repeat and key2margin_repeat[key] != float("inf") else "---"

            align_mitochondria = ';'.join(key2align_mitochondria[key]) if key in key2align_mitochondria and len(key2align_mitochondria[key]) > 0 else "---"
            bscore_mitochondria = str(key2bscore_mitochondria[key]) if key in key2bscore_mitochondria and key2bscore_mitochondria[key] != float("inf") else "---"
            margin_mitochondria = str(key2margin_mitochondria[key]) if key in key2margin_mitochondria and key2margin_mitochondria[key] != float("inf") else "---"

            align_bacteria = ';'.join(key2align_bacteria[key]) if key in key2align_bacteria and len(key2align_bacteria[key]) > 0 else "---"
            bscore_bacteria = str(key2bscore_bacteria[key]) if key in key2bscore_bacteria and key2bscore_bacteria[key] != float("inf") else "---"
            margin_bacteria = str(key2margin_bacteria[key]) if key in key2margin_bacteria and key2margin_bacteria[key] != float("inf") else "---"

            align_human2 = ';'.join(key2align_human2[key]) if key in key2align_human2 and len(key2align_human2[key]) > 0 else "---"
            bscore_human2 = str(key2bscore_human2[key]) if key in key2bscore_human2 and key2bscore_human2[key] != float("inf") else "---"
            margin_human2 = str(key2margin_human2[key]) if key in key2margin_human2 and key2margin_human2[key] != float("inf") else "---"
 
            align_virus2 = ';'.join(key2align_virus2[key]) if key in key2align_virus2 and len(key2align_virus2[key]) > 0 else "---"
            bscore_virus2 = str(key2bscore_virus2[key]) if key in key2bscore_virus2 and key2bscore_virus2[key] != float("inf") else "---"
            margin_virus2 = str(key2margin_virus2[key]) if key in key2margin_virus2 and key2margin_virus2[key] != float("inf") else "---"

            align_repeat2 = ';'.join(key2align_repeat2[key]) if key in key2align_repeat2 and len(key2align_repeat2[key]) > 0 else "---"
            bscore_repeat2 = str(key2bscore_repeat2[key]) if key in key2bscore_repeat2 and key2bscore_repeat2[key] != float("inf") else "---"
            margin_repeat2 = str(key2margin_repeat2[key]) if key in key2margin_repeat2 and key2margin_repeat2[key] != float("inf") else "---"

            align_mitochondria2 = ';'.join(key2align_mitochondria2[key]) if key in key2align_mitochondria2 and len(key2align_mitochondria2[key]) > 0 else "---"
            bscore_mitochondria2 = str(key2bscore_mitochondria2[key]) if key in key2bscore_mitochondria2 and key2bscore_mitochondria2[key] != float("inf") else "---"
            margin_mitochondria2 = str(key2margin_mitochondria2[key]) if key in key2margin_mitochondria2 and key2margin_mitochondria2[key] != float("inf") else "---"

            align_bacteria2 = ';'.join(key2align_bacteria2[key]) if key in key2align_bacteria2 and len(key2align_bacteria2[key]) > 0 else "---"
            bscore_bacteria2 = str(key2bscore_bacteria2[key]) if key in key2bscore_bacteria2 and key2bscore_bacteria2[key] != float("inf") else "---"
            margin_bacteria2 = str(key2margin_bacteria2[key]) if key in key2margin_bacteria2 and key2margin_bacteria2[key] != float("inf") else "---"

            align_human3 = ';'.join(key2align_human3[key]) if key in key2align_human3 and len(key2align_human3[key]) > 0 else "---"
            bscore_human3 = str(key2bscore_human3[key]) if key in key2bscore_human3 and key2bscore_human3[key] != float("inf") else "---"
            margin_human3 = str(key2margin_human3[key]) if key in key2margin_human3 and key2margin_human3[key] != float("inf") else "---"
 
            align_virus3 = ';'.join(key2align_virus3[key]) if key in key2align_virus3 and len(key2align_virus3[key]) > 0 else "---"
            bscore_virus3 = str(key2bscore_virus3[key]) if key in key2bscore_virus3 and key2bscore_virus3[key] != float("inf") else "---"
            margin_virus3 = str(key2margin_virus3[key]) if key in key2margin_virus3 and key2margin_virus3[key] != float("inf") else "---"

            align_repeat3 = ';'.join(key2align_repeat3[key]) if key in key2align_repeat3 and len(key2align_repeat3[key]) > 0 else "---"
            bscore_repeat3 = str(key2bscore_repeat3[key]) if key in key2bscore_repeat3 and key2bscore_repeat3[key] != float("inf") else "---"
            margin_repeat3 = str(key2margin_repeat3[key]) if key in key2margin_repeat3 and key2margin_repeat3[key] != float("inf") else "---"

            align_mitochondria3 = ';'.join(key2align_mitochondria3[key]) if key in key2align_mitochondria3 and len(key2align_mitochondria3[key]) > 0 else "---"
            bscore_mitochondria3 = str(key2bscore_mitochondria3[key]) if key in key2bscore_mitochondria3 and key2bscore_mitochondria3[key] != float("inf") else "---"
            margin_mitochondria3 = str(key2margin_mitochondria3[key]) if key in key2margin_mitochondria3 and key2margin_mitochondria3[key] != float("inf") else "---"

            align_bacteria3 = ';'.join(key2align_bacteria3[key]) if key in key2align_bacteria3 and len(key2align_bacteria3[key]) > 0 else "---"
            bscore_bacteria3 = str(key2bscore_bacteria3[key]) if key in key2bscore_bacteria3 and key2bscore_bacteria3[key] != float("inf") else "---"
            margin_bacteria3 = str(key2margin_bacteria3[key]) if key in key2margin_bacteria3 and key2margin_bacteria3[key] != float("inf") else "---"


            print >> hout, '\t'.join(F) + '\t' + seq + '\t' + junc_seq_consistency + '\t' + align_human + '\t' + bscore_human + '\t' + margin_human + '\t' + \
                           align_virus + '\t' + bscore_virus + '\t' + margin_virus + '\t' + align_repeat + '\t' + bscore_repeat + '\t' + margin_repeat + '\t' + \
                           align_mitochondria + '\t' + bscore_mitochondria + '\t' + margin_mitochondria + '\t' + align_bacteria + '\t' + bscore_bacteria + '\t' + margin_bacteria + '\t' + \
                           pair_seq1 + '\t' + align_human2 + '\t' + bscore_human2 + '\t' + margin_human2 + '\t' + \
                           align_virus2 + '\t' + bscore_virus2 + '\t' + margin_virus2 + '\t' + align_repeat2 + '\t' + bscore_repeat2 + '\t' + margin_repeat2 + '\t' + \
                           align_mitochondria2 + '\t' + bscore_mitochondria2 + '\t' + margin_mitochondria2 + '\t' + align_bacteria2 + '\t' + bscore_bacteria2 + '\t' + margin_bacteria2 + '\t' + \
                           pair_seq2 + '\t' + align_human3 + '\t' + bscore_human3 + '\t' + margin_human3 + '\t' + \
                           align_virus3 + '\t' + bscore_virus3 + '\t' + margin_virus3 + '\t' + align_repeat3 + '\t' + bscore_repeat3 + '\t' + margin_repeat3 + '\t' + \
                           align_mitochondria3 + '\t' + bscore_mitochondria3 + '\t' + margin_mitochondria3 + '\t' + align_bacteria3 + '\t' + bscore_bacteria3 + '\t' + margin_bacteria3

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat.psl"])

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check2.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat2.psl"])

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check3.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat3.psl"])

    
def annotate_break_point(input_file, output_file, genome_id, is_grc):

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, False)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, False)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    hout = open(output_file, 'w')
    header2ind = {}
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        print >> hout, '\t'.join(["Chr", "Pos", "Dir", "Junc_Seq", "Gene", "Exon"] + header[4:])
        for line in hin:
            F = line.rstrip('\n').split('\t')

            ##########
            # check gene annotation
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch(F[header2ind["Chr"]], int(F[header2ind["Pos"]]) - 1, int(F[header2ind["Pos"]]) + 1)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1

            gene = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene.append(record[3])

            gene = list(set(gene))
            if len(gene) == 0: gene.append("---")  

            ##########
            # check gene annotation
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch(F[header2ind["Chr"]], int(F[header2ind["Pos"]]) - 1, int(F[header2ind["Pos"]]) + 1)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1
                
            exon = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon.append(record[3])
                    
            exon = list(set(exon))
            if len(exon) == 0: exon.append("---")

            print >> hout, '\t'.join([F[header2ind[x]] for x in ["Chr", "Pos", "Dir", "Junc_Seq"]]) + '\t' + \
                           ','.join(gene) + '\t' + ';'.join(exon) + '\t' + '\t'.join(F[(header2ind["Junc_Seq"] + 1):])

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])