#! /usr/bin/env python

import sys, os, math, re, subprocess
import pysam
import my_seq
import swalign
import annot_utils.gene, annot_utils.exon


def assemble_seq_fml_asm(readid2seq, junc_seq, tmp_file_path):

    match = 2
    mismatch = -1
    gap_penalty = -1
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
    temp_all_contig = ""
    temp_genome_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_output.fq", 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')
                
                aln_1 = sw.align(tseq, junc_seq)
                aln_2 = sw.align(tseq, my_seq.reverse_complement(junc_seq))
            
                if aln_1.score >= aln_2.score:
                    if aln_1.score > temp_score:
                        temp_contig = tseq[aln_1.r_end:]
                        temp_all_contig = str(tseq)
                        temp_genome_contig = tseq[:aln_1.r_end]
                        temp_score = aln_1.score

                    if aln_1.score == temp_score:
                        if len(ttseq[aln_1.r_end:]) > len(temp_contig): 
                            temp_contig = tseq[aln_1.r_end:]
                            temp_all_contig = str(tseq)
                            temp_genome_contig = tseq[:aln_1.r_end]
                            temp_score = aln_1.score

                if aln_2.score > aln_1.score:
                    if aln_2.score > temp_score:
                        temp_contig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                        temp_all_contig = my_seq.reverse_complement(tseq)
                        temp_genome_contig = my_seq.reverse_complement(tseq[aln_2.r_pos:])
                        temp_score = aln_2.score

                    if aln_2.score == temp_score:
                        if len(ttseq[aln_2.r_end:]) > len(temp_contig): 
                            temp_contig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                            temp_all_contig = my_seq.reverse_complement(tseq)
                            temp_genome_contig = my_seq.reverse_complement(tseq[aln_2.r_pos:])
                            temp_score = aln_2.score
    ###############################

    return temp_contig + '\t' + temp_all_contig + '\t' + temp_genome_contig + '\t' + str(temp_score)
    hin.close()



def assemble_seq_cap3(readid2seq, junc_seq, tmp_file_path, swalign_score, juncseq):

    match = 2
    mismatch = -1
    gap_penalty = -1
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
    temp_all_contig = ""
    temp_genome_contig = ""
    temp_score = 0
    with open(tmp_file_path + ".tmp3.assemble_input.fa.cap.contigs", 'r') as hin:
        # Allhin = hin.read()
        # hin2 = Allhin.replace('\n', '').split(">")
        # tseq = []
        # for i in hin2:
        #     if i != "":
        #         seq =  re.sub(r'\d+','',i[6:])
        #         tseq.append(seq)
        seq = str()
        for line in hin:
            if line[0] == ">":
                F = ('\t')
                seq = seq + F
            else:
                F = (line.strip('\n'))
                seq = seq + F
        tseq = seq.split('\t')
        tseq = filter(lambda str:str != '', tseq)
        
        for i in range(0, len(tseq)):
            ttseq = tseq[i]
            aln_1 = sw.align(ttseq, junc_seq)
            aln_2 = sw.align(ttseq, my_seq.reverse_complement(junc_seq))

#############################
            if aln_1.score >= aln_2.score:
                if aln_1.score > temp_score:
                    temp_contig = ttseq[aln_1.r_end:]
                    temp_all_contig = str(ttseq)
                    temp_genome_contig = ttseq[:aln_1.r_end]
                    temp_score = aln_1.score

                if aln_1.score == temp_score:
                    if len(ttseq[aln_1.r_end:]) > len(temp_contig): 
                        temp_contig = ttseq[aln_1.r_end:]
                        temp_all_contig = str(ttseq)
                        temp_genome_contig = ttseq[:aln_1.r_end]
                        temp_score = aln_1.score

            if aln_2.score > aln_1.score:
                if aln_2.score > temp_score:
                    temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                    temp_all_contig = my_seq.reverse_complement(ttseq)
                    temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                    temp_score = aln_2.score

                if aln_2.score == temp_score:
                    if len(ttseq[aln_2.r_end:]) > len(temp_contig): 
                        temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                        temp_all_contig = my_seq.reverse_complement(ttseq)
                        temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                        temp_score = aln_2.score
###############################

    return temp_contig + '\t' + temp_all_contig + '\t' + temp_genome_contig + '\t' + str(temp_score)
    hin.close()

    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.fa"])
    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.cap*"])


def assemble_seq_sga(readid2seq, junc_seq, tmp_file_path, swalign_score, juncseq):

    match = 2
    mismatch = -1
    gap_penalty = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print >> hout, '>' + tid
        print >> hout, readid2seq[tid]
    hout.close()
    
    os.chdir(os.path.dirname(tmp_file_path))

    sga_input = os.path.basename(tmp_file_path) + ".tmp3.assemble_input.fa"
    

    sret = subprocess.Popen("sga index %s" % (sga_input), shell=True).wait()

    if sret != 0:
        print >> sys.stderr, "sga error, error code: " + str(sret)
        #sys.exit()

    sret2 = subprocess.Popen("sga fm-merge %s" % (sga_input), shell=True).wait()

    if sret2 != 0:
        print >> sys.stderr, "sga error, error code: " + str(sret)
        #sys.exit()

    os.chdir('../')
    line_num = 0
    temp_contig = ""
    temp_all_contig = ""
    temp_genome_contig = ""
    temp_score = 0
    try:
        with open(tmp_file_path + ".tmp3.assemble_input.merged.fa", 'r') as hin:
            tseq = []
            for line in hin:
                F = line.rstrip('\n')
                if F[0] != ">":
                    tseq.append(F)

            for i in range(0, len(tseq)):
                ttseq = tseq[i]
                aln_1 = sw.align(ttseq, junc_seq)
                aln_2 = sw.align(ttseq, my_seq.reverse_complement(junc_seq))

                if aln_1.score >= aln_2.score:
                    if aln_1.score > temp_score:
                        temp_contig = ttseq[aln_1.r_end:]
                        temp_all_contig = str(ttseq)
                        temp_genome_contig = ttseq[:aln_1.r_end]
                        temp_score = aln_1.score

                    if aln_1.score == temp_score:
                        if len(ttseq[aln_1.r_end:]) > len(temp_contig): 
                            temp_contig = ttseq[aln_1.r_end:]
                            temp_all_contig = str(ttseq)
                            temp_genome_contig = ttseq[:aln_1.r_end]
                            temp_score = aln_1.score

                if aln_2.score > aln_1.score:
                    if aln_2.score > temp_score:
                        temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                        temp_all_contig = my_seq.reverse_complement(ttseq)
                        temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                        temp_score = aln_2.score

                    if aln_2.score == temp_score:
                        if len(ttseq[aln_2.r_end:]) > len(temp_contig): 
                            temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                            temp_all_contig = my_seq.reverse_complement(ttseq)
                            temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                            temp_score = aln_2.score

            return temp_contig + '\t' + temp_all_contig + '\t' + temp_genome_contig + '\t' + str(temp_score)

    except:
        return "" + '\t' + "" + '\t' + ""

    hin.close()

    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.merged.fa"])
    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.sai"])
    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.rsai"])
    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.bwt"])
    # subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.rbwt"])


def assemble_seq_velvet(readid2seq, junc_seq, tmp_file_path):

    match = 2
    mismatch = -1
    gap_penalty = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print >> hout, '>' + tid
        print >> hout, readid2seq[tid]
    hout.close()
    
    os.system("mkdir -p velvet_output_dir/")
    velvet_input = tmp_file_path + ".tmp3.assemble_input.fa"
    sret = subprocess.Popen("velveth velvet_output_dir/ 31 %s" % (velvet_input), shell=True).wait()

    if sret != 0:
        print >> sys.stderr, "velveth error, error code: " + str(sret)
        #sys.exit()

    if sret == 0:
        sret2 = subprocess.Popen("velvetg velvet_output_dir/", shell=True).wait()

        if sret2 != 0:
            print >> sys.stderr, "velvetg error, error code: " + str(sret)
            #sys.exit()

        line_num = 0
        temp_contig = ""
        temp_all_contig = ""
        temp_genome_contig = ""
        sys.path.append("velvet_output_dir/")
        with open("velvet_output_dir/contigs.fa", 'r') as hin:
            seq = str()
            for line in hin:
                if line[0] == ">":
                    F = ('\t')
                    seq = seq + F
                else:
                    F = (line.strip('\n'))
                    seq = seq + F
            tseq = seq.split('\t')
            tseq = filter(lambda str:str != '', tseq)

            for i in range(0, len(tseq)):
                ttseq = tseq[i]
                aln_1 = sw.align(ttseq, junc_seq)
                aln_2 = sw.align(ttseq, my_seq.reverse_complement(junc_seq))

#############################
                if aln_1.score >= aln_2.score:
                    if aln_1.score > temp_score:
                        temp_contig = ttseq[aln_1.r_end:]
                        temp_all_contig = str(ttseq)
                        temp_genome_contig = ttseq[:aln_1.r_end]
                        temp_score = aln_1.score

                    if aln_1.score == temp_score:
                        if len(ttseq[aln_1.r_end:]) > len(temp_contig): 
                            temp_contig = ttseq[aln_1.r_end:]
                            temp_all_contig = str(ttseq)
                            temp_genome_contig = ttseq[:aln_1.r_end]
                            temp_score = aln_1.score

                if aln_2.score > aln_1.score:
                    if aln_2.score > temp_score:
                        temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                        temp_all_contig = my_seq.reverse_complement(ttseq)
                        temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                        temp_score = aln_2.score

                    if aln_2.score == temp_score:
                        if len(ttseq[aln_2.r_end:]) > len(temp_contig): 
                            temp_contig = my_seq.reverse_complement(ttseq[:aln_2.r_pos])
                            temp_all_contig = my_seq.reverse_complement(ttseq)
                            temp_genome_contig = my_seq.reverse_complement(ttseq[aln_2.r_pos:])
                            temp_score = aln_2.score
###############################

        return temp_contig + '\t' + temp_all_contig + '\t' + temp_genome_contig + '\t' + str(temp_score)
        hin.close()
        os.system("rm velvet_output_dir/*")



def generate_contig(input_file, output_file, tumor_bp_file, tumor_bam, reference_genome, min_contig_length, swalign_length, swalign_score):

    """
    function for generating contigs
    """

    tumor_bp_db = pysam.TabixFile(tumor_bp_file)

    #readid2key gets paired-reads if either of the pair contains breakpoint
    readid2key = {}
    #readidkey2 gets reads if the read contains breakpoint (pair read is not included)
    readid2key2 = {}
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
                    if record[0] == F[0] and (int(record[1])+1) == int(F[1]) and record[3] == F[2]:
                    #if record[0] == F[0] and (int(record[1])+1) == int(F[1]) and record[3] == F[2] and record[4] == F[3]:
                        for readid in record[5].split(';'):
                            readid2key[re.sub(r'/\d$', '', readid)] = ','.join(F[:4])
                            readid2key2[readid] =  ','.join(F[:4])


    bamfile = pysam.Samfile(tumor_bam, "rb")

    #gets paired-reads if either of the pair contains breakpoint
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

    # #get reads if partner read contains breakpoint
    # hout = open(output_file + ".tmp2.contig3.unsorted", 'w')
    # for read in bamfile.fetch():

    #     flags = format(int(read.flag), "#014b")[:1:-1]

    #     pair_ID = read.qname + ("/2" if flags[6] == "1" else "/1")

    #     if pair_ID in readid2key2:
    #         # skip supplementary alignment
    #         if flags[8] == "1" or flags[11] == "1": continue

    #         # skip duplicated reads
    #         if flags[10] == "1": continue

    #         print >> hout, readid2key2[pair_ID] + '\t' + ("/1" if flags[6] == "1" else "/2") + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence


    # hout.close()

    #sort reads
    hout = open(output_file + ".tmp2.contig.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig.unsorted"], stdout = hout)
    hout.close()

    hout = open(output_file + ".tmp2.contig2.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig2.unsorted"], stdout = hout)
    hout.close()

    # hout = open(output_file + ".tmp2.contig3.sorted", 'w')
    # subprocess.call(["sort", "-k1,1", "-k2,2n", output_file + ".tmp2.contig3.unsorted"], stdout = hout)
    # hout.close()


    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    key2contig_cap3 = {}
    #key2contig_fml_asm = {}
    #key2contig_velvet = {}
    key2contig_sga = {}
    with open(output_file + ".tmp2.contig.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if temp_key != F[0]:
                if len(temp_id2seq) > 0:
                    key2contig_cap3[temp_key] = assemble_seq_cap3(temp_id2seq, temp_junc_seq, output_file, swalign_score, juncseq)
                    #key2contig_fml_asm[temp_key] = assemble_seq_fml_asm(temp_id2seq, temp_junc_seq, output_file)
                    #key2contig_velvet[temp_key] = assemble_seq_velvet(temp_id2seq, temp_junc_seq, output_file)
                    key2contig_sga[temp_key] = assemble_seq_sga(temp_id2seq, temp_junc_seq, output_file, swalign_score, juncseq)

                temp_key = F[0]
                temp_id2seq = {}
                FF = temp_key.split(',')
                if FF[2] == "+":
                    temp_junc_seq = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - swalign_length, int(FF[1]))
                else:
                    temp_junc_seq = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + swalign_length))
                juncseq = FF[3]

            temp_id2seq[F[1]] = F[2]

        if len(temp_id2seq) > 0: 
            key2contig_cap3[temp_key] = assemble_seq_cap3(temp_id2seq, temp_junc_seq, output_file, swalign_score, juncseq)
            #key2contig_fml_asm[temp_key] = assemble_seq_fml_asm(temp_id2seq, temp_junc_seq, output_file)
            #key2contig_velvet[temp_key] = assemble_seq_velvet(temp_id2seq, temp_junc_seq, output_file)
            key2contig_sga[temp_key] = assemble_seq_sga(temp_id2seq, temp_junc_seq, output_file, swalign_score, juncseq)

    temp_key2 = ""
    temp_id2seq2 = {}
    temp_junc_seq2 = ""
    key2contig2_cap3 = {}
    #key2contig2_fml_asm = {}
    #key2contig2_velvet = {}
    #key2contig2_sga = {}
    with open(output_file + ".tmp2.contig2.sorted") as hin2:
        for line in hin2:
            F = line.rstrip('\n').split('\t')
            if temp_key2 != F[0]:
                if len(temp_id2seq2) > 0:
                    key2contig2_cap3[temp_key2] = assemble_seq_cap3(temp_id2seq2, temp_junc_seq2, output_file, swalign_score, juncseq)
                    #key2contig2_fml_asm[temp_key] = assemble_seq_fml_asm(temp_id2seq, temp_junc_seq, output_file)
                    #key2contig2_velvet[temp_key] = assemble_seq_velvet(temp_id2seq, temp_junc_seq, output_file)
                    #key2contig_sga[temp_key] = assemble_seq_sga(temp_id2seq, temp_junc_seq, output_file)

                temp_key2 = F[0]
                temp_id2seq2 = {}
                FF = temp_key2.split(',')
                if FF[2] == "+":
                    temp_junc_seq2 = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - swalign_length, int(FF[1]))
                else:
                    temp_junc_seq2 = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + swalign_length))
                juncseq = FF[3]

            temp_id2seq2[F[1]] = F[2]

        if len(temp_id2seq2) > 0: 
            key2contig2_cap3[temp_key2] = assemble_seq_cap3(temp_id2seq2, temp_junc_seq2, output_file, swalign_score, juncseq)
            #key2contig2_fml_asm[temp_key] = assemble_seq_fml_asm(temp_id2seq, temp_junc_seq, output_file)
            #key2contig2_velvet[temp_key] = assemble_seq_velvet(temp_id2seq, temp_junc_seq, output_file)
            #key2contig_sga[temp_key] = assemble_seq_sga(temp_id2seq, temp_junc_seq, output_file)

    # #readid2key2[pair_ID] + '\t' + ("/1" if flags[6] == "1" else "/2") + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence
    # temp_key3 = ""
    # temp_id2seq3 = {}
    # key2contig3 = {}

    # temp_key4 = ""
    # temp_id2seq4 = {}
    # key2contig4 = {}

    # with open(output_file + ".tmp2.contig3.sorted") as hin3:
    #     for line in hin3:
    #         F = line.rstrip('\n').split('\t')
    #         if F[1] == "/1":
    #             if temp_key3 != F[0]:
    #                 if len(temp_id2seq3) > 0:
    #                     key2contig3[temp_key3] = assemble_seq2(temp_id2seq3, output_file)
    #                 temp_key3 = F[0]
    #                 temp_id2seq3 = {}
    #             temp_id2seq3[F[2]] = F[3]
    #         if F[1] == "/2":
    #             if temp_key4 != F[0]:
    #                 if len(temp_id2seq4) > 0:
    #                     key2contig4[temp_key4] = assemble_seq2(temp_id2seq4, output_file)
    #                 temp_key4 = F[0]
    #                 temp_id2seq4 = {}
    #             temp_id2seq4[F[2]] = F[3]

    #     if len(temp_id2seq3) > 0: 
    #         key2contig3[temp_key3] = assemble_seq2(temp_id2seq3, output_file)

    #     if len(temp_id2seq4) > 0: 
    #         key2contig4[temp_key4] = assemble_seq2(temp_id2seq4, output_file)

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            key = ','.join(F[:4])
            
            long_contig_cap3 = key2contig_cap3[key] if key in key2contig_cap3 else "" + '\t' + "" + '\t' + ""+ '\t' + ""
            short_contig_cap3 = key2contig2_cap3[key] if key in key2contig2_cap3 else "" + '\t' + "" + '\t' + ""+ '\t' + ""
            #if len(contig_cap3) < min_contig_length: continue
            # if contig[:8] != F[3][:8]: continue

            #long_contig_fml_asm = key2contig_fml_asm[key] if key in key2contig_fml_asm else "" + '\t' + "" + '\t' + ""
            #short_contig_fml_asm = key2contig2_fml_asm[key] if key in key2contig2_fml_asm else "" + '\t' + "" + '\t' + ""
            #if len(contig_fml_asm) < min_contig_length: continue

            #long_contig_velvet = key2contig_velvet[key] if key in key2contig_velvet else "" + '\t' + "" + '\t' + ""
            #short_contig_velvet = key2contig2_velvet[key] if key in key2contig2_velvet else "" + '\t' + "" + '\t' + ""
            #if len(contig_velvet) < min_contig_length: continue


            long_contig_sga = key2contig_sga[key] if key in key2contig_sga else "" + '\t' + "" + '\t' + ""+ '\t' + ""
            #short_contig_sga = key2contig2_sga[key] if key in key2contig2_sga else "" + '\t' + "" + '\t' + ""

            #if len(contig_sga) < min_contig_length: continue

            print >> hout, '\t'.join(F) + '\t' + long_contig_cap3 + '\t' + short_contig_cap3 + '\t' + long_contig_sga
            #print >> hout, '\t'.join(F) + '\t' + long_contig_cap3 + '\t' + short_contig_cap3 + '\t' + long_contig_fml_asm + '\t' + short_contig_fml_asm  + '\t' + \
            #               long_contig_velvet + '\t' + short_contig_velvet + '\t' + long_contig_sga + '\t' + short_contig_sga
            #print >> hout, '\t'.join(F) + '\t' + contig + '\t' + long_contig + '\t' + short_contig + '\t' + pair_contig1 + '\t' + pair_contig2

    hout.close()

    subprocess.call(["rm", "-f", output_file + ".tmp2.contig.unsorted"])
    subprocess.call(["rm", "-f", output_file + ".tmp2.contig2.unsorted"])
    subprocess.call(["rm", "-f", output_file + ".tmp2.contig.sorted"])
    subprocess.call(["rm", "-f", output_file + ".tmp2.contig2.sorted"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap3.screen"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.singlets"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.info"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.contigs"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.contigs.qual"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.contigs.links"])
    subprocess.call(["rm", "-f", output_file + ".tmp3.assemble_input.fa.cap.ace"])
    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.merged.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.sai"])
    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.rsai"])
    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.bwt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp3.assemble_input.rbwt"])


def psl_check_for_human(psl_file, key2seq, align_margin = 10000): 

    tempID = ""
    temp_align2score = {}
    key2align = {}
    key2best_score = {}
    key2margin = {}
    with open(psl_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue
            
            if int(F[16]) - int(F[15]) - int(F[10]) > 5: continue
            if int(F[10]) - int(F[0]) > 5: continue
            if int(F[11]) > 5: continue


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



def alignment_contig(input_file, contig_file, output_file, reference_genome, blat_option, virus_db, repeat_db, mitochondria_db, adapter_db):
    
    blat_cmds = ("blat " + blat_option).split(' ')

    # long_contig cap3
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

    #short_contig_cap3
    key2seq2 = {}
    hout = open(output_file + ".tmp4.contig.alignment_check2.fa", 'w')
    with open(contig_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            key2seq2[key] = F[13]
            print >> hout, '>' + key
            print >> hout, F[13]
    hout.close()

    #long_contig_sga = F[27]
    key2seq3 = {}
    hout = open(output_file + ".tmp4.contig.alignment_check3.fa", 'w')
    with open(contig_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            key2seq3[key] = F[17]
            print >> hout, '>' + key
            print >> hout, F[17]

    hout.close()

    # reference genome ##############################
    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check.fa",
                           output_file + ".tmp4.contig.alignment_check.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human, key2bscore_human, key2margin_human = psl_check_for_human(output_file + ".tmp4.contig.alignment_check.psl", key2seq)


    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check2.fa",
                           output_file + ".tmp4.contig.alignment_check2.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human2, key2bscore_human2, key2margin_human2 = psl_check_for_human(output_file + ".tmp4.contig.alignment_check2.psl", key2seq)


    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(blat_cmds + [reference_genome, output_file + ".tmp4.contig.alignment_check3.fa",
                           output_file + ".tmp4.contig.alignment_check3.psl"], stdout = FNULL, stderr = subprocess.STDOUT)

    FNULL.close()
    if sret != 0:
        print >> sys.stderr, "blat error, error code: " + str(sRet)
        sys.exit()

    key2align_human3, key2bscore_human3, key2margin_human3 = psl_check_for_human(output_file + ".tmp4.contig.alignment_check3.psl", key2seq)
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
        key2align_virus, key2bscore_virus, key2margin_virus = {}, {}, {}

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
        key2align_virus2, key2bscore_virus2, key2margin_virus2 = {}, {}, {}

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
        key2align_virus3, key2bscore_virus3, key2margin_virus3 = {}, {}, {}
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
        key2align_repeat, key2bscore_repeat, key2margin_repeat = {}, {}, {} 

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
        key2align_repeat2, key2bscore_repeat2, key2margin_repeat2 = {}, {}, {} 

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
        key2align_repeat3, key2bscore_repeat3, key2margin_repeat3 = {}, {}, {} 
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
        key2align_mitochondria, key2bscore_mitochondria, key2margin_mitochondria = {}, {}, {} 

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
        key2align_mitochondria2, key2bscore_mitochondria2, key2margin_mitochondria2 = {}, {}, {} 

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
        key2align_mitochondria3, key2bscore_mitochondria3, key2margin_mitochondria3 = {}, {}, {} 
    ################################################



    # adapter_genome #################################
    if adapter_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [adapter_db, output_file + ".tmp4.contig.alignment_check.fa",
                        output_file + ".tmp4.contig.alignment_check_adapter.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_adapter, key2bscore_adapter, key2margin_adapter = psl_check(output_file + ".tmp4.contig.alignment_check_adapter.psl", key2seq)
    else:
        key2align_adapter, key2bscore_adapter, key2margin_adapter = {}, {}, {} 

    if adapter_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [adapter_db, output_file + ".tmp4.contig.alignment_check2.fa",
                        output_file + ".tmp4.contig.alignment_check_adapter2.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_adapter2, key2bscore_adapter2, key2margin_adapter2 = psl_check(output_file + ".tmp4.contig.alignment_check_adapter2.psl", key2seq)
    else:
        key2align_adapter2, key2bscore_adapter2, key2margin_adapter2 = {}, {}, {} 

    if adapter_db != "":
        FNULL = open(os.devnull, 'w')
        sret = subprocess.call(blat_cmds + [adapter_db, output_file + ".tmp4.contig.alignment_check3.fa",
                        output_file + ".tmp4.contig.alignment_check_adapter3.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)

        FNULL.close()
        if sret != 0:
            print >> sys.stderr, "blat error, error code: " + str(sRet)
            sys.exit()

        key2align_adapter3, key2bscore_adapter3, key2margin_adapter3 = psl_check(output_file + ".tmp4.contig.alignment_check_adapter3.psl", key2seq)
    else:
        key2align_adapter3, key2bscore_adapter3, key2margin_adapter3 = {}, {}, {} 
    ################################################


    hout = open(output_file, 'w')    
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n')
        print >> hout, header + '\t' + '\t'.join(["Long_Contig_Cap3", "Junc_Seq_Consistency", "Human_Alignment", "Human_Mismatch", "Human_Margin",
                                                  "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                                  "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                                  "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                                                  "Short_Contig_Cap3", "Junc_Seq_Consistency", "Human_Alignment", "Human_Mismatch", "Human_Margin",
                                                  "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                                  "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                                  "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                                                  "Long_Contig_SGA", "Junc_Seq_Consistency", "Human_Alignment", "Human_Mismatch", "Human_Margin",
                                                  "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                                  "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                                  "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                                                  ])

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])

            seq = key2seq[key] if key in key2seq and len(key2seq[key]) > 0 else "---"
            pair_seq1 = key2seq2[key] if key in key2seq2 and len(key2seq2[key]) > 0 else "---"
            pair_seq2 = key2seq3[key] if key in key2seq3 and len(key2seq3[key]) > 0 else "---"

            junc_seq_consistency1 = "TRUE" if seq[:8] == F[3][:8] else "FALSE"
            junc_seq_consistency2 = "TRUE" if pair_seq1[:8] == F[3][:8] else "FALSE"
            junc_seq_consistency3 = "TRUE" if pair_seq2[:8] == F[3][:8] else "FALSE"

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

            align_adapter = ';'.join(key2align_adapter[key]) if key in key2align_adapter and len(key2align_adapter[key]) > 0 else "---"
            bscore_adapter = str(key2bscore_adapter[key]) if key in key2bscore_adapter and key2bscore_adapter[key] != float("inf") else "---"
            margin_adapter = str(key2margin_adapter[key]) if key in key2margin_adapter and key2margin_adapter[key] != float("inf") else "---"

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

            align_adapter2 = ';'.join(key2align_adapter2[key]) if key in key2align_adapter2 and len(key2align_adapter2[key]) > 0 else "---"
            bscore_adapter2 = str(key2bscore_adapter2[key]) if key in key2bscore_adapter2 and key2bscore_adapter2[key] != float("inf") else "---"
            margin_adapter2 = str(key2margin_adapter2[key]) if key in key2margin_adapter2 and key2margin_adapter2[key] != float("inf") else "---"

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

            align_adapter3 = ';'.join(key2align_adapter3[key]) if key in key2align_adapter3 and len(key2align_adapter3[key]) > 0 else "---"
            bscore_adapter3 = str(key2bscore_adapter3[key]) if key in key2bscore_adapter3 and key2bscore_adapter3[key] != float("inf") else "---"
            margin_adapter3 = str(key2margin_adapter3[key]) if key in key2margin_adapter3 and key2margin_adapter3[key] != float("inf") else "---"


            print >> hout, '\t'.join(F) + '\t' + seq + '\t' + junc_seq_consistency1 + '\t' + align_human + '\t' + bscore_human + '\t' + margin_human + '\t' + \
                           align_virus + '\t' + bscore_virus + '\t' + margin_virus + '\t' + align_repeat + '\t' + bscore_repeat + '\t' + margin_repeat + '\t' + \
                           align_mitochondria + '\t' + bscore_mitochondria + '\t' + margin_mitochondria + '\t' + align_adapter + '\t' + bscore_adapter + '\t' + margin_adapter + '\t' + \
                           pair_seq1 + '\t' + junc_seq_consistency2 + '\t' + align_human2 + '\t' + bscore_human2 + '\t' + margin_human2 + '\t' + \
                           align_virus2 + '\t' + bscore_virus2 + '\t' + margin_virus2 + '\t' + align_repeat2 + '\t' + bscore_repeat2 + '\t' + margin_repeat2 + '\t' + \
                           align_mitochondria2 + '\t' + bscore_mitochondria2 + '\t' + margin_mitochondria2 + '\t' + align_adapter2 + '\t' + bscore_adapter2 + '\t' + margin_adapter2 + '\t' + \
                           pair_seq2 + '\t' + junc_seq_consistency3 + '\t' + align_human3 + '\t' + bscore_human3 + '\t' + margin_human3 + '\t' + \
                           align_virus3 + '\t' + bscore_virus3 + '\t' + margin_virus3 + '\t' + align_repeat3 + '\t' + bscore_repeat3 + '\t' + margin_repeat3 + '\t' + \
                           align_mitochondria3 + '\t' + bscore_mitochondria3 + '\t' + margin_mitochondria3 + '\t' + align_adapter3 + '\t' + bscore_adapter3 + '\t' + margin_adapter3

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_mitochondria.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_adapter.psl"])

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check2.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_mitochondria2.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_adapter2.psl"])

    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check3.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_virus3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_repeat3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_mitochondria3.psl"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.contig.alignment_check_adapter3.psl"])


def annotate_break_point(input_file, output_file, genome_id, is_grc):

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, False)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, False)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    hout = open(output_file, 'w')
    header2ind = {}
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        # for (i, cname) in enumerate(header):
        #     header2ind[cname] = i

        print >> hout, '\t'.join(["Chr", "Pos", "Dir", "Junc_Seq", "Gene", "Exon"] + header[4:])
        for line in hin:
            F = line.rstrip('\n').split('\t')

            ##########
            # check gene annotation
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]) + 1)
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
                records = exon_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]) + 1)
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

            print >> hout, '\t'.join(F[:4]) + '\t' + \
                           ','.join(gene) + '\t' + ';'.join(exon) + '\t' + '\t'.join(F[4:])

    hin.close()
    hout.close()
    gene_tb.close()
    exon_tb.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])
    