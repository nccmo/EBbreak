#! /usr/bin/env python

import sys, subprocess
import parse

def parse_main(args):

    parse.parse_bp_from_bam(args.bam_file, args.output_file + ".bp.tmp.txt", args.key_seq_size, args.min_major_clip_size, args.max_minor_clip_size)

    parse.cluster_breakpoint(args.output_file + ".bp.tmp.txt", args.output_file + ".bp.clustered.tmp.txt", args.check_interval)

    hout = open(args.output_file + ".bp.clustered.sorted.tmp.txt", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", args.output_file + ".bp.clustered.tmp.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting merged junction file"
        sys.exit(1)

    hout = open(args.output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", args.output_file + ".bp.clustered.sorted.tmp.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in compression merged junction file"
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "vcf", args.output_file])
    if s_ret != 0:
        print >> sys.stderr, "Error in indexing merged junction file"
        sys.exit(1)

    subprocess.call(["rm", "-f", args.output_file + ".bp.tmp.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".bp.clustered.tmp.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".bp.clustered.sorted.tmp.txt"])
    
