import pysam
import sys
import numpy as np
import getopt


class bucket:

    def __init__(self):
        self.total_depth = 0
        self.matchN = [0, 0, 0.0]
        self.A = [0, 0, 0.0]
        self.G = [0, 0, 0.0]
        self.C = [0, 0, 0.0]
        self.T = [0, 0, 0.0]
        self.insertN = list()
        self.deleteN = list()


def parse_reference(ref_file):
    rstr = ''
    try:
        f = open(ref_file, 'r')
    except FileNotFoundError:
        print('Usage: python3 call_variant_fetch.py AlignmentFile ReferenceFile VariantOutputFilename')
        quit()
    lines = f.readlines()
    for l in lines:
        if l[0] == '>':
            continue
        else:
            rstr += l
    rstr = rstr.replace('\n', '')
    return rstr


def variant():
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'q:t:m:Q:')
    except getopt.GetoptError as err:
        print('Usage: python3 call_variant_pileup.py AlignmentFile ReferenceFile VariantOutputFilename '
              '[-q minimal quality score to count base (Default:20)] '
              '[-t minimal frequency threshold to call variant (Default:0.03)] '
              '[-m minimal number of reads to call variant (Default:0)] '
              '[-Q minimal quality score for a base to be considered during pileup (Default: 0)]')
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    if len(args) < 3:
        print('You should include AlignmentFile, ReferenceFile, and VariantOutputFilename')
        print('Usage: python3 call_variant_pileup.py AlignmentFile ReferenceFile VariantOutputFilename '
              '[-q minimal quality score to count base (Default:20)] '
              '[-t minimal frequency threshold to call variant (Default:0.03)] '
              '[-m minimal number of reads to call variant (Default:0)] '
              '[-Q minimal quality score for a base to be considered during pileup (Default: 0)]')
        quit()

    try:
        # read bam file
        bam = pysam.AlignmentFile(args[0], "rb")
    except FileNotFoundError:
        print("Did not found the alignment file")
        quit()

    ##the path of the reference file
    reference_file = args[1]
    ##reference_sequence: an array of ACGT of the reference sequence itself
    reference_sequence = parse_reference(reference_file)

    ##decide min base quality
    min_qual = 20;
    ##decide min number of reads
    min_depth = 0
    ##decide min threshold for (number of base/number of reads)
    min_threshold = 0.03
    ##decide min base quality for a base to be considered in the pileup
    min_base_quality_pileup = 0

    ##if the user decide not to use the default value, change the parameter according to the command line input options
    try:
        for op, val in optlist:
            if op == '-q':
                min_qual = int(val)
            elif op == '-t':
                min_threshold = float(val)
            elif op == '-m':
                min_depth = int(val)
            elif op == '-Q':
                min_base_quality_pileup = int(val)
    except ValueError:
        print("please give valid inputs for options -q (int) -t (float) -m (int) or -Q (int)")
        quit()

    ##ref_length: the length of the reference sequence
    ref_length = bam.lengths

    ##an array of length ref_length that hold bucket for each reference position
    ##TODO: currently can only works with 1 reference
    buckets_array = []
    for i in range(ref_length[0]):
        buckets_array.append(bucket())

    ##region: the name of the reference sequence
    region = bam.references[0]

    ##write to output file
    outf = open(args[2], 'w')
    outf.write(
        "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\n")

    ref_file = pysam.FastaFile(filename=reference_file)

    ##TODO: in this example, stepper causes no difference, how about other sample that has flags set?
    for pc in bam.pileup(ignore_orphans=False, min_base_quality=min_base_quality_pileup, compute_baq=False,
                         max_depth=600000,
                         stepper='samtools', fastafile=ref_file):
        ##qqual: an array of query bases qualities
        qqual = pc.get_query_qualities()
        ##qseq: an array of query sequences aligned at this pileupcolumn
        qseq = pc.get_query_sequences(mark_ends=True, add_indels=True)
        ##ref_pos: the position on the reference
        ref_pos = pc.reference_pos

        ##mdepth: how many reads are aligned at this ref_pos
        ##len(qseq) != PileupColumn.get_num_aligned(), the latter one is same to the mpileup in samstool
        mdepth = pc.get_num_aligned()

        ##if there is no read aligned, skip this ref_pos and keep everything 0
        if mdepth == 0:
            continue

        ##otherwise add the total_depth of this specific ref_pos by mdepth
        buckets_array[ref_pos].total_depth += mdepth

        ##i: position in qqual
        i = -1
        ##for every base in the qseq
        for b in qseq:
            i += 1

            ##if this base start with ^, erase ^. skip this read if there are no bases following ^
            if b[0] == '^':
                if len(b) == 1:
                    i = i - 1
                    continue
                elif len(b) > 2:
                    b = b[2:]
                ##Testing
                else:
                    print(b)
                    continue
            ##Testing
            if b == '$':
                print("ref_pos: ", ref_pos, " ", b)
                continue

            ##if there is a char $, delete the char
            b = b.replace('$', '')

            ##for whatever reason b goes to '', skip the read
            if b == '':
                continue

            ##put insertion and deletion before the min_qual filter to let +/- all have min_qual

            ##in case of addition(+) or deletion(-)
            ##format: reference_base +/- bases
            if len(b) > 1:
                if b[1] == '+':
                    exist = False
                    ##insertN is an array of array
                    ##subarr is [the specific bases added, the total number this bases showed up, total quality score added, the number of reverse read this pattern has]
                    ##TODO: implement reverse strand count
                    for subarr in buckets_array[ref_pos].insertN:
                        ##if this pattern exists
                        ##Store uppercase string, use subarr[3] to record reverse read (if b[1:] are lower cases then reverse read)
                        if subarr[0] == b[1:].upper():
                            subarr[1] += 1
                            subarr[2] += min_qual
                            if not b[1:].isupper():
                                subarr[3] += 1
                            exist = True
                            break
                    ##if this pattern does not exist
                    if not exist:
                        subarr = list()
                        subarr.append(b[1:].upper())
                        subarr.append(1)
                        subarr.append(min_qual)
                        if not b[1:].isupper():
                            subarr.append(1)
                        else:
                            subarr.append(0)
                        buckets_array[ref_pos].insertN.append(subarr)
                elif b[1] == '-':
                    exist = False
                    for subarr in buckets_array[ref_pos].deleteN:
                        if subarr[0] == b[1:].upper():
                            subarr[1] += 1
                            subarr[2] += min_qual
                            if not b[1:].isupper():
                                subarr[3] += 1
                            exist = True
                            break
                    if not exist:
                        subarr = list()
                        subarr.append(b[1:].upper())
                        subarr.append(1)
                        subarr.append(min_qual)
                        if not b[1:].isupper():
                            subarr.append(1)
                        else:
                            subarr.append(0)
                        buckets_array[ref_pos].deleteN.append(subarr)

            ##filter base with too low quality score
            elif qqual[i] < min_qual:
                continue

            ##each bucket has A,G,C,T array, which is
            ##[total number of this base(A+a), total number of reverse of this base(a), total quality score added]
            ##'ACGTN' are forward read, 'acgtn' are backward read
            elif b == 'A':
                buckets_array[ref_pos].A[0] += 1
                buckets_array[ref_pos].A[2] += qqual[i]
            elif b == 'a':
                buckets_array[ref_pos].A[0] += 1
                buckets_array[ref_pos].A[1] += 1
                buckets_array[ref_pos].A[2] += qqual[i]
            elif b == 'C':
                buckets_array[ref_pos].C[0] += 1
                buckets_array[ref_pos].C[2] += qqual[i]
            elif b == 'c':
                buckets_array[ref_pos].C[0] += 1
                buckets_array[ref_pos].C[1] += 1
                buckets_array[ref_pos].C[2] += qqual[i]
            elif b == 'G':
                buckets_array[ref_pos].G[0] += 1
                buckets_array[ref_pos].G[2] += qqual[i]
            elif b == 'g':
                buckets_array[ref_pos].G[0] += 1
                buckets_array[ref_pos].G[1] += 1
                buckets_array[ref_pos].G[2] += qqual[i]
            elif b == 'T':
                buckets_array[ref_pos].T[0] += 1
                buckets_array[ref_pos].T[2] += qqual[i]
            elif b == 't':
                buckets_array[ref_pos].T[0] += 1
                buckets_array[ref_pos].T[1] += 1
                buckets_array[ref_pos].T[2] += qqual[i]
            ##TODO: what does N present? should I count them?
            elif b == 'N' or b == 'n':
                print("ref_pos: ", ref_pos, " ", b)
            ##TODO: what is reference deletion, what to do in reference deletion?
            elif b[0] == '*':
                continue
            ##Testing
            else:
                print(b)

    # ref_length: reference length of every reference file
    for ref_pos in range(ref_length[0]):

        ##ref: reference base
        ##b: bucket
        ref = reference_sequence[ref_pos]
        b = buckets_array[ref_pos]

        ##calculate pdepth, used in case of insertion/deletion
        pdepth = 0
        pdepth += b.A[0]
        pdepth += b.C[0]
        pdepth += b.G[0]
        pdepth += b.T[0]

        ##check if pdepth pass threshold
        if pdepth < min_depth:
            continue
        ##in case min_depth == 0, avoid error of dividing something by 0
        if pdepth <= 0:
            continue

        ##calculate frequency depth of alternate bases ACGT
        A_freq_depth_0 = float(b.A[0]) / pdepth
        C_freq_depth_0 = float(b.C[0]) / pdepth
        G_freq_depth_0 = float(b.G[0]) / pdepth
        T_freq_depth_0 = float(b.T[0]) / pdepth

        ##ref is the reference base, set matchN to an array that has the same base as the reference base
        matchN = [0, 0, 0.0]
        if ref == 'A':
            matchN = b.A
        elif ref == 'C':
            matchN = b.C
        elif ref == 'G':
            matchN = b.G
        elif ref == 'T':
            matchN = b.T

        ##calculate reference mean quality
        if matchN[0] == 0:
            ref_mean_qual = 0
        else:
            ref_mean_qual = int(matchN[2] / matchN[0])

        # mdepth: number of reads aligned
        mdepth = b.total_depth

        ##if any in/del exceed min_threshold, print them
        for subarr in b.insertN:
            freq_depth_temp = float(subarr[1]) / mdepth
            if freq_depth_temp >= min_threshold:
                outf.write(str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + subarr[0] + '\t' + str(
                    matchN[0]) + '\t' + str(b.matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(
                    subarr[1]) + '\t' + str(subarr[3]) + '\t' + str(int(subarr[2] / subarr[1])) + '\t' + str(
                    round(freq_depth_temp, 6)) + '\t' + str(mdepth) + '\n')
        for subarr in b.deleteN:
            freq_depth_temp = float(subarr[1]) / mdepth
            if freq_depth_temp >= min_threshold:
                outf.write(str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + subarr[0] + '\t' + str(
                    matchN[0]) + '\t' + str(b.matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(
                    subarr[1]) + '\t' + str(subarr[3]) + '\t' + str(int(subarr[2] / subarr[1])) + '\t' + str(
                    round(freq_depth_temp, 6)) + '\t' + str(mdepth) + '\n')
        ##if any ACGT base exceed min_threshold, print them
        if ref != 'A' and A_freq_depth_0 >= min_threshold:
            ##TODO: print in VCF format
            ##outf.write(str(region) + '\t' + str(ref_pos) + '\t' + '.' + '\t' + ref + '\t' + 'A' + '\t' +  str(int(b.A[2]/b.A[0])) + '\n')
            outf.write(
                str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + 'A' + '\t' + str(matchN[0]) + '\t' + str(
                    matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(b.A[0]) + '\t' + str(b.A[1]) + '\t' + str(
                    int(b.A[2] / b.A[0])) + '\t' + str(round(A_freq_depth_0, 6)) + '\t' + str(pdepth) + '\n')
        if ref != 'C' and C_freq_depth_0 >= min_threshold:
            outf.write(
                str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + 'C' + '\t' + str(matchN[0]) + '\t' + str(
                    matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(b.C[0]) + '\t' + str(b.C[1]) + '\t' + str(
                    int(b.C[2] / b.C[0])) + '\t' + str(round(C_freq_depth_0, 6)) + '\t' + str(pdepth) + '\n')
        if ref != 'G' and G_freq_depth_0 >= min_threshold:
            outf.write(
                str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + 'G' + '\t' + str(matchN[0]) + '\t' + str(
                    matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(b.G[0]) + '\t' + str(b.G[1]) + '\t' + str(
                    int(b.G[2] / b.G[0])) + '\t' + str(round(G_freq_depth_0, 6)) + '\t' + str(pdepth) + '\n')
        if ref != 'T' and T_freq_depth_0 >= min_threshold:
            outf.write(
                str(region) + '\t' + str(ref_pos + 1) + '\t' + ref + '\t' + 'T' + '\t' + str(matchN[0]) + '\t' + str(
                    matchN[1]) + '\t' + str(ref_mean_qual) + '\t' + str(b.T[0]) + '\t' + str(b.T[1]) + '\t' + str(
                    int(b.T[2] / b.T[0])) + '\t' + str(round(T_freq_depth_0, 6)) + '\t' + str(pdepth) + '\n')


import time

start_time = time.time()
variant()
print("--- %s seconds ---" % (time.time() - start_time))
