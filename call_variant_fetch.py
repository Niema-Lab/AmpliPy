import pysam
import sys
import numpy as np
import getopt
#import pdb
# 'Usage: python3 call_variant_fetch.py AlignmentFile ReferenceFile VariantOutputFilename'


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
        optlist, args = getopt.getopt(sys.argv[1:], 'q:t:m:')
    except getopt.GetoptError as err:
        print('Usage: python3 call_variant_fetch.py AlignmentFile ReferenceFile VariantOutputFilename '
              '[-q minimal quality score to count base (Default:20)] '
              '[-t minimal frequency threshold to call variant (Default:0.03)] '
              '[-m minimal number of reads to call variant (Default:0)] ')
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    if len(args) < 3:
        print('You should include AlignmentFile, ReferenceFile, and VariantOutputFilename')
        print('Usage: python3 call_variant_fetch.py AlignmentFile ReferenceFile VariantOutputFilename '
              '[-q minimal quality score to count base (Default:20)] '
              '[-t minimal frequency threshold to call variant (Default:0.03)] '
              '[-m minimal number of reads to call variant (Default:0)] ')
        quit()

    try:
        # read bam file
        bam = pysam.AlignmentFile(args[0], "rb")
    except FileNotFoundError:
        print("Did not found the alignment file")
        quit()

    # get reference sequence
    reference_file = args[1]
    reference_sequence = parse_reference(reference_file)
    # decide min base quality
    min_qual = 20;
    # decide min number of reads
    min_depth = 0
    # decide min threshold for (number of base/number of reads)
    min_threshold = 0.03

    ##if the user decide not to use the default value, change the parameter according to the command line input options
    try:
        for op, val in optlist:
            if op == '-q':
                min_qual = int(val)
            elif op == '-t':
                min_threshold = float(val)
            elif op == '-m':
                min_depth = int(val)
    except ValueError:
        print("please give valid inputs for options -q (int) -t (float) -m (int)")
        quit()

    buckets_array = []
    for i in range(len(reference_sequence)):
        buckets_array.append(bucket())

    for alignedSeg in bam.fetch(until_eof=True):
        # query_pos: position in query sequence
        # start with query start position (the first position after soft-clipped bases)
        query_pos = alignedSeg.query_alignment_start
        # ref_pos: position in the reference sequence
        # start with reference start position
        ref_pos = alignedSeg.reference_start
        # reference_end: end position of aligned reference
        reference_end = alignedSeg.reference_end
        # query_alignment_end: end position of aligned query
        query_alignment_end = alignedSeg.query_alignment_end
        # pos: index in aligned_pairs
        pos = query_pos
        # query_sequence: query sequence (including soft-clipped bases)
        query_sequence = alignedSeg.query_sequence
        # aligned_pairs: a list of pairs, each pair is (query_pos, ref_pos)
        # in case of insertion and deletion, one of them is None
        aligned_pairs = alignedSeg.get_aligned_pairs()
        # quals: read sequence base qualities, including soft-clipped bases
        quals = alignedSeg.query_qualities

        # as long as ref_pos and query_pos not reach the ends
        while ref_pos < reference_end and query_pos < query_alignment_end:

            # check first to avoid error
            if pos + 1 < len(aligned_pairs):

                # insertion
                if aligned_pairs[pos + 1][1] is None and aligned_pairs[pos + 1][0] < query_alignment_end:

                    # total_depth is the total number of aligned read
                    buckets_array[ref_pos].total_depth += 1

                    # in case of insertion, ref_pos will be None and query_pos will advance
                    pos += 1
                    query_pos += 1

                    b = ''

                    # collect all bases when reference is None
                    while aligned_pairs[pos][1] is None and query_pos < query_alignment_end:
                        b = b + query_sequence[query_pos]
                        pos += 1
                        query_pos += 1

                    b = '+' + b

                    # put [insertion pattern, num of this insertion happened, sum(qualities), number of reverse]
                    # into bucket_array[ref_pos].insertN
                    exist = False
                    for subarr in buckets_array[ref_pos].insertN:
                        if subarr[0] == b[1:].upper():
                            subarr[1] += 1
                            subarr[2] += min_qual
                            if alignedSeg.is_reverse:
                                subarr[3] += 1
                            exist = True
                            break
                    if not exist:
                        subarr = []
                        subarr.append(b[1:].upper())
                        subarr.append(1)
                        subarr.append(min_qual)
                        if alignedSeg.is_reverse:
                            subarr.append(1)
                        else:
                            subarr.append(0)
                        buckets_array[ref_pos].insertN.append(subarr)

                    # query_pos already ready, no need to change. advance ref_pos by 1
                    ref_pos += 1
                    continue

                ##deletion
                elif aligned_pairs[pos + 1][0] is None and aligned_pairs[pos + 1][1] < reference_end:

                    # total_depth is the total number of aligned read
                    buckets_array[ref_pos].total_depth += 1

                    # in case of deletion, query_pos will be None and ref_pos will advance
                    pos += 1
                    ref_pos_save = ref_pos
                    ref_pos += 1

                    b = ''

                    # same with before
                    while aligned_pairs[pos][0] is None and ref_pos < reference_end:
                        b = b + reference_sequence[ref_pos]
                        pos += 1
                        ref_pos += 1

                    b = '-' + b

                    exist = False
                    for subarr in buckets_array[ref_pos_save].deleteN:
                        if subarr[0] == b[1:].upper():
                            subarr[1] += 1
                            subarr[2] += min_qual
                            if alignedSeg.is_reverse:
                                subarr[3] += 1
                            exist = True
                            break
                    if not exist:
                        subarr = []
                        subarr.append(b[1:].upper())
                        subarr.append(1)
                        subarr.append(min_qual)
                        if alignedSeg.is_reverse:
                            subarr.append(1)
                        else:
                            subarr.append(0)
                        buckets_array[ref_pos_save].deleteN.append(subarr)

                    # ref_pos already ready, no need to change. advance query_pos by 1
                    query_pos += 1
                    continue

            # if the base quality is less than the min_qual, skip it
            if quals[query_pos] < min_qual:

                buckets_array[ref_pos].total_depth += 1
                pos += 1
                query_pos += 1
                ref_pos += 1
                continue

            # if there is no insertion, deletion, then it is ACGT(N)
            else:
                buckets_array[ref_pos].total_depth += 1
                query = query_sequence[query_pos]

                # each bucket has A,G,C,T array, which is
                # [total number of this base(A+a), total number of reverse of this base(a), total quality score added]
                # 'ACGTN' are forward read, 'acgtn' are backward read
                # A,C,G,T [number of this base appears, sum(qualities), number of reverse query]
                if query == 'A':
                    buckets_array[ref_pos].A[0] += 1
                    buckets_array[ref_pos].A[2] += quals[query_pos]
                    if alignedSeg.is_reverse:
                        buckets_array[ref_pos].A[1] += 1
                elif query == 'C':
                    buckets_array[ref_pos].C[0] += 1
                    buckets_array[ref_pos].C[2] += quals[query_pos]
                    if alignedSeg.is_reverse:
                        buckets_array[ref_pos].C[1] += 1
                elif query == 'G':
                    buckets_array[ref_pos].G[0] += 1
                    buckets_array[ref_pos].G[2] += quals[query_pos]
                    if alignedSeg.is_reverse:
                        buckets_array[ref_pos].G[1] += 1
                elif query == 'T':
                    buckets_array[ref_pos].T[0] += 1
                    buckets_array[ref_pos].T[2] += quals[query_pos]
                    if alignedSeg.is_reverse:
                        buckets_array[ref_pos].T[1] += 1
                # TODO: what does N present? should I count them?
                elif query == 'N' or query == 'n':
                    print("ref_pos: ", ref_pos, " ", query)
                else:
                    print("not included query: " + query)
                pos += 1
                query_pos += 1
                ref_pos += 1

    # write to output file
    outf = open(args[2], 'w')
    outf.write(
        "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\n")
    region = bam.references[0]
    reference_sequence = parse_reference(reference_file)
    for ref_pos in range(len(reference_sequence)):

        ref = reference_sequence[ref_pos]
        b = buckets_array[ref_pos]

        pdepth = 0
        pdepth += b.A[0]
        pdepth += b.C[0]
        pdepth += b.G[0]
        pdepth += b.T[0]

        if pdepth < min_depth:
            continue

        if pdepth <= 0:
            continue

        A_freq_depth_0 = float(b.A[0]) / pdepth
        C_freq_depth_0 = float(b.C[0]) / pdepth
        G_freq_depth_0 = float(b.G[0]) / pdepth
        T_freq_depth_0 = float(b.T[0]) / pdepth

        # ref is the reference base, set matchN to an array that has the same base as the reference base
        matchN = [0, 0, 0.0]
        if ref == 'A':
            matchN = b.A
        elif ref == 'C':
            matchN = b.C
        elif ref == 'G':
            matchN = b.G
        elif ref == 'T':
            matchN = b.T

        # calculate reference mean quality
        if matchN[0] == 0:
            ref_mean_qual = 0
        else:
            ref_mean_qual = int(matchN[2] / matchN[0])

        mdepth = b.total_depth

        #every subarr inside b.insertN/deleteN, print any insertion/deletion sequence that exceed min_threshold
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
        # if any ACGT base exceed min_threshold, print them
        if ref != 'A' and A_freq_depth_0 >= min_threshold:
            #outf.write(str(region) + '\t' + str(ref_pos) + '\t' + '.' + '\t' + ref + '\t' + 'A' + '\t' +  str(int(b.A[2]/b.A[0])) + '\n')
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