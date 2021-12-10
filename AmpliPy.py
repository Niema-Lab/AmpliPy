#! /usr/bin/env python3
'''
AmpliPy: Python toolkit for viral amplicon sequencing
'''

# imports
import argparse
import gzip
import pickle
import pysam
from collections import deque
from datetime import datetime
from os.path import isfile
from sys import argv, stderr

# constants
VERSION = '0.0.1'
BUFSIZE = 1048576 # 1 MB
PROGRESS_NUM_READS = 100000

# CIGAR operations
CIGAR = ['BAM_CMATCH', 'BAM_CINS', 'BAM_CDEL', 'BAM_CREF_SKIP', 'BAM_CSOFT_CLIP', 'BAM_CHARD_CLIP', 'BAM_CPAD', 'BAM_CEQUAL', 'BAM_CDIFF', 'BAM_CBACK']
BAM_CMATCH     = 0
BAM_CINS       = 1
BAM_CDEL       = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD       = 6
BAM_CEQUAL     = 7
BAM_CDIFF      = 8
CONSUME_QUERY = [True, True,  False, False, True,  False, False, True, True] # CONSUME_QUERY[i] = True if CIGAR operation i consumes letters from query
CONSUME_REF   = [True, False, True,  True,  False, False, False, True, True] # CONSUME_REF[i]   = True if CIGAR operation i consumes letters from reference

# messages
ERROR_TEXT_EMPTY_BED = "Empty BED file"
ERROR_TEXT_FILE_EXISTS = "File already exists"
ERROR_TEXT_FILE_NOT_FOUND = "File not found"
ERROR_TEXT_INVALID_AMPLIPY_INDEX_EXTENSION = "Invalid AmpliPy index file extension (should be .pkl or .pkl.gz)"
ERROR_TEXT_INVALID_BED_LINE = "Invalid primer BED line"
ERROR_TEXT_INVALID_FASTA = "Invalid FASTA file"
ERROR_TEXT_INVALID_READ_EXTENSION = "Invalid read mapping extension (should be .sam or .bam)"
ERROR_TEXT_MULTIPLE_REF_SEQS = "Multiple sequences in FASTA file"
HELP_TEXT_AMPLIPY_INDEX = "AmpliPy Index (PKL)"
HELP_TEXT_CONSENSUS = "Consensus Sequence (FASTA)"
HELP_TEXT_PRIMER = "Primer File (BED)"
HELP_TEXT_READS_UNTRIMMED = "Untrimmed Reads (SAM/BAM)"
HELP_TEXT_READS_TRIMMED = "Trimmed Reads (SAM/BAM)"
HELP_TEXT_REFERENCE = "Reference Genome (FASTA)"
HELP_TEXT_TRIM_INCLUDE_READS_NO_PRIMER = "Include reads with no primers"
HELP_TEXT_TRIM_MIN_LENGTH = "Minimum length of read to retain after trimming"
HELP_TEXT_TRIM_MIN_QUAL = "Minimum quality threshold for sliding window to pass"
HELP_TEXT_TRIM_PRIMER_POS_OFFSET = "Primer position offset. Reads that occur at the specified offset positions relative to primer positions will also be trimmed"
HELP_TEXT_TRIM_SLIDING_WINDOW_WIDTH = "Width of sliding window"
HELP_TEXT_VARIANTS = "Variant Calls (VCF)"

# print log
def print_log(s='', end='\n'):
    print("[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"),s), end=end, file=stderr); stderr.flush()

# error message
def error(s=None):
    if s is None:
        print_log("ERROR")
    else:
        print_log("ERROR: %s" % s)
    exit(1)

# parse user args
def parse_args():
    # prep arg parser
    if len(argv) == 1:
        argv.append('-h')
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # AmpliPy Trim args
    trim_parser = subparsers.add_parser("trim", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    trim_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_UNTRIMMED)
    trim_parser.add_argument('-p', '--primer', required=True, type=str, help=HELP_TEXT_PRIMER)
    trim_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    trim_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_READS_TRIMMED)
    trim_parser.add_argument('-x', '--primer_pos_offset', required=False, type=int, default=0, help=HELP_TEXT_TRIM_PRIMER_POS_OFFSET)
    trim_parser.add_argument('-m', '--min_length', required=False, type=int, default=30, help=HELP_TEXT_TRIM_MIN_LENGTH)
    trim_parser.add_argument('-q', '--min_quality', required=False, type=int, default=20, help=HELP_TEXT_TRIM_MIN_QUAL)
    trim_parser.add_argument('-s', '--sliding_window_width', required=False, type=int, default=4, help=HELP_TEXT_TRIM_SLIDING_WINDOW_WIDTH)
    trim_parser.add_argument('-e', '--include_no_primer', action='store_true', help=HELP_TEXT_TRIM_INCLUDE_READS_NO_PRIMER)

    # AmpliPy Variants args
    variants_parser = subparsers.add_parser("variants", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    variants_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    variants_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    variants_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_VARIANTS)

    # AmpliPy Consensus args
    consensus_parser = subparsers.add_parser("consensus", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    consensus_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    consensus_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    consensus_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_CONSENSUS)

    # AmpliPy AIO (All-In-One) args
    aio_parser = subparsers.add_parser("aio", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aio_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_UNTRIMMED)
    aio_parser.add_argument('-p', '--primer', required=True, type=str, help=HELP_TEXT_PRIMER)
    aio_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    aio_parser.add_argument('-ot', '--output_trimmed_reads', required=True, type=str, help=HELP_TEXT_READS_TRIMMED)
    aio_parser.add_argument('-ov', '--output_variants', required=True, type=str, help=HELP_TEXT_VARIANTS)
    aio_parser.add_argument('-oc', '--output_consensus', required=True, type=str, help=HELP_TEXT_CONSENSUS)

    # parse args and return
    return parser.parse_args()

# find overlapping primers
def find_overlapping_primers(ref_genome_length, primers):
    '''Find all primers that cover every index of the reference genome

    Args:
        ``ref_genome_length`` (``int``): Length of the reference genome

        ``primers`` (``list`` of ``(int,int)``): List of primer indices, where each primer is represented as a ``(start,end)`` tuple. Note that this uses 0-based indices, ``start`` is **in**clusive, and ``end`` is **ex**clusive. For example, `(0,100)` is the 100-length window starting at index 0 and ending at index 99.

    Returns:
        ``list`` of ``list`` of ``int``: Let ``overlapping_primers`` denote the return value. For each position ``ref_pos`` of the reference genome, ``overlapping_primers[ref_pos]`` is a list of all primers (represented as their index in ``primers``) that cover ``ref_pos``.
    '''
    # set things up
    overlapping_primers = [list() for _ in range(ref_genome_length)]
    num_primers = len(primers); curr_primers = deque(); i = 0 # i = current index of primers; primers[i] is the "current" (start,end) tuple

    # compute overlapping primers for each position of the reference genome
    for ref_pos in range(ref_genome_length):
        # remove primers that have now been passed
        while len(curr_primers) != 0 and ref_pos >= primers[curr_primers[0]][1]:
            curr_primers.popleft()

        # add primers that have now been entered
        while i < num_primers and ref_pos >= primers[i][0]:
            curr_primers.append(i); i += 1

        # the primers in curr_primers must span ref_pos
        if len(curr_primers) != 0:
            overlapping_primers[ref_pos] = list(curr_primers)
        ref_pos += 1
    return overlapping_primers

# load reference genome
def load_ref_genome(reference_fn):
    '''Load reference genome from FASTA

    Args:
        ``reference_fn`` (``str``): Filename of input reference genome FASTA

    Returns:
        ``str``: Reference genome ID

        ``str``: Reference genome sequence
    '''
    if not isfile(reference_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_NOT_FOUND, reference_fn))
    f = open(reference_fn, mode='r', buffering=BUFSIZE); ref_lines = f.read().strip().splitlines(); f.close()
    if len(ref_lines) < 2 or not ref_lines[0].startswith('>'):
        error("%s: %s" % (ERROR_TEXT_INVALID_FASTA, reference_fn))
    ref_genome_ID = ref_lines[0][1:].strip()
    ref_genome_sequence = ''.join(ref_lines[i] for i in range(1,len(ref_lines)))
    if '>' in ref_genome_sequence:
        error("%s: %s" % (ERROR_TEXT_MULTIPLE_REF_SEQS, reference_fn))
    return ref_genome_ID, ref_genome_sequence

# load primers
def load_primers(primer_fn):
    '''Load primer regions from BED

    Args:
        ``primer_fn`` (``str``): Filename of input primer BED

    Returns:
    '''
    if not isfile(primer_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_NOT_FOUND, primer_fn))
    f = open(primer_fn, mode='r', buffering=BUFSIZE); primer_lines = f.read().strip().splitlines(); f.close()
    primers = list()
    for l in primer_lines:
        try:
            curr_ref, curr_start, curr_end, curr_name = l.split('\t')
            primers.append((int(curr_start), int(curr_end)))
        except:
            error("%s: %s" % (ERROR_TEXT_INVALID_BED_LINE, l))
    if len(primers) == 0:
        error("%s: %s" % (ERROR_TEXT_EMPTY_BED, primer_fn))
    primers.sort() # shouldn't be necessary (BED should be sorted), but just in case
    return primers

# write AmpliPy index
def write_amplipy_index(amplipy_index_fn, amplipy_index_tuple):
    if not amplipy_index_fn.lower().endswith('.pkl') and not amplipy_index_fn.lower().endswith('.pkl.gz'):
        error("%s: %s" % (ERROR_TEXT_INVALID_AMPLIPY_INDEX_EXTENSION, amplipy_index_fn))
    if isfile(amplipy_index_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_EXISTS, amplipy_index_fn))
    if amplipy_index_fn.lower().endswith('.gz'):
        f = gzip.open(amplipy_index_fn, mode='wb', compresslevel=9)
    else:
        f = open(amplipy_index_fn, mode='wb', buffering=BUFSIZE)
    pickle.dump(amplipy_index_tuple, f); f.close()

# run AmpliPy Index
def run_index(primer_fn, reference_fn, amplipy_index_fn):
    '''Run AmpliPy Index

    Args:
        ``primer_fn`` (``str``): Filename of input primer BED

        ``reference_fn`` (``str``): Filename of input reference genome FASTA

        ``amplipy_index_fn`` (``str``): Filename of output AmpliPy index PKL
    '''
    print_log("Executing AmpliPy Index (v%s)" % VERSION)
    print_log("Loading reference genome: %s" % reference_fn)
    ref_genome_ID, ref_genome_sequence = load_ref_genome(reference_fn)
    print_log("Loading primers: %s" % primer_fn)
    primers = load_primers(primer_fn)
    print_log("Indexing primers...")
    overlapping_primers = find_overlapping_primers(len(ref_genome_sequence), primers)
    print_log("Writing AmpliPy index to file...")
    amplipy_index_tuple = (ref_genome_ID, ref_genome_sequence, primers, overlapping_primers)
    write_amplipy_index(amplipy_index_fn, amplipy_index_tuple)
    print_log("AmpliPy index successfully written to file: %s" % amplipy_index_fn)
    return amplipy_index_tuple

# create AlignmentFile objects for SAM/BAM input and output files
def create_AlignmentFile_objects(untrimmed_reads_fn, trimmed_reads_fn):
    '''Create ``pysam.AlignmentFile`` objects for the input and output SAM/BAM files

    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of input untrimmed reads SAM/BAM

        ``trimmed_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM

    Returns:
        ``pysam.AlignmentFile``: Stream for reading input untrimmed reads SAM/BAM

        ``pysam.AlignmentFile``: Stream for writing output trimmed reads SAM/BAM
    '''
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    if untrimmed_reads_fn.lower() == 'stdin':
        in_aln = pysam.AlignmentFile('-', 'r') # standard input --> SAM
    elif not isfile(untrimmed_reads_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_NOT_FOUND, untrimmed_reads_fn))
    elif untrimmed_reads_fn.lower().endswith('.sam'):
        in_aln = pysam.AlignmentFile(untrimmed_reads_fn, 'r')
    elif untrimmed_reads_fn.lower().endswith('.bam'):
        in_aln = pysam.AlignmentFile(untrimmed_reads_fn, 'rb')
    else:
        error("%s: %s" % (ERROR_TEXT_INVALID_READ_EXTENSION, untrimmed_reads_fn))
    if trimmed_reads_fn.lower() == 'stdout':
        out_aln = pysam.AlignmentFile('-', 'w', template=in_aln) # standard output --> SAM
    elif isfile(trimmed_reads_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_EXISTS, trimmed_reads_fn))
    elif trimmed_reads_fn.lower().endswith('.sam'):
        out_aln = pysam.AlignmentFile(trimmed_reads_fn, 'w', template=in_aln)
    elif trimmed_reads_fn.lower().endswith('.bam'):
        out_aln = pysam.AlignmentFile(trimmed_reads_fn, 'wb', template=in_aln)
    else:
        error("%s: %s" % (ERROR_TEXT_INVALID_READ_EXTENSION, trimmed_reads_fn))
    pysam.set_verbosity(tmp) # re-enable htslib verbosity
    return in_aln, out_aln

# get the query position that mapped to a given reference position
def ref_pos_to_query_pos(s, ref_pos, mode=0):
    '''Given a reference position and a query alignment, find the position in the query that mapped to that reference position

    Args:
        ``s`` (``pysam.AlignedSegment``): The mapped read

        ``ref_pos`` (``int``): The reference position

        ``mode`` (``int``): How to handle if the reference position doesn't exactly match the query (e.g. deletion). 0 = must match ref_pos, 1 = match ref_pos or first position > ref_pos, 2 = match ref_pos or last position < ref_pos
    '''
    aligned_pair_iter = s.get_aligned_pairs(matches_only=True)
    if mode == 0:
        for q_pos, r_pos in aligned_pair_iter:
            if r_pos == ref_pos:
                return q_pos
    elif mode == 1:
        for q_pos, r_pos in aligned_pair_iter:
            if r_pos >= ref_pos:
                return q_pos
    elif mode == 2:
        prev_r_pos = None
        for q_pos, r_pos in aligned_pair_iter:
            if r_pos >= ref_pos:
                if prev_r_pos is None:
                    error("prev_r_pos is None")
                return prev_r_pos
            prev_r_pos = r_pos
    else:
        error("Invalid mode: %s" % mode)
    error("Reference position did not map to query position: %d\n%s" % (ref_pos, str(s.get_aligned_pairs(matches_only=True))))

# trim a single read
def trim_read(s, overlapping_primers, primers, max_primer_len):
    '''Primer-trim and quality-trim an individual read. In the future, for multiprocessing, I can just change ``s`` argument to be the relevant member variables (e.g. reference_start, reference_end, query_length, is_paired, is_reverse, etc.) and just return the final updated pysam cigartuples (list of tuples)

    Args:
        ``s`` (``pysam.AlignedSegment``): The mapped read to trim (see: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment)
    '''
    # get list of primers that cover the start and end indices (wrt reference) of this alignment (each element in each list is an index in `primers`)
    overlapping_primer_inds_start = overlapping_primers[s.reference_start]
    overlapping_primer_inds_end = overlapping_primers[s.reference_end-1] # "reference_end points to one past the last aligned residue"
    isize_flag = (abs(s.template_length) - max_primer_len) > s.query_length # Why does iVar compare this read's isize (template_length) against the max length of ALL primers (max_primer_len) instead of the max length of primers that cover this read's start position?

    # trim forward primer (paired reads: just forward; unpaired reads: all)
    if len(overlapping_primer_inds_start) > 0 and ((s.is_paired and not s.is_reverse) or not s.is_paired):
        # prepare to trim forward primer
        max_overlap_end = max(primers[i][1] for i in overlapping_primer_inds_start) + 1 # just after max end of overlapping primers
        delete_start = ref_pos_to_query_pos(s, max_overlap_end, mode=1)
        new_cigar = list(); ref_add = 0; del_len = delete_start; pos_start = False; start_pos = 0

        # for each operation in the CIGAR
        for cig, n in s.cigartuples:
            # if nothing left to delete, just append everything
            if del_len == 0 and (pos_start or (CONSUME_QUERY[cig] and CONSUME_REF[cig])):
                pos_start = True; new_cigar.append((cig,n)); continue

            # if our operation consumes the query"
            ref_add = 0
            if CONSUME_QUERY[cig]:
                # How much do we have to delete?
                if del_len >= n: # our entire operation needs to be deleted
                    new_cigar.append((BAM_CSOFT_CLIP,n))
                elif 0 < del_len < n: # We need to delete SOME of our segment, but we still have more later
                    new_cigar.append((BAM_CSOFT_CLIP,del_len))
                elif del_len == 0: # Just keep clipping
                    new_cigar.append((BAM_CSOFT_CLIP,n)); continue

                # Update based on how much we just deleted
                ref_add = min(del_len, n)
                tmp = n
                n = max(n - del_len, 0)
                del_len = max(del_len - tmp, 0)

                # If there is still more left to do, append it
                if n > 0:
                    new_cigar.append((cig,n))

                # If we are done and just consumed, we want to just start appending everything
                if del_len == 0 and CONSUME_QUERY[new_cigar[-1][0]] and CONSUME_REF[new_cigar[-1][0]]:
                    pos_start = True

            # If our trim consumed the reference, we need to move our read's start position forward
            if CONSUME_REF[cig]:
                start_pos += ref_add

        # update our alignment accordingly
        s.cigartuples = new_cigar
        s.reference_start += start_pos

    # trim reverse primer (paired reads: just reverse; unpaired reads: all)
    if len(overlapping_primer_inds_end) > 0 and ((s.is_paired and s.is_reverse) or not s.is_paired):
        # prepare to trim reverse primer
        min_overlap_start = min(primers[i][0] for i in overlapping_primer_inds_end) - 1 # just before min start of overlapping primers
        delete_end = ref_pos_to_query_pos(s, min_overlap_start, mode=2)
        new_cigar = list(); ref_add = 0; del_len = delete_end; pos_start = False

        # for each operation in the CIGAR (reverse order because reverse primer)
        for cig, n in reversed(s.cigartuples):
            # if nothing left to delete, just append everything
            if del_len == 0 and (pos_start or (CONSUME_QUERY[cig] and CONSUME_REF[cig])):
                pos_start = True; new_cigar.append((cig,n)); continue

            # TODO FINISH REVERSE PRIMER

        # update our alignment accordingly
        s.cigartuples = list(reversed(new_cigar)) # I appended to new_cigar backwards, so it needs to be reversed at the end

# run AmpliPy Trim
def run_trim(untrimmed_reads_fn, primer_fn, reference_fn, trimmed_reads_fn, primer_pos_offset, min_length, min_quality, sliding_window_width, include_no_primer):
    '''Run AmpliPy Trim

    * This is where iVar's code does iVar Trim:  https://github.com/TheCrossBoy/ivar/blob/2829ebf359fab633bac5e9f4d19c9c42746629fd/src/trim_primer_quality.cpp#L559-L928
    * This is where it does the actual trimming: https://github.com/TheCrossBoy/ivar/blob/2829ebf359fab633bac5e9f4d19c9c42746629fd/src/trim_primer_quality.cpp#L696-L857


    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of input untrimmed reads SAM/BAM

        ``primer_fn`` (``str``): Filename of input primer BED

        ``reference_fn`` (``str``): Filename of input reference genome FASTA (is this needed?)

        ``trimmed_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM

        ``primer_pos_offset`` (``int``): Primer position offset. Reads that occur at the specified offset positions relative to primer positions will also be trimmed

        ``min_length`` (``int``): Minimum length of read to retain after trimming

        ``min_quality`` (``int``): Minimum quality threshold for sliding window to pass

        ``sliding_window_width`` (``int``): Width of sliding window

        ``include_no_primer`` (``bool``): ``True`` to include reads with no primers, otherwise ``False``
    '''
    # load input files and preprocess
    print_log("Executing AmpliPy Trim (v%s)" % VERSION)
    print_log("Loading reference genome: %s" % reference_fn)
    ref_genome_ID, ref_genome_sequence = load_ref_genome(reference_fn)
    print_log("Loading primers: %s" % primer_fn)
    primers = load_primers(primer_fn)
    max_primer_len = max(end-start for start,end in primers)
    print_log("Precalculate overlapping primers...")
    overlapping_primers = find_overlapping_primers(len(ref_genome_sequence), primers)
    print_log("Input untrimmed SAM/BAM: %s" % untrimmed_reads_fn)
    print_log("Output untrimmed SAM/BAM: %s" % trimmed_reads_fn)
    in_aln, out_aln = create_AlignmentFile_objects(untrimmed_reads_fn, trimmed_reads_fn)

    # initialize counters
    NUM_UNMAPPED = 0
    NUM_NO_CIGAR = 0

    # trim reads
    print_log("Trimming reads...")
    s_i = 0
    for s in in_aln:
        # skip unmapped reads
        if s.is_unmapped:
            NUM_UNMAPPED += 1; continue

        # skip reads without CIGAR
        if s.cigartuples is None:
            NUM_NO_CIGAR += 1; continue

        # trim this read and write the result
        status_trim_primers = trim_read(s, overlapping_primers, primers, max_primer_len)
        out_aln.write(s)

        # print progress update
        s_i += 1
        if s_i % PROGRESS_NUM_READS == 0:
            print_log("Trimmed %d reads..." % s_i)
    print_log("Finished trimming %d reads" % s_i)

	# TODO DELETE WHEN DONE
    error("TRIM NOT IMPLEMENTED\n- untrimmed_reads_fn: %s\n- primer_fn: %s\n- trimmed_reads_fn: %s" % (untrimmed_reads_fn, primer_fn, trimmed_reads_fn)) # TODO

# run AmpliPy Variants
def run_variants(trimmed_reads_fn, variants_fn):
    '''Run AmpliPy Variants

    Args:
        ``trimmed_reads_fn`` (``str``): Filename of input trimmed reads SAM/BAM

        ``reference_fn`` (``str``): Filename of input reference genome FASTA (is this needed?)

        ``variants_fn`` (``str``): Filename of output variants VCF
    '''
    print_log("Executing AmpliPy Variants (v%s)" % VERSION)
    error("VARIANTS NOT IMPLEMENTED\n- trimmed_reads_fn: %s\n- variants_fn: %s" % (trimmed_reads_fn, variants_fn)) # TODO

# run AmpliPy Consensus
def run_consensus(trimmed_reads_fn, consensus_fn):
    '''Run AmpliPy Consensus

    Args:
        ``trimmed_reads_fn`` (``str``): Filename of input trimmed reads SAM/BAM

        ``reference_fn`` (``str``): Filename of input reference genome FASTA (is this needed?)

        ``consensus_fn`` (``str``): Filename of output consensus sequence FASTA
    '''
    print_log("Executing AmpliPy Consensus (v%s)" % VERSION)
    error("CONSENSUS NOT IMPLEMENTED\n- trimmed_reads_fn: %s\n- consensus_fn: %s" % (trimmed_reads_fn, consensus_fn)) # TODO

# run AmpliPy AIO (All-In-One)
def run_aio(untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn, variants_fn, consensus_fn):
    '''Run AmpliPy AIO (All-In-One)

    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of input untrimmed reads SAM/BAM

        ``primer_fn`` (``str``): Filename of input primer BED

        ``reference_fn`` (``str``): Filename of input reference genome FASTA (is this needed?)

        ``trimmed_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM

        ``variants_fn`` (``str``): Filename of output variants VCF

        ``consensus_fn`` (``str``): Filename of output consensus sequence FASTA
    '''
    print_log("Executing AmpliPy AIO (v%s)" % VERSION)
    error("AIO NOT IMPLEMENTED\n- untrimmed_reads_fn: %s\n- amplipy_index_fn: %s\n- trimmed_reads_fn: %s\n- variants_fn: %s\n- consensus_fn: %s" % (untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn, variants_fn, consensus_fn)) # TODO

# main content
if __name__ == "__main__":
    if len(argv) == 1:
        pass # TODO: In the future, run GUI here to fill in argv accordingly (so argparse will run fine)
    args = parse_args()
    if args.command == 'trim':
        run_trim(args.input, args.primer, args.reference, args.output, args.primer_pos_offset, args.min_length, args.min_quality, args.sliding_window_width, args.include_no_primer)
    elif args.command == 'variants':
        run_variants(args.input, args.output)
    elif args.command == 'consensus':
        run_consensus(args.input, args.output)
    elif args.command == 'aio':
        run_aio(args.input, args.amplipy_index, args.output_trimmed_reads, args.output_variants, args.output_consensus)

# OLD STUFF BELOW #
exit() # TODO THIS JUST STOPS EXECUTION BEFORE WE GO LOWER
##### Argument Setup #####
parser = argparse.ArgumentParser()
parser.add_argument("-d", dest="debug", help="Enable debug mode for extra output",
                    action="store_true")  # Currently unused
subparsers = parser.add_subparsers(dest="dest")

# Trim args
trim_parser = subparsers.add_parser("trim", add_help=False)

# Help
trim_help = trim_parser.add_argument_group("help")
trim_help.add_argument("-h", "--help", action="help",
                       help="show this help message and exit")

# Input args
trim_input_args = trim_parser.add_argument_group("required input args")
trim_input_args.add_argument(
    "-i", dest="input_file", required=True, help="Bam file with aligned reads")
trim_input_args.add_argument(
    "-b", dest="primer_file", required=True, help="BED file with primer positions")

# Output args
trim_output_args = trim_parser.add_argument_group("required output args")
trim_output_args.add_argument(
    "-p", dest="out_file", required=True, help="Output file for BAM file")

# Optargs
trim_opt_args = trim_parser.add_argument_group("optional args")
trim_opt_args.add_argument("-m", dest="min_read_len",
                           help="Minimum required length for reads in order to be added to output", type=int, default=30)
trim_opt_args.add_argument("-q", dest="qual_thresh",
                           help="Minimum quality threshold for sliding window to meet in order to pass quality test", type=int, default=20)
trim_opt_args.add_argument("-s", dest="sliding_window",
                           help="Base length of sliding window", type=int, default=4)
trim_opt_args.add_argument("-e", dest="incl_no_primer",
                           help="Include reads with no primers", action="store_true")

# Variant args (todo)
variant_parser = subparsers.add_parser("variant")

args = parser.parse_args()
if (args.dest == "trim"):
    ##### Trim #####
    # Process arguments & verify paths
    if not isfile(args.input_file):
        # No input file
        exit("error: input file %s does not exist" % args.input_file)
    if not isfile(args.primer_file):
        # No primer file
        exit("error: primer file %s does not exist" % args.input_file)

    # Open input & output files
    bam = pysam.AlignmentFile(args.input_file, "r")
    primer_file = open(args.primer_file)
    output = pysam.AlignmentFile(args.out_file, "w", header=bam.header)

    # Call trim
    stats = trim(bam, primer_file, output, min_read_length=args.min_read_len, min_qual_thresh=args.qual_thresh,
                 sliding_window=args.sliding_window, include_no_primer=args.incl_no_primer)

    # Print stats
    print("%i reads were primer trimmed" % stats["primer_trimmed_count"])
    print("%i reads didn't include a primer%s" % (
        stats["no_primer_count"], " (but were included because -e was specified)" if args.incl_no_primer else ""))
    print("%i reads were quality trimmed" % stats["quality_trimmed_count"])
    print("%i reads were excluded due to being shorter than the minimum read length (%i)" % (
        stats["removed_reads"], args.min_read_len))
elif (args.dest == "variant"):
    ##### Call Variant #####
    print("Nothing here yet!")
else:
    # This should never happen! Argparse will handle this, but just in case
    print("This should never happen. Please let us know what you did to get here!")
