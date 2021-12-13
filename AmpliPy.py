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
def find_overlapping_primers(ref_genome_length, primers, primer_pos_offset):
    '''Find all primers that cover every position of the reference genome

    Args:
        ``ref_genome_length`` (``int``): Length of the reference genome

        ``primers`` (``list`` of ``(int,int)``): List of primer indices, where each primer is represented as a ``(start,end)`` tuple. Note that this uses 0-based indices, ``start`` is **in**clusive, and ``end`` is **ex**clusive. For example, `(0,100)` is the 100-length window starting at index 0 and ending at index 99.

        ``primer_pos_offset`` (``int``): Primer position offset. Reads that occur at the specified offset positions relative to primer positions will also be trimmed

    Returns:
        ``list`` of ``int``: The minimum start position (inclusive) of all primers that cover each position of the reference genome

        ``list`` of ``int``: The maximum end position (exclusive) of all primers that cover each position of the reference genome
    '''
    # set things up
    min_primer_start = [None for _ in range(ref_genome_length)]
    max_primer_end = [None for _ in range(ref_genome_length)]
    num_primers = len(primers); curr_primers = deque(); i = 0 # i = current index of primers; primers[i] is the "current" (start,end) tuple

    # compute overlapping primers for each position of the reference genome
    for ref_pos in range(ref_genome_length):
        # remove primers that have now been passed
        while len(curr_primers) != 0 and ref_pos >= (primers[curr_primers[0]][1] + primer_pos_offset):
            curr_primers.popleft()

        # add primers that have now been entered
        while i < num_primers and ref_pos >= (primers[i][0] - primer_pos_offset):
            curr_primers.append(i); i += 1

        # the primers in curr_primers must span ref_pos
        if len(curr_primers) != 0:
            min_primer_start[ref_pos] = min(primers[i][0] for i in curr_primers)
            max_primer_end[ref_pos] = max(primers[i][1] for i in curr_primers)
        ref_pos += 1
    return min_primer_start, max_primer_end

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
        ``list`` of ``tuple``: The primers in the given BED, represented as ``(start, end)`` tuples (start = inclusive, end = exclusive, just like in BED)
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

# NOTE TO NIEMA: I think these get_pos_on_* functions can be sped up by using functions in pysam.AlignedSegment

# get query position on reference
def get_pos_on_ref(cigar, query_pos, ref_start):
    '''Convert a query position to a reference position

    Args:
        ``cigar`` (``list`` of ``tuple``): The CIGAR operations of the alignment

        ``query_pos`` (``int``): The query position

        ``ref_start`` (``int``): The reference start position

    Returns:
        ``int``: The corresponding reference position
    '''
    cur_pos = 0; ref_pos = ref_start
    for cig, n in cigar:
        if CONSUME_QUERY[cig]:
            if query_pos <= cur_pos + n:
                if CONSUME_REF[cig]:
                    ref_pos += (query_pos - cur_pos)
                return ref_pos
            cur_pos += n
        if CONSUME_REF[cig]:
            ref_pos += n
    return ref_pos

# get reference position on query
def get_pos_on_query(cigar, ref_pos, ref_start):
    '''Convert a reference position to a query position

    Args:
        ``cigar`` (``list`` of ``tuple``): The CIGAR operations of the alignment

        ``ref_pos`` (``int``): The reference position

        ``query_start`` (``int``): The query start position

    Returns:
        ``int``: The corresponding query position
    '''
    query_pos = 0; cur_pos = ref_start
    for cig, n in cigar:
        if CONSUME_REF[cig]:
            if ref_pos <= cur_pos + n:
                if CONSUME_QUERY[cig]:
                    query_pos += (ref_pos - cur_pos)
                return query_pos
            cur_pos += n
        if CONSUME_QUERY[cig]:
            query_pos += n
    return query_pos

# fix new CIGAR after trimming
def fix_cigar(new_cigar):
    if not isinstance(new_cigar, list):
        new_cigar = list(new_cigar)
    proper_cigar = list()
    for i in range(len(new_cigar)):
        if i < len(new_cigar)-1 and new_cigar[i][0] == new_cigar[i+1][0]:
            new_cigar[i+1] = (new_cigar[i+1][0], new_cigar[i][1] + new_cigar[i+1][1]); continue
        proper_cigar.append(new_cigar[i])
    return proper_cigar

# trim a single read
def trim_read(s, min_primer_start, max_primer_end, max_primer_len, min_quality, sliding_window_width):
    '''Primer-trim and quality-trim an individual read. In the future, for multiprocessing, I can just change ``s`` argument to be the relevant member variables (e.g. reference_start, reference_end, query_length, is_paired, is_reverse, etc.) and just return the final updated pysam cigartuples (list of tuples)

    Args:
        ``s`` (``pysam.AlignedSegment``): The mapped read to trim (see: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment)

        ``min_primer_start`` (``list`` of ``int``): The minimum start position (inclusive) of all primers that cover each position of the reference genome

        ``max_primer_end`` (``list`` of ``int``): The maximum end position (exclusive) of all primers that cover each position of the reference genome

        ``max_primer_len`` (``int``): Length of longest primer in ``primers``

        ``min_quality`` (``int``): Minimum quality threshold for sliding window to pass

        ``sliding_window_width`` (``int``): Width of sliding window

    Returns:
        ``bool``: ``True`` if primer trimming was done at the start of the read, otherwise ``False``

        ``bool``: ``True`` if primer trimming was done at the end of the read, otherwise ``False``

        ``bool``: ``True`` if quality trimming was done, otherwise ``False``
    '''
    # get list of primers that cover the start and end indices (wrt reference) of this alignment (each element in each list is an index in `primers`)
    left_max_primer_end = max_primer_end[s.reference_start]
    right_min_primer_start = min_primer_start[s.reference_end-1] # "reference_end points to one past the last aligned residue"
    isize_flag = (abs(s.template_length) - max_primer_len) > s.query_length # Why does iVar compare this read's isize (template_length) against the max length of ALL primers (max_primer_len) instead of the max length of primers that cover this read's start position?

    # flags for what happened when trimming this read
    trimmed_primer_start = False
    trimmed_primer_end = False
    trimmed_quality = False

    # trim primer from start of read
    if not (s.is_paired and isize_flag and s.is_reverse) and left_max_primer_end is not None:
        # prepare to trim forward primer
        trimmed_primer_start = True
        delete_start = get_pos_on_query(s.cigartuples, left_max_primer_end + 1, s.reference_start)
        new_cigar = list(); ref_add = 0; del_len = delete_start; pos_start = False; start_pos = 0

        # for each operation in the CIGAR
        for operation in s.cigartuples:
            # If we have nothing left to delete, just append everything; otherwise, separate tuple for convenience
            if del_len == 0 and pos_start:
                new_cigar.append(operation); continue
            cig, n = operation
            
            # If we have nothing left to delete and are consuming both, we want to just append everything
            if del_len == 0 and CONSUME_QUERY[cig] and CONSUME_REF[cig]:
                pos_start = True; new_cigar.append(operation); continue

            # How much our current trim affects our read's start position
            ref_add = 0

            # If our operation consumes the query
            if CONSUME_QUERY[cig]:
                # How much do we have to delete?
                if del_len >= n: # Our entire operation needs to be deleted
                    new_cigar.append((BAM_CSOFT_CLIP, n))
                elif 0 < del_len < n: # We need to delete some of our segment, but we will still have more later
                    new_cigar.append((BAM_CSOFT_CLIP, del_len))
                else: # Since we consume the query, we just need to keep clipping
                    new_cigar.append((BAM_CSOFT_CLIP, n)); continue
                
                # Update based on how much we just deleted
                ref_add = min(del_len, n)
                tmp = n
                n = max(n - del_len, 0)
                del_len = max(del_len - tmp, 0)

                # If there is still more left to do, append it
                if n > 0:
                    new_cigar.append((cig, n))

                # If we are done and just consumed, we want to just start appending everything.
                if del_len == 0 and CONSUME_QUERY[new_cigar[-1][0]] and CONSUME_REF[new_cigar[-1][0]]:
                    pos_start = True

            # If our trim consumed the reference, we need to move our read's start position forwards
            if CONSUME_REF[cig]:
                start_pos += ref_add

        # update our alignment accordingly
        s.cigartuples = fix_cigar(new_cigar)
        s.reference_start += start_pos

    # trim primer from end of read
    if not (s.is_paired and isize_flag and not s.is_reverse) and right_min_primer_start is not None:
        # prepare to trim reverse primer
        trimmed_primer_end = True
        delete_end = s.query_length - get_pos_on_query(s.cigartuples, right_min_primer_start, s.reference_start)
        new_cigar = list(); ref_add = 0; del_len = delete_end; pos_start = False

        # for each operation in the CIGAR (reverse order because reverse primer)
        for operation in reversed(s.cigartuples):
            # if nothing left to delete, just append everything
            if del_len == 0 and pos_start:
                new_cigar.append((operation)); continue
            cig, n = operation
            
            # If we have nothing left to delete and are consuming both, we want to just append everything
            if del_len == 0 and CONSUME_QUERY[cig] and CONSUME_REF[cig]:
                pos_start = True; new_cigar.append(operation); continue

            # If our operation consumes the query
            if CONSUME_QUERY[cig]:
                # How much do we have to delete?
                if del_len >= n: # Our entire operation needs to be deleted
                    new_cigar.append((BAM_CSOFT_CLIP, n))
                elif 0 < del_len < n: # We need to delete some of our segment, but we will still have more later
                    new_cigar.append((BAM_CSOFT_CLIP, del_len))
                else: # Since we consume the query, we just need to keep clipping
                    new_cigar.append((BAM_CSOFT_CLIP, n)); continue

                # Update based on how much we just deleted
                tmp = n
                n = max(n - del_len, 0)
                del_len = max(del_len - tmp, 0)

                # If there is still more left to do, append it
                if n > 0:
                    new_cigar.append((cig, n))

                # If we are done and just consumed, we want to just start appending everything.
                if del_len == 0 and CONSUME_QUERY[new_cigar[-1][0]] and CONSUME_REF[new_cigar[-1][0]]:
                    pos_start = True

        # update our alignment accordingly
        s.cigartuples = fix_cigar(reversed(new_cigar)) # I appended to new_cigar backwards, so it needs to be reversed at the end

    # set things up for quality trimming
    qual = s.query_alignment_qualities
    total = 0; true_start = 0; true_end = len(qual)
    window = min(sliding_window_width, len(qual))

    # quality trim (reverse strand, so trim from beginning)
    if s.is_reverse:
        # build our buffer
        i = true_end
        for offset in range(1, window):
            total += qual[i - offset]

        # Loop through the read, determine when we need to trim
        while i > true_start:
            # We are nearing the end, so we need to shrink our window
            if true_start + window > i:
                window -= 1

            # Still have more to go, so add in our new value
            else:
                total += qual[i - window]

            # Check our current quality score
            if (total / window) < min_quality:
                break

            # Remove the no longer needed quality score, and advance i
            total -= qual[i - 1]; i -= 1

        # Initialization for quality trimming
        new_cigar = list(); del_len = i
        start_pos = get_pos_on_ref(s.cigartuples, del_len + s.query_alignment_start, s.reference_start)

        # Do we need to trim?
        if start_pos > s.reference_start:
            trimmed_quality = True
            # for each operation in the CIGAR
            for operation in s.cigartuples:
                # Nothing left to trim, just append everything
                if del_len == 0:
                    new_cigar.append(operation); continue
                cig, n = operation

                # These are just clips, so it's not part of our quality trim determination
                if cig == BAM_CSOFT_CLIP or cig == BAM_CHARD_CLIP:
                    new_cigar.append(operation); continue

                # We consume the query, so we may need to trim
                if CONSUME_QUERY[cig]:
                    # How much do we need to delete?
                    if del_len >= n: # All of the current operation
                        new_cigar.append((BAM_CSOFT_CLIP, n))
                    else: # Only part of the current operation
                        new_cigar.append((BAM_CSOFT_CLIP, del_len))

                    # Decrease our delete length by how much we've deleted
                    tmp = n
                    n = max(n - del_len, 0)
                    del_len = max(del_len - tmp, 0)

                    # If we ran out of things to delete, we need to append the rest of the operation
                    if n > 0:
                        new_cigar.append((cig, n))

            # update our alignment accordingly
            s.cigartuples = fix_cigar(new_cigar)

    # quality trim (forward strand, so trim from end)
    else:
        # build our buffer
        i = true_start
        for offset in range(window - 1):
            total += qual[i + offset]

        # Loop through the read, determine when we need to trim
        while i < true_end:
            # We are nearing the end, so we need to shrink our window
            if (true_end - window) < i:
                window -= 1

            # Still have more to go, so add in our new value
            else:
                total += qual[i + window - 1]

            # Check our current quality score
            if (total / window) < min_quality:
                break

            # Remove the no longer needed quality score, and advance i
            total -= qual[i]; i += 1

        # Initialization for quality trimming
        new_cigar = list(); del_len = true_end - i
        start_pos = get_pos_on_ref(s.cigartuples, del_len, s.reference_start)

        # Do we need to trim?
        if del_len != 0:
            trimmed_quality = True
            for operation in reversed(s.cigartuples):
                # Nothing left to trim, just append everything
                if del_len == 0:
                    new_cigar.append(operation); continue
                cig, n = operation

                # These are just clips, so it's not part of our quality trim determination
                if cig == BAM_CSOFT_CLIP or cig == BAM_CHARD_CLIP:
                    new_cigar.append(operation); continue

                # We consume the query, so we may need to trim
                if CONSUME_QUERY[cig]:
                    # How much do we need to delete?
                    if del_len >= n: # All of the current operation
                        new_cigar.append((BAM_CSOFT_CLIP, n))
                    else: # Only part of the current operation
                        new_cigar.append((BAM_CSOFT_CLIP, del_len))
                        
                    # Decrease our delete length by how much we've deleted
                    tmp = n
                    n = max(n - del_len, 0)
                    del_len = max(del_len - tmp, 0)

                    # If we ran out of things to delete, we need to append the rest of the operation
                    if n > 0:
                        new_cigar.append((cig, n))

            # update our alignment accordingly
            s.cigartuples = fix_cigar(reversed(new_cigar)) # I appended to new_cigar backwards, so it needs to be reversed at the end
    return trimmed_primer_start, trimmed_primer_end, trimmed_quality

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
    min_primer_start, max_primer_end = find_overlapping_primers(len(ref_genome_sequence), primers, primer_pos_offset)
    print_log("Input untrimmed SAM/BAM: %s" % untrimmed_reads_fn)
    print_log("Output untrimmed SAM/BAM: %s" % trimmed_reads_fn)
    in_aln, out_aln = create_AlignmentFile_objects(untrimmed_reads_fn, trimmed_reads_fn)

    # initialize counters
    NUM_UNMAPPED = 0
    NUM_NO_CIGAR = 0
    NUM_TRIMMED_PRIMER_START = 0
    NUM_TRIMMED_PRIMER_END = 0
    NUM_UNTRIMMED_PRIMER = 0
    NUM_TRIMMED_QUALITY = 0
    NUM_TOO_SHORT = 0
    NUM_WRITTEN = 0

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

        # trim this read
        trimmed_primer_start, trimmed_primer_end, trimmed_quality = trim_read(s, min_primer_start, max_primer_end, max_primer_len, min_quality, sliding_window_width)
        if trimmed_primer_start:
            NUM_TRIMMED_PRIMER_START += 1
        if trimmed_primer_end:
            NUM_TRIMMED_PRIMER_END += 1
        if trimmed_quality:
            NUM_TRIMMED_QUALITY += 1

        # write this read (if applicable)
        write_read = True
        if s.reference_length < min_length:
            NUM_TOO_SHORT += 1; write_read = False
        if not (trimmed_primer_start or trimmed_primer_end):
            NUM_UNTRIMMED_PRIMER += 1
            if not include_no_primer:
                write_read = False
        if write_read:
            out_aln.write(s); NUM_WRITTEN += 1

        # print progress update
        s_i += 1
        if s_i % PROGRESS_NUM_READS == 0:
            print_log("Trimmed %d reads..." % s_i)
    print_log("Finished trimming %d reads" % s_i)
    print_log("- Number of Unmapped Reads: %d" % NUM_UNMAPPED)
    print_log("- Number of Mapped Reads Without CIGAR: %d" % NUM_NO_CIGAR)
    print_log("- Number of Mapped Reads With Primer Trimmed at Start: %d" % NUM_TRIMMED_PRIMER_START)
    print_log("- Number of Mapped Reads With Primer Trimmed at End: %d" % NUM_TRIMMED_PRIMER_END)
    print_log("- Number of Mapped Reads With No Primers Trimmed: %d" % NUM_UNTRIMMED_PRIMER)
    print_log("- Number of Mapped Reads That Were Quality Trimmed: %d" % NUM_TRIMMED_QUALITY)
    print_log("- Number of Mapped Reads Written to Output: %d" % NUM_WRITTEN)

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
