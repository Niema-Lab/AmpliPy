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
VERSION = '0.0.2'
BUFSIZE = 1048576 # 1 MB
PROGRESS_NUM_READS = 50000

# default arguments
DEFAULT_MIN_DEPTH_CONSENSUS = 10
DEFAULT_MIN_DEPTH_VARIANTS = 1
DEFAULT_MIN_FREQ_CONSENSUS = 0
DEFAULT_MIN_FREQ_VARIANTS = 0.03
DEFAULT_MIN_LENGTH = 30
DEFAULT_MIN_QUALITY = 20
DEFAULT_PRIMER_POS_OFFSET = 0
DEFAULT_SLIDING_WINDOW_WIDTH = 4
DEFAULT_UNKNOWN_SYMBOL = 'N'

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
ERROR_TEXT_INVALID_MIN_DEPTH = "Minimum depth must be positive"
ERROR_TEXT_INVALID_MIN_FREQ = "Minimum frequency must be between 0 and 1"
ERROR_TEXT_INVALID_MIN_LENGTH = "Minimum length must be >= 1"
ERROR_TEXT_INVALID_READ_EXTENSION = "Invalid read mapping extension (should be .sam or .bam)"
ERROR_TEXT_INVALID_SLIDING_WINDOW_WIDTH = "Sliding window width must be >= 1"
ERROR_TEXT_INVALID_UNKNOWN_SYMBOL_LENGTH = "Unknown symbol must be exactly 1 character"
ERROR_TEXT_INVALID_VCF_EXTENSION = "Invalid variants extension (should be .vcf, .vcf.gz, or .bcf)"
ERROR_TEXT_MULTIPLE_REF_SEQS = "Multiple sequences in FASTA file"
ERROR_TEXT_NEGATIVE_MIN_QUALITY = "Minimum quality must be non-negative"
ERROR_TEXT_NEGATIVE_PRIMER_POS_OFFSET = "Primer position offset must be non-negative"
HELP_TEXT_AMPLIPY_INDEX = "AmpliPy Index (PKL)"
HELP_TEXT_CONSENSUS = "Consensus Sequence (FASTA)"
HELP_TEXT_MIN_DEPTH_CONSENSUS = "Minimum depth to call consensus"
HELP_TEXT_MIN_DEPTH_VARIANTS = "Minimum depth to call variant"
HELP_TEXT_MIN_FREQ_CONSENSUS = "Minimum frequency threshold (0-1) to call consensus"
HELP_TEXT_MIN_FREQ_VARIANTS = "Minimum frequency threshold (0-1) to call variant"
HELP_TEXT_MIN_QUAL = "Minimum quality threshold"
HELP_TEXT_PRIMER = "Primer File (BED)"
HELP_TEXT_READS_UNTRIMMED = "Untrimmed Reads (SAM/BAM)"
HELP_TEXT_READS_TRIMMED = "Trimmed Reads (SAM/BAM)"
HELP_TEXT_REFERENCE = "Reference Genome (FASTA)"
HELP_TEXT_TRIM_INCLUDE_READS_NO_PRIMER = "Include reads with no primers"
HELP_TEXT_TRIM_MIN_LENGTH = "Minimum length of read to retain after trimming"
HELP_TEXT_TRIM_PRIMER_POS_OFFSET = "Primer position offset. Reads that occur at the specified offset positions relative to primer positions will also be trimmed"
HELP_TEXT_TRIM_SLIDING_WINDOW_WIDTH = "Width of sliding window (average quality of this window must be >= minimum quality threshold)"
HELP_TEXT_UNKNOWN_SYMBOL = "Character to print in regions with less than minimum coverage"
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

# get alignment (for debugging)
def get_alignment(s, ref_genome_sequence=None):
    q = ''; r = ''
    for q_pos, r_pos in s.get_aligned_pairs():
        if q_pos is not None and q_pos < s.query_alignment_start:
            continue
        if q_pos is not None and q_pos >= s.query_alignment_end:
            break
        if q_pos is None:
            q += '-'
        else:
            q += s.query_sequence[q_pos]
        if r_pos is None:
            r += '-'
        elif ref_genome_sequence is None:
            r += '?'
        else:
            r += ref_genome_sequence[r_pos]
    return q, r

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
    trim_parser.add_argument('-x', '--primer_pos_offset', required=False, type=int, default=DEFAULT_PRIMER_POS_OFFSET, help=HELP_TEXT_TRIM_PRIMER_POS_OFFSET)
    trim_parser.add_argument('-ml', '--min_length', required=False, type=int, default=DEFAULT_MIN_LENGTH, help=HELP_TEXT_TRIM_MIN_LENGTH)
    trim_parser.add_argument('-mq', '--min_quality', required=False, type=int, default=DEFAULT_MIN_QUALITY, help=HELP_TEXT_MIN_QUAL)
    trim_parser.add_argument('-s', '--sliding_window_width', required=False, type=int, default=DEFAULT_SLIDING_WINDOW_WIDTH, help=HELP_TEXT_TRIM_SLIDING_WINDOW_WIDTH)
    trim_parser.add_argument('-e', '--include_no_primer', action='store_true', help=HELP_TEXT_TRIM_INCLUDE_READS_NO_PRIMER)

    # AmpliPy Variants args
    variants_parser = subparsers.add_parser("variants", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    variants_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    variants_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    variants_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_VARIANTS)
    variants_parser.add_argument('-mq', '--min_quality', required=False, type=int, default=DEFAULT_MIN_QUALITY, help=HELP_TEXT_MIN_QUAL)
    variants_parser.add_argument('-mf', '--min_freq', required=False, type=float, default=DEFAULT_MIN_FREQ_VARIANTS, help=HELP_TEXT_MIN_FREQ_VARIANTS)
    variants_parser.add_argument('-md', '--min_depth', required=False, type=int, default=DEFAULT_MIN_DEPTH_VARIANTS, help=HELP_TEXT_MIN_DEPTH_VARIANTS)

    # AmpliPy Consensus args
    consensus_parser = subparsers.add_parser("consensus", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    consensus_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    consensus_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    consensus_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_CONSENSUS)
    consensus_parser.add_argument('-mq', '--min_quality', required=False, type=int, default=DEFAULT_MIN_QUALITY, help=HELP_TEXT_MIN_QUAL)
    consensus_parser.add_argument('-mf', '--min_freq', required=False, type=float, default=DEFAULT_MIN_FREQ_CONSENSUS, help=HELP_TEXT_MIN_FREQ_CONSENSUS)
    consensus_parser.add_argument('-md', '--min_depth', required=False, type=int, default=DEFAULT_MIN_DEPTH_CONSENSUS, help=HELP_TEXT_MIN_DEPTH_CONSENSUS)
    consensus_parser.add_argument('-n', '--unknown_symbol', required=False, type=str, default=DEFAULT_UNKNOWN_SYMBOL, help=HELP_TEXT_UNKNOWN_SYMBOL)

    # AmpliPy AIO (All-In-One) args
    aio_parser = subparsers.add_parser("aio", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aio_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_UNTRIMMED)
    aio_parser.add_argument('-p', '--primer', required=True, type=str, help=HELP_TEXT_PRIMER)
    aio_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    aio_parser.add_argument('-ot', '--output_trimmed_reads', required=True, type=str, help=HELP_TEXT_READS_TRIMMED)
    aio_parser.add_argument('-ov', '--output_variants', required=True, type=str, help=HELP_TEXT_VARIANTS)
    aio_parser.add_argument('-oc', '--output_consensus', required=True, type=str, help=HELP_TEXT_CONSENSUS)
    aio_parser.add_argument('-x', '--primer_pos_offset', required=False, type=int, default=DEFAULT_PRIMER_POS_OFFSET, help=HELP_TEXT_TRIM_PRIMER_POS_OFFSET)
    aio_parser.add_argument('-ml', '--min_length', required=False, type=int, default=DEFAULT_MIN_LENGTH, help=HELP_TEXT_TRIM_MIN_LENGTH)
    aio_parser.add_argument('-mq', '--min_quality', required=False, type=int, default=DEFAULT_MIN_QUALITY, help=HELP_TEXT_MIN_QUAL)
    aio_parser.add_argument('-s', '--sliding_window_width', required=False, type=int, default=DEFAULT_SLIDING_WINDOW_WIDTH, help=HELP_TEXT_TRIM_SLIDING_WINDOW_WIDTH)
    aio_parser.add_argument('-mfc', '--min_freq_consensus', required=False, type=float, default=DEFAULT_MIN_FREQ_CONSENSUS, help=HELP_TEXT_MIN_FREQ_CONSENSUS)
    aio_parser.add_argument('-mfv', '--min_freq_variants', required=False, type=float, default=DEFAULT_MIN_FREQ_VARIANTS, help=HELP_TEXT_MIN_FREQ_VARIANTS)
    aio_parser.add_argument('-mdc', '--min_depth_consensus', required=False, type=int, default=DEFAULT_MIN_DEPTH_CONSENSUS, help=HELP_TEXT_MIN_DEPTH_CONSENSUS)
    aio_parser.add_argument('-mdv', '--min_depth_variants', required=False, type=int, default=DEFAULT_MIN_DEPTH_VARIANTS, help=HELP_TEXT_MIN_DEPTH_VARIANTS)
    aio_parser.add_argument('-n', '--unknown_symbol', required=False, type=str, default=DEFAULT_UNKNOWN_SYMBOL, help=HELP_TEXT_UNKNOWN_SYMBOL)
    aio_parser.add_argument('-e', '--include_no_primer', action='store_true', help=HELP_TEXT_TRIM_INCLUDE_READS_NO_PRIMER)

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
    ref_genome_ID = ref_lines[0][1:].split()[0].strip()
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
        header.add_sample('sample')
        error("%s: %s" % (ERROR_TEXT_EMPTY_BED, primer_fn))
    primers.sort() # shouldn't be necessary (BED should be sorted), but just in case
    return primers

# create 
def create_VariantFile_object(output_variants_fn=None, ref_genome_ID=None):
    '''Create ``pysam.VariantFile`` object for the output VCF file

    Args:
        ``output_variants_fn`` (``str``): Filename of output VCF

    Returns:
        ``pysam.VariantFile``: Stream for writing output VCF
    '''
    # build VCF header
    header = pysam.VariantHeader()
    header.add_sample('sample')
    header.add_meta(key='AmpliPyVersion', value=VERSION)
    header.add_meta(key='source', value=' '.join(argv))
    header.add_meta(key='contig', items=[('ID',ref_genome_ID)])
    header.add_meta(key='FORMAT', items=[('ID','GT'), ('Number','1'), ('Type','String'), ('Description',"Genotype")])
    header.add_meta(key='INFO', items=[('ID','DP'), ('Number',1), ('Type','Integer'), ('Description',"Total Depth")])
    header.add_meta(key='INFO', items=[('ID','REF_DP'), ('Number','1'), ('Type','Integer'), ('Description',"Depth of reference base")])
    header.add_meta(key='INFO', items=[('ID','ALT_DP'), ('Number','1'), ('Type','String'), ('Description',"Depth of alternate base")])
    header.add_meta(key='INFO', items=[('ID','REF_FREQ'), ('Number','1'), ('Type','Float'), ('Description',"Frequency of reference base")])
    header.add_meta(key='INFO', items=[('ID','ALT_FREQ'), ('Number','1'), ('Type','String'), ('Description',"Frequency of alternate base")])

    # open VCF file
    if output_variants_fn is None:
        return None
    elif output_variants_fn.lower() == 'stdout':
        return pysam.VariantFile('-', 'w', header=header)
    elif isfile(output_variants_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_EXISTS, output_variants_fn))
    elif output_variants_fn.lower().endswith('.vcf') or output_variants_fn.lower().endswith('.vcf.gz') or output_variants_fn.lower().endswith('.bcf'):
        return pysam.VariantFile(output_variants_fn, 'w', header=header)
    else:
        error("%s: %s" % (ERROR_TEXT_INVALID_VCF_EXTENSION, output_variants_fn))

# create AlignmentFile objects for SAM/BAM input and output files
def create_AlignmentFile_objects(input_reads_fn=None, output_reads_fn=None):
    '''Create ``pysam.AlignmentFile`` objects for the input and output SAM/BAM files

    Args:
        ``input_reads_fn`` (``str``): Filename of input reads SAM/BAM

        ``output_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM

    Returns:
        ``pysam.AlignmentFile``: Stream for reading input reads SAM/BAM

        ``pysam.AlignmentFile``: Stream for writing output trimmed reads SAM/BAM
    '''
    # disable htslib verbosity to avoid "no index file" warning
    tmp = pysam.set_verbosity(0)

    # open input SAM/BAM
    if input_reads_fn is None:
        in_aln = None
    elif input_reads_fn.lower() == 'stdin':
        in_aln = pysam.AlignmentFile('-', 'r') # standard input --> SAM
    elif not isfile(input_reads_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_NOT_FOUND, input_reads_fn))
    elif input_reads_fn.lower().endswith('.sam'):
        in_aln = pysam.AlignmentFile(input_reads_fn, 'r')
    elif input_reads_fn.lower().endswith('.bam'):
        in_aln = pysam.AlignmentFile(input_reads_fn, 'rb')
    else:
        error("%s: %s" % (ERROR_TEXT_INVALID_READ_EXTENSION, input_reads_fn))

    # prep output header (if applicable)
    if in_aln is None:
        error("Input alignment file is None")
    elif output_reads_fn is not None:
        header_dict = in_aln.header.to_dict()
        curr_dict = {
            'PN': 'AmpliPy',
            'PP': header_dict['PG'][-1]['ID'],
            'VN': VERSION,
            'CL': ' '.join(argv),
        }
        num_existing_amplipy = sum(item['PN'] == 'AmpliPy' for item in header_dict['PG'])
        if num_existing_amplipy == 0:
            curr_dict['ID'] = 'AmpliPy'
        else:
            curr_dict['ID'] = 'AmpliPy.%d' % num_existing_amplipy
        header_dict['PG'].append(curr_dict)

    # open output SAM/BAM (if applicable)
    if output_reads_fn is None:
        out_aln = None
    elif output_reads_fn.lower() == 'stdout':
        out_aln = pysam.AlignmentFile('-', 'w', header=header_dict) # standard output --> SAM
    elif isfile(output_reads_fn):
        error("%s: %s" % (ERROR_TEXT_FILE_EXISTS, output_reads_fn))
    elif output_reads_fn.lower().endswith('.sam'):
        out_aln = pysam.AlignmentFile(output_reads_fn, 'w', header=header_dict)
    elif output_reads_fn.lower().endswith('.bam'):
        out_aln = pysam.AlignmentFile(output_reads_fn, 'wb', header=header_dict)
    else:
        error("%s: %s" % (ERROR_TEXT_INVALID_READ_EXTENSION, output_reads_fn))

    # re-enable htslib verbosity and finish up
    pysam.set_verbosity(tmp)
    return in_aln, out_aln

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
    proper_cigar = list(); len_new_cigar = len(new_cigar)
    for i in range(len_new_cigar):
        if i < len_new_cigar-1 and new_cigar[i][0] == new_cigar[i+1][0]:
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

            # If we didn't consume the query but did consume the reference, still need to move start
            elif CONSUME_REF[cig]:
                ref_add += n

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
    window = min(sliding_window_width, true_end)

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
        start_pos = get_pos_on_ref(s.cigartuples, del_len + s.query_alignment_start - 1, s.reference_start)

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

# get the alleles at a given reference position from the symbol counts at that position
def alleles_from_counts(symbol_counts):
    '''Get the alleles at a given reference position from the symbol counts at that position

    Args:
        ``symbol_counts`` (``dict``): The symbol counts at a position of the reference (keys = symbols, values = counts)

    Returns:
        ``int``: The total coverage at the reference position

        ``list`` of ``(int,float,str)``: The alleles at the reference position as a ``list`` of ``(count, frequency, symbol)`` tuples, sorted in descending order of count
    '''
    total_coverage = sum(symbol_counts.values())
    if total_coverage == 0:
        return 0, list()
    else:
        return total_coverage, sorted(((symbol_counts[k], symbol_counts[k]/total_coverage, k) for k in symbol_counts if symbol_counts[k] != 0), reverse=True)

# run AmpliPy
def run_amplipy(
    untrimmed_reads_fn = None,
    primer_fn = None,
    reference_fn = None,
    trimmed_reads_fn = None,
    variants_fn = None,
    consensus_fn = None,
    primer_pos_offset = None,
    min_length = None,
    min_quality = None,
    sliding_window_width = None,
    min_freq_consensus = None,
    min_freq_variants = None,
    min_depth_consensus = None,
    min_depth_variants = None,
    unknown_symbol = None,
    include_no_primer = None,
    run_trim = False,
    run_variants = False,
    run_consensus = False,
):
    '''Run AmpliPy

    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of untrimmed reads SAM/BAM

        ``primer_fn`` (``str``): Filename of input primer BED

        ``reference_fn`` (``str``): Filename of input reference genome FASTA

        ``trimmed_reads_fn`` (``str``): Filename of trimmed reads SAM/BAM

        ``variants_fn`` (``str``): Filename of output variants VCF

        ``consensus_fn`` (``str``): Filename of output consensus sequence FASTA

        ``primer_pos_offset`` (``int``): Trimming: Primer position offset. Reads that occur at the specified offset positions relative to primer positions will also be trimmed

        ``min_length`` (``int``): Trimming: Minimum length of read to retain after trimming

        ``min_quality`` (``int``): Trimming: Minimum quality threshold for sliding window to pass

        ``sliding_window_width`` (``int``): Trimming: Width of sliding window

        ``min_freq_consensus`` (``float``): Consensus: Minimum frequency for consensus calling

        ``min_freq_variants`` (``float``): Variants: Minimum frequency for variant calling

        ``min_depth_consensus`` (``int``): Consensus: Minimum depth for consensus calling

        ``min_depth_variants`` (``int``): Variants: Minimum depth for consensus calling

        ``unknown_symbol`` (``str``): Consensus: Symbol for unknown nucleotides in consensus (e.g. N)

        ``include_no_primer`` (``bool``): Trimming: ``True`` to include reads with no primers trimmed, otherwise ``False``

        ``run_trim`` (``bool``): ``True`` to trim reads, otherwise ``False``

        ``run_variants`` (``bool``): ``True`` to call variants, otherwise ``False``

        ``run_consensus`` (``bool``): ``True`` to call consensus sequence, otherwise ``False``
    '''
    # validity check for non-file arguments
    if primer_pos_offset is not None and primer_pos_offset < 0:
        error("%s: %s" % (ERROR_TEXT_NEGATIVE_PRIMER_POS_OFFSET, primer_pos_offset))
    if min_length is not None and min_length < 1:
        error("%s: %s" % (ERROR_TEXT_INVALID_MIN_LENGTH, min_length))
    if min_quality is not None and min_quality < 0:
        error("%s: %s" % (ERROR_TEXT_NEGATIVE_MIN_QUALITY, min_quality))
    if sliding_window_width is not None and sliding_window_width < 1:
        error("%s: %s" % (ERROR_TEXT_INVALID_SLIDING_WINDOW_WIDTH, sliding_window_width))
    if min_freq_consensus is not None and (min_freq_consensus < 0 or min_freq_consensus > 1):
        error("%s: %s" % (ERROR_TEXT_INVALID_MIN_FREQ, min_freq_consensus))
    if min_freq_variants is not None and (min_freq_variants < 0 or min_freq_variants > 1):
        error("%s: %s" % (ERROR_TEXT_INVALID_MIN_FREQ, min_freq_variants))
    if min_depth_consensus is not None and min_depth_consensus < 0:
        error("%s: %s" % (ERROR_TEXT_INVALID_MIN_DEPTH, min_depth_consensus))
    if min_depth_variants is not None and min_depth_variants < 0:
        error("%s: %s" % (ERROR_TEXT_INVALID_MIN_DEPTH, min_depth_variants))
    if unknown_symbol is not None and len(unknown_symbol) != 1:
        error("%s: %s" % (ERROR_TEXT_INVALID_UNKNOWN_SYMBOL_LENGTH, unknown_symbol))

    # print AmpliPy mode details
    if not (run_trim or run_variants or run_consensus):
        error("Not running any of the AmpliPy operations")
    if run_trim and not (run_variants or run_consensus):
        print_log("Executing AmpliPy Trim (v%s)" % VERSION)
    elif run_variants and not (run_trim or run_consensus):
        print_log("Executing AmpliPy Variants (v%s)" % VERSION)
    elif run_consensus and not (run_trim or run_variants):
        print_log("Executing AmpliPy Consensus (v%s)" % VERSION)
    else:
        print_log("Executing AmpliPy All-In-One (v%s)" % VERSION)

    # load input files and preprocess
    if reference_fn is not None:
        print_log("Loading reference genome: %s" % reference_fn)
        ref_genome_ID, ref_genome_sequence = load_ref_genome(reference_fn)
        ref_genome_len = len(ref_genome_sequence)
    if primer_fn is not None:
        print_log("Loading primers: %s" % primer_fn)
        primers = load_primers(primer_fn)
        max_primer_len = max(end-start for start,end in primers)
        print_log("Precalculating overlapping primers...")
        min_primer_start, max_primer_end = find_overlapping_primers(ref_genome_len, primers, primer_pos_offset)
    if run_trim:
        print_log("Input untrimmed SAM/BAM: %s" % untrimmed_reads_fn)
        print_log("Output trimmed SAM/BAM: %s" % trimmed_reads_fn)
        in_aln, out_aln = create_AlignmentFile_objects(untrimmed_reads_fn, trimmed_reads_fn)
    else:
        print_log("Input trimmed SAM/BAM: %s" % trimmed_reads_fn)
        in_aln, DUMMY = create_AlignmentFile_objects(trimmed_reads_fn, None)
    if variants_fn is not None:
        print_log("Output variants VCF: %s" % variants_fn)
        out_vcf = create_VariantFile_object(variants_fn, ref_genome_ID)

    # if variant/consensus calling, initialize symbol counts
    if run_variants or run_consensus:
        symbol_counts_at_ref_pos = [{'A':0,'C':0,'G':0,'T':0,'N':0,'-':0} for _ in range(ref_genome_len)] # [i] = symbol counts at reference position i

    # process reads
    print_log("Processing reads...")
    s_i = 0
    for s in in_aln:
        # print progress update
        s_i += 1
        if s_i % PROGRESS_NUM_READS == 0:
            print_log("Processed %d reads..." % s_i)

        # skip unmapped reads
        if s.is_unmapped:
            continue

        # skip reads without CIGAR
        if s.cigartuples is None:
            continue

        # trim this read (if applicable)
        if run_trim:
            trimmed_primer_start, trimmed_primer_end, trimmed_quality = trim_read(s, min_primer_start, max_primer_end, max_primer_len, min_quality, sliding_window_width)

            # write this read (if applicable)
            if s.reference_length >= min_length and (trimmed_primer_start or trimmed_primer_end or include_no_primer):
                out_aln.write(s)

        # if variant/consensus calling, update base counts
        if run_variants or run_consensus:
            query_start = s.query_alignment_start # This the index of the first base in seq that is not soft-clipped
            query_end = s.query_alignment_end # This the index JUST PAST the last base in seq that is not soft-clipped
            query_seq = s.query_sequence.upper() # read sequence bases, including soft clipped bases (None if not present)
            query_qual = s.query_qualities
            ref_start = s.reference_start # 0-based leftmost coordinate
            ref_end = s.reference_end # reference_end points to ONE PAST the last aligned residue
            pos_pairs = list(s.get_aligned_pairs()) # list of (q_pos, r_pos) tuples
            len_pos_pairs = len(pos_pairs)
            i = 0
            while i < len_pos_pairs:
                # get this pair and move to next
                q_pos, r_pos = pos_pairs[i]; i += 1

                # deletion
                if q_pos is None:
                    symbol_counts_at_ref_pos[r_pos]['-'] += 1

                # too low of quality (so ignore)
                elif query_qual[q_pos] < min_quality:
                    continue

                # soft-clipped beginning (so ignore)
                elif q_pos < query_start:
                    continue

                # soft-clipped end (so early terminate)
                elif q_pos >= query_end:
                    break

                # insertion
                elif r_pos is None:
                    # search for next reference match
                    q_pos_insertion_start = q_pos
                    while q_pos < query_end and query_qual[q_pos] >= min_quality and r_pos is None:
                        q_pos, r_pos = pos_pairs[i]; i += 1
                    if r_pos == 0:
                        insertion_seq = query_seq[q_pos_insertion_start : q_pos + 1] # if this is an insertion before the reference genome, I need to add 1 letter to the end for the sake of variant calling
                    else:
                        insertion_seq = query_seq[q_pos_insertion_start - 1 : q_pos] # end needs to be exclusive because this is the NEXT match, and I need to include the letter JUST BEFORE the insertion for the sake of variant calling
                    if r_pos is None: # never found another reference match (so this is an insertion at the end of the alignment)
                        ref_insertion_pos = ref_end
                    else: # otherwise, this is a normal insertion before r_pos
                        ref_insertion_pos = r_pos
                        i -= 1 # I need to handle this reference position next loop
                    ref_insertion_pos = max(ref_insertion_pos-1, 0) # move back 1
                    if insertion_seq in symbol_counts_at_ref_pos[ref_insertion_pos]:
                        symbol_counts_at_ref_pos[ref_insertion_pos][insertion_seq] += 1
                    else:
                        symbol_counts_at_ref_pos[ref_insertion_pos][insertion_seq] = 1

                # match/mismatch
                else:
                    if query_qual[q_pos] >= min_quality:
                        symbol_counts_at_ref_pos[r_pos][query_seq[q_pos]] += 1

    # call variants and/or consensus (if applicable)
    if run_variants or run_consensus:
        if run_consensus:
            consensus_symbols = list()
        for ref_pos in range(ref_genome_len):
            # get symbol counts
            ref_symbol = ref_genome_sequence[ref_pos]
            symbol_counts_here = symbol_counts_at_ref_pos[ref_pos]
            total_depth_at_ref_pos, alleles_at_ref_pos = alleles_from_counts(symbol_counts_here)

            # if calling consensus, add current position and insertions after current position
            if run_consensus:
                if len(alleles_at_ref_pos) != 0 and alleles_at_ref_pos[0][0] >= min_depth_consensus and alleles_at_ref_pos[0][1] >= min_freq_consensus:
                    consensus_symbols.append(alleles_at_ref_pos[0][2])
                else:
                    consensus_symbols.append(unknown_symbol)

            # if calling variants, form and write VCF entry
            if run_variants:
                ref_symbol_count = 0; ref_symbol_freq = 0; alt_allele_symbols = list(); alt_allele_counts = list(); alt_allele_freqs = list()
                for count, freq, symbol in alleles_at_ref_pos:
                    if symbol == ref_symbol:
                        ref_symbol_count = count; ref_symbol_freq = freq
                    elif count >= min_depth_variants and freq >= min_freq_variants: # TODO should count be >= min_depth_variants? or total depth at this site?
                        alt_allele_symbols.append(symbol); alt_allele_counts.append(count); alt_allele_freqs.append(freq)
                if len(alt_allele_symbols) != 0:
                    vcf_info = dict()
                    vcf_info['DP'] = total_depth_at_ref_pos
                    vcf_info['REF_DP'] = ref_symbol_count
                    vcf_info['ALT_DP'] = ','.join(str(count) for count in alt_allele_counts)
                    vcf_info['REF_FREQ'] = ref_symbol_freq
                    vcf_info['ALT_FREQ'] = ','.join(str(freq) for freq in alt_allele_freqs)
                    vcf_record = out_vcf.new_record(contig=ref_genome_ID, start=ref_pos, stop=ref_pos+1, alleles=[ref_symbol]+alt_allele_symbols, info=vcf_info, filter='PASS')
                    if ref_symbol_count >= min_depth_variants and ref_symbol_freq >= min_freq_variants:
                        vcf_record.samples['sample']['GT'] = tuple(range(len(alt_allele_symbols)+1))
                    else:
                        vcf_record.samples['sample']['GT'] = tuple(range(1, len(alt_allele_symbols)+1))
                    out_vcf.write(vcf_record)

    # write consensus sequence to disk (if applicable)
    if run_consensus:
        if consensus_fn.lower().endswith('.gz'):
            f = gzip.open(consensus_fn, 'wb')
        else:
            f = open(consensus_fn, 'w')
        f.write('>sample\n%s\n' % ''.join(consensus_symbols)); f.close()

    # finish up
    print_log("Finished Processing %d reads" % s_i)

# main content
if __name__ == "__main__":
    if len(argv) == 1:
        pass # TODO: In the future, run GUI here to fill in argv accordingly (so argparse will run fine)
    args = parse_args()
    if args.command == 'trim':
        run_amplipy(
            untrimmed_reads_fn = args.input,
            primer_fn = args.primer,
            reference_fn = args.reference,
            trimmed_reads_fn = args.output,
            primer_pos_offset = args.primer_pos_offset,
            min_length = args.min_length,
            min_quality = args.min_quality,
            sliding_window_width = args.sliding_window_width,
            include_no_primer = args.include_no_primer,
            run_trim = True,
        )
    elif args.command == 'variants':
        run_amplipy(
            trimmed_reads_fn = args.input,
            reference_fn = args.reference,
            variants_fn = args.output,
            min_quality = args.min_quality,
            min_freq_variants = args.min_freq,
            min_depth_variants = args.min_depth,
            run_variants = True,
        )
    elif args.command == 'consensus':
        run_amplipy(
            trimmed_reads_fn = args.input,
            reference_fn = args.reference,
            consensus_fn = args.output,
            min_quality = args.min_quality,
            min_freq_consensus = args.min_freq,
            min_depth_consensus = args.min_depth,
            unknown_symbol = args.unknown_symbol,
            run_consensus = True,
        )
    elif args.command == 'aio':
        run_amplipy(
            untrimmed_reads_fn = args.input,
            primer_fn = args.primer,
            reference_fn = args.reference,
            trimmed_reads_fn = args.output_trimmed_reads,
            variants_fn = args.output_variants,
            consensus_fn = args.output_consensus,
            primer_pos_offset = args.primer_pos_offset,
            min_length = args.min_length,
            min_quality = args.min_quality,
            sliding_window_width = args.sliding_window_width,
            min_freq_consensus = args.min_freq_consensus,
            min_freq_variants = args.min_freq_variants,
            min_depth_consensus = args.min_depth_consensus,
            min_depth_variants = args.min_depth_variants,
            unknown_symbol = args.unknown_symbol,
            include_no_primer = args.include_no_primer,
            run_trim = True,
            run_variants = True,
            run_consensus = True,
        )
