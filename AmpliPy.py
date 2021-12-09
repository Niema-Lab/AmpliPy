#! /usr/bin/env python3
'''
AmpliPy: Python toolkit for viral amplicon sequencing
'''

# imports
import argparse
import pysam
#from trim import trim
from datetime import datetime
from os.path import isfile
from sys import argv, stderr

# constants
VERSION = '0.0.1'
BUFSIZE = 1048576 # 1 MB

# messages
HELP_TEXT_AMPLIPY_INDEX = "AmpliPy Index (PKL)"
HELP_TEXT_CONSENSUS = "Consensus Sequence (FASTA)"
HELP_TEXT_PRIMER = "Primer File (BED)"
HELP_TEXT_READS_UNTRIMMED = "Untrimmed Reads (SAM/BAM)"
HELP_TEXT_READS_TRIMMED = "Trimmed Reads (SAM/BAM)"
HELP_TEXT_REFERENCE = "Reference Genome (FASTA)"
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

    # AmpliPy Index args
    index_parser = subparsers.add_parser("index", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    index_parser.add_argument('-p', '--primer', required=True, type=str, help=HELP_TEXT_PRIMER)
    index_parser.add_argument('-r', '--reference', required=True, type=str, help=HELP_TEXT_REFERENCE)
    index_parser.add_argument('-o', '--output', required=True, type=str, help=HELP_TEXT_AMPLIPY_INDEX)

    # AmpliPy Trim args
    trim_parser = subparsers.add_parser("trim", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    trim_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_UNTRIMMED)
    trim_parser.add_argument('-a', '--amplipy_index', required=True, type=str, help=HELP_TEXT_AMPLIPY_INDEX)
    trim_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_READS_TRIMMED)

    # AmpliPy Variants args
    variants_parser = subparsers.add_parser("variants", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    variants_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    variants_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_VARIANTS)

    # AmpliPy Consensus args
    consensus_parser = subparsers.add_parser("consensus", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    consensus_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_TRIMMED)
    consensus_parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help=HELP_TEXT_CONSENSUS)

    # AmpliPy AIO (All-In-One) args
    aio_parser = subparsers.add_parser("aio", description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aio_parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help=HELP_TEXT_READS_UNTRIMMED)
    aio_parser.add_argument('-a', '--amplipy_index', required=True, type=str, help=HELP_TEXT_AMPLIPY_INDEX)
    aio_parser.add_argument('-ot', '--output_trimmed_reads', required=True, type=str, help=HELP_TEXT_READS_TRIMMED)
    aio_parser.add_argument('-ov', '--output_variants', required=True, type=str, help=HELP_TEXT_VARIANTS)
    aio_parser.add_argument('-oc', '--output_consensus', required=True, type=str, help=HELP_TEXT_CONSENSUS)

    # parse args and return
    return parser.parse_args()

# run AmpliPy Index
def run_index(primer_fn, reference_fn, amplipy_index_fn):
    '''Run AmpliPy Index

    Args:
        ``primer_fn`` (``str``): Filename of input primer BED

        ``reference_fn`` (``str``): Filename of input reference genome FASTA

        ``amplipy_index_fn`` (``str``): Filename of output AmpliPy index PKL
    '''
    print_log("Executing AmpliPy Index (v%s)" % VERSION)
    f = open(reference_fn, 'r', buffering=BUFSIZE); ref_lines = f.read().splitlines(); f.close()
    error("INDEX NOT IMPLEMENTED\n- primer_fn: %s\n- reference_fn: %s\n- amplipy_index_fn: %s" % (primer_fn, reference_fn, amplipy_index_fn)) # TODO

# run AmpliPy Trim
def run_trim(untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn):
    '''Run AmpliPy Trim

    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of input untrimmed reads SAM/BAM

        ``amplipy_index_fn`` (``str``): Filename of input AmpliPy index PKL

        ``trimmed_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM
    '''
    print_log("Executing AmpliPy Trim (v%s)" % VERSION)
    error("TRIM NOT IMPLEMENTED\n- untrimmed_reads_fn: %s\n- amplipy_index_fn: %s\n- trimmed_reads_fn: %s" % (untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn)) # TODO

# run AmpliPy Variants
def run_variants(trimmed_reads_fn, variants_fn):
    '''Run AmpliPy Variants

    Args:
        ``trimmed_reads_fn`` (``str``): Filename of input trimmed reads SAM/BAM

        ``variants_fn`` (``str``): Filename of output variants VCF
    '''
    print_log("Executing AmpliPy Variants (v%s)" % VERSION)
    error("VARIANTS NOT IMPLEMENTED\n- trimmed_reads_fn: %s\n- variants_fn: %s" % (trimmed_reads_fn, variants_fn)) # TODO

# run AmpliPy Consensus
def run_consensus(trimmed_reads_fn, consensus_fn):
    '''Run AmpliPy Consensus

    Args:
        ``trimmed_reads_fn`` (``str``): Filename of input trimmed reads SAM/BAM

        ``consensus_fn`` (``str``): Filename of output consensus sequence FASTA
    '''
    print_log("Executing AmpliPy Consensus (v%s)" % VERSION)
    error("CONSENSUS NOT IMPLEMENTED\n- trimmed_reads_fn: %s\n- consensus_fn: %s" % (trimmed_reads_fn, consensus_fn)) # TODO

# run AmpliPy AIO (All-In-One)
def run_aio(untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn, variants_fn, consensus_fn):
    '''Run AmpliPy AIO (All-In-One)

    Args:
        ``untrimmed_reads_fn`` (``str``): Filename of input untrimmed reads SAM/BAM

        ``amplipy_index_fn`` (``str``): Filename of input AmpliPy index PKL

        ``trimmed_reads_fn`` (``str``): Filename of output trimmed reads SAM/BAM

        ``variants_fn`` (``str``): Filename of output variants VCF

        ``consensus_fn`` (``str``): Filename of output consensus sequence FASTA
    '''
    print_log("Executing AmpliPy AIO (v%s)" % VERSION)
    error("AIO NOT IMPLEMENTED\n- untrimmed_reads_fn: %s\n- amplipy_index_fn: %s\n- trimmed_reads_fn: %s\n- variants_fn: %s\n- consensus_fn: %s" % (untrimmed_reads_fn, amplipy_index_fn, trimmed_reads_fn, variants_fn, consensus_fn)) # TODO

# main content
if __name__ == "__main__":
    args = parse_args()
    if args.command == 'index':
        run_index(args.primer, args.reference, args.output)
    elif args.command == 'trim':
        run_trim(args.input, args.amplipy_index, args.output)
    elif args.command == 'variants':
        run_variants(args.input, args.output)
    elif args.command == 'consensus':
        run_consensus(args.input, args.output)
    elif args.command == 'aio':
        run_aio(args.input, args.amplipy_index, args.output_trimmed_reads, args.output_variants, args.output_consensus)
    exit(1)

# OLD STUFF BELOW #

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
