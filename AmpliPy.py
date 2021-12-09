#! /usr/bin/env python3
'''
AmpliPy: Python toolkit for viral amplicon sequencing
'''

# imports
import argparse
import pysam
from trim import trim
from sys import argv, stderr
from os.path import isfile

# messages
COMMANDS = '''COMMANDS:
- index     = Index primer BED + reference FASTA
- trim      = Trim SAM/BAM
- variants  = Call variants
- consensus = Call consensus sequence'''

# index primer data
def run_index_primers():
    print("INDEX NOT IMPLEMENTED", file=stderr); exit(1) # TODO

# main content
if __name__ == "__main__":
    if len(argv) != 1:
        command = argv[1].strip().lower()
        if command == 'index':
            run_index_primers()
        elif command == 'trim':
            print("NOT IMPLEMENTED", file=stderr); exit(1) # TODO
        elif command == 'variants':
            print("NOT IMPLEMENTED", file=stderr); exit(1) # TODO
        elif command == 'consensus':
            print("NOT IMPLEMENTED", file=stderr); exit(1) # TODO
    print("USAGE: %s <command>\n\n%s" % (argv[0], COMMANDS)); exit(1)

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
