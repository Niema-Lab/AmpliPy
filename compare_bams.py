import pysam
import sys
import os.path

file1 = sys.argv[1]
file2 = sys.argv[2]

if not os.path.isfile(file1):
	print("'%s' does not exist"%file1)
	sys.exit()

if not os.path.isfile(file2):
	print("'%s' does not exist"%file2)
	sys.exit()

bam1 = pysam.AlignmentFile(file1, "r")
bam2 = pysam.AlignmentFile(file2, "r")

reads = {}
reverses = {}
incorrects = []

for r in bam1:
	if r.is_reverse:
		reverses[r.query_name] = r
	else:
		reads[r.query_name] = r

for r in bam2:
	if r.is_reverse:
		if r.query_name in reverses:
			# compare them
			if r.cigarstring != reverses[r.query_name].cigarstring:
				incorrects.append([r, reverses[r.query_name]])
			reverses.pop(r.query_name)
		#else:
			# in right but not left
	else:
		if r.query_name in reads:
			# compare them
			if r.cigarstring != reads[r.query_name].cigarstring:
				incorrects.append([r, reads[r.query_name]])
			reads.pop(r.query_name)
		#else:
			# in right but not left
			

for [read, other] in incorrects:
	print("===================")
	print(read.query_name)
	print(read.cigarstring)
	print(other.cigarstring)

print(len(incorrects), len(reads), len(reverses))