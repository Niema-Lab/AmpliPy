import pysam
import sys
import numpy as np
import os, psutil

# Constants, so we can easily see what type of operation we want
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9

MIN_READ_LENGTH = 30

# File opening and handling, currently hard coded for ease of development
bam = pysam.AlignmentFile("Trimmed_Primer_Only.sam", "r")
output = pysam.AlignmentFile("Trimmed_Both.sam", "w", header=bam.header)
primerFile = open("/home/josh/ivar-examples/primers/swift/sarscov2_v2_primers.bed")

# Information about cigar operations
cigarMap = ["M", "I", "D", "N", "S", "H", "P", "=", "X"]
consumeQuery = [True, True, False, False, True, False, False, True, True]
consumeReference = [True, False, True, True, False, False, False, True, True]

# Build our primer list
primers = []
for primer in primerFile:
	data = primer.split()

	primers.append((int(data[1]), int(data[2])))

# Sort primers, this may help with optimizations later?
# primers.sort(key=lambda e : int(e[1]))

# Determine which primers overlap with a given start of a read
def getOverlappingPrimers(start):
	overlapping = []
	for primer in primers:
		# If the start is between the primer's start and end, it's an overlap
		# [-----primer[0]-----start-----primer[1]-----] as an example, is an overlap
		if start >= primer[0] and start <= primer[1]:
			overlapping.append(primer)
	
	return overlapping

# Convert a position on a reference to the position on the current read
def getPosOnQuery(cigar, pos, seg_start):
	# Initializations
	queryPos = 0
	curPos = seg_start

	# For each operation in our cigar
	for operation in cigar:
		# If we consume the reference
		if consumeReference[operation[0]]:
			# If our desired position is found somewhere inside the current operation
			if pos <= curPos + operation[1]:
				# If we consume the query, we want to advance the position on the query by as much as we need to go
				if consumeQuery[operation[0]]:
					queryPos += (pos - curPos)

				# Return, as our desired position is found in this operation
				return queryPos
			
			# Advance past the current operation, as it exists outside of the current operation
			curPos += operation[1]
		
		# If we consume the query, we want to advance our query's pointer
		if consumeQuery[operation[0]]:
			queryPos += operation[1]
		
	# We ran out of operations, so it must be here (or outside our query)
	return queryPos

def cigarToRefLen(cigar):
	len = 0
	for operation in cigar:
		if consumeReference[operation[0]]:
			len += operation[1]
	
	return len

# Convert a position on the read to a position on the reference
def getPosOnReference(cigar, pos, ref_start):
	curPos = 0
	referencePos = ref_start

	for operation in cigar:
		if consumeQuery[operation[0]]:
			if pos <= curPos + operation[1]:
				if consumeReference[operation[0]]:
					referencePos += (pos - curPos)
			
				return referencePos
		
			curPos += operation[1]
	
		if consumeReference[operation[0]]:
			referencePos += operation[1]
	
	return referencePos

QUALITY_THRESHOLD = 20
SLIDING_WINDOW = 4
itr = 0
# Iterating through the reads
reads = []
for r in bam:
	# itr = itr + 1
	# print(itr)
	# print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
	# If we don't have a CIGAR, we should just skip.
	if r.cigartuples == None:
		continue

	if r.is_reverse:
		sum = 0
		qual = r.query_alignment_qualities
		window = SLIDING_WINDOW if SLIDING_WINDOW <= len(qual) else len(qual)
		truestart = 0
		trueend = len(qual)

		i = trueend
		for offset in range(1, window):
			sum += qual[i - offset]

		while i > truestart:
			if truestart + window > i:
				window -= 1
			else:
				sum += qual[i - window]

			if sum / window < QUALITY_THRESHOLD:
				break
						
			sum -= qual[i - 1]
			i -= 1

		newcigar = []
		del_len = i
		start_pos = getPosOnReference(r.cigartuples, del_len + r.qstart, r.reference_start)

		if start_pos > r.reference_start:
			# # print("Need to delete", del_len)

			for operation in r.cigartuples:
				if (del_len == 0):
					newcigar.append(operation)
					continue
				
				cig = operation[0]
				n = operation[1]

				if cig == BAM_CSOFT_CLIP or cig == BAM_CHARD_CLIP:
					newcigar.append(operation)
					continue

				if consumeQuery[cig]:
					if del_len >= n:
						newcigar.append((BAM_CSOFT_CLIP, n))
					elif del_len < n:
						newcigar.append((BAM_CSOFT_CLIP, del_len))
					
					temp = n
					n = max(n - del_len, 0)
					del_len = max(del_len - temp, 0)

					if n > 0:
						newcigar.append((cig, n))
			
			# Update our cigar string, since that's what will be written
			cigarstr = ""
			propercigar = []
			for i in range(0, len(newcigar)):
				if i < len(newcigar)-1 and newcigar[i][0] == newcigar[i+1][0]:
					newcigar[i+1] = (newcigar[i+1][0], newcigar[i][1] + newcigar[i+1][1])
					continue
				
				cigarstr = cigarstr + str(newcigar[i][1]) + cigarMap[newcigar[i][0]]
				propercigar.append(newcigar[i])

			r.cigarstring = cigarstr
			r.cigartuples = propercigar

			# Move our position on the reference forward, if needed
			r.reference_start = start_pos
	else:
		sum = 0
		qual = r.query_alignment_qualities
		window = SLIDING_WINDOW if SLIDING_WINDOW <= len(qual) else len(qual)
		# # print(window)
		# # print(len(qual), qual)

		# Figure out where we should start and end in our quality array, since it includes clips
		truestart = 0
		trueend = len(qual)	

		i = truestart
		for offset in range(0, window - 1):
			sum += qual[i + offset]
		
		while i < trueend:
			if trueend - window < i:
				window -= 1
			else:
				sum += qual[i + window - 1]
			
			if sum / window < QUALITY_THRESHOLD:
				break
			
			sum -= qual[i]
			i += 1
		
		newcigar = []
		del_len = trueend - i
		start_pos = getPosOnReference(r.cigartuples, del_len, r.reference_start)

		for operation in reversed(r.cigartuples):
			if (del_len == 0):
				newcigar.append(operation)
				continue
			
			cig = operation[0]
			n = operation[1]

			if cig == BAM_CSOFT_CLIP or cig == BAM_CHARD_CLIP:
				newcigar.append(operation)
				continue

			if consumeQuery[cig]:
				if del_len >= n:
					newcigar.append((BAM_CSOFT_CLIP, n))
				elif del_len < n:
					newcigar.append((BAM_CSOFT_CLIP, del_len))
				
				temp = n
				n = max(n - del_len, 0)
				del_len = max(del_len - temp, 0)

				if n > 0:
					newcigar.append((cig, n))
		
		# Update our cigar string, since that's what will be written
		cigarstr = ""
		propercigar = []
		for i in reversed(range(0, len(newcigar))):
			if i > 0 and newcigar[i][0] == newcigar[i-1][0]:
				newcigar[i-1] = (newcigar[i-1][0], newcigar[i][1] + newcigar[i-1][1])
				continue

			cigarstr = cigarstr + str(newcigar[i][1]) + cigarMap[newcigar[i][0]]
			propercigar.append(newcigar[i])

		r.cigarstring = cigarstr
		r.cigartuples = propercigar

		# Move our position on the reference forward, if needed
		r.reference_start = start_pos

	# Append the trimmed read
	reads.append(r)

# Write our reads
for r in reads:
	if cigarToRefLen(r.cigartuples) >= MIN_READ_LENGTH:
		# print("not", r)
		output.write(r)
	# else:
		# print("removed", r)