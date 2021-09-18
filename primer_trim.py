import pysam
import sys
import numpy as np

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
bam = pysam.AlignmentFile("Mapped.bam", "r")
output = pysam.AlignmentFile("Trimmed_Primer_Only.sam", "w", header=bam.header)
primerFile = open("/home/josh/ivar-examples/primers/swift/sarscov2_v2_primers.bed")

# Information about cigar operations
cigarMap = ["M", "I", "D", "N", "S", "H", "P", "=", "X"]
consumeQuery = [True, True, False, False, True, False, False, True, True]
consumeReference = [True, False, True, True, False, False, False, True, True]

max_primer_len = 0

# Build our primer list
primers = []
for primer in primerFile:
	data = primer.split()
	start = int(data[1])
	# End isn't 0 based in bed
	end = int(data[2]) - 1

	primers.append((start, end))
	# Determine the longest primer
	if end - start + 1 > max_primer_len:
		max_primer_len = end - start + 1

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

def cigarToQLen(cigar):
	len = 0
	for operation in cigar:
		if consumeQuery[operation[0]]:
			len += operation[1]
	
	return len

reads = []
for r in bam:
	overlapping_start = getOverlappingPrimers(r.reference_start)
	overlapping_end = getOverlappingPrimers(r.reference_end - 1)

	isize_flag = abs(r.template_length) - max_primer_len > abs(cigarToQLen(r.cigartuples))

	ref_start_offset = 0

	if not (r.is_paired and isize_flag and r.is_reverse) and len(overlapping_start) > 0:
		overlap_end = 0
		for primer in overlapping_start:
			if primer[1] > overlap_end:
				overlap_end = primer[1]
		max_delete_start = getPosOnQuery(r.cigartuples, overlap_end + 1, r.reference_start)
		max_delete_start = max(max_delete_start, 0)

		# Initializations
		newcigar = []
		ref_add = 0
		delLen = max_delete_start
		pos_start = False
		start_pos = 0

		# For each operation in our Cigar
		for operation in r.cigartuples:
			# If we have nothing left to delete, just append everything
			if (delLen == 0 and pos_start):
				newcigar.append(operation)
				continue
			
			# For convenience, this is the data about our current cigar operation
			cig = operation[0]
			n = operation[1]

			# If we have nothing left to delete and are consuming both, we want to just append everything
			if (delLen == 0 and consumeQuery[cig] and consumeReference[cig]):
				pos_start = True
				newcigar.append(operation) # Remember we need to include this one!
				continue
			
			# How much our current trim affects our read's start position
			ref_add = 0

			# If our operation consumes the query
			if consumeQuery[cig]:
				# How much do we have to delete?
				if delLen >= n:
					# Our entire operation needs to be deleted
					newcigar.append((BAM_CSOFT_CLIP, n))
				elif delLen < n and delLen > 0:
					# We need to delete some of our segment, but we will still have more later
					newcigar.append((BAM_CSOFT_CLIP, delLen))
				elif delLen == 0:
					# Since we consume the query, we just need to keep clipping
					newcigar.append((BAM_CSOFT_CLIP, n))
					continue
				
				# Update based on how much we just deleted
				ref_add = min(delLen, n)
				temp = n
				n = max(n - delLen, 0)
				delLen = max(delLen - temp, 0)

				# If there is still more left to do, append it
				if n > 0:
					newcigar.append((cig, n))

				# If we are done and just consumed, we want to just start appending everything.
				if delLen == 0 and consumeQuery[newcigar[-1][0]] and consumeReference[newcigar[-1][0]]:
					pos_start = True
			
			# If our trim consumed the reference, we need to move our read's start position forwards
			if consumeReference[cig]:
				start_pos += ref_add

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
		r.reference_start += start_pos
		ref_start_offset += start_pos

	if not (r.is_paired and isize_flag and not r.is_reverse) and len(overlapping_end) > 0:
		overlap_start = float("inf")
		for primer in overlapping_end:
			if primer[0] < overlap_start:
				overlap_start = primer[0]
		max_delete_end = cigarToQLen(r.cigartuples) - getPosOnQuery(r.cigartuples, overlap_start, r.reference_start)

		# Initializations
		newcigar = []
		ref_add = 0
		delLen = max_delete_end
		pos_start = False

		# For each operation in our Cigar
		for operation in reversed(r.cigartuples):
			# If we have nothing left to delete, just append everything
			if (delLen == 0 and pos_start):
				newcigar.append(operation)
				continue
			
			# For convenience, this is the data about our current cigar operation
			cig = operation[0]
			n = operation[1]

			# If we have nothing left to delete and are consuming both, we want to just append everything
			if (delLen == 0 and consumeQuery[cig] and consumeReference[cig]):
				pos_start = True
				newcigar.append(operation) # Remember we need to include this one!
				continue
			
			# If our operation consumes the query
			if consumeQuery[cig]:
				# How much do we have to delete?
				if delLen >= n:
					# Our entire operation needs to be deleted
					newcigar.append((BAM_CSOFT_CLIP, n))
				elif delLen < n and delLen > 0:
					# We need to delete some of our segment, but we will still have more later
					newcigar.append((BAM_CSOFT_CLIP, delLen))
				elif delLen == 0:
					# Since we consume the query, we just need to keep clipping
					newcigar.append((BAM_CSOFT_CLIP, n))
					continue
				
				# Update based on how much we just deleted
				temp = n
				n = max(n - delLen, 0)
				delLen = max(delLen - temp, 0)

				# If there is still more left to do, append it
				if n > 0:
					newcigar.append((cig, n))

				# If we are done and just consumed, we want to just start appending everything.
				if delLen == 0 and consumeQuery[newcigar[-1][0]] and consumeReference[newcigar[-1][0]]:
					pos_start = True

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
	
	reads.append(r)

# Write our reads
for r in reads:
	if cigarToRefLen(r.cigartuples) >= MIN_READ_LENGTH:
		# print(r)
		output.write(r)
	# else:
		# print("removed", r)