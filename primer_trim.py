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

# File opening and handling, currently hard coded for ease of development
bam = pysam.AlignmentFile("Mapped1.bam", "rb")
output = pysam.AlignmentFile("Trimmed1.sam", "w", header=bam.header)
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

# Iterating through the reads
reads = []
for r in bam:
	# If we don't have a CIGAR, we should just skip.
	if r.cigartuples == None:
		continue

	# If we have a reverse read, we need to go backwards
	if (r.is_reverse):
		# Check what is overlapping from the start of the reverse (the end of the read)
		overlapping = getOverlappingPrimers(r.reference_end - 1)

		# If nothing overlaps, we don't need to continue trimming
		if len(overlapping) > 0:
			# Find the minimum, or the largest we'd need to trim
			overlap_start = float("inf")
			for primer in overlapping:
				if primer[0] < overlap_start:
					overlap_start = primer[0]
			
			# Determine how much we need to delete
			maxDeleteLen = r.infer_query_length() - getPosOnQuery(r.cigartuples, overlap_start, r.reference_start)
			maxDeleteLen = maxDeleteLen if maxDeleteLen > 0 else 0

			# Initializations
			newcigar = []
			ref_add = 0
			delLen = maxDeleteLen
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

			for i in reversed(range(0, len(newcigar))):
				cigarstr = cigarstr + str(newcigar[i][1]) + cigarMap[newcigar[i][0]]
			
			r.cigarstring = cigarstr
	else: # We are a forward read, so we want to continue forwards
		# Check what is overlapping from the start of the reverse (the end of the read)
		overlapping = getOverlappingPrimers(r.reference_start)

		# If nothing overlaps, we don't need to continue trimming
		if len(overlapping) > 0:
			# Find the maximum, or the largest we'd need to trim
			overlap_end = 0
			for primer in overlapping:
				if primer[1] > overlap_end:
					overlap_end = primer[1]

			# Determine how much we need to delete
			maxDeleteLen = getPosOnQuery(r.cigartuples, overlap_end, r.reference_start)
			maxDeleteLen = maxDeleteLen if maxDeleteLen > 0 else 0

			# Initializations
			newcigar = []
			ref_add = 0
			delLen = maxDeleteLen
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

			for i in range(0, len(newcigar)):
				cigarstr = cigarstr + str(newcigar[i][1]) + cigarMap[newcigar[i][0]]

			r.cigarstring = cigarstr

			# Move our position on the reference forward, if needed
			r.reference_start += start_pos

	
	# Append the trimmed read
	reads.append(r)

# Write our reads
for r in reads:
	output.write(r)