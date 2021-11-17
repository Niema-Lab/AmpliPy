import sys
import time

# Determine which primers overlap with a given start of a read
def get_overlapping_primers_original(start, primers):
	overlapping = []
	for primer in primers:
		# If the start is between the primer's start and end, it's an overlap
		# [-----primer[0]-----start-----primer[1]-----] as an example, is an overlap
		if start >= primer[0] and start <= primer[1]:
			overlapping.append(primer)
	
	return overlapping

def find_surrounding_primers(primers, pos, mid, backwards):
	if mid < 0 or mid >= len(primers): return []
	direct = -1 if backwards else 1
	if (pos >= primers[mid][0] and pos <= primers[mid][1]):
		return [primers[mid]] + find_surrounding_primers(primers, pos, mid+direct, backwards)
	return []


def find_primer_binary(primers, pos, low, high):
	low = max(0, low)
	high = min(high, len(primers)-1)
	# print(low, high)
	if low <= high:
		mid = (low + high) // 2
		# print(mid)
		# print(primers[mid])
		if (pos >= primers[mid][0] and pos <= primers[mid][1]):
			# print("FOUND",mid)
			return [primers[mid]] + find_surrounding_primers(primers, pos, mid-1, True) + find_surrounding_primers(primers, pos, mid+1, False)
		elif pos < primers[mid][0]:
			# print("less")
			return find_primer_binary(primers, pos, low, mid-1)
		else:
			# print("greater")
			return find_primer_binary(primers, pos, mid+1, high)
	# print("HERE")
	return []


def get_overlapping_primers_new(start, primers):
	overlapping = find_primer_binary(primers, start, 0, len(primers)-1)
	return overlapping





# Test the implementations
##### Build our primer list #####
primer_file = open("sarscov2_v2_primers.bed")
max_primer_len = 0
primers = []
for primer in primer_file:
	data = primer.split()
	start = int(data[1])
	# End isn't 0 based in bed
	end = int(data[2]) - 1

	primers.append((start, end))
	# Determine the longest primer
	if end - start + 1 > max_primer_len:
		max_primer_len = end - start + 1

sorted_primers = primers
sorted_primers.sort(key=lambda e : e[0])

print("CHECKING CORRECTNESS")
for i in range(30000):
	original = get_overlapping_primers_original(i, primers)
	new = get_overlapping_primers_new(i, sorted_primers)
	if original != new:
		print("FAILED!")
		print("i:", i)
		print("Original:", original)
		print("New:", new)
		sys.exit()

print("PASSED CORRECTNESS")
print("CHECKING TIME")

start_time = time.time()

for i in range(30000):
	original = get_overlapping_primers_original(i, primers)
	# new = get_overlapping_primers_new(i, primers)

finish_original = time.time()

for i in range(30000):
	# original = get_overlapping_primers_original(i, primers)
	new = get_overlapping_primers_new(i, primers)

finish_all = time.time()

print("FINISHED TIMING")
print("Original:", finish_original - start_time)
print("New:", finish_all - finish_original)