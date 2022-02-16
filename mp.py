import threading
import multiprocessing as mp
from time import time
import numpy as np

sem = threading.Semaphore

def process_read(read):
	print("completed!")
	print(str(read))

def process_distributor(queue, signal_queue):
	# if we have added all jobs and we have no jobs left to do, break
	while not (signal_queue.empty() and not queue.empty()):
		sem.acquire()
		read = queue.get()
		process_read(read)

def main():
	bam = np.random.rand(1000)
	time_start = time()
	m = mp.Manager()
	queued_jobs = m.Queue()
	signal_queue = m.Queue()

	distributor1 = mp.Process(target=process_distributor,
							  args=(queued_jobs, signal_queue))
	distributor2 = mp.Process(target=process_distributor,
							  args=(queued_jobs, signal_queue))
	distributor3 = mp.Process(target=process_distributor,
							  args=(queued_jobs, signal_queue))
	
	distributor1.start()
	distributor2.start()
	distributor3.start()

	for read in bam:
		queued_jobs.put(read)
		sem.release()
	distributor1.join()
	distributor2.join()
	distributor3.join()

	time_end = time.time()
	print("mp took " + str(time_end - time_start))

	# linear approach
	time_start = time.time()
	for read in bam:
		process_read(read)
	
	time_end = time.time()
	print("linear took " + str(time_end - time_start))

main()