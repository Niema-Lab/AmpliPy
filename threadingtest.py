import asyncio
import numpy as np
import time
import multiprocessing as mp


data = []


def output(result):
    global data
    data.append(result)


async def process_read_async(read):
    sum = 0
    for data in read:
        sum += data
    output(sum)


def process_read_linear(read):
    sum = 0
    for data in read:
        sum += data
    output(sum)


async def main_threading_approach(reads):
    start_time = time.time()
    tasks = []

    for read in reads:
        tasks.append(process_read_async(read))

    for task in tasks:
        await task

    end_time = time.time()
    return end_time - start_time


def process_read_mp(read):
    sum = 0
    for data in read:
        sum += data
    output(sum)


def process_reads(queue, signal_queue):
    # if we have added all jobs and we have no jobs left to do, break
    while not (signal_queue.empty() and not queue.empty()):
        # check da queue
        # if something on da queue
        if not queue.empty():
            read = queue.get()
            process_read_mp(read)


def main_mp_approach(reads):
    start_time = time.time()

    m = mp.Manager()
    queued_jobs = m.Queue()
    signal_queue = m.Queue()

    distributor1 = mp.Process(target=process_reads,
                              args=(queued_jobs, signal_queue))
    distributor2 = mp.Process(target=process_reads,
                              args=(queued_jobs, signal_queue))
    distributor3 = mp.Process(target=process_reads,
                              args=(queued_jobs, signal_queue))

    distributor1.start()
    distributor2.start()
    distributor3.start()
    for read in reads:
        queued_jobs.put(read)

    signal_queue.put(True)  # signals to all that we're done
    distributor1.join()
    distributor2.join()
    distributor3.join()

    end_time = time.time()
    return end_time - start_time


def main_linear_approach(reads):
    start_time = time.time()
    for read in reads:
        process_read_linear(read)

    end_time = time.time()
    return end_time - start_time


reads = np.random.randint(10, 90, (1000, 10000))

loop = asyncio.get_event_loop()

sample_size = 10
avg_time_threading = 0
avg_time_linear = 0
avg_time_mp = 0
for i in range(sample_size):
    print("sample run", str(i))
    avg_time_threading += loop.run_until_complete(
        main_threading_approach(reads))
    avg_time_linear += main_linear_approach(reads)
    avg_time_mp += main_mp_approach(reads)


avg_time_threading /= sample_size
avg_time_linear /= sample_size
avg_time_mp /= sample_size

print("avg_time_threading", str(avg_time_threading))
print("avg_time_linear", str(avg_time_linear))
print("avg_time_mp", str(avg_time_mp))
