from multiprocessing import Queue, Process
import time


class Trimmer:
    def __init__(self, bam, output):
        queue = Queue()
        idx = 0
        self.processes = []

        for read in bam:
            process = Process(target=self.process_read,
                              args=([read, idx, queue]))
            self.processes.append(process)
            process.start()
            idx += 1

        for process in self.processes:
            process.join()

        results = [queue.get() for p in self.processes]
        for p in self.processes:
            r = queue.get()
            output.write(r)

    def process_read(self, r, idx, queue):
        queue.put((idx, r))
