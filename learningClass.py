import random
import time
from multiprocessing import Process, Queue
from threading import Thread


class Runner:
    def __init__(self, payload=10):
        self.data = []
        self.queue = Queue()

        start_time = time.time()

        for x in range(payload):
            self.rand_num(x, "linear")
        print("[linear] elapsed time: " + str(time.time() - start_time))

        self.processes = [
            Process(target=self.rand_num, args=([x, "process", self.queue])) for x in range(payload)]

        self.threads = [
            Thread(target=self.rand_num, args=([x, "thread"])) for x in range(payload)
        ]

    def rand_num(self, amt, type, queue=3):
        # print("started " + str(type) + "with amt " + str(amt))
        num = random.random()
        entry = {
            "id": amt,
            "num": num
        }
        # print("finished " + str(type) + "with amt " + str(amt))
        if (type == "process"):
            queue.put(entry)

    def run_process(self):
        start_time = time.time()

        for p in self.processes:
            p.start()
        for p in self.processes:
            p.join()

        results = [self.queue.get() for p in self.processes]
        print(len(results))
        print("[process] elapsed time: " + str(time.time() - start_time))

    def run_threads(self):
        start_time = time.time()
        for t in self.threads:
            t.start()
        for t in self.threads:
            t.join()

        print("[threads] elapsed time: " + str(time.time() - start_time))


usain = Runner(3000)
# usain.run_threads()
usain.run_process()
