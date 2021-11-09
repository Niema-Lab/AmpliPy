from multiprocessing import Process, Queue
import random
import time

data = []


def rand_num(amt):
    print("process " + str(amt))
    num = random.random()
    time.sleep(amt)
    print(data)
    data.append({
        "id": amt,
        "num": num
    })
    print(data)
    print(num)


if __name__ == "__main__":
    queue = Queue()

    processes = [Process(target=rand_num, args=([x])) for x in range(4)]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    print(data)
