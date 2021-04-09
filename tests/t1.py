from dask.distributed import Client, wait
import time

def inc(x):
    return x + 1

def add(x, y):
    return x + y


if __name__ == '__main__':
    client = Client(processes=True)
    print(client)
    a = client.submit(inc, 10)
    b = client.submit(inc, 20)
    c = client.submit(add, a, b)
    wait([c])
    print(c.result())
    futures = client.map(inc, range(1000))
    wait(futures)
    # print(list(map(lambda x: x.result(), futures)))
