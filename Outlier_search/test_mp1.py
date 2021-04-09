from multiprocessing import Pool
import math
import random
import os
import pickle
from pprint import pprint

#os.environ["MKL_NUM_THREADS"] = "1" 
#os.environ["NUMEXPR_NUM_THREADS"] = "1" 
#os.environ["OMP_NUM_THREADS"] = "1" 

pprint(dict(os.environ))


import numpy as np

def compute(x):
    y = x

    for i in range(1000):
        y = np.dot(y, y)
        m = np.max(y)
        if(m != 0.0):
            y /= m
    return y

def run():
    N = 2
    n = 1000
    a = []
    for k in range(100):
        a.append(np.random.random((n, n)))


    with open('env.pickle', 'wb') as f:
        pickle.dump(dict(os.environ), f)


    j = 0
    while(True):
        print(f"Iteration {j}")
        j += 1
        with Pool(N) as p:
            print(p.map(compute, a))


run()
