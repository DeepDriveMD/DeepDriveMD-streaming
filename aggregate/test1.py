from utils import *

import numpy as np

n=10

b = []
for i in range(3):
    a = np.random.random(n)
    a=a.reshape([1,n])
    at = a.T
    b.append(at)

print("="*5 + " b = " + "="*5)
print(b)

cm_all = np.hstack(b)
print("="*5 + " cm_all = " + "="*5)
print(cm_all)

c = np.array([triu_to_full(d) for d in cm_all.T])

print("="*5 + " c = " + "="*5)
print(c)
print(c.shape)

pad_f = lambda x: (0,0) if x%2 == 0 else (0,1)
padding_buffer=[(0,0)]
for x in c.shape[1:]:
    padding_buffer.append(pad_f(x))

#print(padding_buffer)

cm_data_full = np.pad(c, padding_buffer, mode='constant')
print("="*5 + " cm_data_full = " + "="*5)
print(cm_data_full)
print(cm_data_full.shape)



cc = cm_data_full.reshape(cm_data_full.shape + (1,))
print("="*5 + " cc = " + "="*5)
print(cc)
print(cc.shape)
