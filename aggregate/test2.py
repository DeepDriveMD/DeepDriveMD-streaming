import numpy as np
import adios2

with adios2.open("aggregator.bp", "r") as fr:
    n = fr.steps()
    vars = fr.available_variables()
    print("vars = ", vars)
    print("n = ", n)
    z = fr.read('md5')
    print(z)

'''

    for v in vars:
        print(v)
        shape = list(map(int, vars[v]['Shape'].split(",")))
        zs = list(np.zeros(len(shape), dtype='int'))
        z = fr.read(v, zs, shape, 0, n)
        print(z.shape)
'''

'''
(OpenMM5) [iyakushin@login3.summit aggregate]$ python test.py
{'contact_map': {'AvailableStepsCount': '2538', 'Max': '1', 'Min': '0', 'Shape': '210', 'SingleValue': 'false', 'Type': 'double'}, 'positions': {'AvailableStepsCount': '2538', 'Max': '21.9658', 'Min': '-20.5932', 'Shape': '264, 3', 'SingleValue': 'false', 'Type': 'double'}}
(2538, 210)
'''
