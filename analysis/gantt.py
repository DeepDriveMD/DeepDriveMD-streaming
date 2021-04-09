#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

fn = sys.argv[1]
n_md = int(sys.argv[2])

pf = pd.read_csv(fn)

print(pf.date.min())
print(pf.date.max())


pfu = {}
pfd = {}

for u in pf['unit'].unique():
    tmp = pf[pf.unit == u].reset_index()
    pfu[u] = {}
    pfd[u] = {}
    for l in tmp.label.unique():
        pfu[u][l] = tmp[tmp.label == l].reset_index()
        starts = []
        ends = []
        durations = []
        units = []
        labels = []
        for i in range(0,len(pfu[u][l])//2*2,2):
            start = pfu[u][l].gps[i]
            end = pfu[u][l].gps[i+1]
            duration = end - start
            unit = pfu[u][l].unit[i]
            label = pfu[u][l].label[i]
            starts.append(start)
            ends.append(end)
            durations.append(duration)
            units.append(unit)
            labels.append(label)
        pfd[u][l] = pd.DataFrame(columns=["start","end","duration","unit","label"])
        pfd[u][l].start = starts
        pfd[u][l].end = ends
        pfd[u][l].duration = durations
        pfd[u][l].unit = units
        pfd[u][l].label = labels


merged_simulation_step = pfd[0]['simulation.step'].copy()
for u in range(1,n_md):
    merged_simulation_step = merged_simulation_step.append(pfd[u]['simulation.step'])

print(merged_simulation_step.duration.mean())

print(merged_simulation_step.describe())

