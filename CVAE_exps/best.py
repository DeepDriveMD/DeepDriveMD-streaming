import glob
import numpy as np
import subprocess

losses = glob.glob("cvae_runs_*/loss.npy")

mylosses = {}

for loss in losses:
    a = np.load(loss)
    mylosses[loss]=a[-1]


best_loss = 10000
best_dir = ''

for k in mylosses.keys():
    if(mylosses[k] < best_loss):
        best_loss = mylosses[k]
        best_dir = k.split("/")[0]


print(f"best_loss = {best_loss}, best_dir = {best_dir}")

print(mylosses)

print(subprocess.getstatusoutput(f"ln -sf {best_dir} best_model"))




