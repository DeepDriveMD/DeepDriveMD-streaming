import time
from inspect import currentframe, getframeinfo
import sys
import os
import subprocess
import glob

def mytime_label(label, start=1): #start = 1 - start, start = -1 - stop, start = 0 - neither
    t = time.localtime()
    gps = time.mktime(t)
    readable = time.asctime(t)
    frameinfo = getframeinfo(currentframe().f_back)
    fractions = time.perf_counter()
    print(f'TLaBeL|{label}|{start}|{gps}|{readable}|{frameinfo.filename}|{frameinfo.lineno}|{fractions}')
    sys.stdout.flush()

def mytop():
    print("="*10 + " top " + "="*10)
    user = os.getenv("USER")
    print(subprocess.getstatusoutput(f"top -U {user} -b -n 1")[1])
    print("="*5)
    print(subprocess.getstatusoutput("nvidia-smi")[1])
    print("="*25)

def get_bp_files(dir, n):
    flist = glob.glob(dir + "/*.bp")
    if(len(flist) != n):
        return None
    else:
        return flist

def bp_steps(flist):
    steps = 0
    for fn in flist:
        steps += int(subprocess.getstatusoutput(f"bpls {fn}")[1].split('\n')[0].split('*')[0].split(" ")[-1])
    return steps


