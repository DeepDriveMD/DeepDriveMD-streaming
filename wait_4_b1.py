import os
import time
import subprocess
import datetime
import sys

def get_now():
    return datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")

barrier_fn = sys.argv[1]

time.sleep(30)

while(not os.path.exists(barrier_fn)):
    t = get_now()
    print(f'wait {t}')
    time.sleep(30)

t = get_now()
subprocess.getstatusoutput(f'rm {barrier_fn}')
print(f"go {t}")

