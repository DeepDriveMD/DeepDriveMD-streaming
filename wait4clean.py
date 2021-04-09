import sys
import os
import time
import subprocess

fn = sys.argv[1]

sleeptime = 30

while(True):
    if(os.path.exists(fn)):
        print(f"{fn} found. Exiting...")
        break
    print(f"Sleeping for {sleeptime} seconds")
    time.sleep(sleeptime)




