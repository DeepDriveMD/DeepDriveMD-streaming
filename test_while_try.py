import sys
import time

def f(a,b):
    return a/b

i=3


while(True):
    try:
        print(f"i={i}, f={f(11,i)}"); sys.stdout.flush()
    except Exception  as e:
        print(e)
        print(f"i={i}")
        time.sleep(5)
        i = 5
        continue
    i -= 1


print("After loop")
