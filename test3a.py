def f1(a):
    return a*a*a

while(True):
    z = 0
    for i in range(1000):
        z += f1(i)
    print(z)
